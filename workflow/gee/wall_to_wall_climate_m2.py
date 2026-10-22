"""
wall_to_wall_climate_m2.py -- Model 2 (CMIP6) wall-to-wall climate, ONE
GeoTIFF per (scenario, year) on the fixed grid, holding the equal-weight mean
of the model ensemble for each requested variable:

    winTmin     minimum over Jan-Mar of the monthly mean of tasmin (deg C)
    annualPrec  annual sum of daily pr x 86400 (mm)
    (meanTemp   optional, only if requested with --vars)

DEFAULT BANDS ARE winTmin AND annualPrec ONLY. The locked-in Model 2 formula
(M2_06) does not use meanTemp, and dropping it removes 12 monthly means per
model from every export task. If the formula ever changes to use meanTemp,
re-export with:  --vars winTmin annualPrec meanTemp  (and update R accordingly).

ENSEMBLE: 13 models (13MME), reduced from the 16 in Wu et al. (2023) because
CESM2, CESM2-WACCM and IITM-ESM have no tasmax/tasmin bands in
NASA/GDDP-CMIP6 (confirmed by listing bandNames() on the live catalog). That
is a data-availability constraint, not a modelling choice -- flagged for Tong.
Equal-weight mean; the paper's own weighting was not recoverable.

FAIL FAST: before submitting anything, every model is checked client-side for
images and the needed bands in each requested scenario. The previous version
only discovered CESM2-WACCM's missing band after the task had been queued.

WHY THIS IS CHEAP NOW: GDDP-CMIP6 is 0.25 deg (~25 km). Exporting it on a
1 km vector grid repeated each value ~600 times. Here the ensemble is
computed on the native grid and interpolated (bilinear by default) to the
4 km raster, so each task only has to produce ~1 million cells.

YEARS: export EVERY YEAR that falls inside a prediction window. R averages the
yearly rasters (window = label year -5 .. +4, clipped to the scenario's years) and
predicts once per decade label (1950, 1960, ..., 2100). That needs 1950-2014 for
"historical" and 2015-2100 for each ssp scenario: 65 + 2 x 86 = 237 small exports.
In R, w2w_missing() lists exactly which years are still absent.

Usage:
  python gee/wall_to_wall_climate_m2.py --check --years 2010 --scenarios historical
  python gee/wall_to_wall_climate_m2.py --years 2009 2010 --scenarios historical --test   # pipeline test
  python gee/wall_to_wall_climate_m2.py --years $(seq 1950 2014) --scenarios historical
  python gee/wall_to_wall_climate_m2.py --years $(seq 2015 2100) --scenarios ssp245
  python gee/wall_to_wall_climate_m2.py --years $(seq 2015 2100) --scenarios ssp585
"""
import argparse

import ee

from wall_to_wall_common import (add_common_args, init_ee, interpolate_to_grid, make_grid,
                                 stem_m2, submit_export, existing_task_states)

COLLECTION = "NASA/GDDP-CMIP6"
WINTER_MONTHS = [1, 2, 3]
# 13MME -- see header for why this is not the original 16.
DEFAULT_MODELS = [
    "ACCESS-CM2", "ACCESS-ESM1-5", "BCC-CSM2-MR",
    "CMCC-CM2-SR5", "CMCC-ESM2", "CNRM-CM6-1",
    "MIROC6",
    "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0",
    "NorESM2-LM", "NorESM2-MM", "TaiESM1",
]
NEEDED_BANDS = {"winTmin": ["tasmin"], "annualPrec": ["pr"], "meanTemp": ["tasmax", "tasmin"]}


def base_collection(year, scenario, models):
    return (ee.ImageCollection(COLLECTION)
            .filter(ee.Filter.eq("scenario", scenario))
            .filter(ee.Filter.calendarRange(year, year, "year"))
            .filter(ee.Filter.inList("model", models)))


def validate_models(models, scenario, year, variables):
    """Raise before anything is submitted if a model lacks images or bands."""
    need = sorted({b for v in variables for b in NEEDED_BANDS[v]})
    base = (ee.ImageCollection(COLLECTION).filter(ee.Filter.eq("scenario", scenario))
            .filter(ee.Filter.calendarRange(year, year, "year")))
    bad = {}
    for m in models:
        ic = base.filter(ee.Filter.eq("model", m))
        if ic.size().getInfo() == 0:
            bad[m] = f"no {scenario} images for {year}"
            continue
        # check EVERY image of the year, not just the first: a model is unusable only if
        # no image has the bands; a few defective days are tolerated (and reported)
        flags = (ic.map(lambda im: im.set("ok", ee.Algorithms.If(im.bandNames().containsAll(need), 1, 0)))
                   .aggregate_array("ok").getInfo())
        n_ok = sum(flags)
        if n_ok == 0:
            bad[m] = f"no image has bands {need}"
        elif n_ok < len(flags):
            print(f"  WARNING {m} {scenario} {year}: {len(flags) - n_ok} of {len(flags)} images lack {need} (skipped)")
    if bad:
        lines = "\n".join(f"  {m}: {why}" for m, why in bad.items())
        raise SystemExit(f"ERROR: {len(bad)} model(s) unusable for {scenario} {year}:\n{lines}\n"
                         "Remove them with --models (and change --ensemble_label to match).")
    print(f"models OK for {scenario} {year}: {len(models)} models have {need}")


def one_model(ic, variables):
    """Per-model value of each requested variable (same definitions as the
    plot-level fit side)."""
    bands = []
    if "winTmin" in variables:
        # a few catalog images lack tasmin (e.g. ssp585 TaiESM1 20930202); skip just those days
        ic_t = (ic.map(lambda im: im.set("hasTmin", ee.Algorithms.If(im.bandNames().contains("tasmin"), 1, 0)))
                  .filter(ee.Filter.eq("hasTmin", 1)))
        monthly = [ic_t.filter(ee.Filter.calendarRange(m, m, "month")).select("tasmin").mean().subtract(273.15)
                   for m in WINTER_MONTHS]
        bands.append(ee.Image.cat(monthly).reduce(ee.Reducer.min()).rename("winTmin"))
    if "annualPrec" in variables:
        bands.append(ic.select("pr").map(lambda im: im.multiply(86400)).sum().rename("annualPrec"))
    if "meanTemp" in variables:
        mt = []
        for m in range(1, 13):
            mo = ic.filter(ee.Filter.calendarRange(m, m, "month"))
            mt.append(mo.select("tasmax").mean().subtract(273.15)
                        .add(mo.select("tasmin").mean().subtract(273.15)).divide(2))
        bands.append(ee.Image.cat(mt).reduce(ee.Reducer.mean()).rename("meanTemp"))
    return ee.Image.cat(bands)


def build_m2(year, scenario, models, variables, g, resample):
    if scenario == "historical" and year > 2014:
        raise SystemExit(f"scenario=historical has no data after 2014 (got {year}); use ssp245/ssp585")
    base = base_collection(year, scenario, models)
    per_model = [one_model(base.filter(ee.Filter.eq("model", m)), variables) for m in models]
    ens = ee.ImageCollection.fromImages(per_model).mean()
    proj = ee.Image(base.first()).select(NEEDED_BANDS[variables[0]][0]).projection()
    ens = ens.setDefaultProjection(proj)          # compute on the native 0.25 deg grid
    return interpolate_to_grid(ens, g, resample).select(variables)


if __name__ == "__main__":
    ap = add_common_args(argparse.ArgumentParser())
    ap.add_argument("--years", type=int, nargs="+", required=True)
    ap.add_argument("--scenarios", nargs="+", default=["historical"])
    ap.add_argument("--models", nargs="+", default=DEFAULT_MODELS)
    ap.add_argument("--ensemble_label", default="13MME",
                    help="goes into file names; must match GCM_ENSEMBLE_LABEL in R/00_config.R")
    ap.add_argument("--vars", nargs="+", default=["winTmin", "annualPrec"],
                    choices=list(NEEDED_BANDS))
    ap.add_argument("--resample", default="bilinear", choices=["bilinear", "nearest"])
    args = ap.parse_args()

    init_ee()
    for s in args.scenarios:
        validate_models(args.models, s, args.years[0], args.vars)
    if args.check:
        raise SystemExit(0)

    g = make_grid(args.res, args.test)
    states = existing_task_states()
    for s in args.scenarios:
        for y in args.years:
            submit_export(build_m2(y, s, args.models, args.vars, g, args.resample),
                          stem_m2(g, args.ensemble_label, s, y), g, states=states, force=args.force)
    print("bands in file order:", args.vars)
