"""
wall_to_wall_climate_m1.py -- Model 1 (remote sensing) wall-to-wall climate,
ONE three-band GeoTIFF per calendar year on the fixed grid:

    1 precipGPM   GPM IMERG monthly precipitation, summed over the year (mm)
    2 LST         MODIS daytime land surface temperature, annual mean (deg C)
    3 airT_TC     TerraClimate air temperature ((tmmx+tmmn)/2), annual mean (deg C)

diffLST_airT is NOT exported: R computes it after download as
LST - airT_TC (Surface minus Air, the locked-in convention). That is the
same arithmetic as before and one less band to ship.

Why this should not hit "Computed value is too large" any more: the old
script reduced three products over ~15 million vector cells in one request.
This one never creates vector cells; each variable is resampled straight
onto a raster grid.

HOW EACH SOURCE IS PUT ON THE GRID
  GPM (0.1 deg) and TerraClimate (1/24 deg) are COARSER than a 4 km cell, so
  they are interpolated (bilinear by default; --resample nearest keeps native
  blocks). MODIS LST (1 km) is FINER than a 4 km cell, so it is area-averaged.

--lst_product MOD11A2 swaps the daily product for the 8-day composite: about
8x less data to read, annual means are very close, but the fitting side used
daily MOD11A1, so keep the default unless the export is too slow.

I could not run this against live Earth Engine when writing it. Use --test first.

YEARS: export EVERY YEAR that falls inside a prediction window. R averages the
yearly rasters (window = label year -5 .. +4, so the 2010 prediction uses 2005-2014
and the extended 2020 prediction uses 2015-2020), then predicts once per decade
label. The default core window therefore needs 2005-2014 (10 files); the extended
window adds 2015-2020. In R, w2w_missing() lists exactly which years are still absent.

Usage:
  python gee/wall_to_wall_climate_m1.py --years 2009 2010 --test        # small pipeline test
  python gee/wall_to_wall_climate_m1.py --years $(seq 2005 2014)        # core window
  python gee/wall_to_wall_climate_m1.py --years $(seq 2015 2020)        # extended window (optional)
"""
import argparse

import ee

from wall_to_wall_common import (add_common_args, init_ee, interpolate_to_grid, make_grid,
                                 reduce_to_grid, stem_m1, submit_export, existing_task_states)

M1_BANDS = ["precipGPM", "LST", "airT_TC"]
MODIS_NATIVE_M = 926.6   # MODIS "1 km" sinusoidal pixel


def annual_precip_gpm(year):
    ic = (ee.ImageCollection("NASA/GPM_L3/IMERG_MONTHLY_V07")
          .filterDate(f"{year}-01-01", f"{year + 1}-01-01"))

    def monthly_total(img):
        date = ee.Date(img.get("system:time_start"))
        hours = date.advance(1, "month").difference(date, "day").multiply(24)
        return img.select("precipitation").multiply(ee.Image.constant(hours)).rename("pr_gpm_month")

    proj = ee.Image(ic.first()).select("precipitation").projection()
    return ic.map(monthly_total).sum().rename("precipGPM").setDefaultProjection(proj)


def annual_air_t(year):
    ic = (ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE")
          .filterDate(f"{year}-01-01", f"{year + 1}-01-01"))

    def monthly_mean(img):
        t = img.select(["tmmx", "tmmn"]).multiply(0.1)
        return t.select("tmmx").add(t.select("tmmn")).divide(2).rename("airT_month")

    proj = ee.Image(ic.first()).select("tmmx").projection()
    return ic.map(monthly_mean).mean().rename("airT_TC").setDefaultProjection(proj)


def annual_lst(year, product):
    ic = (ee.ImageCollection(f"MODIS/061/{product}")
          .filterDate(f"{year}-01-01", f"{year + 1}-01-01").select("LST_Day_1km"))
    proj = ee.Image(ic.first()).select("LST_Day_1km").projection()
    return (ic.mean().multiply(0.02).subtract(273.15).rename("LST").setDefaultProjection(proj))


def build_m1(year, g, lst_product, resample):
    precip_g = interpolate_to_grid(annual_precip_gpm(year), g, resample)
    air_g = interpolate_to_grid(annual_air_t(year), g, resample)
    lst_g = reduce_to_grid(annual_lst(year, lst_product), g,
                           max_pixels=4 * (g.res_m / MODIS_NATIVE_M) ** 2 + 64)
    return ee.Image.cat([precip_g, lst_g, air_g]).select(M1_BANDS)


if __name__ == "__main__":
    ap = add_common_args(argparse.ArgumentParser())
    ap.add_argument("--years", type=int, nargs="+", required=True)
    ap.add_argument("--lst_product", default="MOD11A1", choices=["MOD11A1", "MOD11A2"])
    ap.add_argument("--resample", default="bilinear", choices=["bilinear", "nearest"])
    args = ap.parse_args()

    init_ee()
    if args.check:
        for coll, band in [("NASA/GPM_L3/IMERG_MONTHLY_V07", "precipitation"),
                           ("IDAHO_EPSCOR/TERRACLIMATE", "tmmx"),
                           (f"MODIS/061/{args.lst_product}", "LST_Day_1km")]:
            have = ee.Image(ee.ImageCollection(coll).first()).bandNames().getInfo()
            print(f"{coll:34s} {band:14s}", "OK" if band in have else f"MISSING (has {have})")
        raise SystemExit(0)

    g = make_grid(args.res, args.test)
    states = existing_task_states()
    for y in args.years:
        submit_export(build_m1(y, g, args.lst_product, args.resample),
                      stem_m1(g, y), g, states=states, force=args.force)
    print("bands in file order:", M1_BANDS)
