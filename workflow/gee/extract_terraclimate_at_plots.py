"""
extract_terraclimate_at_plots.py -- Model 2's fitting climate (meanTemp,
winTmin, annualPrec) from TerraClimate, at FIA plot locations, 2001-2020.

WHY THIS EXISTS: Model 2 ("climate-only") is the model that has to survive
into CMIP6 for the historical-baseline and future wall-to-wall stages, so
it is fit ONLY on variables that exist in both TerraClimate (here) and
NASA/GDDP-CMIP6 (see gee/wall_to_wall_climate_m2.py) -- mean/min
temperature and precipitation. This restores the project's ORIGINAL
pre-migration formula definitions exactly:
  meanTemp   = mean of monthly (tmmx + tmmn) / 2, all months, all years
               (NOT the `tas`-only convention used by the now-abandoned
               NEX-GDDP fitting approach -- Model 2 restores the older,
               tmmx/tmmn-based definition on purpose)
  winTmin    = coldest MONTHLY tmmn value, Jan/Feb/Mar only (months 1-3,
               NOT December -- matches WINTER_MONTHS in R/00_config.R and
               the project's original clusterFIA.R convention)
  annualPrec = sum of monthly pr, all months, all years, / number of years

SCALE FACTORS, verified against the GEE catalog page for
IDAHO_EPSCOR/TERRACLIMATE -- do not skip this: tmmx/tmmn are stored with a
0.1 scale factor (raw pixel value x 0.1 = degrees C); pr has no scale
factor (already mm). Both are applied below before any reduction.

COMPUTATIONAL-COMPLEXITY NOTE: unlike the old NEX-GDDP extraction (which
exported 720 raw monthly bands and did the meanTemp/winTmin/annualPrec
reduction in R), this script does the temporal reduction SERVER-SIDE in
GEE and exports only 3 final columns per plot. Much smaller download,
much lighter reduceRegions call, no need for the old script's tileScale=16
tuning.

Prerequisites:
  pip install earthengine-api
  earthengine authenticate

Usage:
  python extract_terraclimate_at_plots.py --asset users/<you>/fia_plot_points \
      --out terraclimate_at_plots
"""
import argparse
import ee

FIT_YEARS     = list(range(2001, 2021))   # 2001-2020, matches R/00_config.R FIT_YEARS
WINTER_MONTHS = [1, 2, 3]                 # matches WINTER_MONTHS in R/00_config.R

def build_climate_image() -> ee.Image:
    tc = (ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE")
            .filterDate(f"{FIT_YEARS[0]}-01-01", f"{FIT_YEARS[-1] + 1}-01-01"))

    # apply TerraClimate's documented scale factors: tmmx/tmmn x0.1 -> degrees
    # C; pr has no scale factor (already mm).
    def scale(img):
        temp = img.select(["tmmx", "tmmn"]).multiply(0.1)
        return temp.addBands(img.select("pr")).copyProperties(img, ["system:time_start"])
    tc_scaled = tc.map(scale)

    # meanTemp: mean of monthly (tmmx+tmmn)/2, ALL months, ALL years -- the
    # ORIGINAL clusterFIA.R definition, restored (not the tas-only
    # convention the abandoned NEX-GDDP approach used).
    def monthly_mean_temp(img):
        return img.select("tmmx").add(img.select("tmmn")).divide(2).rename("meanTemp_month")
    meanTemp = tc_scaled.map(monthly_mean_temp).mean().rename("meanTemp")

    # winTmin: coldest MONTHLY tmmn, Jan/Feb/Mar only, across all years.
    tc_winter = tc_scaled.filter(ee.Filter.calendarRange(WINTER_MONTHS[0], WINTER_MONTHS[-1], "month"))
    winTmin = tc_winter.select("tmmn").min().rename("winTmin")

    # annualPrec: sum of monthly pr, all months, all years, / number of years.
    annualPrec = tc_scaled.select("pr").sum().divide(len(FIT_YEARS)).rename("annualPrec")

    return meanTemp.addBands(winTmin).addBands(annualPrec)

def extract(plots_fc: ee.FeatureCollection, out_name: str) -> ee.batch.Task:
    img = build_climate_image()
    sampled = img.reduceRegions(
        collection=plots_fc,
        reducer=ee.Reducer.mean(),   # trivial at a point -- one value already
        scale=4638,                  # TerraClimate native resolution (~4 km)
        tileScale=4,                 # only 3 output bands -- far lighter than
                                     # the old 720-band NEX-GDDP export
    )
    task = ee.batch.Export.table.toDrive(
        collection=sampled,
        description=out_name,
        fileNamePrefix=out_name,
        fileFormat="CSV",
        selectors=["rowid", "meanTemp", "winTmin", "annualPrec"],
    )
    task.start()
    return task

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--asset", required=True,
                    help="EE FeatureCollection asset ID of FIA plot points "
                         "(property 'rowid') -- from R/01d_export_plot_points.R")
    ap.add_argument("--out", default="terraclimate_at_plots")
    args = ap.parse_args()

    ee.Initialize()
    fc = ee.FeatureCollection(args.asset)
    t = extract(fc, args.out)
    print(f"submitted export task: {t.status()['id']}")
    print(f"result will appear at Google Drive/{args.out}.csv (3 columns + rowid --")
    print("small enough this should finish quickly; no need for a subset test run.")
