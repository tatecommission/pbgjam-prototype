"""
extract_remote_sensing_at_plots.py -- Model 1's THREE fitting predictors at
FIA plot locations: GPM precipitation, MODIS LST, and TerraClimate air
temperature (the reference for the difference term).

REVISED: VegDRI REMOVED entirely. Tong's call -- VegDRI and diffLST_airT
"measure a similar thing" (both are moisture/thermal-stress proxies), so
Model 1 keeps diffLST_airT only rather than carrying two overlapping
drought terms. This also removes the old sub-window problem: all three
remaining predictors now use the SAME FIT_YEARS window (2001-2020) --
there is no more VegDRI-driven 2009 floor forcing one term onto a shorter
period than the others.

THREE PREDICTORS, all FIT_YEARS (2001-2020):
  precipGPM     GPM IMERG monthly precipitation, summed and averaged to an
                annual rate.
  LST           MODIS MOD11A1 daytime land surface temperature, mean.
  airT_TC       TerraClimate mean air temperature ((tmmx+tmmn)/2), mean --
                exported as its own column (not just folded into the
                difference term) so the two components can be sanity-
                checked separately in R.
  diffLST_airT  LST - airT_TC (Surface minus Air) -- computed here
                directly from the same underlying images, exported
                alongside its two components so R can verify
                diffLST_airT == LST - airT_TC exactly.

SCALE FACTORS (unchanged from the prior version of this script):
  MOD11A1 LST_Day_1km: raw x 0.02 = Kelvin; subtract 273.15 for Celsius.
  TerraClimate tmmx/tmmn: raw x 0.1 = degrees C.
  GPM IMERG_MONTHLY precipitation band: mm/hr, average rate for the month
  -- multiplied by hours-in-that-month below to get a monthly total.

Prerequisites:
  pip install earthengine-api
  earthengine authenticate

Usage:
  python extract_remote_sensing_at_plots.py --asset users/<you>/fia_plot_points \
      --out remote_sensing_at_plots
"""
import argparse
import ee

FIT_YEARS = list(range(2001, 2021))   # 2001-2020

def gpm_precip() -> ee.Image:
    ic = (ee.ImageCollection("NASA/GPM_L3/IMERG_MONTHLY_V07")
            .filterDate(f"{FIT_YEARS[0]}-01-01", f"{FIT_YEARS[-1] + 1}-01-01"))

    def monthly_total(img):
        date = ee.Date(img.get("system:time_start"))
        hours_in_month = date.advance(1, "month").difference(date, "day").multiply(24)
        total = img.select("precipitation").multiply(ee.Image.constant(hours_in_month))
        return total.rename("pr_gpm_month").copyProperties(img, ["system:time_start"])

    return ic.map(monthly_total).sum().divide(len(FIT_YEARS)).rename("precipGPM")

def modis_lst() -> ee.Image:
    ic = (ee.ImageCollection("MODIS/061/MOD11A1")
            .filterDate(f"{FIT_YEARS[0]}-01-01", f"{FIT_YEARS[-1] + 1}-01-01")
            .select("LST_Day_1km"))
    return ic.mean().multiply(0.02).subtract(273.15).rename("LST")

def terraclimate_air_t() -> ee.Image:
    tc = (ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE")
            .filterDate(f"{FIT_YEARS[0]}-01-01", f"{FIT_YEARS[-1] + 1}-01-01"))

    def monthly_mean_temp(img):
        temp = img.select(["tmmx", "tmmn"]).multiply(0.1)
        return temp.select("tmmx").add(temp.select("tmmn")).divide(2).rename("airT_month")

    return tc.map(monthly_mean_temp).mean().rename("airT_TC")

def build_image() -> ee.Image:
    lst  = modis_lst()
    airt = terraclimate_air_t()
    precip = gpm_precip()
    diff = lst.subtract(airt).rename("diffLST_airT")   # Surface minus Air
    return precip.addBands(lst).addBands(airt).addBands(diff)

def extract(plots_fc: ee.FeatureCollection, out_name: str) -> ee.batch.Task:
    img = build_image()
    sampled = img.reduceRegions(
        collection=plots_fc,
        reducer=ee.Reducer.mean(),
        scale=1000,
        tileScale=4,
    )
    task = ee.batch.Export.table.toDrive(
        collection=sampled,
        description=out_name,
        fileNamePrefix=out_name,
        fileFormat="CSV",
        selectors=["rowid", "precipGPM", "LST", "airT_TC", "diffLST_airT"],
    )
    task.start()
    return task

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--asset", required=True,
                    help="EE FeatureCollection asset ID of FIA plot points "
                         "(property 'rowid') -- from R/01d_export_plot_points.R")
    ap.add_argument("--out", default="remote_sensing_at_plots")
    args = ap.parse_args()

    ee.Initialize()
    fc = ee.FeatureCollection(args.asset)
    t = extract(fc, args.out)
    print(f"submitted export task: {t.status()['id']}")
    print(f"result will appear at Google Drive/{args.out}.csv (4 columns + rowid)")
