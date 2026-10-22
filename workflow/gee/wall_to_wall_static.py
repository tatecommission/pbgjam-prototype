"""
wall_to_wall_static.py -- static wall-to-wall layers on the fixed EPSG:5070
grid (see wall_to_wall_common.py), as ONE seven-band GeoTIFF:

    1 north_component   mean of sin(slope)*cos(aspect) over land pixels
    2 east_component    mean of sin(slope)*sin(aspect) over land pixels
    3 cec30             SoilGrids CEC, mean of 0-5, 5-15, 15-30 cm
    4 clay30            SoilGrids clay, same depths
    5 sand30            SoilGrids sand, same depths
    6 forest_frac       fraction of the cell that is forest (see below)
    7 land_frac         fraction of the cell that is NLCD land (valid and not open water)

ECOLOGICAL / EFFICIENCY CHOICES (each is a decision to confirm with Tong)

 1. Domain: WHOLE CONUS LAND (changed from forest-only). Topography is averaged
    over all NLCD land pixels (--topo_mask land, the default), so every land
    cell has values. land_frac (NLCD valid and not open water) defines which
    cells are US land; R keeps cells with land_frac >= LAND_FRAC_MIN. Predicting
    outside forest is extrapolation beyond the fitted plots -- forest_frac is
    still exported so maps can be masked or the run restricted to forest later
    (FOREST_FRAC_MIN in R; --topo_mask forest restores forest-weighted
    topography). "Forest" = NLCD deciduous (41), evergreen (42), mixed (43)
    and woody wetland (90) -- woody wetland because FIA counts bottomland and
    swamp forest as forest.

 2. Topography is computed at --topo_scale (default 100 m), not 1 m or 10 m.
    FIA slope/aspect are field measurements on a plot of a few tens of
    metres, so ~100 m hillslope aspect is the right ecological scale, and
    Earth Engine can build it from the DEM's coarse pyramid levels instead
    of reading every 10 m pixel. 100 m is chosen because it divides 1, 2 and 4 km
    exactly. The north/east components are then averaged over
    LAND pixels (or FOREST pixels with --topo_mask forest). Both terms
    enter both formulas linearly, so a cell mean is the
    correct summary for the linear predictor.

 3. DEM = USGS/3DEP/10m, a seamless CONUS product. The previous
    USGS/3DEP/1m is a patchwork of project tiles and does not cover CONUS.

 4. Soil: SoilGrids at 250 m, area-averaged into each cell, reading the same
    three layers as before (unweighted mean). UNITS: SoilGrids stores CEC and
    clay/sand at 10x the conventional values (mmol(c)/kg and g/kg). They are
    multiplied by SOIL_SCALE = 0.1 here, giving cmol(c)/kg and percent, because
    the check on the North Carolina test window showed the fit-side cec30,
    clay30 and sand30 are in those units (grid ranges after /10 match the
    fitted ranges). R/07a_build_static_table.R still compares grid and fit.

TOPOGRAPHY UNIT NOTE (unchanged): slope here is converted from DEGREES
(ee.Terrain.slope). On the fitting side slope is FIA percent grade and uses
atan(pct/100). Do not swap the two conversions.

I could not run this against live Earth Engine when writing it. Run
`--check` first (seconds), then `--test` (minutes), then the full export.

Usage:
  python gee/wall_to_wall_static.py --check
  python gee/wall_to_wall_static.py --test
  python gee/wall_to_wall_static.py            # full CONUS at 4 km, whole-US land domain
  python gee/wall_to_wall_static.py --topo_mask forest   # old forest-weighted topography
"""
import argparse
import math

import ee

from wall_to_wall_common import (add_common_args, init_ee, make_grid, reduce_to_grid,
                                 stem_static, submit_export, pick_shard, CRS)

DEM_ID = "USGS/3DEP/10m"
DEM_BAND = "elevation"
NLCD_ID = "USGS/NLCD_RELEASES/2019_REL/NLCD"   # most recent image in the collection is used
FOREST_CLASSES = [41, 42, 43, 90]
WATER_CLASS = 11
SOIL_ASSETS = {"cec": "projects/soilgrids-isric/cec_mean",
               "clay": "projects/soilgrids-isric/clay_mean",
               "sand": "projects/soilgrids-isric/sand_mean"}
SOIL_DEPTHS = ["0-5cm", "5-15cm", "15-30cm"]
SOIL_NATIVE_M = 250
SOIL_SCALE = 0.1      # SoilGrids mapped units -> cmol(c)/kg (cec) and % (clay, sand)
STATIC_BANDS = ["north_component", "east_component", "cec30", "clay30", "sand30",
                "forest_frac", "land_frac"]


def check_assets():
    print("DEM      ", DEM_ID, ee.Image(DEM_ID).bandNames().getInfo())
    ic = ee.ImageCollection(NLCD_ID)
    print("NLCD     ", NLCD_ID, "images:", ic.aggregate_array("system:index").getInfo())
    print("         bands of latest:", ee.Image(ic.sort("system:time_start", False).first()).bandNames().getInfo())
    for k, aid in SOIL_ASSETS.items():
        bands = ee.Image(aid).bandNames().getInfo()
        need = [f"{k}_{d}_mean" for d in SOIL_DEPTHS]
        print(f"SoilGrids {k:5s}", aid, "OK" if all(b in bands for b in need) else f"MISSING {need} in {bands}")


def nlcd_landcover():
    ic = ee.ImageCollection(NLCD_ID)
    return ee.Image(ic.sort("system:time_start", False).first()).select("landcover")


def forest_binary():
    lc = nlcd_landcover()
    f = lc.eq(FOREST_CLASSES[0])
    for c in FOREST_CLASSES[1:]:
        f = f.Or(lc.eq(c))
    return f.unmask(0).rename("forest")     # outside the NLCD footprint = not forest


def land_binary():
    """1 where NLCD has a valid class other than open water (11); 0 over water
    and outside the NLCD (US) footprint."""
    return nlcd_landcover().neq(WATER_CLASS).unmask(0).rename("land")


def top30(asset, name):
    img = ee.Image(asset)
    picks = [f"{name}_{d}_mean" for d in SOIL_DEPTHS]
    return img.select(picks).reduce(ee.Reducer.mean()).multiply(SOIL_SCALE).rename(name + "30")


def build_static(g, topo_scale, topo_mask="land"):
    t_topo = g.subgrid_transform(topo_scale)
    per_cell_topo = (g.res_m / topo_scale) ** 2

    # forest and land masks at topo_scale (NLCD is categorical, so Earth Engine
    # reads its coarse pyramid level here: a majority-class mask at ~100 m)
    forest_t = forest_binary().reproject(crs=CRS, crsTransform=t_topo)
    land_t = land_binary().reproject(crs=CRS, crsTransform=t_topo)
    mask_t = land_t if topo_mask == "land" else forest_t

    # slope/aspect from the DEM at topo_scale (coarse pyramid levels, cheap)
    dem_t = ee.Image(DEM_ID).select(DEM_BAND).reproject(crs=CRS, crsTransform=t_topo)
    slope_r = ee.Terrain.slope(dem_t).multiply(math.pi / 180)
    aspect_r = ee.Terrain.aspect(dem_t).multiply(math.pi / 180)
    north = slope_r.sin().multiply(aspect_r.cos()).rename("north_component")
    east = slope_r.sin().multiply(aspect_r.sin()).rename("east_component")
    topo = north.addBands(east).updateMask(mask_t)

    topo_g = reduce_to_grid(topo, g, max_pixels=2 * per_cell_topo + 64)
    forest_g = reduce_to_grid(forest_t.rename("forest_frac"), g, max_pixels=2 * per_cell_topo + 64)
    land_g = reduce_to_grid(land_t.rename("land_frac"), g, max_pixels=2 * per_cell_topo + 64)

    soil = top30(SOIL_ASSETS["cec"], "cec").addBands(top30(SOIL_ASSETS["clay"], "clay")) \
                                           .addBands(top30(SOIL_ASSETS["sand"], "sand"))
    soil_g = reduce_to_grid(soil, g, max_pixels=4 * (g.res_m / SOIL_NATIVE_M) ** 2 + 64)

    return ee.Image.cat([topo_g, soil_g, forest_g, land_g]).select(STATIC_BANDS)


if __name__ == "__main__":
    ap = add_common_args(argparse.ArgumentParser())
    ap.add_argument("--topo_scale", type=int, default=100, choices=[50, 100, 200],
                    help="scale (m) at which slope/aspect are computed before cell averaging; "
                         "all choices divide 1, 2 and 4 km exactly")
    ap.add_argument("--topo_mask", default="land", choices=["land", "forest"],
                    help="average slope/aspect over all land pixels (default, whole US) "
                         "or over forest pixels only")
    args = ap.parse_args()

    init_ee()
    if args.check:
        check_assets()
        raise SystemExit(0)

    g = make_grid(args.res, args.test)
    shard = pick_shard(g.res_m, min(args.topo_scale, SOIL_NATIVE_M))
    print(f"export tile size: {shard} px (intermediate {shard * g.res_m / min(args.topo_scale, SOIL_NATIVE_M):.0f} px)")
    submit_export(build_static(g, args.topo_scale, args.topo_mask), stem_static(g), g,
                  force=args.force, shard_size=shard)
    print("bands in file order:", STATIC_BANDS)
    print(f"output: Drive/pbgjam_w2w/{stem_static(g)}.tif")
