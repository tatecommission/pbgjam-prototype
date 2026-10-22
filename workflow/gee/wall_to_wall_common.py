"""
wall_to_wall_common.py -- shared pieces for the raster wall-to-wall pipeline
(wall_to_wall_static.py, wall_to_wall_climate_m1.py, wall_to_wall_climate_m2.py).

WHY RASTERS INSTEAD OF CSV TABLES
  The previous version exported one CSV row per 1 km cell of the whole
  lon/lat bounding box (~15.3 million cells, including ocean, Canada and
  Mexico) via reduceRegions over a FeatureCollection. That is what produced
  "Computed value is too large", multi-GB CSVs, and runs that never finished.
  Here every product is exported as a GeoTIFF on ONE fixed EPSG:5070 grid, so
  all layers line up cell-for-cell and are joined in R by cell index.

THE GRID (identical for every export)
  CRS EPSG:5070 (CONUS Albers, same as NLCD). Top-left corner
  (-2,500,000, 3,320,000) m, bottom-right (2,356,000, 168,000) m. Every edge
  is a multiple of 4 km, so 1, 2 and 4 km grids all nest exactly. Verified
  to cover CONUS with >= 97 km of margin on every side. At 4 km this is
  1214 x 788 = 956,632 cells (about 3.8 MB per Float32 band); at 1 km it is
  4856 x 3152 = 15,306,112 cells (about 61 MB per band).

NO ee CALLS AT IMPORT TIME. (Creating an ee.Geometry before ee.Initialize()
is what crashed the previous versions of these scripts.)
"""
import csv
import datetime
import os
from dataclasses import dataclass

import ee

CRS = "EPSG:5070"
X0, Y0 = -2_500_000, 3_320_000      # top-left corner of the full grid (m)
X1, Y1 = 2_356_000, 168_000         # bottom-right corner (m)
# Small North Carolina window (Blue Ridge -> Piedmont -> coast) for cheap
# end-to-end tests: xmin, ymin, xmax, ymax in metres. 960 x 480 km.
TEST_BOX = (1_000_000, 1_280_000, 1_960_000, 1_760_000)

DEFAULT_RES_M = 4000
ALLOWED_RES_M = (1000, 2000, 4000)
NODATA = -9999.0                    # sentinel written into every GeoTIFF
DRIVE_FOLDER = os.environ.get("W2W_DRIVE_FOLDER", "pbgjam_w2w")         # all exports land in this Drive folder
MANIFEST = os.path.join(os.path.dirname(os.path.abspath(__file__)), "task_manifest.csv")
ACTIVE_STATES = {"READY", "PENDING", "RUNNING", "COMPLETED"}


def init_ee():
    """Initialize Earth Engine with an explicit project (no more 'no project
    found' errors in a fresh terminal)."""
    ee.Initialize(project=os.environ.get("GOOGLE_CLOUD_PROJECT") or "pbgjam-v2")


@dataclass
class Grid:
    res_m: int
    test: bool
    xmin: int
    ymin: int
    xmax: int
    ymax: int

    @property
    def transform(self):
        return [self.res_m, 0, self.xmin, 0, -self.res_m, self.ymax]

    @property
    def ncol(self):
        return (self.xmax - self.xmin) // self.res_m

    @property
    def nrow(self):
        return (self.ymax - self.ymin) // self.res_m

    @property
    def tag(self):
        return f"{self.res_m // 1000}km"

    @property
    def suffix(self):
        return "_test" if self.test else ""

    def subgrid_transform(self, res_m):
        """Transform of a finer/coarser lattice with the same origin (used for
        the intermediate topography grid)."""
        return [res_m, 0, self.xmin, 0, -res_m, self.ymax]

    def region(self):
        """Only call after init_ee()."""
        return ee.Geometry.Rectangle([self.xmin, self.ymin, self.xmax, self.ymax],
                                     proj=CRS, geodesic=False)


def make_grid(res_m=DEFAULT_RES_M, test=False) -> Grid:
    if res_m not in ALLOWED_RES_M:
        raise SystemExit(f"--res must be one of {ALLOWED_RES_M} so the grid nests exactly")
    xmin, ymin, xmax, ymax = TEST_BOX if test else (X0, Y1, X1, Y0)
    g = Grid(res_m, test, xmin, ymin, xmax, ymax)
    for edge in (xmin, ymin, xmax, ymax):
        if edge % res_m:
            raise SystemExit(f"grid edge {edge} is not a multiple of {res_m} m")
    return g


# ---- file-name stems (mirrored in R/w2w_helpers.R -- keep in sync) -----------
def stem_static(g: Grid):
    return f"gee_static_conus_{g.tag}{g.suffix}"


def stem_m1(g: Grid, year):
    return f"gee_climate_m1_{g.tag}_{year}{g.suffix}"


def stem_m2(g: Grid, label, scenario, year):
    return f"gee_climate_m2_{label}_{g.tag}_{scenario}_{year}{g.suffix}"


# ---- resampling helpers ---------------------------------------------------------
def reduce_to_grid(image, g: Grid, max_pixels):
    """Area-average a FINER image into grid cells (masked pixels are ignored; a
    cell with no valid pixels stays masked). The image must have a fixed
    projection (assets do; reprojected/setDefaultProjection images do)."""
    return (image.reduceResolution(reducer=ee.Reducer.mean(), maxPixels=int(max_pixels))
                 .reproject(crs=CRS, crsTransform=g.transform))


def interpolate_to_grid(image, g: Grid, mode="bilinear"):
    """Resample a COARSER image (GCM/GPM/TerraClimate pixels bigger than the
    grid cells) onto the grid. bilinear avoids 25 km blocks in the maps; use
    'nearest' to keep the native blocks. The image must carry its native
    projection (setDefaultProjection)."""
    img = image if mode == "nearest" else image.resample(mode)
    return img.reproject(crs=CRS, crsTransform=g.transform)


def pick_shard(res_m, finest_m, target_px=2560):
    """Largest export tile (output px per side) whose intermediate render stays
    under EE's reprojection limit (~10,000 px). Intermediate = shard * res_m/finest_m."""
    ratio = res_m / finest_m
    for s in (256, 128, 64, 32, 16):
        if s * ratio <= target_px:
            return s
    return 16


# ---- export guard + manifest ------------------------------------------------------
def existing_task_states():
    """{description: set(states)} from Earth Engine's task history. History
    expires after a while, so this guards against duplicate submissions, it is
    not a record of what is in Drive."""
    states = {}
    try:
        for t in ee.batch.Task.list():
            desc = (getattr(t, "config", None) or {}).get("description")
            st = getattr(t, "state", None)
            if desc and st:
                states.setdefault(desc, set()).add(st)
    except Exception as exc:  # listing is a convenience; never block a submit on it
        print(f"WARNING: could not read task history ({exc}); duplicate guard is off")
    return states


def log_manifest(name, task_id, g: Grid):
    new = not os.path.exists(MANIFEST)
    with open(MANIFEST, "a", newline="") as fh:
        w = csv.writer(fh)
        if new:
            w.writerow(["submitted_utc", "export_name", "task_id", "res_m", "test"])
        w.writerow([datetime.datetime.utcnow().isoformat(timespec="seconds"),
                    name, task_id, g.res_m, int(g.test)])


def submit_export(image, name, g: Grid, states=None, force=False, shard_size=None):
    """Export `image` as a Float32 GeoTIFF on the fixed grid into the Drive
    folder. Masked pixels are written as NODATA (-9999) because Earth Engine
    otherwise writes them as 0, which is a real value for several predictors.
    Skips a name that is already queued, running or completed unless force."""
    states = existing_task_states() if states is None else states
    active = states.get(name, set()) & ACTIVE_STATES
    if active and not force:
        print(f"skip {name}: already {sorted(active)} (use --force to resubmit)")
        return None
    img = image.toFloat().unmask(NODATA)
    opts = {"shardSize": int(shard_size)} if shard_size else {}
    task = ee.batch.Export.image.toDrive(
        image=img, description=name, folder=DRIVE_FOLDER, fileNamePrefix=name,
        region=g.region(), crs=CRS, crsTransform=g.transform,
        maxPixels=int(1e10), fileFormat="GeoTIFF", **opts)
    task.start()
    tid = task.status().get("id")
    log_manifest(name, tid, g)
    print(f"submitted {name}  ({g.ncol} x {g.nrow} cells): {tid}")
    return task


def add_common_args(ap):
    ap.add_argument("--res", type=int, default=DEFAULT_RES_M, choices=ALLOWED_RES_M,
                    help="grid cell size in metres (default 4000). Must match "
                         "GRID_RES_M in R/00_config.R")
    ap.add_argument("--test", action="store_true",
                    help="export only the small North Carolina test window "
                         "(files get a _test suffix) -- run this first")
    ap.add_argument("--force", action="store_true",
                    help="resubmit even if a task with this name is queued/running/completed")
    ap.add_argument("--check", action="store_true",
                    help="only check that the assets/bands this script needs exist, then exit")
    return ap
