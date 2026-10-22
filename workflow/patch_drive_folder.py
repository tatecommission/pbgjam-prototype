# Lets the GEE exports go to a different Drive folder via the W2W_DRIVE_FOLDER environment variable.
# Default is unchanged ("pbgjam_w2w"), so 4 km behaviour is identical. Fails loudly if the line is not found exactly once.
p = "gee/wall_to_wall_common.py"
s = open(p).read()
old = 'DRIVE_FOLDER = "pbgjam_w2w"'
new = 'DRIVE_FOLDER = os.environ.get("W2W_DRIVE_FOLDER", "pbgjam_w2w")'
if s.count(old) != 1:
    raise SystemExit(f"expected exactly one match for {old!r}, found {s.count(old)} -- nothing changed")
open(p, "w").write(s.replace(old, new))
print("patched OK")
