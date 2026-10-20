"""Build data/basemap_naturalearth_10m.geojson, the base map of Figure 11, from Natural Earth (public domain).

Run once; analyse.py only reads the GeoJSON it writes. Needs pyshp (pip install pyshp==2.3.1) and the two Natural
Earth 1:10m zips, which are downloaded if they are not in the folder given as the first argument:
  ne_10m_admin_0_countries.zip   https://naciscdn.org/naturalearth/10m/cultural/ne_10m_admin_0_countries.zip
  ne_10m_lakes.zip               https://naciscdn.org/naturalearth/10m/physical/ne_10m_lakes.zip

Rings that do not reach the map window are dropped, the rest are simplified (Douglas-Peucker) and very small islands
and lakes, which would be invisible at the printed scale, are removed.
"""
import io
import json
import sys
import urllib.request
import zipfile
from pathlib import Path

import numpy as np
import shapefile

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "data" / "basemap_naturalearth_10m.geojson"
URL = "https://naciscdn.org/naturalearth/10m/{}/{}.zip"
WINDOW = (92.0, 145.0, -14.0, 9.5)  # lon min, lon max, lat min, lat max (a margin around the plotted extent)
TOL = {"IDN": 0.008, "other": 0.02, "lake": 0.005}  # simplification tolerance, degrees
MIN_AREA = 0.0004  # square degrees (about 5 km2)


def load(folder, theme, name):
    z = Path(folder) / f"{name}.zip" if folder else None
    data = z.read_bytes() if z and z.exists() else urllib.request.urlopen(URL.format(theme, name)).read()
    zf = zipfile.ZipFile(io.BytesIO(data))
    part = {ext: io.BytesIO(zf.read(f"{name}.{ext}")) for ext in ("shp", "shx", "dbf")}
    return shapefile.Reader(shp=part["shp"], shx=part["shx"], dbf=part["dbf"], encoding="utf-8")


def douglas_peucker(pts, tol):
    keep = np.zeros(len(pts), bool)
    keep[[0, -1]] = True
    stack = [(0, len(pts) - 1)]
    while stack:
        a, b = stack.pop()
        if b <= a + 1:
            continue
        seg, rel = pts[b] - pts[a], pts[a + 1:b] - pts[a]
        norm = np.hypot(*seg)
        d = np.abs(seg[0] * rel[:, 1] - seg[1] * rel[:, 0]) / norm if norm else np.hypot(rel[:, 0], rel[:, 1])
        i = int(np.argmax(d))
        if d[i] > tol:
            keep[a + 1 + i] = True
            stack += [(a, a + 1 + i), (a + 1 + i, b)]
    return pts[keep]


def area(pts):
    x, y = pts[:, 0], pts[:, 1]
    return 0.5 * (np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def rings(shape, tol):
    bounds = list(shape.parts) + [len(shape.points)]
    for a, b in zip(bounds[:-1], bounds[1:]):
        pts = np.asarray(shape.points[a:b], float)
        lo, hi = pts.min(0), pts.max(0)
        if hi[0] < WINDOW[0] or lo[0] > WINDOW[1] or hi[1] < WINDOW[2] or lo[1] > WINDOW[3]:
            continue
        if abs(area(pts)) < MIN_AREA:
            continue
        s = douglas_peucker(pts, tol)
        if len(s) >= 4:
            yield np.round(s, 3).tolist(), bool(area(pts) > 0)  # shapefile: outer rings clockwise (negative area)


def main(folder=None):
    features = []
    countries = load(folder, "cultural", "ne_10m_admin_0_countries")
    fields = [f[0] for f in countries.fields[1:]]
    for sr in countries.iterShapeRecords():
        rec = dict(zip(fields, sr.record))
        iso = rec["ADM0_A3"]
        rr = list(rings(sr.shape, TOL["IDN"] if iso == "IDN" else TOL["other"]))
        if rr:
            features.append(dict(type="Feature",
                                 properties=dict(kind="country", iso=iso, name=rec["NAME"], holes=[h for _, h in rr]),
                                 geometry=dict(type="MultiLineString", coordinates=[r for r, _ in rr])))
    lakes = load(folder, "physical", "ne_10m_lakes")
    fields = [f[0] for f in lakes.fields[1:]]
    for sr in lakes.iterShapeRecords():
        rec = dict(zip(fields, sr.record))
        rr = list(rings(sr.shape, TOL["lake"]))
        if rr:
            features.append(dict(type="Feature",
                                 properties=dict(kind="lake", name=rec.get("name") or "", holes=[h for _, h in rr]),
                                 geometry=dict(type="MultiLineString", coordinates=[r for r, _ in rr])))
    gj = dict(type="FeatureCollection", features=features,
              source="Natural Earth 1:10m Admin 0 Countries and Lakes (public domain), naturalearthdata.com; "
                     "rings outside the map window removed, simplified, islands and lakes under 0.0004 sq. deg removed")
    OUT.write_text(json.dumps(gj, separators=(",", ":")))
    n = sum(len(r) for f in features for r in f["geometry"]["coordinates"])
    print(f"wrote {OUT.name}: {len(features)} features, {n} points, {OUT.stat().st_size / 1e6:.2f} MB")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)
