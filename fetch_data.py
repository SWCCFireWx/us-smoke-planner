#!/usr/bin/env python3
"""
National Smoke Dispersion & Burn-Window Planner - data builder.

Runs inside GitHub Actions (.github/workflows/update-data.yml) and writes the
site's data folder:

  manifest.json                 run time, valid times, fields, domains, time zones
  <domain>/<field>/<NNN>.png    one 8-bit grayscale image per field per forecast time
  <domain>/blk/<row>_<col>.png  "point" files: every field and time for a 16x16 block
  <domain>/tz.png               index of the local time zone of every cell
  gacc.geojson, psa.geojson, perimeters.geojson     NIFC layers
  alerts.geojson                NWS Red Flag Warnings and Fire Weather Watches
  airnow.geojson                latest hourly PM2.5 at AirNow monitors
  class1.geojson                Mandatory Class I federal areas (EPA)
  nonattain.geojson             PM2.5 nonattainment areas (EPA Green Book)
  erc.geojson                   ERC / BI percentiles by PSA (NIFC), as published
  meta.json                     when each source was fetched and whether it worked

Forecast sources
  * NOAA National Blend of Models (NBM) "core" GRIB2 for CONUS, Alaska, Hawaii and
    Puerto Rico: mixing height, transport wind, 10-m wind/gust (+ spread), 2-m
    temperature and RH, PoP and QPF. Out to 7 days (hourly to 48 h, 3-hourly to
    96 h, 6-hourly to 168 h) when the NBM run reaches that far.
  * NOAA HRRR (CONUS): near-surface and vertically integrated smoke, out to 48 h.
Only the needed GRIB messages are downloaded (HTTP byte ranges from .idx files).

Requires: numpy, scipy, pygrib, timezonefinder  (installed by the workflow).
Usage:   python fetch_data.py --out _site/data [--domains conus,ak,hi,pr]
"""
import argparse
import collections
import concurrent.futures as cf
import datetime as dt
import json
import math
import os
import re
import struct
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import warnings
import zlib

import numpy as np

warnings.filterwarnings("ignore", category=RuntimeWarning)

# ---------------------------------------------------------------------------
# Domains. Longitudes run west-negative and may go below -180 so the western
# Aleutians (east of the date line) stay next to the rest of Alaska.
#   dlat/dlon = output grid spacing (deg); k = NBM points averaged per cell
# ---------------------------------------------------------------------------
DOMAINS = collections.OrderedDict([
    ("conus", {"file": "co", "label": "CONUS", "west": -125.0, "east": -66.5, "south": 24.0, "north": 50.0,
               "dlat": 0.125, "dlon": 0.125, "k": 24}),
    ("ak",    {"file": "ak", "label": "Alaska", "west": -188.0, "east": -129.0, "south": 51.0, "north": 71.5,
               "dlat": 0.1, "dlon": 0.25, "k": 12}),
    ("hi",    {"file": "hi", "label": "Hawaii", "west": -160.5, "east": -154.5, "south": 18.75, "north": 22.5,
               "dlat": 0.05, "dlon": 0.05, "k": 4}),
    ("pr",    {"file": "pr", "label": "Puerto Rico / USVI", "west": -67.5, "east": -64.5, "south": 17.6, "north": 18.6,
               "dlat": 0.025, "dlon": 0.025, "k": 4}),
])
BLOCK = 16                                            # cells per side of a point block

BUCKET = "https://noaa-nbm-grib2-pds.s3.amazonaws.com"
FHOURS = list(range(1, 49)) + list(range(51, 97, 3)) + list(range(102, 169, 6))
HRRR_BUCKET = "https://noaa-hrrr-bdp-pds.s3.amazonaws.com"

NIFC = "https://services3.arcgis.com/T4QMspbfLg3qTGWY/ArcGIS/rest/services"
GACC_URL = NIFC + "/DMP_NationalGACCBoundaries_Public/FeatureServer/0/query"
PSA_URL = NIFC + "/DMP_Predictive_Service_Area__PSA_Boundaries_Public/FeatureServer/0/query"
PERIM_URL = NIFC + "/WFIGS_Interagency_Perimeters_Current/FeatureServer/0/query"
ERC_URLS = [NIFC + f"/PSA_ERC_and_BI_Percentiles_and_Trends/FeatureServer/{i}/query" for i in (0, 1, 2)]
WWA_URL = "https://mapservices.weather.noaa.gov/eventdriven/rest/services/WWA/watch_warn_adv/MapServer/1/query"
EPA = "https://services.arcgis.com/cJ9YHowT8TU7DUyn/ArcGIS/rest/services"
CLASS1_URL = EPA + "/Mandatory_Class1_Federal_Areas/FeatureServer/0/query"
NAA_URLS = [(EPA + "/Nonattainment_Areas_and_Designations/FeatureServer/5/query", "PM2.5 24-hr (2006)"),
            (EPA + "/Nonattainment_Areas_and_Designations/FeatureServer/7/query", "PM2.5 annual (2012)")]
AIRNOW_BASES = ["https://files.airnowtech.org/airnow", "https://s3-us-west-1.amazonaws.com/files.airnowtech.org/airnow"]

UA = "National-Smoke-Dispersion-Planner (GitHub Actions; non-commercial)"
M_TO_FT = 3.28084
MS_TO_KT = 1.943844
MS_TO_MPH = 2.236936
M2S_TO_KTFT = MS_TO_KT * M_TO_FT                      # 1 m2/s = 6.3775 kt-ft
MISSING = 255

# Output fields in a fixed order: name -> (scale, offset, units, encoding).
#   lin:  pixel = round((v - offset) / scale)        v = pixel * scale + offset
#   sqrt: pixel = round(sqrt(v) * scale)             v = (pixel / scale) ** 2
FIELDS = collections.OrderedDict([
    ("vi",   (2000.0, 0.0,  "kt-ft",  "lin")),   # ventilation rate (to 508,000)
    ("mh",   (100.0,  0.0,  "ft AGL", "lin")),   # mixing height
    ("tw",   (1.0,    0.0,  "kt",     "lin")),   # transport wind speed
    ("twd",  (2.0,    0.0,  "deg",    "lin")),   # transport wind direction (from)
    ("pop",  (1.0,    0.0,  "%",      "lin")),   # chance of >=0.01 in, period ending at this time
    ("qpf",  (0.01,   0.0,  "in",     "lin")),   # precip amount, period ending at this time
    ("ws",   (1.0,    0.0,  "mph",    "lin")),   # 10-m wind speed
    ("wd",   (2.0,    0.0,  "deg",    "lin")),   # 10-m wind direction (from)
    ("wg",   (1.0,    0.0,  "mph",    "lin")),   # 10-m gust
    ("t",    (1.0,  -60.0,  "F",      "lin")),   # 2-m temperature
    ("rh",   (1.0,    0.0,  "%",      "lin")),   # 2-m relative humidity
    ("vsd",  (2000.0, 0.0,  "kt-ft",  "lin")),   # NBM spread (1 std dev) of ventilation
    ("wsd",  (0.5,    0.0,  "mph",    "lin")),   # NBM spread (1 std dev) of 10-m wind
    ("smk",  (10.0,   0.0,  "ug/m3",  "sqrt")),  # HRRR near-surface smoke (to ~645 ug/m3)
    ("csmk", (4.0,    0.0,  "mg/m2",  "sqrt")),  # HRRR vertically integrated smoke (to ~4,000 mg/m2)
])


# ---------------------------------------------------------------------------
# HTTP helpers
# ---------------------------------------------------------------------------
def http_get(url, params=None, headers=None, retries=5, timeout=120):
    full = url + ("?" + urllib.parse.urlencode(params) if params else "")
    last = None
    for attempt in range(1, retries + 1):
        try:
            h = {"User-Agent": UA}
            h.update(headers or {})
            with urllib.request.urlopen(urllib.request.Request(full, headers=h), timeout=timeout) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            last = e
            if e.code in (403, 404):
                raise
            if e.code in (429, 500, 502, 503, 504):
                time.sleep(min(60, 5 * attempt))
                continue
            raise
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            last = e
            time.sleep(5 * attempt)
    raise RuntimeError(f"giving up on {url}: {last}")


def http_get_json(url, params):
    return json.loads(http_get(url, params).decode("utf-8"))


def exists(url):
    try:
        req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=30) as r:
            return r.status == 200
    except Exception:
        return False


# ---------------------------------------------------------------------------
# File helpers
# ---------------------------------------------------------------------------
def write_json(out, name, obj):
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, name)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, separators=(",", ":"))
    print(f"  wrote {name} ({os.path.getsize(path) / 1024:.0f} KB)", flush=True)


def write_png_gray(path, arr):
    """Write a uint8 2-D array as an 8-bit grayscale PNG (PNG 'Up' filter for small files)."""
    arr = np.ascontiguousarray(arr, dtype=np.uint8)
    h, w = arr.shape
    up = np.empty_like(arr)
    up[0] = arr[0]
    up[1:] = (arr[1:].astype(np.int16) - arr[:-1].astype(np.int16)) & 0xFF
    raw = np.concatenate([np.full((h, 1), 2, np.uint8), up], axis=1).tobytes()

    def chunk(tag, data):
        c = struct.pack(">I", len(data)) + tag + data
        return c + struct.pack(">I", zlib.crc32(tag + data) & 0xFFFFFFFF)

    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as f:
        f.write(b"\x89PNG\r\n\x1a\n"
                + chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 0, 0, 0, 0))
                + chunk(b"IDAT", zlib.compress(raw, 6))
                + chunk(b"IEND", b""))


def quantize(vals, field):
    scale, offset, _, enc = FIELDS[field]
    if enc == "sqrt":
        q = np.round(np.sqrt(np.maximum(vals, 0)) * scale)
    else:
        q = np.round((vals - offset) / scale)
    return np.where(np.isfinite(q), np.clip(q, 0, MISSING - 1), MISSING).astype(np.uint8)


# ---------------------------------------------------------------------------
# ArcGIS layers
# ---------------------------------------------------------------------------
def arcgis_geojson(url, extra, page=True):
    feats, offset = [], 0
    while True:
        params = {"where": "1=1", "outFields": "*", "returnGeometry": "true", "outSR": "4326", "f": "geojson"}
        if page:
            params.update({"resultOffset": offset, "resultRecordCount": 1000})
        params.update(extra)
        j = http_get_json(url, params)
        if "error" in j:
            raise RuntimeError(f"ArcGIS error: {j['error']}")
        batch = j.get("features") or []
        feats.extend(batch)
        exceeded = j.get("exceededTransferLimit") or (j.get("properties") or {}).get("exceededTransferLimit")
        if not page or not exceeded or not batch:
            break
        offset += len(batch)
    return feats


def round_coords(geom, nd):
    def r(c):
        return round(c, nd) if isinstance(c, (int, float)) else [r(x) for x in c]
    if geom and "coordinates" in geom:
        geom["coordinates"] = r(geom["coordinates"])
    return geom


def pick(props, *names):
    for n in names:
        v = props.get(n)
        if v not in (None, ""):
            return v
    return None


def name_like(props):
    """First property whose key looks like a name."""
    for k, v in props.items():
        if isinstance(v, str) and v and re.search(r"name", k, re.I):
            return v
    return ""


def fetch_gacc(out):
    print("Fetching GACC boundaries ...", flush=True)
    feats = arcgis_geojson(GACC_URL, {"outFields": "GACCName,GACCAbbreviation,GACCUnitID",
                                      "maxAllowableOffset": "0.01", "geometryPrecision": "3"})
    res = [{"type": "Feature", "geometry": round_coords(f.get("geometry"), 3), "properties": {
        "name": pick(f.get("properties") or {}, "GACCName") or "",
        "abbr": pick(f.get("properties") or {}, "GACCAbbreviation", "GACCUnitID") or ""}}
        for f in feats if f.get("geometry")]
    if not res:
        raise RuntimeError("no GACC features returned")
    write_json(out, "gacc.geojson", {"type": "FeatureCollection", "features": res})
    return len(res)


def fetch_psa(out):
    print("Fetching Predictive Service Areas ...", flush=True)
    feats = arcgis_geojson(PSA_URL, {"maxAllowableOffset": "0.01", "geometryPrecision": "3"})
    res = []
    for f in feats:
        p = f.get("properties") or {}
        if not f.get("geometry"):
            continue
        res.append({"type": "Feature", "geometry": round_coords(f["geometry"], 3), "properties": {
            "code": pick(p, "PSANationalCode", "PSANATIONALCODE", "PSA_ID", "PSAID") or "",
            "name": pick(p, "PSANAME", "PSAName", "PSA_NAME") or "",
            "gacc": pick(p, "GACCUnitID", "GACCName") or ""}})
    write_json(out, "psa.geojson", {"type": "FeatureCollection", "features": res})
    return len(res)


def fetch_perimeters(out):
    print("Fetching WFIGS current perimeters ...", flush=True)
    feats = arcgis_geojson(PERIM_URL, {"maxAllowableOffset": "0.001", "geometryPrecision": "4"})
    res = []
    for f in feats:
        p = f.get("properties") or {}
        if not f.get("geometry"):
            continue
        typ = (pick(p, "attr_IncidentTypeCategory", "IncidentTypeCategory", "poly_IncidentTypeCategory") or "").upper()
        res.append({"type": "Feature", "geometry": round_coords(f["geometry"], 4), "properties": {
            "name": pick(p, "poly_IncidentName", "attr_IncidentName", "IncidentName") or "Unnamed",
            "type": typ or "WF",
            "acres": pick(p, "poly_GISAcres", "attr_IncidentSize", "GISAcres"),
            "contained": pick(p, "attr_PercentContained", "PercentContained"),
            "discovered": pick(p, "attr_FireDiscoveryDateTime", "FireDiscoveryDateTime"),
            "updated": pick(p, "poly_DateCurrent", "poly_PolygonDateTime", "attr_ModifiedOnDateTime_dt"),
            "state": pick(p, "attr_POOState", "POOState"),
            "id": pick(p, "attr_UniqueFireIdentifier", "UniqueFireIdentifier", "poly_IRWINID"),
        }})
    print(f"  {len(res)} perimeters ({sum(1 for x in res if x['properties']['type'] == 'RX')} prescribed)")
    write_json(out, "perimeters.geojson", {"type": "FeatureCollection", "features": res})
    return len(res)


def parse_time(v):
    """Alert times arrive as epoch ms, ISO text or compact text; return epoch seconds or None."""
    if v in (None, ""):
        return None
    if isinstance(v, (int, float)):
        return int(v / 1000) if v > 1e11 else int(v)
    s = str(v).strip()
    for fmt in ("%Y-%m-%dT%H:%M:%S%z", "%Y-%m-%dT%H:%M%z", "%Y-%m-%d %H:%M:%S%z", "%Y-%m-%dT%H:%M:%SZ",
                "%Y-%m-%d %H:%M:%S", "%Y-%m-%dT%H:%M:%S", "%Y%m%d%H%M", "%m/%d/%Y %H:%M:%S", "%m/%d/%Y %H:%M"):
        try:
            d = dt.datetime.strptime(s.replace("Z", "+0000") if "%z" in fmt else s, fmt)
            if d.tzinfo is None:
                d = d.replace(tzinfo=dt.timezone.utc)
            return int(d.timestamp())
        except ValueError:
            continue
    try:
        d = dt.datetime.fromisoformat(s.replace("Z", "+00:00"))
        return int((d if d.tzinfo else d.replace(tzinfo=dt.timezone.utc)).timestamp())
    except ValueError:
        return None


def fetch_alerts(out):
    print("Fetching NWS Red Flag Warnings / Fire Weather Watches ...", flush=True)
    feats = arcgis_geojson(WWA_URL, {
        "where": "prod_type IN ('Red Flag Warning','Fire Weather Watch')",
        "outFields": "prod_type,onset,ends,expiration,issuance,wfo,url,event",
        "maxAllowableOffset": "0.005", "geometryPrecision": "3"}, page=False)
    res = []
    for f in feats:
        p = f.get("properties") or {}
        if not f.get("geometry"):
            continue
        start = parse_time(p.get("onset")) or parse_time(p.get("issuance"))
        end = parse_time(p.get("ends")) or parse_time(p.get("expiration"))
        res.append({"type": "Feature", "geometry": round_coords(f["geometry"], 3), "properties": {
            "type": p.get("prod_type") or "", "start": start, "end": end,
            "wfo": p.get("wfo") or "", "url": p.get("url") or ""}})
    write_json(out, "alerts.geojson", {"type": "FeatureCollection", "features": res})
    return len(res)


def fetch_class1(out):
    print("Fetching Class I areas ...", flush=True)
    feats = arcgis_geojson(CLASS1_URL, {"outFields": "NAME,STATE,AGENCY", "maxAllowableOffset": "0.005",
                                        "geometryPrecision": "3"})
    res = [{"type": "Feature", "geometry": round_coords(f["geometry"], 3), "properties": {
        "name": pick(f.get("properties") or {}, "NAME") or "Class I area",
        "state": pick(f.get("properties") or {}, "STATE") or "",
        "agency": pick(f.get("properties") or {}, "AGENCY") or ""}} for f in feats if f.get("geometry")]
    write_json(out, "class1.geojson", {"type": "FeatureCollection", "features": res})
    return len(res)


def fetch_nonattainment(out):
    print("Fetching PM2.5 nonattainment areas ...", flush=True)
    res = []
    for url, label in NAA_URLS:
        try:
            feats = arcgis_geojson(url, {"maxAllowableOffset": "0.005", "geometryPrecision": "3"})
        except Exception as e:
            print(f"  {label}: {e}", flush=True)
            continue
        for f in feats:
            if not f.get("geometry"):
                continue
            p = f.get("properties") or {}
            status = next((v for k, v in p.items() if isinstance(v, str) and re.search(r"status|class", k, re.I)), "")
            res.append({"type": "Feature", "geometry": round_coords(f["geometry"], 3),
                        "properties": {"name": name_like(p) or "Nonattainment area", "standard": label, "status": status}})
    if not res:
        raise RuntimeError("no nonattainment features")
    write_json(out, "nonattain.geojson", {"type": "FeatureCollection", "features": res})
    return len(res)


def fetch_erc(out):
    """ERC/BI percentiles by PSA as NIFC publishes them; field names are passed through for the page to detect."""
    print("Fetching ERC / BI percentiles by PSA ...", flush=True)
    last = None
    for url in ERC_URLS:
        try:
            feats = arcgis_geojson(url, {"maxAllowableOffset": "0.01", "geometryPrecision": "3"})
        except Exception as e:
            last = e
            continue
        if not feats:
            continue
        res = []
        for f in feats:
            p = {k: v for k, v in (f.get("properties") or {}).items()
                 if isinstance(v, (int, float, str)) and k.lower() not in ("shape__area", "shape__length", "globalid")}
            res.append({"type": "Feature", "geometry": round_coords(f.get("geometry"), 3) if f.get("geometry") else None,
                        "properties": p})
        print(f"  {len(res)} features from {url.split('/services/')[1]}; fields: {', '.join(list(res[0]['properties'])[:40])}")
        write_json(out, "erc.geojson", {"type": "FeatureCollection", "features": res, "source": url})
        return len(res)
    raise RuntimeError(f"ERC/BI layer not reachable: {last}")


def fetch_airnow(out):
    """Latest hourly PM2.5 at AirNow monitors (public files, no key)."""
    print("Fetching AirNow PM2.5 ...", flush=True)

    def get_text(path):
        for base in AIRNOW_BASES:
            try:
                return http_get(f"{base}/{path}", retries=2, timeout=60).decode("latin-1")
            except Exception:
                continue
        return None

    sites_txt = get_text("today/monitoring_site_locations.dat")
    if not sites_txt:
        raise RuntimeError("AirNow site file not reachable")
    sites = {}
    for line in sites_txt.splitlines():
        c = line.split("|")
        if len(c) >= 10 and c[1].strip().upper() in ("PM2.5", "PM25"):
            try:
                sites[c[0].strip()] = (float(c[8]), float(c[9]), c[3].strip(), c[6].strip())
            except ValueError:
                pass
    now = dt.datetime.now(dt.timezone.utc).replace(minute=0, second=0, microsecond=0)
    for back in range(1, 5):
        h = now - dt.timedelta(hours=back)
        txt = get_text(f"today/HourlyData_{h:%Y%m%d%H}.dat") or get_text(f"{h:%Y}/{h:%Y%m%d}/HourlyData_{h:%Y%m%d%H}.dat")
        if not txt:
            continue
        res = []
        for line in txt.splitlines():
            c = line.split("|")
            if len(c) < 8 or c[5].strip().upper() not in ("PM2.5", "PM25"):
                continue
            s = sites.get(c[2].strip())
            try:
                val = float(c[7])
            except ValueError:
                continue
            if not s or val < -5:
                continue
            res.append({"type": "Feature", "geometry": {"type": "Point", "coordinates": [round(s[1], 4), round(s[0], 4)]},
                        "properties": {"name": s[2], "agency": s[3], "pm25": round(val, 1)}})
        if res:
            write_json(out, "airnow.geojson", {"type": "FeatureCollection", "time": int(h.timestamp()), "features": res})
            return len(res)
    raise RuntimeError("no recent AirNow hourly file")


# ---------------------------------------------------------------------------
# GRIB helpers
# ---------------------------------------------------------------------------
def parse_idx(text):
    """Return (start, end_or_None, line) byte ranges from a GRIB .idx file."""
    rows = []
    for line in text.splitlines():
        parts = line.split(":")
        if len(parts) > 2 and parts[1].isdigit():
            rows.append((int(parts[1]), line))
    return [(s, rows[i + 1][0] - 1 if i + 1 < len(rows) else None, l) for i, (s, l) in enumerate(rows)]


def fetch_range(url, start, end):
    return http_get(url, headers={"Range": f"bytes={start}-" + ("" if end is None else str(end))}, timeout=180)


def decode_grib(msg, latlons=False):
    """Decode one GRIB2 message -> (2-D float32 values with NaN, lats, lons). Replaced in tests."""
    import pygrib
    g = pygrib.fromstring(msg)
    vals = np.ma.filled(np.ma.masked_invalid(g.values).astype(np.float32), np.nan)
    if latlons:
        lats, lons = g.latlons()
        return vals, lats, lons
    return vals, None, None


# ---------------------------------------------------------------------------
# NBM message selection
# ---------------------------------------------------------------------------
def nbm_url(cycle, fh, ftype):
    return f"{BUCKET}/blend.{cycle:%Y%m%d}/{cycle:%H}/core/blend.t{cycle:%H}z.core.f{fh:03d}.{ftype}.grib2"


ACC_RE = re.compile(r":APCP:surface:(\d+)-(\d+) hour acc fcst:")


def _plain(line):
    return "std dev" not in line and "prob" not in line


def _sd(line):
    return "std dev" in line and "prob" not in line


def select_messages(idx_rows, fh, prev_fh):
    """Pick the byte range for each needed NBM field in one forecast file."""
    want = {
        "MIXHT": lambda l: ":MIXHT:" in l and _plain(l),
        "TRWSPD": lambda l: ":TRWSPD:" in l and _plain(l),
        "TRWDIR": lambda l: ":TRWDIR:" in l and _plain(l),
        "WIND": lambda l: ":WIND:10 m above ground:" in l and _plain(l),
        "WDIR": lambda l: ":WDIR:10 m above ground:" in l and _plain(l),
        "GUST": lambda l: ":GUST:10 m above ground:" in l and _plain(l),
        "TMP": lambda l: ":TMP:2 m above ground:" in l and _plain(l),
        "RH": lambda l: ":RH:2 m above ground:" in l and _plain(l),
        # spread (ensemble standard deviation), used for the confidence layers when NBM provides it
        "VRATE_SD": lambda l: ":VRATE:" in l and _sd(l),
        "MIXHT_SD": lambda l: ":MIXHT:" in l and _sd(l),
        "TRWSPD_SD": lambda l: ":TRWSPD:" in l and _sd(l),
        "WIND_SD": lambda l: ":WIND:10 m above ground:" in l and _sd(l),
    }
    inst = f":{fh} hour fcst:"
    sel = {}
    for key, ok in want.items():
        cands = [(s, e, l) for s, e, l in idx_rows if ok(l)]
        if cands:   # prefer the record valid at exactly this hour
            exact = [c for c in cands if inst in c[2]]
            sel[key] = (exact or cands)[0]
    gap = fh - prev_fh
    for key, is_prob in (("POP", True), ("QPF", False)):
        best = None
        for s, e, l in idx_rows:
            m = ACC_RE.search(l)
            if not m or int(m.group(2)) != fh or "std dev" in l:
                continue
            if is_prob != (">0.254" in l) or (not is_prob and "prob" in l):
                continue
            win = int(m.group(2)) - int(m.group(1))
            score = (0 if win == gap else 1, win)   # match the time step first, then the shortest window
            if best is None or score < best[0]:
                best = (score, (s, e, l), win)
        if best:
            sel[key] = best[1]
            sel[key + "_WIN"] = best[2]
    return sel


# ---------------------------------------------------------------------------
# Regridding to a regular lat/lon grid (nearest k points on the sphere, averaged)
# ---------------------------------------------------------------------------
def _xyz(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.column_stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])


def domain_shape(c):
    nx = int(round((c["east"] - c["west"]) / c["dlon"])) + 1
    ny = int(round((c["north"] - c["south"]) / c["dlat"])) + 1
    return ny, nx


def domain_latlon(c):
    ny, nx = domain_shape(c)
    lat = c["north"] - np.arange(ny) * c["dlat"]     # row 0 = north
    lon = c["west"] + np.arange(nx) * c["dlon"]
    LON, LAT = np.meshgrid(lon, lat)
    return LAT, LON


class Regridder:
    def __init__(self, c, lats, lons, k=None):
        from scipy.spatial import cKDTree
        k = k or c["k"]
        lonw = np.mod(lons, 360.0) - 360.0               # west-negative, may go below -180
        mask = ((lats >= c["south"] - 0.6) & (lats <= c["north"] + 0.6) &
                (lonw >= c["west"] - 1.5) & (lonw <= c["east"] + 1.5))
        self.mask = mask.ravel()
        if not self.mask.any():
            raise RuntimeError("source grid does not cover this domain")
        self.ny, self.nx = domain_shape(c)
        LAT, LON = domain_latlon(c)
        tree = cKDTree(_xyz(lats.ravel()[self.mask], lonw.ravel()[self.mask]))
        dist, self.idx = tree.query(_xyz(LAT.ravel(), LON.ravel()), k=k)
        if k == 1:
            dist, self.idx = dist[:, None], self.idx[:, None]
        midlat = np.radians((c["north"] + c["south"]) / 2)
        cell = np.radians(max(c["dlat"], c["dlon"] * np.cos(midlat)))
        self.outside = dist[:, 0] > 2.0 * cell            # beyond the edge of the source grid

    def scalar(self, vals):
        out = np.nanmean(vals.ravel()[self.mask][self.idx], axis=1)
        out[self.outside] = np.nan
        return out.reshape(self.ny, self.nx)

    def vector(self, speed, direction):
        """Average wind as u/v; returns (mean speed, vector-mean direction it blows FROM)."""
        s = speed.ravel()[self.mask][self.idx]
        d = np.radians(direction.ravel()[self.mask][self.idx])
        u = np.nanmean(-s * np.sin(d), axis=1)
        v = np.nanmean(-s * np.cos(d), axis=1)
        spd = np.nanmean(s, axis=1)
        dirn = (np.degrees(np.arctan2(-u, -v)) + 360) % 360
        spd[self.outside] = np.nan
        dirn[self.outside] = np.nan
        return spd.reshape(self.ny, self.nx), dirn.reshape(self.ny, self.nx)


# ---------------------------------------------------------------------------
# Local time zone of every cell
# ---------------------------------------------------------------------------
def _fallback_zone(dkey, lat, lon):
    if dkey == "ak":
        return "America/Adak" if lon < -169.5 else "America/Anchorage"
    if dkey == "hi":
        return "Pacific/Honolulu"
    if dkey == "pr":
        return "America/Puerto_Rico"
    if -114.8 <= lon <= -109.05 and 31.3 <= lat <= 37.0:
        return "America/Phoenix"
    if lon < -114.5:
        return "America/Los_Angeles"
    if lon < -101.5:
        return "America/Denver"
    if lon < -86.5:
        return "America/Chicago"
    return "America/New_York"


_TF = None


def _zone_finder():
    global _TF
    if _TF is None:
        try:
            from timezonefinder import TimezoneFinder
            try:
                _TF = TimezoneFinder(in_memory=True)
            except TypeError:
                _TF = TimezoneFinder()
        except Exception as e:  # keep going with a rough rule if the package is unavailable
            print(f"  timezonefinder unavailable ({e}); using approximate zones", flush=True)
            _TF = False
    return _TF


def zone_grid(dkey, c, zones):
    """uint8 grid of indexes into `zones` (shared list of IANA names, grown here).
    Looks up every other cell, then only the cells between lookups that disagree."""
    tf = _zone_finder()
    LAT, LON = domain_latlon(c)
    ny, nx = LAT.shape
    lookup = {z: i for i, z in enumerate(zones)}

    def zone_at(r, col):
        lat, lon = float(LAT[r, col]), float(LON[r, col])
        z = None
        if tf:
            try:
                z = tf.timezone_at(lng=((lon + 180) % 360) - 180, lat=lat)
            except Exception:
                z = None
        z = z or _fallback_zone(dkey, lat, lon)
        if z not in lookup:
            if len(zones) >= 254:
                z = _fallback_zone(dkey, lat, lon)
            if z not in lookup:
                lookup[z] = len(zones)
                zones.append(z)
        return lookup[z]

    out = np.full((ny, nx), 255, np.uint8)
    rows, cols = list(range(0, ny, 2)), list(range(0, nx, 2))
    if rows[-1] != ny - 1:
        rows.append(ny - 1)
    if cols[-1] != nx - 1:
        cols.append(nx - 1)
    rowset, colset = set(rows), set(cols)
    for r in rows:
        for col in cols:
            out[r, col] = zone_at(r, col)
    ri = np.searchsorted(rows, np.arange(ny))
    ci = np.searchsorted(cols, np.arange(nx))
    for r in range(ny):
        r1 = rows[min(ri[r], len(rows) - 1)]
        r0 = r1 if r in rowset else rows[max(ri[r] - 1, 0)]
        for col in range(nx):
            if out[r, col] != 255:
                continue
            c1 = cols[min(ci[col], len(cols) - 1)]
            c0 = c1 if col in colset else cols[max(ci[col] - 1, 0)]
            corners = {out[r0, c0], out[r0, c1], out[r1, c0], out[r1, c1]}
            out[r, col] = corners.pop() if len(corners) == 1 else zone_at(r, col)
    return out


# ---------------------------------------------------------------------------
# NBM forecast build
# ---------------------------------------------------------------------------
def find_cycle(now=None):
    """Newest NBM run that reaches 7 days; otherwise the newest that reaches 96 h."""
    now = now or dt.datetime.now(dt.timezone.utc).replace(minute=0, second=0, microsecond=0)
    fallback = None
    for back in range(1, 14):
        cyc = now - dt.timedelta(hours=back)
        if not exists(nbm_url(cyc, FHOURS[0], "co") + ".idx"):
            continue
        if exists(nbm_url(cyc, FHOURS[-1], "co") + ".idx"):
            return cyc
        if fallback is None and exists(nbm_url(cyc, 96, "co") + ".idx"):
            fallback = cyc
    if fallback:
        print("  no recent 7-day NBM run; using the newest 4-day run", flush=True)
        return fallback
    raise RuntimeError("no complete NBM cycle found in the last 13 hours")


def get_idx(url):
    try:
        return parse_idx(http_get(url + ".idx").decode("utf-8", "ignore"))
    except urllib.error.HTTPError as e:
        if e.code in (403, 404):
            return None
        raise


def download(url, sel):
    return {k: fetch_range(url, s, e) for k, (s, e, _) in ((k, v) for k, v in sel.items() if not k.endswith("_WIN"))}


def add_hrrr_smoke(c, times, Q):
    """Fill Q['smk'] and Q['csmk'] for CONUS from the newest HRRR run that covers the first 48 h."""
    now = dt.datetime.now(dt.timezone.utc).replace(minute=0, second=0, microsecond=0)
    url = lambda cy, f: f"{HRRR_BUCKET}/hrrr.{cy:%Y%m%d}/conus/hrrr.t{cy:%H}z.wrfsfcf{f:02d}.grib2"
    cyc, last, short = None, 0, None
    for back in range(1, 16):           # prefer the newest 48-h run (00/06/12/18Z); else the newest 18-h run
        cy = now - dt.timedelta(hours=back)
        if cy.hour % 6 == 0 and exists(url(cy, 48) + ".idx"):
            cyc, last = cy, 48
            break
        if short is None and exists(url(cy, 18) + ".idx"):
            short = cy
    if cyc is None and short is not None:
        cyc, last = short, 18
    if cyc is None:
        raise RuntimeError("no recent HRRR run found")
    print(f"[conus] HRRR smoke from {cyc:%Y-%m-%d %HZ} (to {last} h)", flush=True)
    jobs = []
    for i, e in enumerate(times):
        fh = int(round((e - cyc.timestamp()) / 3600))
        if 0 <= fh <= last:
            jobs.append((i, fh))

    def one(job):
        i, fh = job
        rows = get_idx(url(cyc, fh))
        if not rows:
            return i, {}
        sel = {}
        for s, e, l in rows:
            if ":MASSDEN:8 m above ground:" in l and "MASSDEN" not in sel:
                sel["MASSDEN"] = (s, e, l)
            if ":COLMD:entire atmosphere" in l and "COLMD" not in sel:
                sel["COLMD"] = (s, e, l)
        return i, download(url(cyc, fh), sel)

    rg = None
    n = 0
    with cf.ThreadPoolExecutor(max_workers=8) as ex:
        for i, raw in ex.map(one, jobs):
            for key, msg in raw.items():
                vals, lats, lons = decode_grib(msg, latlons=rg is None)
                if rg is None:
                    rg = Regridder(c, lats, lons, k=16)
                if key == "MASSDEN":
                    Q["smk"][i] = quantize(rg.scalar(vals) * 1e9, "smk")      # kg/m3 -> ug/m3
                else:
                    Q["csmk"][i] = quantize(rg.scalar(vals) * 1e6, "csmk")    # kg/m2 -> mg/m2
            n += bool(raw)
    print(f"[conus] HRRR smoke: {n} times", flush=True)
    return {"cycle": int(cyc.timestamp()), "times": n}


def build_domain(dkey, c, cycle, times_fh, out, zones):
    """Fetch, regrid and write one domain. times_fh = global list of forecast hours (None = decide here)."""
    t0 = time.time()
    print(f"[{dkey}] reading inventories ...", flush=True)
    wanted = times_fh or FHOURS
    with cf.ThreadPoolExecutor(max_workers=8) as ex:
        idxs = {fh: rows for fh, rows in zip(wanted, ex.map(lambda f: get_idx(nbm_url(cycle, f, c["file"])), wanted)) if rows}
    hours = sorted(idxs)
    if not hours:
        raise RuntimeError(f"no NBM {c['label']} files in cycle {cycle:%Y%m%d %HZ}")
    times_fh = times_fh or hours
    gidx = {fh: i for i, fh in enumerate(times_fh)}
    missing = [f for f in times_fh if f not in idxs]
    if missing:
        print(f"[{dkey}] forecast hours not available (left blank): {missing}", flush=True)
    sels = {fh: select_messages(idxs[fh], fh, hours[i - 1] if i else 0) for i, fh in enumerate(hours)}
    sd_found = sorted({k for s in sels.values() for k in s if k.endswith("_SD")})
    print(f"[{dkey}] NBM spread fields found: {', '.join(sd_found) or 'none'}", flush=True)

    ny, nx = domain_shape(c)
    nt = len(times_fh)
    Q = {f: np.full((nt, ny, nx), MISSING, np.uint8) for f in FIELDS}
    avail = [False] * nt
    qwin = [None] * nt
    rg = None

    def process(fh, raw, win):
        nonlocal rg
        dec = {}
        for key, msg in raw.items():
            vals, lats, lons = decode_grib(msg, latlons=rg is None)
            if rg is None:
                rg = Regridder(c, lats, lons)
            dec[key] = vals
        if rg is None:
            return
        nan = np.full((ny, nx), np.nan, np.float32)
        g = lambda k: rg.scalar(dec[k]) if k in dec else nan
        mh_m = g("MIXHT")
        mh = mh_m * M_TO_FT
        tw, twd = rg.vector(dec["TRWSPD"], dec["TRWDIR"]) if "TRWSPD" in dec and "TRWDIR" in dec else (nan, nan)
        tw_ms = tw
        tw = tw * MS_TO_KT
        ws, wd = rg.vector(dec["WIND"], dec["WDIR"]) if "WIND" in dec and "WDIR" in dec else (nan, nan)
        vi = mh * tw
        if "VRATE_SD" in dec:
            vsd = g("VRATE_SD") * M2S_TO_KTFT
        elif "MIXHT_SD" in dec and "TRWSPD_SD" in dec:   # spread of a product, from relative spreads
            vsd = np.abs(vi) * np.sqrt((g("MIXHT_SD") / np.maximum(mh_m, 1)) ** 2 + (g("TRWSPD_SD") / np.maximum(tw_ms, 0.1)) ** 2)
        else:
            vsd = nan
        vals = {"vi": vi, "mh": mh, "tw": tw, "twd": twd, "pop": g("POP"), "qpf": g("QPF") / 25.4,
                "ws": ws * MS_TO_MPH, "wd": wd, "wg": g("GUST") * MS_TO_MPH,
                "t": (g("TMP") - 273.15) * 9 / 5 + 32, "rh": g("RH"),
                "vsd": vsd, "wsd": g("WIND_SD") * MS_TO_MPH}
        i = gidx[fh]
        for f, v in vals.items():
            Q[f][i] = quantize(v, f)
        avail[i] = True
        qwin[i] = win

    with cf.ThreadPoolExecutor(max_workers=8) as ex:
        dq = collections.deque()
        for fh in hours:
            url = nbm_url(cycle, fh, c["file"])
            dq.append((fh, sels[fh].get("QPF_WIN"), ex.submit(download, url, sels[fh])))
            if len(dq) >= 8:
                fh0, w0, fut = dq.popleft()
                process(fh0, fut.result(), w0)
        while dq:
            fh0, w0, fut = dq.popleft()
            process(fh0, fut.result(), w0)
    print(f"[{dkey}] decoded {sum(avail)} times in {time.time() - t0:.0f} s", flush=True)

    extra = {}
    if dkey == "conus":
        times = [int((cycle + dt.timedelta(hours=h)).timestamp()) for h in times_fh]
        try:
            extra["smoke"] = add_hrrr_smoke(c, times, Q)
        except Exception as e:
            print(f"!! HRRR smoke failed: {e}", file=sys.stderr, flush=True)

    ddir = os.path.join(out, dkey)
    fav = {}
    for f in FIELDS:
        bits = []
        for i in range(nt):
            ok = avail[i] and bool((Q[f][i] != MISSING).any())
            if ok:
                write_png_gray(os.path.join(ddir, f, f"{i:03d}.png"), Q[f][i])
            bits.append("1" if ok else "0")
        fav[f] = "".join(bits)

    nf, B = len(FIELDS), BLOCK
    nby, nbx = -(-ny // B), -(-nx // B)
    allq = np.stack([Q[f] for f in FIELDS])
    nblk = 0
    for by in range(nby):
        for bx in range(nbx):
            blk = np.full((nf, nt, B, B), MISSING, np.uint8)
            sub = allq[:, :, by * B:(by + 1) * B, bx * B:(bx + 1) * B]
            blk[:, :, :sub.shape[2], :sub.shape[3]] = sub
            if (blk[0] == MISSING).all():
                continue
            write_png_gray(os.path.join(ddir, "blk", f"{by}_{bx}.png"), blk.reshape(nf * nt, B * B))
            nblk += 1

    print(f"[{dkey}] computing local time zones ...", flush=True)
    tz = zone_grid(dkey, c, zones)
    write_png_gray(os.path.join(ddir, "tz.png"), tz)
    present = [f for f in FIELDS if (Q[f] != MISSING).any()]
    print(f"[{dkey}] done: {ny}x{nx} grid, {sum(avail)} times, {nblk} blocks, fields {present}, {time.time() - t0:.0f} s", flush=True)
    return times_fh, dict({
        "label": c["label"], "west": c["west"], "east": c["east"], "south": c["south"], "north": c["north"],
        "dlat": c["dlat"], "dlon": c["dlon"], "nx": nx, "ny": ny, "avail": avail, "qpf_window_h": qwin,
        "zones": sorted(set(int(z) for z in np.unique(tz))), "fields": present, "fav": fav,
    }, **extra)


def fetch_forecast(out, domain_keys):
    cycle = find_cycle()
    print(f"Using NBM cycle {cycle:%Y-%m-%d %HZ}", flush=True)
    zones, doms, times_fh = [], {}, None
    for dkey in domain_keys:
        try:
            times_fh, doms[dkey] = build_domain(dkey, DOMAINS[dkey], cycle, times_fh, out, zones)
        except Exception as e:
            if dkey == "conus":
                raise            # CONUS sets the time axis; without it there is no site
            print(f"!! {dkey} failed: {e}", file=sys.stderr, flush=True)
    write_json(out, "manifest.json", {
        "generated": int(time.time()),
        "source": "NOAA National Blend of Models (NBM) core; HRRR smoke",
        "cycle": cycle.strftime("%Y-%m-%dT%H:00Z"),
        "cycle_epoch": int(cycle.timestamp()),
        "times": [int((cycle + dt.timedelta(hours=h)).timestamp()) for h in times_fh],
        "fields": {f: {"scale": s, "offset": o, "units": u, "enc": e} for f, (s, o, u, e) in FIELDS.items()},
        "field_order": list(FIELDS),
        "block": BLOCK,
        "missing": MISSING,
        "zones": zones,
        "domains": doms,
    })
    return len(doms)


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="_site/data")
    ap.add_argument("--domains", default=",".join(DOMAINS))
    args = ap.parse_args()
    out = args.out
    keys = [k for k in args.domains.split(",") if k in DOMAINS]
    if "conus" in keys:
        keys = ["conus"] + [k for k in keys if k != "conus"]
    os.makedirs(out, exist_ok=True)
    meta = {"sources": {}}

    def run(key, fn, label, required=False):
        try:
            n = fn(out)
            meta["sources"][key] = {"label": label, "updated": int(time.time()), "count": n, "ok": True}
        except Exception as e:
            print(f"!! {key} failed: {e}", file=sys.stderr, flush=True)
            meta["sources"][key] = {"label": label, "ok": False, "error": str(e)[:300]}
            if required:
                raise
            write_json(out, f"{key}.geojson", {"type": "FeatureCollection", "features": []})

    run("gacc", fetch_gacc, "NIFC National GACC Boundaries")
    run("psa", fetch_psa, "NIFC Predictive Service Areas")
    run("perimeters", fetch_perimeters, "NIFC WFIGS Current Interagency Fire Perimeters")
    run("alerts", fetch_alerts, "NWS Red Flag Warnings / Fire Weather Watches")
    run("class1", fetch_class1, "EPA Mandatory Class I Federal Areas")
    run("nonattain", fetch_nonattainment, "EPA PM2.5 Nonattainment Areas")
    run("erc", fetch_erc, "NIFC ERC / BI percentiles by PSA")
    run("airnow", fetch_airnow, "AirNow hourly PM2.5")
    # The forecast is required: if NBM CONUS fails the job fails and the previous site stays live.
    run("forecast", lambda o: fetch_forecast(o, keys), "NOAA NBM + HRRR smoke", required=True)
    meta["last_run"] = int(time.time())
    write_json(out, "meta.json", meta)


if __name__ == "__main__":
    main()
