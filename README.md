# National Smoke Dispersion & Burn-Window Planner

A free web page, hosted on GitHub Pages, for smoke dispersion and burn-window planning
anywhere in the United States, built for Predictive Services at every GACC.

Live site (once set up): **https://swccfirewx.github.io/us-smoke-planner/**

## What it shows

**Forecast (NOAA National Blend of Models, NBM)** on four grids: CONUS (2.5 km), Alaska (3 km,
including the Aleutians across the date line), Hawaii, and Puerto Rico / USVI. Out to **7 days**:
hourly to 48 h, 3-hourly to 4 days, 6-hourly to 7 days (when the NBM run reaches that far).
Smooth, interpolated maps at about 10–14 km.

| Map layer | What it is |
|---|---|
| Ventilation rate | Mixing height × transport wind, in categories (default NWS 40k/60k/100k/150k kt-ft; editable) |
| Burn-window hours (selected day) | Hours meeting *your* criteria, in each location's own local time |
| Meets burn criteria (this hour) | Yes / no at the selected time |
| Confidence: ventilation / wind | Chance ventilation meets your minimum, or 10-m wind stays under your max, from NBM's ensemble spread |
| Overnight smoke settling risk | Screening index for the night starting that evening (see below) |
| Near-surface smoke / vertically integrated smoke | NOAA HRRR smoke (CONUS, first 48 h) |
| Mixing height, transport wind, precip chance, 24-h precip, 10-m wind, gusts, RH, temperature | NBM |

**Overlays**
- **NWS Red Flag Warnings and Fire Weather Watches.** Checked live from NWS when the page opens, with a
  cached copy as backup. Alerts in effect at the selected time are drawn bold, others faint. By default,
  hours under a warning or watch are excluded from burn windows.
- **AirNow PM2.5 monitors** (latest hourly value, colored by EPA PM2.5 breakpoints; not the 24-hour AQI)
- **Class I areas** and **PM2.5 nonattainment areas** (EPA)
- **PSAs colored by ERC percentile** (NIFC ERC/BI percentiles by PSA)
- GACCs, PSAs, current wildfire and prescribed-fire perimeters (NIFC)

**Click any point** for charts (dispersion with NBM spread, precip and wind, HRRR smoke), burn windows
in local time + UTC + your main time, and flags for:
- Red Flag Warnings / Fire Weather Watches at or near the window
- smoke drifting toward a Class I or nonattainment area within 40 miles (from the window's mean transport wind)
- the PSA's ERC and BI percentiles
- the nearest AirNow monitor

There's a daily table (including overnight settling risk) and a **CSV download** of the full point forecast.

**Other tools**
- **Area menu**: Lower 48, any GACC, Alaska, Hawaii, Puerto Rico. Sets the usual time zones for that area.
- **Clocks**: main zone, optional second zone, and UTC.
- **Copy link**: the URL captures area, layer, time, clicked point, criteria, ventilation scale and units.
- **Print / PDF**: a landscape page with title, times, criteria and legend (use "Save as PDF" in the print dialog).
- **My burn units**: save units with their own prescriptions; each shows its next window. Export/import
  as a JSON file to share with others; download all units' windows as CSV.
- **Display**: ventilation in kt-ft, m²/s or mph-ft; wind in mph, kt or m/s; °F/°C; ft/m; in/mm.
- **Colorblind-safe colors** and **Lite mode** (less data; on automatically on phones).
- On phones, the controls open from **☰ Controls**.

### Overnight smoke settling risk (this site's screening index)
For each night (6 PM–8 AM local), from the NBM forecast:
- **Moderate:** lowest mixing height under 1,500 ft *and* lowest transport wind under 10 kt
- **High:** under 1,000 ft *and* under 6 kt
- **Very high:** High plus RH ≥ 85% (smoke + fog)

This is not LVORI. Lavdas ADI and LVORI can be added once their published lookup tables are
provided (they are only published as chart images).

## How it works
A GitHub Action (`.github/workflows/update-data.yml`) runs every 3 hours. It downloads only the needed
fields from NOAA's public NBM and HRRR data on AWS, plus NWS alerts, AirNow, EPA and NIFC layers. It finds
each map cell's time zone, packages everything with `index.html`, and publishes to GitHub Pages.
Nothing is committed to the repository.

| File | What it does |
|---|---|
| `index.html` | The web page |
| `fetch_data.py` | Downloads and processes all data into `data/` |
| `.github/workflows/update-data.yml` | Runs the script every 3 hours and publishes the site |

---

## Step-by-step setup (about 15 minutes, all in a web browser)

1. **Sign in** at https://github.com as **SWCCFireWx**.
2. **Create a repository:** **+ → New repository**, name `us-smoke-planner`, **Public**,
   check **Add a README file** → **Create repository**.
3. **Upload files:** **Add file → Upload files**, drag in `index.html`, `fetch_data.py`, `README.md`
   (the files, not the zip) → **Commit changes**.
4. **Workflow file:** **Add file → Create new file**, name it exactly `.github/workflows/update-data.yml`,
   paste in the contents of `update-data.yml` → **Commit changes**.
5. **Turn on Pages:** **Settings → Pages → Source: GitHub Actions**.
6. **First run:** **Actions → Update data → Run workflow**. Wait for **build** and **deploy** to show
   green checks (about 25–45 minutes the first time).
7. **Open** https://swccfirewx.github.io/us-smoke-planner/ and share the link.

**Updating an existing copy:** upload the new `index.html`, `fetch_data.py` and `README.md` over the old
ones, then open `.github/workflows/update-data.yml`, click the pencil, replace everything with the new
version, and commit. The commit starts a new run automatically.

## Checking the first real run
Open the run log (**Actions → the run → build → "Build site and data"**) and look for:
- `NBM spread fields found: …`: which spread fields NBM provided. If `none`, the two confidence layers
  will say spread isn't available.
- `features from PSA_ERC_and_BI_Percentiles_and_Trends …; fields: …`: the ERC/BI field names found.
  If the page doesn't color PSAs by ERC, send these lines so the field match can be adjusted.
- `HRRR smoke from …`: the HRRR run used.
- Any line starting with `!!`: that source failed this run. Every source except NBM CONUS is optional,
  and the site still publishes without it.

## Limitations
- NBM and HRRR are guidance. They under-resolve valley inversions and drainage flow; overnight smoke in
  valleys can be worse than shown. **Always get a NWS spot forecast** and follow the burn plan and smoke permit.
- Confidence uses NBM's ensemble standard deviation with a normal approximation.
- HRRR smoke covers CONUS only, about 48 hours.
- AirNow values are the latest hour, not NowCast or 24-hour AQI.
- Default ventilation categories are the common NWS scale; set your state's breakpoints on the page.
- GitHub Pages has a soft limit of about 100 GB of downloads a month (thousands of visits). For heavy
  agency-wide use, host the same files on an agency web server.

## Troubleshooting
| Problem | Fix |
|---|---|
| Site shows 404 | Settings → Pages → Source must be **GitHub Actions**, then run the workflow |
| "No forecast data yet" | Run the workflow and wait for both green checks, then reload |
| Red ✗ on a run | Open the failed step and read the last lines; NOAA is sometimes briefly unavailable. The previous site stays up |
| Header says "not in this run: Alaska" (etc.) | That grid failed this run (log line `!! ak failed`); usually fixes itself |
| Legend says alerts are "as of the last data update" | The live NWS check didn't respond; the cached copy (up to 3 h old) is shown |
| Header shows "⚠ … h old" | Check the Actions tab; rerun **Update data** |
| Run cancelled at 120 minutes | NOAA downloads slow; rerun. If frequent, raise `timeout-minutes` |
| Scheduled runs stopped | GitHub pauses schedules after 60 days with no repo activity: **Actions → Update data → Enable workflow** |

## Data sources
- NOAA National Blend of Models (core; CONUS, Alaska, Hawaii, Puerto Rico), NOAA Open Data on AWS `noaa-nbm-grib2-pds`
- NOAA HRRR smoke (MASSDEN 8 m, COLMD), NOAA Open Data on AWS `noaa-hrrr-bdp-pds`
- NWS watches/warnings map service (Red Flag Warnings, Fire Weather Watches)
- AirNow public hourly data files (files.airnowtech.org)
- EPA Mandatory Class I Federal Areas; EPA Nonattainment Areas (PM2.5 24-hr 2006, annual 2012)
- NIFC: WFIGS current perimeters, GACC boundaries, PSA boundaries, ERC/BI percentiles by PSA
- Time zones: timezone-boundary-builder via the `timezonefinder` Python package
