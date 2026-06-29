# Uncertainties — Python / Jupyter port

This folder's MATLAB analyses have been ported to Python Jupyter notebooks. The
translation is **faithful**: every computation reproduces the MATLAB analysis;
only the plotting back-end (matplotlib / Cartopy) and the I/O format differ.

## Environment

A dedicated conda env `geomip` (conda-forge) holds everything:

```
conda create -n geomip -c conda-forge python=3.11 numpy scipy matplotlib pandas \
    xarray netcdf4 cartopy pyshp palettable jupyter nbconvert ipykernel jupytext
python -m ipykernel install --user --name geomip --display-name "Python (geomip)"
```

The notebooks live in **`j_notebooks/`**. Open them with the **Python (geomip)**
kernel. Their first cell automatically re-anchors the working directory to the
project root (it walks up until it finds `python/`), so all the relative data
paths resolve regardless of where the notebook is launched from.

## Layout

| Path | Contents |
|---|---|
| `j_notebooks/*.ipynb` | the executed notebooks, one per analysis |
| `python/src/*.py` | the notebook sources (jupytext "percent" format — edit these) |
| `python/mat_utils.py` | MATLAB-compatibility helpers (see below) |
| `python/build.sh` | build one notebook: `bash python/build.sh <name>` |
| `python/run_all.py` | rebuild + execute **all** notebooks in dependency order, with a pass/fail table (`--clean` wipes generated data first) |
| `data_py/*.npz` | datasets **generated** by notebooks in this folder (Python-native) |
| `py_outputs/*.png` | regenerated figures (high-resolution PNG, 300 dpi) |

The original `.m`, `.mat`, `.eps` files are left **untouched**.

## `mat_utils` helpers

Faithful ports of the MATLAB-specific behaviour used by the scripts:

- `brewermap(N, scheme)` — Stephen Cobeldick's ColorBrewer function (exact node
  tables, Lab-space cubic-spline interpolation, `*`-reversal, qualitative repeat).
- `brighten(rgb, beta)` — MATLAB `brighten` (`gamma = (1-beta)/(1+beta)`).
- `movmean(x, k)` — MATLAB moving mean (centred window, even windows lean left,
  shrinking edge windows, NaN propagation).
- `ncread(file, var)` — MATLAB dimension order (full axis reversal), fill → NaN.
- `shaperead` / `inpolygon` — MATLAB shapefile records (NaN-separated `.X`/`.Y`)
  and point-in-polygon (multipart polygons, boundary included).
- `interp1` / `interp2` — MATLAB-compatible linear interpolation with an explicit
  extrapolation value.
- `worldmap` / `pcolorm` / `contourfm` / `geoshow_coastlines` — Cartopy stand-ins
  for the MATLAB Mapping Toolbox calls.
- `save_npz(name, …)` / `load_npz(name)` — read/write `data_py/<name>.npz`.

## Data format policy (`.mat` → `.npz`)

- **External inputs** (no generator in this folder — e.g. `../GeoMIP_G6Sulfur/*.mat`,
  `borderdata.mat`, the provided `*_regions{1,2,3}.mat` / `Ndeg_*.mat`
  intermediates) are read as `.mat` (unchanged).
- **Files created by a script in this folder** are written and read in the
  Python-native `.npz` format under `data_py/` (the original `.mat` is never
  overwritten). These are: `GeoMIP_tas`, `GeoMIP_pr`, `data_all_SOCOL`,
  `map_plots_GeoMIP_SOCOL`, `aod_SOCOL`, `tas_timeseries_SOCOL`, `G615kSAI`,
  `tas_regions1_std`, `pr_regions1_std`, `2deg_pr`, `AODVIS`.

The data-prep notebooks (`read_GeoMIP_models`, `save_SOCOL`,
`save_SOCOL_alternatw`) regenerate their datasets from the raw NetCDF on
`/Volumes/D_Visioni2/…` and **verify** them against the existing `.mat`
references (the GeoMIP `tas`/`pr` datasets reproduce bit-identically; the SOCOL
ones match to float32 round-off).

## Notebook ↔ MATLAB mapping

| Notebook | MATLAB source(s) | Notes |
|---|---|---|
| `read_GeoMIP_models` | `read_GeoMIP_models.m` | runs for both `tas` and `pr` → `data_py/GeoMIP_{tas,pr}.npz` |
| `save_SOCOL` | `save_SOCOL.m` | → `tas_timeseries_SOCOL`, `map_plots_GeoMIP_SOCOL`, `aod_SOCOL` |
| `save_SOCOL_alternatw` | `save_SOCOL_alternatw.m` | → `data_all_SOCOL` |
| `include_SOCOL_HS` | `include_SOCOL_HS.m` | GeoMIP + SOCOL global/regional uncertainty; global time series, tas/pr partition, and the `..._tas_Africa` / `..._pr_Africa` 7-region panels |
| `plot_HS_GeoMIP` | `plot_HS_GeoMIP.m` | also writes the 3-var `tas_regions1_std` |
| `plot_HSpr_GeoMIP` | `plot_HSpr_GeoMIP.m` | writes `pr_regions1_std` |
| `plot_HS_countries` | `plot_HS_countries.m` | writes `tas_regions1_std_countries` (see name-collision note) |
| `plot_tas_vs_pr` | `plot_tas_vs_pr.m` | |
| `plot_g6vsSSPunc` | `plot_g6vsSSPunc.m` | |
| `plot_g6_uncertainties_aerosols_vs_solar` | `plot_g6_uncertainties_aerosols_vs_solar.m` | |
| `plot_spatial_aggregation_tas` | `plot_spatial_aggregation_tas.m` | writes `2deg_pr` |
| `plot_all_spatial_aggregation` | `plot_all_spatial_aggregation.m` | |
| `plot_maps` | `plot_maps.m` | Cartopy |
| `plot_maps_15k` | `plot_maps_15k.m` | Cartopy; writes `G615kSAI` |
| `compare_g6s` | `compare_g6s.m` | Cartopy |
| `plot_uncertainty` | `plot_uncertainty.m` | |
| `plot_AOD_vs_tas` | `plot_AOD_vs_tas.m` | |
| `plot_HS_TREFHT` | `plot_HS.m` + `plot_regions_T.m` | merged (shared MATLAB workspace) |
| `plot_HS_PRECT` | `plot_HS_p.m` + `plot_regions_P.m` | merged |
| `plot_AOD` | `plot_AOD.m` + `plot_AOD_regions.m` | merged; writes `AODVIS` |
| `plot_HS_regions_map` | (new) | map "dashboard": 8 WGII regions coloured on a world map with per-region line + partition insets and leader lines, for `tas` and `pr` |

`brewermap.m` became `mat_utils.brewermap` (a helper, not a notebook). The three
`plot_regions_*`/`plot_AOD_regions` scripts reused another script's workspace, so
they were merged into the corresponding loader notebook.

## Faithfulness notes / documented deviations

- **`inpolygon` boundary handling** differs slightly between MATLAB and
  matplotlib's `Path.contains_points`, so region-edge grid cells can be assigned
  marginally differently. This produces relative differences of ~0.1–0.3 % in the
  regional variance values vs the shipped `*_std.mat` — not a translation bug.
- **SOCOL fields** in the reference `.mat` were stored as float32; the
  regenerated values match to float32 round-off (≤ ~2×10⁻⁴ on ~250 K fields).
- **`tas_regions1_std` name collision:** in MATLAB both `plot_HS_GeoMIP.m`
  (3 variables, via a commented `save`) and `plot_HS_countries.m` (2 variables)
  write `tas_regions1_std.mat`. The downstream scatter scripts need the 3-variable
  version, which matches the shipped file, so `plot_HS_GeoMIP` produces
  `tas_regions1_std` and `plot_HS_countries` writes `tas_regions1_std_countries`.
- **`plot_AOD.m`** references undefined baseline variables (`tab_zm`/`tab_gm`) — a
  latent MATLAB bug that does not affect any AOD output; those lines are omitted.
- **`plot_AOD_vs_tas.m`** ends with a non-functional `figure(3)` loop (indexes an
  81-row colormap up to 91 and scatters mismatched lengths); the real outputs are
  reproduced and that dead fragment is omitted.
- **`plot_single`** map blocks (false in the MATLAB) and unused `T0/T1/T2`
  zonal-gradient terms are omitted where they produce no output.
- Cartopy maps use the **Robinson** projection (matching MATLAB `worldmap`'s
  default oval) and reproduce the same panels, colormaps, data ranges and
  sign-agreement stippling; coastlines are black (vs blue) and the manual MATLAB
  subplot/colorbar position nudges are not pixel-matched.
- A parity pass against the original MATLAB `.pdf`/`.eps` outputs confirmed every
  figure matches. Fixes applied during that pass: corrected the `include_SOCOL_HS`
  partition legend order; set the `plot_HS_*`/`plot_AOD` variance-partition areas
  to MATLAB's default blue/orange/yellow; fixed the `plot_HS_PRECT` baseline
  (zero- vs NaN-fill) so the ensemble-mean lines render; consolidated
  `plot_all_spatial_aggregation` into a single cell (split cells + the inline
  backend had produced an empty figure); switched maps to Robinson.

## Reproducing everything

```
cd "<this folder>"
/opt/anaconda3/envs/geomip/bin/python python/run_all.py --clean
```

This wipes `data_py/` and `py_outputs/`, regenerates the in-folder datasets from
the raw NetCDF, executes every notebook in dependency order, and prints a
pass/fail table. Requires `/Volumes/D_Visioni2/…` (raw model output) and the
parent `../GeoMIP_G6Sulfur/` folder to be present.
