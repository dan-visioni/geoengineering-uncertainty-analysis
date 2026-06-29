# Geoengineering Uncertainty Analysis

Code repository accompanying the paper on uncertainty decomposition in SAI (Stratospheric Aerosol Injection) geoengineering scenarios, comparing GeoMIP G6sulfur / G6solar multi-model ensembles with the CESM2 MA-GAUSS and SOCOL injection-strategy ensembles.

---

## Repository structure

```
python/src/       Jupytext percent-format notebook sources (.py)
j_notebooks/      Built Jupyter notebooks (.ipynb, with outputs)
data_py/          Small pre-computed intermediates (committed here)
                  Large files: download from Zenodo (see below)
py_outputs/       300-dpi PNG figures (and companion PDFs)
matlab/           Legacy MATLAB scripts (reference only)
IPCC-WGI-*        IPCC AR6 WGI reference regions shapefile (58 regions)
IPCC_WGII_*       IPCC AR6 WGII continental regions shapefile (8 regions)
*.mat             Static spatial-resolution inputs read by notebooks
```

---

## Setup

### 1 — Create the conda environment

```bash
conda env create -f environment.yml
conda activate geomip
```

### 2 — Download large data files from Zenodo

Three large intermediate files are too big for git and must be downloaded separately. Place them in `data_py/`:

| File | Size | Zenodo DOI |
|------|------|-----------|
| `GeoMIP_tas.npz` | 1.3 GB | *(to be added after upload)* |
| `GeoMIP_pr.npz` | 1.3 GB | *(to be added after upload)* |
| `data_all_SOCOL.npz` | 354 MB | *(to be added after upload)* |

For `plot_HS_regions_map`, the CESM2 MA-GAUSS raw NetCDF data is also needed (21 GB). Download from Zenodo and place at `/Volumes/D_Visioni2/CESM2_MA/`, or update the path variable at the top of `python/src/plot_HS_TREFHT.py` and `python/src/plot_HS_PRECT.py`.

### 3 — Build notebooks from source (optional)

The `j_notebooks/` directory already contains built `.ipynb` files with outputs. To rebuild from the `.py` sources:

```bash
cd python && bash build.sh
```

---

## Notebook dependency order

If you need to regenerate the large intermediates from raw NetCDF data:

| Step | Notebook | Output |
|------|----------|--------|
| 1 | `read_GeoMIP_models` | `data_py/GeoMIP_tas.npz`, `data_py/GeoMIP_pr.npz` |
| 2 | `save_SOCOL_alternatw` | `data_py/data_all_SOCOL.npz` |
| 3 | `plot_HS_GeoMIP` | `data_py/tas_regions1_std.npz` |
| 4 | `plot_HSpr_GeoMIP` | `data_py/pr_regions1_std.npz` |
| 5 | `plot_spatial_aggregation_tas` | `data_py/2deg_pr.npz` |

Once `data_py/` is populated from Zenodo, you can run the four manuscript-figure notebooks directly without step 1–5.

---

## Manuscript figures

| Figure file | Producing notebook |
|---|---|
| `GeoMIP_uncertainties_aerosols.png` | `plot_g6_uncertainties_aerosols_vs_solar` |
| `GeoMIP_spatial_uncertainties.png` | `plot_all_spatial_aggregation` |
| `Africa_SOCOL.png` | `include_SOCOL_HS` |
| `GeoMIP_uncertainties_aerosols_SOCOL_tas_pr.png` | `include_SOCOL_HS` |
| `HS_plot_TREFHT_regions_map.png` | `plot_HS_regions_map` |
| `HS_plot_PRECT_regions_map.png` | `plot_HS_regions_map` |

---

## Raw data provenance

| Dataset | Used by | How to obtain |
|---------|---------|---------------|
| GeoMIP CMIP6 (G6sulfur, G6solar, SSP5-8.5, SSP2-4.5), 6 models | `read_GeoMIP_models` | [ESGF node](https://esgf-node.llnl.gov/search/cmip6/) |
| SOCOL G6sulfur ensemble (3 members) | `save_SOCOL` | Zenodo *(DOI to be added)* |
| SOCOL alternative-strategy ensemble (G6, al, ca, di × 3 members) | `save_SOCOL_alternatw` | Zenodo *(DOI to be added)* |
| CESM2 MA-GAUSS 7-strategy ensemble (7 × 3 members) | `plot_HS_TREFHT`, `plot_HS_PRECT`, `plot_HS_regions_map` | Zenodo *(DOI to be added)* |
