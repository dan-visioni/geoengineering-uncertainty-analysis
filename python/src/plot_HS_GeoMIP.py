# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#   kernelspec:
#     display_name: Python (geomip)
#     language: python
#     name: geomip
# ---

# %% [markdown]
# # plot_HS_GeoMIP — regional GeoMIP surface-temperature response
#
# Loads the 6-model GeoMIP `tas` ensemble and computes cosine-latitude-weighted
# global and IPCC-WGI-region mean temperature responses for G6sulfur, G6solar,
# and SSP5-8.5 (each relative to SSP2-4.5, 2020–2100). Plots per-region time
# series with per-model thin lines and the multi-model mean. Also computes
# end-of-century (2080–2100) multi-model variances saved to `data_py/tas_regions1_std.npz`.

# %%
import sys
import os as _os
_r = _os.getcwd()
while not _os.path.isdir(_os.path.join(_r, "python")) and _os.path.dirname(_r) != _r:
    _r = _os.path.dirname(_r)
_os.chdir(_r)   # notebooks live in j_notebooks/; anchor to project root
sys.path.insert(0, "python")
import numpy as np
import matplotlib.pyplot as plt
import mat_utils as mu
from mat_utils import brewermap, movmean, shaperead, inpolygon

# %% [markdown]
# ## Load data
# `GeoMIP_tas` contains the multi-model ensemble: `x_e[model, experiment, member]`
# of `tas` fields (lon × lat × year) plus the per-model `LAT`/`LON` coordinate arrays.

# %%
m = mu.load_npz("GeoMIP_tas")     # data_py/GeoMIP_tas.npz (regenerated from raw; see read_GeoMIP_models)
x_e = m["x_e"]          # object array, shape (6, 4, 3)
LAT = m["LAT"]          # object array, shape (1, 6)
LON = m["LON"]          # object array, shape (1, 6)

# %% [markdown]
# ## Global means
# For each model/experiment, take the zonal mean (over longitude) and the
# cosine-latitude weighted mean (over latitude). Models whose record starts in
# 2015 (86 years) are trimmed to 2020–2100.

# %%
x_e_gm = np.full((6, 4, 3, 81), np.nan)

for im in range(6):
    lat = LAT[0, im].ravel()
    ww = np.cos(lat / 180 * np.pi)
    for ie in range(4):
        for inn in range(1):                       # MATLAB: for in=1:1
            field = x_e[im, ie, inn]
            if np.size(field) > 0 and field.shape[0] > 0:
                zm = np.mean(field, axis=0)         # mean over lon -> (lat, time)
                y = np.sum(zm * ww[:, None], axis=0) / np.sum(ww)
                if y.shape[0] == 86:
                    x_e_gm[im, ie, inn, 0:81] = y[5:86]
                else:
                    x_e_gm[im, ie, inn, :] = y

# %% [markdown]
# ## Regional means
# Regions come from the IPCC-WGI reference regions shapefile (the
# WGII file is read first then immediately overwritten by the WGI file, so only
# the WGI regions matter). For each model the area-weighted regional mean of
# the (experiment − SSP2-4.5) difference is accumulated for each year.

# %%
regions = shaperead("IPCC-WGI-reference-regions-v4.shx")
nR = len(regions)

ta_reg_sulf = np.full((6, 81, nR), np.nan)
ta_reg_solar = np.full((6, 81, nR), np.nan)
ta_reg_ssp = np.full((6, 81, nR), np.nan)

for im in range(6):
    lat = LAT[0, im].ravel()
    lon = LON[0, im].ravel()
    ww = np.cos(lat / 180 * np.pi)

    # MATLAB meshgrid(LAT,LON-180) and meshgrid(ww,LON) -> arrays (nlon, nlat)
    latM, lonM = np.meshgrid(lat, lon - 180)       # both (nlon, nlat)
    AA = np.meshgrid(ww, lon)[0]                    # (nlon, nlat), AA[i,j]=ww[j]

    ta_1 = (x_e[im, 0, 0][:, :, 0:81] - x_e[im, 3, 0][:, :, 0:81]) * AA[:, :, None]
    ta_2 = (x_e[im, 1, 0][:, :, 0:81] - x_e[im, 3, 0][:, :, 0:81]) * AA[:, :, None]
    ta_3 = (x_e[im, 2, 0][:, :, 0:81] - x_e[im, 3, 0][:, :, 0:81]) * AA[:, :, None]

    for i in range(nR):
        print(regions[i].Name)
        ind = inpolygon(latM, lonM, regions[i].Y, regions[i].X)
        wS = np.sum(AA[ind])
        for j in range(81):
            ta_reg_sulf[im, j, i] = np.nansum(ta_1[:, :, j][ind]) / wS
            ta_reg_solar[im, j, i] = np.nansum(ta_2[:, :, j][ind]) / wS
            ta_reg_ssp[im, j, i] = np.nansum(ta_3[:, :, j][ind]) / wS

# %% [markdown]
# ## Plot the first region
# Per-model thin lines plus the multi-model mean (thick), each smoothed with a
# 10-year moving mean. Black = G6sulfur, blue = G6solar, red = SSP5-8.5.

# %%
i = 0                       # MATLAB region index 1
time = np.arange(2020, 2101)
mvm = 10

fig = plt.figure()
ax = fig.add_subplot(2, 1, 1)
for im in range(6):
    ax.plot(time, movmean(ta_reg_sulf[im, :, i], mvm), "k")
    ax.plot(time, movmean(ta_reg_solar[im, :, i], mvm), "b")
    ax.plot(time, movmean(ta_reg_ssp[im, :, i], mvm), "r")

ax.plot(time, movmean(np.mean(ta_reg_sulf[:, :, i], axis=0), mvm), "k", linewidth=2)
ax.plot(time, movmean(np.mean(ta_reg_solar[:, :, i], axis=0), mvm), "b", linewidth=2)
ax.plot(time, movmean(np.mean(ta_reg_ssp[:, :, i], axis=0), mvm), "r", linewidth=2)

ax.set_title(regions[i].Name, fontsize=18, fontweight="bold", color="k")
ax.set_ylabel("K")
for s in ax.spines.values():
    s.set_linewidth(2)
ax.tick_params(labelsize=14)
fig.set_facecolor("w")
ax.set_facecolor("w")
plt.show()

# %% [markdown]
# ## End-of-century variances
# Multi-model variance of the 2080–2100 mean response, using the sample
# standard deviation (ddof=1).

# %%
ta_var_sulf = np.full(nR, np.nan)
ta_var_solar = np.full(nR, np.nan)
ta_var_ssp = np.full(nR, np.nan)

for i in range(nR):
    x = np.mean(ta_reg_sulf[:, 60:81, i], axis=1)
    y = np.mean(ta_reg_solar[:, 60:81, i], axis=1)
    z = np.mean(ta_reg_ssp[:, 60:81, i], axis=1)
    ta_var_sulf[i] = np.std(x, ddof=1) ** 2
    ta_var_solar[i] = np.std(y, ddof=1) ** 2
    ta_var_ssp[i] = np.std(z, ddof=1) ** 2

print("region 0 variances (sulf, solar, ssp):",
      ta_var_sulf[0], ta_var_solar[0], ta_var_ssp[0])

# %% [markdown]
# ## Save outputs
# The three end-of-century variance arrays (`ta_var_sulf`, `ta_var_solar`,
# `ta_var_ssp`) are saved to `data_py/tas_regions1_std.npz` for use by the
# downstream scatter analyses (`plot_g6vsSSPunc`,
# `plot_g6_uncertainties_aerosols_vs_solar`).

# %%
mu.save_npz("tas_regions1_std",
            ta_var_sulf=ta_var_sulf,
            ta_var_solar=ta_var_solar,
            ta_var_ssp=ta_var_ssp)
print("saved data_py/tas_regions1_std.npz")
