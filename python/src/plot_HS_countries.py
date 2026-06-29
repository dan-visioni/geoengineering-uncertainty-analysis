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
# # plot_HS_countries — regional GeoMIP surface-temperature response
#
# Loads the 6-model GeoMIP `tas` ensemble and computes NaN-omitting
# cosine-latitude-weighted global and IPCC-WGI-region mean temperature responses
# for G6sulfur and G6solar (relative to SSP2-4.5, 2020–2100). Computes the
# end-of-century (2080–2100) multi-model variances for both experiments, saved
# to `data_py/tas_regions1_std_countries.npz` for use in downstream scatter analyses.

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
from mat_utils import brewermap, brighten, movmean, ncread, shaperead, inpolygon, interp1, interp2

# %% [markdown]
# ## Load data
# `GeoMIP_tas` contains `x_e[model, experiment, member]` of `tas` fields
# (lon × lat × year) plus the per-model `LAT`/`LON` coordinate arrays.
# `borderdata.mat` is loaded for consistency but is not used in any computation.

# %%
m = mu.load_npz("GeoMIP_tas")     # data_py/GeoMIP_tas.npz (regenerated from raw; see read_GeoMIP_models)
x_e = m["x_e"]          # object array, shape (6, 4, 3)
LAT = m["LAT"]          # object array, shape (1, 6)
LON = m["LON"]          # object array, shape (1, 6)

borderdata = mu.loadmat("borderdata.mat")   # loaded for fidelity; unused

# %% [markdown]
# ## Global means
# For each model/experiment/member, take the (NaN-omitting) zonal mean over
# longitude and the cosine-latitude weighted mean over latitude. Models whose
# record starts in 2015 (86 years) are trimmed to 2020–2100. All three
# ensemble members are used.

# %%
x_e_gm = np.full((6, 4, 3, 81), np.nan)

for im in range(6):
    lat = LAT[0, im].ravel()
    ww = np.cos(lat / 180 * np.pi)
    for ie in range(4):
        for inn in range(3):                        # MATLAB: for in=1:3
            field = x_e[im, ie, inn]
            if np.size(field) > 0 and field.shape[0] > 0:
                zm = np.nanmean(field, axis=0)      # mean over lon (omitnan)
                y = np.sum(zm * ww[:, None], axis=0) / np.sum(ww)
                if y.shape[0] == 86:
                    x_e_gm[im, ie, inn, 0:81] = y[5:86]
                else:
                    x_e_gm[im, ie, inn, :] = y

# %% [markdown]
# ## Regional means
# Regions come from the IPCC-WGI reference regions shapefile. For each model the
# area-weighted regional mean of the (experiment − SSP2-4.5) difference is
# accumulated for every year and **every** region. Only the two SRM experiments
# are carried (G6sulfur, G6solar).

# %%
regions = shaperead("IPCC-WGI-reference-regions-v4.shx")
nR = len(regions)

ta_reg_sulf = np.full((6, 81, nR), np.nan)
ta_reg_solar = np.full((6, 81, nR), np.nan)

for im in range(6):
    lat = LAT[0, im].ravel()
    lon = LON[0, im].ravel()
    ww = np.cos(lat / 180 * np.pi)

    # MATLAB meshgrid(LAT,LON-180) and meshgrid(ww,LON) -> arrays (nlon, nlat)
    latM, lonM = np.meshgrid(lat, lon - 180)       # both (nlon, nlat)
    AA = np.meshgrid(ww, lon)[0]                    # (nlon, nlat), AA[i,j]=ww[j]

    ta_1 = (x_e[im, 0, 0][:, :, 0:81] - x_e[im, 3, 0][:, :, 0:81]) * AA[:, :, None]
    ta_2 = (x_e[im, 1, 0][:, :, 0:81] - x_e[im, 3, 0][:, :, 0:81]) * AA[:, :, None]

    for i in range(nR):
        print(regions[i].Name)
        ind = inpolygon(latM, lonM, regions[i].Y, regions[i].X)
        wS = np.sum(AA[ind])
        for j in range(81):
            ta_reg_sulf[im, j, i] = np.nansum(ta_1[:, :, j][ind]) / wS
            ta_reg_solar[im, j, i] = np.nansum(ta_2[:, :, j][ind]) / wS

# %% [markdown]
# ## Plot the first region
# Per-region time-series plot for the first IPCC-WGI region: per-model thin
# lines plus the multi-model mean (thick), each smoothed with a 10-year moving
# mean. Black = G6sulfur, blue = G6solar.

# %%
i = 0                       # MATLAB region index 1
time = np.arange(2020, 2101)
mvm = 10

fig = plt.figure()
ax = fig.add_subplot(2, 1, 1)
for im in range(6):
    ax.plot(time, movmean(ta_reg_sulf[im, :, i], mvm), "k")
    ax.plot(time, movmean(ta_reg_solar[im, :, i], mvm), "b")

ax.plot(time, movmean(np.nanmean(ta_reg_sulf[:, :, i], axis=0), mvm), "k", linewidth=2)
ax.plot(time, movmean(np.nanmean(ta_reg_solar[:, :, i], axis=0), mvm), "b", linewidth=2)

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
# Multi-model variance of the 2080–2100 mean response, computed NaN-omitting
# with the sample standard deviation (ddof=1).

# %%
ta_var_sulf = np.full(nR, np.nan)
ta_var_solar = np.full(nR, np.nan)

for i in range(nR):
    x = np.nanmean(ta_reg_sulf[:, 60:81, i], axis=1)
    y = np.nanmean(ta_reg_solar[:, 60:81, i], axis=1)
    ta_var_sulf[i] = np.std(x, ddof=1) ** 2
    ta_var_solar[i] = np.std(y, ddof=1) ** 2

print("region 0 variances (sulf, solar):",
      ta_var_sulf[0], ta_var_solar[0])

# %% [markdown]
# ## Save outputs
# The two end-of-century variance arrays (`ta_var_sulf`, `ta_var_solar`) are
# saved to `data_py/tas_regions1_std_countries.npz`.
# Note: `plot_HS_GeoMIP` saves a 3-variable file called `tas_regions1_std`
# (which also includes `ta_var_ssp`) that the downstream scatter scripts depend
# on. This script's distinct 2-variable output is saved under a separate name
# to avoid any conflict.

# %%
mu.save_npz("tas_regions1_std_countries",
            ta_var_sulf=ta_var_sulf,
            ta_var_solar=ta_var_solar)
print("saved data_py/tas_regions1_std_countries.npz")
