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
# # plot_AOD_vs_tas — zonal-mean stratospheric AOD across GeoMIP models
#
# Plots the end-of-century (2090–2100) zonal-mean stratospheric AOD for each
# GeoMIP model under G6sulfur and SSP5-8.5, interpolates them onto a common 2°
# latitude grid, and computes the inter-model standard deviation as a measure
# of the structural spread across models.

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
from mat_utils import brewermap, interp1

a = mu.loadmat("../GeoMIP_G6Sulfur/AOD_GeoMIP.mat")
lat = np.asarray(a["lat"], float).ravel()
lat2 = np.asarray(a["lat2"], float).ravel()
lat3 = np.asarray(a["lat3"], float).ravel()
x2_y_CNRM = np.asarray(a["x2_y_CNRM"], float).ravel()
x2_y_MPI = np.asarray(a["x2_y_MPI"], float).ravel()
x2_y_MPIH = np.asarray(a["x2_y_MPIH"], float).ravel()
LAT = a["LAT"]; x_lat = a["x_lat"]

# per-model zonal AOD curves (G6sulfur - SSP, as in the MATLAB)
def lat_of(i):      # LAT{i} (1-based) -> 0-based
    return np.asarray(LAT[0, i - 1], float).ravel()
def xl(i, j):       # x_lat{i,j} (1-based)
    return np.asarray(x_lat[i - 1, j - 1], float).ravel()

# %% [markdown]
# ## Zonal-mean stratospheric AOD per model

# %%
fc1 = brewermap(7, "Set2")
fig, ax = plt.subplots(figsize=(6, 7))
ax.plot(lat, x2_y_CNRM, color=fc1[0], linewidth=2)
ax.plot(lat_of(2), xl(2, 1) - xl(2, 3) * 0, color=fc1[1], linewidth=2)
ax.plot(lat_of(3), xl(3, 1) - xl(3, 3), color=fc1[2], linewidth=2)
ax.plot(lat_of(4), xl(4, 1) - xl(4, 3), color=fc1[3], linewidth=2)
ax.plot(lat2, x2_y_MPI, color=fc1[4], linewidth=2)
ax.plot(lat3, x2_y_MPIH, color=fc1[5], linewidth=2)
ax.set_title("a) 2090-2100 zonal mean stratospheric AOD")
ax.set_xticks(np.arange(-90, 91, 30)); ax.set_xlim(-90, 90); ax.set_ylim(0, .9)
for s in ax.spines.values():
    s.set_linewidth(2)
plt.show()

# %% [markdown]
# ## Interpolate onto a common 2° grid and take the inter-model spread

# %%
lat_int1 = np.arange(-90, 91, 2.0)
aod_all = np.zeros((lat_int1.size, 6))
aod_all[:, 0] = interp1(lat, x2_y_CNRM, lat_int1, "linear", np.nan)
aod_all[:, 1] = interp1(lat_of(2), xl(2, 1) - xl(2, 3) * 0, lat_int1, "linear", np.nan)
aod_all[:, 2] = interp1(lat_of(3), xl(3, 1) - xl(3, 3), lat_int1, "linear", np.nan)
aod_all[:, 3] = interp1(lat_of(4), xl(4, 1) - xl(4, 3), lat_int1, "linear", np.nan)
aod_all[:, 4] = interp1(lat2, x2_y_MPI, lat_int1, "linear", np.nan)
aod_all[:, 5] = interp1(lat3, x2_y_MPIH, lat_int1, "linear", np.nan)

aod_all_std = np.std(aod_all, axis=1, ddof=1)
print("aod_all shape:", aod_all.shape, "| max inter-model std:", np.nanmax(aod_all_std))

# 2deg_tas.mat (provided input) is loaded by the MATLAB script for the
# subsequent (non-functional) figure(3); we load it for completeness.
td = mu.loadmat("2deg_tas.mat")
print("2deg_tas sulf_delta_std1 shape:", np.asarray(td["sulf_delta_std1"]).shape)
