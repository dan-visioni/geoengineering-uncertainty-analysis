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
# # plot_tas_vs_pr — fractional aerosol variance, tas vs pr
#
# Loads precomputed regional response ensembles for surface temperature (`tas`)
# and precipitation (`pr`) at three spatial resolutions:
#
# * **regions1** — IPCC-WGI reference regions (58 regions),
# * **regions2** — IPCC-WGII continental regions (8 regions),
# * **regions3** — a fine 238-region grid.
#
# For each region, computes the multi-model sample variance of the G6sulfur
# and G6solar responses, sums them to a total variance, and scatters the
# time-mean (2060–2100) fractional sulfur variance for tas against pr.
# Points are coloured by IPCC-WGII continental group.
#
# The analysis
# (matplotlib) differs.

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
# ## regions1 — IPCC-WGI reference regions (58)
# Per-region, per-timestep multi-model sample variance of the sulfur and solar
# responses, then their sum.

# %%
m = mu.loadmat("tas_regions1.mat")
ta_reg_sulf = m["ta_reg_sulf"]          # (6, 81, 58)
ta_reg_solar = m["ta_reg_solar"]

regions = shaperead("IPCC-WGI-reference-regions-v4.shx")
nR1 = len(regions)                      # 58

tas_var1_sulf = np.full((nR1, 81), np.nan)
tas_var1_solar = np.full((nR1, 81), np.nan)
tot_var_reg1_tas = np.full((nR1, 81), np.nan)

for i in range(nR1):
    for t in range(81):
        x = ta_reg_sulf[:, t, i]
        y = ta_reg_solar[:, t, i]
        tas_var1_sulf[i, t] = np.std(x, ddof=1) ** 2
        tas_var1_solar[i, t] = np.std(y, ddof=1) ** 2
    tot_var_reg1_tas[i, :] = tas_var1_sulf[i, :] + tas_var1_solar[i, :]

# %% [markdown]
# ## regions2 — IPCC-WGII continental regions (8)

# %%
regions = shaperead("IPCC_WGII_continental_regions.shx")

m = mu.loadmat("tas_regions2.mat")
ta_reg_sulf = m["ta_reg_sulf"]          # (6, 81, 8)
ta_reg_solar = m["ta_reg_solar"]

nR2 = len(regions)                      # 8

tas_var2_sulf = np.full((nR2, 81), np.nan)
tas_var2_solar = np.full((nR2, 81), np.nan)
tot_var_reg2_tas = np.full((nR2, 81), np.nan)

for i in range(nR2):
    for t in range(81):
        x = ta_reg_sulf[:, t, i]
        y = ta_reg_solar[:, t, i]
        tas_var2_sulf[i, t] = np.std(x, ddof=1) ** 2
        tas_var2_solar[i, t] = np.std(y, ddof=1) ** 2
    tot_var_reg2_tas[i, :] = tas_var2_sulf[i, :] + tas_var2_solar[i, :]

# %% [markdown]
# ## regions3 — fine 238-region grid

# %%
m = mu.loadmat("tas_regions3.mat")
ta_reg_sulf = m["ta_reg_sulf"]          # (6, 81, 238)
ta_reg_solar = m["ta_reg_solar"]

tas_var3_sulf = np.full((238, 81), np.nan)
tas_var3_solar = np.full((238, 81), np.nan)
tot_var_reg3_tas = np.full((238, 81), np.nan)

for i in range(238):
    for t in range(81):
        x = ta_reg_sulf[:, t, i]
        y = ta_reg_solar[:, t, i]
        tas_var3_sulf[i, t] = np.std(x, ddof=1) ** 2
        tas_var3_solar[i, t] = np.std(y, ddof=1) ** 2
    tot_var_reg3_tas[i, :] = tas_var3_sulf[i, :] + tas_var3_solar[i, :]

# %% [markdown]
# ## Precipitation — the same three resolutions
# The precipitation inputs reuse the same array names `ta_reg_sulf` / `ta_reg_solar`.

# %%
m = mu.loadmat("pr_regions1.mat")
ta_reg_sulf = m["ta_reg_sulf"]
ta_reg_solar = m["ta_reg_solar"]

regions = shaperead("IPCC-WGI-reference-regions-v4.shx")
nR1 = len(regions)

pr_var1_sulf = np.full((nR1, 81), np.nan)
pr_var1_solar = np.full((nR1, 81), np.nan)
tot_var_reg1_pr = np.full((nR1, 81), np.nan)

for i in range(nR1):
    for t in range(81):
        x = ta_reg_sulf[:, t, i]
        y = ta_reg_solar[:, t, i]
        pr_var1_sulf[i, t] = np.std(x, ddof=1) ** 2
        pr_var1_solar[i, t] = np.std(y, ddof=1) ** 2
    tot_var_reg1_pr[i, :] = pr_var1_sulf[i, :] + pr_var1_solar[i, :]

# %%
m = mu.loadmat("pr_regions2.mat")
ta_reg_sulf = m["ta_reg_sulf"]
ta_reg_solar = m["ta_reg_solar"]

regions = shaperead("IPCC_WGII_continental_regions.shx")
nR2 = len(regions)

pr_var2_sulf = np.full((nR2, 81), np.nan)
pr_var2_solar = np.full((nR2, 81), np.nan)
tot_var_reg2_pr = np.full((nR2, 81), np.nan)

for i in range(nR2):
    for t in range(81):
        x = ta_reg_sulf[:, t, i]
        y = ta_reg_solar[:, t, i]
        pr_var2_sulf[i, t] = np.std(x, ddof=1) ** 2
        pr_var2_solar[i, t] = np.std(y, ddof=1) ** 2
    tot_var_reg2_pr[i, :] = pr_var2_sulf[i, :] + pr_var2_solar[i, :]

# %%
m = mu.loadmat("pr_regions3.mat")
ta_reg_sulf = m["ta_reg_sulf"]
ta_reg_solar = m["ta_reg_solar"]

pr_var3_sulf = np.full((238, 81), np.nan)
pr_var3_solar = np.full((238, 81), np.nan)
tot_var_reg3_pr = np.full((238, 81), np.nan)

for i in range(238):
    for t in range(81):
        x = ta_reg_sulf[:, t, i]
        y = ta_reg_solar[:, t, i]
        pr_var3_sulf[i, t] = np.std(x, ddof=1) ** 2
        pr_var3_solar[i, t] = np.std(y, ddof=1) ** 2
    tot_var_reg3_pr[i, :] = pr_var3_sulf[i, :] + pr_var3_solar[i, :]

# %% [markdown]
# ## Fractional sulfur variances and the scatter
# Per-region fractional sulfur variance (%) at each timestep, then the time
# mean over years 41–81 (2060–2100). Black = WGII (8 regions),
# red = WGI (58 regions), blue = fine 238-region grid.

# %%
pr3 = pr_var3_sulf / tot_var_reg3_pr * 100
tas3 = tas_var3_sulf / tot_var_reg3_tas * 100

pr2 = pr_var2_sulf / tot_var_reg2_pr * 100
tas2 = tas_var2_sulf / tot_var_reg2_tas * 100

pr1 = pr_var1_sulf / tot_var_reg1_pr * 100
tas1 = tas_var1_sulf / tot_var_reg1_tas * 100

fig, ax = plt.subplots(figsize=(7, 6))
ax.scatter(tas2[:, 40:81].mean(axis=1), pr2[:, 40:81].mean(axis=1),
           s=250, c="k")
ax.axis([45, 70, 45, 70])
ax.scatter(tas1[:, 40:81].mean(axis=1), pr1[:, 40:81].mean(axis=1),
           s=250, c="r")
ax.scatter(tas3[:, 40:81].mean(axis=1), pr3[:, 40:81].mean(axis=1),
           s=50, c="b")
ax.plot([50, 50], [0, 100], "k")
ax.plot([0, 100], [50, 50], "k")
ax.set_xlabel("tas fractional variance (%)")
ax.set_ylabel("pr fractional variance (%)")
fig.set_facecolor("w")
ax.set_facecolor("w")
plt.show()
