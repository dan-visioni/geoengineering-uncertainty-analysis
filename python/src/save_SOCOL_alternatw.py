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
# # save_SOCOL_alternatw — SOCOL alternative-aerosol-injection ensemble
#
# Reads the SOCOL "alternative injection strategy" surface output from
# `/Volumes/D_Visioni2/SOCOL/GEOMIP_SOCOL/`: 2 m temperature (`var167`) and
# precipitation (large-scale `var142` + convective `var143`) for four G6-style
# injection strategies (G6, equatorial "al", cross-equatorial "ca", deep-tropics
# "di"), each with three members, plus the SSP2-4.5 and SSP5-8.5 references
# (three members each). Monthly fields are collapsed to annual means and saved
# as `data_py/data_all_SOCOL.npz`.

# %%
import sys
import os as _os
_r = _os.getcwd()
while not _os.path.isdir(_os.path.join(_r, "python")) and _os.path.dirname(_r) != _r:
    _r = _os.path.dirname(_r)
_os.chdir(_r)   # notebooks live in j_notebooks/; anchor to project root
sys.path.insert(0, "python")
import numpy as np
import mat_utils as mu

# %% [markdown]
# ## Read the raw monthly fields
# Builds 5-D arrays `ta_g6[lon, lat, time, member, strategy]` where the strategy
# axis order is G6=0, "al"=1, "ca"=2, "di"=3. Temperature uses `var167`;
# precipitation is large-scale `var142` + convective `var143`.

# %%
fold = "/Volumes/D_Visioni2/SOCOL/GEOMIP_SOCOL/"


def rd(tag, var):
    return mu.ncread(f"{fold}GEOMIP_socol4_{tag}_surface.nc", var)


# temperature: ta_g6[:,:,:, member, strategy]  (strategy: g6,al,ca,di -> 0..3)
nlon, nlat, nt = rd("g61", "var167").shape
ta_g6 = np.empty((nlon, nlat, nt, 3, 4))
for m, tag in enumerate(["al1", "al2", "al3"]):
    ta_g6[:, :, :, m, 1] = rd(tag, "var167")
for m, tag in enumerate(["ca1", "ca2", "ca3"]):
    ta_g6[:, :, :, m, 2] = rd(tag, "var167")
for m, tag in enumerate(["di1", "di2", "di3"]):
    ta_g6[:, :, :, m, 3] = rd(tag, "var167")
for m, tag in enumerate(["g61", "g62", "g63"]):
    ta_g6[:, :, :, m, 0] = rd(tag, "var167")

ta_ssp2 = np.stack([rd(t, "var167") for t in ["ssp21", "ssp22", "ssp23"]], axis=3)
ta_ssp5 = np.stack([rd(t, "var167") for t in ["ssp51", "ssp52", "ssp53"]], axis=3)

# %% [markdown]
# Precipitation = large-scale (`var142`) + convective (`var143`). Note: for the
# "al" strategy, all three members share `al1`'s convective term; this is
# preserved as-is from the original experiment setup.

# %%
pr_g6 = np.empty((nlon, nlat, nt, 3, 4))
al143 = rd("al1", "var143")
pr_g6[:, :, :, 0, 1] = rd("al1", "var142") + al143
pr_g6[:, :, :, 1, 1] = rd("al2", "var142") + al143   # MATLAB: + al1 var143
pr_g6[:, :, :, 2, 1] = rd("al3", "var142") + al143   # MATLAB: + al1 var143
for m, tag in enumerate(["ca1", "ca2", "ca3"]):
    pr_g6[:, :, :, m, 2] = rd(tag, "var142") + rd(tag, "var143")
for m, tag in enumerate(["di1", "di2", "di3"]):
    pr_g6[:, :, :, m, 3] = rd(tag, "var142") + rd(tag, "var143")
for m, tag in enumerate(["g61", "g62", "g63"]):
    pr_g6[:, :, :, m, 0] = rd(tag, "var142") + rd(tag, "var143")

pr_ssp2 = np.stack([rd(t, "var142") + rd(t, "var143")
                    for t in ["ssp21", "ssp22", "ssp23"]], axis=3)
pr_ssp5 = np.stack([rd(t, "var142") + rd(t, "var143")
                    for t in ["ssp51", "ssp52", "ssp53"]], axis=3)

# %% [markdown]
# ## Coordinates

# %%
lat_SOCOL = rd("ssp21", "lat")
lon_SOCOL = rd("ssp21", "lon")

# %% [markdown]
# ## Annual means
# Reshape the month axis into (12, 70 years) and average the 12 months,
# using Fortran-order reshape to preserve the original time-axis layout.

# %%
def annual(X, extra):
    """Reshape (192,96,840,*extra) -> (192,96,12,70,*extra) and mean over months."""
    Y = np.reshape(X, (192, 96, 12, 70, *extra), order="F")
    return np.squeeze(np.mean(Y, axis=2))


pr_g6_y = annual(pr_g6, (3, 4))
ta_g6_y = annual(ta_g6, (3, 4))
pr_ssp5_y = annual(pr_ssp5, (3,))
ta_ssp5_y = annual(ta_ssp5, (3,))
pr_ssp2_y = annual(pr_ssp2, (3,))
ta_ssp2_y = annual(ta_ssp2, (3,))

mu.save_npz("data_all_SOCOL",
            pr_g6_y=pr_g6_y, ta_g6_y=ta_g6_y,
            ta_ssp5_y=ta_ssp5_y, ta_ssp2_y=ta_ssp2_y,
            pr_ssp5_y=pr_ssp5_y, pr_ssp2_y=pr_ssp2_y,
            lat_SOCOL=lat_SOCOL, lon_SOCOL=lon_SOCOL)

# %% [markdown]
# ## Verification
# The reference arrays are stored as float32, so float32-eps round-off is
# expected (precipitation/temperature magnitudes give ~1e-3 absolute on float32).

# %%
ref = mu.loadmat("data_all_SOCOL.mat")
arrays = {
    "pr_g6_y": pr_g6_y, "ta_g6_y": ta_g6_y,
    "ta_ssp5_y": ta_ssp5_y, "ta_ssp2_y": ta_ssp2_y,
    "pr_ssp5_y": pr_ssp5_y, "pr_ssp2_y": pr_ssp2_y,
    "lat_SOCOL": lat_SOCOL, "lon_SOCOL": lon_SOCOL,
}
ok = True
for name, a in arrays.items():
    a = np.asarray(a, dtype=float)
    b = np.asarray(ref[name], dtype=float)
    if name in ("lat_SOCOL", "lon_SOCOL"):
        a = a.ravel(); b = b.ravel()
    if a.shape != b.shape:
        print(f"  SHAPE MISMATCH {name}: py={a.shape} ref={b.shape}")
        ok = False
        continue
    d = np.nanmax(np.abs(a - b))
    # relative scale for the tolerance (float32-eps * field magnitude)
    scale = max(1.0, np.nanmax(np.abs(b)))
    tol = 1e-5 * scale
    status = "ok" if d < tol else "DIFF"
    if d >= tol:
        ok = False
    print(f"  {name}: shape {a.shape}, max|Δ| = {d:.3e}  ({status}, tol {tol:.1e})")

print()
print("SUMMARY: data_all_SOCOL", "PASS" if ok else "FAIL")
assert ok, "verification failed"
