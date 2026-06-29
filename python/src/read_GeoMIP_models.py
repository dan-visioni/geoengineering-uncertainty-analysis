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
# # read_GeoMIP_models — build the GeoMIP multi-model ensemble datasets
#
# Walks the raw CMIP6/GeoMIP NetCDF archive on `/Volumes/D_Visioni2/GeoMIP/`,
# reads monthly fields for the four experiments (`G6sulfur`, `G6solar`, `ssp585`,
# `ssp245`) across six models (CNRM-ESM2-1, IPSL-CM6A-LR, CESM2-WACCM,
# UKESM1-0-LL, MPI-ESM1-2-LR, MPI-ESM1-2-HR), applies model-specific time-axis
# fix-ups, and collapses each record to **annual** means.
#
# The loader runs twice — once for `tas` (K) and once for `pr`
# (kg m⁻² s⁻¹ → mm day⁻¹, factor 86400) — and writes `data_py/GeoMIP_tas.npz`
# and `data_py/GeoMIP_pr.npz`. These datasets feed all downstream regional and
# global analyses.

# %%
import sys
import os as _os
_r = _os.getcwd()
while not _os.path.isdir(_os.path.join(_r, "python")) and _os.path.dirname(_r) != _r:
    _r = _os.path.dirname(_r)
_os.chdir(_r)   # notebooks live in j_notebooks/; anchor to project root
sys.path.insert(0, "python")
import glob
import numpy as np
import mat_utils as mu

# %% [markdown]
# ## Set-up
# Experiment/model lists, the ensemble-member table `n_ens` (rows per model,
# columns per experiment), and the unit-conversion factors (tas: ×1, pr: ×86400).

# %%
fold = "/Volumes/D_Visioni2/GeoMIP/"
n_exp = ["G6sulfur", "G6solar", "ssp585", "ssp245"]
n_mod = ["CNRM-ESM2-1", "IPSL-CM6A-LR", "CESM2-WACCM",
         "UKESM1-0-LL", "MPI-ESM1-2-LR", "MPI-ESM1-2-HR"]

# ensemble members available per (model, experiment); MATLAB n_ens{im,ie}
n_ens = [
    [[1, 2, 3], [1],       [1, 2, 3], [1, 2, 3]],   # CNRM-ESM2-1
    [[1],       [1],       [1],       [1]],          # IPSL-CM6A-LR
    [[1, 2],    [1, 2],    [1],       [1]],          # CESM2-WACCM
    [[1, 4, 8], [1, 4, 8], [1, 4, 8], [1, 4, 8]],    # UKESM1-0-LL
    [[1, 2, 3], [1, 2, 3], [1, 2, 3], [1, 2, 3]],    # MPI-ESM1-2-LR
    [[1, 2, 3], [1, 2, 3], [1, 2],    [1, 2]],       # MPI-ESM1-2-HR
]

# variable names and the unit-conversion factors (tas: 1, pr: 86400)
var1 = ["tas", "pr"]
corr = [1, 86400]


# %% [markdown]
# ## The loader
# A single function runs the double loop (models × experiments) for one variable. It
# returns the `(6, 4, 3)` object array `x_e` and the length-6 `LAT`/`LON` object
# arrays. Each model is processed independently and its raw fields are released
# before moving on, so the (large) `pr` run never holds more than one model at a
# time.

# %%
def load_geomip(iv):
    """Build x_e / LAT / LON for variable index ``iv`` (0=tas, 1=pr).

    Faithful translation of the MATLAB loop in read_GeoMIP_models.m."""
    vname = var1[iv]
    cfac = corr[iv]

    x_e = np.empty((6, 4, 3), dtype=object)
    for idx in np.ndindex(x_e.shape):
        x_e[idx] = np.zeros((0, 0))          # MATLAB [] (unassigned cell)
    LAT = np.empty(6, dtype=object)
    LON = np.empty(6, dtype=object)

    print("Model\tens size")
    for im in range(6):                                   # MATLAB im=1:6
        print(n_mod[im], end="\t(\t)\t")
        last_list = []
        for ie in range(4):                               # MATLAB ie=1:4
            print(n_exp[ie], end=" ")
            nne = n_ens[im][ie]
            for i in range(len(nne)):                     # MATLAB i=1:length(nne)
                pat = (f"{fold}{n_exp[ie]}/{n_mod[im]}/{vname}_*"
                       f"{n_mod[im]}*{n_exp[ie]}_r{nne[i]}*.nc")
                lst = sorted(glob.glob(pat))              # MATLAB dir(...) sorted
                if len(lst) == 0:
                    x = np.array(0.0)
                    if i == 0:
                        print("(N)", end=" ")
                else:
                    if i == 0:
                        print("(Y)", end=" ")
                    x = None
                    for j in range(len(lst)):             # MATLAB j=1:length(list)
                        x2 = mu.ncread(lst[j], vname) * cfac
                        # CESM2-WACCM: weird beginning length, fix first file.
                        # MATLAB: im==3 && j==1 && ie<3  ->  im==2, j==0, ie<=2
                        if im == 2 and j == 0 and ie <= 2:
                            if ie == 1 and i == 0:        # MATLAB ie==2 && i==1
                                x2 = x2[:, :, 12:]        # x2(:,:,13:end)
                            if i == 1:                    # MATLAB i==2
                                x2 = x2[:, :, 60:]        # x2(:,:,61:end)
                        if len(lst) > 1 and j > 0:
                            x = np.concatenate([x, x2], axis=2)
                        else:
                            x = x2

                # MPI-ESM1-2-LR (im==5, ie<3): pad a missing final year
                if im == 4 and ie <= 2:
                    if x.shape[2] == 1020:
                        x = np.concatenate([x, x[:, :, 1008:1020]], axis=2)
                # MPI-ESM1-2-HR (im==6)
                if im == 5 and ie <= 1:                   # MATLAB im==6 && ie<2
                    if x.shape[2] == 960:
                        x = np.concatenate([x, x[:, :, 948:960]], axis=2)
                if im == 5 and ie >= 1:                   # MATLAB im==6 && ie>=2
                    if x.shape[2] == 1020:
                        x = np.concatenate([x, x[:, :, 1008:1020]], axis=2)

                # reshape (lon,lat,12,nyears), mean over the months axis
                nlon, nlat, nt = x.shape
                x = np.reshape(x, (nlon, nlat, 12, nt // 12), order="F")
                x = np.squeeze(np.mean(x, axis=2))
                x_e[im, ie, i] = x

            if len(lst) > 0:
                last_list = lst
        if last_list:
            LAT[im] = mu.ncread(last_list[0], "lat")
            LON[im] = mu.ncread(last_list[0], "lon")
        print()
    return x_e, LAT, LON


# %% [markdown]
# ## Run for `tas` and save `GeoMIP_tas`

# %%
x_e_tas, LAT_tas, LON_tas = load_geomip(0)
# save LAT/LON as (1,6) object arrays to mirror the MATLAB cell layout that the
# consumer notebooks index as LAT[0, im]
mu.save_npz("GeoMIP_tas", x_e=x_e_tas, LAT=LAT_tas.reshape(1, 6), LON=LON_tas.reshape(1, 6))

# %% [markdown]
# ## Run for `pr` and save `GeoMIP_pr`
# `corr = 86400` converts kg m⁻² s⁻¹ to mm day⁻¹.

# %%
x_e_pr, LAT_pr, LON_pr = load_geomip(1)
mu.save_npz("GeoMIP_pr", x_e=x_e_pr, LAT=LAT_pr.reshape(1, 6), LON=LON_pr.reshape(1, 6))

# %% [markdown]
# ## Verification
# For every non-empty array `x_e[im, ie, i]`, asserts the shape matches the
# reference and reports the maximum absolute difference (expected ~0).
# The `LAT`/`LON` coordinate arrays are checked the same way.

# %%
def verify_geomip(name, x_e, LAT, LON):
    ref = mu.loadmat(name + ".mat")
    xe_r = ref["x_e"]
    LAT_r = ref["LAT"]
    LON_r = ref["LON"]
    max_diff = 0.0
    shape_ok = True
    ncells = 0
    for im in range(6):
        for ie in range(4):
            for i in range(3):
                a = np.asarray(x_e[im, ie, i])
                b = np.asarray(xe_r[im, ie, i])
                a_empty = a.size == 0
                b_empty = b.size == 0
                if a_empty != b_empty:
                    shape_ok = False
                    print(f"  EMPTY MISMATCH at [{im},{ie},{i}]: "
                          f"py={a.shape} ref={b.shape}")
                    continue
                if a_empty:
                    continue
                ncells += 1
                if a.shape != b.shape:
                    shape_ok = False
                    print(f"  SHAPE MISMATCH at [{im},{ie},{i}]: "
                          f"py={a.shape} ref={b.shape}")
                    continue
                d = np.nanmax(np.abs(a - b))
                max_diff = max(max_diff, d)
    # LAT / LON
    coord_diff = 0.0
    for im in range(6):
        for arr, ref_arr in [(LAT, LAT_r), (LON, LON_r)]:
            a = np.asarray(arr[im]).ravel()
            b = np.asarray(ref_arr[0, im]).ravel()
            if a.shape != b.shape:
                shape_ok = False
                print(f"  COORD SHAPE MISMATCH model {im}: {a.shape} vs {b.shape}")
                continue
            coord_diff = max(coord_diff, np.nanmax(np.abs(a - b)))
    status = "PASS" if (shape_ok and max_diff < 1e-6 and coord_diff < 1e-6) else "FAIL"
    print(f"[{name}] {status}: {ncells} non-empty cells, "
          f"max|Δ x_e| = {max_diff:.3e}, max|Δ LAT/LON| = {coord_diff:.3e}")
    return status == "PASS"


ok_tas = verify_geomip("GeoMIP_tas", x_e_tas, LAT_tas, LON_tas)
ok_pr = verify_geomip("GeoMIP_pr", x_e_pr, LAT_pr, LON_pr)
print()
print("SUMMARY:",
      "GeoMIP_tas", "PASS" if ok_tas else "FAIL",
      "| GeoMIP_pr", "PASS" if ok_pr else "FAIL")
assert ok_tas and ok_pr, "verification failed"
