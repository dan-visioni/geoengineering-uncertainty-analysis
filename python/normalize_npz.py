"""Normalize data_py/*.npz so each variable's shape mirrors the original .mat,
so that consumer notebooks can switch loadmat->load_npz with NO indexing changes.
Re-saves in place and re-verifies numeric equality against the .mat."""
import sys, os
sys.path.insert(0, "python")
import numpy as np
import mat_utils as mu

PAIRS = [
    ("GeoMIP_tas.mat", "GeoMIP_tas"),
    ("GeoMIP_pr.mat", "GeoMIP_pr"),
    ("data_all_SOCOL.mat", "data_all_SOCOL"),
    ("map_plots_GeoMIP_SOCOL.mat", "map_plots_GeoMIP_SOCOL"),
    ("aod_SOCOL.mat", "aod_SOCOL"),
    ("tas_timeseries_SOCOL.mat", "tas_timeseries_SOCOL"),
]


def reshape_like(npz_val, mat_val):
    """Reshape npz_val to mat_val.shape when element counts match."""
    mv = np.asarray(mat_val)
    nv = np.asarray(npz_val)
    if nv.dtype == object or mv.dtype == object:
        # object/cell arrays: reshape the object container to the .mat layout
        if nv.size == mv.size and nv.shape != mv.shape:
            nv = nv.reshape(mv.shape)
        return nv
    if nv.shape != mv.shape and nv.size == mv.size:
        nv = nv.reshape(mv.shape)
    return nv


for matf, name in PAIRS:
    if not os.path.exists(matf):
        print(f"skip {name}: {matf} absent"); continue
    m = mu.loadmat(matf)
    d = mu.load_npz(name)
    keys = [k for k in m if not k.startswith("__") and k in d]
    out = {}
    changed = []
    for k in keys:
        before = np.asarray(d[k]).shape
        out[k] = reshape_like(d[k], m[k])
        after = np.asarray(out[k]).shape
        if before != after:
            changed.append(f"{k}: {before}->{after}")
    # keep any npz-only keys untouched
    for k in d:
        if k not in out:
            out[k] = d[k]
    mu.save_npz(name, **out)
    print(f"{name}: reshaped {changed if changed else 'none'}")

print("\n=== re-verify numeric max-diff vs .mat ===")
for matf, name in PAIRS:
    if not os.path.exists(matf):
        continue
    m = mu.loadmat(matf)
    d = mu.load_npz(name)
    worst = 0.0
    for k in [k for k in m if not k.startswith("__") and k in d]:
        mv, nv = np.asarray(m[k]), np.asarray(d[k])
        if mv.dtype == object:
            for idx in np.ndindex(mv.shape):
                a, b = np.asarray(mv[idx], float), np.asarray(nv[idx], float)
                if a.size and a.shape == b.shape:
                    worst = max(worst, float(np.nanmax(np.abs(a - b))))
        elif np.issubdtype(mv.dtype, np.number):
            a, b = mv.astype(float), nv.astype(float)
            if a.shape == b.shape:
                worst = max(worst, float(np.nanmax(np.abs(a - b))))
    print(f"  {name}: max|diff| = {worst:.3e}")
