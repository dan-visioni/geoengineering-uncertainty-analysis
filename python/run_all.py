#!/usr/bin/env python
"""
End-to-end acceptance runner for the Python/Jupyter port.

Executes every notebook in DEPENDENCY ORDER with the `geomip` kernel, captures
pass/fail + runtime + first error for each, and prints a markdown table.

Usage (from the project root):
    /opt/anaconda3/envs/geomip/bin/python python/run_all.py [--clean]

--clean  : delete data_py/ and py_outputs/ first, so the run proves the pipeline
           is self-contained (in-folder-generated datasets are rebuilt from raw
           NetCDF; only true external .mat inputs are read from disk).

The notebooks are executed in place (outputs/figures embedded). Order respects
the producer -> consumer .npz dependencies:

  L0 data-prep (raw NetCDF -> data_py/*.npz)
  L1 analyses that consume L0 datasets / external .mat / raw data
  L2 analyses that consume L1-produced .npz (region_std, 2deg_pr, ...)
"""
import os
import sys
import time
import shutil
import subprocess

ENVBIN = "/opt/anaconda3/envs/geomip/bin"
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Dependency-ordered notebook list. (Edit as notebooks are finalized.)
ORDER = [
    # --- L0: data-prep, raw NetCDF -> data_py/*.npz ---
    ("read_GeoMIP_models",      "L0 data-prep"),
    ("save_SOCOL_alternatw",    "L0 data-prep"),
    ("save_SOCOL",              "L0 data-prep"),
    ("plot_maps_15k",           "L0 data-prep (G615kSAI)"),
    # --- L1: consume L0 npz + external .mat + raw ---
    ("read? include_SOCOL_HS",  "L1"),   # placeholder names corrected below
    ("plot_HS_GeoMIP",          "L1"),
    ("plot_HSpr_GeoMIP",        "L1 -> pr_regions1_std"),
    ("plot_HS_countries",       "L1 -> tas_regions1_std"),
    ("plot_maps",               "L1"),
    ("compare_g6s",             "L1 (needs G615kSAI)"),
    ("plot_spatial_aggregation_tas", "L1 -> 2deg_pr"),
    ("plot_HS_TREFHT",          "L1 CESM"),
    ("plot_HS_PRECT",           "L1 CESM"),
    ("plot_AOD",                "L1 CESM -> AODVIS"),
    ("plot_HS_regions_map",     "L1 CESM (dashboard)"),
    ("plot_uncertainty",        "L1"),
    ("plot_AOD_vs_tas",         "L1"),
    # --- L2: consume L1-produced npz ---
    ("plot_tas_vs_pr",          "L2"),
    ("plot_g6vsSSPunc",         "L2 (needs *_regions1_std)"),
    ("plot_g6_uncertainties_aerosols_vs_solar", "L2 (needs *_regions1_std)"),
    ("plot_all_spatial_aggregation", "L2 (needs 2deg_pr)"),
]
# correct the include_SOCOL_HS placeholder
ORDER[4] = ("include_SOCOL_HS", "L1")


def count_errors(path):
    import nbformat
    nb = nbformat.read(path, as_version=4)
    figs = 0
    for c in nb.cells:
        for o in c.get("outputs", []):
            if o.get("output_type") == "error":
                return False, f"{o.get('ename')}: {o.get('evalue')}", figs
            if "image/png" in o.get("data", {}):
                figs += 1
    return True, "", figs


def run_one(name):
    src = os.path.join(ROOT, "python", "src", name + ".py")
    os.makedirs(os.path.join(ROOT, "j_notebooks"), exist_ok=True)
    nb = os.path.join(ROOT, "j_notebooks", name + ".ipynb")
    if not os.path.exists(src):
        return ("MISSING", 0.0, "no source", 0)
    t0 = time.time()
    # rebuild the .ipynb from its percent source (sources are the source of truth)
    c = subprocess.run([f"{ENVBIN}/jupytext", "--to", "ipynb", src, "-o", nb, "--quiet"],
                       capture_output=True, text=True, cwd=ROOT)
    if c.returncode != 0:
        return ("FAIL", time.time() - t0, "jupytext: " + (c.stderr.strip()[-120:]), 0)
    p = subprocess.run(
        [f"{ENVBIN}/jupyter", "nbconvert", "--to", "notebook", "--execute",
         "--inplace", nb, "--ExecutePreprocessor.timeout=2400",
         "--ExecutePreprocessor.kernel_name=geomip"],
        capture_output=True, text=True, cwd=ROOT)
    dt = time.time() - t0
    if p.returncode != 0:
        tail = (p.stderr.strip().splitlines() or ["nbconvert failed"])[-1]
        return ("FAIL", dt, tail[:160], 0)
    ok, err, figs = count_errors(nb)
    return ("PASS" if ok else "FAIL", dt, err[:160], figs)


def main():
    clean = "--clean" in sys.argv
    if clean:
        for d in ("data_py", "py_outputs"):
            shutil.rmtree(os.path.join(ROOT, d), ignore_errors=True)
        print("Cleaned data_py/ and py_outputs/ (in-folder datasets will be rebuilt).\n")

    rows = []
    for name, lvl in ORDER:
        status, dt, err, figs = run_one(name)
        rows.append((name, lvl, status, dt, figs, err))
        print(f"  {status:7s} {name:42s} {dt:6.1f}s  figs={figs}  {err}")

    print("\n| Notebook | Stage | Status | Time (s) | Figures | Error |")
    print("|---|---|---|---:|---:|---|")
    for name, lvl, status, dt, figs, err in rows:
        print(f"| `{name}` | {lvl} | {status} | {dt:.0f} | {figs} | {err} |")

    npass = sum(1 for r in rows if r[2] == "PASS")
    print(f"\n**{npass}/{len(rows)} notebooks PASS.**")
    return 0 if npass == len(rows) else 1


if __name__ == "__main__":
    sys.exit(main())
