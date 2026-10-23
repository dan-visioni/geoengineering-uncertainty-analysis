"""
Plot inter-model variance σ²(t) of temperature anomalies (relative to SSP2-4.5)
for G6sulfur, G6solar, and SSP5-8.5 — global mean and selected IPCC regions.
"""
import sys, numpy as np, matplotlib.pyplot as plt
sys.path.insert(0, '/Users/danielevisioni/Desktop/Articles_and_Research/Articles WIP/Uncertainties/python')
import mat_utils as mu
from mat_utils import movmean

# ── load ──────────────────────────────────────────────────────────────────────
mt  = mu.load_npz('GeoMIP_tas')
x_e = mt['x_e']          # object (6, 4, 3);  inner: (lon, lat, time)
LAT = mt['LAT']           # object (1, 6);     inner: (nlat,)
LON = mt['LON']           # object (1, 6);     inner: (nlon,)

NMOD = 6
YEARS = np.arange(2020, 2101)   # 81 years common to all experiments
NYR = 81

# region bounding boxes [lat_min, lat_max, lon_min, lon_max]  (lon 0-360 convention)
REGIONS = {
    'Global':        None,
    'Africa':        (-35,  37,  340, 360+55),   # lon handled below
    'Europe':        ( 36,  72,  348, 360+42),
    'N. America':    ( 25,  70,  220, 305),
    'SE Asia':       (-10,  30,   65, 150),
}

EXP_IDX  = {'G6sulfur': 0, 'G6solar': 1, 'SSP5-8.5': 2}
BASE_IDX  = 3   # SSP2-4.5

COLORS  = {'G6sulfur': '#1a6faf', 'G6solar': '#e87d0d', 'SSP5-8.5': '#c0392b'}
LABELS  = {'G6sulfur': 'G6sulfur', 'G6solar': 'G6solar', 'SSP5-8.5': 'SSP5-8.5'}

# ── regional weighted mean helper ─────────────────────────────────────────────
def regmean(data3d, lat1d, lon1d, bbox):
    """Cosine-weighted mean of (lon, lat, time) over a bounding box."""
    if bbox is None:
        wt = np.cos(np.deg2rad(lat1d))
        return np.average(data3d.mean(axis=0), axis=0, weights=wt)  # (time,)
    lat_min, lat_max, lon_min, lon_max = bbox
    ln = lon1d % 360                          # normalise to 0-360
    a, b = lon_min % 360, lon_max % 360      # both in 0-360
    lmask = (ln >= a) & (ln <= b) if a <= b else (ln >= a) | (ln <= b)
    jmask = (lat1d >= lat_min) & (lat1d <= lat_max)
    ii, jj = np.where(lmask)[0], np.where(jmask)[0]
    if len(ii) == 0 or len(jj) == 0:
        return np.full(data3d.shape[-1], np.nan)
    sub  = data3d[np.ix_(ii, jj)]             # (nlon_r, nlat_r, time)
    wt   = np.cos(np.deg2rad(lat1d[jj]))
    return np.average(sub.mean(axis=0), axis=0, weights=wt)

# ── build per-model, per-experiment time series (81 yr) ───────────────────────
def get_ts(im, ie, bbox):
    """Ensemble-mean regional/global mean for model im, experiment ie."""
    lat1d = np.array(LAT[0, im])
    lon1d = np.array(LON[0, im])
    ens_ts = []
    for ie2 in range(3):
        arr = np.array(x_e[im, ie, ie2])
        if arr.size == 0:
            continue
        arr = arr[..., -NYR:]           # trim to last 81 years
        ens_ts.append(regmean(arr, lat1d, lon1d, bbox))
    return np.nanmean(ens_ts, axis=0)  # (81,)

# ── compute σ²(t) for each experiment in each region ─────────────────────────
def compute_sigma2(exp_name, bbox):
    ie = EXP_IDX[exp_name]
    anom = []
    for im in range(NMOD):
        ts_exp  = get_ts(im, ie,       bbox)
        ts_base = get_ts(im, BASE_IDX, bbox)
        anom.append(ts_exp - ts_base)
    return np.nanvar(anom, axis=0, ddof=1)   # (81,)

MVM = 10   # 10-year centred moving mean — consistent with all other paper figures

# ── plot ──────────────────────────────────────────────────────────────────────
nreg = len(REGIONS)
fig, axes = plt.subplots(1, nreg, figsize=(3.5*nreg, 3.8), sharey=False)

for ax, (rname, bbox) in zip(axes, REGIONS.items()):
    for exp in ('G6sulfur', 'G6solar', 'SSP5-8.5'):
        s2 = compute_sigma2(exp, bbox)
        ax.plot(YEARS, movmean(s2, MVM), color=COLORS[exp], lw=1.8, label=LABELS[exp])
    ax.set_title(rname, fontsize=10, fontweight='bold')
    ax.set_xlabel('Year', fontsize=9)
    ax.set_xlim(2020, 2100)
    ax.tick_params(labelsize=8)
    ax.grid(alpha=0.25, lw=0.5)
    ax.spines[['top','right']].set_visible(False)

axes[0].set_ylabel('Inter-model variance σ² [K²]', fontsize=9)
axes[-1].legend(fontsize=8, framealpha=0.8)
fig.suptitle('Inter-model temperature variance σ²(t) relative to SSP2-4.5', fontsize=11, y=1.02)
plt.tight_layout()

outpath = '/Users/danielevisioni/Desktop/Articles_and_Research/Articles WIP/Uncertainties/py_outputs/sigma2_timeseries_2.png'
plt.savefig(outpath, dpi=150, bbox_inches='tight')
print(f'Saved: {outpath}')
