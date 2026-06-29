"""
mat_utils.py
============
Helper functions that reproduce the MATLAB-specific behaviour used by the
*Uncertainties* analysis scripts, so the Python/Jupyter port stays faithful to
the original MATLAB code.

Implemented here:

* ``brewermap(N, scheme)``      -- faithful port of Stephen Cobeldick's
                                   brewermap.m (ColorBrewer schemes, Lab-space
                                   cubic-spline interpolation, ``*`` reversal,
                                   qualitative repetition).
* ``brighten(rgb, beta)``       -- MATLAB ``brighten`` (gamma = (1-b)/(1+b)).
* ``movmean(x, k, ...)``        -- MATLAB ``movmean`` (centred window, even
                                   windows lean left, shrinking edge windows).
* ``ncread(file, var)``         -- MATLAB ``ncread`` dimension order (full axis
                                   reversal) with fill values -> NaN.
* ``loadmat(file)``             -- thin scipy.io.loadmat wrapper.
* ``cell(obj)`` / ``mlcell``    -- index MATLAB cell arrays loaded by scipy.
* ``shaperead(file)``           -- MATLAB ``shaperead`` (records with NaN-
                                   separated .X/.Y plus every dbf attribute).
* ``inpolygon(xq, yq, xv, yv)`` -- MATLAB ``inpolygon`` (NaN-separated multipart
                                   polygons, boundary points included).
* ``interp1`` / ``interp2``     -- MATLAB-compatible linear interpolation with
                                   an explicit extrapolation value.
* cartopy map helpers:
    ``worldmap`` / ``pcolorm`` / ``contourfm`` / ``geoshow_coastlines``.

The numeric helpers (brewermap, brighten, movmean, ncread, interp*) carry no
heavyweight dependencies; the mapping and shapefile helpers import cartopy /
pyshp lazily so the module is usable even where those are absent.
"""

import os
import numpy as np
import scipy.io as sio
from scipy.interpolate import CubicSpline, interp1d, RegularGridInterpolator


# ---------------------------------------------------------------------------
# ColorBrewer node tables (RGB 0-255), copied verbatim from brewermap.m
# ---------------------------------------------------------------------------
_BREWER = {
    # Diverging (num = 11)
    'BrBG': ('Diverging', [[84,48,5],[140,81,10],[166,97,26],[191,129,45],[216,179,101],[223,194,125],[246,232,195],[245,245,245],[199,234,229],[128,205,193],[90,180,172],[53,151,143],[1,133,113],[1,102,94],[0,60,48]]),
    'PiYG': ('Diverging', [[142,1,82],[197,27,125],[208,28,139],[222,119,174],[233,163,201],[241,182,218],[253,224,239],[247,247,247],[230,245,208],[184,225,134],[161,215,106],[127,188,65],[77,172,38],[77,146,33],[39,100,25]]),
    'PRGn': ('Diverging', [[64,0,75],[118,42,131],[123,50,148],[153,112,171],[175,141,195],[194,165,207],[231,212,232],[247,247,247],[217,240,211],[166,219,160],[127,191,123],[90,174,97],[0,136,55],[27,120,55],[0,68,27]]),
    'PuOr': ('Diverging', [[127,59,8],[179,88,6],[230,97,1],[224,130,20],[241,163,64],[253,184,99],[254,224,182],[247,247,247],[216,218,235],[178,171,210],[153,142,195],[128,115,172],[94,60,153],[84,39,136],[45,0,75]]),
    'RdBu': ('Diverging', [[103,0,31],[178,24,43],[202,0,32],[214,96,77],[239,138,98],[244,165,130],[253,219,199],[247,247,247],[209,229,240],[146,197,222],[103,169,207],[67,147,195],[5,113,176],[33,102,172],[5,48,97]]),
    'RdGy': ('Diverging', [[103,0,31],[178,24,43],[202,0,32],[214,96,77],[239,138,98],[244,165,130],[253,219,199],[255,255,255],[224,224,224],[186,186,186],[153,153,153],[135,135,135],[64,64,64],[77,77,77],[26,26,26]]),
    'RdYlBu': ('Diverging', [[165,0,38],[215,48,39],[215,25,28],[244,109,67],[252,141,89],[253,174,97],[254,224,144],[255,255,191],[224,243,248],[171,217,233],[145,191,219],[116,173,209],[44,123,182],[69,117,180],[49,54,149]]),
    'RdYlGn': ('Diverging', [[165,0,38],[215,48,39],[215,25,28],[244,109,67],[252,141,89],[253,174,97],[254,224,139],[255,255,191],[217,239,139],[166,217,106],[145,207,96],[102,189,99],[26,150,65],[26,152,80],[0,104,55]]),
    'Spectral': ('Diverging', [[158,1,66],[213,62,79],[215,25,28],[244,109,67],[252,141,89],[253,174,97],[254,224,139],[255,255,191],[230,245,152],[171,221,164],[153,213,148],[102,194,165],[43,131,186],[50,136,189],[94,79,162]]),
    # Qualitative (num = #nodes)
    'Accent': ('Qualitative', [[127,201,127],[190,174,212],[253,192,134],[255,255,153],[56,108,176],[240,2,127],[191,91,23],[102,102,102]]),
    'Dark2': ('Qualitative', [[27,158,119],[217,95,2],[117,112,179],[231,41,138],[102,166,30],[230,171,2],[166,118,29],[102,102,102]]),
    'Paired': ('Qualitative', [[166,206,227],[31,120,180],[178,223,138],[51,160,44],[251,154,153],[227,26,28],[253,191,111],[255,127,0],[202,178,214],[106,61,154],[255,255,153],[177,89,40]]),
    'Pastel1': ('Qualitative', [[251,180,174],[179,205,227],[204,235,197],[222,203,228],[254,217,166],[255,255,204],[229,216,189],[253,218,236],[242,242,242]]),
    'Pastel2': ('Qualitative', [[179,226,205],[253,205,172],[203,213,232],[244,202,228],[230,245,201],[255,242,174],[241,226,204],[204,204,204]]),
    'Set1': ('Qualitative', [[228,26,28],[55,126,184],[77,175,74],[152,78,163],[255,127,0],[255,255,51],[166,86,40],[247,129,191],[153,153,153]]),
    'Set2': ('Qualitative', [[102,194,165],[252,141,98],[141,160,203],[231,138,195],[166,216,84],[255,217,47],[229,196,148],[179,179,179]]),
    'Set3': ('Qualitative', [[141,211,199],[255,255,179],[190,186,218],[251,128,114],[128,177,211],[253,180,98],[179,222,105],[252,205,229],[217,217,217],[188,128,189],[204,235,197],[255,237,111]]),
    # Sequential (num = 9)
    'Blues': ('Sequential', [[247,251,255],[239,243,255],[222,235,247],[198,219,239],[189,215,231],[158,202,225],[107,174,214],[66,146,198],[49,130,189],[33,113,181],[8,81,156],[8,69,148],[8,48,107]]),
    'BuGn': ('Sequential', [[247,252,253],[237,248,251],[229,245,249],[204,236,230],[178,226,226],[153,216,201],[102,194,164],[65,174,118],[44,162,95],[35,139,69],[0,109,44],[0,88,36],[0,68,27]]),
    'BuPu': ('Sequential', [[247,252,253],[237,248,251],[224,236,244],[191,211,230],[179,205,227],[158,188,218],[140,150,198],[140,107,177],[136,86,167],[136,65,157],[129,15,124],[110,1,107],[77,0,75]]),
    'GnBu': ('Sequential', [[247,252,240],[240,249,232],[224,243,219],[204,235,197],[186,228,188],[168,221,181],[123,204,196],[78,179,211],[67,162,202],[43,140,190],[8,104,172],[8,88,158],[8,64,129]]),
    'Greens': ('Sequential', [[247,252,245],[237,248,233],[229,245,224],[199,233,192],[186,228,179],[161,217,155],[116,196,118],[65,171,93],[49,163,84],[35,139,69],[0,109,44],[0,90,50],[0,68,27]]),
    'Greys': ('Sequential', [[255,255,255],[247,247,247],[240,240,240],[217,217,217],[204,204,204],[189,189,189],[150,150,150],[115,115,115],[99,99,99],[82,82,82],[37,37,37],[37,37,37],[0,0,0]]),
    'OrRd': ('Sequential', [[255,247,236],[254,240,217],[254,232,200],[253,212,158],[253,204,138],[253,187,132],[252,141,89],[239,101,72],[227,74,51],[215,48,31],[179,0,0],[153,0,0],[127,0,0]]),
    'Oranges': ('Sequential', [[255,245,235],[254,237,222],[254,230,206],[253,208,162],[253,190,133],[253,174,107],[253,141,60],[241,105,19],[230,85,13],[217,72,1],[166,54,3],[140,45,4],[127,39,4]]),
    'PuBu': ('Sequential', [[255,247,251],[241,238,246],[236,231,242],[208,209,230],[189,201,225],[166,189,219],[116,169,207],[54,144,192],[43,140,190],[5,112,176],[4,90,141],[3,78,123],[2,56,88]]),
    'PuBuGn': ('Sequential', [[255,247,251],[246,239,247],[236,226,240],[208,209,230],[189,201,225],[166,189,219],[103,169,207],[54,144,192],[28,144,153],[2,129,138],[1,108,89],[1,100,80],[1,70,54]]),
    'PuRd': ('Sequential', [[247,244,249],[241,238,246],[231,225,239],[212,185,218],[215,181,216],[201,148,199],[223,101,176],[231,41,138],[221,28,119],[206,18,86],[152,0,67],[145,0,63],[103,0,31]]),
    'Purples': ('Sequential', [[252,251,253],[242,240,247],[239,237,245],[218,218,235],[203,201,226],[188,189,220],[158,154,200],[128,125,186],[117,107,177],[106,81,163],[84,39,143],[74,20,134],[63,0,125]]),
    'RdPu': ('Sequential', [[255,247,243],[254,235,226],[253,224,221],[252,197,192],[251,180,185],[250,159,181],[247,104,161],[221,52,151],[197,27,138],[174,1,126],[122,1,119],[122,1,119],[73,0,106]]),
    'Reds': ('Sequential', [[255,245,240],[254,229,217],[254,224,210],[252,187,161],[252,174,145],[252,146,114],[251,106,74],[239,59,44],[222,45,38],[203,24,29],[165,15,21],[153,0,13],[103,0,13]]),
    'YlGn': ('Sequential', [[255,255,229],[255,255,204],[247,252,185],[217,240,163],[194,230,153],[173,221,142],[120,198,121],[65,171,93],[49,163,84],[35,132,67],[0,104,55],[0,90,50],[0,69,41]]),
    'YlGnBu': ('Sequential', [[255,255,217],[255,255,204],[237,248,177],[199,233,180],[161,218,180],[127,205,187],[65,182,196],[29,145,192],[44,127,184],[34,94,168],[37,52,148],[12,44,132],[8,29,88]]),
    'YlOrBr': ('Sequential', [[255,255,229],[255,255,212],[255,247,188],[254,227,145],[254,217,142],[254,196,79],[254,153,41],[236,112,20],[217,95,14],[204,76,2],[153,52,4],[140,45,4],[102,37,6]]),
    'YlOrRd': ('Sequential', [[255,255,204],[255,255,178],[255,237,160],[254,217,118],[254,204,92],[254,178,76],[253,141,60],[252,78,42],[240,59,32],[227,26,28],[189,0,38],[177,0,38],[128,0,38]]),
}

# number of defining nodes by type (brewermap.m bmColors trailer)
_NUM = {'Diverging': 11, 'Qualitative': None, 'Sequential': 9}


def _bm_index(N, num, typ):
    """Port of brewermap.m::bmIndex -> (idx 0-based, interpolate flag)."""
    itp = N > num
    if typ == 'Qualitative':
        return (np.mod(np.arange(N), num)).astype(int), False
    if typ == 'Diverging':
        table = {1: [8], 2: [4, 12], 3: [5, 8, 11], 4: [3, 6, 10, 13],
                 5: [3, 6, 8, 10, 13], 6: [2, 5, 7, 9, 11, 14],
                 7: [2, 5, 7, 8, 9, 11, 14], 8: [2, 4, 6, 7, 9, 10, 12, 14],
                 9: [2, 4, 6, 7, 8, 9, 10, 12, 14],
                 10: [1, 2, 4, 6, 7, 9, 10, 12, 14, 15]}
        idx = table.get(N, [1, 2, 4, 6, 7, 8, 9, 10, 12, 14, 15])
    elif typ == 'Sequential':
        table = {1: [6], 2: [4, 8], 3: [3, 6, 9], 4: [2, 5, 7, 10],
                 5: [2, 5, 7, 9, 11], 6: [2, 4, 6, 7, 9, 11],
                 7: [2, 4, 6, 7, 8, 10, 12], 8: [1, 3, 4, 6, 7, 8, 10, 12]}
        idx = table.get(N, [1, 3, 4, 6, 7, 8, 10, 11, 13])
    else:
        raise ValueError(f'unknown colorscheme type {typ}')
    return np.asarray(idx, dtype=int) - 1, itp   # MATLAB 1-based -> 0-based


# sRGB <-> Lab transforms, identical to brewermap.m
_M = np.array([[3.2406255, -1.5372080, -0.4986286],
               [-0.9689307, 1.8757561, 0.0415175],
               [0.0557101, -0.2040211, 1.0569959]])
_WPT = np.array([0.95047, 1.0, 1.08883])  # D65


def _gamma_inv(rgb):
    out = np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
    return np.real(out)


def _gamma_cor(rgb):
    out = np.where(rgb <= 0.0031308, 12.92 * rgb, 1.055 * np.abs(rgb) ** (1 / 2.4) - 0.055)
    return np.real(out)


def _rgb2lab(rgb):
    xyz = _gamma_inv(rgb) @ np.linalg.inv(_M.T)
    xyz = xyz / _WPT
    idx = xyz > (6 / 29) ** 3
    F = np.where(idx, np.cbrt(np.where(idx, xyz, 1.0)),
                 xyz * (29 / 6) ** 2 / 3 + 4 / 29)
    lab = np.empty_like(rgb)
    lab[:, 1] = 500 * (F[:, 0] - F[:, 1])
    lab[:, 2] = 200 * (F[:, 1] - F[:, 2])
    lab[:, 0] = 116 * F[:, 1] - 16
    return lab


def _lab2rgb(lab):
    tmp = lab[:, [1, 0, 2]] / np.array([500.0, np.inf, -200.0])
    tmp = tmp + (lab[:, 0:1] + 16) / 116
    idx = tmp > (6 / 29)
    tmp = np.where(idx, tmp ** 3, 3 * (6 / 29) ** 2 * (tmp - 4 / 29))
    xyz = tmp * _WPT
    rgb = _gamma_cor(xyz @ _M.T)
    return np.clip(rgb, 0.0, 1.0)


def brewermap(N, scheme):
    """Faithful port of brewermap.m (RGB colormap, values in [0, 1]).

    Prefix ``scheme`` with ``*`` to reverse, e.g. ``brewermap(11, '*RdBu')``.
    """
    reverse = scheme.startswith('*')
    name = scheme[1:] if reverse else scheme
    typ, rgb255 = _BREWER[name]
    rgb_nodes = np.asarray(rgb255, dtype=float)
    num = _NUM[typ] if _NUM[typ] is not None else rgb_nodes.shape[0]
    N = int(N)
    if N == 0:
        return np.zeros((0, 3))
    idx, itp = _bm_index(N, num, typ)
    cmap = rgb_nodes[idx, :] / 255.0
    if itp:
        lab = _rgb2lab(cmap)
        xs = np.arange(1, num + 1)
        xq = np.linspace(1, num, N)
        # MATLAB interp1(...,'spline') == not-a-knot cubic spline
        out = np.empty((N, 3))
        for c in range(3):
            out[:, c] = CubicSpline(xs, lab[:, c], bc_type='not-a-knot')(xq)
        cmap = _lab2rgb(out)
    if reverse:
        cmap = cmap[::-1, :]
    return cmap


def brighten(rgb, beta):
    """MATLAB ``brighten``: newmap = map .^ gamma, gamma = (1-beta)/(1+beta)."""
    rgb = np.asarray(rgb, dtype=float)
    beta = float(np.clip(beta, -1 + np.finfo(float).eps, 1 - np.finfo(float).eps))
    gamma = (1.0 - beta) / (1.0 + beta)
    return np.clip(rgb, 0.0, 1.0) ** gamma


# ---------------------------------------------------------------------------
# movmean -- MATLAB moving mean
# ---------------------------------------------------------------------------
def _movmean_1d(x, k, omitnan):
    x = np.asarray(x, dtype=float)
    n = x.size
    left = k // 2            # even windows lean left (more trailing elements)
    right = k - 1 - left
    out = np.empty(n, dtype=float)
    for i in range(n):
        lo = max(0, i - left)
        hi = min(n, i + right + 1)
        seg = x[lo:hi]
        if omitnan:
            m = seg[~np.isnan(seg)]
            out[i] = np.nan if m.size == 0 else m.mean()
        else:
            out[i] = seg.mean()      # NaN propagates, matching MATLAB default
    return out


def movmean(x, k, omitnan=False, axis=None):
    """MATLAB ``movmean(x, k)``.

    Centred moving mean of window length ``k``.  For even ``k`` the window
    includes ``k/2`` trailing and ``k/2 - 1`` leading neighbours (MATLAB
    convention).  Edge windows shrink.  NaNs propagate unless ``omitnan``.
    By default operates along the first non-singleton dimension (as MATLAB).
    """
    x = np.asarray(x, dtype=float)
    if x.ndim == 1:
        return _movmean_1d(x, k, omitnan)
    if axis is None:
        axis = next((d for d, s in enumerate(x.shape) if s > 1), 0)
    return np.apply_along_axis(lambda v: _movmean_1d(v, k, omitnan), axis, x)


# ---------------------------------------------------------------------------
# NetCDF reading with MATLAB dimension order
# ---------------------------------------------------------------------------
def ncread(filename, varname):
    """MATLAB ``ncread``: returns the variable with axes fully reversed
    relative to the CDL/Python order (MATLAB is column-major), scale/offset
    applied, and fill values converted to NaN."""
    import netCDF4
    with netCDF4.Dataset(filename) as ds:
        var = ds.variables[varname]
        var.set_auto_maskandscale(True)
        data = var[:]
    if np.ma.isMaskedArray(data):
        data = data.filled(np.nan).astype(float)
    else:
        data = np.asarray(data, dtype=float)
    return np.transpose(data)            # reverse all axes -> MATLAB order


# ---------------------------------------------------------------------------
# .mat loading + cell-array access
# ---------------------------------------------------------------------------
def loadmat(filename, **kw):
    """scipy.io.loadmat wrapper (struct_as_record=False, squeeze_me=False by
    default so MATLAB shapes such as cell arrays load as object ndarrays)."""
    return sio.loadmat(filename, **kw)


DATA_DIR = "data_py"   # Python-native datasets generated by notebooks in this folder


def save_npz(name, **arrays):
    """Save a dataset generated by an in-folder notebook to ``data_py/<name>.npz``.

    Replaces MATLAB ``save('<name>.mat', ...)`` for files that are *created by a
    program in this folder* (the original ``.mat`` files are left untouched).
    Ragged/cell data may be passed as object arrays (stored with pickling)."""
    os.makedirs(DATA_DIR, exist_ok=True)
    path = os.path.join(DATA_DIR, name + ".npz")
    np.savez(path, **arrays)
    return path


def load_npz(name):
    """Load a dataset previously written with :func:`save_npz`.

    Returns a dict of arrays.  Replaces MATLAB ``load('<name>.mat')`` /
    ``mu.loadmat`` for files created by a program in this folder."""
    path = os.path.join(DATA_DIR, name + ".npz")
    with np.load(path, allow_pickle=True) as d:
        return {k: d[k] for k in d.files}


def cell(obj, *idx):
    """Index a MATLAB cell array (object ndarray) the MATLAB way.

    ``cell(LAT, im)`` for a 1xN cell, ``cell(x_e, im, ie, inn)`` for ND.
    Indices are 0-based.
    """
    arr = np.asarray(obj, dtype=object)
    if len(idx) == 1 and arr.ndim == 2 and 1 in arr.shape:
        flat = arr.ravel()
        return flat[idx[0]]
    return arr[idx]


# ---------------------------------------------------------------------------
# Shapefile reading + point-in-polygon, mirroring shaperead / inpolygon
# ---------------------------------------------------------------------------
class _Record(dict):
    """Attribute-accessible shapefile record (mirrors MATLAB struct).

    ``.X`` / ``.Y`` are NaN-separated coordinate vectors (lon in X, lat in Y),
    exactly as returned by MATLAB ``shaperead``."""
    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError as e:
            raise AttributeError(name) from e


def shaperead(filename):
    """MATLAB ``shaperead``: list of records with NaN-separated .X/.Y plus all
    dbf attributes as fields."""
    import shapefile  # pyshp
    sf = shapefile.Reader(filename)
    fields = [f[0] for f in sf.fields[1:]]
    records = []
    for sr in sf.shapeRecords():
        shp = sr.shape
        pts = np.asarray(shp.points, dtype=float)
        parts = list(shp.parts) + [len(pts)]
        xs, ys = [], []
        for p in range(len(parts) - 1):
            seg = pts[parts[p]:parts[p + 1]]
            xs.append(seg[:, 0]); xs.append([np.nan])
            ys.append(seg[:, 1]); ys.append([np.nan])
        rec = _Record()
        rec['X'] = np.concatenate(xs) if xs else np.array([])
        rec['Y'] = np.concatenate(ys) if ys else np.array([])
        rec['Geometry'] = shp.shapeTypeName
        for fname, val in zip(fields, sr.record):
            rec[fname] = val
        records.append(rec)
    return records


def inpolygon(xq, yq, xv, yv):
    """MATLAB ``inpolygon``: boolean mask of query points (xq, yq) inside the
    polygon (xv, yv).  Handles NaN-separated multipart polygons; boundary
    points count as inside."""
    from matplotlib.path import Path
    xq = np.asarray(xq, dtype=float)
    yq = np.asarray(yq, dtype=float)
    pts = np.column_stack([xq.ravel(), yq.ravel()])
    xv = np.asarray(xv, dtype=float).ravel()
    yv = np.asarray(yv, dtype=float).ravel()

    verts, codes = [], []
    start = 0
    n = xv.size
    for i in range(n + 1):
        if i == n or np.isnan(xv[i]) or np.isnan(yv[i]):
            if i > start:
                ring = np.column_stack([xv[start:i], yv[start:i]])
                verts.append(ring)
                codes.append(np.full(len(ring), Path.LINETO))
                codes[-1][0] = Path.MOVETO
            start = i + 1
    if not verts:
        return np.zeros(xq.shape, dtype=bool)
    path = Path(np.vstack(verts), np.concatenate(codes))
    # radius=tiny so points exactly on the boundary count as inside (MATLAB)
    inside = path.contains_points(pts, radius=1e-9) | path.contains_points(pts, radius=-1e-9)
    return inside.reshape(xq.shape)


# ---------------------------------------------------------------------------
# interp1 / interp2 (MATLAB-compatible linear interpolation)
# ---------------------------------------------------------------------------
def interp1(x, v, xq, kind='linear', extrap=np.nan):
    """MATLAB ``interp1(x, v, xq, kind, extrapval)``."""
    x = np.asarray(x, dtype=float).ravel()
    v = np.asarray(v, dtype=float).ravel()
    order = np.argsort(x)
    x, v = x[order], v[order]
    f = interp1d(x, v, kind=kind, bounds_error=False, fill_value=extrap)
    return f(np.asarray(xq, dtype=float))


def interp2(X, Y, V, Xq, Yq, kind='linear', extrap=np.nan):
    """MATLAB ``interp2(X, Y, V, Xq, Yq, 'linear', extrapval)``.

    ``V`` has shape (len(Y), len(X)).  When ``Xq`` and ``Yq`` are vectors the
    result is the grid meshgrid(Xq, Yq) with shape (len(Yq), len(Xq)), matching
    MATLAB's behaviour for vector query inputs of differing orientation.
    """
    X = np.asarray(X, dtype=float).ravel()
    Y = np.asarray(Y, dtype=float).ravel()
    V = np.asarray(V, dtype=float)
    # RegularGridInterpolator needs strictly increasing axes
    xo = np.argsort(X); yo = np.argsort(Y)
    Xs, Ys = X[xo], Y[yo]
    Vs = V[np.ix_(yo, xo)]
    rgi = RegularGridInterpolator((Ys, Xs), Vs, method=kind,
                                  bounds_error=False, fill_value=extrap)
    Xq = np.asarray(Xq, dtype=float).ravel()
    Yq = np.asarray(Yq, dtype=float).ravel()
    XX, YY = np.meshgrid(Xq, Yq)        # -> shape (len(Yq), len(Xq))
    return rgi(np.column_stack([YY.ravel(), XX.ravel()])).reshape(XX.shape)


# ---------------------------------------------------------------------------
# Cartopy mapping helpers (approximate MATLAB Mapping Toolbox calls)
# ---------------------------------------------------------------------------
def worldmap(ax=None, latlim=(-90, 90), lonlim=(-180, 180), figsize=None):
    """Return a cartopy GeoAxes set up like MATLAB ``worldmap``."""
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
    if latlim == (-90, 90) and lonlim == (-180, 180):
        ax.set_global()                      # full globe (MATLAB worldmap default = Robinson oval)
    else:
        ax.set_extent([lonlim[0], lonlim[1], latlim[0], latlim[1]],
                      crs=ccrs.PlateCarree())
    return ax


def _to_pm180(lon, data_lonaxis0):
    """Shift a lon axis given in 0..360 to -180..180 and roll data to match.
    data is assumed (lat, lon)."""
    lon = np.asarray(lon, dtype=float)
    if lon.max() > 180:
        shift = lon > 180
        lon2 = np.where(shift, lon - 360, lon)
        order = np.argsort(lon2)
        return lon2[order], data_lonaxis0[:, order]
    return lon, data_lonaxis0


def pcolorm(ax, lat, lon, data, cmap=None, vmin=None, vmax=None, shift=True):
    """MATLAB ``pcolorm(lat, lon, data)`` -> pcolormesh on a GeoAxes.
    ``data`` has shape (len(lat), len(lon))."""
    import cartopy.crs as ccrs
    from matplotlib.colors import ListedColormap
    lat = np.asarray(lat, dtype=float).ravel()
    lon = np.asarray(lon, dtype=float).ravel()
    data = np.asarray(data, dtype=float)
    if shift:
        lon, data = _to_pm180(lon, data)
    if cmap is not None and not hasattr(cmap, 'N'):
        cmap = ListedColormap(cmap)
    return ax.pcolormesh(lon, lat, data, transform=ccrs.PlateCarree(),
                         cmap=cmap, vmin=vmin, vmax=vmax, shading='auto')


def contourfm(ax, lat, lon, data, levels, cmap=None, shift=True):
    """MATLAB ``contourfm(lat, lon, data, levels)`` -> filled contours."""
    import cartopy.crs as ccrs
    from matplotlib.colors import ListedColormap
    lat = np.asarray(lat, dtype=float).ravel()
    lon = np.asarray(lon, dtype=float).ravel()
    data = np.asarray(data, dtype=float)
    if shift:
        lon, data = _to_pm180(lon, data)
    if cmap is not None and not hasattr(cmap, 'N'):
        cmap = ListedColormap(cmap)
    return ax.contourf(lon, lat, data, levels=levels, transform=ccrs.PlateCarree(),
                       cmap=cmap, extend='both')


def geoshow_coastlines(ax, **kw):
    """MATLAB ``load coastlines; geoshow(coastlat, coastlon)`` -> coastlines."""
    kw.setdefault('color', 'k')
    kw.setdefault('linewidth', 0.8)
    ax.coastlines(**kw)
