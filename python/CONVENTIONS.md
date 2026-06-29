# MATLAB → Python conversion conventions

Goal: **faithful** translation of the *Uncertainties* MATLAB scripts into
Jupyter notebooks. The analysis (every number/array) must be identical; figures
should be close (matplotlib), not pixel-perfect.

## Mechanics
- Author each notebook as a **percent-format `.py`** in `python/src/<name>.py`
  (same `<name>` as the `.m` file). Use `# %%` for code cells and
  `# %% [markdown]` for explanation cells. Copy the header block (jupytext +
  `kernelspec` `geomip`) from `python/src/plot_HS_GeoMIP.py`.
- Build & execute with: `bash python/build.sh <name>` (run from project root).
  It converts to `<name>.ipynb` in the project root and executes it with the
  `geomip` kernel. **The notebook must run with no errors and embed the
  figures.** Iterate until clean.
- First code cell of every notebook:
  ```python
  import sys; sys.path.insert(0, "python")
  import numpy as np
  import matplotlib.pyplot as plt
  import mat_utils as mu
  from mat_utils import brewermap, brighten, movmean, ncread, shaperead, inpolygon, interp1, interp2
  ```
- Notebooks execute with **cwd = project root**, so use the *same relative data
  paths as the MATLAB script* (`'GeoMIP_tas.mat'`, `'../GeoMIP_G6Sulfur/...'`,
  `'/Volumes/D_Visioni2/...'`, `'AREA.nc'`, shapefiles, etc.).

## Indexing & numerics (critical for fidelity)
- MATLAB is **1-based**, Python **0-based**. `x(6:86)` → `x[5:86]`;
  `(:,71:81)` → `[:, 70:81]`; `x(:,:,1:81)` → `x[:, :, 0:81]`.
- MATLAB column-major. `.mat` numeric arrays keep their MATLAB shape via
  `mu.loadmat`. Cell arrays load as **object ndarrays**: a `1×N` cell `LAT` →
  `LAT[0, k]`; an `M×N×P` cell `x_e{im,ie,in}` → `x_e[im-1, ie-1, in-1]`.
- `std(X)` and `std(X,0)` use **N-1** → `np.std(x, ddof=1)`. `std(X,1)` →
  `np.std(x, ddof=0)`. `var` likewise. **Default MATLAB std is ddof=1.**
- `mean`/`sum` over a MATLAB dim `d` → numpy `axis=d-1`. `squeeze` → `np.squeeze`.
- `nanmean`/`nansum`/`mean(...,'omitnan')` → `np.nanmean`/`np.nansum`.
- `movmean(x,k)` → `mu.movmean(x,k)` (handles MATLAB edge/even-window rules and
  first-non-singleton-dim). `movmean(...,k)` with `'omitnan'` → add `omitnan=True`.
- `reshape`: MATLAB reshape is column-major (Fortran order). Use
  `np.reshape(a, shape, order='F')`. e.g. `reshape(x,12,N,...)` then
  `mean(...,1)` (mean over the 12 months) → `np.reshape(x,(12,...),order='F').mean(axis=0)`.
- `cos(lat/180*pi)` weights: `ww = np.cos(lat/180*np.pi)`.

## Grids, polygons, NetCDF
- `[A,B]=meshgrid(u,v)` (MATLAB) → `A,B = np.meshgrid(u, v)` (numpy default
  `indexing='xy'` matches: output shape `(len(v), len(u))`, `A[i,j]=u[j]`,
  `B[i,j]=v[i]`).
- `AA = meshgrid(ww, lon)` (single output) → `AA = np.meshgrid(ww, lon)[0]`
  (shape `(nlon, nlat)`, `AA[i,j]=ww[j]`).
- `regions = shaperead('...shx')` → `regions = mu.shaperead('...shx')`; access
  `regions[i].Name` / `.Region` / `.X` (lon, NaN-separated) / `.Y` (lat).
  `size(regions,1)` → `len(regions)`.
- `inpolygon(latM,lonM,regions(i).Y,regions(i).X)` →
  `mu.inpolygon(latM, lonM, regions[i].Y, regions[i].X)` (returns boolean array
  same shape as `latM`). `field(ind)` → `field[ind]` (1-D); sums:
  `np.nansum(field[ind])`.
- `ncread(file,var)` → `mu.ncread(file, var)` (already gives MATLAB lon×lat×time
  order with fill→NaN). `ncdisp` → drop or replace with a comment.
- `interp1(x,v,xq,'linear',NaN)` → `mu.interp1(x,v,xq,'linear',np.nan)`.
  `interp2(X,Y,V,Xq,Yq,'linear',NaN)` → `mu.interp2(X,Y,V,Xq,Yq,'linear',np.nan)`
  (`mu.interp2` returns the meshgrid result shape `(len(Yq), len(Xq))`).
- `a:step:b` → `np.arange(a, b+step, step)` (include endpoint when b is hit).

## Colours & plotting
- `brewermap(N,'Scheme')` → `brewermap(N,'Scheme')` (returns N×3 in [0,1]).
  `'*Scheme'` reversal supported. `brighten(map,beta)` → `brighten(map,beta)`.
- A colormap row `fc(k,:)` → `fc[k-1]`. Pass as `color=fc[k]`.
- `figure; subplot(m,n,p)` → `fig, axes = plt.subplots(...)` or
  `fig.add_subplot(m,n,p)`. `hold on`/`box on` are default in matplotlib.
- `plot(t,y,'Color',fc(1,:),'Linewidth',2)` →
  `ax.plot(t, y, color=fc[0], linewidth=2)`.
- `area(...)` → `ax.fill_between(x, y, 0, ...)` or `ax.stackplot`/overlapping
  `fill_between` to mimic the layered MATLAB `area` calls.
- `errorbar(x,y,e,...)` → `ax.errorbar(x, y, yerr=e, ...)`.
- `set(gca,'XTick',...)`, `axis([..])`, `title`, `xlabel`, `ylabel`, `legend`
  → matplotlib equivalents. `movmean`-smoothed lines exactly as MATLAB.
- Maps: `worldmap([-90 90],[-180 180])` → `ax = mu.worldmap(figsize=...)` (pass
  an existing GeoAxes via `mu.worldmap(ax=ax)` when using subplots created with
  `subplot_kw={'projection': ccrs.PlateCarree()}`). `pcolorm(lat,lon,Z)` →
  `mu.pcolorm(ax, lat, lon, Z, cmap=fc, vmin=mmin, vmax=mmax)`.
  `contourfm(lat,lon,Z,v)` → `mu.contourfm(ax, lat, lon, Z, v, cmap=fc)`.
  `load coastlines; geoshow(coastlat,coastlon)` → `mu.geoshow_coastlines(ax)`.
  `caxis([a b])` → pass `vmin=a, vmax=b` to pcolorm or `mappable.set_clim(a,b)`.
  `colorbar(...)` → `fig.colorbar(mappable, ax=ax, ...)`.

## Saving outputs
- `save('foo.mat', 'a','b')` → `from scipy.io import savemat`. **Do NOT overwrite
  existing `.mat` inputs.** For data-prep scripts that regenerate an existing
  file, save to a `py_outputs/` subfolder (create it) with the *same base name*,
  e.g. `savemat('py_outputs/GeoMIP_pr.mat', {...})`, and note this in markdown.
- `print(gcf,'-depsc2','name.eps')` → `fig.savefig('py_outputs/name.eps')`
  (also create `py_outputs/`; don't clobber the existing `.eps`).

## When the MATLAB script is not standalone
Some scripts reuse workspace variables left by another script, or reference
variables that are undefined on their own. Per the agreed plan, **merge** such
pairs into one notebook (loader first, then the dependent section). If a
reference is genuinely undefined even in MATLAB (a latent bug), reproduce the
*intended* computation by loading the needed data and add a short markdown note
explaining the deviation. Each notebook must stand alone and execute.
