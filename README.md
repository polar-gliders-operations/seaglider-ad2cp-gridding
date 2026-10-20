# seaglider-ad2cp-gridding

Tools for loading and gridding velocity data from the **Nortek Signature 1000 (1 MHz) AD2CP** mounted on a Seaglider. The code reads the per-dive NetCDF files from the Seaglider basestation (`pGGGDDDD.nc`) and builds 2D (depth × time) arrays, so you can look at how currents evolve over a whole mission.

Developed by the Polar Gliders group at the University of Gothenburg. The first application is Seaglider **SG267** on a Cape Basin transect.

## What it does

The repository handles two AD2CP products found in each dive file:

| Product | Variables | Function | Output shape |
|---|---|---|---|
| **Raw ADCP velocities** (relative to the glider, 15 cells) | `ad2cp_velX`, `ad2cp_velY` | `load_adcp_files` + `grid_adcp` | `(cell, pressure, time)` |
| **Ocean velocities** from the inverse solution (eastward/northward, 5 m grid, 0–1000 m) | `ad2cp_inv_profile_uocn`, `ad2cp_inv_profile_vocn` | `load_inv_profiles` | `(depth, time)` |

In both cases the `time` axis has **one column per half-dive**: each dive gives a descent column and an ascent column. Every column carries `dive` and `cast` coordinates (`cast` = 0 for descent, 1 for ascent), so you can always trace it back to the source file.

## Repository layout

```
code/
  adcp_gridding.py           # helper functions for reading and gridding
  explore_adcp_data.ipynb    # example notebook: load SG267 ocean velocities and plot them
documentation/
  Seaglider_ADCP_processing.pdf   # description of the APL-UW AD2CP processing / inverse method
```

## Requirements

- Python 3.10+
- `numpy`
- `xarray`
- `netCDF4` (NetCDF backend for xarray)
- `matplotlib` (for the notebook)
- `jupyter` (to run the notebook)

For example, with conda:

```bash
conda create -n sg-adcp python numpy xarray netcdf4 matplotlib jupyter
conda activate sg-adcp
```

## Getting started

```bash
git clone git@github.com:polar-gliders-operations/seaglider-ad2cp-gridding.git
cd seaglider-ad2cp-gridding/code
jupyter lab explore_adcp_data.ipynb
```

The notebook expects the dive files in `../data/sg267/` relative to `code/`, i.e.

```
seaglider-ad2cp-gridding/
  code/
  data/
    sg267/
      p2670001.nc
      p2670002.nc
      ...
```

Data files are **not** included in the repository. Change `path_to_files` in the notebook if your data live elsewhere.

## Usage

### Ocean velocities (inverse solution)

```python
from adcp_gridding import load_inv_profiles

ds_ocn = load_inv_profiles("../data/sg267/")
uocn = ds_ocn["ad2cp_inv_profile_uocn"]   # eastward velocity, (depth, time), m/s
vocn = ds_ocn["ad2cp_inv_profile_vocn"]   # northward velocity, (depth, time), m/s
```

These profiles are already on a common depth grid in each file, so no binning is needed; they are placed side by side in time order. The time of each column is the mean glider time of that half-dive, with the dive split at its deepest point. Files without inverse profiles are skipped (for SG267, only dives 1–244 contain them).

### Raw ADCP velocities, binned in pressure

```python
from adcp_gridding import load_adcp_files, grid_adcp

ds_raw = load_adcp_files("../data/sg267/")
grid = grid_adcp(ds_raw, dp=5.0)          # 5 dbar bins from the surface to the deepest ping
velX = grid["ad2cp_velX"]                 # (cell, pressure, time)
```

`load_adcp_files` concatenates all pings and labels each one with its dive number and cast. `grid_adcp` then averages all pings that fall into the same (pressure bin, half-dive) box, separately for each of the 15 cells. Empty boxes are NaN.

### Caching

Reading several hundred dive files (especially from a network or Google Drive mount) can take minutes. The notebook saves the result to `../data/sg267_adcp_inv_profiles.nc` the first time and re-opens it afterwards. Delete that file to force a rebuild when the raw data change.

## Notes and assumptions

- Dive numbers are taken from the filename (`pGGGDDDD.nc` → characters 5–8). Files that don't follow this naming will break the loaders.
- Files are opened with `decode_times=False` because some contain unrelated, undecodable time variables. Only the variables that are kept are decoded.
- Descent and ascent are separated at the point of maximum pressure (raw data) or maximum glider depth (inverse profiles).
- If a file's inverse-profile depth grid differs from the first file's, its profiles are linearly interpolated onto the common grid.

## Background

The inverse solution for ocean velocity is produced by the APL-UW Seaglider AD2CP processing, which follows Todd et al. (2017) and the LADCP approach of Visbeck (2002). See `documentation/Seaglider_ADCP_processing.pdf` (Rainville, Shilling, Gobat & Lee, APL-UW, 2024) and the [iop-apl-uw/ad2cp_python](https://github.com/iop-apl-uw/ad2cp_python) repository.

- Todd, R. E., Rudnick, D. L., Sherman, J. T., Owens, W. B., & George, L. (2017). Absolute velocity estimates from autonomous underwater gliders equipped with Doppler current profilers. *J. Atmos. Oceanic Technol.*, 34, 309–333.
- Visbeck, M. (2002). Deep velocity profiling using lowered acoustic Doppler current profilers: Bottom track and inverse solutions. *J. Atmos. Oceanic Technol.*, 19, 794–807.

## Contact

Marcel du Plessis — Polar Gliders, University of Gothenburg
