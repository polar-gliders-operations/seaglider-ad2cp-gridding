"""Helper functions to read Seaglider AD2CP data from many dive files and build 2D arrays.

- load_adcp_files + grid_adcp: raw ADCP velocities (velX, velY) for each of the 15 cells,
  binned into (pressure, time) arrays.
- load_inv_profiles: processed ocean velocities (uocn, vocn) as (depth, time) arrays.

"time" columns are half-dives: every dive gives one descent column and one ascent column.
"""

import glob
import os

import numpy as np
import xarray as xr

# The only variables we need from each file (reading just these is much faster)
ADCP_VARS = ["ad2cp_time", "ad2cp_pressure", "ad2cp_velX", "ad2cp_velY"]


def load_adcp_files(path_to_files, pattern="*.nc"):
    """Read all dive files and concatenate the AD2CP data along ad2cp_data_data_point.

    Adds `dive` (dive number from the filename) and `cast` (0 = descent, 1 = ascent) per data point.
    Files without AD2CP data are skipped.
    """
    # List all dive files in the folder, sorted so dives are in order
    files = sorted(glob.glob(os.path.join(path_to_files, pattern)))
    pieces, skipped = [], []
    for f in files:
        # Seaglider names: pGGGDDDD.nc -> glider GGG, dive DDDD
        dive_idx = int(os.path.basename(f)[4:8])
        # decode_times=False: some files contain unrelated, undecodable time variables,
        # so we only convert times for the variables we actually keep (decode_cf below)
        with xr.open_dataset(f, decode_times=False) as ds:
            # Some dives have no ADCP data at all -> skip them
            if not all(v in ds for v in ADCP_VARS) or ds.sizes.get("ad2cp_data_data_point", 0) == 0:
                skipped.append(os.path.basename(f))
                continue
            sub = xr.decode_cf(ds[ADCP_VARS]).load()

        # Split the dive at its deepest point: pings before = descent (0), after = ascent (1)
        p = sub["ad2cp_pressure"].values
        i_max = np.nanargmax(p) if np.isfinite(p).any() else len(p)
        cast = (np.arange(len(p)) > i_max).astype(int)
        # Label every ping with its dive number and descent/ascent flag
        sub = sub.assign(
            dive=("ad2cp_data_data_point", np.full(len(p), dive_idx)),
            cast=("ad2cp_data_data_point", cast),
        )
        pieces.append(sub)

    print(f"loaded {len(pieces)} files, skipped {len(skipped)} without AD2CP data")
    # Stack all dives into one long dataset of pings
    return xr.concat(pieces, dim="ad2cp_data_data_point")


def grid_adcp(ds, dp=5.0, p_max=None):
    """Bin each profile (dive + cast) onto a regular pressure grid.

    Returns a Dataset with ad2cp_velX / ad2cp_velY of shape
    (ad2cp_cell_data_point, pressure, time); `time` is the mean time of each profile.
    """
    # --- Pressure axis: bins of width dp from the surface to the deepest ping ---
    p = ds["ad2cp_pressure"].values
    if p_max is None:
        p_max = np.nanmax(p)
    edges = np.arange(0, p_max + dp, dp)          # bin edges, e.g. 0, 5, 10, ...
    centres = 0.5 * (edges[:-1] + edges[1:])      # bin centres, e.g. 2.5, 7.5, ...
    n_p = len(centres)

    # --- Time axis: one column per half-dive ---
    # Unique number per half-dive: dive 5 descent -> 10, dive 5 ascent -> 11
    profile_id = ds["dive"].values * 2 + ds["cast"].values
    profiles = np.unique(profile_id)
    n_t = len(profiles)
    n_cells = ds.sizes["ad2cp_cell_data_point"]

    # For every ping: which pressure row (p_bin) and which time column (t_idx) it belongs to
    p_bin = np.digitize(p, edges) - 1
    valid_p = np.isfinite(p) & (p_bin >= 0) & (p_bin < n_p)
    t_idx = np.searchsorted(profiles, profile_id)

    # Time of each column = average time of all pings in that half-dive
    times = np.array(
        [ds["ad2cp_time"].values[profile_id == pid].astype("int64").mean() for pid in profiles]
    ).astype("datetime64[ns]")

    # --- Fill the grids: average all pings that fall into the same (pressure, time) box ---
    out = {}
    for name in ["ad2cp_velX", "ad2cp_velY"]:
        v = ds[name].values  # (n_points, n_cells)
        grid = np.full((n_cells, n_p, n_t), np.nan)  # NaN where there is no data
        for c in range(n_cells):
            ok = valid_p & np.isfinite(v[:, c])
            s = np.zeros((n_p, n_t))  # sum of velocities in each box
            n = np.zeros((n_p, n_t))  # number of pings in each box
            np.add.at(s, (p_bin[ok], t_idx[ok]), v[ok, c])
            np.add.at(n, (p_bin[ok], t_idx[ok]), 1)
            with np.errstate(invalid="ignore"):
                grid[c] = s / n  # mean = sum / count (0/0 -> NaN for empty boxes)
        out[name] = (("ad2cp_cell_data_point", "pressure", "time"), grid, ds[name].attrs)

    return xr.Dataset(
        out,
        coords={
            "ad2cp_cell_data_point": np.arange(n_cells),
            "pressure": ("pressure", centres, {"units": "dbar"}),
            "time": times,
            "dive": ("time", profiles // 2),
            "cast": ("time", profiles % 2, {"description": "0 = descent, 1 = ascent"}),
        },
    )


INV_VARS = ["ad2cp_inv_profile_depth", "ad2cp_inv_profile_uocn", "ad2cp_inv_profile_vocn",
            "ad2cp_inv_glider_time", "ad2cp_inv_glider_depth"]


def load_inv_profiles(path_to_files, pattern="*.nc"):
    """Read the inverse-solution ocean velocity profiles (uocn, vocn) from all dive files.

    Each file holds two profiles on a common depth grid: index 0 = descent, 1 = ascent.
    Returns a Dataset with ad2cp_inv_profile_uocn / ad2cp_inv_profile_vocn of shape
    (depth, time), one time column per descent/ascent. `time` is the mean
    ad2cp_inv_glider_time of that half-dive (split at the deepest glider point).
    """
    files = sorted(glob.glob(os.path.join(path_to_files, pattern)))
    depth = None
    cols_u, cols_v, times, dives, casts = [], [], [], [], []
    skipped = 0
    for f in files:
        dive = int(os.path.basename(f)[4:8])  # dive number from filename pGGGDDDD.nc
        with xr.open_dataset(f, decode_times=False) as ds:
            # Later dives have no inverse solution -> skip them
            if not all(v in ds for v in INV_VARS) or ds.sizes.get("ad2cp_inv_profile", 0) != 2:
                skipped += 1
                continue
            sub = xr.decode_cf(ds[INV_VARS]).load()

        # Depth grid of the profiles (0-1000 m every 5 m); the first file defines the common grid
        z = sub["ad2cp_inv_profile_depth"].values
        if depth is None:
            depth = z
        # Shape (depth, 2): column 0 = descent profile, column 1 = ascent profile
        u = sub["ad2cp_inv_profile_uocn"].transpose("ad2cp_inv_profile_point", "ad2cp_inv_profile").values
        v = sub["ad2cp_inv_profile_vocn"].transpose("ad2cp_inv_profile_point", "ad2cp_inv_profile").values
        if not np.array_equal(z, depth):
            # put profiles onto the common depth grid if a file uses a different one
            u = np.column_stack([np.interp(depth, z, u[:, k], left=np.nan, right=np.nan) for k in range(2)])
            v = np.column_stack([np.interp(depth, z, v[:, k], left=np.nan, right=np.nan) for k in range(2)])

        # Give each profile a time: split the glider track at its deepest point and
        # take the mean time of the descent half and of the ascent half
        t = sub["ad2cp_inv_glider_time"].values.astype("int64")
        gd = sub["ad2cp_inv_glider_depth"].values
        i_max = np.nanargmax(gd) if np.isfinite(gd).any() else len(gd) // 2
        halves = [t[: i_max + 1], t[i_max + 1:] if i_max + 1 < len(t) else t[i_max:]]

        # Store descent then ascent as two new time columns
        for k in range(2):
            cols_u.append(u[:, k])
            cols_v.append(v[:, k])
            times.append(int(halves[k].mean()))
            dives.append(dive)
            casts.append(k)

    print(f"loaded {len(cols_u) // 2} files, skipped {skipped} without inverse profiles")
    # Glue all columns side by side into (depth, time) arrays
    attrs_u = {"units": "m s-1", "description": "Eastward ocean velocity (inverse solution)"}
    attrs_v = {"units": "m s-1", "description": "Northward ocean velocity (inverse solution)"}
    return xr.Dataset(
        {
            "ad2cp_inv_profile_uocn": (("depth", "time"), np.column_stack(cols_u), attrs_u),
            "ad2cp_inv_profile_vocn": (("depth", "time"), np.column_stack(cols_v), attrs_v),
        },
        coords={
            "depth": ("depth", depth, {"units": "m"}),
            "time": np.array(times).astype("datetime64[ns]"),
            "dive": ("time", np.array(dives)),
            "cast": ("time", np.array(casts), {"description": "0 = descent, 1 = ascent"}),
        },
    )
