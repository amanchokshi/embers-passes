import concurrent.futures
import json
from itertools import repeat
from pathlib import Path

import h5py
import healpy as hp
import numpy as np
from scipy.stats import beta

from embers.rf_tools.rf_data import tile_names
from embers.sat_utils.sat_channels import time_tree
from embers.tile_maps import tile_maps
from embers.tile_maps.beam_utils import chisq_fit_test, rotate_map


def get_Imatr(x, Nc=10):
    """
    Get a design matrix based on Bernstein polynomials where the first coefficient is a constant offset,
    and the next Nc coefficients parametrize the derivatives. If the coefficients are sorted and positive,
    the resulting polynomial will have a strictly increasing positive derivative (i.e. monotonic and concave up).

    Parameters:
        x (array):
            X coordinates of samples, RESCALED to live on [0, 1]. Do not need to be evenly spaced
        Nc (int):
            Degree of polynomial (e.g. Nc = 3 makes a cubic).

    Returns:
        Imatr (array):
            Design matrix.
    """
    Nx = len(x)
    Imatr = np.zeros([Nx, Nc + 1])
    Imatr[:, 0] = np.ones(Nx)
    for k in range(1, Nc + 1):
        Imatr[:, k] = beta(k, Nc + 1 - k).cdf(x)

    return Imatr


def scale_to_unit(
    x: np.ndarray,
    fit_min: float,
    fit_max: float,
) -> np.ndarray:
    """
    Linearly scale data from [fit_min, fit_max] to [0, 1].

    Parameters
    ----------
    x
        Input data.
    fit_min
        Lower bound of the fitting interval.
    fit_max
        Upper bound of the fitting interval.

    Returns
    -------
    ndarray
        Data scaled to the unit interval.
    """
    x = np.asarray(x, dtype=float)
    return (x - fit_min) / (fit_max - fit_min)


def scale_from_unit(
    u: np.ndarray,
    fit_min: float,
    fit_max: float,
) -> np.ndarray:
    """
    Convert data from the unit interval [0, 1] back to the original scale.

    Parameters
    ----------
    u
        Data in the unit interval.
    fit_min
        Lower bound of the fitting interval.
    fit_max
        Upper bound of the fitting interval.

    Returns
    -------
    ndarray
        Data on the original scale.
    """
    u = np.asarray(u, dtype=float)
    return fit_min + u * (fit_max - fit_min)


def load_rfe_fit(
    fit_file: str | Path,
    fit_key: str,
) -> dict:
    """Load one tile/reference compression fit."""
    with Path(fit_file).open("r") as stream:
        payload = json.load(stream)

    try:
        return payload["fits"][fit_key]
    except KeyError as exc:
        available = ", ".join(sorted(payload.get("fits", {})))
        raise KeyError(
            f"No fit named {fit_key!r}. Available fits: {available}"
        ) from exc


def evaluate_rfe_residual(
    observed_power: np.ndarray,
    fit: dict,
    *,
    below_fit: str = "zero",
    above_fit: str = "clip",
) -> np.ndarray:
    """
    Evaluate the missing-power correction for observed RF Explorer power.

    Parameters
    ----------
    observed_power
        Measured tile power in dBm.
    fit
        Record returned by `load_rfe_fit`.
    below_fit
        Behaviour below `fit_min_dbm`:

        - ``"zero"``: return zero correction;
        - ``"extrapolate"``: evaluate outside the fitted domain.

    above_fit
        Behaviour above `fit_max_dbm`:

        - ``"clip"``: use the correction at the upper boundary;
        - ``"extrapolate"``: evaluate outside the fitted domain;
        - ``"nan"``: return NaN.

    Returns
    -------
    ndarray
        Additive correction in dB.
    """
    power = np.asarray(observed_power, dtype=float)

    fit_min = float(fit["fit_min_dbm"])
    fit_max = float(fit["fit_max_dbm"])
    degree = int(fit["degree"])

    coeffs = np.asarray(
        fit["derivative_coefficients"],
        dtype=float,
    )

    if coeffs.size != degree:
        raise ValueError(f"Expected {degree} coefficients, found {coeffs.size}.")

    u = scale_to_unit(power, fit_min, fit_max)
    u_eval = u.copy()

    if above_fit == "clip":
        u_eval = np.minimum(u_eval, 1.0)
    elif above_fit == "nan":
        pass
    elif above_fit != "extrapolate":
        raise ValueError(f"Unknown above_fit option: {above_fit!r}")

    if below_fit == "zero":
        u_eval = np.maximum(u_eval, 0.0)
    elif below_fit != "extrapolate":
        raise ValueError(f"Unknown below_fit option: {below_fit!r}")

    design = get_Imatr(u_eval, Nc=degree)
    correction = design[:, 1:] @ coeffs

    if below_fit == "zero":
        correction = np.where(power < fit_min, 0.0, correction)

    if above_fit == "nan":
        correction = np.where(power > fit_max, np.nan, correction)

    correction = np.where(np.isfinite(power), correction, np.nan)

    return correction


def fit_offset(data, model):
    """Fit an additive offset between data and model in dB."""
    data = np.asarray(data, dtype=float)
    model = np.asarray(model, dtype=float)

    mask = np.isfinite(data) & np.isfinite(model)
    if not np.any(mask):
        return np.nan

    return np.median(data[mask] - model[mask])


def project_tile_passes(
    start_date,
    stop_date,
    tile_pair,
    sat_thresh,
    noi_thresh,
    pow_thresh,
    ref_model,
    fee_map,
    rfe_cali,
    obs_point_json,
    align_dir,
    chrono_dir,
    chan_map_dir,
    out_dir,
    rfe_cali_bool=True,
):
    """
    Extract calibrated satellite passes and save to HDF5.

    For each satellite pass:
    - Apply thresholding
    - Interpolate beam models to native alt/az
    - Construct calibrated pass
    - Fit gain offset
    - Keep passes with pval >= 0.8

    Output:
        HDF5 file with passes grouped by pointing and sorted by mean unix time.
    """

    good_sats = {
        25338,
        25982,
        25984,
        25985,
        28654,
        40086,
        40087,
        40091,
        41179,
        41180,
        41182,
        41183,
        41184,
        41185,
        41187,
        41188,
        41189,
        44387,
    }

    valid_pointings = {0, 2, 4, 41}

    ref, tile = tile_pair

    rfe_fit = None

    if rfe_cali_bool:
        fit_key = f"{tile}_{ref}"

        rfe_fit = load_rfe_fit(
            fit_file=rfe_cali,
            fit_key=fit_key,
        )

    dates, timestamps = time_tree(start_date, stop_date)

    out_path = Path(out_dir) / f"{tile}_{ref}_passes.h5"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    ref_fee_model = np.load(ref_model, allow_pickle=True)
    fee_m = np.load(fee_map, allow_pickle=True)

    if "XX" in tile:
        ref_fee = ref_fee_model["XX"]
    else:
        ref_fee = ref_fee_model["YY"]

    rotated_fee = rotate_map(
        hp.get_nside(ref_fee),
        angle=-(np.pi / 2),
        healpix_array=ref_fee,
    )

    accepted_passes = {p: [] for p in sorted(valid_pointings)}

    for day in range(len(dates)):
        for window in range(len(timestamps[day])):
            timestamp = timestamps[day][window]

            point = tile_maps.check_pointing(timestamp, obs_point_json)
            if point not in valid_pointings:
                continue

            if "XX" in tile:
                mwa_fee = fee_m[str(point)][0]
            else:
                mwa_fee = fee_m[str(point)][1]

            ali_file = Path(
                f"{align_dir}/{dates[day]}/{timestamp}/{ref}_{tile}_{timestamp}_aligned.npz"
            )

            if not ali_file.is_file():
                continue

            chrono_file = Path(f"{chrono_dir}/{timestamp}.json")
            channel_map_file = Path(f"{chan_map_dir}/{timestamp}.json")

            if not channel_map_file.is_file():
                continue

            with open(chrono_file) as f:
                chrono_ephem = json.load(f)

            if not chrono_ephem:
                continue

            with open(channel_map_file) as f:
                chan_map = json.load(f)

            for sat in map(int, chan_map.keys()):
                if sat not in good_sats:
                    continue

                chan = chan_map[str(sat)]

                sat_data = tile_maps.rf_apply_thresholds(
                    ali_file,
                    chrono_file,
                    sat,
                    chan,
                    sat_thresh,
                    noi_thresh,
                    pow_thresh,
                    point,
                    plots=False,
                    out_dir=f"{out_dir}/tile_maps_raw",
                )

                if sat_data == 0:
                    continue

                (
                    ref_power,
                    tile_power,
                    alt_deg,
                    az_rad,
                    times,
                    mwa_sigma_db,
                ) = sat_data

                ref_pass = np.asarray(ref_power)
                tile_pass = np.asarray(tile_power)
                alt_deg = np.asarray(alt_deg)
                az_rad = np.asarray(az_rad)
                times_pass = np.asarray(times)

                alt_rad = np.radians(alt_deg)
                az_deg = np.degrees(az_rad)
                za = np.pi / 2 - alt_rad

                ref_fee_pass = hp.get_interp_val(rotated_fee, za, az_rad)
                mwa_fee_pass = hp.get_interp_val(mwa_fee, za, az_rad)

                if rfe_cali_bool:
                    correction_db = evaluate_rfe_residual(
                        observed_power=tile_pass,
                        fit=rfe_fit,
                        below_fit="zero",
                        above_fit="clip",
                    )

                    tile_pass_rfe = tile_pass + correction_db
                else:
                    correction_db = np.zeros_like(tile_pass, dtype=float)
                    tile_pass_rfe = tile_pass.copy()

                mwa_pass = tile_pass_rfe - ref_pass + ref_fee_pass

                # offset = chisq_fit_gain(data=mwa_pass, model=mwa_fee_pass)
                # mwa_pass_fit = mwa_pass - offset[0]
                offset = fit_offset(data=mwa_pass, model=mwa_fee_pass)
                mwa_pass_fit = mwa_pass - offset

                if mwa_pass_fit.size == 0 or not np.isfinite(mwa_pass_fit).all():
                    continue

                pval = chisq_fit_test(
                    data=mwa_pass_fit,
                    model=mwa_fee_pass,
                )

                if not np.isfinite(pval) or pval < 0.8:
                    continue

                accepted_passes[point].append(
                    {
                        "norad_id": int(sat),
                        "unix_time": float(np.mean(times_pass)),
                        "alt_deg": alt_deg,
                        "az_deg": az_deg,
                        "power_db": mwa_pass_fit,
                        "rfe_correction_db": correction_db,
                        "noise_db": float(mwa_sigma_db),
                    }
                )

    for point in accepted_passes:
        accepted_passes[point].sort(key=lambda x: x["unix_time"])

    with h5py.File(out_path, "w") as h5:
        h5.attrs["tile"] = tile
        h5.attrs["ref"] = ref
        h5.attrs["start_date"] = start_date
        h5.attrs["stop_date"] = stop_date
        h5.attrs["pointings"] = np.array(sorted(valid_pointings), dtype=int)
        h5.attrs["n_passes_total"] = int(
            sum(len(accepted_passes[p]) for p in accepted_passes)
        )

        pointings_grp = h5.create_group("pointings")

        for point in sorted(valid_pointings):
            point_grp = pointings_grp.create_group(str(point))
            point_grp.attrs["pointing"] = int(point)
            point_grp.attrs["n_passes"] = len(accepted_passes[point])

            passes_grp = point_grp.create_group("passes")

            for i, p in enumerate(accepted_passes[point]):
                grp = passes_grp.create_group(f"{i:06d}")

                grp.attrs["norad_id"] = p["norad_id"]
                grp.attrs["unix_time"] = p["unix_time"]
                grp.attrs["noise_db"] = p["noise_db"]

                grp.create_dataset("alt_deg", data=p["alt_deg"], compression="gzip")
                grp.create_dataset("az_deg", data=p["az_deg"], compression="gzip")
                grp.create_dataset("power_db", data=p["power_db"], compression="gzip")
                grp.create_dataset(
                    "rfe_correction_db", data=p["rfe_correction_db"], compression="gzip"
                )


def tile_passes_batch(
    start_date,
    stop_date,
    sat_thresh,
    noi_thresh,
    pow_thresh,
    ref_model,
    fee_map,
    rfe_cali,
    obs_point_json,
    align_dir,
    chrono_dir,
    chan_map_dir,
    out_dir,
    rfe_cali_bool=True,
    max_cores=None,
):
    """Batch process satellite RF data to create calibrated pass files.

    For each valid reference/tile pair, extract calibrated satellite passes and
    save them to an HDF5 file using :func:`project_tile_passes`.

    Parameters
    ----------
    start_date : str
        Start date in ``YYYY-MM-DD-HH:MM`` format.
    stop_date : str
        Stop date in ``YYYY-MM-DD-HH:MM`` format.
    sat_thresh : float
        Sigma threshold used in noise-floor estimation.
    noi_thresh : float
        Noise threshold multiplier.
    pow_thresh : float
        Minimum peak-above-noise threshold required for a pass.
    ref_model : str or Path
        Path to the reference beam model ``.npz`` file.
    fee_map : str or Path
        Path to the MWA FEE model ``.npz`` file.
    rfe_cali : str or Path
        Path to the RFE calibration solution.
    obs_point_json : str or Path
        Path to observation pointing metadata JSON.
    align_dir : str or Path
        Directory containing aligned RF data.
    chrono_dir : str or Path
        Directory containing chronological ephemeris JSON files.
    chan_map_dir : str or Path
        Directory containing satellite channel map JSON files.
    out_dir : str or Path
        Output directory where pass HDF5 files will be written.
    rfe_cali_bool : bool, optional
        If True, apply RFE calibration. Default is True.
    max_cores : int or None, optional
        Maximum number of worker processes. Default is None, meaning use the
        executor default.

    Returns
    -------
    None
    """

    refs = tile_names()[:4]
    tiles = tile_names()[4:]

    tile_pairs = []
    for ref in refs:
        if "XX" in ref:
            for tile in (t for t in tiles if "XX" in t):
                tile_pairs.append([ref, tile])
        else:
            for tile in (t for t in tiles if "YY" in t):
                tile_pairs.append([ref, tile])

    Path(out_dir).mkdir(parents=True, exist_ok=True)

    with concurrent.futures.ProcessPoolExecutor(max_workers=max_cores) as executor:
        list(
            executor.map(
                project_tile_passes,
                repeat(start_date),
                repeat(stop_date),
                tile_pairs,
                repeat(sat_thresh),
                repeat(noi_thresh),
                repeat(pow_thresh),
                repeat(ref_model),
                repeat(fee_map),
                repeat(rfe_cali),
                repeat(obs_point_json),
                repeat(align_dir),
                repeat(chrono_dir),
                repeat(chan_map_dir),
                repeat(out_dir),
                repeat(rfe_cali_bool),
            )
        )


if __name__ == "__main__":
    tile_passes_batch(
        start_date="2019-09-10",
        stop_date="2020-03-18",
        sat_thresh=1,
        noi_thresh=3,
        pow_thresh=5,
        ref_model="./embers_out/tile_maps/ref_models/ref_dipole_models.npz",
        fee_map="./embers_out/mwa_utils/mwa_fee/mwa_fee_beam.npz",
        rfe_cali="./embers_out/tile_maps/rfe_calibration/20260714/rfe_compression_monotonic_concave_fits.json",
        obs_point_json="./embers_out/mwa_utils/obs_pointings.json",
        align_dir="./embers_out/rf_tools/align_data",
        chrono_dir="./embers_out/sat_utils/ephem_chrono",
        chan_map_dir="./embers_out//sat_utils/sat_channels/window_maps",
        out_dir="./embers_out/tile_maps/passes_concave_monotonic",
        rfe_cali_bool=True,
        max_cores=None,
    )
