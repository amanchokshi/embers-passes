from __future__ import annotations

import argparse
import json
import os
import pickle
from pathlib import Path
from typing import Any

import healpy as hp
import mwa_hyperbeam
import numpy as np
import pymc as pm
import pytensor.tensor as pt
from pytensor.graph.op import Op
from scipy.stats import median_abs_deviation as mad

from embers_passes import PassFile


def pass_hyperbeam_model(
    p,
    beam,
    freq_hz: float = 137e6,
    delays: list[int] | None = None,
    amps: list[float] | None = None,
    norm_to_zenith: bool = True,
    parallactic: bool = False,
) -> np.ndarray:
    """
    Evaluate the MWA hyperbeam model along a single pass.

    Parameters
    ----------
    p : PassRecord
        Pass record with alt_deg and az_deg arrays.
    beam : mwa_hyperbeam.FEEBeam
        Hyperbeam beam object.
    freq_hz : float, optional
        Frequency in Hz.
    delays : list[int], optional
        Beamformer delays.
    amps : list[float], optional
        Dipole amplitudes.
    norm_to_zenith : bool, optional
        Whether to normalize to zenith.
    parallactic : bool, optional
        Whether to include parallactic rotation.

    Returns
    -------
    np.ndarray
        Beam model in dB sampled at the pass coordinates.
    """
    if delays is None:
        delays = [0] * 16
    if amps is None:
        amps = [1.0] * 16

    az = np.radians(np.asarray(p.az_deg, dtype=float))
    za = np.radians(90.0 - np.asarray(p.alt_deg, dtype=float))

    jones = beam.calc_jones_array(
        az,
        za,
        freq_hz,
        delays,
        amps,
        norm_to_zenith,
        parallactic,
    )

    unpol_beam = make_unpol_instrumental_response(jones, jones)

    pol = p.tile[-2:]  # "XX" or "YY"

    if pol == "XX":
        model = 10.0 * np.log10(np.real(unpol_beam[:, 0]))
    elif pol == "YY":
        model = 10.0 * np.log10(np.real(unpol_beam[:, 3]))
    else:
        raise ValueError(f"Could not infer polarization from tile name: {p.tile}")

    return model


def passes_to_healpix(
    passes,
    nside: int = 32,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Bin pass samples onto a HEALPix map.

    Parameters
    ----------
    passes
        Iterable of pass records. Each record must provide ``alt_deg``,
        ``az_deg``, and ``power_db`` attributes.
    nside
        HEALPix NSIDE parameter controlling angular resolution.

    Returns
    -------
    hp_map
        HEALPix map containing the mean power per pixel in dB.
        Pixels with no contributing samples are NaN.
    counts
        Number of samples contributing to each pixel.
    mad_map
        HEALPix map of the median absolute deviation of samples in each pixel.
        Pixels with no contributing samples are NaN.

    Notes
    -----
    The median absolute deviation is scaled using ``scale="normal"`` so that
    it estimates the equivalent Gaussian standard deviation.
    """
    npix = hp.nside2npix(nside)

    sums = np.zeros(npix, dtype=float)
    counts = np.zeros(npix, dtype=int)
    values_by_pix: list[list[float]] = [[] for _ in range(npix)]

    for p in passes:
        alt_deg = np.asarray(p.alt_deg, dtype=float)
        az_deg = np.asarray(p.az_deg, dtype=float)
        power_db = np.asarray(p.power_db, dtype=float)

        mask = np.isfinite(alt_deg) & np.isfinite(az_deg) & np.isfinite(power_db)
        if not np.any(mask):
            continue

        theta = np.radians(90.0 - alt_deg[mask])
        phi = np.radians(az_deg[mask])

        pix = hp.ang2pix(nside, theta, phi)
        vals = power_db[mask]

        np.add.at(sums, pix, vals)
        np.add.at(counts, pix, 1)

        for pix_i, val_i in zip(pix, vals, strict=False):
            values_by_pix[pix_i].append(val_i)

    hp_map = np.full(npix, np.nan, dtype=float)
    mad_map = np.full(npix, np.nan, dtype=float)

    valid = counts > 0
    hp_map[valid] = sums[valid] / counts[valid]

    for pix_i, values in enumerate(values_by_pix):
        if values:
            mad_map[pix_i] = mad(values, scale="normal")

    return hp_map, counts, mad_map


def passes_to_healpix_concave_monotonic(
    passes,
    beam_file: str | Path,
    mwa_model_threshold_db: float = -10.0,
    residual_thresh: float = 2.0,
    nside: int = 32,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Bin monotonic-concave calibrated pass samples onto a HEALPix map.

    Samples from accepted satellite passes are assigned to HEALPix pixels
    using their altitude and azimuth coordinates. The output power map contains
    the median of all samples contributing to each pixel.

    Passes with a non-zero RFE calibration correction are compared with the
    MWA FEE beam model. A pass is rejected when its measured power exceeds the
    model by more than ``residual_thresh`` within the region where the model
    response is greater than ``mwa_model_threshold_db``.

    Parameters
    ----------
    passes
        Iterable of satellite-pass records. Each record must provide
        ``alt_deg``, ``az_deg``, ``power_db``, and ``cal_db`` attributes and
        must be compatible with :func:`pass_hyperbeam_model`.
    mwa_model_threshold_db
        Minimum MWA model response used when testing calibrated passes, in dB
        relative to the beam maximum. Default is -10 dB.
    residual_thresh
        Maximum permitted positive residual between the measured power and
        MWA beam model within the selected model region, in dB. Default is
        2 dB.
    nside
        HEALPix NSIDE parameter controlling the angular resolution. Default
        is 32.

    Returns
    -------
    power_map
        HEALPix map containing the median power per pixel in dB. Pixels with
        no contributing samples are NaN.
    count_map
        Number of samples contributing to each pixel.
    mad_map
        HEALPix map containing the median absolute deviation of the power
        samples in each pixel. Pixels with no contributing samples are NaN.

    Notes
    -----
    The residual threshold is one-sided: passes are rejected only when the
    measured power lies more than ``residual_thresh`` above the MWA model.

    The median absolute deviation is scaled using ``scale="normal"`` so that
    it estimates the equivalent Gaussian standard deviation.
    """
    npix = hp.nside2npix(nside)
    beam = mwa_hyperbeam.FEEBeam(str(beam_file))

    values_by_pixel: list[list[float]] = [[] for _ in range(npix)]

    for pass_record in passes:
        alt_deg = np.asarray(pass_record.alt_deg, dtype=float)
        az_deg = np.asarray(pass_record.az_deg, dtype=float)
        power_db = np.asarray(pass_record.power_db, dtype=float)
        calibration_db = np.asarray(pass_record.cal_db, dtype=float)

        mwa_model_db = np.asarray(
            pass_hyperbeam_model(pass_record, beam),
            dtype=float,
        )

        if not (
            alt_deg.shape
            == az_deg.shape
            == power_db.shape
            == calibration_db.shape
            == mwa_model_db.shape
        ):
            raise ValueError(
                "Pass altitude, azimuth, power, calibration, and model "
                "arrays must have matching shapes."
            )

        finite = (
            np.isfinite(alt_deg)
            & np.isfinite(az_deg)
            & np.isfinite(power_db)
            & np.isfinite(calibration_db)
            & np.isfinite(mwa_model_db)
        )

        if not np.any(finite):
            continue

        alt_deg = alt_deg[finite]
        az_deg = az_deg[finite]
        power_db = power_db[finite]
        calibration_db = calibration_db[finite]
        mwa_model_db = mwa_model_db[finite]

        calibration_applied = np.any(calibration_db > 0.0)

        if calibration_applied:
            model_region = mwa_model_db > mwa_model_threshold_db

            if np.any(model_region):
                residual_db = power_db[model_region] - mwa_model_db[model_region]

                if np.max(residual_db) > residual_thresh:
                    continue

        theta = np.radians(90.0 - alt_deg)
        phi = np.radians(az_deg)
        pixels = hp.ang2pix(nside, theta, phi)

        for pixel, value in zip(pixels, power_db, strict=False):
            values_by_pixel[pixel].append(float(value))

    power_map = np.full(npix, np.nan, dtype=float)
    mad_map = np.full(npix, np.nan, dtype=float)
    count_map = np.zeros(npix, dtype=int)

    for pixel, values in enumerate(values_by_pixel):
        if not values:
            continue

        pixel_values = np.asarray(values, dtype=float)

        count_map[pixel] = pixel_values.size
        power_map[pixel] = np.median(pixel_values)
        mad_map[pixel] = mad(pixel_values, scale="normal")

    return power_map, count_map, mad_map


def make_unpol_instrumental_response(j1: np.ndarray, j2: np.ndarray) -> np.ndarray:
    """
    Convert Jones matrices to unpolarized instrumental beam responses.

    Parameters
    ----------
    j1, j2 : np.ndarray
        Arrays of Jones terms with shape (N, 4), where the Jones matrix is
        flattened as [j00, j01, j10, j11].

    Returns
    -------
    np.ndarray
        Array of shape (N, 4) containing [XX, XY, YX, YY] instrumental
        responses for an unpolarized source.
    """
    result = np.empty_like(j1)

    result[:, 0] = j1[:, 0] * j2[:, 0].conjugate() + j1[:, 1] * j2[:, 1].conjugate()
    result[:, 1] = j1[:, 0] * j2[:, 2].conjugate() + j1[:, 1] * j2[:, 3].conjugate()
    result[:, 2] = j1[:, 2] * j2[:, 0].conjugate() + j1[:, 3] * j2[:, 1].conjugate()
    result[:, 3] = j1[:, 2] * j2[:, 2].conjugate() + j1[:, 3] * j2[:, 3].conjugate()

    return result


def hyperbeam_healpix(
    beam,
    nside: int = 32,
    freq_hz: float = 137e6,
    delays: np.ndarray | None = None,
    amps: np.ndarray | None = None,
    norm_to_zenith: bool = True,
    parallactic: bool = False,
) -> np.ndarray:
    """Evaluate Hyperbeam on an above-horizon HEALPix grid.

    Parameters
    ----------
    beam
        Initialized ``mwa_hyperbeam.FEEBeam`` instance.
    nside
        HEALPix NSIDE parameter controlling angular resolution.
    freq_hz
        Frequency at which to evaluate the beam in Hz.
    delays
        Beamformer delays for the 16 dipoles.
        If ``None``, all delays are set to zero.
    amps
        Dipole amplitudes for the 16 dipoles.
        If ``None``, all amplitudes are set to unity.
    norm_to_zenith
        If ``True``, normalize the beam to unity at zenith.
    parallactic
        If ``True``, apply parallactic-angle correction.

    Returns
    -------
    beam_db
        Array of shape ``(npix // 2, 2)`` containing above-horizon beam
        responses in dB. Column 0 contains XX polarization and column 1
        contains YY polarization.

    Notes
    -----
    The beam response is computed from the Jones matrices using
    ``makeUnpolInstrumentalResponse`` and converted to decibels via
    ``10 * log10(power)``.
    Any pixels with non-positive power will yield ``-inf``.
    """
    if delays is None:
        delays = np.zeros(16, dtype=int)

    if amps is None:
        amps = np.ones(16, dtype=float)

    delays = np.asarray(delays, dtype=int)
    amps = np.asarray(amps, dtype=float)

    if delays.shape != (16,):
        raise ValueError(f"Expected delays with shape (16,), got {delays.shape}.")

    if amps.shape != (16,):
        raise ValueError(f"Expected amps with shape (16,), got {amps.shape}.")

    npix = hp.nside2npix(nside)
    above_horizon = np.arange(npix // 2)

    za, az = hp.pix2ang(nside, above_horizon)

    jones = beam.calc_jones_array(
        az,
        za,
        freq_hz,
        delays,
        amps,
        norm_to_zenith,
        parallactic,
    )

    unpol_beam = make_unpol_instrumental_response(jones, jones)

    beam_xx = np.real(unpol_beam[:, 0])
    beam_yy = np.real(unpol_beam[:, 3])

    beam_db_xx = 10.0 * np.log10(beam_xx)
    beam_db_yy = 10.0 * np.log10(beam_yy)

    return np.column_stack([beam_db_xx, beam_db_yy])


class HyperbeamLogLike(Op):
    """PyTensor Op wrapping a Hyperbeam Gaussian log-likelihood.

    Parameters
    ----------
    beam_file
        Path to the Hyperbeam HDF5 beam model file.
    data_db
        Observed beam measurements in dB.
    sigma_db
        Per-pixel Gaussian uncertainties in dB.
    fit_mask
        Boolean mask selecting valid above-horizon pixels.
    nside
        HEALPix NSIDE parameter used for beam evaluation.
    freq_hz
        Frequency at which the beam is evaluated in Hz.
    pol
        Polarization to fit. Must be either ``"xx"`` or ``"yy"``.

    Notes
    -----
    This class wraps a black-box likelihood based on the Hyperbeam Rust
    library. Since Hyperbeam is not differentiable within PyTensor/JAX, the
    likelihood is evaluated directly inside ``perform``.

    The likelihood assumes independent Gaussian errors
    """

    itypes = [pt.dvector]
    otypes = [pt.dscalar]

    def __init__(
        self,
        beam_file: str | Path,
        data_db: np.ndarray,
        sigma_db: np.ndarray,
        fit_mask: np.ndarray,
        nside: int = 32,
        freq_hz: float = 137e6,
        pol: str = "xx",
    ) -> None:
        self.beam_file = str(beam_file)
        self.data_db = np.asarray(data_db, dtype=float)
        self.sigma_db = np.asarray(sigma_db, dtype=float)
        self.fit_mask = np.asarray(fit_mask, dtype=bool)
        self.nside = int(nside)
        self.freq_hz = float(freq_hz)
        self.pol = pol.lower()

        if self.pol not in {"xx", "yy"}:
            raise ValueError("pol must be either 'xx' or 'yy'.")

        if self.data_db.shape != self.sigma_db.shape:
            raise ValueError(
                "data_db and sigma_db must have the same shape. "
                f"Got {self.data_db.shape} and {self.sigma_db.shape}."
            )

    def perform(self, node: Any, inputs: list[np.ndarray], outputs: list[Any]) -> None:
        """Evaluate the Gaussian log-likelihood."""
        beam = mwa_hyperbeam.FEEBeam(self.beam_file)

        (theta,) = inputs

        model = hyperbeam_healpix(
            beam,
            nside=self.nside,
            freq_hz=self.freq_hz,
            amps=theta,
        )

        pol_idx = 0 if self.pol == "xx" else 1
        model_db = model[:, pol_idx][self.fit_mask]

        resid = self.data_db - model_db

        loglike = -0.5 * np.sum(
            (resid / self.sigma_db) ** 2 + np.log(2.0 * np.pi * self.sigma_db**2)
        )

        outputs[0][0] = np.asarray(loglike, dtype=float)


def clean_datatree_for_save(dt):
    """Remove object-dtype variables from a DataTree.

    Parameters
    ----------
    dt
        xarray DataTree containing PyMC/ArviZ inference outputs.

    Returns
    -------
    dt_clean
        Copy of the input DataTree with object-dtype variables removed.

    Notes
    -----
    Some PyMC SMC bookkeeping variables such as ``beta`` are stored as object
    arrays containing Python lists of varying lengths. These cannot be written
    cleanly to NetCDF or Zarr and are therefore removed.
    """
    dt = dt.copy()

    for path, node in list(dt.subtree_with_keys):
        ds = node.ds

        drop_vars = [name for name, da in ds.data_vars.items() if da.dtype == object]

        if drop_vars:
            print(f"Dropping from {path}: {drop_vars}")
            dt[path].ds = ds.drop_vars(drop_vars)

    return dt


def save_smc_stats(idata, path: Path) -> None:
    """Save useful SMC bookkeeping variables separately.

    Parameters
    ----------
    idata
        ArviZ inference data object returned by ``pm.sample_smc``.
    path
        Output pickle path.

    Notes
    -----
    This preserves useful SMC diagnostics such as:
    - beta annealing schedules
    - acceptance rates
    - log marginal likelihood estimates

    even when they are removed from NetCDF/Zarr outputs.
    """
    stats = idata["/sample_stats"].ds

    names = ["beta", "accept_rate", "log_marginal_likelihood"]
    smc_stats = {name: stats[name].values for name in names if name in stats}

    with open(path, "wb") as f:
        pickle.dump(smc_stats, f)


def build_fit_data(
    pass_file: str | Path,
    beam_file: str | Path,
    pointing: int,
    nside: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Construct above-horizon fitting vectors from pass data.

    Parameters
    ----------
    pass_file
        Path to the EMBERS pass HDF5 file.
    pointing
        Beamformer pointing index.
    nside
        HEALPix NSIDE parameter.

    Returns
    -------
    data_db
        Above-horizon observed beam values in dB.
    sigma_db
        Per-pixel MAD-based uncertainties in dB.
    fit_mask
        Boolean mask selecting valid above-horizon pixels.
    data_map
        Full-sky HEALPix data map.
    count_map
        Number of contributing samples per HEALPix pixel.

    Notes
    -----
    Pixels are considered valid if:
    - the beam value is finite,
    - the MAD estimate is finite,
    - the MAD estimate is positive.
    """
    pf = PassFile(pass_file)
    passes = pf.read_passes(pointing=pointing)

    # data_map, count_map, mad_map = passes_to_healpix(passes, nside=nside)
    data_map, count_map, mad_map = passes_to_healpix_concave_monotonic(
        passes, beam_file=beam_file, nside=nside
    )

    npix = hp.nside2npix(nside)
    above_horizon = np.arange(npix // 2)

    data_above = data_map[above_horizon]
    sigma_above = mad_map[above_horizon]

    fit_mask = np.isfinite(data_above) & np.isfinite(sigma_above) & (sigma_above > 0)

    data_db = np.asarray(data_above[fit_mask], dtype=float)
    sigma_db = np.asarray(sigma_above[fit_mask], dtype=float)
    fit_mask = np.asarray(fit_mask, dtype=bool)

    return data_db, sigma_db, fit_mask, data_map, count_map


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Fit MWA dipole gains using Hyperbeam and PyMC SMC."
    )

    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--draws", type=int, default=2048)
    parser.add_argument("--outdir", type=str, default="data/mcmc")

    parser.add_argument(
        "--pass-file",
        type=str,
        default="../passes/rf0/S06XX_rf0XX_passes.h5",
    )
    parser.add_argument(
        "--beam-file",
        type=str,
        default=os.environ.get("MWA_BEAM_FILE"),
    )
    parser.add_argument("--nside", type=int, default=32)
    parser.add_argument("--freq-hz", type=float, default=137e6)
    parser.add_argument("--pointing", type=int, default=0)
    parser.add_argument("--pol", choices=["xx", "yy"], default="xx")

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.beam_file is None:
        raise ValueError(
            "No beam file supplied. Pass --beam-file or set MWA_BEAM_FILE."
        )

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    metadata = {
        "seed": args.seed,
        "draws": args.draws,
        "pass_file": str(Path(args.pass_file).resolve()),
        "beam_file": str(Path(args.beam_file).resolve()),
        "nside": args.nside,
        "freq_hz": args.freq_hz,
        "pointing": args.pointing,
        "pol": args.pol,
        # "prior": "theta = 16 * Dirichlet(ones(16))",
        "prior": "theta = Truncated Gaussian, mu=1, sigma=0.3",
        "chains": 1,
        "cores": 1,
    }

    pass_stem = Path(args.pass_file).stem

    stem = f"{pass_stem[:11]}_{metadata['pointing']}_S{args.seed}_D{args.draws}"

    with open(outdir / f"{stem}_metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)

    data_db, sigma_db, fit_mask, _, _ = build_fit_data(
        pass_file=args.pass_file,
        beam_file=args.beam_file,
        pointing=args.pointing,
        nside=args.nside,
    )

    loglike_op = HyperbeamLogLike(
        beam_file=args.beam_file,
        data_db=data_db,
        sigma_db=sigma_db,
        fit_mask=fit_mask,
        nside=args.nside,
        freq_hz=args.freq_hz,
        pol=args.pol,
    )

    # with pm.Model() as model:
    #     p = pm.Dirichlet("p", a=np.ones(16))
    #     theta = pm.Deterministic("theta", 16.0 * p)
    #
    #     pm.Potential("beam_likelihood", loglike_op(theta))
    #
    #     idata = pm.sample_smc(
    #         draws=args.draws,
    #         chains=1,
    #         cores=1,
    #         random_seed=args.seed,
    #         progressbar=False,
    #     )
    #
    # with open(outdir / f"{stem}.pkl", "wb") as f:
    #     pickle.dump(idata, f)
    #
    # save_smc_stats(idata, outdir / f"{stem}_smc_stats.pkl")
    #
    # idata_save = clean_datatree_for_save(idata)
    # idata_save.to_netcdf(outdir / f"{stem}.nc")

    with pm.Model() as model:
        theta = pm.TruncatedNormal(
            "theta",
            mu=1.0,
            sigma=0.3,
            lower=0.0,
            shape=16,
        )

        theta_sum = pm.Deterministic("theta_sum", pt.sum(theta))
        theta_mean = pm.Deterministic("theta_mean", pt.mean(theta))
        theta_std = pm.Deterministic("theta_std", pt.std(theta))

        pm.Potential("beam_likelihood", loglike_op(theta))

        idata = pm.sample_smc(
            draws=args.draws,
            chains=1,
            cores=1,
            random_seed=args.seed,
            progressbar=False,
        )

    with open(outdir / f"{stem}.pkl", "wb") as f:
        pickle.dump(idata, f)

    save_smc_stats(idata, outdir / f"{stem}_smc_stats.pkl")

    idata_save = clean_datatree_for_save(idata)
    idata_save.to_netcdf(outdir / f"{stem}.nc")


if __name__ == "__main__":
    main()
