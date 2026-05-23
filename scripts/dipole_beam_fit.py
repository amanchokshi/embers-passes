import argparse
import os
import pickle
from pathlib import Path

import pymc as pm
import numpy as np
import healpy as hp
import mwa_hyperbeam
import pytensor.tensor as pt
from pytensor.graph.op import Op
from scipy.stats import median_abs_deviation as mad

from embers_passes import PassFile


def passes_to_healpix(
    passes,
    nside: int = 32,
    nest: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Bin pass samples onto a HEALPix map.

    Parameters
    ----------
    passes
        Iterable of pass records. Each record must have ``alt_deg``,
        ``az_deg``, and ``power_db`` attributes.
    nside
        HEALPix NSIDE parameter. Default is 32.
    nest
        Use nested ordering if ``True``. Default is ``False``, corresponding
        to ring ordering.

    Returns
    -------
    hp_map
        HEALPix map of mean binned power values. Pixels with no samples are
        NaN.
    counts
        Number of samples contributing to each pixel.
    mad_map
        HEALPix map of the median absolute deviation of samples in each pixel.
        Pixels with no samples are NaN.
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

        pix = hp.ang2pix(nside, theta, phi, nest=nest)
        vals = power_db[mask]

        np.add.at(sums, pix, vals)
        np.add.at(counts, pix, 1)

        for pix_i, val_i in zip(pix, vals):
            values_by_pix[pix_i].append(val_i)

    hp_map = np.full(npix, np.nan, dtype=float)
    mad_map = np.full(npix, np.nan, dtype=float)

    valid = counts > 0
    hp_map[valid] = sums[valid] / counts[valid]

    for pix_i, values in enumerate(values_by_pix):
        if values:
            mad_map[pix_i] = mad(values, scale="normal")

    return hp_map, counts, mad_map


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
    delays=np.zeros(16, dtype=int),
    amps=np.ones(16),
    norm_to_zenith: bool = True,
    parallactic: bool = False,
) -> np.ndarray:
    """Evaluate the MWA Hyperbeam model on an above-horizon HEALPix grid.
    The beam is evaluated for all HEALPix pixels with zenith angle
    0 <= ZA <= pi/2, corresponding to the upper hemisphere. For a HEALPix
    map with npix = 12 * nside**2 pixels, this includes npix // 2 pixels.

    Parameters
    ----------
    beam
        Initialized mwa_hyperbeam.FEEBeam object.
    nside
        HEALPix NSIDE parameter defining the angular resolution of the grid.
    freq_hz
        Frequency at which to evaluate the beam, in Hz.
    delays
        List of 16 beamformer delays, one per dipole. If None, all
        delays are set to zero.
    amps
        List of 16 dipole amplitudes. If ``None``, all amplitudes are set
        to unity.
    norm_to_zenith
        If True, normalize the beam response to unity at zenith.
    parallactic
        If True, include parallactic-angle rotation.

    Returns
    -------
    beam_db
        Array of shape ``(npix // 2, 2)`` containing the unpolarized beam
        response in dB for all above-horizon pixels. Column 0 contains the
        XX polarization and column 1 contains the YY polarization.

    Notes
    -----
    The beam response is computed from the Jones matrices using
    ``makeUnpolInstrumentalResponse`` and converted to decibels via
    ``10 * log10(power)``.
    Any pixels with non-positive power will yield ``-inf``.
    """

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
    itypes = [pt.dvector]
    otypes = [pt.dscalar]

    def __init__(self, beam_file, data_db, sigma_db, valid_mask):
        self.beam_file = str(beam_file)
        self.data_db = np.asarray(data_db, dtype=float)
        self.sigma_db = np.asarray(sigma_db, dtype=float)
        self.valid_mask = np.asarray(valid_mask, dtype=bool)

    def perform(self, node, inputs, outputs):
        beam = mwa_hyperbeam.FEEBeam(self.beam_file)

        (theta,) = inputs

        model_db = hyperbeam_healpix(
            beam,
            nside=32,
            freq_hz=137e6,
            amps=theta,
        )[:, 0]

        model_db = model_db[self.valid_mask]
        resid = self.data_db - model_db

        loglike = -0.5 * np.sum(
            (resid / self.sigma_db) ** 2 + np.log(2.0 * np.pi * self.sigma_db**2)
        )

        outputs[0][0] = np.asarray(loglike, dtype=float)


def clean_datatree_for_save(dt):
    """Drop object-dtype variables from every DataTree node."""
    dt = dt.copy()

    for path, node in list(dt.subtree_with_keys):
        ds = node.ds

        drop_vars = [name for name, da in ds.data_vars.items() if da.dtype == object]

        if drop_vars:
            print(f"Dropping from {path}: {drop_vars}")
            dt[path].ds = ds.drop_vars(drop_vars)

    return dt


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--draws", type=int, default=2048)
    parser.add_argument("--outdir", type=str, default="data/mcmc")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    pf = PassFile("../passes/rf0/S06XX_rf0XX_passes.h5")
    passes = pf.read_passes(pointing=0)
    S06XX_map, count_map, mad_map = passes_to_healpix(passes, nside=32)

    # Restrict to above horizon pix
    npix = hp.nside2npix(32)
    above_horizon = np.arange(npix // 2)

    data_db = S06XX_map[above_horizon]
    sigma_db = mad_map[above_horizon]

    # Mask bad pix
    valid_mask = np.isfinite(data_db) & np.isfinite(sigma_db) & (sigma_db > 0)

    data_db = np.asarray(data_db[valid_mask])
    sigma_db = np.asarray(sigma_db[valid_mask])
    valid_mask = np.asarray(valid_mask)

    loglike_op = HyperbeamLogLike(
        beam_file=os.environ["MWA_BEAM_FILE"],
        data_db=data_db,
        sigma_db=sigma_db,
        valid_mask=valid_mask,
    )

    with pm.Model() as model:
        p = pm.Dirichlet("p", a=np.ones(16))
        theta = pm.Deterministic("theta", 16.0 * p)

        pm.Potential("beam_likelihood", loglike_op(theta))

        idata = pm.sample_smc(
            draws=args.draws,
            chains=1,
            cores=1,
            random_seed=args.seed,
            progressbar=False,
        )

    with open(outdir / f"smc_dirichlet_seed{args.seed}.pkl", "wb") as f:
        pickle.dump(idata, f)

    idata_save = clean_datatree_for_save(idata)
    idata_save.to_netcdf(outdir / f"smc_dirichlet_seed{args.seed}.nc")
