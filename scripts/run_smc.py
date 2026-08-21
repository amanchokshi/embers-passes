"""Run amplitude-only and amplitude+phase SMC fits for one pass file."""

from __future__ import annotations

import argparse
import pickle
import re
import sys
import time
from dataclasses import asdict
from datetime import UTC, datetime
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("-o", "--output", type=Path)
    parser.add_argument("--chunk-index", type=int, default=0)
    parser.add_argument("--chunk-days", type=float, default=200.0)
    parser.add_argument("--pointing", type=int, default=0)
    parser.add_argument("--nside", type=int, default=32)
    parser.add_argument("--n-devices", type=int, default=12)
    parser.add_argument("--particles-per-device", type=int, default=192)
    parser.add_argument("--dirichlet-alpha", type=float, default=10.0)
    parser.add_argument("--phase-prior-std-deg", type=float, default=10.0)
    parser.add_argument("--prior-seed", type=int, default=2030)
    parser.add_argument("--smc-seed", type=int, default=2040)
    parser.add_argument("--final-mutation-seed", type=int, default=2041)
    parser.add_argument("--quiet", action="store_true")
    return parser.parse_args()


def infer_tile_pol(path):
    match = re.search(r"(S\d{2})(XX|YY)", path.name, re.IGNORECASE)
    if match is None:
        raise ValueError(f"Could not infer tile/pol from {path.name!r}.")
    tile, pol = match.group(1).upper(), match.group(2).upper()
    return tile, pol, {"XX": 0, "YY": 1}[pol]


def gaussian_log_likelihood(data, model, sigma):
    import numpy as np
    return float(
        -0.5 * np.sum(((data - model) / sigma) ** 2 + np.log(2 * np.pi * sigma**2))
    )


def package_result(smc, model, result, n_data):
    import numpy as np

    initial_state = smc.state_to_numpy(result["initial_state"])
    posterior_state = smc.state_to_numpy(result["posterior_state"])
    posterior_sites = smc.physical_sites_to_numpy(
        smc.physical_sites(model, result["posterior_state"].particles)
    )

    best_log_likelihood = float(np.max(posterior_state["log_likelihood"]))
    n_parameters = int(model["latent_dim"])

    return {
        "latent_dim": n_parameters,
        "initial_state": initial_state,
        "posterior_state": posterior_state,
        "posterior_sites": posterior_sites,
        "stage_diagnostics": smc.diagnostics_to_numpy(result["stage_diagnostics"]),
        "final_mutation_diagnostics": smc.diagnostics_to_numpy(
            result["final_mutation_diagnostics"]
        ),
        "cumulative_log_evidence": np.asarray(result["cumulative_log_evidence"]),
        "log_evidence": float(result["log_evidence"]),
        "final_nuts_step_size": float(result["final_nuts_step_size"]),
        "posterior_nuts_step_size": float(result["posterior_nuts_step_size"]),
        "best_log_likelihood": best_log_likelihood,
        "bic": float(n_parameters * np.log(n_data) - 2 * best_log_likelihood),
    }


def main():
    args = parse_args()

    input_path = args.input.expanduser().resolve()
    if not input_path.is_file():
        raise FileNotFoundError(input_path)
    if args.chunk_index < 0 or args.chunk_days <= 0:
        raise ValueError("Chunk index must be non-negative and chunk days positive.")

    output_path = (
        args.output.expanduser().resolve()
        if args.output
        else input_path.with_name(f"{input_path.stem}_smc.pkl")
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # This must happen before JAX initializes its CPU backend.
    import numpyro
    numpyro.set_host_device_count(args.n_devices)

    import healpy as hp
    import jax
    import jax.numpy as jnp
    import numpy as np
    from mwa_jaxbeam import aee_137mhz as aee
    from embers_passes import PassFile, chunk_passes, chunks_to_healpix_counts
    from embers_passes import smc

    tile, pol_name, pol = infer_tile_pol(input_path)
    n_particles = args.n_devices * args.particles_per_device
    phase_std_rad = np.deg2rad(args.phase_prior_std_deg)
    phase_kappa = float(1 / phase_std_rad**2)

    if jax.local_device_count() < args.n_devices:
        raise RuntimeError(
            f"Requested {args.n_devices} devices; only {jax.local_device_count()} available."
        )

    print(f"{tile}{pol_name}: {n_particles} particles on {args.n_devices} devices")

    # ------------------------------------------------------------------
    # Data and likelihood pixels
    # ------------------------------------------------------------------

    passes = PassFile(input_path).read_passes(pointing=args.pointing)
    chunks, bin_edges, bin_centers = chunk_passes(
        passes,
        chunk_sec=args.chunk_days * 86400,
    )
    if args.chunk_index >= len(chunks):
        raise IndexError(f"Chunk {args.chunk_index} requested; only {len(chunks)} exist.")

    median_maps, count_maps, mad_maps, median_sem_maps = chunks_to_healpix_counts(
        passes,
        chunks,
        decimate=1,
        nside=args.nside,
    )

    data_map = np.asarray(
        median_maps[args.chunk_index],
        dtype=np.float32,
    )
    count_map = np.asarray(
        count_maps[args.chunk_index],
    )
    error_map = np.sqrt(
        np.asarray(mad_maps[args.chunk_index], dtype=np.float32) ** 2
        + np.asarray(median_sem_maps[args.chunk_index], dtype=np.float32) ** 2
    ).astype(np.float32)

    npix = hp.nside2npix(args.nside)
    pixel_indices = np.arange(npix)
    za_hpx, az_hpx = hp.pix2ang(args.nside, pixel_indices)

    visible = (
        (za_hpx <= np.pi / 2)
        & (count_map >= 3)
    )

    visible_indices = np.flatnonzero(visible)

    reference_power = np.full(npix, np.nan, dtype=np.float32)
    reference_power[visible] = np.asarray(
        aee.power(
            az_rad=az_hpx[visible],
            za_rad=za_hpx[visible],
            normalize=True,
        )[pol],
        dtype=np.float32,
    )

    valid = (
        visible
        & np.isfinite(data_map)
        & np.isfinite(error_map)
        & (error_map > 0)
        & np.isfinite(reference_power)
        & (reference_power > 0)
    )
    comparison_indices = np.flatnonzero(valid)
    if comparison_indices.size == 0:
        raise ValueError("No valid comparison pixels.")

    visible_lookup = np.full(npix, -1, dtype=np.int32)
    visible_lookup[visible_indices] = np.arange(visible_indices.size, dtype=np.int32)
    comparison_visible_indices = visible_lookup[comparison_indices]

    az_visible = jnp.asarray(az_hpx[visible], dtype=jnp.float32)
    za_visible = jnp.asarray(za_hpx[visible], dtype=jnp.float32)
    comparison_visible_indices = jnp.asarray(comparison_visible_indices, dtype=jnp.int32)
    observed_db = jnp.asarray(data_map[comparison_indices], dtype=jnp.float32)
    observed_sigma_db = jnp.asarray(error_map[comparison_indices], dtype=jnp.float32)

    # ------------------------------------------------------------------
    # Shared beam forward model
    # ------------------------------------------------------------------

    @jax.jit
    def predict_beam_db(dipole_gains_pol):
        fixed = jnp.ones(16, dtype=jnp.complex64)
        dipole_gains_pol = dipole_gains_pol.astype(jnp.complex64)
        dipole_gains = (
            jnp.stack([dipole_gains_pol, fixed])
            if pol == 0
            else jnp.stack([fixed, dipole_gains_pol])
        )
        power = aee.power(
            az_rad=az_visible,
            za_rad=za_visible,
            dipole_gains=dipole_gains,
            normalize=True,
        )[pol, comparison_visible_indices]
        power = jnp.maximum(power, jnp.finfo(power.dtype).tiny)
        return 10 * jnp.log10(power)

    nominal_model_db = np.asarray(
        predict_beam_db(jnp.ones(16, dtype=jnp.complex64))
    )
    observed_db_np = np.asarray(observed_db)
    observed_sigma_np = np.asarray(observed_sigma_db)
    nominal_residual = observed_db_np - nominal_model_db
    nominal_log_likelihood = gaussian_log_likelihood(
        observed_db_np, nominal_model_db, observed_sigma_np
    )
    nominal_chi2 = float(np.sum((nominal_residual / observed_sigma_np) ** 2))

    # ------------------------------------------------------------------
    # SMC fits
    # ------------------------------------------------------------------

    config = smc.SMCConfig(n_particles=n_particles, n_devices=args.n_devices)
    model_kwargs = dict(
        n_dipoles=16,
        dirichlet_alpha=args.dirichlet_alpha,
        interface_seed=2026,
        dtype=jnp.float32,
    )
    run_kwargs = dict(
        prior_seed=args.prior_seed,
        smc_seed=args.smc_seed,
        final_mutation_seed=args.final_mutation_seed,
        print_progress=not args.quiet,
    )

    fits = {}
    for kind in ("amplitude", "complex"):
        print(f"\nRunning {kind} SMC")
        model = smc.build_model_interface(
            predict_beam_db,
            observed_db,
            observed_sigma_db,
            kind=kind,
            phase_kappa=phase_kappa if kind == "complex" else None,
            **model_kwargs,
        )
        start = time.perf_counter()
        raw_result = smc.run_smc(model, config, **run_kwargs)
        fits[kind] = package_result(smc, model, raw_result, observed_db.size)
        fits[kind]["runtime_seconds"] = time.perf_counter() - start
        del raw_result

    # ------------------------------------------------------------------
    # Pickle only the expensive inference products and the canonical data.
    # ------------------------------------------------------------------

    config_dict = asdict(config)
    config_dict["dtype"] = "float32"

    result = {
        "metadata": {
            "created_utc": datetime.now(UTC).isoformat(),
            "input_path": str(input_path),
            "tile": tile,
            "polarization": pol_name,
            "pol_index": pol,
            "pointing": args.pointing,
            "chunk_index": args.chunk_index,
            "chunk_days": args.chunk_days,
            "chunk_start_unix": float(bin_edges[args.chunk_index]),
            "chunk_stop_unix": float(bin_edges[args.chunk_index + 1]),
            "chunk_center_unix": float(bin_centers[args.chunk_index]),
            "nside": args.nside,
            "dirichlet_alpha": args.dirichlet_alpha,
            "phase_prior_std_deg": args.phase_prior_std_deg,
            "phase_kappa": phase_kappa,
            "prior_seed": args.prior_seed,
            "smc_seed": args.smc_seed,
            "final_mutation_seed": args.final_mutation_seed,
            "smc_config": config_dict,
            "command": sys.argv,
        },
        "data": {
            "data_map": data_map,
            "error_map": error_map,
            "comparison_indices": comparison_indices.astype(np.int32),
        },
        "nominal": {
            "log_likelihood": nominal_log_likelihood,
            "chi2": nominal_chi2,
            "reduced_chi2": nominal_chi2 / observed_db.size,
            "residual_rms_db": float(np.sqrt(np.mean(nominal_residual**2))),
            "bic": float(-2 * nominal_log_likelihood),
        },
        **fits,
    }

    with output_path.open("wb") as handle:
        pickle.dump(result, handle, protocol=pickle.HIGHEST_PROTOCOL)

    print(f"\nSaved: {output_path}")
    print(f"Comparison pixels: {comparison_indices.size}")
    print(f"BIC nominal / amplitude / complex: "
          f"{result['nominal']['bic']:.1f} / {fits['amplitude']['bic']:.1f} / "
          f"{fits['complex']['bic']:.1f}")


if __name__ == "__main__":
    main()
