"""Run BlackJAX amplitude-only and amplitude+phase SMC fits for one pass file.

This is the project-facing launcher for :mod:`embers_passes.smc_blackjax`.
Machine-specific execution choices (logical CPU device count, particles per
logical device, seeds, paths) live here; the inference algorithm lives in
``smc_blackjax.py``.

The default run mirrors the validated notebook configuration:

- BlackJAX adaptive tempered SMC
- NUTS rejuvenation
- BlackJAX SMC inner-kernel step-size tuning
- identity inverse mass matrix
- max_num_doublings = 9
- pmap(vmap(NUTS)) across the requested logical CPU devices
"""

from __future__ import annotations

import argparse
import os
import pickle
import re
import sys
from dataclasses import asdict
from datetime import UTC, datetime
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("-o", "--output", type=Path)
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Directory for the default output filename (ignored with --output).",
    )
    parser.add_argument(
        "--tag",
        type=str,
        help="Optional label inserted into the default output filename.",
    )

    parser.add_argument("--chunk-index", type=int, default=0)
    parser.add_argument("--chunk-days", type=float, default=200.0)
    parser.add_argument("--pointing", type=int, default=0)
    parser.add_argument("--nside", type=int, default=32)
    parser.add_argument(
        "--single-pass-z-cut",
        type=float,
        default=None,
        help=(
            "Optionally mask anomalously low-uncertainty pixels with exactly "
            "one contributing pass. The cut is applied to log10(total sigma) "
            "using a pooled robust median/MAD across all chunks. For example, "
            "--single-pass-z-cut -3 applies a 3-sigma lower-tail cut. "
            "Default: disabled."
        ),
    )

    # Execution configuration.  These are deliberately NOT baked into the
    # scientific SMC module.
    parser.add_argument(
        "--devices", "--n-devices",
        dest="n_devices",
        type=int,
        default=12,
        help="Number of logical CPU devices used by pmap (default: 12).",
    )
    parser.add_argument(
        "--particles",
        type=int,
        help=(
            "Total SMC particle count. Must be divisible by --devices. "
            "If omitted, uses --devices * --particles-per-device."
        ),
    )
    parser.add_argument(
        "--particles-per-device",
        type=int,
        default=32,
        help="Particles per logical device when --particles is omitted (default: 32).",
    )

    # Model priors.
    parser.add_argument("--dirichlet-alpha", type=float, default=10.0)
    parser.add_argument("--phase-prior-std-deg", type=float, default=10.0)

    # BlackJAX SMC / NUTS settings.
    parser.add_argument("--target-ess", type=float, default=0.80)
    parser.add_argument("--target-accept", type=float, default=0.80)
    parser.add_argument("--max-doublings", type=int, default=9)
    parser.add_argument("--max-stages", type=int, default=256)
    parser.add_argument("--step-search-initial-scale", type=float, default=1.0)

    # A single seed defines the complete BlackJAX run.  The module splits it
    # internally into prior, step-size-search, and SMC substreams.
    parser.add_argument("--seed", type=int, default=314159)

    # Usually run both models for direct model comparison, but allow focused
    # regression / scaling runs when desired.
    parser.add_argument(
        "--kind",
        choices=("both", "amplitude", "complex"),
        default="both",
    )
    parser.add_argument("--quiet", action="store_true")
    return parser.parse_args()


def configure_xla_cpu_devices(n_devices: int) -> None:
    """Expose ``n_devices`` logical CPU devices before JAX is imported."""

    if n_devices < 1:
        raise ValueError("n_devices must be positive.")

    flag = f"--xla_force_host_platform_device_count={n_devices}"
    existing = os.environ.get("XLA_FLAGS", "").strip()

    # Avoid silently retaining a conflicting device-count flag inherited from
    # the shell / notebook environment.
    pattern = r"--xla_force_host_platform_device_count=\d+"
    if re.search(pattern, existing):
        existing = re.sub(pattern, flag, existing)
        os.environ["XLA_FLAGS"] = existing
    elif existing:
        os.environ["XLA_FLAGS"] = f"{existing} {flag}"
    else:
        os.environ["XLA_FLAGS"] = flag


def infer_tile_pol(path: Path):
    match = re.search(r"(S\d{2})(XX|YY)", path.name, re.IGNORECASE)
    if match is None:
        raise ValueError(f"Could not infer tile/pol from {path.name!r}.")
    tile, pol = match.group(1).upper(), match.group(2).upper()
    return tile, pol, {"XX": 0, "YY": 1}[pol]


def gaussian_log_likelihood(data, model, sigma):
    import numpy as np

    return float(
        -0.5
        * np.sum(
            ((data - model) / sigma) ** 2
            + np.log(2 * np.pi * sigma**2)
        )
    )


def package_fit(smc, model, result, n_data: int) -> dict:
    """Convert one ``SMCResult`` into long-lived analysis products."""

    import numpy as np

    packed = smc.result_to_dict(result)

    posterior_sites = smc.physical_sites_to_numpy(
        smc.physical_sites(model, result.particles)
    )

    best_log_likelihood = float(np.max(result.log_likelihood))
    n_parameters = int(model["latent_dim"])

    packed.update(
        {
            "latent_dim": n_parameters,
            "posterior_sites": posterior_sites,
            "best_log_likelihood": best_log_likelihood,
            "bic": float(
                n_parameters * np.log(n_data)
                - 2.0 * best_log_likelihood
            ),
        }
    )
    return packed


def contributing_pass_count_map(passes, chunk, nside: int):
    """Return number of independent passes contributing to each HEALPix pixel."""

    import healpy as hp
    import numpy as np

    npix = hp.nside2npix(nside)
    n_pass_map = np.zeros(npix, dtype=np.int32)

    for pass_index in chunk:
        pass_record = passes[pass_index]

        alt_deg = np.asarray(pass_record.alt_deg, dtype=float)
        az_deg = np.asarray(pass_record.az_deg, dtype=float)
        power_db = np.asarray(pass_record.power_db, dtype=float)

        if not (
            alt_deg.shape
            == az_deg.shape
            == power_db.shape
        ):
            raise ValueError(
                "alt_deg, az_deg, and power_db must have matching shapes."
            )

        valid = (
            np.isfinite(alt_deg)
            & np.isfinite(az_deg)
            & np.isfinite(power_db)
            & (alt_deg > 0.0)
        )

        if not np.any(valid):
            continue

        theta = np.radians(90.0 - alt_deg[valid])
        phi = np.radians(az_deg[valid])

        pixels = hp.ang2pix(
            nside,
            theta,
            phi,
        )

        n_pass_map[np.unique(pixels)] += 1

    return n_pass_map


def main():
    args = parse_args()

    input_path = args.input.expanduser().resolve()
    if not input_path.is_file():
        raise FileNotFoundError(input_path)
    if args.chunk_index < 0:
        raise ValueError("chunk_index must be non-negative.")
    if args.chunk_days <= 0:
        raise ValueError("chunk_days must be positive.")
    if args.n_devices < 1:
        raise ValueError("devices must be positive.")
    if args.particles_per_device < 1:
        raise ValueError("particles_per_device must be positive.")
    if args.single_pass_z_cut is not None and args.single_pass_z_cut >= 0:
        raise ValueError(
            "--single-pass-z-cut must be negative; "
            "for example, --single-pass-z-cut -3."
        )

    n_particles = (
        args.particles
        if args.particles is not None
        else args.n_devices * args.particles_per_device
    )
    if n_particles < 2:
        raise ValueError("particles must be at least two.")
    if n_particles % args.n_devices != 0:
        raise ValueError(
            f"{n_particles} particles cannot be divided evenly across "
            f"{args.n_devices} devices."
        )
    particles_per_device = n_particles // args.n_devices

    if args.output is not None:
        output_path = args.output.expanduser().resolve()
    else:
        output_dir = (
            args.output_dir.expanduser().resolve()
            if args.output_dir is not None
            else input_path.parent
        )
        tag_suffix = f"_{args.tag}" if args.tag else ""
        output_path = output_dir / (
            f"{input_path.stem}{tag_suffix}_smc_blackjax.pkl"
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # IMPORTANT: configure logical CPU devices before importing JAX,
    # BlackJAX, NumPyro, mwa_jaxbeam, or smc_blackjax.
    # ------------------------------------------------------------------

    configure_xla_cpu_devices(args.n_devices)

    import healpy as hp
    import jax
    import jax.numpy as jnp
    import numpy as np
    from mwa_jaxbeam import aee_137mhz as aee
    from embers_passes import PassFile, chunk_passes, chunks_to_healpix_counts
    from embers_passes import smc_blackjax as smc

    tile, pol_name, pol = infer_tile_pol(input_path)

    phase_std_rad = np.deg2rad(args.phase_prior_std_deg)
    if phase_std_rad <= 0:
        raise ValueError("phase_prior_std_deg must be positive.")
    phase_kappa = float(1.0 / phase_std_rad**2)

    cpu_devices = jax.devices("cpu")
    if len(cpu_devices) < args.n_devices:
        raise RuntimeError(
            f"Requested {args.n_devices} CPU devices; "
            f"JAX exposes only {len(cpu_devices)}."
        )

    print(
        f"{tile}{pol_name}: {n_particles} particles "
        f"({particles_per_device}/device) on {args.n_devices} CPU devices"
    )
    print(
        f"BlackJAX settings: ESS={args.target_ess:.2f}, "
        f"target_accept={args.target_accept:.2f}, "
        f"max_doublings={args.max_doublings} "
        f"(ceiling={2**args.max_doublings - 1})"
    )

    # ------------------------------------------------------------------
    # Data and likelihood pixels
    # ------------------------------------------------------------------

    passes = PassFile(input_path).read_passes(pointing=args.pointing)
    chunks, bin_edges, bin_centers = chunk_passes(
        passes,
        chunk_sec=args.chunk_days * 86400,
    )
    if args.chunk_index >= len(chunks):
        raise IndexError(
            f"Chunk {args.chunk_index} requested; only {len(chunks)} exist."
        )

    median_maps, count_maps, mad_maps, median_sem_maps = chunks_to_healpix_counts(
        passes,
        chunks,
        decimate=1,
        nside=args.nside,
    )

    n_pass_maps = [
        contributing_pass_count_map(
            passes,
            chunk,
            args.nside,
        )
        for chunk in chunks
    ]

    single_pass_cut_metadata = {
        "enabled": args.single_pass_z_cut is not None,
        "robust_z_cut": args.single_pass_z_cut,
        "pooled_n_single_pass_pixels": 0,
        "median_log10_sigma": None,
        "robust_sigma_log10_sigma": None,
        "sigma_cut_db": None,
    }

    if args.single_pass_z_cut is not None:
        pooled_log_sigma = []

        for i in range(len(chunks)):
            counts_i = np.asarray(count_maps[i])
            n_pass_i = np.asarray(n_pass_maps[i])

            sigma_i = np.sqrt(
                np.asarray(mad_maps[i], dtype=np.float64) ** 2
                + np.asarray(median_sem_maps[i], dtype=np.float64) ** 2
            )

            valid_single = (
                (counts_i >= 3)
                & (n_pass_i == 1)
                & np.isfinite(sigma_i)
                & (sigma_i > 0)
            )

            pooled_log_sigma.append(np.log10(sigma_i[valid_single]))

        pooled_log_sigma = np.concatenate(pooled_log_sigma)

        if pooled_log_sigma.size < 10:
            raise ValueError(
                "Too few valid single-pass pixels to define the robust "
                f"uncertainty cut ({pooled_log_sigma.size} found)."
            )

        median_log_sigma = float(np.median(pooled_log_sigma))
        robust_sigma_log = float(
            1.482602218505602
            * np.median(
                np.abs(
                    pooled_log_sigma
                    - median_log_sigma
                )
            )
        )

        if not np.isfinite(robust_sigma_log) or robust_sigma_log <= 0:
            raise ValueError(
                "Could not estimate a finite positive robust scatter for "
                "single-pass log uncertainties."
            )

        log_sigma_cut = (
            median_log_sigma
            + args.single_pass_z_cut * robust_sigma_log
        )
        sigma_cut_db = float(10.0 ** log_sigma_cut)

        single_pass_cut_metadata.update(
            {
                "pooled_n_single_pass_pixels": int(pooled_log_sigma.size),
                "median_log10_sigma": median_log_sigma,
                "robust_sigma_log10_sigma": robust_sigma_log,
                "sigma_cut_db": sigma_cut_db,
            }
        )

        print("\nSingle-pass low-uncertainty quality cut")
        print("---------------------------------------")
        print(f"Pooled single-pass pixels: {pooled_log_sigma.size}")
        print(f"Median log10(sigma):       {median_log_sigma:.6f}")
        print(f"Robust scatter:            {robust_sigma_log:.6f}")
        print(f"Robust z cut:              {args.single_pass_z_cut:.2f}")
        print(f"Equivalent sigma cut:      {sigma_cut_db:.6f} dB")

    data_map = np.asarray(
        median_maps[args.chunk_index],
        dtype=np.float32,
    )
    count_map = np.asarray(count_maps[args.chunk_index])
    error_map = np.sqrt(
        np.asarray(mad_maps[args.chunk_index], dtype=np.float32) ** 2
        + np.asarray(median_sem_maps[args.chunk_index], dtype=np.float32) ** 2
    ).astype(np.float32)

    n_pass_map = np.asarray(
        n_pass_maps[args.chunk_index],
        dtype=np.int32,
    )

    single_pass_robust_z = np.full(
        count_map.shape,
        np.nan,
        dtype=np.float64,
    )
    single_pass_mask = np.zeros(
        count_map.shape,
        dtype=bool,
    )

    if args.single_pass_z_cut is not None:
        one_pass_valid = (
            (count_map >= 3)
            & (n_pass_map == 1)
            & np.isfinite(error_map)
            & (error_map > 0)
        )

        single_pass_robust_z[one_pass_valid] = (
            (
                np.log10(error_map[one_pass_valid])
                - single_pass_cut_metadata["median_log10_sigma"]
            )
            / single_pass_cut_metadata[
                "robust_sigma_log10_sigma"
            ]
        )

        single_pass_mask = (
            one_pass_valid
            & (single_pass_robust_z < args.single_pass_z_cut)
        )

    masked_indices = np.flatnonzero(single_pass_mask)

    print(
        f"\nSingle-pass uncertainty mask: "
        f"{masked_indices.size} pixels removed "
        f"from chunk {args.chunk_index}"
    )

    if masked_indices.size:
        print(
            "pixel       count  passes       sigma[dB]       robust_z"
        )
        print(
            "---------------------------------------------------------"
        )

        for pixel in masked_indices:
            print(
                f"{pixel:5d} "
                f"{int(count_map[pixel]):11d} "
                f"{int(n_pass_map[pixel]):7d} "
                f"{float(error_map[pixel]):15.6e} "
                f"{float(single_pass_robust_z[pixel]):12.3f}"
            )

    npix = hp.nside2npix(args.nside)
    pixel_indices = np.arange(npix)
    za_hpx, az_hpx = hp.pix2ang(args.nside, pixel_indices)

    visible = (za_hpx <= np.pi / 2) & (count_map >= 3)
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
        & ~single_pass_mask
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
    visible_lookup[visible_indices] = np.arange(
        visible_indices.size,
        dtype=np.int32,
    )
    comparison_visible_indices = visible_lookup[comparison_indices]

    az_visible = jnp.asarray(az_hpx[visible], dtype=jnp.float32)
    za_visible = jnp.asarray(za_hpx[visible], dtype=jnp.float32)
    comparison_visible_indices_jax = jnp.asarray(
        comparison_visible_indices,
        dtype=jnp.int32,
    )
    observed_db = jnp.asarray(
        data_map[comparison_indices],
        dtype=jnp.float32,
    )
    observed_sigma_db = jnp.asarray(
        error_map[comparison_indices],
        dtype=jnp.float32,
    )


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
        )[pol, comparison_visible_indices_jax]
        power = jnp.maximum(power, jnp.finfo(power.dtype).tiny)
        return 10.0 * jnp.log10(power)

    # ------------------------------------------------------------------
    # Nominal model diagnostics
    # ------------------------------------------------------------------

    nominal_model_db = np.asarray(
        predict_beam_db(jnp.ones(16, dtype=jnp.complex64))
    )
    observed_db_np = np.asarray(observed_db)
    observed_sigma_np = np.asarray(observed_sigma_db)
    nominal_residual = observed_db_np - nominal_model_db
    nominal_log_likelihood = gaussian_log_likelihood(
        observed_db_np,
        nominal_model_db,
        observed_sigma_np,
    )
    nominal_chi2 = float(
        np.sum((nominal_residual / observed_sigma_np) ** 2)
    )

    # ------------------------------------------------------------------
    # BlackJAX SMC configuration
    # ------------------------------------------------------------------

    config = smc.SMCConfig(
        n_particles=n_particles,
        target_ess=args.target_ess,
        target_accept=args.target_accept,
        max_num_doublings=args.max_doublings,
        max_stages=args.max_stages,
        num_mcmc_steps=1,
        step_search_initial_scale=args.step_search_initial_scale,
        dtype=jnp.float32,
    )
    execution = smc.ExecutionConfig(
        n_devices=args.n_devices,
        platform="cpu",
    )

    model_kwargs = dict(
        n_dipoles=16,
        dirichlet_alpha=args.dirichlet_alpha,
        interface_seed=2026,
        dtype=jnp.float32,
    )

    kinds = (
        ("amplitude", "complex")
        if args.kind == "both"
        else (args.kind,)
    )

    fits = {}
    for kind in kinds:
        print(f"\nRunning {kind} BlackJAX SMC")

        model = smc.build_model_interface(
            predict_beam_db,
            observed_db,
            observed_sigma_db,
            kind=kind,
            phase_kappa=phase_kappa if kind == "complex" else None,
            **model_kwargs,
        )

        raw_result = smc.run_smc(
            model,
            config,
            seed=args.seed,
            execution=execution,
            print_progress=not args.quiet,
        )
        fits[kind] = package_fit(
            smc,
            model,
            raw_result,
            observed_db.size,
        )
        del raw_result, model

    # ------------------------------------------------------------------
    # Pickle canonical data plus expensive inference products only.
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
            "n_comparison_pixels": int(comparison_indices.size),
            "single_pass_uncertainty_cut": single_pass_cut_metadata,
            "n_single_pass_masked": int(masked_indices.size),
            "dirichlet_alpha": args.dirichlet_alpha,
            "phase_prior_std_deg": args.phase_prior_std_deg,
            "phase_kappa": phase_kappa,
            "seed": args.seed,
            "kind": args.kind,
            "tag": args.tag,
            "n_particles": n_particles,
            "n_devices": args.n_devices,
            "particles_per_device": particles_per_device,
            "smc_config": config_dict,
            "execution_config": asdict(execution),
            "command": sys.argv,
        },
        "data": {
            "data_map": data_map,
            "error_map": error_map,
            "count_map": count_map,
            "n_pass_map": n_pass_map,
            "single_pass_robust_z": single_pass_robust_z.astype(np.float32),
            "single_pass_mask": single_pass_mask,
            "masked_indices": masked_indices.astype(np.int32),
            "comparison_indices": comparison_indices.astype(np.int32),
        },
        "nominal": {
            "log_likelihood": nominal_log_likelihood,
            "chi2": nominal_chi2,
            "reduced_chi2": nominal_chi2 / observed_db.size,
            "residual_rms_db": float(
                np.sqrt(np.mean(nominal_residual**2))
            ),
            # Nominal has no fitted beam-gain parameters.
            "bic": float(-2.0 * nominal_log_likelihood),
        },
        **fits,
    }

    with output_path.open("wb") as handle:
        pickle.dump(
            result,
            handle,
            protocol=pickle.HIGHEST_PROTOCOL,
        )

    print(f"\nSaved: {output_path}")
    print(f"Comparison pixels: {comparison_indices.size}")

    if args.single_pass_z_cut is not None:
        print(
            f"Single-pass uncertainty mask: "
            f"{masked_indices.size} pixels "
            f"(z < {args.single_pass_z_cut:g}, "
            f"sigma_cut={single_pass_cut_metadata['sigma_cut_db']:.4f} dB)"
        )

    bic_parts = [f"nominal={result['nominal']['bic']:.1f}"]
    for kind in kinds:
        bic_parts.append(f"{kind}={fits[kind]['bic']:.1f}")
    print("BIC: " + " / ".join(bic_parts))

    for kind in kinds:
        fit = fits[kind]
        final_steps = np.asarray(fit["final_step_sizes"])
        print(
            f"{kind}: beta={fit['final_beta']:.8f}, "
            f"stages={len(fit['stage_diagnostics']['stage'])}, "
            f"runtime={fit['runtime_seconds']:.1f}s, "
            f"logZ={fit['log_evidence']:.3f}, "
            f"step_med={np.median(final_steps):.3e}"
        )


if __name__ == "__main__":
    main()
