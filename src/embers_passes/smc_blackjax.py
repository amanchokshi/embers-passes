"""BlackJAX adaptive-tempered SMC + NUTS for MWA beam inference.

This module contains the inference machinery only.  Machine-specific launch
configuration (for example XLA_FLAGS / logical CPU device count), data loading,
file naming, and analysis belong in a runner script or notebook.

Core algorithm
--------------
- BlackJAX adaptive tempered SMC
- systematic resampling
- BlackJAX NUTS rejuvenation
- one-time ``find_reasonable_step_size`` at the first non-zero temperature
- BlackJAX SMC inner-kernel step-size tuning via
  ``update_scale_from_acceptance_rate``
- fixed identity inverse mass matrix
- optional multi-device ``pmap(vmap(...))`` execution of the independent NUTS
  particle mutations

The statistical SMC logic remains BlackJAX-native.  The custom code in this
module only specializes model construction, prior initialization, diagnostics,
and the device-level execution of independent NUTS transitions.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Callable, Literal, NamedTuple
import platform as _platform
import time

import blackjax
import jax
import jax.numpy as jnp
import numpy as np
from jax.scipy.special import gammaln, i0e, logsumexp

from blackjax.adaptation.step_size import find_reasonable_step_size
from blackjax.smc import ess as smc_ess
from blackjax.smc import extend_params, resampling, solver
from blackjax.smc.base import SMCInfo
from blackjax.smc.inner_kernel_tuning import as_top_level_api as inner_kernel_tuning
from blackjax.smc.tuning.from_kernel_info import update_scale_from_acceptance_rate


Array = jax.Array
ModelKind = Literal["amplitude", "complex"]


# =============================================================================
# Configuration / result containers
# =============================================================================


@dataclass(frozen=True)
class SMCConfig:
    """Scientific/numerical configuration for BlackJAX SMC.

    ``n_devices`` is deliberately *not* part of this configuration.  Device
    count is an execution concern and is supplied separately to ``run_smc``.
    """

    n_particles: int = 384
    target_ess: float = 0.80
    target_accept: float = 0.80
    max_num_doublings: int = 9
    max_stages: int = 256
    num_mcmc_steps: int = 1
    step_search_initial_scale: float = 1.0
    dtype: Any = jnp.float32
    print_every: int = 1

    def validate(self) -> None:
        if self.n_particles < 2:
            raise ValueError("n_particles must be at least two.")
        if not 0.0 < self.target_ess <= 1.0:
            raise ValueError("target_ess must lie in (0, 1].")
        if not 0.0 < self.target_accept < 1.0:
            raise ValueError("target_accept must lie in (0, 1).")
        if self.max_num_doublings < 1:
            raise ValueError("max_num_doublings must be positive.")
        if self.max_stages < 1:
            raise ValueError("max_stages must be positive.")
        if self.num_mcmc_steps != 1:
            raise ValueError(
                "The validated parallel backend currently supports exactly "
                "one NUTS transition per SMC stage."
            )
        if self.step_search_initial_scale <= 0.0:
            raise ValueError("step_search_initial_scale must be positive.")
        if self.print_every < 1:
            raise ValueError("print_every must be positive.")


@dataclass(frozen=True)
class ExecutionConfig:
    """Execution-only settings.

    The runner should arrange for JAX to expose the desired logical devices
    *before importing this module*.  For CPU runs this normally means setting
    ``XLA_FLAGS=--xla_force_host_platform_device_count=N`` before importing JAX.
    """

    n_devices: int = 1
    platform: str = "cpu"

    def validate(self, n_particles: int) -> None:
        if self.n_devices < 1:
            raise ValueError("n_devices must be positive.")
        devices = jax.devices(self.platform)
        if len(devices) < self.n_devices:
            raise ValueError(
                f"Requested {self.n_devices} {self.platform!r} devices, but JAX "
                f"exposes only {len(devices)}."
            )
        if n_particles % self.n_devices:
            raise ValueError(
                f"{n_particles} particles cannot be divided evenly across "
                f"{self.n_devices} devices."
            )


@dataclass
class SMCResult:
    """Pickle-friendly result of one complete SMC run."""

    model_kind: str
    seed: int
    particles: np.ndarray
    weights: np.ndarray
    log_prior: np.ndarray
    log_likelihood: np.ndarray
    final_beta: float
    final_step_sizes: np.ndarray
    stage_diagnostics: dict[str, np.ndarray]
    cumulative_log_evidence: np.ndarray
    log_evidence: float
    runtime_seconds: float
    config: dict[str, Any]
    execution: dict[str, Any]
    metadata: dict[str, Any]


class TunedNUTSInfo(NamedTuple):
    """NUTS info plus the step size that produced the transition."""

    nuts_info: Any
    step_size: Array


# =============================================================================
# Explicit beam prior / coordinate interface
# =============================================================================


def _alr_to_fraction(amplitude_coordinates: Array) -> Array:
    """Invert the additive-log-ratio (ALR) transform.

    For K simplex fractions ``f`` we use the final component as reference,

        z_i = log(f_i / f_{K-1}),  i=0,...,K-2.

    The inverse is ``softmax([z, 0])``.
    """

    augmented = jnp.concatenate(
        [amplitude_coordinates, jnp.zeros(1, dtype=amplitude_coordinates.dtype)]
    )
    return jax.nn.softmax(augmented)


def _fraction_to_alr(fraction: Array) -> Array:
    """Map a simplex fraction vector to additive-log-ratio coordinates."""

    return jnp.log(fraction[:-1] / fraction[-1])


def _log_dirichlet_in_alr(
    amplitude_coordinates: Array,
    concentration: Array,
) -> Array:
    """Dirichlet log density expressed in ALR coordinates.

    The ordinary Dirichlet density is defined with respect to simplex
    coordinates.  NUTS samples the unconstrained ALR coordinates, so the
    change-of-variables Jacobian is included explicitly:

        log |df/dz| = sum_i log(f_i).
    """

    fraction = _alr_to_fraction(amplitude_coordinates)
    log_fraction = jnp.log(fraction)
    log_normalizer = gammaln(jnp.sum(concentration)) - jnp.sum(
        gammaln(concentration)
    )
    log_dirichlet = log_normalizer + jnp.sum(
        (concentration - 1.0) * log_fraction
    )
    log_jacobian = jnp.sum(log_fraction)
    return log_dirichlet + log_jacobian


def _log_von_mises(relative_phase: Array, phase_kappa: Array) -> Array:
    """Independent zero-centred Von Mises log density on [-pi, pi).

    ``i0e`` is used for a stable log normalizer:

        log I0(kappa) = log(i0e(kappa)) + |kappa|.

    We deliberately keep one fundamental angular interval.  This avoids the
    infinitely repeated density that would result from extending cos(phi) over
    all of R while preserving the explicit coordinates used by the validated
    notebook implementation.
    """

    kappa = jnp.asarray(phase_kappa, dtype=relative_phase.dtype)
    log_i0 = jnp.log(i0e(kappa)) + jnp.abs(kappa)
    per_phase = (
        kappa * jnp.cos(relative_phase)
        - jnp.log(jnp.asarray(2.0 * jnp.pi, dtype=relative_phase.dtype))
        - log_i0
    )
    in_fundamental_interval = jnp.all(
        (relative_phase >= -jnp.pi) & (relative_phase < jnp.pi)
    )
    return jnp.where(
        in_fundamental_interval,
        jnp.sum(per_phase),
        -jnp.inf,
    )


def _sample_dirichlet(
    rng_key: Array,
    concentration: Array,
    n_samples: int,
    dtype,
) -> Array:
    """Sample a Dirichlet exactly by normalizing independent Gamma draws."""

    gamma = jax.random.gamma(
        rng_key,
        concentration,
        shape=(n_samples, concentration.shape[0]),
        dtype=dtype,
    )
    return gamma / jnp.sum(gamma, axis=-1, keepdims=True)


def _sample_von_mises_one(rng_key: Array, kappa: Array, dtype) -> Array:
    """Sample Von Mises(0, kappa) with the Best--Fisher rejection method."""

    kappa = jnp.asarray(kappa, dtype=dtype)

    def sample_concentrated(key):
        a = 1.0 + jnp.sqrt(1.0 + 4.0 * kappa**2)
        b = (a - jnp.sqrt(2.0 * a)) / (2.0 * kappa)
        r = (1.0 + b**2) / (2.0 * b)

        def cond_fn(carry):
            _, _, accepted = carry
            return ~accepted

        def body_fn(carry):
            key, _, _ = carry
            key, draw_key = jax.random.split(key)
            u = jax.random.uniform(
                draw_key,
                shape=(3,),
                minval=jnp.asarray(0.0, dtype=dtype),
                maxval=jnp.asarray(1.0, dtype=dtype),
                dtype=dtype,
            )
            z = jnp.cos(jnp.pi * u[0])
            f = (1.0 + r * z) / (r + z)
            c = kappa * (r - f)
            accept = (u[1] < c * (2.0 - c)) | (
                jnp.log(c / u[1]) + 1.0 - c >= 0.0
            )
            theta = jnp.where(
                u[2] > 0.5,
                jnp.arccos(f),
                -jnp.arccos(f),
            )
            return key, theta, accept

        _, theta, _ = jax.lax.while_loop(
            cond_fn,
            body_fn,
            (key, jnp.asarray(0.0, dtype=dtype), jnp.asarray(False)),
        )
        return theta

    def sample_uniform(key):
        return jax.random.uniform(
            key,
            shape=(),
            minval=-jnp.pi,
            maxval=jnp.pi,
            dtype=dtype,
        )

    return jax.lax.cond(kappa < 1e-6, sample_uniform, sample_concentrated, rng_key)


def _sample_relative_phases(
    rng_key: Array,
    n_samples: int,
    n_phase_coords: int,
    phase_kappa: float,
    dtype,
) -> Array:
    """Draw independent zero-centred Von Mises phase coordinates."""

    keys = jax.random.split(rng_key, n_samples * n_phase_coords)
    draws = jax.vmap(
        lambda key: _sample_von_mises_one(key, phase_kappa, dtype)
    )(keys)
    return draws.reshape(n_samples, n_phase_coords)


def build_model_interface(
    predict_beam_db: Callable[[Array], Array],
    observed_db,
    observed_sigma_db,
    *,
    kind: ModelKind,
    n_dipoles: int = 16,
    dirichlet_alpha: float = 10.0,
    phase_kappa: float | None = None,
    interface_seed: int | None = None,
    dtype=jnp.float32,
) -> dict[str, Any]:
    """Build the explicit BlackJAX target used by the validated notebook.

    Coordinates
    -----------
    ``amplitude``
        ``n_dipoles - 1`` ALR coordinates for a Dirichlet-distributed gain
        fraction vector.

    ``complex``
        The same ALR amplitude coordinates followed by ``n_dipoles - 1``
        direct relative phases in the fundamental interval [-pi, pi).

    ``interface_seed`` is retained only for backwards-compatible runner calls;
    the explicit interface itself has no random initialization.
    """

    del interface_seed

    if n_dipoles < 2:
        raise ValueError("n_dipoles must be at least two.")
    if dirichlet_alpha <= 0.0:
        raise ValueError("dirichlet_alpha must be positive.")
    if kind not in ("amplitude", "complex"):
        raise ValueError(f"Unknown model kind: {kind!r}.")
    if kind == "complex" and (phase_kappa is None or phase_kappa <= 0.0):
        raise ValueError("A positive phase_kappa is required for the complex model.")

    observed_db = jnp.asarray(observed_db, dtype=dtype)
    observed_sigma_db = jnp.asarray(observed_sigma_db, dtype=dtype)
    if observed_db.ndim != 1 or observed_sigma_db.shape != observed_db.shape:
        raise ValueError("observed_db and observed_sigma_db must be matching 1-D arrays.")
    if not bool(jnp.all(jnp.isfinite(observed_db))):
        raise ValueError("observed_db contains non-finite values.")
    if not bool(jnp.all(jnp.isfinite(observed_sigma_db) & (observed_sigma_db > 0))):
        raise ValueError("observed_sigma_db must be finite and strictly positive.")

    n_amp_coords = n_dipoles - 1
    n_phase_coords = n_dipoles - 1 if kind == "complex" else 0
    latent_dim = n_amp_coords + n_phase_coords
    concentration = jnp.full(n_dipoles, dirichlet_alpha, dtype=dtype)
    phase_kappa_array = (
        jnp.asarray(phase_kappa, dtype=dtype) if kind == "complex" else None
    )

    @jax.jit
    def position_to_physical(position):
        amplitude_coordinates = position[:n_amp_coords]
        fraction = _alr_to_fraction(amplitude_coordinates)
        amplitude = n_dipoles * fraction

        if kind == "complex":
            relative_phase = position[n_amp_coords:]
            phase = jnp.concatenate(
                [jnp.zeros(1, dtype=dtype), relative_phase]
            )
            gain = (amplitude * jnp.exp(1j * phase)).astype(jnp.complex64)
        else:
            relative_phase = jnp.zeros(0, dtype=dtype)
            phase = jnp.zeros(n_dipoles, dtype=dtype)
            gain = amplitude.astype(jnp.complex64)

        return fraction, amplitude, relative_phase, phase, gain

    @jax.jit
    def log_prior_flat(position):
        amplitude_coordinates = position[:n_amp_coords]
        logp = _log_dirichlet_in_alr(amplitude_coordinates, concentration)
        if kind == "complex":
            relative_phase = position[n_amp_coords:]
            logp = logp + _log_von_mises(relative_phase, phase_kappa_array)
        return logp

    @jax.jit
    def log_likelihood_flat(position):
        _, _, _, _, gains = position_to_physical(position)
        model_db = predict_beam_db(gains)
        residual = (observed_db - model_db) / observed_sigma_db
        return -0.5 * jnp.sum(
            residual**2
            + jnp.log(
                jnp.asarray(2.0 * jnp.pi, dtype=dtype)
                * observed_sigma_db**2
            )
        )

    @jax.jit
    def log_posterior_flat(position, beta=1.0):
        return log_prior_flat(position) + beta * log_likelihood_flat(position)

    @jax.jit
    def physical_sites_from_flat(position):
        fraction, amplitude, relative_phase, phase, gain = position_to_physical(
            position
        )
        sites = {
            "dipole_gain_fraction": fraction,
            "dipole_gain_amplitude": amplitude,
            "dipole_gain_phase": phase,
            "dipole_gain_complex": gain,
        }
        if kind == "complex":
            sites["relative_phase"] = relative_phase
        return sites

    return {
        "kind": kind,
        "n_dipoles": n_dipoles,
        "n_amplitude_coordinates": n_amp_coords,
        "n_phase_coordinates": n_phase_coords,
        "latent_dim": latent_dim,
        "dirichlet_alpha": float(dirichlet_alpha),
        "dirichlet_concentration": concentration,
        "phase_kappa": None if phase_kappa is None else float(phase_kappa),
        "log_prior_flat": log_prior_flat,
        "log_likelihood_flat": log_likelihood_flat,
        "log_posterior_flat": log_posterior_flat,
        "position_to_physical": position_to_physical,
        "physical_sites_from_flat": physical_sites_from_flat,
        "log_prior_particles": jax.jit(jax.vmap(log_prior_flat, in_axes=0)),
        "log_likelihood_particles": jax.jit(
            jax.vmap(log_likelihood_flat, in_axes=0)
        ),
        "physical_sites": jax.jit(jax.vmap(physical_sites_from_flat, in_axes=0)),
    }


# =============================================================================
# Exact prior initialization
# =============================================================================


def sample_prior_positions(
    model: dict[str, Any],
    rng_key: Array,
    n_particles: int,
    *,
    dtype=jnp.float32,
) -> Array:
    """Draw exact physical-prior samples in the explicit BlackJAX coordinates."""

    amplitude_key, phase_key = jax.random.split(rng_key)
    concentration = jnp.asarray(model["dirichlet_concentration"], dtype=dtype)
    fractions = _sample_dirichlet(
        amplitude_key, concentration, n_particles, dtype
    )
    amplitude_coordinates = jax.vmap(_fraction_to_alr)(fractions)

    if model["kind"] == "complex":
        relative_phase = _sample_relative_phases(
            phase_key,
            n_particles,
            int(model["n_phase_coordinates"]),
            float(model["phase_kappa"]),
            dtype,
        )
        positions = jnp.concatenate(
            [amplitude_coordinates, relative_phase],
            axis=1,
        )
    else:
        positions = amplitude_coordinates

    positions = jnp.asarray(positions, dtype=dtype)
    expected_shape = (n_particles, int(model["latent_dim"]))
    if positions.shape != expected_shape:
        raise RuntimeError(
            f"Unexpected prior position shape {positions.shape}; expected {expected_shape}."
        )
    if not bool(jnp.all(jnp.isfinite(positions))):
        raise FloatingPointError("Initial SMC population contains non-finite values.")
    return positions


# =============================================================================
# Parallel mutation backend
# =============================================================================


def _parallel_shard(array: Array, n_devices: int) -> Array:
    array = jnp.asarray(array)
    n_particles = array.shape[0]
    if n_particles % n_devices:
        raise ValueError(
            f"{n_particles} particles cannot be divided evenly across "
            f"{n_devices} devices."
        )
    per_device = n_particles // n_devices
    return array.reshape((n_devices, per_device) + array.shape[1:])


def _parallel_flatten_pytree(pytree):
    return jax.tree_util.tree_map(
        lambda x: x.reshape((x.shape[0] * x.shape[1],) + x.shape[2:]), pytree
    )


def _expand_particle_parameter(parameter: Array, n_particles: int) -> Array:
    """Expand BlackJAX shared/per-particle parameter conventions explicitly."""

    # Explicit host round-trip removes committed pmap sharding inherited from
    # the previous stage.  This is intentional and mirrors the validated
    # notebook implementation.
    parameter = jnp.asarray(jax.device_get(parameter))

    if parameter.shape[0] == 1:
        parameter = jnp.broadcast_to(
            parameter, (n_particles,) + parameter.shape[1:]
        )
    elif parameter.shape[0] != n_particles:
        raise ValueError(
            f"Unexpected parameter leading dimension {parameter.shape[0]}; "
            f"expected 1 or {n_particles}."
        )
    return jnp.asarray(jax.device_get(parameter))


def build_parallel_update_particles_fn(
    model: dict[str, Any],
    config: SMCConfig,
    execution: ExecutionConfig,
    *,
    nuts_kernel,
):
    """Build BlackJAX ``update_particles_fn`` with pmap(vmap(NUTS))."""

    config.validate()
    execution.validate(config.n_particles)

    devices = jax.devices(execution.platform)[: execution.n_devices]
    log_prior_fn = model["log_prior_flat"]
    log_likelihood_fn = model["log_likelihood_flat"]

    def one_particle(
        rng_key,
        position,
        beta,
        step_size,
        inverse_mass_matrix,
    ):
        def logdensity_fn(x):
            return log_prior_fn(x) + beta * log_likelihood_fn(x)

        state = blackjax.nuts.init(position, logdensity_fn)
        new_state, nuts_info = nuts_kernel(
            rng_key,
            state,
            logdensity_fn,
            step_size,
            inverse_mass_matrix,
            max_num_doublings=config.max_num_doublings,
        )
        return (
            new_state.position,
            TunedNUTSInfo(nuts_info=nuts_info, step_size=step_size),
        )

    population_kernel = jax.pmap(
        jax.vmap(
            one_particle,
            in_axes=(0, 0, None, 0, 0),
        ),
        in_axes=(0, 0, None, 0, 0),
        devices=devices,
    )

    def update_particles(
        rng_key,
        state,
        num_mcmc_steps,
        mcmc_parameters,
        tempered_logposterior_fn,
        log_weights_fn,
    ):
        if int(num_mcmc_steps) != 1:
            raise NotImplementedError(
                "The validated parallel backend currently supports "
                "num_mcmc_steps == 1 only."
            )
        del tempered_logposterior_fn

        updating_key, resampling_key = jax.random.split(rng_key, 2)
        n_particles = int(state.weights.shape[0])
        if n_particles % execution.n_devices:
            raise ValueError(
                f"{n_particles} particles must be divisible by "
                f"{execution.n_devices} devices."
            )

        ancestors = resampling.systematic(
            resampling_key,
            state.weights,
            n_particles,
        )
        resampled_particles = state.particles[ancestors]

        step_sizes = _expand_particle_parameter(
            mcmc_parameters["step_size"], n_particles
        )
        inverse_mass_matrices = _expand_particle_parameter(
            mcmc_parameters["inverse_mass_matrix"], n_particles
        )
        particle_keys = jax.random.split(updating_key, n_particles)

        particles_sharded = _parallel_shard(
            resampled_particles, execution.n_devices
        )
        keys_sharded = _parallel_shard(particle_keys, execution.n_devices)
        steps_sharded = _parallel_shard(step_sizes, execution.n_devices)
        mass_sharded = _parallel_shard(
            inverse_mass_matrices, execution.n_devices
        )

        particles_sharded_new, info_sharded = population_kernel(
            keys_sharded,
            particles_sharded,
            state.tempering_param,
            steps_sharded,
            mass_sharded,
        )

        new_particles_sharded = _parallel_flatten_pytree(particles_sharded_new)
        update_info_sharded = _parallel_flatten_pytree(info_sharded)

        # Explicitly cross the pmap boundary.  Without this, the returned arrays
        # remain committed to the device mesh and the following SMC stage can
        # fail when the same values are reshaped / remapped.
        new_particles = jnp.asarray(
            jax.device_get(new_particles_sharded),
            dtype=state.particles.dtype,
        )
        update_info = jax.tree_util.tree_map(
            lambda x: jnp.asarray(jax.device_get(x)),
            update_info_sharded,
        )

        log_weights = jax.vmap(log_weights_fn)(new_particles)
        logsum_weights = logsumexp(log_weights)
        normalizing_constant = logsum_weights - jnp.log(n_particles)
        new_weights = jnp.exp(log_weights - logsum_weights)

        new_state = state._replace(
            particles=new_particles,
            weights=new_weights,
        )
        info = SMCInfo(
            ancestors,
            normalizing_constant,
            update_info,
        )
        return new_state, info

    return update_particles


# =============================================================================
# Sampler construction / initialization
# =============================================================================


def _compute_next_beta(state, log_likelihood_fn, target_ess: float):
    beta_previous = state.tempering_param
    max_delta = 1.0 - beta_previous
    delta = smc_ess.ess_solver(
        jax.vmap(log_likelihood_fn),
        state.particles,
        target_ess,
        max_delta,
        solver.dichotomy,
    )
    delta = jnp.clip(delta, 0.0, max_delta)
    return beta_previous + delta


def _find_initial_step_size(
    model: dict[str, Any],
    particles: Array,
    first_beta: Array,
    inverse_mass_matrix: Array,
    rng_key: Array,
    config: SMCConfig,
    *,
    nuts_kernel,
) -> tuple[Array, int]:
    """Run BlackJAX's one-time reasonable-step-size search."""

    log_likelihoods = model["log_likelihood_particles"](particles)
    order = jnp.argsort(log_likelihoods)
    reference_index = int(order[particles.shape[0] // 2])
    reference_position = particles[reference_index]

    def logdensity_fn(position):
        return (
            model["log_prior_flat"](position)
            + first_beta * model["log_likelihood_flat"](position)
        )

    reference_state = blackjax.nuts.init(reference_position, logdensity_fn)

    def kernel_generator(step_size):
        def kernel(key, state):
            return nuts_kernel(
                key,
                state,
                logdensity_fn,
                step_size,
                inverse_mass_matrix,
                max_num_doublings=config.max_num_doublings,
            )

        return kernel

    step_size = find_reasonable_step_size(
        rng_key,
        kernel_generator,
        reference_state,
        jnp.asarray(config.step_search_initial_scale, dtype=config.dtype),
        target_accept=config.target_accept,
    )
    return jnp.asarray(step_size, dtype=config.dtype), reference_index


def build_sampler(
    model: dict[str, Any],
    config: SMCConfig,
    execution: ExecutionConfig,
):
    """Build the tuned BlackJAX SMC algorithm for a model/config pair."""

    config.validate()
    execution.validate(config.n_particles)

    latent_dim = int(model["latent_dim"])
    inverse_mass_matrix = jnp.ones(latent_dim, dtype=config.dtype)
    nuts_kernel = blackjax.nuts.build_kernel()

    def tuned_nuts_step(
        rng_key,
        state,
        logdensity_fn,
        step_size,
        inverse_mass_matrix,
    ):
        new_state, nuts_info = nuts_kernel(
            rng_key,
            state,
            logdensity_fn,
            step_size,
            inverse_mass_matrix,
            max_num_doublings=config.max_num_doublings,
        )
        return new_state, TunedNUTSInfo(nuts_info, step_size)

    parallel_update_particles = build_parallel_update_particles_fn(
        model,
        config,
        execution,
        nuts_kernel=nuts_kernel,
    )

    return {
        "inverse_mass_matrix": inverse_mass_matrix,
        "nuts_kernel": nuts_kernel,
        "tuned_nuts_step": tuned_nuts_step,
        "parallel_update_particles": parallel_update_particles,
    }


# =============================================================================
# Full run
# =============================================================================


def run_smc(
    model: dict[str, Any],
    config: SMCConfig,
    *,
    seed: int,
    execution: ExecutionConfig | None = None,
    print_progress: bool = True,
) -> SMCResult:
    """Run BlackJAX adaptive-tempered SMC from the prior to beta=1."""

    config.validate()
    if execution is None:
        execution = ExecutionConfig()
    execution.validate(config.n_particles)

    backend = build_sampler(model, config, execution)
    inverse_mass_matrix = backend["inverse_mass_matrix"]
    nuts_kernel = backend["nuts_kernel"]
    tuned_nuts_step = backend["tuned_nuts_step"]
    parallel_update_particles = backend["parallel_update_particles"]

    # One master key for a reproducible run; split into independent logical
    # substreams for prior initialization, step-size search, and SMC stages.
    master_key = jax.random.key(seed)
    prior_key, stepfind_key, smc_key = jax.random.split(master_key, 3)

    particles = sample_prior_positions(
        model,
        prior_key,
        config.n_particles,
        dtype=config.dtype,
    )

    log_likelihoods_initial = model["log_likelihood_particles"](particles)
    initial_state = blackjax.adaptive_tempered_smc.init(particles)
    first_beta = _compute_next_beta(
        initial_state,
        model["log_likelihood_flat"],
        config.target_ess,
    )

    initial_step_size, reference_index = _find_initial_step_size(
        model,
        particles,
        first_beta,
        inverse_mass_matrix,
        stepfind_key,
        config,
        nuts_kernel=nuts_kernel,
    )

    n_particles = config.n_particles

    def mcmc_parameter_update_fn(rng_key, state, info):
        del rng_key, state
        acceptance_rates = jnp.ravel(info.update_info.nuts_info.acceptance_rate)
        used_step_sizes = jnp.ravel(info.update_info.step_size)
        # One transition per particle in the validated implementation.
        acceptance_rates = acceptance_rates[:n_particles]
        used_step_sizes = used_step_sizes[:n_particles]
        next_step_sizes = update_scale_from_acceptance_rate(
            used_step_sizes,
            acceptance_rates,
            target_acceptance_rate=config.target_accept,
        )
        return {
            "step_size": next_step_sizes,
            "inverse_mass_matrix": inverse_mass_matrix[None, :],
        }

    initial_parameters = extend_params(
        {
            "step_size": initial_step_size,
            "inverse_mass_matrix": inverse_mass_matrix,
        }
    )

    tuned_smc = inner_kernel_tuning(
        logprior_fn=model["log_prior_flat"],
        loglikelihood_fn=model["log_likelihood_flat"],
        mcmc_step_fn=tuned_nuts_step,
        mcmc_init_fn=blackjax.nuts.init,
        resampling_fn=resampling.systematic,
        smc_algorithm=blackjax.adaptive_tempered_smc,
        mcmc_parameter_update_fn=mcmc_parameter_update_fn,
        initial_parameter_value=initial_parameters,
        target_ess=config.target_ess,
        root_solver=solver.dichotomy,
        num_mcmc_steps=config.num_mcmc_steps,
        update_particles_fn=parallel_update_particles,
    )

    state = tuned_smc.init(particles)

    history: dict[str, list[Any]] = {
        "stage": [],
        "beta_previous": [],
        "beta": [],
        "delta_beta": [],
        "step_used_median": [],
        "step_used_min": [],
        "step_used_max": [],
        "step_next_median": [],
        "step_next_min": [],
        "step_next_max": [],
        "mean_acceptance": [],
        "min_acceptance": [],
        "median_integration_steps": [],
        "max_integration_steps": [],
        "ceiling_fraction": [],
        "divergence_fraction": [],
        "n_unique_ancestors": [],
        "log_evidence_increment": [],
        "runtime_seconds": [],
    }

    if print_progress:
        print(
            f"Initial SMC population: {particles.shape}\n"
            f"Initial logL range: {float(jnp.min(log_likelihoods_initial)):.3f} "
            f"to {float(jnp.max(log_likelihoods_initial)):.3f}\n"
            f"Initial beta: {float(initial_state.tempering_param):.8f}\n"
            f"First ESS-selected beta: {float(first_beta):.6e}\n"
            f"Reference particle index: {reference_index}\n"
            f"Reference particle logL: "
            f"{float(log_likelihoods_initial[reference_index]):.3f}\n"
            f"BlackJAX reasonable initial step size: "
            f"{float(initial_step_size):.6e}\n"
            f"Devices: {execution.n_devices} x {execution.platform}"
        )
        print("\nBlackJAX adaptive tempered SMC + NUTS + inner-kernel tuning")
        print("---------------------------------------------------------")

    max_steps = 2**config.max_num_doublings - 1
    cumulative_log_evidence = [0.0]
    run_start = time.perf_counter()

    for stage_index in range(1, config.max_stages + 1):
        sampler_state = state.sampler_state
        beta_previous = float(sampler_state.tempering_param)
        if beta_previous >= 1.0 - 1e-7:
            break

        used_steps = jnp.ravel(state.parameter_override["step_size"])
        if used_steps.size == 1:
            used_steps_diag = jnp.repeat(used_steps, n_particles)
        else:
            used_steps_diag = used_steps

        smc_key, stage_key = jax.random.split(smc_key)
        stage_start = time.perf_counter()
        state, smc_info = tuned_smc.step(stage_key, state)
        jax.block_until_ready(state.sampler_state.particles)
        stage_runtime = time.perf_counter() - stage_start

        beta = float(state.sampler_state.tempering_param)
        delta_beta = beta - beta_previous

        nuts_info = smc_info.update_info.nuts_info
        acceptance = jnp.ravel(nuts_info.acceptance_rate)
        integration_steps = jnp.ravel(nuts_info.num_integration_steps)
        divergent = jnp.ravel(nuts_info.is_divergent)
        next_steps = jnp.ravel(state.parameter_override["step_size"])

        ceiling_fraction = jnp.mean(integration_steps >= max_steps)
        divergence_fraction = jnp.mean(divergent.astype(jnp.float32))
        unique_ancestors = jnp.unique(
            smc_info.ancestors,
            size=n_particles,
            fill_value=-1,
        )
        n_unique_ancestors = jnp.sum(unique_ancestors >= 0)
        log_z_increment = float(smc_info.log_likelihood_increment)
        cumulative_log_evidence.append(
            cumulative_log_evidence[-1] + log_z_increment
        )

        history["stage"].append(stage_index)
        history["beta_previous"].append(beta_previous)
        history["beta"].append(beta)
        history["delta_beta"].append(delta_beta)
        history["step_used_median"].append(float(jnp.median(used_steps_diag)))
        history["step_used_min"].append(float(jnp.min(used_steps_diag)))
        history["step_used_max"].append(float(jnp.max(used_steps_diag)))
        history["step_next_median"].append(float(jnp.median(next_steps)))
        history["step_next_min"].append(float(jnp.min(next_steps)))
        history["step_next_max"].append(float(jnp.max(next_steps)))
        history["mean_acceptance"].append(float(jnp.mean(acceptance)))
        history["min_acceptance"].append(float(jnp.min(acceptance)))
        history["median_integration_steps"].append(
            float(jnp.median(integration_steps))
        )
        history["max_integration_steps"].append(int(jnp.max(integration_steps)))
        history["ceiling_fraction"].append(float(ceiling_fraction))
        history["divergence_fraction"].append(float(divergence_fraction))
        history["n_unique_ancestors"].append(int(n_unique_ancestors))
        history["log_evidence_increment"].append(log_z_increment)
        history["runtime_seconds"].append(stage_runtime)

        if print_progress and (
            stage_index % config.print_every == 0 or beta >= 1.0 - 1e-7
        ):
            print(
                f"stage={stage_index:03d} "
                f"beta={beta:.6e} "
                f"dBeta={delta_beta:.2e} "
                f"step={float(jnp.median(used_steps_diag)):.3e}"
                f"[{float(jnp.min(used_steps_diag)):.2e},"
                f"{float(jnp.max(used_steps_diag)):.2e}] -> "
                f"{float(jnp.median(next_steps)):.3e} "
                f"acc={float(jnp.mean(acceptance)):.3f} "
                f"steps={float(jnp.median(integration_steps)):.0f}/"
                f"{int(jnp.max(integration_steps))} "
                f"ceiling={100.0 * float(ceiling_fraction):4.1f}% "
                f"div={int(jnp.sum(divergent))} "
                f"anc={int(n_unique_ancestors)}/{n_particles} "
                f"time={stage_runtime:.2f}s"
            )
    else:
        raise RuntimeError(
            f"Did not reach beta=1 within max_stages={config.max_stages}."
        )

    runtime_seconds = time.perf_counter() - run_start
    final_sampler_state = state.sampler_state
    final_particles = final_sampler_state.particles
    final_weights = final_sampler_state.weights
    final_step_sizes = jnp.ravel(state.parameter_override["step_size"])
    final_log_prior = model["log_prior_particles"](final_particles)
    final_log_likelihood = model["log_likelihood_particles"](final_particles)
    final_log_likelihood.block_until_ready()

    history_np = {
        name: np.asarray(values)
        for name, values in history.items()
    }

    metadata = {
        "jax_version": getattr(jax, "__version__", "unknown"),
        "blackjax_version": getattr(blackjax, "__version__", "unknown"),
        "python_version": _platform.python_version(),
        "platform": _platform.platform(),
        "jax_devices": [str(d) for d in jax.devices(execution.platform)],
        "latent_dim": int(model["latent_dim"]),
        "max_integration_steps": int(max_steps),
        "reference_particle_index": int(reference_index),
        "initial_step_size": float(initial_step_size),
        "first_beta": float(first_beta),
    }

    return SMCResult(
        model_kind=str(model["kind"]),
        seed=int(seed),
        particles=np.asarray(final_particles),
        weights=np.asarray(final_weights),
        log_prior=np.asarray(final_log_prior),
        log_likelihood=np.asarray(final_log_likelihood),
        final_beta=float(final_sampler_state.tempering_param),
        final_step_sizes=np.asarray(final_step_sizes),
        stage_diagnostics=history_np,
        cumulative_log_evidence=np.asarray(cumulative_log_evidence),
        log_evidence=float(cumulative_log_evidence[-1]),
        runtime_seconds=float(runtime_seconds),
        config=_jsonable_dataclass_dict(config),
        execution=asdict(execution),
        metadata=metadata,
    )


# =============================================================================
# Analysis / serialization helpers
# =============================================================================


def physical_sites(model: dict[str, Any], particles: Array):
    """Transform explicit BlackJAX coordinates back to physical gain sites."""

    return model["physical_sites"](particles)


def physical_sites_to_numpy(sites) -> dict[str, np.ndarray]:
    return {name: np.asarray(value) for name, value in sites.items()}


def result_to_dict(result: SMCResult) -> dict[str, Any]:
    """Return a plain pickle-friendly dictionary."""

    return {
        "model_kind": result.model_kind,
        "seed": result.seed,
        "particles": result.particles,
        "weights": result.weights,
        "log_prior": result.log_prior,
        "log_likelihood": result.log_likelihood,
        "final_beta": result.final_beta,
        "final_step_sizes": result.final_step_sizes,
        "stage_diagnostics": result.stage_diagnostics,
        "cumulative_log_evidence": result.cumulative_log_evidence,
        "log_evidence": result.log_evidence,
        "runtime_seconds": result.runtime_seconds,
        "config": result.config,
        "execution": result.execution,
        "metadata": result.metadata,
    }


def _jsonable_dataclass_dict(obj) -> dict[str, Any]:
    raw = asdict(obj)
    output = {}
    for key, value in raw.items():
        if value in (jnp.float32, np.float32):
            output[key] = "float32"
        elif value in (jnp.float64, np.float64):
            output[key] = "float64"
        else:
            output[key] = value
    return output


__all__ = [
    "ExecutionConfig",
    "ModelKind",
    "SMCConfig",
    "SMCResult",
    "build_model_interface",
    "build_parallel_update_particles_fn",
    "build_sampler",
    "physical_sites",
    "physical_sites_to_numpy",
    "result_to_dict",
    "run_smc",
    "sample_prior_positions",
]
