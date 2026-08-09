"""Minimal annealed SMC utilities for MWA beam inference.

The module deliberately keeps the implementation procedural.  The only
structured types are the numerical configuration, particle-population state,
and per-stage diagnostics.  File I/O, HEALPix handling, posterior beam
evaluation, plotting, mode finding, and model comparison belong in the CLI or
analysis notebook.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Literal, NamedTuple

import blackjax
import jax
import jax.numpy as jnp
import numpy as np
import numpyro
import numpyro.distributions as dist
from jax.flatten_util import ravel_pytree
from numpyro.handlers import reparam
from numpyro.infer.initialization import init_to_sample
from numpyro.infer.reparam import CircularReparam
from numpyro.infer.util import initialize_model

Array = jax.Array
ModelKind = Literal["amplitude", "complex"]


# -----------------------------------------------------------------------------
# Small structured containers
# -----------------------------------------------------------------------------


@dataclass(frozen=True)
class SMCConfig:
    """Numerical settings for annealed SMC."""

    n_particles: int = 2304
    n_devices: int = 12
    ess_target_fraction: float = 0.8
    beta_tolerance: float = 1e-10
    max_bisection_steps: int = 80
    max_stages: int = 256
    min_beta_increment: float = 1e-12
    print_every: int = 1

    nuts_initial_step_size: float = 0.20
    nuts_max_num_doublings: int = 7
    nuts_divergence_threshold: float = 1_000.0
    nuts_target_acceptance_low: float = 0.75
    nuts_target_acceptance_high: float = 0.95
    nuts_step_size_decrease_factor: float = 0.5
    nuts_step_size_increase_factor: float = 1.25
    nuts_min_step_size: float = 1e-5
    nuts_max_step_size: float = 0.50
    nuts_max_mutation_attempts: int = 8
    nuts_max_trial_divergence_fraction: float = 0.0

    n_final_mutation_steps: int = 3
    dtype: object = jnp.float32

    def validate(self) -> None:
        if self.n_particles < 2:
            raise ValueError("n_particles must be at least two.")
        if self.n_devices < 1:
            raise ValueError("n_devices must be positive.")
        if self.n_particles % self.n_devices != 0:
            raise ValueError(
                f"{self.n_particles} particles cannot be divided evenly "
                f"across {self.n_devices} devices."
            )
        if not 0.0 < self.ess_target_fraction <= 1.0:
            raise ValueError("ess_target_fraction must lie in (0, 1].")
        if self.nuts_initial_step_size <= 0.0:
            raise ValueError("nuts_initial_step_size must be positive.")
        if self.n_final_mutation_steps < 0:
            raise ValueError("n_final_mutation_steps must be non-negative.")


class SMCState(NamedTuple):
    """State of the annealed SMC particle population."""

    particles: Array
    log_prior: Array
    log_likelihood: Array
    log_weights: Array
    beta: Array
    ancestor_indices: Array


class SMCStageDiagnostics(NamedTuple):
    """Diagnostics recorded for one tempering/resampling/mutation stage."""

    beta_previous: Array
    beta_next: Array
    delta_beta: Array
    ess_before_resampling: Array
    log_evidence_increment: Array
    n_unique_ancestors: Array
    nuts_step_size: Array
    mutation_attempts: Array
    mean_acceptance_rate: Array
    minimum_acceptance_rate: Array
    n_divergent: Array
    median_integration_steps: Array
    maximum_integration_steps: Array
    median_displacement: Array


# -----------------------------------------------------------------------------
# NumPyro beam models
# -----------------------------------------------------------------------------


def build_amplitude_models(
    predict_beam_db: Callable[[Array], Array],
    *,
    n_dipoles: int = 16,
    dirichlet_alpha: float = 10.0,
    dtype=jnp.float32,
) -> tuple[Callable, Callable]:
    """Return prior-only and full amplitude-only NumPyro models."""

    concentration = jnp.full(n_dipoles, dirichlet_alpha, dtype=dtype)

    def sample_excitations():
        fraction = numpyro.sample("excitation_fraction", dist.Dirichlet(concentration))
        amplitude = n_dipoles * fraction
        phase = jnp.zeros(n_dipoles, dtype=dtype)
        excitation = amplitude.astype(jnp.complex64)
        numpyro.deterministic("excitation_amplitude", amplitude)
        numpyro.deterministic("excitation_phase", phase)
        numpyro.deterministic("excitation_complex", excitation)
        return excitation

    def prior_model():
        sample_excitations()

    def full_model(observed_db, observed_sigma_db):
        excitation = sample_excitations()
        model_db = predict_beam_db(excitation)
        numpyro.deterministic("model_db", model_db)
        numpyro.sample(
            "data",
            dist.Normal(model_db, observed_sigma_db).to_event(1),
            obs=observed_db,
        )

    return prior_model, full_model


def build_complex_models(
    predict_beam_db: Callable[[Array], Array],
    *,
    n_dipoles: int = 16,
    dirichlet_alpha: float = 10.0,
    phase_kappa: float,
    dtype=jnp.float32,
) -> tuple[Callable, Callable]:
    """Return prior-only and full amplitude+phase NumPyro models."""

    concentration = jnp.full(n_dipoles, dirichlet_alpha, dtype=dtype)
    phase_loc = jnp.zeros(n_dipoles - 1, dtype=dtype)
    phase_concentration = jnp.full(n_dipoles - 1, phase_kappa, dtype=dtype)

    def sample_excitations():
        fraction = numpyro.sample("excitation_fraction", dist.Dirichlet(concentration))
        relative_phase = numpyro.sample(
            "relative_phase",
            dist.VonMises(phase_loc, phase_concentration).to_event(1),
        )
        amplitude = n_dipoles * fraction
        phase = jnp.concatenate([jnp.zeros(1, dtype=dtype), relative_phase])
        excitation = (amplitude * jnp.exp(1j * phase)).astype(jnp.complex64)
        numpyro.deterministic("excitation_amplitude", amplitude)
        numpyro.deterministic("excitation_phase", phase)
        numpyro.deterministic("excitation_complex", excitation)
        return excitation

    def prior_model_base():
        sample_excitations()

    def full_model_base(observed_db, observed_sigma_db):
        excitation = sample_excitations()
        model_db = predict_beam_db(excitation)
        numpyro.deterministic("model_db", model_db)
        numpyro.sample(
            "data",
            dist.Normal(model_db, observed_sigma_db).to_event(1),
            obs=observed_db,
        )

    config = {"relative_phase": CircularReparam()}
    return (
        reparam(prior_model_base, config=config),
        reparam(full_model_base, config=config),
    )


def build_model_interface(
    predict_beam_db: Callable[[Array], Array],
    observed_db,
    observed_sigma_db,
    *,
    kind: ModelKind,
    n_dipoles: int = 16,
    dirichlet_alpha: float = 10.0,
    phase_kappa: float | None = None,
    interface_seed: int = 2026,
    dtype=jnp.float32,
) -> dict:
    """Build the flat unconstrained functions needed by SMC.

    A plain dictionary is returned intentionally: there is no stateful model
    object, only a small bundle of functions and metadata.
    """

    observed_db = jnp.asarray(observed_db, dtype=dtype)
    observed_sigma_db = jnp.asarray(observed_sigma_db, dtype=dtype)
    if observed_db.ndim != 1 or observed_sigma_db.shape != observed_db.shape:
        raise ValueError("observed_db and observed_sigma_db must be matching 1-D arrays.")
    if not bool(jnp.all(jnp.isfinite(observed_db))):
        raise ValueError("observed_db contains non-finite values.")
    if not bool(jnp.all(jnp.isfinite(observed_sigma_db) & (observed_sigma_db > 0))):
        raise ValueError("observed_sigma_db must be finite and strictly positive.")

    if kind == "amplitude":
        prior_model, full_model = build_amplitude_models(
            predict_beam_db,
            n_dipoles=n_dipoles,
            dirichlet_alpha=dirichlet_alpha,
            dtype=dtype,
        )
    elif kind == "complex":
        if phase_kappa is None:
            raise ValueError("phase_kappa is required for the complex model.")
        prior_model, full_model = build_complex_models(
            predict_beam_db,
            n_dipoles=n_dipoles,
            dirichlet_alpha=dirichlet_alpha,
            phase_kappa=phase_kappa,
            dtype=dtype,
        )
    else:
        raise ValueError(f"Unknown model kind: {kind!r}.")

    key = jax.random.key(interface_seed)
    prior_info = initialize_model(key, prior_model, init_strategy=init_to_sample())
    full_info = initialize_model(
        key,
        full_model,
        init_strategy=init_to_sample(),
        model_kwargs={
            "observed_db": observed_db,
            "observed_sigma_db": observed_sigma_db,
        },
    )

    prior_template = prior_info.param_info.z
    full_template = full_info.param_info.z
    if prior_template.keys() != full_template.keys():
        raise RuntimeError("Prior and full models expose different latent sites.")
    for name in prior_template:
        if np.shape(prior_template[name]) != np.shape(full_template[name]):
            raise RuntimeError(f"Latent site {name!r} has inconsistent shapes.")

    flat_template, unravel_latent = ravel_pytree(prior_template)
    prior_potential_fn = prior_info.potential_fn
    full_potential_fn = full_info.potential_fn
    postprocess_fn = prior_info.postprocess_fn

    @jax.jit
    def log_prior_flat(z_flat):
        return -prior_potential_fn(unravel_latent(z_flat))

    @jax.jit
    def log_likelihood_flat(z_flat):
        z = unravel_latent(z_flat)
        return prior_potential_fn(z) - full_potential_fn(z)

    @jax.jit
    def log_posterior_flat(z_flat, beta=1.0):
        z = unravel_latent(z_flat)
        prior_potential = prior_potential_fn(z)
        full_potential = full_potential_fn(z)
        return -(1.0 - beta) * prior_potential - beta * full_potential

    @jax.jit
    def physical_sites_from_flat(z_flat):
        return postprocess_fn(unravel_latent(z_flat))

    return {
        "kind": kind,
        "n_dipoles": n_dipoles,
        "prior_model": prior_model,
        "latent_dim": int(flat_template.size),
        "log_prior_flat": log_prior_flat,
        "log_likelihood_flat": log_likelihood_flat,
        "log_posterior_flat": log_posterior_flat,
        "physical_sites_from_flat": physical_sites_from_flat,
        "log_prior_particles": jax.jit(jax.vmap(log_prior_flat, in_axes=0)),
        "log_likelihood_particles": jax.jit(
            jax.vmap(log_likelihood_flat, in_axes=0)
        ),
        "physical_sites": jax.jit(
            jax.vmap(physical_sites_from_flat, in_axes=0)
        ),
    }


# -----------------------------------------------------------------------------
# Weighting, tempering, and resampling
# -----------------------------------------------------------------------------


@jax.jit
def normalize_log_weights(log_weights):
    return log_weights - jax.scipy.special.logsumexp(log_weights)


@jax.jit
def effective_sample_size_from_log_weights(log_weights):
    log_weights = normalize_log_weights(log_weights)
    return jnp.exp(-jax.scipy.special.logsumexp(2.0 * log_weights))


@jax.jit
def incremental_log_weights(log_likelihood, beta_previous, beta_candidate):
    return (beta_candidate - beta_previous) * log_likelihood


@jax.jit
def updated_log_weights(log_weights, log_likelihood, beta_previous, beta_candidate):
    return log_weights + incremental_log_weights(
        log_likelihood, beta_previous, beta_candidate
    )


@jax.jit
def candidate_effective_sample_size(
    log_weights, log_likelihood, beta_previous, beta_candidate
):
    return effective_sample_size_from_log_weights(
        updated_log_weights(log_weights, log_likelihood, beta_previous, beta_candidate)
    )


def find_next_beta(
    log_weights,
    log_likelihood,
    beta_previous,
    target_ess,
    *,
    beta_tolerance: float = 1e-10,
    max_bisection_steps: int = 80,
):
    """Choose the largest next beta retaining the requested ESS."""

    beta_previous = float(beta_previous)
    target_ess = float(target_ess)
    n_particles = int(log_likelihood.shape[0])
    if not 0.0 <= beta_previous <= 1.0:
        raise ValueError("beta_previous must lie in [0, 1].")
    if not 1.0 <= target_ess <= n_particles:
        raise ValueError("target_ess must lie between 1 and n_particles.")

    if float(
        candidate_effective_sample_size(
            log_weights, log_likelihood, beta_previous, 1.0
        )
    ) >= target_ess:
        return 1.0

    lower, upper = beta_previous, 1.0
    for _ in range(max_bisection_steps):
        midpoint = 0.5 * (lower + upper)
        ess = float(
            candidate_effective_sample_size(
                log_weights, log_likelihood, beta_previous, midpoint
            )
        )
        if ess >= target_ess:
            lower = midpoint
        else:
            upper = midpoint
        if upper - lower <= beta_tolerance:
            break
    return lower


@jax.jit
def systematic_resample_indices(rng_key, normalized_log_weights):
    n_particles = normalized_log_weights.shape[0]
    weights = jnp.exp(normalized_log_weights)
    cumulative = jnp.cumsum(weights).at[-1].set(1.0)
    offset = jax.random.uniform(
        rng_key,
        shape=(),
        minval=0.0,
        maxval=1.0 / n_particles,
        dtype=weights.dtype,
    )
    positions = offset + jnp.arange(n_particles, dtype=weights.dtype) / n_particles
    return jnp.searchsorted(cumulative, positions, side="right").astype(jnp.int32)


@jax.jit
def systematic_resample_state(state: SMCState, rng_key) -> SMCState:
    ancestors = systematic_resample_indices(rng_key, state.log_weights)
    n_particles = state.particles.shape[0]
    return SMCState(
        particles=state.particles[ancestors],
        log_prior=state.log_prior[ancestors],
        log_likelihood=state.log_likelihood[ancestors],
        log_weights=jnp.full(
            n_particles, -jnp.log(n_particles), dtype=state.log_weights.dtype
        ),
        beta=state.beta,
        ancestor_indices=ancestors,
    )


# -----------------------------------------------------------------------------
# Parallel BlackJAX NUTS
# -----------------------------------------------------------------------------


def build_parallel_nuts_functions(model: dict, config: SMCConfig):
    """Build reusable ``pmap(vmap(...))`` NUTS initialization/step functions."""

    config.validate()
    devices = jax.devices("cpu")[: config.n_devices]
    if len(devices) < config.n_devices:
        raise ValueError(
            f"Requested {config.n_devices} CPU devices, but JAX exposes "
            f"only {len(jax.devices('cpu'))}."
        )

    log_posterior_flat = model["log_posterior_flat"]
    nuts_kernel = blackjax.nuts.build_kernel(
        divergence_threshold=config.nuts_divergence_threshold
    )

    def initialize_one(particle, beta):
        def logdensity_fn(z_flat):
            return log_posterior_flat(z_flat, beta=beta)

        return blackjax.nuts.init(particle, logdensity_fn)

    def step_one(rng_key, state, beta, step_size, inverse_mass_matrix):
        def logdensity_fn(z_flat):
            return log_posterior_flat(z_flat, beta=beta)

        return nuts_kernel(
            rng_key,
            state,
            logdensity_fn,
            step_size,
            inverse_mass_matrix,
            max_num_doublings=config.nuts_max_num_doublings,
        )

    initialize_parallel = jax.pmap(
        jax.vmap(initialize_one, in_axes=(0, None)),
        in_axes=(0, None),
        devices=devices,
    )
    step_parallel = jax.pmap(
        jax.vmap(step_one, in_axes=(0, 0, None, None, None)),
        in_axes=(0, 0, None, None, None),
        devices=devices,
    )
    return initialize_parallel, step_parallel


def _shard_particles(particles, n_devices):
    n_particles = particles.shape[0]
    if n_particles % n_devices:
        raise ValueError(
            f"{n_particles} particles cannot be divided evenly across {n_devices} devices."
        )
    per_device = n_particles // n_devices
    return particles.reshape((n_devices, per_device) + particles.shape[1:])


def _flatten_sharded_pytree(pytree):
    return jax.tree_util.tree_map(
        lambda a: a.reshape((a.shape[0] * a.shape[1],) + a.shape[2:]), pytree
    )


def mutate_particles_with_adaptive_nuts(
    particles,
    rng_key,
    beta,
    *,
    initial_step_size,
    inverse_mass_matrix,
    initialize_parallel,
    step_parallel,
    config: SMCConfig,
):
    """Apply one adaptively tuned population-wide parallel NUTS mutation."""

    n_particles = particles.shape[0]
    initial_states = initialize_parallel(
        _shard_particles(particles, config.n_devices),
        jnp.asarray(beta, dtype=config.dtype),
    )
    step_size = float(initial_step_size)

    for attempt in range(1, config.nuts_max_mutation_attempts + 1):
        rng_key, attempt_key = jax.random.split(rng_key)
        per_device = n_particles // config.n_devices
        keys = jax.random.split(attempt_key, n_particles)
        keys = keys.reshape(
            config.n_devices,
            per_device,
            *keys.shape[1:],
        )
        states_sharded, info_sharded = step_parallel(
            keys,
            initial_states,
            jnp.asarray(beta, dtype=config.dtype),
            jnp.asarray(step_size, dtype=config.dtype),
            jnp.asarray(inverse_mass_matrix, dtype=config.dtype),
        )
        states = _flatten_sharded_pytree(states_sharded)
        info = _flatten_sharded_pytree(info_sharded)

        acceptance = np.asarray(info.acceptance_rate)
        divergent = np.asarray(info.is_divergent)
        mean_acceptance = float(acceptance.mean())
        divergence_fraction = float(divergent.mean())
        all_finite = bool(np.asarray(jnp.all(jnp.isfinite(states.position))))

        print(
            f"    mutation attempt {attempt}: step={step_size:.3e}, "
            f"acceptance={mean_acceptance:.3f}, "
            f"divergences={divergent.sum()}/{n_particles}"
        )

        if (
            all_finite
            and divergence_fraction <= config.nuts_max_trial_divergence_fraction
            and mean_acceptance >= config.nuts_target_acceptance_low
        ):
            # Remove committed pmap sharding before carrying particles into the
            # next SMC stage.
            accepted_particles = jnp.asarray(
                jax.device_get(states.position), dtype=particles.dtype
            )
            accepted_info = info
            accepted_attempt = attempt
            break

        step_size *= config.nuts_step_size_decrease_factor
        if step_size < config.nuts_min_step_size:
            raise RuntimeError(
                "NUTS step size fell below the configured minimum "
                f"({config.nuts_min_step_size:.3e})."
            )
    else:
        raise RuntimeError(
            "No acceptable NUTS mutation was found after "
            f"{config.nuts_max_mutation_attempts} attempts."
        )

    if (
        mean_acceptance > config.nuts_target_acceptance_high
        and divergence_fraction == 0.0
    ):
        next_step_size = min(
            step_size * config.nuts_step_size_increase_factor,
            config.nuts_max_step_size,
        )
    else:
        next_step_size = step_size

    return (
        accepted_particles,
        accepted_info,
        step_size,
        next_step_size,
        accepted_attempt,
    )


# -----------------------------------------------------------------------------
# Annealed SMC
# -----------------------------------------------------------------------------


def initialize_smc_state(model: dict, config: SMCConfig, seed: int) -> SMCState:
    """Draw the initial equally weighted particle population from the prior."""

    keys = jax.random.split(jax.random.key(seed), config.n_particles)
    particles = []
    for key in keys:
        info = initialize_model(
            key, model["prior_model"], init_strategy=init_to_sample()
        )
        flat, _ = ravel_pytree(info.param_info.z)
        if flat.shape != (model["latent_dim"],):
            raise RuntimeError(
                f"Unexpected latent shape {flat.shape}; expected "
                f"({model['latent_dim']},)."
            )
        particles.append(flat)

    particles = jnp.stack(particles)
    log_prior = model["log_prior_particles"](particles)
    log_likelihood = model["log_likelihood_particles"](particles)
    log_likelihood.block_until_ready()

    if not bool(
        jnp.all(jnp.isfinite(particles))
        & jnp.all(jnp.isfinite(log_prior))
        & jnp.all(jnp.isfinite(log_likelihood))
    ):
        raise FloatingPointError("Initial SMC population contains non-finite values.")

    n = config.n_particles
    return SMCState(
        particles=particles,
        log_prior=log_prior,
        log_likelihood=log_likelihood,
        log_weights=jnp.full(n, -jnp.log(n), dtype=config.dtype),
        beta=jnp.asarray(0.0, dtype=config.dtype),
        ancestor_indices=jnp.arange(n, dtype=jnp.int32),
    )


def smc_transition(
    state: SMCState,
    rng_key,
    *,
    model: dict,
    config: SMCConfig,
    nuts_step_size,
    inverse_mass_matrix,
    initialize_parallel,
    step_parallel,
):
    """Advance the particle population through one adaptive tempering stage."""

    n_particles = state.particles.shape[0]
    target_ess = config.ess_target_fraction * n_particles
    beta_next = find_next_beta(
        state.log_weights,
        state.log_likelihood,
        state.beta,
        target_ess,
        beta_tolerance=config.beta_tolerance,
        max_bisection_steps=config.max_bisection_steps,
    )

    incremental = incremental_log_weights(
        state.log_likelihood, state.beta, beta_next
    )
    updated = updated_log_weights(
        state.log_weights, state.log_likelihood, state.beta, beta_next
    )
    normalized = normalize_log_weights(updated)
    ess = effective_sample_size_from_log_weights(normalized)
    log_z_increment = jax.scipy.special.logsumexp(
        state.log_weights + incremental
    )

    weighted_state = state._replace(
        log_weights=normalized,
        beta=jnp.asarray(beta_next, dtype=config.dtype),
    )
    resampling_key, mutation_key = jax.random.split(rng_key)
    resampled_state = systematic_resample_state(weighted_state, resampling_key)

    (
        mutated_particles,
        mutation_info,
        accepted_step_size,
        next_step_size,
        mutation_attempts,
    ) = mutate_particles_with_adaptive_nuts(
        resampled_state.particles,
        mutation_key,
        resampled_state.beta,
        initial_step_size=nuts_step_size,
        inverse_mass_matrix=inverse_mass_matrix,
        initialize_parallel=initialize_parallel,
        step_parallel=step_parallel,
        config=config,
    )

    log_prior = model["log_prior_particles"](mutated_particles)
    log_likelihood = model["log_likelihood_particles"](mutated_particles)

    new_state = SMCState(
        particles=mutated_particles,
        log_prior=log_prior,
        log_likelihood=log_likelihood,
        log_weights=resampled_state.log_weights,
        beta=resampled_state.beta,
        ancestor_indices=resampled_state.ancestor_indices,
    )

    unique_ancestors = jnp.unique(
        resampled_state.ancestor_indices,
        size=n_particles,
        fill_value=-1,
    )
    displacement = jnp.linalg.norm(
        mutated_particles - resampled_state.particles, axis=1
    )
    diagnostics = SMCStageDiagnostics(
        beta_previous=state.beta,
        beta_next=resampled_state.beta,
        delta_beta=resampled_state.beta - state.beta,
        ess_before_resampling=ess,
        log_evidence_increment=log_z_increment,
        n_unique_ancestors=jnp.sum(unique_ancestors >= 0),
        nuts_step_size=jnp.asarray(accepted_step_size, dtype=config.dtype),
        mutation_attempts=jnp.asarray(mutation_attempts, dtype=jnp.int32),
        mean_acceptance_rate=jnp.mean(mutation_info.acceptance_rate),
        minimum_acceptance_rate=jnp.min(mutation_info.acceptance_rate),
        n_divergent=jnp.sum(mutation_info.is_divergent),
        median_integration_steps=jnp.median(mutation_info.num_integration_steps),
        maximum_integration_steps=jnp.max(mutation_info.num_integration_steps),
        median_displacement=jnp.median(displacement),
    )
    return new_state, diagnostics, next_step_size


def run_annealed_smc(
    initial_state: SMCState,
    rng_key,
    *,
    model: dict,
    config: SMCConfig,
    nuts_step_size=None,
    inverse_mass_matrix=None,
    initialize_parallel,
    step_parallel,
    print_progress: bool = True,
):
    """Temper an SMC population from beta=0 to beta=1."""

    state = initial_state
    current_step_size = float(
        config.nuts_initial_step_size if nuts_step_size is None else nuts_step_size
    )
    if inverse_mass_matrix is None:
        inverse_mass_matrix = jnp.ones(model["latent_dim"], dtype=config.dtype)

    diagnostics = []
    cumulative_log_evidence = [0.0]
    n_particles = state.particles.shape[0]

    if print_progress:
        print("Annealed SMC run")
        print("----------------")

    for stage_index in range(config.max_stages):
        if float(state.beta) >= 1.0 - 1e-7:
            break

        rng_key, stage_key = jax.random.split(rng_key)
        new_state, stage_diag, next_step_size = smc_transition(
            state,
            stage_key,
            model=model,
            config=config,
            nuts_step_size=current_step_size,
            inverse_mass_matrix=inverse_mass_matrix,
            initialize_parallel=initialize_parallel,
            step_parallel=step_parallel,
        )
        new_state.log_likelihood.block_until_ready()

        beta_previous = float(stage_diag.beta_previous)
        beta_next = float(stage_diag.beta_next)
        delta_beta = float(stage_diag.delta_beta)
        if not bool(
            jnp.all(jnp.isfinite(new_state.particles))
            & jnp.all(jnp.isfinite(new_state.log_prior))
            & jnp.all(jnp.isfinite(new_state.log_likelihood))
            & jnp.all(jnp.isfinite(new_state.log_weights))
        ):
            raise FloatingPointError(
                f"Stage {stage_index + 1} produced non-finite values."
            )
        if beta_next <= beta_previous:
            raise RuntimeError(
                f"Stage {stage_index + 1} did not advance beta: "
                f"{beta_previous:.12e} -> {beta_next:.12e}."
            )
        if delta_beta < config.min_beta_increment and beta_next < 1.0 - 1e-7:
            raise RuntimeError(
                f"Annealing stalled at stage {stage_index + 1}: "
                f"delta beta={delta_beta:.3e}."
            )

        diagnostics.append(stage_diag)
        cumulative_log_evidence.append(
            cumulative_log_evidence[-1] + float(stage_diag.log_evidence_increment)
        )
        state = new_state
        current_step_size = float(next_step_size)

        if print_progress and (
            stage_index % config.print_every == 0 or beta_next >= 1.0 - 1e-7
        ):
            print(
                f"Stage {stage_index + 1:03d}: beta={beta_next:.6e}, "
                f"delta={delta_beta:.2e}, "
                f"ESS={float(stage_diag.ess_before_resampling):5.1f}, "
                f"ancestors={int(stage_diag.n_unique_ancestors):4d}/{n_particles}, "
                f"step={float(stage_diag.nuts_step_size):.2e}, "
                f"attempts={int(stage_diag.mutation_attempts)}, "
                f"acc={float(stage_diag.mean_acceptance_rate):.3f}, "
                f"div={int(stage_diag.n_divergent):2d}, "
                f"steps={float(stage_diag.median_integration_steps):4.1f}, "
                f"logZ={cumulative_log_evidence[-1]:.3f}"
            )
    else:
        raise RuntimeError(f"Reached max_stages={config.max_stages} before beta=1.")

    return (
        state,
        diagnostics,
        np.asarray(cumulative_log_evidence),
        current_step_size,
    )


def run_final_mutations(
    state: SMCState,
    rng_key,
    *,
    model: dict,
    config: SMCConfig,
    nuts_step_size,
    inverse_mass_matrix,
    initialize_parallel,
    step_parallel,
    n_mutation_steps=None,
    print_progress=True,
):
    """Apply pure beta=1 NUTS rejuvenation after annealing."""

    if not np.isclose(float(state.beta), 1.0, rtol=0.0, atol=1e-6):
        raise ValueError("Final posterior mutations require beta=1.")
    if n_mutation_steps is None:
        n_mutation_steps = config.n_final_mutation_steps

    current_state = state
    current_step_size = float(nuts_step_size)
    diagnostics = []

    for mutation_index in range(1, n_mutation_steps + 1):
        rng_key, mutation_key = jax.random.split(rng_key)
        (
            particles,
            info,
            accepted_step_size,
            next_step_size,
            attempts,
        ) = mutate_particles_with_adaptive_nuts(
            current_state.particles,
            mutation_key,
            current_state.beta,
            initial_step_size=current_step_size,
            inverse_mass_matrix=inverse_mass_matrix,
            initialize_parallel=initialize_parallel,
            step_parallel=step_parallel,
            config=config,
        )
        log_prior = model["log_prior_particles"](particles)
        log_likelihood = model["log_likelihood_particles"](particles)
        log_likelihood.block_until_ready()
        displacement = jnp.linalg.norm(particles - current_state.particles, axis=1)

        diag = {
            "mutation_index": mutation_index,
            "nuts_step_size": float(accepted_step_size),
            "mutation_attempts": int(attempts),
            "mean_acceptance_rate": float(jnp.mean(info.acceptance_rate)),
            "minimum_acceptance_rate": float(jnp.min(info.acceptance_rate)),
            "n_divergent": int(jnp.sum(info.is_divergent)),
            "median_integration_steps": float(jnp.median(info.num_integration_steps)),
            "maximum_integration_steps": int(jnp.max(info.num_integration_steps)),
            "median_displacement": float(jnp.median(displacement)),
        }
        diagnostics.append(diag)
        current_state = current_state._replace(
            particles=particles,
            log_prior=log_prior,
            log_likelihood=log_likelihood,
        )
        current_step_size = float(next_step_size)

        if print_progress:
            print(
                f"Mutation {mutation_index:02d}: step={diag['nuts_step_size']:.3e}, "
                f"attempts={diag['mutation_attempts']}, "
                f"acc={diag['mean_acceptance_rate']:.3f}, "
                f"div={diag['n_divergent']}, "
                f"steps={diag['median_integration_steps']:.1f}"
            )

    return current_state, diagnostics, current_step_size


def run_smc(
    model: dict,
    config: SMCConfig,
    *,
    prior_seed: int,
    smc_seed: int,
    final_mutation_seed: int,
    inverse_mass_matrix=None,
    print_progress: bool = True,
) -> dict:
    """Run prior initialization, annealing, and final posterior rejuvenation."""

    config.validate()
    if inverse_mass_matrix is None:
        inverse_mass_matrix = jnp.ones(model["latent_dim"], dtype=config.dtype)

    initialize_parallel, step_parallel = build_parallel_nuts_functions(model, config)
    initial_state = initialize_smc_state(model, config, prior_seed)
    annealed_state, stage_diagnostics, log_evidence, final_step_size = run_annealed_smc(
        initial_state,
        jax.random.key(smc_seed),
        model=model,
        config=config,
        inverse_mass_matrix=inverse_mass_matrix,
        initialize_parallel=initialize_parallel,
        step_parallel=step_parallel,
        print_progress=print_progress,
    )
    posterior_state, final_diagnostics, posterior_step_size = run_final_mutations(
        annealed_state,
        jax.random.key(final_mutation_seed),
        model=model,
        config=config,
        nuts_step_size=final_step_size,
        inverse_mass_matrix=inverse_mass_matrix,
        initialize_parallel=initialize_parallel,
        step_parallel=step_parallel,
        print_progress=print_progress,
    )

    return {
        "initial_state": initial_state,
        "annealed_state": annealed_state,
        "posterior_state": posterior_state,
        "stage_diagnostics": stage_diagnostics,
        "final_mutation_diagnostics": final_diagnostics,
        "cumulative_log_evidence": log_evidence,
        "log_evidence": float(log_evidence[-1]),
        "final_nuts_step_size": final_step_size,
        "posterior_nuts_step_size": posterior_step_size,
    }


def physical_sites(model: dict, particles):
    """Transform flat posterior particles back to physical NumPyro sites."""

    return model["physical_sites"](particles)


# -----------------------------------------------------------------------------
# Pickle-friendly helpers
# -----------------------------------------------------------------------------


def state_to_numpy(state: SMCState) -> dict:
    return {
        "particles": np.asarray(state.particles),
        "log_prior": np.asarray(state.log_prior),
        "log_likelihood": np.asarray(state.log_likelihood),
        "log_weights": np.asarray(state.log_weights),
        "beta": float(state.beta),
        "ancestor_indices": np.asarray(state.ancestor_indices),
    }


def diagnostics_to_numpy(diagnostics) -> list[dict]:
    output = []
    for item in diagnostics:
        if isinstance(item, dict):
            output.append(item)
        else:
            output.append(
                {
                    name: np.asarray(value).item()
                    if np.asarray(value).ndim == 0
                    else np.asarray(value)
                    for name, value in item._asdict().items()
                }
            )
    return output


def physical_sites_to_numpy(sites) -> dict[str, np.ndarray]:
    return {name: np.asarray(value) for name, value in sites.items()}


__all__ = [
    "SMCConfig",
    "SMCState",
    "SMCStageDiagnostics",
    "build_amplitude_models",
    "build_complex_models",
    "build_model_interface",
    "build_parallel_nuts_functions",
    "candidate_effective_sample_size",
    "diagnostics_to_numpy",
    "effective_sample_size_from_log_weights",
    "find_next_beta",
    "incremental_log_weights",
    "initialize_smc_state",
    "mutate_particles_with_adaptive_nuts",
    "normalize_log_weights",
    "physical_sites",
    "physical_sites_to_numpy",
    "run_annealed_smc",
    "run_final_mutations",
    "run_smc",
    "smc_transition",
    "state_to_numpy",
    "systematic_resample_indices",
    "systematic_resample_state",
    "updated_log_weights",
]
