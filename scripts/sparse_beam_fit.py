from pathlib import Path
from embers_passes import PassFile, SphericalSpline, eval_spline_batch, \
    make_knots_from_grid, load_config, run_diagnostics, scatter_coeffs, \
    make_disk_mask, make_flat_index, make_constraint_info, inject_pivot, \
    chunk_passes, healpix_to_pyuvdata, bin_chunk, prep_mwa_beam, get_pol_ind,\
    get_res
    

import numpy as np
import healpy as hp

import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm


from pyuvdata import UVBeam

from functools import partial
import pickle

import arviz as az


def prep_spline_params(beam, ortho_knots=False, spice_variant=0):
    if ortho_knots:
        if spice_variant:
            if spice_variant > 3: # Increasing density only at horizon
                interval_1 = np.arange(-1, -0.6, 0.05 / (2**(spice_variant - 3)))
                interval_2 = np.arange(-0.6, -0.4, 0.025)
                interval_3 = np.arange(-0.4, 0., 0.1)
            elif spice_variant == 3: # Twice as dense except in main lobe
                interval_1 = np.arange(-1., -0.6, 0.025)
                interval_2 = np.arange(-0.6, -0.4, 0.0125)
                interval_3 = np.arange(-0.4, 0., 0.1)
            elif spice_variant == 2: # Twice as dense except at horizon
                interval_1 = np.arange(-1., -0.6, 0.05)
                interval_2 = np.arange(-0.6, -0.4, 0.0125)
                interval_3 = np.arange(-0.4, 0., 0.05)
            else: # Twice as dense everywhere
                interval_1 = np.arange(-1., -0.6, 0.025)
                interval_2 = np.arange(-0.6, -0.4, 0.0125)
                interval_3 = np.arange(-0.4, 0., 0.05)
        else:
            interval_1 = np.arange(-1., -0.6, 0.05)
            interval_2 = np.arange(-0.6, -0.4, 0.025)
            interval_3 = np.arange(-0.4, 0., 0.1)

        base_knots = np.concatenate([
            interval_1, 
            interval_2, 
            interval_3,
            np.atleast_1d([0.]),
            -np.flip(interval_3),
            -np.flip(interval_2),
            -np.flip(interval_1)
        ])

        t_x, t_y, p, q, n_x, n_y = make_knots_from_grid(base_knots, base_knots)
        Nparam = n_x * n_y

        return t_x, t_y, p, q, n_x, n_y, Nparam
    else:
        theta_param = beam.axis2_array[1::10] # avoid the origin by a quarter degree
        phi_param = beam.axis1_array[::10]

        # ONCE, outside the model:
        t_theta, t_phi, p, q, n_th, n_ph = make_knots_from_grid(theta_param, phi_param)
        Nparam = n_th * n_ph
        return t_theta, t_phi, p, q, n_th, n_ph, Nparam

def contributing_pass_count_map(passes, chunk, nside):
    """
    Number of independent passes contributing to each HEALPix pixel.

    A pass contributes at most once to a given pixel, regardless of how many
    samples from that pass land in the pixel.
    """
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

        # Each independent pass counts only once per pixel.
        n_pass_map[np.unique(pixels)] += 1

    return n_pass_map

def fit_single_pass_uncertainty_cut(
    count_maps,
    mad_maps,
    median_sem_maps,
    n_pass_maps,
    min_count=3,
    robust_z_cut=-3.0,
):
    """
    Estimate a robust lower-tail uncertainty cut for single-pass pixels.

    The distribution is defined in log10(sigma), pooling valid single-pass
    pixels across all chunks.

    Returns
    -------
    metadata : dict
        Median/scatter of log10(sigma), robust-z threshold, and equivalent
        sigma threshold in dB.
    """
    pooled_log_sigma = []

    for count_map, mad_map, sem_map, n_pass_map in zip(
        count_maps,
        mad_maps,
        median_sem_maps,
        n_pass_maps,
    ):
        count_map = np.asarray(count_map)
        n_pass_map = np.asarray(n_pass_map)

        sigma = np.sqrt(
            np.asarray(mad_map, dtype=np.float64) ** 2
            + np.asarray(sem_map, dtype=np.float64) ** 2
        )

        valid_single = (
            (count_map >= min_count)
            & (n_pass_map == 1)
            & np.isfinite(sigma)
            & (sigma > 0)
        )

        pooled_log_sigma.append(
            np.log10(sigma[valid_single])
        )

    pooled_log_sigma = np.concatenate(
        pooled_log_sigma
    )

    if pooled_log_sigma.size < 10:
        raise ValueError(
            "Too few valid single-pass pixels to define a robust "
            f"uncertainty cut ({pooled_log_sigma.size} found)."
        )

    median_log_sigma = np.median(
        pooled_log_sigma
    )

    # Gaussian-equivalent scatter estimated from MAD.
    robust_sigma_log = (
        1.482602218505602
        * np.median(
            np.abs(
                pooled_log_sigma
                - median_log_sigma
            )
        )
    )

    if (
        not np.isfinite(robust_sigma_log)
        or robust_sigma_log <= 0
    ):
        raise ValueError(
            "Could not estimate a finite positive robust scatter."
        )

    log_sigma_cut = (
        median_log_sigma
        + robust_z_cut * robust_sigma_log
    )

    sigma_cut_db = 10.0 ** log_sigma_cut

    return {
        "robust_z_cut": float(robust_z_cut),
        "pooled_n_single_pass_pixels": int(
            pooled_log_sigma.size
        ),
        "median_log10_sigma": float(
            median_log_sigma
        ),
        "robust_sigma_log10_sigma": float(
            robust_sigma_log
        ),
        "sigma_cut_db": float(
            sigma_cut_db
        ),
    }

def single_pass_low_uncertainty_mask(
    count_map,
    mad_map,
    median_sem_map,
    n_pass_map,
    cut_metadata,
    min_count=3,
):
    """
    Flag anomalously low-uncertainty single-pass pixels.

    Returns
    -------
    mask : bool array
        True for pixels to exclude.
    robust_z : float array
        Robust z-score in log10(sigma); NaN elsewhere.
    sigma : float array
        Combined per-pixel uncertainty.
    """
    count_map = np.asarray(count_map)
    n_pass_map = np.asarray(n_pass_map)

    sigma = np.sqrt(
        np.asarray(mad_map, dtype=np.float64) ** 2
        + np.asarray(median_sem_map, dtype=np.float64) ** 2
    )

    one_pass_valid = (
        (count_map >= min_count)
        & (n_pass_map == 1)
        & np.isfinite(sigma)
        & (sigma > 0)
    )

    robust_z = np.full(
        sigma.shape,
        np.nan,
        dtype=np.float64,
    )

    robust_z[one_pass_valid] = (
        (
            np.log10(
                sigma[one_pass_valid]
            )
            - cut_metadata[
                "median_log10_sigma"
            ]
        )
        / cut_metadata[
            "robust_sigma_log10_sigma"
        ]
    )

    mask = (
        one_pass_valid
        & (
            robust_z
            < cut_metadata[
                "robust_z_cut"
            ]
        )
    )

    return mask, robust_z, sigma


if __name__ == "__main__":

    import sys
    args = load_config(sys.argv[1])
    from shutil import copy
    import random as pyrandom
    import os

    from scipy.stats import binned_statistic_2d
    
    import numpyro
    numpyro.set_host_device_count(args.ndevice)
    from numpyro import distributions as dist
    from numpyro.infer import MCMC, NUTS, SVI, Trace_ELBO, init_to_value
    from numpyro.infer.autoguide import AutoDelta, AutoLaplaceApproximation
    
    from jax import random
    import jax.numpy as jnp
    import jax
    jax.config.update("jax_debug_nans", True)

    import arviz as az

    # FIXME: Hardcode -- only zenith
    beam = prep_mwa_beam(args.mwa_beamfile)
    pol_ind = get_pol_ind(beam, args.pol)

    subdir = f"rf{args.rf_num}/S{args.tile}/{args.pol}"
    if args.jackknife:
        subdir += f"/{args.jackknife_mode}_jackknife/chunk{args.jackknife_chunk}"
    outdir = f"{args.outdir}/{subdir}"
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    # Bring config over so it's easy to find later
    copy(sys.argv[1], f"{outdir}/")

    path = Path(f"{args.indir}/rf{args.rf_num}/S{args.tile}{args.pol}_rf{args.rf_num}{args.pol}_passes.h5")
    pf = PassFile(path)
    # FIXME: Only 0 pointing
    passes = pf.read_passes(pointing=0)
    Npass_total = len(passes)

    if args.jackknife:
        if args.jackknife_mode == "random": # shuffle the passes first
            pyrandom.seed(args.jackknife_key)
            pyrandom.shuffle(passes)
            if args.jackknife_chunk: # implicitly second half
                chunk = list(range(Npass_total // 2, Npass_total))
            else: # otherwise first half
                chunk = list(range(Npass_total // 2))
        elif args.jackknife_mode == "time":
            chunks, bins_edges, bin_centers = chunk_passes(passes, args.chunk_sec)
            chunk = chunks[args.jackknife_chunk]
        else:
            raise ValueError("Invalid jackknife mode. Should be 'time' or 'random'. Check config file." )
    elif args.num_pass is None:
        chunk = list(range(Npass_total))
    else:
        chunk = list(range(args.num_pass))

    median_map, count_map, mad_map, median_sem_map = bin_chunk(
        passes, chunk
    )

    n_pass_map = contributing_pass_count_map(
        passes, chunk, nside=32
    )

    cut_metadata = fit_single_pass_uncertainty_cut(
        [count_map], [mad_map], [median_sem_map], [n_pass_map]
    )

    mask, robust_z, sigma = single_pass_low_uncertainty_mask(
        count_map, mad_map, median_sem_map, n_pass_map, cut_metadata
    )

    # No great way to pass this mask through healpix_to_pyuvdata -- set counts to 0
    count_map[mask] = 0

    median_values, za_rad, az_rad = healpix_to_pyuvdata(median_map)
    count_values, _, _ = healpix_to_pyuvdata(count_map)
    mad_values, _, _ = healpix_to_pyuvdata(mad_map)
    median_sem_values, _, _ = healpix_to_pyuvdata(median_sem_map)

    count_gt_2 = count_values > 2 # MAD estimates completely untrustworthy for counts less than 3
    median_values = median_values[count_gt_2]
    za_rad = za_rad[count_gt_2]
    az_rad = az_rad[count_gt_2]
    mad_values = mad_values[count_gt_2]
    count_values = count_values[count_gt_2]
    median_sem_values = median_sem_values[count_gt_2]

    # Do some Bayesian inference with MCMC (NUTS) on some short chains
    if args.ortho_knots:
        # Basic knot stuff
        t_x, t_y, p, q, n_x, n_y, _= prep_spline_params(beam, ortho_knots=True)
        knot_mask = make_disk_mask(t_x, t_y, p, q)           # (n_x, n_y) bool, numpy
        knot_ij, Nparam = make_flat_index(knot_mask)            # static arrays
        if args.enforce_boresight: # Sample one less parameter to satisfy constraint
            Nparam -= 1


        (
            pivot_flat_idx, pivot_weight, contrib_flat_idxs, contrib_weights
        ) = make_constraint_info(
            t_x, t_y, p, q, u0=0, v0=0, ij=knot_ij # FIXME: Hardcode zenith-pointed
        )
        dat_coord_1 = jnp.sin(za_rad) * jnp.cos(az_rad)
        dat_coord_2 = jnp.sin(za_rad) * jnp.sin(az_rad)
    else:
        t_theta, t_phi, p, q, n_th, n_ph, Nparam = prep_spline_params(beam)
        dat_coord_1 = za_rad
        dat_coord_2 = az_rad

    res, model_beam = get_res(az_rad, za_rad, median_values, pol_ind, beam)

    class SkewLogistic(dist.Distribution, metaclass=dist.distribution.DistributionMeta):
        """
        Type I Generlized Logistic Distribution.
        """
        arg_constraints = {
            "loc": dist.constraints.real,
            "scale": dist.constraints.positive,
            "shape": dist.constraints.positive
        }
        support = dist.constraints.real
        has_rsample = False

        def __init__(self, loc=0., scale=1., shape=1., *, validate_args=None):
            self.loc = loc
            self.scale = scale
            self.shape = shape
            

            batch_shape = jax.lax.broadcast_shapes(
                jnp.shape(loc),
                jnp.shape(scale),
                jnp.shape(shape),
            )
            super(SkewLogistic, self).__init__(batch_shape=batch_shape, validate_args=validate_args)

        def __new__(cls, *args, **kwargs):
            return object.__new__(cls)

        def log_prob(self, x):
            centered = (x - self.loc) / self.scale
            return -jnp.log(self.scale) + jnp.log(self.shape) - centered - (self.shape + 1) * jnp.log1p(jnp.exp(-centered))

        def sample(self, key, sample_shape=()):
            final_sample_shape = self.batch_shape + sample_shape
            beta_samp = random.beta(key, self.shape, jnp.ones_like(self.shape), shape=final_sample_shape)
            centered_samp = jnp.log(beta_samp / (1 - beta_samp))
            return self.scale * centered_samp + self.loc

    def vanilla_model(
            dat_coord_1, 
            dat_coord_2, 
            Ndat,
            noise_scale,
            data=None, 
            ortho_knots=False,
            enforce_boresight=False,
    ):
        scale_vals = numpyro.sample("scale_vals", dist.Uniform(low=1, high=1e2))
        shape_vals = numpyro.sample("shape_vals", dist.Uniform(low=1e-2, high=1e2))
        scale_grad = numpyro.sample("scale_grad", dist.Uniform(low=1, high=1e2))
        with numpyro.plate("Nparam", Nparam): # Need to enforce peak-normalization

            surf_vals = numpyro.sample(
                "surf_vals",
                SkewLogistic(
                    loc=0,
                    scale=scale_vals,
                    shape=shape_vals
                )
            )
        spl, grad = get_spl_grad(ortho_knots, enforce_boresight, surf_vals)
        model_vals = eval_spline_batch(spl, dat_coord_1, dat_coord_2) # may be x_data, y_data
        
        lp1 = dist.SoftLaplace(loc=0, scale=scale_grad).log_prob(grad[0]).sum()
        lp2 = dist.SoftLaplace(loc=0, scale=scale_grad).log_prob(grad[1]).sum()
        numpyro.factor("grad_pot", lp1 + lp2)
        
        with numpyro.plate("Ndat", Ndat):
            obs = numpyro.sample(
                "obs",
                dist.SoftLaplace(
                    loc=model_vals, 
                    scale=noise_scale * 2 / jnp.pi
                ),
                obs=data
            )

    def get_spl_grad(ortho_knots, enforce_boresight, surf_vals):
        if enforce_boresight:
            # Hardcoded 0 since log-beam 0 at peak
            surf_vals_model = inject_pivot(
                surf_vals, 
                pivot_flat_idx, 
                pivot_weight,
                contrib_flat_idxs, 
                contrib_weights, 
                0
            ) 
        else:
            surf_vals_model = surf_vals
            
        if ortho_knots:
            surf_vals_model = scatter_coeffs(surf_vals_model, knot_mask, knot_ij)
            spl = SphericalSpline(t_x, t_y, surf_vals_model, p, q)

            grad = jnp.gradient(surf_vals_model)
            grad = (grad[0][knot_mask], grad[1][knot_mask])
        else:
            spl = SphericalSpline(t_theta, t_phi, surf_vals_model.reshape(n_th, n_ph), p, q)
            grad = jnp.gradient(surf_vals_model.reshape(n_th, n_ph))
        return spl,grad
            
    def mixture_model(za_data, az_data, Ndat, data=None):
        scale_vals = numpyro.sample("scale_vals", dist.Uniform(low=1e-2, high=1e2))
        skewlog_scale_vals = numpyro.sample("skewlog_scale_vals", dist.Uniform(low=1e-2, high=1e2))
        skewlog_shape_vals = numpyro.sample("skewlog_shape_vals", dist.Uniform(low=1e-2, high=1e2))
        scale_grad = numpyro.sample("scale_grad", dist.Uniform(low=1e-2, high=1e2))
        scale_noise = numpyro.sample("scale_noise", dist.Uniform(low=1e-2, high=1e2))
        mix_weight = numpyro.sample("mix_weight", dist.Beta(11, 1))
        with numpyro.plate("Nparam", Nparam):
            surf_vals = numpyro.sample("surf_vals",
                                    dist.MixtureGeneral(
                                        dist.Categorical(
                                            jnp.array([
                                                mix_weight,
                                                1 - mix_weight
                                            ])
                                        ), 
                                        [
                                            dist.SoftLaplace(
                                            loc=0, scale=scale_vals
                                        ), 
                                            SkewLogistic(
                                                loc=0,
                                                scale=skewlog_scale_vals,
                                                shape=skewlog_shape_vals
                                        )
                                        ]
                                        )
                                    )
        spl = SphericalSpline(t_theta, t_phi, surf_vals.reshape(n_th, n_ph), p, q)
        model_vals = eval_spline_batch(spl, za_data, az_data)
        grad = jnp.gradient(surf_vals.reshape(n_th, n_ph))
        pot1 = dist.SoftLaplace(loc=0, scale=scale_grad).log_prob(grad[0]).sum()
        pot2 = dist.SoftLaplace(loc=0, scale=scale_grad).log_prob(grad[1]).sum()

        numpyro.factor("grad_pot", pot1 + pot2)
        
        
        with numpyro.plate("Ndat", Ndat):
            obs = numpyro.sample("obs",
                                dist.SoftLaplace(loc=model_vals, scale=scale_noise),
                                obs=data)
            
    model = mixture_model if args.mix_model else vanilla_model

    model_args = (
        dat_coord_1, 
        dat_coord_2, 
        len(median_values),
        np.sqrt(mad_values**2 + median_sem_values**2),
    )
    model_kwargs = {
        "data": res.real, 
        "ortho_knots": args.ortho_knots, 
        "enforce_boresight": args.enforce_boresight
    }

    # dat_for_inference = mean_res.real
    # count_gt_0 = count_res > 0
    # dat_for_inference = mean_res.real[count_gt_0]
    # model_args = (
    #     X[count_gt_0],
    #     Y[count_gt_0],
    #     len(dat_for_inference),
    #     1/np.sqrt(count_res)[count_gt_0]
    # )
    # model_kwargs = {
    #     "data": dat_for_inference,
    #     "ortho_knots": args.ortho_knots,
    #     "enforce_boresight": args.enforce_boresight
    # }
    if args.inference: # Need to do inference
        key = random.key(args.key)
        if args.svi:
            #guide = AutoDelta(model)
            """
            Train an AutoLaplace guide to extract the correlation structure
            of the spline coefficients. Use this in a multivariate normal guide
            to debias the posterior of the scale parameters. 
            """
            guide = AutoLaplaceApproximation(model)
            svi = SVI(model, guide, numpyro.optim.Adam(args.learning_rate), loss=Trace_ELBO())
            svi_result = svi.run(
                random.key(args.svi_key), 
                args.num_warmup, 
                *model_args, 
                **model_kwargs
            )
            params = svi_result.params
            plt.plot(svi_result.losses)
            plt.ylabel("Loss")
            plt.xlabel("Iteration")
            plt.savefig(f"{outdir}/losses.pdf")

            svi_samps = guide.sample_posterior(
                key, 
                params,
                *model_args,
                sample_shape=(args.num_sample,),
                **model_kwargs,
            )

            with open(f"{outdir}/svi_laplace_params.pkl", "wb") as svi_params_file:
                pickle.dump(params, svi_params_file)
            with open(f"{outdir}/svi_laplace_samps.pkl", "wb") as svi_samps_file:
                pickle.dump(svi_samps, svi_samps_file)

            if args.double_svi:

                def get_index_map(guide):
                    """Map each latent site name -> array of positions in the flattened
                    unconstrained vector, using the guide's own unpack logic (so it's
                    correct regardless of internal ordering)."""
                    latent_dim = guide.latent_dim
                    idx = jnp.arange(latent_dim, dtype=jnp.result_type(float))
                    unpacked = guide._unpack_latent(idx)
                    return {name: np.asarray(v).reshape(-1).astype(int) for name, v in unpacked.items()}


                def subset_scale_tril(guide, params, names):
                    """Lower-Cholesky factor of the approximate posterior restricted to `names`."""
                    posterior = guide.get_posterior(params)
                    scale_tril_full = np.asarray(posterior.scale_tril)
                    cov_full = scale_tril_full @ scale_tril_full.T

                    idx_map = get_index_map(guide)
                    idx = np.concatenate([idx_map[n] for n in names])

                    cov_sub = cov_full[np.ix_(idx, idx)]
                    scale_tril_sub = np.linalg.cholesky(cov_sub)
                    return scale_tril_sub, idx_map


                def extract_scale_tril(auto_laplace_guide, params, site_name):
                    """Full (not just correlation) Cholesky factor of the Laplace covariance
                    for one site, to be frozen and reused as-is in the custom guide."""
                    L_sub, idx_map = subset_scale_tril(auto_laplace_guide, params, [site_name])
                    cov = L_sub @ L_sub.T
                    cov = (cov + cov.T) / 2  # symmetrize away fp noise

                    # nudge onto the PD cone if numerically borderline, same as before
                    eigval_floor = 1e-6
                    w, v = np.linalg.eigh(cov)
                    w = np.clip(w, eigval_floor, None)
                    cov = (v * w) @ v.T

                    L_cov = np.linalg.cholesky(cov)
                    return jnp.asarray(L_cov)


                def make_custom_guide(Nparam, L_cov_fixed):
                    def guide(dat_coord_1, dat_coord_2, Ndat, noise_scale, data=None,
                            ortho_knots=False, enforce_boresight=False):

                        # --- surf_vals: fixed full covariance from Laplace, only loc is fit ---
                        loc_surf = numpyro.param("surf_vals_loc", jnp.zeros(Nparam))

                        surf_vals = numpyro.sample(
                            "surf_vals_aux",
                            dist.MultivariateNormal(loc=loc_surf, scale_tril=L_cov_fixed),
                            infer={"is_auxiliary": True},
                        )

                        with numpyro.plate("Nparam", Nparam):
                            numpyro.sample(
                                "surf_vals",
                                dist.Delta(surf_vals, event_dim=0).mask(False),
                            )

                        # --- bounded scalar params: diagonal normal pushed through biject_to ---
                        for name, low, high in [
                            ("scale_vals", 1.0, 1e2),
                            ("shape_vals", 1e-2, 1e2),
                            ("scale_grad", 1.0, 1e2),
                        ]:
                            loc = numpyro.param(f"{name}_loc", 0.0)
                            log_scale = numpyro.param(f"{name}_log_scale", 0.0)
                            base = dist.Normal(loc, jnp.exp(log_scale))
                            numpyro.sample(
                                name,
                                dist.TransformedDistribution(
                                    base, dist.transforms.biject_to(dist.constraints.interval(low, high))
                                ),
                            )

                    return guide


                L_cov_fixed = extract_scale_tril(guide, params, "surf_vals")
                custom_guide = make_custom_guide(Nparam, L_cov_fixed)
                svi = SVI(vanilla_model, custom_guide, numpyro.optim.Adam(1e-1), loss=Trace_ELBO(num_particles=8))
                svi_result = svi.run(random.key(2356), 2 * args.num_sample, *model_args, **model_kwargs)

                params = svi_result.params
                plt.plot(svi_result.losses)
                plt.ylabel("Loss")
                plt.xlabel("Iteration")
                plt.savefig(f"{outdir}/custom_guide_losses.pdf")

                from numpyro.infer import Predictive

                predictive = Predictive(custom_guide, params=svi_result.params, num_samples=1000)
                svi_samps = predictive(
                    random.key(562), *model_args,
                    ortho_knots=model_kwargs["ortho_knots"], 
                    enforce_boresight=model_kwargs["enforce_boresight"],
                    # data=None  -- guide doesn't use data, but pass same non-data args as training
                )

                with open(f"{outdir}/svi_params.pkl", "wb") as svi_params_file:
                    pickle.dump(params, svi_params_file)
                with open(f"{outdir}/svi_samps.pkl", "wb") as svi_samps_file:
                    pickle.dump(svi_samps, svi_samps_file)

            
        else:
            kernel = NUTS(model, dense_mass=args.dense_mass)
            mcmc = MCMC(
                kernel, 
                num_warmup=args.num_warmup, 
                num_samples=args.num_sample, 
                num_chains=args.ndevice
            )
            mcmc.run(
                key, 
                *model_args, 
                **model_kwargs
            )
            idata = az.from_numpyro(mcmc)
            idata.to_netcdf(f"{outdir}/mcmc_out.nc")

            run_diagnostics(idata)
    