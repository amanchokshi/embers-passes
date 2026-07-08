from pathlib import Path
from embers_passes import PassFile, SphericalSpline, eval_spline_batch, \
    make_knots_from_grid, load_config, run_diagnostics, scatter_coeffs, \
    make_disk_mask, make_flat_index, make_constraint_info, inject_pivot, \
    eval_spline, eval_spline_samples

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm


from pyuvdata import UVBeam

from functools import partial
import pickle

import arviz as az

def to_lin(db):
    """
    Convert a quantity in decibels to linear units.
    """
    return 10**(db/10)

def get_bin_cent(edges):
    """
    Utility function for histogram bins that given edges returns bin centers.
    """
    return (edges[:-1] + edges[1:]) / 2

def altaz_to_rad(az_deg, alt_deg):
    """
    Convert azimuth and altitude in degrees to azimuth and zenith angle in radians.
    Switches from 'North to East' convention to 'East to North' convention to
    be consistent with pyuvdata.
    """

    az_rad = (np.deg2rad(90. - az_deg)) % (2 * np.pi) # Difference in convention b/w data and UVBeam
    za_rad = np.deg2rad(90. - alt_deg)
    
    return az_rad, za_rad

def get_res(az_rad, za_rad, power_db):
    """
    Get the residual beam in dB (so the log of the ratio) between the measured
    power and the simulated model power at the az_deg and alt_deg in question.
    """
    # FIXME: Uses beam as global
    model_beam = beam.interp(az_array=az_rad, za_array=za_rad)[0][0,pol_ind,0].real # It's a copol power beam
    model_beam_db = 10 * np.log10(model_beam) 
    assert np.all(np.isfinite(model_beam_db))
    
    res = power_db - model_beam_db
    return res, model_beam

# Good candidate for some util functions for the module -- should get rid of global variable dependence
def get_pass_subset(passes, slc=slice(None), db_cut=3):
    """
    Get a pass subset including the altitudes, azimuths, and power measurements in dB.

    Args:
        passes (list):
            List of satellite pass records.
        slc (slice):
            A slice object for slicing into passes.
        db_cut (float):
            A lower cutoff such that if any power measurements are above it, 
            the corresponding pass data are excluded from the returned items.
    Returns:
        alt_deg (array):
            Stacked array of altitudes in degrees for selected passes.
        az_deg (array):
            Stacked array of azimuths (North to East) in degrees for selected passes.
        power_db (array):
            Stacked array of power measurements in dB.
    """
    ret = []
    for attr in ["alt_deg", "az_deg", "power_db"]:
        ret.append(
            stack_pass_attr(
                attr, 
                passes, 
                slc, 
                db_cut=db_cut
            )
        )
    
    return ret

def get_pass_cond(p, db_cut=3):
        # Sometimes there are outlier passes with crazy data -- 3 dB cut gets rid of them for at least a couple instances
        # Some points come in below the horizon
    return (max(p.power_db) < db_cut) and (all(p.alt_deg > 0))

def stack_pass_attr(attr, passes, slc=slice(None), db_cut=3):
    filter_fn = partial(get_pass_cond, db_cut=db_cut)
    filtered_passes = filter(filter_fn, passes[slc])
    return np.concatenate([getattr(p, attr) for p in filtered_passes])

def plot_pass_subset(alt_deg, az_deg, power_db):
    az_rad, za_rad = altaz_to_rad(az_deg, alt_deg)
    
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})

    im = ax.scatter(az_rad, np.sin(za_rad), c=power_db, s=0.1)
    fig.colorbar(im)
    
    return

def plot_pass_ratio(outfile, az_rad, za_rad, res):
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})

    im = plt.scatter(
            az_rad, 
            np.sin(za_rad), 
            c=res, 
            s=0.1, 
            vmin=-30, 
            vmax=30, 
            cmap="coolwarm"
        )
    fig.colorbar(im, label="Data/Sim (dB)")
    ax.grid(False)

    fig.savefig(outfile)
    plt.close(fig)

def prep_spline_params(beam, ortho_knots=False):
    if ortho_knots:
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

def plot_hist(idata, ax, param_name, to_hist, xlabel):
    scale_samples = idata.posterior[param_name]
    scales_to_plot = az.extract(
        idata,
        combined=True,
        var_names=param_name,
        num_samples=30,
    )
    counts, _, _ = ax.hist(
            to_hist,
            histtype="step",
            bins="auto",
            density=True
        )
    mean_scale = scale_samples.mean().values
    xplot = 10 * mean_scale * np.linspace(-1, 1, num=100)
    ax.plot(xplot, 
            np.exp(dist.SoftLaplace(
                loc=0, 
                scale=jnp.array(scales_to_plot)[None, :]
            ).log_prob(xplot[:, None])),
            color="tab:orange",
            alpha=0.1)
    
    ax.set_yscale("log")
    ax.set_ylim([0.5*np.amin(counts[counts > 0]), 2 * np.amax(counts)])
    ax.set_xlabel(xlabel)

    return

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
    from numpyro.infer import MCMC, NUTS, SVI, Trace_ELBO
    from numpyro.infer.autoguide import AutoDelta   
    
    from jax import random
    import jax.numpy as jnp
    import jax
    #jax.config.update("jax_debug_nans", True)

    import arviz as az

    # FIXME: Hardcode -- only zenith
    delays = np.zeros(16, dtype=int)
    delays = np.array([delays, delays])

    beam = UVBeam.from_file(
        args.mwa_beamfile,
        freq_range=[136e6, 138e6],
        delays=delays,
        za_range=(0,90)
        )
    beam.peak_normalize()
    beam.efield_to_power()


    pol_ind_dict = {
        "XX": np.where(beam.polarization_array == -5)[0][0],
        "YY": np.where(beam.polarization_array == -6)[0][0]
    }
    assert args.pol in pol_ind_dict.keys(), "Invalid pol; must be XX or YY"
    pol_ind = pol_ind_dict[args.pol]

    subdir = f"rf{args.rf_num}/S{args.tile}/{args.pol}"
    if args.jackknife:
        subdir += f"/{args.jackknife_mode}_jackknife/half{args.jackknife_half}"
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
        if args.jackknife_half: # implicitly second half
            pass_slice = slice(
                Npass_total // 2, 
                Npass_total
            )
        else: # otherwise first half
            pass_slice = slice(
                0, 
                Npass_total // 2
            )
    elif args.num_pass is None:
        pass_slice = slice(
            0,
            Npass_total
        )
    else:
        pass_slice = slice(
            0,
            args.num_pass
        )

    alt_deg, az_deg, power_db = get_pass_subset(passes, slc=pass_slice)
    az_rad, za_rad = altaz_to_rad(az_deg, alt_deg)
    res, model_beam = get_res(az_rad, za_rad, power_db)


    # Do some Bayesian inference with MCMC (NUTS) on some short chains
    if args.ortho_knots:
        # Basic knot stuff
        t_x, t_y, p, q, n_x, n_y, _= prep_spline_params(beam, ortho_knots=True)
        knot_mask = make_disk_mask(t_x, t_y, p, q)           # (n_x, n_y) bool, numpy
        knot_ij, Nparam = make_flat_index(knot_mask)            # static arrays
        if args.enforce_boresight: # Sample one less parameter to satisfy constraint
            Nparam -= 1
        
        # Gradient stuff
        xgrad_base = np.linspace(-1, 1, num=100)
        ygrad_base = np.copy(xgrad_base)
        Xgrad, Ygrad = np.meshgrid(xgrad_base, ygrad_base, indexing="ij")

        # This mask is not dealing with spline knots, so don't use the helper
        grad_mask = Xgrad**2 + Ygrad**2 <= 1.0
        #Xgrad = Xgrad[grad_mask]
        #Ygrad = Ygrad[grad_mask]


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

    class SkewLogistic(dist.Distribution, metaclass=dist.distribution.DistributionMeta):
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
            data=None, 
            add_noise_bias=False,
            noise_bias=1e-6,
            ortho_knots=False,
            enforce_boresight=False,
    ):
        scale_vals = numpyro.sample("scale_vals", dist.Uniform(low=1e-2, high=1e2))
        scale_grad = numpyro.sample("scale_grad", dist.Uniform(low=1e-2, high=1e2))
        scale_noise = numpyro.sample("scale_noise", dist.Uniform(low=1e-2, high=1e2))
        with numpyro.plate("Nparam", Nparam): # Need to enforce peak-normalization
            surf_vals = numpyro.sample("surf_vals",
                                       dist.SoftLaplace(loc=0, scale=scale_vals))
        spl, grad = get_spl_grad(ortho_knots, enforce_boresight, surf_vals)
        model_vals = eval_spline_batch(spl, dat_coord_1, dat_coord_2) # may be x_data, y_data
        if add_noise_bias:
            model_vals += 10 * jnp.log10(1 + noise_bias / model_beam / to_lin(model_vals))
        
        
        lp1 = dist.SoftLaplace(loc=0, scale=scale_grad).log_prob(grad[0]).sum()
        lp2 = dist.SoftLaplace(loc=0, scale=scale_grad).log_prob(grad[1]).sum()
        numpyro.factor("grad_pot", lp1 + lp2)
        
        
        with numpyro.plate("Ndat", Ndat):
            obs = numpyro.sample("obs",
                                dist.SoftLaplace(loc=model_vals, scale=scale_noise),
                                obs=data)

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
        dat_coord_1[::args.decimation_factor], 
        dat_coord_2[::args.decimation_factor], 
        len(res[::args.decimation_factor])
    )
    model_kwargs = {
        "data": res.real[::args.decimation_factor], 
        "ortho_knots": args.ortho_knots, 
        "enforce_boresight": args.enforce_boresight
    }
    if args.inference: # Need to do inference
        key = random.key(args.key)
        if args.svi:
            guide = AutoDelta(model)
            svi = SVI(model, guide, numpyro.optim.Adam(1e-1), loss=Trace_ELBO())
            svi_result = svi.run(
                key, 
                args.num_sample, 
                *model_args, 
                **model_kwargs
            )
            params = svi_result.params
            plt.plot(svi_result.losses)
            plt.ylabel("Loss")
            plt.xlabel("Iteration")
            plt.savefig(f"{outdir}/losses.pdf")

            with open(f"{outdir}/svi_params.pkl", "wb") as svi_params_file:
                pickle.dump(params, svi_params_file)
        else:
            kernel = NUTS(model)
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
    if args.postprocess: # Already have samples, make some plots
        mean_res, xedges, yedges, bn = binned_statistic_2d(
            dat_coord_1[::5], 
            dat_coord_2[::5], 
            res[::5], 
            bins=np.linspace(-1, 1, num=101), 
            expand_binnumbers=True
        )
        xcent = get_bin_cent(xedges)
        ycent = get_bin_cent(yedges)
        X, Y = np.meshgrid(xcent, ycent, indexing="ij")

        idata = az.from_netcdf(f"{outdir}/mcmc_out.nc")
        num_samp_total = args.ndevice * args.num_sample
        c_samps = jnp.array(idata.posterior["surf_vals"]).reshape(num_samp_total, -1)
        coeffs_for_mod_samps = jnp.zeros([num_samp_total, n_x, n_y])
        for k in range(num_samp_total):
            coeffs_for_mod_samps = coeffs_for_mod_samps.at[k].set(
                scatter_coeffs(
                    c_samps[k],
                    knot_mask, 
                    knot_ij
                )
            )
        spl = SphericalSpline(
            t_x, t_y, coeffs_for_mod_samps, p , q
        )


        res_samples = eval_spline_samples(
            spl, 
            X.flatten(), 
            Y.flatten()
        ).reshape(
            num_samp_total, 
            *X.shape
        )

        
        r_interp = np.sqrt(X**2 + Y**2)
        za_interp = np.arcsin(r_interp)
        az_interp = np.arctan2(Y, X) % (2 * np.pi)
        interp_mask = r_interp <= 1.0

        bm_interp = np.zeros(za_interp.shape)
        bm_interp[interp_mask] = beam.interp(
            za_array=za_interp[interp_mask], 
            az_array=az_interp[interp_mask]
        )[0][0, pol_ind, 0]
        bm_interp[~interp_mask] = np.nan

        # Sample beams in linear units
        beam_samps = to_lin(res_samples) * bm_interp
        beam_from_mean = to_lin(mean_res) * bm_interp

        mean_beam = beam_samps.mean(axis=0)
        fig, ax = plt.subplots(
            nrows=4, sharex=True, sharey=True, figsize=(6, 12)
        )
        mean_im = ax[0].pcolormesh(
            X, Y, mean_beam, cmap="plasma", vmin=0, vmax=1
        )
        fig.colorbar(mean_im, ax=ax[0], label="Beam Value")

        mb_minus_interp = mean_beam - bm_interp
        lin_res_im = ax[1].pcolormesh(
            X,
            Y,
            mb_minus_interp,
            cmap="Spectral",
            norm=SymLogNorm(vmin=-1, vmax=1, linthresh=1e-3)
        )
        fig.colorbar(lin_res_im, ax=ax[1], label="Residual Beam")

        lin_z_im = ax[2].pcolormesh(
            X,
            Y,
            mb_minus_interp / beam_samps.std(axis=0),
            cmap="coolwarm",
            vmin=-30,
            vmax=30,
        )
        fig.colorbar(lin_z_im, ax=ax[2], label="Residual Beam z-score")

        mean_res_post = res_samples.mean(axis=0)
        mean_res_post = mean_res_post.at[~interp_mask].set(jnp.nan)
        ratio_im = ax[3].pcolormesh(
            X, 
            Y, 
            mean_res_post, 
            vmin=-25, 
            vmax=25, 
            cmap="coolwarm"
        )
        fig.colorbar(ratio_im, ax=ax[3], label="Fit Ratio (dB)")

        ax[-1].set_xlabel("EW Direction Cosine")
        for ax_ob in ax:
            ax_ob.set_ylabel("NS Direction Cosine")

        
        fig.tight_layout()
        fig.savefig(f"{outdir}/mean_beam_check.pdf", bbox_inches="tight")
        plt.close(fig)

        # Compare to inferred distributions using the scale parameters
        fig, ax = plt.subplots(nrows=3, figsize=(6, 4))
        # Look at inferred coefficient distributions
        mean_surf_vals = jnp.array(idata.posterior["surf_vals"]).mean(axis=(0,1))
        plot_hist(
            idata, 
            ax[0], 
            "scale_vals",
            mean_surf_vals,
            "Coefficient Values"
        )
        mean_spl, mean_grad = get_spl_grad(
            args.ortho_knots,
            args.enforce_boresight,
            mean_surf_vals
        )
        plot_hist(
            idata, 
            ax[1], 
            "scale_grad",
            jnp.concatenate(mean_grad, axis=0).flatten(),
            "Coefficient Finite Differences"
        )
        mean_mod_for_data = eval_spline_batch(
            mean_spl, 
            *model_args[:2]
        )
        plot_hist(
            idata, 
            ax[2], 
            "scale_noise",
            model_kwargs["data"] - mean_mod_for_data,
            "Data Residuals"
        )
        fig.savefig(f"{outdir}/dist_plot.pdf")
        plt.close(fig)

        # Compare to 10 random passes
        # FIXME: hardcode db_cut
        fig, ax = plt.subplots()
        spl_for_passes = SphericalSpline(
            t_x, t_y, coeffs_for_mod_samps[::num_samp_total//30], p , q
        )
        
        lengths = [len(p.power_db) for p in passes if np.amax(p.power_db) < 3]
        pass_boundaries = np.cumsum(lengths)
        pass_boundaries = np.append(0, pass_boundaries)
        pass_inds = np.random.choice(len(lengths) - 1, 10, replace=False)
        xstart = 0
        for trial_ind, pass_ind in enumerate(pass_inds):
            pass_slc = slice(
                pass_boundaries[pass_ind], 
                pass_boundaries[pass_ind + 1]
            )
            xvals = np.arange(xstart, xstart + lengths[pass_ind])
            xstart += lengths[pass_ind]
            ax.plot(
                xvals,
                res.real[pass_slc], marker=".", linestyle="none", color=f"C{trial_ind}"
            )

            pass_mod = eval_spline_samples(
                spl_for_passes,
                dat_coord_1[pass_slc],
                dat_coord_2[pass_slc]
            )

            ax.plot(xvals, pass_mod.T, color=f"C{trial_ind}", alpha=0.2)
        ax.set_ylabel("Residual (dB)")
        ax.set_xlabel("Time Index")

        # 10 random linear beam samples

        # Slice plots


        



