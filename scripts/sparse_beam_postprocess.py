import argparse
from itertools import product

import numpy as np
import matplotlib.pyplot as plt

from embers_passes import SphericalSpline, eval_spline_batch, make_knots_from_grid, eval_spline_samples
from sparse_beam_fit import prep_spline_params
from pyuvdata import UVBeam
import arviz as az

import jax
import jax.numpy as jnp

parser = argparse.ArgumentParser()
parser.add_argument("--indir", help="Where the MCMC outputs are stored.",
                    type=str, required=True)
parser.add_argument("--tilefile", help="Path to list of tiles to load in.",
                    type=str, required=True)
parser.add_argument("--mwa-beamfile", dest="mwa_beamfile", type=str,
                    help="Path to MWA beam file.")
args = parser.parse_args()


with open(args.tilefile, "r") as tilefile:
    tilelist = tilefile.read().splitlines().zfill(2) # Add leading zeros if necessary


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

t_theta, t_phi, p, q, n_th, n_ph, Nparam = prep_spline_params(beam)

##### THINGS TO CHECK #####
##--- Peak normalization
##--- Tension b/w jackknives, full run
##--- Smoothness
##--- Hierarchical structure
##### THINGS TO CHECK #####

pols = ["XX", "YY"]
jk_modes = ["random", "time"]
jk_halves = [0, 1]
product_args = (tilelist, pols, jk_modes, jk_halves)

def samps_to_arr(path):
    idata = az.from_netcdf(path)
    return jnp.array(idata.posterior["surf_vals"]).reshape(-1, n_th, n_ph)

# Load in all the outputs
outputs = {}
for tile, pol, jk_mode, jk_half in product(*product_args):

    tiledir = f"{args.indir}/rf0/S{tile}"
    outputs[tile] = {}
    outputs[tile][pol] = {}
    poldir = f"{tiledir}/{pol}"
    outputs[tile][pol]["full_mcmc"] = samps_to_arr(f"{poldir}/mcmc_out.nc")

    jk_dir = f"{poldir}/{jk_mode}_jackknife"
    outputs[tile][pol][jk_mode] = []

    half_dir = f"{jk_dir}/half{jk_half}"
    outputs[tile][pol][jk_mode].append(samps_to_arr(f"{half_dir}/mcmc_out.nc"))

def eval_peak_norm(coeffs):
    spl = SphericalSpline(t_theta, t_phi, coeffs, p, q)
    val = eval_spline_samples(spl, jnp.atleast_1d(0), jnp.atleast_1d(0))
    return val

# Gives a pytree as above but with peak values
peak_norms = jax.tree.map(eval_peak_norm, outputs) 




