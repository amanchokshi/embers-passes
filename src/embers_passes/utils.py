import yaml
import argparse
import random
import itertools

import arviz as az
import numpy as np
import healpy as hp

def load_config(config_path: str) -> argparse.Namespace:
    """
    Load configuration from a YAML file and return a Namespace object.
    """
    defaults = {
        "mwa_beamfile": None,   # required
        "outdir": None,         # required
        "pol": "XX",
        "tile": "08",
        "rf_num": 0,
        "plotting": False,
        "num_pass": None,
        "ndevice": 4,
        "key": 338422580,
        "num_warmup": 1000,
        "num_sample": 2000,
        "mix_model": False,
        "jackknife": False,
        "jackknife_mode": "time",
        "jackknife_half": 0,
        "jackknife_key": 224007541,
        "ortho_knots": False
    }

    with open(config_path, "r") as f:
        user_config = yaml.safe_load(f)

    config = {**defaults, **user_config}

    # Validate required fields
    required = ["mwa_beamfile", "outdir"]
    missing = [field for field in required if config.get(field) is None]
    if missing:
        raise ValueError(f"Missing required config fields: {missing}")
    
    # Check for valid jackknife modes
    valid_jk_modes = ["time", "random"]
    selected_mode = config.get("jackknife_mode")
    if selected_mode not in valid_jk_modes:
        raise ValueError(f"jackknife_mode {selected_mode} is invalid; must be in {valid_jk_modes}")

    return argparse.Namespace(**config)

def write_configs(
    base_config_path: str,
    outdir: str,
    tiles: list = None,
    pols: list = None,
    jackknives: list = None,
) -> None:
    """
    Write one YAML config file per (tile, pol) combination.
    
    Args:
        base_config_path: path to yaml of shared settings 
            (existing values for tile, pol, key will be overwritten in new yaml)
        tiles:       list of tile strings e.g. ["08", "09", "10"]
        pols:        list of polarisations e.g. ["XX", "YY"]
        outdir:      directory to write config files into
        jackknives:  list of jackknives to do ("time" or "random").
    """
    
    JAX_KEY_MAX = 2**32 - 1
    
    with open(base_config_path, "r") as base_config_file:
        base_config = yaml.safe_load(base_config_file)

    if tiles is None:
        # All the tiles in the rf0 directory, listed as ints
        tiles = list(range(6, 11)) + [12] + list(range(29, 37))
    # convert to strings with leading 0
    tiles = [f"{t:02d}" if isinstance(t, int) else t for t in tiles]
    pols = pols if pols is not None else ["XX", "YY"]

    product_args = [tiles, pols]
    if jackknives is not None:
        product_args.append(jackknives)
        product_args.append([0, 1])

    for param_tuple in itertools.product(*product_args):
        key = random.randint(0, JAX_KEY_MAX)
        pass_plot_seed = random.randint(0, JAX_KEY_MAX) # Same max for numpy
        tile, pol = param_tuple[:2]
        config = {
            **base_config,
            "tile": tile,
            "pol": pol,
            "key": key,
            "pass_plot_seed": pass_plot_seed
        }
        filename = f"{outdir}/config_tile{tile}_{pol}"
        if jackknives is not None:
            jk_mode, jk_half = param_tuple[2:]
            jk_key =  random.randint(0, JAX_KEY_MAX)
            config["jackknife"] = True
            config["jackknife_mode"] = jk_mode
            config["jackknife_half"] = jk_half
            config["jackknife_key"] = jk_key

            filename += f"jk{jk_mode}_half{jk_half}"

        filename += ".yaml"
        with open(filename, "w") as f:
            yaml.dump(config, f, default_flow_style=False)

def run_diagnostics(idata):
    """
    Run arviz diagnostics on an InferenceData object and print useful info.
    """
    has_errors, diagnostics = az.diagnose(idata, return_diagnostics=True)

    if has_errors:
        print("Some tests failed!")
        print(diagnostics)
    else:
        print("All tests passed!")

def chunk_passes(
    passes: list,
    chunk_sec: float,
) -> tuple[list[list[int]], np.ndarray, np.ndarray]:
    """Group pass indices into fixed-duration time chunks.

    The first time bin begins at the earliest ``unix_time`` in ``passes``.
    Passes are assigned to consecutive bins of width ``chunk_sec``.

    Parameters
    ----------
    passes
        List of pass records. Each record must contain ``unix_time`` attributes.
    chunk_sec
        Width of each time bin in seconds. Must be positive.

    Returns
    -------
    chunks
        Nested list of pass indices. Each inner list contains the indices of
        passes belonging to the corresponding time bin.
    bin_edges
        Unix timestamps defining the edges of the time bins. Has length
        ``len(chunks) + 1``.
    bin_centers
        Unix timestamps at the centers of the time bins. Has length
        ``len(chunks)``.
    """
    if chunk_sec <= 0:
        raise ValueError("chunk_sec must be positive.")

    if len(passes) == 0:
        return [], np.array([], dtype=float), np.array([], dtype=float)

    unix_times = np.asarray(
        [pass_record.unix_time for pass_record in passes],
        dtype=float,
    )

    if not np.all(np.isfinite(unix_times)):
        raise ValueError("All pass unix_time values must be finite.")

    start_time = unix_times.min()

    chunk_indices = np.floor(
        (unix_times - start_time) / chunk_sec
    ).astype(int)

    n_chunks = int(chunk_indices.max()) + 1

    bin_edges = start_time + np.arange(n_chunks + 1) * chunk_sec
    bin_centers = bin_edges[:-1] + 0.5 * chunk_sec

    chunks: list[list[int]] = [[] for _ in range(n_chunks)]

    for pass_index, chunk_index in enumerate(chunk_indices):
        pass_record = passes[pass_index]

        chunks[chunk_index].append(pass_index)

    return chunks, bin_edges, bin_centers

def chunks_to_healpix_counts(
    passes,
    chunks: list[list[int]],
    decimate: int = 1,
    nside: int = 32,
) -> tuple[list[np.ndarray], list[np.ndarray]]:
    """Bin pass samples into per-chunk HEALPix median and count maps.

    Each element of ``chunks`` contains indices selecting pass records from
    ``passes``. A separate median power map and sample-count map are
    constructed for each chunk.

    Only finite samples with altitude strictly greater than 0 degrees are
    retained. The remaining samples are optionally decimated before being
    assigned to HEALPix pixels.

    Parameters
    ----------
    passes
        Sequence of pass records. Each record must provide ``alt_deg``,
        ``az_deg``, and ``power_db`` attributes.
    chunks
        Nested list of indices into ``passes``. Each inner list defines one
        time chunk.
    decimate
        Keep every ``decimate``-th valid sample from each pass. A value of 1
        retains all valid samples. Default is 1.
    nside
        HEALPix NSIDE parameter. Default is 32.

    Returns
    -------
    median_maps
        List containing one median power map per chunk. Unpopulated pixels
        are set to NaN.
    count_maps
        List containing one sample-count map per chunk. Each map has shape
        ``(hp.nside2npix(nside),)``.
    """
    if not hp.isnsideok(nside):
        raise ValueError(f"Invalid HEALPix NSIDE: {nside}.")

    if decimate < 1:
        raise ValueError("decimate must be >= 1.")

    npix = hp.nside2npix(nside)

    median_maps: list[np.ndarray] = []
    count_maps: list[np.ndarray] = []

    for chunk in chunks:
        power_values: list[list[float]] = [[] for _ in range(npix)]

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

            alt_valid = alt_deg[valid]
            az_valid = az_deg[valid]
            power_valid = power_db[valid]

            if decimate > 1:
                alt_valid = alt_valid[::decimate]
                az_valid = az_valid[::decimate]
                power_valid = power_valid[::decimate]

            theta = np.radians(90.0 - alt_valid)
            phi = np.radians(az_valid)

            pixels = hp.ang2pix(
                nside,
                theta,
                phi,
            )

            for pixel, power in zip(
                pixels,
                power_valid,
            ):
                power_values[pixel].append(float(power))

        median_map = np.full(npix, np.nan, dtype=float)
        count_map = np.zeros(npix, dtype=int)

        for pixel, values in enumerate(power_values):
            count_map[pixel] = len(values)

            if values:
                median_map[pixel] = np.median(values)

        median_maps.append(median_map)
        count_maps.append(count_map)

    return median_maps, count_maps

def healpix_to_pyuvdata(
    healpix_map: np.ndarray,
    nside: int = 32,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert a HEALPix map to PyUVData coordinates.

    Return the values of all above-horizon pixels together with their
    zenith angles and azimuths in the PyUVData East-to-North convention.
    """
    healpix_map = np.asarray(healpix_map, dtype=float)
    npix = hp.nside2npix(nside)

    pixel_indices = np.arange(npix, dtype=int)

    # Convert HEALPix pixel indices to zenith angle and azimuth in the
    # original North-to-East convention.
    za_rad_ne, az_rad_ne = hp.pix2ang(
        nside,
        pixel_indices,
        nest=False,
    )

    # Retain only pixels whose centres lie above the horizon.
    above_horizon = za_rad_ne < np.pi / 2.0

    pixel_values = healpix_map[above_horizon]
    za_rad_ne = za_rad_ne[above_horizon]
    az_rad_ne = az_rad_ne[above_horizon]

    # Convert azimuth from the North-to-East convention used by the
    # HEALPix map to the East-to-North convention expected by PyUVData.
    az_rad_pyuvdata = (
        np.pi / 2.0 - az_rad_ne
    ) % (2.0 * np.pi)

    # Zenith angle is unchanged by the azimuth-convention conversion.
    za_rad_pyuvdata = za_rad_ne

    return (
        pixel_values,
        za_rad_pyuvdata,
        az_rad_pyuvdata,
    )