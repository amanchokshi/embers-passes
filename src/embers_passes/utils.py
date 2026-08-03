import yaml
import argparse
import random
import itertools

import arviz as az
import numpy as np
from scipy.stats import median_abs_deviation as mad
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
        "jackknife_chunk": 0,
        "jackknife_key": 224007541,
        "ortho_knots": False,
        "svi": False,
        "enforce_boresight": False,
        "ds_factor": 1,
        "inference": True,
        "postprocess": False,
        "pass_plot_seed": 321906282,
        "chunk_sec": 1728000, # 20 days,
        "dense_mass": False
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
    Nchunk: int = 4,
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
        product_args.append(list(range(Nchunk)))

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
            jk_mode, jk_chunk = param_tuple[2:]
            jk_key =  random.randint(0, JAX_KEY_MAX)
            config["jackknife"] = True
            config["jackknife_mode"] = jk_mode
            config["jackknife_chunk"] = jk_chunk
            config["jackknife_key"] = jk_key

            filename += f"jk{jk_mode}_half{jk_chunk}"

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
    mad_maps
        Count-weighted RMS of per-pass MADs per pixel (i.e. sqrt of the
        weighted average of per-pass MAD-squared values, weighted by each
        pass's sample count in that pixel), scaled to match standard
        deviation of normal distribution, then divided by additional
        sqrt(count_map) to match expected rms of averaged i.i.d. noise. 
        Unpopulated pixels set to nan.
    median_sem_maps
        MAD of the per-pass medians per pixel, divided by sqrt(n), where n
        is the number of passes contributing to that pixel. Unpopulated
        pixels set to nan. Captures variation of per-pass correlated systematics 
        in final map.
    """

    if decimate < 1:
        raise ValueError("decimate must be >= 1.")


    median_maps: list[np.ndarray] = []
    count_maps: list[np.ndarray] = []
    mad_maps: list[np.ndarray] = []
    median_sem_maps: list[np.ndarray] = []

    for chunk in chunks:
        median_map, count_map, mad_map, median_sem_map = bin_chunk(
            passes, 
            chunk, 
            decimate=decimate, 
            nside=nside
        )

        median_maps.append(median_map)
        count_maps.append(count_map)
        mad_maps.append(mad_map)
        median_sem_maps.append(median_sem_map)

    return median_maps, count_maps, mad_maps, median_sem_maps

def bin_chunk(
        passes, 
        chunk,
        decimate: int = 1, 
        nside: int = 32, 
        filter: bool = False,
        outlier_threshold: float = 3., 
):
    """
    Bin a particular chunk from passes.

    Parameters:
        passes
            List of pass records.
        chunk
            List of indices into passes indicating the chunk
        decimate
            Keep every ``decimate``-th valid sample from each pass. A value of 1
            retains all valid samples. Default is 1.
        nside
            HEALPix NSIDE parameter. Default is 32.
        filter
            Whether to trim outliers before binning (default False)
        outlier_threshold:
            Multiple of MAD to throw away
    Returns:
        median_map
            Median power per pixel, pooled over all raw samples from
            contributing passes, with unpopulated pixels set to nan
        count_map
            Total number of samples (summed across contributing passes) per pixel
        mad_map
            Count-weighted RMS of per-pass MADs per pixel (i.e. sqrt of the
            weighted average of per-pass MAD-squared values, weighted by each
            pass's sample count in that pixel), scaled to match standard
            deviation of normal distribution, then divided by additional
            sqrt(count_map) to match expected rms of averaged i.i.d. noise. 
            Unpopulated pixels set to nan.
        median_sem_map
            MAD of the per-pass medians per pixel, divided by sqrt(n), where n
            is the number of passes contributing to that pixel. Unpopulated
            pixels set to nan. Captures variation of per-pass correlated 
            systematics in final map.
    """
    if not hp.isnsideok(nside):
        raise ValueError(f"Invalid HEALPix NSIDE: {nside}.")
    npix = hp.nside2npix(nside)

    # Per-pixel lists, one entry per pass that contributed samples to that
    # pixel: raw sample arrays (for the pooled global median), plus each
    # pass's median, MAD, and sample count.
    pass_raw_values: list[list[np.ndarray]] = [[] for _ in range(npix)]
    pass_median_values: list[list[float]] = [[] for _ in range(npix)]
    pass_mad_values: list[list[float]] = [[] for _ in range(npix)]
    pass_count_values: list[list[int]] = [[] for _ in range(npix)]

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

        # Accumulate this pass's samples per-pixel before reducing to a
        # single median/MAD/count per pixel for this pass.
        pass_power_values: dict[int, list[float]] = {}
        for pixel, power in zip(
                pixels,
                power_valid,
            ):
            pass_power_values.setdefault(pixel, []).append(float(power))

        for pixel, values in pass_power_values.items():
            values_arr = np.asarray(values, dtype=float)
            pass_median = np.median(values_arr)
            pass_mad = mad(
                values_arr,
                scale="normal",
                nan_policy="omit",
            )
            pass_raw_values[pixel].append(values_arr)
            pass_median_values[pixel].append(pass_median)
            pass_mad_values[pixel].append(pass_mad)
            pass_count_values[pixel].append(values_arr.size)

    median_map = np.full(npix, np.nan, dtype=float)
    count_map = np.zeros(npix, dtype=int)
    mad_map = np.full(npix, np.nan, dtype=float)
    median_sem_map = np.full(npix, np.nan, dtype=float)

    for pixel, values in enumerate(pass_median_values):
        raw_arrays = pass_raw_values[pixel]
        mads = pass_mad_values[pixel]
        counts = pass_count_values[pixel]

        if filter:
            if not values:
                continue

            values = np.asarray(values, dtype=float)
            mads = np.asarray(mads, dtype=float)
            counts = np.asarray(counts, dtype=int)

            pixel_median = np.median(values)
            pixel_mad = mad(
                values,
                scale="normal",
                nan_policy="omit",
            )

            if pixel_mad == 0.0:
                keep = values == pixel_median
            else:
                keep = (
                    np.abs(values - pixel_median)
                    <= outlier_threshold * pixel_mad
                )

            filtered_values = values[keep]
            filtered_mads = mads[keep]
            filtered_counts = counts[keep]
            filtered_raw_arrays = [
                arr for arr, k in zip(raw_arrays, keep) if k
            ]

            if filtered_values.size == 0:
                continue

            filtered_pooled_raw = np.concatenate(filtered_raw_arrays)
            filtered_global_median = np.median(filtered_pooled_raw)
            filtered_mad = np.sqrt(
                np.average(filtered_mads**2, weights=filtered_counts)
            )
            filtered_count = int(filtered_counts.sum())
            filtered_median_mad = mad(
                filtered_values,
                scale="normal",
                nan_policy="omit",
            )

            median_map[pixel] = filtered_global_median
            mad_map[pixel] = filtered_mad
            count_map[pixel] = filtered_count
            median_sem_map[pixel] = filtered_median_mad / np.sqrt(filtered_values.size)
        else:
            counts = np.asarray(counts, dtype=int)
            count_map[pixel] = int(counts.sum())

            if values:
                mads = np.asarray(mads, dtype=float)
                pooled_raw = np.concatenate(raw_arrays)
                median_map[pixel] = np.median(pooled_raw)
                # Extra counts factor in denominator since noise is independent _per sample_
                mad_map[pixel] = np.sqrt(
                    np.average(mads**2, weights=counts) / count_map[pixel]
                )
                # Variation expected from different constant _per pass_
                median_of_medians_mad = mad(
                    values,
                    scale="normal",
                    nan_policy="omit",
                )
                median_sem_map[pixel] = median_of_medians_mad / np.sqrt(len(values))
    return median_map, count_map, mad_map, median_sem_map

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