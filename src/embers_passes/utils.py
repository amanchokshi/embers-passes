import yaml
import argparse
import random

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
        "num_pass": 1000,
        "ndevice": 4,
        "key": 338422580,
        "num_warmup": 1000,
        "num_sample": 2000,
        "mix_model": False,
    }

    with open(config_path, "r") as f:
        user_config = yaml.safe_load(f)

    config = {**defaults, **user_config}

    # Validate required fields
    required = ["mwa_beamfile", "outdir"]
    missing = [field for field in required if config.get(field) is None]
    if missing:
        raise ValueError(f"Missing required config fields: {missing}")

    return argparse.Namespace(**config)

def write_configs(
    base_config: dict,
    outdir: str,
    tiles: list = None,
    pols: list = None,
) -> None:
    """
    Write one YAML config file per (tile, pol, key) combination.
    
    Args:
        base_config: dict of shared settings 
            (existing values for tile, pol, key will be overwritten in new yaml)
        tiles:       list of tile strings e.g. ["08", "09", "10"]
        pols:        list of polarisations e.g. ["XX", "YY"]
        num_keys:    how many random keys to generate per (tile, pol) pair
        outdir:      directory to write config files into
    """
    if tiles is None:
        # All the tiles in the rf0 directory, listed as ints
        tiles = list(range(6, 11)) + [12] + list(range(29, 37))
    # convert to strings with leading 0
    tiles = [f"{t:02d}" if isinstance(t, int) else t for t in tiles]
    pols = pols if pols is not None else ["XX", "YY"]

    num_keys = len(tiles) * len(pols)
    
    JAX_KEY_MAX = 2**32 - 1
    keys = [random.randint(0, JAX_KEY_MAX) for _ in range(num_keys)]

    for tile in tiles:
        for pol in pols:
            for key in keys:
                config = {
                    **base_config,
                    "tile": tile,
                    "pol": pol,
                    "key": key,
                }
                filename = f"{outdir}/config_tile{tile}_{pol}_key{key}.yaml"
                with open(filename, "w") as f:
                    yaml.dump(config, f, default_flow_style=False)