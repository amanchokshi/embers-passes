import yaml
import argparse
import random
import arviz as az

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
        "jackknife_key": 224007541
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
    valid_jk_modes = ["time", "random", "full"]
    selected_mode = config.get("jackknife_mode")
    if selected_mode not in valid_jk_modes:
        raise ValueError(f"jackknife_mode {selected_mode} is invalid; must be in {valid_jk_modes}")

    return argparse.Namespace(**config)

def write_configs(
    base_config_path: str,
    outdir: str,
    tiles: list = None,
    pols: list = None,
) -> None:
    """
    Write one YAML config file per (tile, pol) combination.
    
    Args:
        base_config_path: path to yaml of shared settings 
            (existing values for tile, pol, key will be overwritten in new yaml)
        tiles:       list of tile strings e.g. ["08", "09", "10"]
        pols:        list of polarisations e.g. ["XX", "YY"]
        outdir:      directory to write config files into
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

    for tile in tiles:
        for pol in pols:
            key = random.randint(0, JAX_KEY_MAX) 
            config = {
                **base_config,
                "tile": tile,
                "pol": pol,
                "key": key,
            }
            filename = f"{outdir}/config_tile{tile}_{pol}.yaml"
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