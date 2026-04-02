# embers-passes

A catalog of calibrated MWA satellite beam measurements at 137 MHz from
[Chokshi et. al., 2021](https://ui.adsabs.harvard.edu/abs/2021MNRAS.502.1990C/abstract),
and [Chokshi et. al., 2020](https://ui.adsabs.harvard.edu/abs/2020JOSS....5.2629C/abstract)

---

## Data size

This repository contains ~2.3 GB of HDF5 data products.

If you only need the code, you can clone without data:

```bash
git clone --filter=blob:none https://github.com/amanchokshi/embers-passes.git
```

## Installation

#### Basic install: Clone the repository and install locally

```bash
git clone https://github.com/amanchokshi/embers-passes.git
cd embers-passes

python -m venv .venv
source .venv/bin/activate

pip install .
```

#### Using uv (recommended)

```bash
uv venv
source .venv/bin/activate

uv pip install .
```

### Development / exploration / fun (recommended)

This gives you access to fun plotting stuff, the ability to compare to MWA beams,
interactive notebooks, healpix maps and more.

```bash
git clone https://github.com/amanchokshi/embers-passes.git
cd embers-passes

pip install .[dev]

# or

uv pip install .[dev]
```

This includes:

- Jupyter / JupyterLab
- Matplotlib + cmasher colormaps
- mwa-hyperbeam (for beam modelling)
- healpy

## Accessing the data

```python
>>> from embers_passes import PassFile

>>> pf = PassFile("passes/rf0/S06XX_rf0XX_passes.h5")
>>> pf
    PassFile(path=passes/rf0/S06XX_rf0XX_passes.h5, tile='S06XX', ref='rf0XX', pointings=(0, 2, 4, 41))

>>> passes = pf.read_passes(pointing=0)
>>> len(passes)
    1823

>>> p0 = passes[0]
>>> p0
    PassRecord(norad_id=41188, pointing=0, n_samples=490, peak_alt=17.17°)

>>> print(p0.norad_id, p0.pointing, p0.unix_time)
    41188 0 1568247947.5

>>> print(p0.alt_deg.shape, p0.az_deg.shape, p0.power_db.shape)
    (490,) (490,) (490,)
```


## HDF5 file structure

Each file corresponds to a single reference/tile pair (e.g.
S06XX_rf0XX_passes.h5) and contains all extracted satellite passes for that
pair.

The internal layout is:

```code
<file>.h5
├── attrs
│   ├── tile           (str)
│   ├── ref            (str)
│   ├── start_date     (str)
│   ├── stop_date      (str)
│   ├── pointings      (int[]) e.g. [0, 2, 4, 41]
│   └── n_passes_total (int)
│
└── pointings/
    ├── 0/
    │   └── passes/
    │       ├── 000000/
    │       │   ├── attrs
    │       │   │   ├── norad_id (int)
    │       │   │   ├── unix_time (float)
    │       │   │   └── noise_db (float)
    │       │   │
    │       │   ├── alt_deg   (float[N])
    │       │   ├── az_deg    (float[N])
    │       │   └── power_db  (float[N])
    │       │
    │       └── ...
    │
    ├── 2/
    ├── 4/
    └── 41/
```

### Notes

- Passes are grouped by MWA pointing (0, 2, 4, 41)
- Within each pointing, passes are sorted by mean unix time
- All arrays within a pass have equal length N
- Units:
  - alt_deg → degrees
  - az_deg → degrees
  - power_db → calibrated power (dB)
  - noise_db is a per-pass noise estimate

## Reproducing the dataset

See:

```code
embers_extract/
```
