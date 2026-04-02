## A guide to reproduce this data

Resurrecting `EMBERS` has been an archaeological challenge, and required
dragging the codebase kicking and screaming into the modern era. None of the
original `py38` code worked on my M4 MacBook Pro, and I could not even install
the legacy HDF5 stack cleanly.

The [`EMBERS:1.0.1`](https://github.com/amanchokshi/EMBERS/tree/1.0.1) branch
now supports `py313`, and the new `Dockerfile` builds a containerised,
semi-modern version available at
[`amanchokshi/embers:1.0.1-py313`](https://hub.docker.com/repository/docker/amanchokshi/embers/general),
which can also be used on HPC systems.

Assuming a typical `EMBERS` pipeline has already been run and the usual data
products are available
(see [EMBERS by Example](https://embers.readthedocs.io/en/latest/embersbyexample.html)),
the Python and Slurm scripts in this directory should let you reproduce the
data products archived in this repository.

In particular:

- `extract_tile_passes.py`
- `sbatch_extract_tile_passes.sh`
