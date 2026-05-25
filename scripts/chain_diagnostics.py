import arviz as az
import sys
from embers_passes import run_diagnostics

# First arg to python script is a path to a netcdf file
idata = az.from_netcdf(sys.argv[1])

run_diagnostics(idata)

