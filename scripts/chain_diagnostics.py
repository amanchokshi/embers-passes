import arviz as az
import sys

# First arg to python script is a path to a netcdf file
idata = az.from_netcdf(sys.argv[1])

has_errors, diagnostics = az.diagnose(idata, return_diagnostics=True)

if has_errors:
    print("Some tests failed!")
    print(diagnostics)
else:
    print("All tests passed!")

