# faerie
Gaussian process regression code in Fortran

Capablities:

- Object orientated allowing for multiple GP instances
- Complex input & output warpings
- Custom mean functions
- Exponential, Matern, and exponentiated quadratic kernels, and summations thereof
- Easily extensible with user custom kernels, input/output warpings, and mean functions
- Mean & variance predictions for single outputs, with necessary conversion and reversion between original and warped spaces
- All necessary hyperparameters and datasets specified via netCDF input file

Lacks:

- Optimisation of hyperparameters; this is best handled by external python libraries (such as PyMC) and interfaced with faerie via the netCDF input file
- Multi-output predictions

## netCDF file format

Global attributes:

- kernel (GP kernel selection; string)
- noise (GP Gaussian noise; float)
- yconrev (output transform selection; string)
- mean (mean function selection; string)

Dimensions:

- inputs (number of input dimensions; int)
- lscales (number of kernel lengthscales; int)
- varns (number of kernel variances; int)
- outputs (number of output dimensions; int) -- ***must currently be set to 1***
- samples (number of data samples; int)
- xcargs (number of input warping hyperparameters; int)
- ycargs (number of output warping hyperparameters; int)
- mcoeffs (number of mean function hyperparameters; int)
- maxstr (upper bound string length; int)

Variables:

- lengthscales (kernel lengthscales; dimensions[lscales]; floats)
- variances (kernel variances; dimensions[varns]; floats)
- input_samples (input data; dimensions[samples,inputs]; floats)
- output_samples (output data; dimensions[samples,outputs]; floats)
- xconrevs (list of input transformations; dimensions[inputs, maxstr]; strings)
- xnorms (list of logicals determining whether to initially normalise x data to the range [0,1] with a safety factor at the bounds, by an affine transform based on the max and min of the data; dimensions[inputs]; ints (0 or 1))
- xconrevargs (input warping hyperparameters; dimensions[xcargs]; floats)
- yconrevargs (output warping hyperparameters; dimensions[ycargs]; floats)
- meancoeffs (mean function hyperparameters; dimensions[mcoeffs]; floats) 
