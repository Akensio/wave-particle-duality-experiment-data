# Common Utilities

This directory contains shared utility functions used by both the simulation and experiment components of the wave-particle duality project.

## Contents

### ODR Fitting Utilities

`odr_fit.py` - Provides functions for Orthogonal Distance Regression (ODR) fitting and visualization.

Key features:
- `perform_odr_fit()` - Performs linear regression accounting for errors in both x and y variables
- `plot_odr_fit()` - Creates visualization of fitted data with error bars and confidence intervals

ODR fitting is especially important for experimental physics analysis where both the independent and dependent variables have measurement uncertainties.

## Usage

Import these utilities in your analysis scripts:

```python
# For experiment scripts
from common.odr_fit import perform_odr_fit, plot_odr_fit

# Usage example
slope, slope_error, intercept, intercept_error = perform_odr_fit(
    x_data, y_data, x_error, y_error
)

fig = plot_odr_fit(
    x_data, y_data, slope, intercept, slope_error, intercept_error,
    x_err=x_error, y_err=y_error,
    x_label="X Label", y_label="Y Label", 
    title="Plot Title"
)
```
