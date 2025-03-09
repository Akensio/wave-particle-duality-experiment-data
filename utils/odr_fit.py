"""
Simple ODR Fit Utilities

A clear, straightforward implementation of Orthogonal Distance Regression (ODR)
for fitting data with uncertainties in both x and y variables.
"""

from typing import Tuple, Optional
import numpy as np
from numpy.typing import NDArray
from scipy import odr
import matplotlib.pyplot as plt
from matplotlib.figure import Figure


def linear_function(
    params: NDArray[np.float64], x: NDArray[np.float64]
) -> NDArray[np.float64]:
    """
    Simple linear function: y = m*x + b

    Args:
        params: List containing [slope, intercept]
        x: Independent variable values

    Returns:
        Dependent variable values calculated from the model
    """
    slope, intercept = params
    return slope * x + intercept


def perform_odr_fit(
    x: NDArray[np.float64],
    y: NDArray[np.float64],
    x_err: Optional[NDArray[np.float64]] = None,
    y_err: Optional[NDArray[np.float64]] = None,
) -> Tuple[float, float, float, float]:
    """
    Perform Orthogonal Distance Regression to fit a linear model.

    Args:
        x: Independent variable data as numpy array
        y: Dependent variable data as numpy array
        x_err: Uncertainties in x as numpy array (optional)
        y_err: Uncertainties in y as numpy array (optional)

    Returns:
        tuple: (slope, slope_error, intercept, intercept_error)
    """
    # Create ODR model and data objects
    linear_model = odr.Model(linear_function)

    # Handle the case when errors are not provided
    if x_err is None:
        x_err = np.zeros_like(x)  # Default to zero uncertainty

    if y_err is None:
        y_err = np.zeros_like(y)  # Default to zero uncertainty

    data = odr.RealData(x, y, sx=x_err, sy=y_err)

    # Initial parameter guess: [slope, intercept]
    initial_slope = (y.max() - y.min()) / (x.max() - x.min())
    initial_intercept = np.mean(y) - initial_slope * np.mean(x)

    # Create and run ODR
    odr_instance = odr.ODR(data, linear_model, beta0=[initial_slope, initial_intercept])
    result = odr_instance.run()

    # Extract results
    slope: float = result.beta[0]
    intercept: float = result.beta[1]
    slope_error: float = result.sd_beta[0]
    intercept_error: float = result.sd_beta[1]

    return slope, slope_error, intercept, intercept_error


def plot_odr_fit(
    x: NDArray[np.float64],
    y: NDArray[np.float64],
    slope: float,
    intercept: float,
    slope_error: float,
    intercept_error: float,
    x_err: Optional[NDArray[np.float64]] = None,
    y_err: Optional[NDArray[np.float64]] = None,
    x_label: str = "X",
    y_label: str = "Y",
    title: str = "ODR Fit",
    save_path: Optional[str] = None,
) -> Figure:
    """
    Plot data with ODR fit line.

    Args:
        x: Independent variable data as numpy array
        y: Dependent variable data as numpy array
        slope: Fitted slope
        intercept: Fitted intercept
        slope_error: Error in slope
        intercept_error: Error in intercept
        x_err: Uncertainties in x as numpy array (optional)
        y_err: Uncertainties in y as numpy array (optional)
        x_label: Label for x-axis
        y_label: Label for y-axis
        title: Plot title
        save_path: Path to save the plot (if None, plot is displayed)

    Returns:
        matplotlib.figure.Figure: The created figure
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot data points with error bars
    ax.errorbar(
        x, y, xerr=x_err, yerr=y_err, fmt="o", label="Data", color="blue", markersize=5
    )

    # Plot fit line
    x_range = np.linspace(min(x), max(x), 100)
    y_fit = slope * x_range + intercept
    ax.plot(x_range, y_fit, "r-", label=f"Fit: {y_label} = {slope:.4f} * {x_label} + {intercept:.4f}")

    # Add labels and title
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.grid(True, linestyle="--", alpha=0.7)
    ax.legend()

    # Add fit parameters as text
    text = f"Slope = {slope:.4f} ± {slope_error:.4f}\nIntercept = {intercept:.4f} ± {intercept_error:.4f}"
    ax.text(
        0.05,
        0.95,
        text,
        transform=ax.transAxes,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )

    # Save or show plot
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")

    return fig
