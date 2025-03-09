#!/usr/bin/env python3
"""
Black Body Radiation Spectrum Analysis

This script analyzes experimental data from a black body radiation spectrum experiment
and fits it to Planck's law to determine the temperature of the black body.
"""

from typing import Optional, Tuple, List, Dict, Any
import numpy as np
from numpy.typing import NDArray
import pandas as pd
import matplotlib.pyplot as plt
from scipy import constants
import os
from utils.odr_fit import perform_odr_fit, plot_odr_fit

datafile = os.path.join(
    os.path.dirname(__file__), "data", "black_body_radiation_spectrum.xlsx"
)


def main() -> None:
    """Main function to run the analysis"""
    ############################################################
    # ROTATION RATIO
    ############################################################
    # 1) Load the sheet
    rotation_df = pd.read_excel(datafile, sheet_name="rotation_ratio")

    # 2) Extract the single column (dropping any NaNs, if they exist)
    col_name = "small_circle_angle_change_per_large_circle_full_rotation-rad"
    rotation_measurements_in_rad: NDArray[np.float64] = rotation_df[col_name].values

    # 3) Compute mean + error (with std dev and SEM)
    angle_change_per_large_circle_full_rotation_in_rad: float = np.mean(rotation_measurements_in_rad)
    angle_change_per_large_circle_full_rotation_error_in_rad: float = np.std(rotation_measurements_in_rad, ddof=1) / np.sqrt(len(rotation_measurements_in_rad))

    # 4) Convert to dimensionless ratio by dividing by 2π
    ANGLE_RATIO: float = angle_change_per_large_circle_full_rotation_in_rad / (2 * np.pi)  # final ratio (dimensionless)
    ANGLE_RATIO_ERROR: float = angle_change_per_large_circle_full_rotation_error_in_rad / (2 * np.pi)  # uncertainty (SEM) on that ratio

    print("\nRotation Ratio Analysis:")
    print(f"    Mean small-circle angle per large-circle revolution = {angle_change_per_large_circle_full_rotation_in_rad:.3f} rad")
    print(f"    Error on mean = {angle_change_per_large_circle_full_rotation_error_in_rad:.3f} rad")
    print(f"    Dimensionless ratio = {ANGLE_RATIO:.5f} ± {ANGLE_RATIO_ERROR:.5f}")

    ############################################################
    # STARTING ANGLE
    ############################################################
    # 1) Load the sheet
    starting_angle_df: pd.DataFrame = pd.read_excel(datafile, sheet_name="starting_angle")
    
    # 2) Extract the single column
    col_name: str = "starting_angle-rad"
    unadjusted_starting_angle_measurements_in_rad: NDArray[np.float64] = starting_angle_df[col_name].values
    
    # 3) Compute mean + error (with std dev and SEM)
    unadjusted_starting_angle_in_rad: float = np.mean(unadjusted_starting_angle_measurements_in_rad)
    unadjusted_starting_angle_error_in_rad: float = \
        np.std(unadjusted_starting_angle_measurements_in_rad, ddof=1) / np.sqrt(len(unadjusted_starting_angle_measurements_in_rad))
    
    # 4) Multiply by the rotation ratio to get the adjusted angle
    STARTING_ANGLE_IN_RAD: float = unadjusted_starting_angle_in_rad / ANGLE_RATIO
    
    # 5) Propagate errors (using the formula for division of independent variables)
    # For z = x / y, the relative error is: (Δz/z)² = (Δx/x)² + (Δy/y)²
    # So Δz = z * sqrt((Δx/x)² + (Δy/y)²)
    relative_error_squared: float = (unadjusted_starting_angle_error_in_rad/unadjusted_starting_angle_in_rad)**2 + (ANGLE_RATIO_ERROR/ANGLE_RATIO)**2
    STARTING_ANGLE_ERROR_IN_RAD: float = STARTING_ANGLE_IN_RAD * np.sqrt(relative_error_squared)
    
    print("\nStarting Angle Analysis:")
    print(f"    Mean starting angle (unadjusted - measured in units of the small circle) = {unadjusted_starting_angle_in_rad:.3f} rad")
    print(f"    Standard error (unadjusted) = {unadjusted_starting_angle_error_in_rad:.3f} rad")
    print(f"    Starting angle (unadjusted) = {unadjusted_starting_angle_in_rad:.3f} ± {unadjusted_starting_angle_error_in_rad:.3f} rad")
    print(f"    Starting angle (adjusted) = {STARTING_ANGLE_IN_RAD:.3f} ± {STARTING_ANGLE_ERROR_IN_RAD:.3f} rad")

    ############################################################
    # RESISTANCE MEASUREMENTS
    ############################################################
    # 1) Load the sheet
    resistance_df: pd.DataFrame = pd.read_excel(datafile, sheet_name="lightbulb_with_wires")
    
    # 2) Extract the columns
    voltage_in_mV: NDArray[np.float64] = resistance_df["voltage-mV"].values
    current_in_mA: NDArray[np.float64] = resistance_df["current-mA"].values
    
    # 3) Compute resistance and error using ODR fit
    
    # We assume 0.01 mA/mV uncertainty in current and voltage measurements.
    voltage_error_in_mV: NDArray[np.float64] = np.ones_like(voltage_in_mV) * 0.01
    current_error_in_mA: NDArray[np.float64] = np.ones_like(current_in_mA) * 0.01

    # Perform ODR fit (slope = resistance in Ohms since V/I = R).
    R_total_in_ohms, R_total_error_in_ohms, voltage_offset_in_mV, offset_error_in_mV = perform_odr_fit(
        current_in_mA, voltage_in_mV, current_error_in_mA, voltage_error_in_mV
    )

    # 4) Compute the resistance of the wire
    R_LIGHTBULB_IN_OHMS = 1.28
    R_LIGHTBULB_ERROR_IN_OHMS = 0.01
    R_WIRES_IN_OHMS = R_total_in_ohms - R_LIGHTBULB_IN_OHMS
    R_WIRES_ERROR_IN_OHMS = np.sqrt(R_total_error_in_ohms**2 + R_LIGHTBULB_ERROR_IN_OHMS**2)

    # 4) Print results
    print("\nResistance Analysis:")
    print(f"    Resistance = {R_total_in_ohms:.3f} ± {R_total_error_in_ohms:.3f} Ω")
    if abs(voltage_offset_in_mV) > 1e-3:  # Only show offset if it's non-zero
        print(f"    Offset voltage = {voltage_offset_in_mV:.3f} ± {offset_error_in_mV:.3f} mV")
    print(f"    Resistance of the lightbulb = {R_LIGHTBULB_IN_OHMS:.3f} ± {R_LIGHTBULB_ERROR_IN_OHMS:.3f} Ω")
    print(f"    Resistance of the wire = {R_WIRES_IN_OHMS:.3f} ± {R_WIRES_ERROR_IN_OHMS:.3f} Ω")
    
    # 5) Create plot
    fig = plot_odr_fit(
        current_in_mA, voltage_in_mV, 
        R_total_in_ohms, voltage_offset_in_mV, 
        R_total_error_in_ohms, offset_error_in_mV,
        x_label="Current (mA)", 
        y_label="Voltage (mV)",
        title="Resistance Measurement (Four-Wire Method)"
    )
    
    # Save the plot
    plt.savefig("resistance_fit.png", dpi=300, bbox_inches='tight')
    plt.close(fig)


if __name__ == "__main__":
    main()
