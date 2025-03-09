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
    rotation_measurements_in_rad = rotation_df[col_name].values

    # 3) Compute mean + error (with std dev and SEM)
    angle_change_per_large_circle_full_rotation_in_rad = np.mean(rotation_measurements_in_rad)
    angle_change_per_large_circle_full_rotation_error_in_rad = np.std(rotation_measurements_in_rad, ddof=1) / np.sqrt(len(rotation_measurements_in_rad))

    # 4) Convert to dimensionless ratio by dividing by 2π
    angle_ratio = angle_change_per_large_circle_full_rotation_in_rad / (2 * np.pi)  # final ratio (dimensionless)
    angle_ratio_error = angle_change_per_large_circle_full_rotation_error_in_rad / (2 * np.pi)  # uncertainty (SEM) on that ratio

    print("\nRotation Ratio Analysis:")
    print(f"    Mean small-circle angle per large-circle revolution = {angle_change_per_large_circle_full_rotation_in_rad:.3f} rad")
    print(f"    Std dev = {angle_change_per_large_circle_full_rotation_error_in_rad:.3f} rad")
    print(f"    Dimensionless ratio = {angle_ratio:.5f} ± {angle_ratio_error:.5f}")

    ############################################################
    # STARTING ANGLE
    ############################################################
    # 1) Load the sheet
    starting_angle_df = pd.read_excel(datafile, sheet_name="starting_angle")
    
    # 2) Extract the single column
    col_name = "starting_angle-rad"
    unadjusted_starting_angle_measurements_in_rad = starting_angle_df[col_name].values
    
    # 3) Compute mean + error (with std dev and SEM)
    unadjusted_starting_angle_in_rad = np.mean(unadjusted_starting_angle_measurements_in_rad)
    unadjusted_starting_angle_error_in_rad = np.std(unadjusted_starting_angle_measurements_in_rad, ddof=1) / np.sqrt(len(unadjusted_starting_angle_measurements_in_rad))
    
    # 4) Multiply by the rotation ratio to get the adjusted angle
    starting_angle_in_rad = unadjusted_starting_angle_in_rad * angle_ratio
    
    # 5) Propagate errors (using the formula for multiplication of independent variables)
    # For z = x * y, the relative error is: (Δz/z)² = (Δx/x)² + (Δy/y)²
    # So Δz = z * sqrt((Δx/x)² + (Δy/y)²)
    relative_error_squared = (unadjusted_starting_angle_error_in_rad/unadjusted_starting_angle_in_rad)**2 + (angle_ratio_error/angle_ratio)**2
    starting_angle_error_in_rad = starting_angle_in_rad * np.sqrt(relative_error_squared)
    
    print("\nStarting Angle Analysis:")
    print(f"    Mean starting angle (unadjusted - measured in units of the small circle) = {unadjusted_starting_angle_in_rad:.3f} rad")
    print(f"    Standard error (unadjusted) = {unadjusted_starting_angle_error_in_rad:.3f} rad")
    print(f"    Starting angle (unadjusted) = {unadjusted_starting_angle_in_rad:.3f} ± {unadjusted_starting_angle_error_in_rad:.3f} rad")
    print(f"    Starting angle (adjusted) = {starting_angle_in_rad:.3f} ± {starting_angle_error_in_rad:.3f} rad")

    ############################################################
    # RESISTANCE MEASUREMENTS
    ############################################################
    # 1) Load the sheet
    resistance_df = pd.read_excel(datafile, sheet_name="lightbulb_four_wires")
    
    # 2) Extract the columns
    voltage_in_mV: NDArray[np.float64] = resistance_df["voltage-mV"].values
    current_in_mA: NDArray[np.float64] = resistance_df["current-mA"].values
    
    # 3) Compute resistance and error using ODR fit
    # We assume 1% uncertainty in both current and voltage measurements
    current_err_in_mA: NDArray[np.float64] = current_in_mA * 0.01
    voltage_err_in_mV: NDArray[np.float64] = voltage_in_mV * 0.01
    
    # Perform ODR fit (slope = resistance in Ohms since V/I = R)
    resistance_in_ohms, resistance_error_in_ohms, offset_in_mV, offset_error_in_mV = perform_odr_fit(
        current_in_mA, voltage_in_mV, current_err_in_mA, voltage_err_in_mV
    )
    
    # 4) Print results
    print("\nResistance Analysis:")
    print(f"    Resistance = {resistance_in_ohms:.3f} ± {resistance_error_in_ohms:.3f} Ω")
    if abs(offset_in_mV) > 1e-3:  # Only show offset if it's non-zero
        print(f"    Offset voltage = {offset_in_mV:.3f} ± {offset_error_in_mV:.3f} mV")
    
    # 5) Create plot
    fig = plot_odr_fit(
        current_in_mA, voltage_in_mV, 
        resistance_in_ohms, offset_in_mV, 
        resistance_error_in_ohms, offset_error_in_mV,
        current_err_in_mA, voltage_err_in_mV,
        x_label="Current (mA)", 
        y_label="Voltage (mV)",
        title="Resistance Measurement (Four-Wire Method)"
    )
    
    # Save the plot
    plt.savefig("resistance_fit.png", dpi=300, bbox_inches='tight')
    plt.close(fig)


if __name__ == "__main__":
    main()
