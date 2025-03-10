#!/usr/bin/env python3
"""
Black Body Radiation Spectrum Analysis

This script analyzes experimental data from a black body radiation spectrum experiment
and fits it to Planck's law to determine the temperature of the black body.
"""

import os
from typing import Any, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from numpy.typing import NDArray
from scipy import constants

from utils.odr_fit import perform_odr_fit, plot_odr_fit

datafile = os.path.join(os.path.dirname(__file__), "data", "black_body_radiation_spectrum.xlsx")


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
    angle_change_per_large_circle_full_rotation_error_in_rad: float = np.std(rotation_measurements_in_rad, ddof=1) / np.sqrt(
        len(rotation_measurements_in_rad)
    )

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
    unadjusted_starting_angle_error_in_rad: float = np.std(unadjusted_starting_angle_measurements_in_rad, ddof=1) / np.sqrt(
        len(unadjusted_starting_angle_measurements_in_rad)
    )

    # 4) Multiply by the rotation ratio to get the adjusted angle
    STARTING_ANGLE_IN_RAD: float = unadjusted_starting_angle_in_rad / ANGLE_RATIO

    # 5) Propagate errors (using the formula for division of independent variables)
    # For z = x / y, the relative error is: (Δz/z)² = (Δx/x)² + (Δy/y)²
    # So Δz = z * sqrt((Δx/x)² + (Δy/y)²)
    relative_error_squared: float = (unadjusted_starting_angle_error_in_rad / unadjusted_starting_angle_in_rad) ** 2 + (
        ANGLE_RATIO_ERROR / ANGLE_RATIO
    ) ** 2
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
        current_in_mA,
        voltage_in_mV,
        R_total_in_ohms,
        voltage_offset_in_mV,
        R_total_error_in_ohms,
        offset_error_in_mV,
        x_label="Current (mA)",
        y_label="Voltage (mV)",
        title="Resistance Measurement (Four-Wire Method)",
    )

    plt.savefig("resistance_fit.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


    ############################################################
    # MAXIMUM INTENSITY WAVELENGTH
    ############################################################
    # 1) Load the sheet
    max_intensity_df = pd.read_excel(datafile, sheet_name="max_intensity_wavelength")

    # 2) Extract the columns
    voltage_in_V: NDArray[np.float64] = max_intensity_df["voltage-V"].values
    voltage_error_in_V: NDArray[np.float64] = max_intensity_df["voltage_error-V"].values
    current_in_A: NDArray[np.float64] = max_intensity_df["current-A"].values
    current_error_in_A: NDArray[np.float64] = max_intensity_df["current_error-A"].values
    
    # 3) Calculate resistance and its error for each measurement
    resistance_in_ohms: NDArray[np.float64] = voltage_in_V / current_in_A
    
    # Error propagation for R = V/I
    # (ΔR/R)² = (ΔV/V)² + (ΔI/I)²
    relative_R_error_squared: NDArray[np.float64] = (voltage_error_in_V / voltage_in_V)**2 + (current_error_in_A / current_in_A)**2
    resistance_error_in_ohms: NDArray[np.float64] = resistance_in_ohms * np.sqrt(relative_R_error_squared)
    
    # 4) Calculate resistivity using rho = rho0 * ((R - Rw) / R0)
    rho0 = 5.65e-8  # Resistivity at room temperature in ohm*meters
    R0 = R_LIGHTBULB_IN_OHMS  # Resistance of lightbulb at room temperature
    Rw = R_WIRES_IN_OHMS  # Resistance of the wires
    
    resistivity = rho0 * ((resistance_in_ohms - Rw) / R0)
    
    # Error propagation for resistivity
    # For rho = rho0 * ((R - Rw) / R0)
    # We need to propagate errors from R, Rw, and R0
    term_error_squared = (resistance_error_in_ohms / (resistance_in_ohms - Rw))**2 + (R_WIRES_ERROR_IN_OHMS / (resistance_in_ohms - Rw))**2 + (R_LIGHTBULB_ERROR_IN_OHMS / R0)**2
    resistivity_error = resistivity * np.sqrt(term_error_squared)
    
    # 5) Calculate temperature using the formula from the image
    # T = (103 + 38.1*rho - 0.095*rho^2 + 2.48e-4*rho^3)
    # Where rho is in 10^-8 ohm*m, so we need to convert our resistivity
    rho_scaled = resistivity * 1e8  # Convert to units of 10^-8 ohm*m
    rho_scaled_error = resistivity_error * 1e8
    
    temperature_in_K = 103 + 38.1 * rho_scaled - 0.095 * rho_scaled**2 + 2.48e-4 * rho_scaled**3
    
    # Error propagation for temperature
    # For T = 103 + 38.1*rho - 0.095*rho^2 + 2.48e-4*rho^3
    # dT/drho = 38.1 - 2*0.095*rho + 3*2.48e-4*rho^2
    dT_drho = 38.1 - 2 * 0.095 * rho_scaled + 3 * 2.48e-4 * rho_scaled**2
    temperature_error_in_K = np.abs(dT_drho * rho_scaled_error)
    
    # 6) Extract the max intensity angle measurements
    max_intensity_angle_1 = max_intensity_df["max_intensity_angle_measurement_1-rad"].values
    max_intensity_angle_2 = max_intensity_df["max_intensity_angle_measurement_2-rad"].values
    max_intensity_angle_3 = max_intensity_df["max_intensity_angle_measurement_3-rad"].values
    
    # Calculate the mean angle for each measurement
    max_intensity_angles = np.column_stack([max_intensity_angle_1, max_intensity_angle_2, max_intensity_angle_3])
    mean_max_intensity_angle = np.mean(max_intensity_angles, axis=1)
    
    # Calculate the standard error of the mean for each measurement
    max_intensity_angle_error = np.std(max_intensity_angles, axis=1, ddof=1) / np.sqrt(3)
    
    # 7) Print results
    print("\nMaximum Intensity Wavelength Analysis:")
    print("Measurement | Voltage (V) | Current (A) | Resistance (Ω) | Resistivity (10^-8 Ω·m) | Temperature (K) | Max Intensity Angle (rad)")
    print("-" * 110)
    
    for i in range(len(voltage_in_V)):
        print(f"{i+1:^11} | {voltage_in_V[i]:^11.3f} | {current_in_A[i]:^11.3f} | {resistance_in_ohms[i]:^14.3f} | {rho_scaled[i]:^23.3f} | {temperature_in_K[i]:^15.1f} | {mean_max_intensity_angle[i]:^24.3f}")
    
    print("\nDetailed Results with Errors:")
    for i in range(len(voltage_in_V)):
        print(f"\nMeasurement {i+1}:")
        print(f"    Voltage = {voltage_in_V[i]:.3f} ± {voltage_error_in_V[i]:.3f} V")
        print(f"    Current = {current_in_A[i]:.3f} ± {current_error_in_A[i]:.3f} A")
        print(f"    Resistance = {resistance_in_ohms[i]:.3f} ± {resistance_error_in_ohms[i]:.3f} Ω")
        print(f"    Resistivity = {rho_scaled[i]:.3f} ± {rho_scaled_error[i]:.3f} × 10^-8 Ω·m")
        print(f"    Temperature = {temperature_in_K[i]:.1f} ± {temperature_error_in_K[i]:.1f} K")
        print(f"    Max Intensity Angle = {mean_max_intensity_angle[i]:.3f} ± {max_intensity_angle_error[i]:.3f} rad")
    
    # 8) Create a plot of temperature vs. max intensity angle
    plt.figure(figsize=(10, 6))
    plt.errorbar(temperature_in_K, mean_max_intensity_angle, 
                xerr=temperature_error_in_K, yerr=max_intensity_angle_error,
                fmt='o', capsize=5, markersize=8, elinewidth=1, label='Measurements')
    
    plt.xlabel('Temperature (K)')
    plt.ylabel('Maximum Intensity Angle (rad)')
    plt.title('Maximum Intensity Angle vs. Temperature')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    plt.savefig("max_intensity_vs_temperature.png", dpi=300, bbox_inches="tight")
    plt.close()




if __name__ == "__main__":
    main()
