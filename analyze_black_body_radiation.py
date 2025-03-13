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
    starting_angle_measurements_in_rad: NDArray[np.float64] = starting_angle_df[col_name].values

    # 3) Compute mean + error (with std dev and SEM)
    STARTING_ANGLE_IN_RAD: float = np.mean(starting_angle_measurements_in_rad)
    STARTING_ANGLE_ERROR_IN_RAD: float = np.std(starting_angle_measurements_in_rad, ddof=1) / np.sqrt(
        len(starting_angle_measurements_in_rad)
    )

    print("\nStarting Angle Analysis:")
    print(f"    Starting angle = {STARTING_ANGLE_IN_RAD:.3f} ± {STARTING_ANGLE_ERROR_IN_RAD:.3f} rad")

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
    rho0_in_ohm_m = 5.65e-8  # Resistivity at room temperature in ohm*meters
    R0 = R_LIGHTBULB_IN_OHMS  # Resistance of lightbulb at room temperature
    Rw = R_WIRES_IN_OHMS  # Resistance of the wires
    
    resistivity_in_ohm_m: NDArray[np.float64] = rho0_in_ohm_m * ((resistance_in_ohms - Rw) / R0)
    
    # Error propagation for resistivity
    # For rho = rho0 * ((R - Rw) / R0)
    # We need to propagate errors from R, Rw, and R0
    # Calculate the error term for the resistivity:
    # Δρ/ρ = Δ[(R-Rw)/R₀]/[(R-Rw)/R₀] = √( [Δ(R-Rw)/(R-Rw)]² + [ΔR₀/R₀]² ) = √( [ΔR/(R-Rw)]² + [ΔRw/(R-Rw)]² + [ΔR₀/R₀]² )
    term_error_squared: NDArray[np.float64] = (resistance_error_in_ohms / (resistance_in_ohms - Rw))**2 + (R_WIRES_ERROR_IN_OHMS / (resistance_in_ohms - Rw))**2 + (R_LIGHTBULB_ERROR_IN_OHMS / R0)**2
    resistivity_error_in_ohm_m: NDArray[np.float64] = resistivity_in_ohm_m * np.sqrt(term_error_squared)
    
    # 5) Calculate temperature using the formula from the lab manual
    # T = (103 + 38.1*rho - 0.095*rho^2 + 2.48e-4*rho^3)
    # Where rho is in 10^-8 ohm*m, so we need to convert our resistivity
    rho_scaled: NDArray[np.float64] = resistivity_in_ohm_m * 1e8  # Convert to units of 10^-8 ohm*m
    rho_scaled_error: NDArray[np.float64] = resistivity_error_in_ohm_m * 1e8
    
    TEMPERATURE_IN_K: NDArray[np.float64] = 103 + 38.1 * rho_scaled - 0.095 * rho_scaled**2 + 2.48e-4 * rho_scaled**3
    
    # Error propagation for temperature
    # For T = 103 + 38.1*rho - 0.095*rho^2 + 2.48e-4*rho^3
    # dT/drho = 38.1 - 2*0.095*rho + 3*2.48e-4*rho^2
    dT_drho: NDArray[np.float64] = 38.1 - 2 * 0.095 * rho_scaled + 3 * 2.48e-4 * rho_scaled**2
    TEMPERATURE_ERROR_IN_K: NDArray[np.float64] = np.abs(dT_drho * rho_scaled_error)

    print("\nTemperature Analysis:")
    for i in range(len(TEMPERATURE_IN_K)):
        print(f"    Temperature = {TEMPERATURE_IN_K[i]:.6f} ± {TEMPERATURE_ERROR_IN_K[i]:.6f} K")
    
    # 6) Extract the max intensity angle measurements
    max_intensity_angle_1: NDArray[np.float64] = max_intensity_df["max_intensity_angle_measurement_1-rad"].values
    max_intensity_angle_2: NDArray[np.float64] = max_intensity_df["max_intensity_angle_measurement_2-rad"].values
    max_intensity_angle_3: NDArray[np.float64] = max_intensity_df["max_intensity_angle_measurement_3-rad"].values
    
    # Calculate the mean angle for each measurement
    max_intensity_angles: NDArray[np.float64] = np.column_stack([max_intensity_angle_1, max_intensity_angle_2, max_intensity_angle_3])
    mean_max_intensity_angle_in_rad: NDArray[np.float64] = np.mean(max_intensity_angles, axis=1)
    
    # Calculate the standard error of the mean for each measurement
    max_intensity_angle_error_in_rad: NDArray[np.float64] = np.std(max_intensity_angles, axis=1, ddof=1) / np.sqrt(3)

    # 7) Calculate the refraction index using the formula from the lab manual
    # n = sqrt((1.1547*sin((Init - Angle)/Ratio) + 0.5)^2 + 0.75)
    # Ignoring the filter part as requested
    
    # Calculate the term inside the sin function: (Init - Angle)/Ratio
    # Note: We're using the unadjusted angles directly as they are already in the correct form
    inner_term: NDArray[np.float64] = (STARTING_ANGLE_IN_RAD - mean_max_intensity_angle_in_rad) / ANGLE_RATIO
    
    # Calculate the sin term: 1.1547*sin(inner_term) + 0.5
    sin_term: NDArray[np.float64] = 1.1547 * np.sin(inner_term) + 0.5
    
    # Calculate the refraction index: sqrt(sin_term^2 + 0.75)
    refraction_index: NDArray[np.float64] = np.sqrt(sin_term**2 + 0.75)
    
    # Error propagation
    # For n = sqrt((1.1547*sin((Init - Angle)/Ratio) + 0.5)^2 + 0.75)
    
    # For inner_term = (STARTING_ANGLE_IN_RAD - mean_max_intensity_angle_in_rad) / ANGLE_RATIO
    # The error comes from STARTING_ANGLE_IN_RAD, mean_max_intensity_angle_in_rad, and ANGLE_RATIO
    inner_term_error: NDArray[np.float64] = np.sqrt(
        (STARTING_ANGLE_ERROR_IN_RAD / ANGLE_RATIO)**2 + 
        (max_intensity_angle_error_in_rad / ANGLE_RATIO)**2 + 
        ((STARTING_ANGLE_IN_RAD - mean_max_intensity_angle_in_rad) * ANGLE_RATIO_ERROR / ANGLE_RATIO**2)**2
    )
    
    # For sin_term = 1.1547 * np.sin(inner_term) + 0.5
    # The derivative of sin(x) is cos(x)
    sin_term_error: NDArray[np.float64] = 1.1547 * np.abs(np.cos(inner_term) * inner_term_error)
    
    # For refraction_index = np.sqrt(sin_term**2 + 0.75)
    # The derivative of sqrt(x) is 1/(2*sqrt(x))
    # The derivative of sin_term**2 + 0.75 with respect to sin_term is 2*sin_term
    # Therefore Δn = 1/(2*sqrt(sin_term**2 + 0.75)) * (2*sin_term) * Δsin_term = sin_term * Δsin_term / refraction_index
    refraction_index_error: NDArray[np.float64] = np.abs(sin_term * sin_term_error / refraction_index)

    print("\nRefraction Index Analysis:")
    for i in range(len(refraction_index)):
        print(f"    Refraction index = {refraction_index[i]:.6f} ± {refraction_index_error[i]:.6f}")

    # 8) Calculate wavelengths using the formula from equation (2.12)
    # λ = 3000/√(A + Bn + Cn² + Dn³ + En⁴ + Fn⁵ + Gn⁶ + Hn⁷ + In⁸)
    
    # Define the coefficients from the table
    A = -49852133
    B = 86092018.9
    C = -29983328.35
    D = -14354236.56
    E = 835425.05
    F = 5647432.02
    G = 1863438.86
    H = -2719226.18
    I = 574967.82
    
    # Calculate the denominator polynomial for each refraction index
    n = refraction_index
    denominator: NDArray[np.float64] = np.sqrt(
        A + B*n + C*n**2 + D*n**3 + E*n**4 + F*n**5 + G*n**6 + H*n**7 + I*n**8
    )
    
    # Calculate wavelength in nm
    WAVELENGTH_IN_NM: NDArray[np.float64] = 3000 / denominator
    
    # Error propagation for wavelength
    # For λ = 3000/√(polynomial(n))
    # Δλ = dλ/dn * Δn
    # We need dλ/dn
    
    # First, calculate d(denominator)/dn
    # d/dn[√(A + Bn + Cn² + ...)] = (1/2)(A + Bn + Cn² + ...)^(-1/2) * (B + 2Cn + 3Dn² + ...) = 
    # (B + 2Cn + 3Dn² + ...)/(2*sqrt(A + Bn + Cn² + ...)) = (B + 2Cn + 3Dn² + ...)/(2*denominator)
    polynomial_derivative: NDArray[np.float64] = (
        B + 2*C*n + 3*D*n**2 + 4*E*n**3 + 5*F*n**4 + 6*G*n**5 + 7*H*n**6 + 8*I*n**7
    )
    
    denominator_derivative: NDArray[np.float64] = 0.5 * polynomial_derivative / denominator
    
    # Now calculate dλ/dn = -3000 * (denominator)^(-2) * d(denominator)/dn
    wavelength_derivative: NDArray[np.float64] = -3000 * denominator_derivative / denominator**2
    
    # Finally, calculate the wavelength error
    WAVELENGTH_ERROR_IN_NM: NDArray[np.float64] = np.abs(wavelength_derivative * refraction_index_error)

    print("\nWavelength Analysis:")
    for i in range(len(WAVELENGTH_IN_NM)):
        print(f"    Wavelength = {WAVELENGTH_IN_NM[i]:.6f} ± {WAVELENGTH_ERROR_IN_NM[i]:.6f} nm")

    ############################################################
    # INVERSE TEMPERATURE AND WIEN'S LAW ANALYSIS
    ############################################################
    # Calculate 1/T and its error
    INVERSE_TEMPERATURE_IN_K_INVERSE: NDArray[np.float64] = 1.0 / TEMPERATURE_IN_K
    
    # Error propagation for 1/T
    # For y = 1/x, the error is Δy = |dy/dx| * Δx = |-1/x²| * Δx = Δx/x²
    INVERSE_TEMPERATURE_ERROR_IN_K_INVERSE: NDArray[np.float64] = TEMPERATURE_ERROR_IN_K / (TEMPERATURE_IN_K**2)
    
    print("\nInverse Temperature Analysis:")
    for i in range(len(INVERSE_TEMPERATURE_IN_K_INVERSE)):
        print(f"    1/T = {INVERSE_TEMPERATURE_IN_K_INVERSE[i]:.8f} ± {INVERSE_TEMPERATURE_ERROR_IN_K_INVERSE[i]:.8f} K⁻¹")
    
    # Perform ODR fit for wavelength vs 1/T (Wien's displacement law)
    # Wien's law states that λ_max * T = constant, so λ_max = constant * (1/T)
    # This means we expect a linear relationship between wavelength and 1/T with zero intercept
    
    # Convert wavelength from nm to m for physical constants
    wavelength_in_m = WAVELENGTH_IN_NM * 1e-9
    wavelength_error_in_m = WAVELENGTH_ERROR_IN_NM * 1e-9
    
    # Perform the ODR fit
    wien_constant, wien_constant_error, intercept, intercept_error = perform_odr_fit(
        INVERSE_TEMPERATURE_IN_K_INVERSE,
        wavelength_in_m,
        INVERSE_TEMPERATURE_ERROR_IN_K_INVERSE,
        wavelength_error_in_m
    )
    
    # Convert the Wien constant to mm·K for comparison with the theoretical value
    wien_constant_in_mm_K = wien_constant * 1000
    wien_constant_error_in_mm_K = wien_constant_error * 1000
    
    # Theoretical Wien's law constant (in mm·K)
    theoretical_wien_constant = 2.9 # mm·K
    
    # Calculate the relative difference
    relative_difference = abs(wien_constant_in_mm_K - theoretical_wien_constant) / theoretical_wien_constant * 100
    
    print("\nWien's Law Analysis:")
    print(f"    Experimental Wien's constant = {wien_constant_in_mm_K:.6f} ± {wien_constant_error_in_mm_K:.6f} mm·K")
    print(f"    Theoretical Wien's constant = {theoretical_wien_constant:.6f} mm·K")
    print(f"    Relative difference = {relative_difference:.2f}%")
    
    # Create a plot of wavelength vs 1/T
    fig = plot_odr_fit(
        INVERSE_TEMPERATURE_IN_K_INVERSE,
        wavelength_in_m,
        wien_constant,
        intercept,
        wien_constant_error,
        intercept_error,
        x_err=INVERSE_TEMPERATURE_ERROR_IN_K_INVERSE,
        y_err=wavelength_error_in_m,
        x_label="1/T (K⁻¹)",
        y_label="Wavelength (m)",
        title="Wien's Displacement Law: λ_max vs 1/T",
    )
    
    plt.savefig("wien_law_fit.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
