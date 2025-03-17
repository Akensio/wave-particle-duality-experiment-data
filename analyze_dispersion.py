#!/usr/bin/env python3
"""
Dispersion Analysis

This script analyzes experimental data from a voltage-diameter relationship experiment.
"""

import os
import numpy as np
import pandas as pd
from utils import odr_fit
from numpy.typing import NDArray
import matplotlib.pyplot as plt

datafile = os.path.join(os.path.dirname(__file__), "data", "dispersion.xlsx")


def main():
    # Load data
    df = pd.read_excel("data/dispersion.xlsx")

    #########################################################
    # Wavelength calculations
    #########################################################
    # 1) Load the sheet
    dispersion_df: pd.DataFrame = pd.read_excel(datafile, sheet_name="dispersion")

    # 2) Extract the columns
    voltage_in_V: NDArray[np.float64] = dispersion_df["voltage-kv"].values * 1000

    # 3) Calculate f(V)=1/√V and its error
    # Δf(V) = |df/dV| × ΔV = |-1/2 × V^(-3/2)| × ΔV = (1/2) × V^(-3/2) × ΔV = (0.5/√V) × (ΔV/V)
    # ΔV = 0.1 kV (Measurement device error)
    reciprocal_sqrt_voltage_in_v = 1 / np.sqrt(voltage_in_V)
    reciprocal_sqrt_voltage_error_in_v = (
        0.5 * reciprocal_sqrt_voltage_in_v * (100 / voltage_in_V)
    )

    # 4) Calculate the de Broglie wavelength
    h_in_J_s = 6.62607015e-34
    m_in_kg = 9.1093837015e-31
    e_in_C = 1.602176634e-19
    wavelength_in_m = (h_in_J_s / np.sqrt(2 * m_in_kg * e_in_C)) * reciprocal_sqrt_voltage_in_v
    wavelength_error_in_m = wavelength_in_m * reciprocal_sqrt_voltage_error_in_v / reciprocal_sqrt_voltage_in_v

    wavelength_in_pm = wavelength_in_m * 1e12
    wavelength_error_in_pm = wavelength_error_in_m * 1e12

    print("\nWavelength Analysis:")
    for i in range(len(voltage_in_V)):
        print(f"    Voltage = {voltage_in_V[i]:.3f} V, Wavelength = {wavelength_in_pm[i]:.3f} ± {wavelength_error_in_pm[i]:.3f} pm")

    #########################################################
    # Radius calculations
    #########################################################

    # 1) Extract the measurements of ring1 diameter
    # These are the two measurements of the inner part of ring1
    diameter_in_internal_mm = dispersion_df["diameter_in_internal-mm"].values
    diameter_in_internal_meas2_mm = dispersion_df["diameter2_in_internal-mm"].values
    # These are the two measurements of the outer part of ring1
    diameter_in_external_mm = dispersion_df["diameter_in_external-mm"].values
    diameter_in_external_meas2_mm = dispersion_df["diameter2_in_external-mm"].values

    # 2) Calculate the mean and error of the first measurement of ring1
    diameter1_meas1_mm = (diameter_in_internal_mm + diameter_in_external_mm) / 2
    diameter1_meas1_error_mm = np.abs(diameter_in_external_mm - diameter_in_internal_mm) / 2

    # 3) Calculate the mean and error of the second measurement of ring1
    diameter1_meas2_mm = (diameter_in_internal_meas2_mm + diameter_in_external_meas2_mm) / 2
    diameter1_meas2_error_mm = np.abs(diameter_in_external_meas2_mm - diameter_in_internal_meas2_mm) / 2

    # 4) Calculate the mean and error of d1
    diameter1_mm = (diameter1_meas1_mm + diameter1_meas2_mm) / 2
    diameter1_error_mm = np.sqrt(diameter1_meas1_error_mm**2 + diameter1_meas2_error_mm**2) / 2

    r1_mm = diameter1_mm / 2
    r1_error_mm = diameter1_error_mm / 2

    # Now we do the same calculations for ring2.

    # 5) Extract the measurements of ring2 diameter
    # These are the two measurements of the inner part of ring2
    diameter_ex_internal_mm = dispersion_df["diameter_ex_internal-mm"].values
    diameter_ex_internal_meas2_mm = dispersion_df["diameter2_ex_internal-mm"].values
    # These are the two measurements of the outer part of ring2
    diameter_ex_external_mm = dispersion_df["diameter_ex_external-mm"].values
    diameter_ex_external_meas2_mm = dispersion_df["diameter2_ex_external-mm"].values

    # 6) Calculate the mean and error of the first measurement of ring2
    diameter2_meas1_mm = (diameter_ex_internal_mm + diameter_ex_external_mm) / 2
    diameter2_meas1_error_mm = np.abs(diameter_ex_external_mm - diameter_ex_internal_mm) / 2

    # 7) Calculate the mean and error of the second measurement of ring2
    diameter2_meas2_mm = (diameter_ex_internal_meas2_mm + diameter_ex_external_meas2_mm) / 2
    diameter2_meas2_error_mm = np.abs(diameter_ex_external_meas2_mm - diameter_ex_internal_meas2_mm) / 2

    # 8) Calculate the mean and error of d2
    diameter2_mm = (diameter2_meas1_mm + diameter2_meas2_mm) / 2
    diameter2_error_mm = np.sqrt(diameter2_meas1_error_mm**2 + diameter2_meas2_error_mm**2) / 2

    r2_mm = diameter2_mm / 2
    r2_error_mm = diameter2_error_mm / 2

    print("\nRadius Analysis of ring 1:")
    for i in range(len(r1_mm)):
        print(f"    Ring 1: r = {r1_mm[i]:.3f} ± {r1_error_mm[i]:.3f} mm")

    print("\nRadius Analysis of ring 2:")
    for i in range(len(r2_mm)):
        print(f"    Ring 2: r = {r2_mm[i]:.3f} ± {r2_error_mm[i]:.3f} mm")

    #########################################################
    # Dispersion analysis
    #########################################################

    # 1) ODR fit of r1 vs wavelength and r2 vs wavelength
    slope1, slope1_error, intercept1, intercet1_error = odr_fit.perform_odr_fit(
        wavelength_in_pm, r1_mm, wavelength_error_in_pm, r1_error_mm
    )
    slope2, slope2_error, intercept2, intercet2_error = odr_fit.perform_odr_fit(
        wavelength_in_pm, r2_mm, wavelength_error_in_pm, r2_error_mm
    )

    # 2) Build graphs of r1 vs wavelength and r2 vs wavelength
    fig = odr_fit.plot_odr_fit(
        wavelength_in_pm,
        r1_mm,
        slope1,
        intercept1,
        slope1_error,
        intercet1_error,
        wavelength_error_in_pm,
        r1_error_mm,
        "Wavelength (pm)",
        "Radius (mm)",
        "r1 vs wavelength",
        "r1_vs_wavelength.png",
    )
    fig.savefig("r1_vs_wavelength.png")
    plt.close(fig)

    fig = odr_fit.plot_odr_fit(
        wavelength_in_pm,
        r2_mm,
        slope2,
        intercept2,
        slope2_error,
        intercet2_error,
        wavelength_error_in_pm,
        r2_error_mm,
        "Wavelength (pm)",
        "Radius (mm)",
        "r2 vs wavelength",
        "r2_vs_wavelength.png",
    )
    fig.savefig("r2_vs_wavelength.png")
    plt.close(fig)

    print("\nSlopes:")
    print(f"    r1_slope: {slope1:.3f} ± {slope1_error:.3f}")
    print(f"    r2_slope: {slope2:.3f} ± {slope2_error:.3f}")

    # 3) Calculate the first two interplanar distances in the crystal
    DIAMETER_OF_BULB_IN_MM = 135
    DIAMETER_OF_BULB_ERROR_IN_MM = 3
    
    RADIUS_OF_BULB_IN_MM = DIAMETER_OF_BULB_IN_MM / 2
    RADIUS_OF_BULB_ERROR_IN_MM = DIAMETER_OF_BULB_ERROR_IN_MM / 2

    


if __name__ == "__main__":
    main()
