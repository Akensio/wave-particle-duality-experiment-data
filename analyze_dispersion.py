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
    reciprocal_sqrt_voltage = voltage_in_V.apply(lambda v: 1 / np.sqrt(v))
    reciprocal_sqrt_voltage_error = 0.5 * reciprocal_sqrt_voltage * (0.1 / voltage_in_V) # voltage_in_V.apply(lambda v: (0.5 / np.sqrt(v)) * (0.1 / v))

    # 4) Calculate the de Broglie wavelength
    h_in_J_s = 6.62607015e-34
    m_in_kg = 9.1093837015e-31
    e_in_C = 1.602176634e-19
    wavelength_in_m = (h_in_J_s / np.sqrt(2 * m_in_kg * e_in_C)) * reciprocal_sqrt_voltage
    wavelength_error_in_m = wavelength_in_m * reciprocal_sqrt_voltage_error / reciprocal_sqrt_voltage

    wavelength_in_nm = wavelength_in_m * 1e9
    wavelength_error_in_nm = wavelength_error_in_m * 1e9

    print("Wavelength Analysis:")
    for i in range(len(voltage_in_V)):
        print(f"    Voltage = {voltage_in_V[i]:.3f} V, Wavelength = {wavelength_in_nm[i]:.3f} ± {wavelength_error_in_nm[i]:.3f} nm")

    #########################################################
    # Dispersion Analysis
    #########################################################

    # 1) Extract the columns for r1 - the measurements of the internal ring diameter
    # These are the two measurements of the inner part of the ring
    diameter_in_internal_mm = dispersion_df["diameter_in_internal-mm"].values
    diameter2_in_internal_mm = dispersion_df["diameter2_in_internal-mm"].values
    # These are the two measurements of the outer part of the ring
    diameter_in_external_mm = dispersion_df["diameter_in_external-mm"].values
    diameter2_in_external_mm = dispersion_df["diameter2_in_external-mm"].values

    # 2) Calculate r1 and its error
    r1s = (diameter_in_internal_mm + diameter2_in_internal_mm) / 4
    dr1s = (diameter_in_internal_mm - diameter2_in_internal_mm) / 8



if __name__ == "__main__":
    main() 