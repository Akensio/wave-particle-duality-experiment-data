#!/usr/bin/env python3
"""
Black Body Radiation Spectrum Analysis

This script analyzes experimental data from a black body radiation spectrum experiment
and fits it to Planck's law to determine the temperature of the black body.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import constants
from scipy.optimize import curve_fit
import os

datafile = os.path.join(
    os.path.dirname(__file__), "data", "black_body_radiation_spectrum.xlsx"
)


def main():
    """Main function to run the analysis"""
    ############################################################
    # ROTATION RATIO
    ############################################################
    # 1) Load the sheet
    rotation_df = pd.read_excel(datafile, sheet_name="rotation_ratio")

    # 2) Extract the single column
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
    starting_angle_measurements_in_rad = starting_angle_df[col_name].values
    
    # 3) Compute mean + error (with std dev and SEM)
    starting_angle_in_rad = np.mean(starting_angle_measurements_in_rad)
    starting_angle_error_in_rad = np.std(starting_angle_measurements_in_rad, ddof=1) / np.sqrt(len(starting_angle_measurements_in_rad))
    
    print("\nStarting Angle Analysis:")
    print(f"    Mean starting angle = {starting_angle_in_rad:.3f} rad")
    print(f"    Standard error = {starting_angle_error_in_rad:.3f} rad")
    print(f"    Starting angle = {starting_angle_in_rad:.3f} ± {starting_angle_error_in_rad:.3f} rad")

    ############################################################
    # RESISTANCE MEASUREMENTS
    ############################################################
    # 1) Load the sheet
    resistance_df = pd.read_excel(datafile, sheet_name="lightbulb_four_wires")
    
    # 2) Extract the columns
    voltage_in_mV = resistance_df["voltage-mV"].values
    current_in_mA = resistance_df["current-mA"].values
    



if __name__ == "__main__":
    main()
