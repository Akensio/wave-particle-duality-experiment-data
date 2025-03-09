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

# Set up matplotlib for high-quality plots
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 12

def planck_law(wavelength, T, A):
    """
    Planck's law for spectral radiance.
    
    Parameters:
    wavelength : array-like, wavelength in meters
    T : float, temperature in Kelvin
    A : float, amplitude scaling factor
    
    Returns:
    array-like, spectral radiance
    """
    h = constants.h  # Planck's constant
    c = constants.c  # Speed of light
    k = constants.k  # Boltzmann constant
    
    # Convert wavelength from nm to m if needed
    if np.mean(wavelength) > 1e-6:  # If wavelength is in nm
        wavelength = wavelength * 1e-9  # Convert to meters
    
    # Calculate the exponent term
    exponent = h * c / (wavelength * k * T)
    
    # Handle potential numerical issues with large exponents
    with np.errstate(over='ignore'):
        spectral_radiance = A * (2.0 * h * c**2) / (wavelength**5) / (np.exp(exponent) - 1.0)
    
    # Replace any NaN or inf values with 0
    spectral_radiance = np.nan_to_num(spectral_radiance)
    
    return spectral_radiance

def main():
    """Main function to run the analysis"""
    print("Black Body Radiation Spectrum Analysis")
    print("-" * 40)
    
    # Load the data from Excel file
    file_path = 'data/black_body_radiation_spectrum.xlsx'
    print(f"Loading data from {file_path}...")
    
    try:
        # First, check the Excel file structure
        xls = pd.ExcelFile(file_path)
        print(f"Excel file contains {len(xls.sheet_names)} sheets: {xls.sheet_names}")
        
        # Read data from each sheet
        for sheet_name in xls.sheet_names:
            print(f"\nAnalyzing sheet: {sheet_name}")
            df = pd.read_excel(file_path, sheet_name=sheet_name)
            
            print(f"Sheet dimensions: {df.shape[0]} rows x {df.shape[1]} columns")
            print("Column names:", df.columns.tolist())
            print("\nFirst few rows:")
            print(df.head())
            
            # Check if this sheet has wavelength and intensity data
            wavelength_cols = [col for col in df.columns if 'wavelength' in col.lower() or 'lambda' in col.lower() or 'nm' in col.lower()]
            intensity_cols = [col for col in df.columns if 'intensity' in col.lower() or 'radiance' in col.lower() or 'power' in col.lower()]
            
            if wavelength_cols and intensity_cols:
                print(f"\nFound potential wavelength column(s): {wavelength_cols}")
                print(f"Found potential intensity column(s): {intensity_cols}")
                
                # Use the first identified columns
                wavelength_col = wavelength_cols[0]
                intensity_col = intensity_cols[0]
                
                # Extract data
                wavelengths = df[wavelength_col].values
                intensities = df[intensity_col].values
                
                # Remove any NaN values
                valid_indices = ~(np.isnan(wavelengths) | np.isnan(intensities))
                wavelengths = wavelengths[valid_indices]
                intensities = intensities[valid_indices]
                
                if len(wavelengths) < 10:
                    print("Not enough valid data points for analysis.")
                    continue
                
                print(f"\nAnalyzing {len(wavelengths)} data points...")
                
                # Plot the raw data
                plt.figure(figsize=(12, 6))
                plt.scatter(wavelengths, intensities, color='blue', label='Experimental Data')
                plt.xlabel(wavelength_col)
                plt.ylabel(intensity_col)
                plt.title(f'Black Body Radiation Spectrum - {sheet_name}')
                plt.grid(True, alpha=0.3)
                plt.legend()
                plt.savefig(f'black_body_raw_data_{sheet_name}.png', dpi=300, bbox_inches='tight')
                plt.close()
                
                # Fit to Planck's law
                try:
                    # Initial guess for parameters
                    T_initial = 3000  # K (typical for a light bulb)
                    A_initial = 1e12  # Arbitrary scaling factor
                    initial_params = [T_initial, A_initial]
                    
                    # Perform the curve fitting
                    params, covariance = curve_fit(
                        planck_law, wavelengths, intensities, 
                        p0=initial_params, 
                        bounds=([500, 1e8], [10000, 1e15]),
                        maxfev=10000
                    )
                    
                    # Extract the fitted parameters
                    T_fit, A_fit = params
                    T_err, A_err = np.sqrt(np.diag(covariance))
                    
                    print(f"Fitted temperature: {T_fit:.2f} ± {T_err:.2f} K")
                    print(f"Scaling factor: {A_fit:.3e} ± {A_err:.3e}")
                    
                    # Calculate the fitted curve
                    wavelength_fine = np.linspace(min(wavelengths), max(wavelengths), 1000)
                    intensity_fit = planck_law(wavelength_fine, T_fit, A_fit)
                    
                    # Plot the results
                    plt.figure(figsize=(12, 8))
                    plt.scatter(wavelengths, intensities, color='blue', label='Experimental Data')
                    plt.plot(wavelength_fine, intensity_fit, 'r-', label=f'Planck\'s Law Fit (T = {T_fit:.2f} K)')
                    plt.xlabel(wavelength_col)
                    plt.ylabel(intensity_col)
                    plt.title(f'Black Body Radiation Spectrum with Planck\'s Law Fit - {sheet_name}')
                    plt.grid(True, alpha=0.3)
                    plt.legend()
                    plt.savefig(f'black_body_fit_{sheet_name}.png', dpi=300, bbox_inches='tight')
                    plt.close()
                    
                    # Calculate residuals
                    residuals = intensities - planck_law(wavelengths, T_fit, A_fit)
                    
                    # Plot residuals
                    plt.figure(figsize=(12, 4))
                    plt.scatter(wavelengths, residuals)
                    plt.axhline(y=0, color='r', linestyle='-')
                    plt.xlabel(wavelength_col)
                    plt.ylabel('Residuals')
                    plt.title(f'Fit Residuals - {sheet_name}')
                    plt.grid(True, alpha=0.3)
                    plt.savefig(f'black_body_residuals_{sheet_name}.png', dpi=300, bbox_inches='tight')
                    plt.close()
                    
                    # Find the wavelength of maximum intensity in the fitted curve
                    max_intensity_idx = np.argmax(intensity_fit)
                    lambda_max = wavelength_fine[max_intensity_idx]
                    print(f"Wavelength of maximum emission: {lambda_max:.2f} nm")
                    
                    # Calculate temperature using Wien's displacement law
                    b = 2.898e-3  # Wien's displacement constant in m·K
                    lambda_max_m = lambda_max * 1e-9  # Convert nm to m
                    T_wien = b / lambda_max_m
                    
                    print(f"Temperature from Wien's law: {T_wien:.2f} K")
                    print(f"Difference from fitted temperature: {abs(T_wien - T_fit):.2f} K ({abs(T_wien - T_fit)/T_fit*100:.2f}%)")
                    
                    # Calculate total power using Stefan-Boltzmann law
                    sigma = constants.sigma  # Stefan-Boltzmann constant
                    total_power = sigma * (T_fit ** 4)  # W/m²
                    
                    print(f"Total power emitted per unit area: {total_power:.2e} W/m²")
                    
                except Exception as e:
                    print(f"Error during curve fitting: {e}")
            else:
                print("Could not identify wavelength and intensity columns in this sheet.")
    
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main() 