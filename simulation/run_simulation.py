#!/usr/bin/env python3
"""
Main script to run the diffraction grating simulation with various options.
"""

import sys
import os
import argparse
import numpy as np

# Add parent directory to path for imports to work
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from simulation.ui.interactive import run_interactive_simulation
from simulation.core.simulator import DiffractionSimulator


def parse_wavelength(wavelength_str):
    """Parse wavelength string to a value in meters."""
    try:
        # Try parsing as a direct value (in meters)
        return float(wavelength_str)
    except ValueError:
        # Try parsing as a value with units
        if wavelength_str.endswith('nm'):
            return float(wavelength_str.rstrip('nm')) * 1e-9
        elif wavelength_str.endswith('um') or wavelength_str.endswith('µm'):
            return float(wavelength_str.rstrip('umµ')) * 1e-6
        else:
            raise ValueError(f"Unrecognized wavelength format: {wavelength_str}. Use e.g. '500nm' or '0.5e-6'.")


def gaussian_spectrum(center, width):
    """
    Create a Gaussian spectrum intensity function.
    
    Args:
        center: Center wavelength in meters
        width: Width (standard deviation) in meters
        
    Returns:
        Function that takes wavelengths and returns intensities
    """
    def intensity_function(wavelengths):
        return np.exp(-((wavelengths - center) / width) ** 2 / 2)
    
    return intensity_function


def double_peak_spectrum(center1, center2, width):
    """
    Create a double-peaked spectrum intensity function.
    
    Args:
        center1: First peak center wavelength in meters
        center2: Second peak center wavelength in meters
        width: Width (standard deviation) in meters
        
    Returns:
        Function that takes wavelengths and returns intensities
    """
    def intensity_function(wavelengths):
        return (np.exp(-((wavelengths - center1) / width) ** 2 / 2) + 
                np.exp(-((wavelengths - center2) / width) ** 2 / 2))
    
    return intensity_function


def main():
    """Main function to run the simulation based on command line arguments."""
    parser = argparse.ArgumentParser(description='Diffraction Grating Simulation')
    
    # Define simulation modes
    parser.add_argument('--mode', type=str, default='interactive',
                       choices=['interactive', 'interactive-infinite', 'spectrum', 'spectrum-infinite', 
                               'monochromatic', 'monochromatic-infinite', 'custom', 'custom-infinite',
                               'compare'],
                       help='Simulation mode (default: interactive)')
    
    # Define general parameters
    parser.add_argument('--spacing', type=float, default=1.67e-6,
                       help='Grating spacing in meters (default: 1.67e-6 for 600 lines/mm)')
    parser.add_argument('--distance', type=float, default=1.0,
                       help='Distance to screen in meters (default: 1.0)')
    parser.add_argument('--width', type=float, default=1.0,
                       help='Screen width in meters (default: 1.0)')
    parser.add_argument('--slits', type=int, default=100,
                       help='Number of slits (default: 100)')
    
    # Define wavelength parameters (for monochromatic modes)
    parser.add_argument('--wavelength', type=str, default='532nm',
                       help='Wavelength (for monochromatic modes) in nm or m (default: 532nm)')
    
    # Define spectrum parameters (for custom modes)
    spectrum_group = parser.add_argument_group('Custom spectrum options')
    spectrum_group.add_argument('--spectrum-type', type=str, default='gaussian',
                              choices=['gaussian', 'double-peak'],
                              help='Type of spectrum to generate (default: gaussian)')
    spectrum_group.add_argument('--center', type=str, default='550nm',
                              help='Center wavelength for Gaussian spectrum (default: 550nm)')
    spectrum_group.add_argument('--spectrum-width', type=str, default='50nm',
                              help='Width (standard deviation) for spectrum (default: 50nm)')
    spectrum_group.add_argument('--center2', type=str, default='650nm',
                              help='Second center wavelength for double-peak spectrum (default: 650nm)')
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Convert spacing from lines/mm if needed
    if args.spacing > 0 and args.spacing < 0.1:
        # If the value is small, assume it's directly in meters
        grating_spacing = args.spacing
    else:
        # Otherwise, assume it's lines/mm and convert
        grating_spacing = 1 / (args.spacing * 1000)  # Convert to meters
    
    # Interactive UI modes
    if args.mode == 'interactive' or args.mode == 'interactive-infinite':
        use_infinite_slits = (args.mode == 'interactive-infinite')
        run_interactive_simulation(use_infinite_slits=use_infinite_slits)
        return
    
    # Create simulation instance for all other modes
    simulator = DiffractionSimulator(
        grating_spacing=grating_spacing,
        distance_to_screen=args.distance,
        screen_width=args.width,
        num_slits=args.slits
    )
    
    # Run the appropriate simulation mode
    if args.mode == 'spectrum' or args.mode == 'spectrum-infinite':
        use_infinite_slits = (args.mode == 'spectrum-infinite')
        simulator.simulate_visible_spectrum(
            show_maxima=True,
            use_infinite_slits=use_infinite_slits
        )
    
    elif args.mode == 'monochromatic' or args.mode == 'monochromatic-infinite':
        wavelength = parse_wavelength(args.wavelength)
        use_infinite_slits = (args.mode == 'monochromatic-infinite')
        simulator.simulate_monochromatic(
            wavelength=wavelength,
            show_maxima=True,
            use_infinite_slits=use_infinite_slits
        )
    
    elif args.mode == 'custom' or args.mode == 'custom-infinite':
        # Parse spectrum parameters
        center = parse_wavelength(args.center)
        width = parse_wavelength(args.spectrum_width)
        use_infinite_slits = (args.mode == 'custom-infinite')
        
        if args.spectrum_type == 'gaussian':
            intensity_function = gaussian_spectrum(center, width)
        elif args.spectrum_type == 'double-peak':
            center2 = parse_wavelength(args.center2)
            intensity_function = double_peak_spectrum(center, center2, width)
        
        simulator.simulate_custom_spectrum(
            intensity_function=intensity_function,
            wavelength_range=(400e-9, 700e-9),  # Visible range
            num_wavelengths=20,
            use_infinite_slits=use_infinite_slits
        )
    
    elif args.mode == 'compare':
        wavelength = parse_wavelength(args.wavelength)
        slit_numbers = [2, 5, 10, 50, 100, 500]
        
        simulator.compare_slit_numbers(
            wavelength=wavelength,
            slit_numbers=slit_numbers,
            show_infinite=True
        )


if __name__ == "__main__":
    main() 