#!/usr/bin/env python3
"""Run script for the wave packet diffraction simulation."""

import argparse
import matplotlib.pyplot as plt
import numpy as np

from simulation.core.physics import GaussianWavePacketPhysics
from simulation.visualization import WavePacketVisualizer
from simulation.core.config import (
    DEFAULT_GRATING_SPACING, DEFAULT_DISTANCE_TO_SCREEN,
    DEFAULT_SCREEN_WIDTH, DEFAULT_NUM_SLITS, WAVELENGTH_OPTIONS
)


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Wave packet diffraction simulation")
    
    parser.add_argument(
        "--center-wavelength", type=str, default="green",
        choices=list(WAVELENGTH_OPTIONS.keys()),
        help="Center wavelength of the wave packet (color name)"
    )
    
    parser.add_argument(
        "--wavelength-width", type=float, default=5e-9,
        help="Width of the Gaussian wavelength distribution in nanometers"
    )
    
    parser.add_argument(
        "--grating-spacing", type=float, default=DEFAULT_GRATING_SPACING,
        help="Spacing between slits in the diffraction grating (in meters)"
    )
    
    parser.add_argument(
        "--distance", type=float, default=DEFAULT_DISTANCE_TO_SCREEN,
        help="Distance from grating to screen (in meters)"
    )
    
    parser.add_argument(
        "--screen-width", type=float, default=DEFAULT_SCREEN_WIDTH,
        help="Width of the screen (in meters)"
    )
    
    parser.add_argument(
        "--num-slits", type=int, default=DEFAULT_NUM_SLITS,
        help="Number of slits in the diffraction grating"
    )
    
    parser.add_argument(
        "--analysis", action="store_true",
        help="Show full analysis visualization with multiple plots"
    )
    
    return parser.parse_args()


def main():
    """Run the wave packet diffraction simulation."""
    args = parse_arguments()
    
    # Get the center wavelength from the color name
    center_wavelength = WAVELENGTH_OPTIONS[args.center_wavelength]
    
    # Convert the width from nanometers to meters
    wavelength_width = args.wavelength_width * 1e-9
    
    # Create the physics model
    physics = GaussianWavePacketPhysics(
        grating_spacing=args.grating_spacing,
        distance_to_screen=args.distance,
        screen_width=args.screen_width
    )
    
    # Create the visualizer
    visualizer = WavePacketVisualizer()
    
    # Generate the visualizations
    if args.analysis:
        # Full analysis with multiple plots
        fig = visualizer.plot_wave_packet_full_analysis(
            physics, center_wavelength, wavelength_width, args.num_slits
        )
        plt.suptitle(f"Gaussian Wave Packet Analysis - Center: {args.center_wavelength}, "
                    f"Width: {args.wavelength_width} nm", fontsize=16)
    else:
        # Show the diffraction pattern comparison
        fig, ax = visualizer.plot_wave_packet_diffraction(
            physics, center_wavelength, wavelength_width, args.num_slits
        )
        
        # Add a small inset with the wavelength spectrum
        ax_inset = fig.add_axes([0.15, 0.15, 0.3, 0.2])
        visualizer.plot_wave_packet_spectrum(
            physics, center_wavelength, wavelength_width, ax=ax_inset
        )
    
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main() 