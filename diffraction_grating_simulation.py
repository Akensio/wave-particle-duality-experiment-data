#!/usr/bin/env python3
# Diffraction Grating Simulation
# This program simulates light passing through a diffraction grating and displays the resulting pattern

import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple, List


def calculate_intensity_pattern(
    wavelength: float,
    grating_spacing: float, 
    distance_to_screen: float, 
    screen_width: float, 
    num_points: int = 1000,
    num_slits: int = 5
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate the intensity pattern on a screen from a diffraction grating.
    
    Args:
        wavelength: Light wavelength in meters
        grating_spacing: Distance between slits in meters
        distance_to_screen: Distance from grating to screen in meters
        screen_width: Width of the screen in meters
        num_points: Number of points to calculate on the screen
        num_slits: Number of slits in the grating
        
    Returns:
        Tuple of (screen positions, intensity values)
    """
    # Calculate positions on the screen
    screen_positions = np.linspace(-screen_width/2, screen_width/2, num_points)
    
    # Calculate angles to each position on the screen
    angles = np.arctan(screen_positions / distance_to_screen)
    
    # Calculate the phase difference between adjacent slits
    delta = (2 * np.pi / wavelength) * grating_spacing * np.sin(angles)
    
    # Calculate the multiple slit interference factor
    intensity_factor = np.sin(num_slits * delta / 2)**2 / np.sin(delta / 2)**2
    # Replace NaN values that occur when sin(delta/2) is zero
    intensity_factor = np.nan_to_num(intensity_factor, nan=num_slits**2)
    
    # Normalize the intensity
    intensity = intensity_factor / np.max(intensity_factor)
    
    return screen_positions, intensity


def plot_diffraction_pattern(
    wavelengths: List[float],
    colors: List[str] = None,
    labels: List[str] = None,
    grating_spacing: float = 2e-6,  # 2 micrometers
    distance_to_screen: float = 1.0,  # 1 meter
    screen_width: float = 0.5,  # 0.5 meters
    num_slits: int = 5,
    show_maxima: bool = True,
    max_order: int = 3
) -> None:
    """
    Plot diffraction patterns for multiple wavelengths on the same graph.
    
    Args:
        wavelengths: List of light wavelengths in meters
        colors: List of colors for each wavelength (optional)
        labels: List of labels for each wavelength (optional)
        grating_spacing: Distance between slits in meters
        distance_to_screen: Distance from grating to screen in meters
        screen_width: Width of the screen in meters
        num_slits: Number of slits in the grating
        show_maxima: Whether to show diffraction maxima positions
        max_order: Maximum diffraction order to mark
    """
    if colors is None:
        # Default colors if not provided
        colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k'] * (len(wavelengths) // 7 + 1)
    
    if labels is None:
        # Default labels if not provided
        labels = [f"{w*1e9:.1f} nm" for w in wavelengths]
    
    # Set up the figure
    plt.figure(figsize=(12, 8))
    
    # Plot each wavelength
    for i, wavelength in enumerate(wavelengths):
        screen_positions, intensity = calculate_intensity_pattern(
            wavelength, grating_spacing, distance_to_screen, screen_width, num_slits=num_slits
        )
        
        plt.plot(screen_positions, intensity, color=colors[i], label=labels[i], linewidth=2)
        
        # Add markers for diffraction maxima angles if requested
        if show_maxima:
            for m in range(-max_order, max_order + 1):
                # Calculate position of mth maximum using grating equation
                if abs(m * wavelength / grating_spacing) < 1:  # Only if physically possible
                    angle = np.arcsin(m * wavelength / grating_spacing)
                    pos = distance_to_screen * np.tan(angle)
                    if abs(pos) <= screen_width/2:  # Only if on screen
                        plt.axvline(x=pos, color=colors[i], linestyle='--', alpha=0.3)
                        # Only label maxima for the first wavelength to avoid clutter
                        if i == 0:
                            plt.text(pos, 0.1, f"m={m}", ha='center', color='black')
    
    # Add title and labels
    plt.title(f"Diffraction Pattern from {num_slits} Slits\n"
             f"Grating spacing: {grating_spacing*1e6:.1f} μm, Distance to screen: {distance_to_screen:.1f} m")
    plt.xlabel("Position on Screen (m)")
    plt.ylabel("Relative Intensity")
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Add a 2D representation of the wall with diffraction pattern
    ax_inset = plt.axes([0.35, 0.02, 0.3, 0.25])  # [left, bottom, width, height]
    
    # Create a 2D image representing the pattern on the wall
    wall_height = screen_width  # Make it square for simplicity
    wall_y = np.linspace(-wall_height/2, wall_height/2, 500)
    wall_x, wall_y = np.meshgrid(screen_positions, wall_y)
    
    # For simplicity, assuming the pattern is the same in vertical direction
    # For each wavelength, add to the wall pattern with appropriate color and alpha
    wall_pattern = np.zeros((len(wall_y), len(screen_positions), 3))  # RGB array
    
    for i, wavelength in enumerate(wavelengths):
        _, intensity = calculate_intensity_pattern(
            wavelength, grating_spacing, distance_to_screen, screen_width, num_slits=num_slits
        )
        
        # Convert color string to RGB
        color_map = {'r': [1, 0, 0], 'g': [0, 1, 0], 'b': [0, 0, 1], 
                    'c': [0, 1, 1], 'm': [1, 0, 1], 'y': [1, 1, 0], 'k': [0, 0, 0]}
        
        # Get RGB values for the current color
        if colors[i] in color_map:
            rgb = color_map[colors[i]]
        else:
            # Default to white if color not recognized
            rgb = [1, 1, 1]
        
        # Add this wavelength's contribution to the wall pattern
        # Create a 2D pattern (same in vertical direction)
        pattern_2d = np.tile(intensity, (len(wall_y), 1))
        
        # Add to the RGB channels with appropriate color
        for j in range(3):  # RGB channels
            wall_pattern[:, :, j] += pattern_2d * rgb[j] / len(wavelengths)
    
    # Clip values to valid range [0, 1]
    wall_pattern = np.clip(wall_pattern, 0, 1)
    
    # Display the 2D pattern
    ax_inset.imshow(wall_pattern, extent=[-screen_width/2, screen_width/2, -wall_height/2, wall_height/2], 
                   aspect='auto', interpolation='bilinear')
    ax_inset.set_title("Pattern on Wall")
    ax_inset.set_xlabel("Horizontal Position (m)")
    ax_inset.set_ylabel("Vertical Position (m)")
    
    # Use subplots_adjust instead of tight_layout to avoid warning
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.3)
    plt.show()


def main():
    # Default parameters
    grating_spacing = 2e-6  # 2 micrometers between slits
    distance_to_screen = 1.0  # 1 meter
    screen_width = 0.5  # 0.5 meters
    num_slits = 5
    
    # Example 1: Compare different wavelengths - visible spectrum
    wavelengths = [
        700e-9,  # Red (700 nm)
        620e-9,  # Orange (620 nm)
        580e-9,  # Yellow (580 nm)
        530e-9,  # Green (530 nm)
        470e-9,  # Blue (470 nm)
        420e-9   # Violet (420 nm)
    ]
    
    colors = ['r', 'darkorange', 'gold', 'g', 'b', 'violet']
    labels = ["Red (700 nm)", "Orange (620 nm)", "Yellow (580 nm)", 
             "Green (530 nm)", "Blue (470 nm)", "Violet (420 nm)"]
    
    print("Plotting diffraction patterns for visible spectrum...")
    plot_diffraction_pattern(
        wavelengths=wavelengths,
        colors=colors,
        labels=labels,
        grating_spacing=grating_spacing,
        distance_to_screen=distance_to_screen,
        screen_width=screen_width,
        num_slits=num_slits
    )
    
    # Example 2: Compare different wavelengths - RGB
    rgb_wavelengths = [
        700e-9,  # Red (700 nm)
        530e-9,  # Green (530 nm)
        470e-9   # Blue (470 nm)
    ]
    
    rgb_colors = ['r', 'g', 'b']
    rgb_labels = ["Red (700 nm)", "Green (530 nm)", "Blue (470 nm)"]
    
    print("Plotting diffraction patterns for RGB wavelengths...")
    plot_diffraction_pattern(
        wavelengths=rgb_wavelengths,
        colors=rgb_colors,
        labels=rgb_labels,
        grating_spacing=grating_spacing,
        distance_to_screen=distance_to_screen,
        screen_width=screen_width,
        num_slits=num_slits
    )
    
    # Example 3: Single wavelength with more detail
    print("Plotting detailed diffraction pattern for a single wavelength...")
    plot_diffraction_pattern(
        wavelengths=[550e-9],  # Green (550 nm)
        colors=['g'],
        labels=["Green (550 nm)"],
        grating_spacing=grating_spacing,
        distance_to_screen=distance_to_screen,
        screen_width=screen_width,
        num_slits=num_slits,
        show_maxima=True,
        max_order=5
    )


if __name__ == "__main__":
    main() 