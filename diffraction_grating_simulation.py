#!/usr/bin/env python3
# Diffraction Grating Simulation
# This program simulates light passing through a diffraction grating and displays the resulting pattern

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, CheckButtons, RadioButtons
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
                    'c': [0, 1, 1], 'm': [1, 0, 1], 'y': [1, 1, 0], 'k': [0, 0, 0],
                    'darkorange': [1, 0.55, 0], 'gold': [1, 0.84, 0], 'violet': [0.93, 0.51, 0.93]}
        
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


def interactive_diffraction_simulation():
    """
    Create an interactive simulation with sliders to control parameters.
    """
    # Set up the initial parameters
    initial_grating_spacing = 2e-6  # 2 micrometers
    initial_distance = 1.0  # 1 meter
    initial_slits = 5
    initial_screen_width = 1.0  # 1 meter (wider to see more peaks)
    
    # Predefined wavelengths for visible spectrum
    wavelength_options = {
        'Red': 700e-9,
        'Orange': 620e-9,
        'Yellow': 580e-9,
        'Green': 530e-9,
        'Blue': 470e-9,
        'Violet': 420e-9
    }
    
    # Color mapping for each wavelength
    color_map = {
        'Red': 'r',
        'Orange': 'darkorange',
        'Yellow': 'gold',
        'Green': 'g',
        'Blue': 'b',
        'Violet': 'violet'
    }
    
    # Set up the figure and axes
    fig, (ax1, ax_total, ax2) = plt.subplots(3, 1, figsize=(12, 14), 
                                           gridspec_kw={'height_ratios': [3, 1, 1]})
    plt.subplots_adjust(left=0.15, bottom=0.25, right=0.85, top=0.95, hspace=0.4)
    
    # Make space for the sliders
    fig.subplots_adjust(bottom=0.25)
    
    # Initial calculation and plot
    selected_wavelengths = ['Red', 'Green', 'Blue']  # Default selection
    wavelength_values = [wavelength_options[color] for color in selected_wavelengths]
    colors = [color_map[color] for color in selected_wavelengths]
    
    # Calculate initial patterns
    lines = []
    maxima_lines = []
    
    for i, wavelength in enumerate(wavelength_values):
        screen_positions, intensity = calculate_intensity_pattern(
            wavelength, initial_grating_spacing, initial_distance, initial_screen_width, num_slits=initial_slits
        )
        line, = ax1.plot(screen_positions, intensity, color=colors[i], label=f"{selected_wavelengths[i]} ({wavelength*1e9:.0f} nm)", linewidth=2)
        lines.append(line)
    
    # Set up the axes
    ax1.set_xlim(-initial_screen_width/2, initial_screen_width/2)
    ax1.set_ylim(0, 1.05)
    ax1.set_xlabel("Position on Screen (m)")
    ax1.set_ylabel("Relative Intensity")
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='upper right')
    
    # Calculate total intensity (sum of all wavelengths)
    screen_positions, _ = calculate_intensity_pattern(
        wavelength_values[0], initial_grating_spacing, initial_distance, 
        initial_screen_width, num_slits=initial_slits
    )
    
    total_intensity = np.zeros_like(screen_positions)
    for wavelength in wavelength_values:
        _, intensity = calculate_intensity_pattern(
            wavelength, initial_grating_spacing, initial_distance, 
            initial_screen_width, num_slits=initial_slits
        )
        total_intensity += intensity
    
    # Normalize the total intensity
    if np.max(total_intensity) > 0:
        total_intensity = total_intensity / np.max(total_intensity)
    
    # Plot total intensity
    ax_total.plot(screen_positions, total_intensity, 'k-', linewidth=2)
    ax_total.set_xlim(-initial_screen_width/2, initial_screen_width/2)
    ax_total.set_ylim(0, 1.05)
    ax_total.set_ylabel("Total Intensity")
    ax_total.grid(True, alpha=0.3)
    ax_total.set_title("Combined Intensity (Sum of All Selected Wavelengths)")
    
    # Create the wall pattern visualization
    wall_height = initial_screen_width/2  # Make it rectangle for better visibility
    wall_positions = np.linspace(-initial_screen_width/2, initial_screen_width/2, 1000)
    wall_y = np.linspace(-wall_height, wall_height, 500)
    
    # Initialize the wall pattern
    wall_pattern = np.zeros((len(wall_y), len(wall_positions), 3))
    
    for i, wavelength in enumerate(wavelength_values):
        _, intensity = calculate_intensity_pattern(
            wavelength, initial_grating_spacing, initial_distance, initial_screen_width, num_slits=initial_slits
        )
        
        # Interpolate intensity to match wall_positions
        from scipy.interpolate import interp1d
        screen_positions, _ = calculate_intensity_pattern(
            wavelength, initial_grating_spacing, initial_distance, initial_screen_width, num_slits=initial_slits
        )
        intensity_interp = interp1d(screen_positions, intensity, kind='linear', bounds_error=False, fill_value=0)
        intensity_wall = intensity_interp(wall_positions)
        
        # Create 2D pattern (same in vertical direction)
        pattern_2d = np.tile(intensity_wall, (len(wall_y), 1))
        
        # Get RGB values
        color_rgb_map = {'r': [1, 0, 0], 'g': [0, 1, 0], 'b': [0, 0, 1], 
                        'darkorange': [1, 0.55, 0], 'gold': [1, 0.84, 0], 'violet': [0.93, 0.51, 0.93]}
        
        rgb = color_rgb_map[colors[i]]
        
        # Add to RGB channels
        for j in range(3):
            wall_pattern[:, :, j] += pattern_2d * rgb[j] / len(wavelength_values)
    
    # Clip values
    wall_pattern = np.clip(wall_pattern, 0, 1)
    
    # Display the wall pattern
    wall_img = ax2.imshow(wall_pattern, extent=[-initial_screen_width/2, initial_screen_width/2, -wall_height, wall_height],
                         aspect='auto', interpolation='bilinear')
    ax2.set_title("Pattern on Wall")
    ax2.set_xlabel("Horizontal Position (m)")
    ax2.set_ylabel("Vertical Position (m)")
    
    # Title for the entire figure
    fig.suptitle(f"Interactive Diffraction Grating Simulation\nGrating spacing: {initial_grating_spacing*1e6:.1f} μm, Number of slits: {initial_slits}", fontsize=16)
    
    # Create sliders
    ax_grating = plt.axes([0.15, 0.17, 0.7, 0.02])
    ax_distance = plt.axes([0.15, 0.13, 0.7, 0.02])
    ax_slits = plt.axes([0.15, 0.09, 0.7, 0.02])
    ax_width = plt.axes([0.15, 0.05, 0.7, 0.02])
    
    # Define slider ranges
    grating_slider = Slider(ax_grating, 'Grating Spacing (μm)', 0.5, 10.0, valinit=initial_grating_spacing*1e6)
    distance_slider = Slider(ax_distance, 'Distance to Screen (m)', 0.1, 5.0, valinit=initial_distance)
    slits_slider = Slider(ax_slits, 'Number of Slits', 2, 20, valinit=initial_slits, valstep=1)
    width_slider = Slider(ax_width, 'Screen Width (m)', 0.2, 3.0, valinit=initial_screen_width)
    
    # Create checkboxes for wavelength selection
    ax_check = plt.axes([0.02, 0.45, 0.1, 0.2])
    check_labels = list(wavelength_options.keys())
    
    # Initialize with RGB colors checked
    initial_selection = [check_labels.index(color) for color in selected_wavelengths]
    check = CheckButtons(ax_check, check_labels, [l in selected_wavelengths for l in check_labels])
    
    # Update function for the sliders and checkboxes
    def update(val):
        # Get current values from sliders
        grating_spacing = grating_slider.val * 1e-6  # Convert μm to m
        distance = distance_slider.val
        num_slits = int(slits_slider.val)
        screen_width = width_slider.val
        
        # Get selected wavelengths
        selected_colors = []
        for i, label in enumerate(check_labels):
            if check.get_status()[i]:
                selected_colors.append(label)
        
        # Make sure at least one color is selected
        if not selected_colors:
            selected_colors = ['Green']  # Default to green if nothing selected
        
        # Get wavelength values and colors
        wavelength_values = [wavelength_options[color] for color in selected_colors]
        colors = [color_map[color] for color in selected_colors]
        
        # Clear existing plots
        ax1.clear()
        
        # Recalculate and plot for each wavelength
        for i, wavelength in enumerate(wavelength_values):
            screen_positions, intensity = calculate_intensity_pattern(
                wavelength, grating_spacing, distance, screen_width, num_slits=num_slits
            )
            ax1.plot(screen_positions, intensity, color=colors[i], 
                    label=f"{selected_colors[i]} ({wavelength*1e9:.0f} nm)", linewidth=2)
            
            # Add diffraction maxima markers
            max_order = 5  # Show more orders
            for m in range(-max_order, max_order + 1):
                if abs(m * wavelength / grating_spacing) < 1:  # Only if physically possible
                    angle = np.arcsin(m * wavelength / grating_spacing)
                    pos = distance * np.tan(angle)
                    if abs(pos) <= screen_width/2:  # Only if on screen
                        ax1.axvline(x=pos, color=colors[i], linestyle='--', alpha=0.3)
                        # Only label maxima for the first wavelength to avoid clutter
                        if i == 0:
                            ax1.text(pos, 1.02, f"m={m}", ha='center', color='black', fontsize=8)
        
        # Update axis limits
        ax1.set_xlim(-screen_width/2, screen_width/2)
        ax1.set_ylim(0, 1.05)
        ax1.set_xlabel("Position on Screen (m)")
        ax1.set_ylabel("Relative Intensity")
        ax1.grid(True, alpha=0.3)
        ax1.legend(loc='upper right')
        
        # Calculate and plot total intensity
        if wavelength_values:
            # Get screen positions from the first wavelength
            screen_positions, _ = calculate_intensity_pattern(
                wavelength_values[0], grating_spacing, distance, screen_width, num_slits=num_slits
            )
            
            # Calculate total intensity (sum of all wavelengths)
            total_intensity = np.zeros_like(screen_positions)
            for wavelength in wavelength_values:
                _, intensity = calculate_intensity_pattern(
                    wavelength, grating_spacing, distance, screen_width, num_slits=num_slits
                )
                total_intensity += intensity
            
            # Normalize the total intensity
            if np.max(total_intensity) > 0:
                total_intensity = total_intensity / np.max(total_intensity)
        else:
            # If no wavelengths are selected, show flat line
            screen_positions = np.linspace(-screen_width/2, screen_width/2, 1000)
            total_intensity = np.zeros_like(screen_positions)
        
        # Update total intensity plot
        ax_total.clear()
        ax_total.plot(screen_positions, total_intensity, 'k-', linewidth=2)
        ax_total.set_xlim(-screen_width/2, screen_width/2)
        ax_total.set_ylim(0, 1.05)
        ax_total.set_ylabel("Total Intensity")
        ax_total.grid(True, alpha=0.3)
        ax_total.set_title("Combined Intensity (Sum of All Selected Wavelengths)")
        
        # Update the wall pattern
        wall_height = screen_width/2
        wall_positions = np.linspace(-screen_width/2, screen_width/2, 1000)
        wall_y = np.linspace(-wall_height, wall_height, 500)
        
        wall_pattern = np.zeros((len(wall_y), len(wall_positions), 3))
        
        for i, wavelength in enumerate(wavelength_values):
            _, intensity = calculate_intensity_pattern(
                wavelength, grating_spacing, distance, screen_width, num_slits=num_slits
            )
            
            # Interpolate intensity to match wall_positions
            from scipy.interpolate import interp1d
            screen_positions, _ = calculate_intensity_pattern(
                wavelength, grating_spacing, distance, screen_width, num_slits=num_slits
            )
            intensity_interp = interp1d(screen_positions, intensity, kind='linear', bounds_error=False, fill_value=0)
            intensity_wall = intensity_interp(wall_positions)
            
            # Create 2D pattern
            pattern_2d = np.tile(intensity_wall, (len(wall_y), 1))
            
            # Get RGB values
            color_rgb_map = {'r': [1, 0, 0], 'g': [0, 1, 0], 'b': [0, 0, 1], 
                           'darkorange': [1, 0.55, 0], 'gold': [1, 0.84, 0], 'violet': [0.93, 0.51, 0.93]}
            
            rgb = color_rgb_map[colors[i]]
            
            # Add to RGB channels
            for j in range(3):
                wall_pattern[:, :, j] += pattern_2d * rgb[j] / len(wavelength_values)
        
        # Clip values
        wall_pattern = np.clip(wall_pattern, 0, 1)
        
        # Update the wall display
        ax2.clear()
        ax2.imshow(wall_pattern, extent=[-screen_width/2, screen_width/2, -wall_height, wall_height],
                  aspect='auto', interpolation='bilinear')
        ax2.set_title("Pattern on Wall")
        ax2.set_xlabel("Horizontal Position (m)")
        ax2.set_ylabel("Vertical Position (m)")
        
        # Update the main title
        fig.suptitle(f"Interactive Diffraction Grating Simulation\n"
                   f"Grating spacing: {grating_spacing*1e6:.1f} μm, Number of slits: {num_slits}", fontsize=16)
        
        fig.canvas.draw_idle()
    
    # Connect the update function to the sliders and checkbox
    grating_slider.on_changed(update)
    distance_slider.on_changed(update)
    slits_slider.on_changed(update)
    width_slider.on_changed(update)
    check.on_clicked(lambda label: update(0))  # Pass a dummy value
    
    plt.show()


def main():
    # Run the interactive version
    print("Starting interactive diffraction grating simulation...")
    interactive_diffraction_simulation()


if __name__ == "__main__":
    main() 