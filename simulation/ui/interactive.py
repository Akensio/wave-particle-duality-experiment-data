#!/usr/bin/env python3
# Interactive UI components for the diffraction grating simulation

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, CheckButtons, RadioButtons
from typing import List, Dict, Callable, Optional, Union, Tuple

from simulation.core.physics import DiffractionModel
from simulation.core.config import (
    DEFAULT_GRATING_SPACING, DEFAULT_DISTANCE_TO_SCREEN, DEFAULT_NUM_SLITS,
    DEFAULT_SCREEN_WIDTH, WAVELENGTH_OPTIONS, COLOR_MAP, GRATING_PRESETS,
    DEFAULT_WAVELENGTHS, COLOR_RGB_MAP
)
from simulation.visualization.plotter import DiffractionVisualizer

class InteractiveSimulation:
    """Class for creating interactive diffraction grating simulations."""
    
    def __init__(self, use_infinite_slits: bool = False):
        """
        Initialize the interactive simulation.
        
        Args:
            use_infinite_slits: Whether to use the infinite slits model
        """
        self.use_infinite_slits = use_infinite_slits
        
        # Initialize parameters
        self.grating_spacing = DEFAULT_GRATING_SPACING
        self.distance_to_screen = DEFAULT_DISTANCE_TO_SCREEN
        self.num_slits = DEFAULT_NUM_SLITS
        self.screen_width = DEFAULT_SCREEN_WIDTH
        
        # Initialize selected wavelengths
        self.selected_wavelengths = DEFAULT_WAVELENGTHS[:]
        
        # Create the figure and axes
        if self.use_infinite_slits:
            self.fig, (self.ax1, self.ax_total, self.ax2) = plt.subplots(
                3, 1, figsize=(12, 14), gridspec_kw={'height_ratios': [3, 1, 1]}
            )
            self.fig.suptitle(
                f"Interactive Diffraction Grating Simulation (Infinite Slits)\n"
                f"Grating spacing: {self.grating_spacing*1e6:.1f} μm", 
                fontsize=16
            )
        else:
            self.fig, (self.ax1, self.ax_total, self.ax2) = plt.subplots(
                3, 1, figsize=(12, 14), gridspec_kw={'height_ratios': [3, 1, 1]}
            )
            self.fig.suptitle(
                f"Interactive Diffraction Grating Simulation\n"
                f"Grating spacing: {self.grating_spacing*1e6:.1f} μm, Number of slits: {self.num_slits}", 
                fontsize=16
            )
        
        plt.subplots_adjust(left=0.15, bottom=0.25, right=0.85, top=0.95, hspace=0.4)
        
        # Initial plot
        self._update_plots()
        
        # Create sliders and controls
        self._create_controls()
        
    def _create_controls(self) -> None:
        """Create sliders, checkboxes, and radio buttons for interaction."""
        # Create sliders
        ax_grating = plt.axes([0.15, 0.17, 0.7, 0.02])
        ax_distance = plt.axes([0.15, 0.13, 0.7, 0.02])
        ax_width = plt.axes([0.15, 0.05, 0.7, 0.02])
        
        self.grating_slider = Slider(
            ax_grating, 'Grating Spacing (µm)', 
            0.5, 10.0, valinit=self.grating_spacing*1e6
        )
        self.distance_slider = Slider(
            ax_distance, 'Distance to Screen (m)', 
            0.1, 5.0, valinit=self.distance_to_screen
        )
        self.width_slider = Slider(
            ax_width, 'Screen Width (m)', 
            0.2, 3.0, valinit=self.screen_width
        )
        
        # Add slits slider only for finite slits model
        if not self.use_infinite_slits:
            ax_slits = plt.axes([0.15, 0.09, 0.7, 0.02])
            self.slits_slider = Slider(
                ax_slits, 'Number of Slits', 
                2, 1000, valinit=self.num_slits, valstep=1
            )
            self.slits_slider.on_changed(self._update)
        
        # Create checkboxes for wavelength selection
        ax_check = plt.axes([0.02, 0.45, 0.1, 0.2])
        check_labels = list(WAVELENGTH_OPTIONS.keys())
        
        # Initialize with default colors checked
        self.check = CheckButtons(
            ax_check, check_labels, 
            [label in self.selected_wavelengths for label in check_labels]
        )
        
        # Add realistic presets for common lab gratings
        ax_preset = plt.axes([0.02, 0.25, 0.1, 0.15])
        preset_buttons = RadioButtons(
            ax_preset, list(GRATING_PRESETS.keys()), active=0
        )
        
        # Connect event handlers
        self.grating_slider.on_changed(self._update)
        self.distance_slider.on_changed(self._update)
        self.width_slider.on_changed(self._update)
        self.check.on_clicked(lambda label: self._update(0))  # Dummy value
        preset_buttons.on_clicked(self._handle_preset)
    
    def _handle_preset(self, label: str) -> None:
        """
        Handle preset button clicks.
        
        Args:
            label: The preset label that was clicked
        """
        if label in GRATING_PRESETS:
            preset = GRATING_PRESETS[label]
            self.grating_slider.set_val(preset['spacing'] * 1e6)  # Convert to µm
            
            if not self.use_infinite_slits and hasattr(self, 'slits_slider'):
                self.slits_slider.set_val(preset['slits'])
        
        self._update(0)  # Update with dummy value
    
    def _update(self, val: Union[float, str]) -> None:
        """
        Update the simulation when controls are changed.
        
        Args:
            val: The value from the control that triggered the update (not used)
        """
        # Get current values from sliders
        self.grating_spacing = self.grating_slider.val * 1e-6  # Convert μm to m
        self.distance_to_screen = self.distance_slider.val
        self.screen_width = self.width_slider.val
        
        if not self.use_infinite_slits and hasattr(self, 'slits_slider'):
            self.num_slits = int(self.slits_slider.val)
        
        # Get selected wavelengths
        check_labels = list(WAVELENGTH_OPTIONS.keys())
        self.selected_wavelengths = [
            label for i, label in enumerate(check_labels) 
            if self.check.get_status()[i]
        ]
        
        # Make sure at least one color is selected
        if not self.selected_wavelengths:
            self.selected_wavelengths = ['Green']  # Default to green if nothing selected
        
        # Update the plots
        self._update_plots()
    
    def _update_plots(self) -> None:
        """Update all plots with current parameter values."""
        # Clear the axes
        self.ax1.clear()
        self.ax_total.clear()
        self.ax2.clear()
        
        # Get wavelength values and colors
        wavelength_values = [WAVELENGTH_OPTIONS[color] for color in self.selected_wavelengths]
        colors = [COLOR_MAP[color] for color in self.selected_wavelengths]
        labels = [f"{color} ({WAVELENGTH_OPTIONS[color]*1e9:.0f} nm)" for color in self.selected_wavelengths]
        
        # Plot diffraction patterns
        if self.use_infinite_slits:
            # Plot for infinite slits
            for i, wavelength in enumerate(wavelength_values):
                screen_positions, intensity = DiffractionModel.calculate_infinite_slit_pattern(
                    wavelength, self.grating_spacing, self.distance_to_screen, self.screen_width
                )
                self.ax1.plot(screen_positions, intensity, color=colors[i], label=labels[i], linewidth=2)
        else:
            # Plot for finite slits
            for i, wavelength in enumerate(wavelength_values):
                screen_positions, intensity = DiffractionModel.calculate_intensity_pattern(
                    wavelength, self.grating_spacing, self.distance_to_screen, 
                    self.screen_width, num_slits=self.num_slits
                )
                self.ax1.plot(screen_positions, intensity, color=colors[i], label=labels[i], linewidth=2)
                
                # Add diffraction maxima markers
                max_order = 5  # Show more orders
                maxima = DiffractionModel.calculate_maxima_positions(
                    wavelength, self.grating_spacing, self.distance_to_screen, 
                    self.screen_width, max_order
                )
                
                for m, pos in maxima:
                    self.ax1.axvline(x=pos, color=colors[i], linestyle='--', alpha=0.3)
                    # Only label maxima for the first wavelength to avoid clutter
                    if i == 0:
                        self.ax1.text(pos, 1.02, f"m={m}", ha='center', color='black', fontsize=8)
        
        # Set up the axes
        self.ax1.set_xlim(-self.screen_width/2, self.screen_width/2)
        self.ax1.set_ylim(0, 1.05)
        self.ax1.set_xlabel("Position on Screen (m)")
        self.ax1.set_ylabel("Relative Intensity")
        self.ax1.grid(True, alpha=0.3)
        self.ax1.legend(loc='upper right')
        
        # Plot total intensity
        if self.use_infinite_slits:
            # For infinite slits, we need to calculate differently
            screen_positions = np.linspace(-self.screen_width/2, self.screen_width/2, 1000)
            total_intensity = np.zeros_like(screen_positions)
            
            for wavelength in wavelength_values:
                # Get screen positions and intensity for this wavelength
                screen_pos, intensity = DiffractionModel.calculate_infinite_slit_pattern(
                    wavelength, self.grating_spacing, self.distance_to_screen, self.screen_width
                )
                
                # Interpolate to match the desired positions
                from scipy.interpolate import interp1d
                intensity_interp = interp1d(screen_pos, intensity, kind='linear', 
                                         bounds_error=False, fill_value=0)
                total_intensity += intensity_interp(screen_positions)
            
            # Normalize
            if np.max(total_intensity) > 0:
                total_intensity = total_intensity / np.max(total_intensity)
        else:
            # For finite slits, use the built-in function
            screen_positions, total_intensity = DiffractionModel.calculate_total_intensity(
                wavelength_values, self.grating_spacing, self.distance_to_screen, 
                self.screen_width, self.num_slits
            )
        
        # Plot total intensity
        self.ax_total.plot(screen_positions, total_intensity, 'k-', linewidth=2)
        self.ax_total.set_xlim(-self.screen_width/2, self.screen_width/2)
        self.ax_total.set_ylim(0, 1.05)
        self.ax_total.set_ylabel("Total Intensity")
        self.ax_total.grid(True, alpha=0.3)
        self.ax_total.set_title("Combined Intensity (Sum of All Selected Wavelengths)")
        
        # Update wall pattern visualization
        self._update_wall_pattern(wavelength_values, colors)
        
        # Update title
        if self.use_infinite_slits:
            self.fig.suptitle(
                f"Interactive Diffraction Grating Simulation (Infinite Slits)\n"
                f"Grating spacing: {self.grating_spacing*1e6:.1f} μm", 
                fontsize=16
            )
        else:
            self.fig.suptitle(
                f"Interactive Diffraction Grating Simulation\n"
                f"Grating spacing: {self.grating_spacing*1e6:.1f} μm, Number of slits: {self.num_slits}", 
                fontsize=16
            )
        
        self.fig.canvas.draw_idle()
    
    def _update_wall_pattern(self, wavelength_values: List[float], colors: List[str]) -> None:
        """
        Update the wall pattern visualization.
        
        Args:
            wavelength_values: List of wavelengths in meters
            colors: List of colors for each wavelength
        """
        wall_height = self.screen_width/2
        wall_positions = np.linspace(-self.screen_width/2, self.screen_width/2, 1000)
        wall_y = np.linspace(-wall_height, wall_height, 500)
        
        wall_pattern = np.zeros((len(wall_y), len(wall_positions), 3))
        
        for i, wavelength in enumerate(wavelength_values):
            if self.use_infinite_slits:
                screen_positions, intensity = DiffractionModel.calculate_infinite_slit_pattern(
                    wavelength, self.grating_spacing, self.distance_to_screen, self.screen_width
                )
            else:
                screen_positions, intensity = DiffractionModel.calculate_intensity_pattern(
                    wavelength, self.grating_spacing, self.distance_to_screen, 
                    self.screen_width, num_slits=self.num_slits
                )
            
            # Interpolate intensity to match wall_positions
            from scipy.interpolate import interp1d
            intensity_interp = interp1d(screen_positions, intensity, kind='linear', 
                                      bounds_error=False, fill_value=0)
            intensity_wall = intensity_interp(wall_positions)
            
            # Create 2D pattern
            pattern_2d = np.tile(intensity_wall, (len(wall_y), 1))
            
            # Get RGB values
            if colors[i] in COLOR_RGB_MAP:
                rgb = COLOR_RGB_MAP[colors[i]]
            else:
                rgb = [1, 1, 1]  # Default to white
            
            # Add to RGB channels
            for j in range(3):
                wall_pattern[:, :, j] += pattern_2d * rgb[j] / len(wavelength_values)
        
        # Clip values
        wall_pattern = np.clip(wall_pattern, 0, 1)
        
        # Update the wall display
        self.ax2.imshow(wall_pattern, extent=[-self.screen_width/2, self.screen_width/2, 
                                           -wall_height, wall_height],
                     aspect='auto', interpolation='bilinear')
        self.ax2.set_title("Pattern on Wall")
        self.ax2.set_xlabel("Horizontal Position (m)")
        self.ax2.set_ylabel("Vertical Position (m)")
    
    def show(self) -> None:
        """Display the interactive simulation."""
        plt.show()


def run_interactive_simulation(use_infinite_slits: bool = False) -> None:
    """
    Run the interactive diffraction grating simulation.
    
    Args:
        use_infinite_slits: Whether to use the infinite slits model
    """
    print(f"Starting interactive diffraction grating simulation...")
    if use_infinite_slits:
        print("Using infinite slits model (perfect grating)")
    
    sim = InteractiveSimulation(use_infinite_slits=use_infinite_slits)
    sim.show() 