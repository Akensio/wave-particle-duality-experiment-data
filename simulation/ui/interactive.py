#!/usr/bin/env python3
"""
Interactive UI components for the diffraction grating simulation.
Provides an interactive interface for experimenting with diffraction parameters.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, CheckButtons, RadioButtons
from typing import List, Dict, Callable, Optional, Union, Tuple, Set

from simulation.core.physics import DiffractionPhysics, InfiniteSlitDiffractionPhysics
from simulation.core.config import (
    DEFAULT_GRATING_SPACING, DEFAULT_DISTANCE_TO_SCREEN, DEFAULT_NUM_SLITS,
    DEFAULT_SCREEN_WIDTH, WAVELENGTH_OPTIONS, COLOR_MAP, GRATING_PRESETS,
    DEFAULT_WAVELENGTHS
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
        self.selected_wavelengths = DEFAULT_WAVELENGTHS.copy()
        
        # Create physics model based on mode
        self._init_physics_model()
        
        # Create visualizer
        self.visualizer = DiffractionVisualizer()
        
        # Create the figure and axes
        self._create_figure()
        
        # Initial plot
        self._update_plots()
        
        # Create sliders and controls
        self._create_controls()
    
    def _init_physics_model(self) -> None:
        """Initialize physics model based on simulation mode."""
        if self.use_infinite_slits:
            self.physics = InfiniteSlitDiffractionPhysics(
                grating_spacing=self.grating_spacing,
                distance_to_screen=self.distance_to_screen,
                screen_width=self.screen_width
            )
        else:
            self.physics = DiffractionPhysics(
                grating_spacing=self.grating_spacing,
                distance_to_screen=self.distance_to_screen,
                screen_width=self.screen_width
            )
    
    def _create_figure(self) -> None:
        """Create the figure and axes for the interactive simulation."""
        self.fig, (self.ax1, self.ax_total, self.ax2) = plt.subplots(
            3, 1, figsize=(12, 14), gridspec_kw={'height_ratios': [3, 1, 1]}
        )
        
        suptitle = f"Interactive Diffraction Grating Simulation"
        if self.use_infinite_slits:
            suptitle += " (Infinite Slits)"
        suptitle += f"\nGrating spacing: {self.grating_spacing*1e6:.1f} μm"
        
        if not self.use_infinite_slits:
            suptitle += f", Number of slits: {self.num_slits}"
            
        self.fig.suptitle(suptitle, fontsize=16)
        
        plt.subplots_adjust(left=0.15, bottom=0.25, right=0.85, top=0.95, hspace=0.4)
    
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
        active = [i for i, color in enumerate(check_labels) if color in self.selected_wavelengths]
        
        self.check = CheckButtons(
            ax_check, check_labels, 
            [i in active for i in range(len(check_labels))]
        )
        
        # Create radio buttons for grating presets
        ax_radio = plt.axes([0.02, 0.15, 0.1, 0.1])
        self.radio = RadioButtons(
            ax_radio, list(GRATING_PRESETS.keys()), active=0
        )
        
        # Connect callback functions for all controls
        self.grating_slider.on_changed(self._update)
        self.distance_slider.on_changed(self._update)
        self.width_slider.on_changed(self._update)
        self.check.on_clicked(self._update_wavelengths)
        self.radio.on_clicked(self._update_grating_preset)
    
    def _update_plots(self) -> None:
        """Update all plots based on current parameters."""
        # Clear axes
        self.ax1.clear()
        self.ax_total.clear()
        self.ax2.clear()
        
        # Get selected wavelengths
        wavelengths = [WAVELENGTH_OPTIONS[color] for color in self.selected_wavelengths]
        colors = [COLOR_MAP[color] for color in self.selected_wavelengths]
        labels = [f"{color} ({wl*1e9:.0f} nm)" for color, wl in zip(self.selected_wavelengths, wavelengths)]
        
        if not wavelengths:
            # If no wavelengths selected, show empty plot
            self.ax1.set_title("Select at least one wavelength")
            return
        
        # Plot individual wavelengths
        for i, (wavelength, color, label) in enumerate(zip(wavelengths, colors, labels)):
            if self.use_infinite_slits:
                screen_positions, intensity = self.physics.calculate_intensity_pattern(wavelength)
            else:
                screen_positions, intensity = self.physics.calculate_intensity_pattern(
                    wavelength, num_slits=self.num_slits
                )
            
            self.ax1.plot(screen_positions, intensity, color=color, label=label)
            
            # Add diffraction maxima markers
            maxima = self.physics.calculate_maxima_positions(wavelength)
            for m, pos in maxima:
                self.ax1.axvline(x=pos, color=color, linestyle='--', alpha=0.3)
        
        # Plot total intensity
        if len(wavelengths) > 1:
            if self.use_infinite_slits:
                self._plot_combined_intensity(self.ax_total, wavelengths)
            else:
                # For finite slits, calculate combined intensity
                self._plot_combined_intensity(self.ax_total, wavelengths, self.num_slits)
        
        # Set up axis labels and title
        self.ax1.set_xlim(-self.screen_width/2, self.screen_width/2)
        self.ax1.set_ylim(0, 1.05)
        self.ax1.set_xlabel("Position on Screen (m)")
        self.ax1.set_ylabel("Relative Intensity")
        self.ax1.grid(True, alpha=0.3)
        self.ax1.legend()
        
        # Set up title
        title = f"Diffraction Pattern"
        if self.use_infinite_slits:
            title += " (Infinite Slits)"
        else:
            title += f" ({self.num_slits} Slits)"
        title += f"\nGrating spacing: {self.grating_spacing*1e6:.1f} µm, Distance to screen: {self.distance_to_screen:.1f} m"
        self.ax1.set_title(title)
        
        # Set up total intensity axis
        self.ax_total.set_xlim(-self.screen_width/2, self.screen_width/2)
        self.ax_total.set_ylim(0, 1.05)
        self.ax_total.set_xlabel("Position on Screen (m)")
        self.ax_total.set_ylabel("Total Intensity")
        self.ax_total.grid(True, alpha=0.3)
        self.ax_total.set_title("Combined Intensity (Sum of All Selected Wavelengths)")
        
        # Create a wall pattern visualization
        self._create_wall_pattern(wavelengths, colors)
        
        self.fig.canvas.draw_idle()
    
    def _plot_combined_intensity(self, ax: plt.Axes, wavelengths: List[float], 
                               num_slits: Optional[int] = None) -> None:
        """
        Plot combined intensity of multiple wavelengths.
        
        Args:
            ax: Axes to plot on
            wavelengths: List of wavelengths
            num_slits: Number of slits (for finite slits model)
        """
        # Calculate total intensity
        if self.use_infinite_slits:
            positions, total = self._calculate_combined_intensity(wavelengths)
        else:
            positions, total = self._calculate_combined_intensity(wavelengths, num_slits)
            
        # Plot total intensity
        ax.plot(positions, total, 'k-', linewidth=2)
    
    def _calculate_combined_intensity(self, wavelengths: List[float], 
                                    num_slits: Optional[int] = None) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate total intensity from multiple wavelengths.
        
        Args:
            wavelengths: List of wavelengths
            num_slits: Number of slits (for finite slits model)
            
        Returns:
            Tuple of (positions, total intensity)
        """
        # Get positions from first wavelength
        if self.use_infinite_slits or num_slits is None:
            positions, _ = self.physics.calculate_intensity_pattern(wavelengths[0])
        else:
            positions, _ = self.physics.calculate_intensity_pattern(wavelengths[0], num_slits=num_slits)
        
        # Initialize total intensity
        total = np.zeros_like(positions)
        
        # Add contribution from each wavelength
        for wl in wavelengths:
            if self.use_infinite_slits or num_slits is None:
                _, intensity = self.physics.calculate_intensity_pattern(wl)
            else:
                _, intensity = self.physics.calculate_intensity_pattern(wl, num_slits=num_slits)
            
            total += intensity
        
        # Normalize
        if np.max(total) > 0:
            total = total / np.max(total)
            
        return positions, total
    
    def _create_wall_pattern(self, wavelengths: List[float], colors: List[str]) -> None:
        """
        Create the 2D wall pattern visualization.
        
        Args:
            wavelengths: List of wavelengths
            colors: List of colors
        """
        # Create visualization in the bottom axes
        if self.use_infinite_slits:
            self.visualizer._add_wall_pattern_visualization(
                self.fig, self.physics, wavelengths, colors, ax=self.ax2
            )
        else:
            self.visualizer._add_wall_pattern_visualization(
                self.fig, self.physics, wavelengths, colors, num_slits=self.num_slits, ax=self.ax2
            )
            
        self.ax2.set_title("Pattern on Wall")
    
    def _update(self, val: float) -> None:
        """
        Update callback for sliders.
        
        Args:
            val: New value from slider (not used directly)
        """
        # Get values from sliders
        self.grating_spacing = self.grating_slider.val * 1e-6  # Convert from µm to m
        self.distance_to_screen = self.distance_slider.val
        self.screen_width = self.width_slider.val
        
        if not self.use_infinite_slits:
            self.num_slits = int(self.slits_slider.val)
        
        # Update physics model parameters
        self.physics.grating_spacing = self.grating_spacing
        self.physics.distance_to_screen = self.distance_to_screen
        self.physics.screen_width = self.screen_width
        
        # Update the figure title
        suptitle = f"Interactive Diffraction Grating Simulation"
        if self.use_infinite_slits:
            suptitle += " (Infinite Slits)"
        suptitle += f"\nGrating spacing: {self.grating_spacing*1e6:.1f} μm"
        
        if not self.use_infinite_slits:
            suptitle += f", Number of slits: {self.num_slits}"
            
        self.fig.suptitle(suptitle, fontsize=16)
        
        # Update all plots
        self._update_plots()
    
    def _update_wavelengths(self, label: str) -> None:
        """
        Update callback for wavelength checkboxes.
        
        Args:
            label: Wavelength label that was toggled
        """
        # Toggle the selected wavelength
        if label in self.selected_wavelengths:
            self.selected_wavelengths.remove(label)
        else:
            self.selected_wavelengths.append(label)
        
        # Update plots
        self._update_plots()
    
    def _update_grating_preset(self, label: str) -> None:
        """
        Update callback for grating presets.
        
        Args:
            label: Selected preset label
        """
        # Get preset values
        preset = GRATING_PRESETS[label]
        self.grating_spacing = preset['spacing']
        
        if not self.use_infinite_slits:
            self.num_slits = preset['slits']
        
        # Update slider values
        self.grating_slider.set_val(self.grating_spacing * 1e6)  # Convert to µm
        
        if not self.use_infinite_slits:
            self.slits_slider.set_val(self.num_slits)
        
        # Update physics model
        self.physics.grating_spacing = self.grating_spacing
        
        # Update plots
        self._update_plots()


def run_interactive_simulation(use_infinite_slits: bool = False) -> None:
    """
    Run the interactive simulation.
    
    Args:
        use_infinite_slits: Whether to use the infinite slits model
    """
    # Create and show the simulation
    simulation = InteractiveSimulation(use_infinite_slits=use_infinite_slits)
    plt.show() 