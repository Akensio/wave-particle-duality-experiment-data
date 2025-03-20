#!/usr/bin/env python3
# Visualization module for the diffraction grating simulation

import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple, Dict, Optional
from scipy.interpolate import interp1d

from simulation.core.physics import DiffractionModel
from simulation.core.config import COLOR_RGB_MAP

class DiffractionVisualizer:
    """Class responsible for visualizing diffraction patterns."""
    
    @staticmethod
    def plot_diffraction_pattern(
        wavelengths: List[float],
        colors: List[str],
        labels: List[str],
        grating_spacing: float,
        distance_to_screen: float,
        screen_width: float,
        num_slits: int,
        show_maxima: bool = True,
        max_order: int = 3,
        ax: Optional[plt.Axes] = None,
        show_wall_pattern: bool = True
    ) -> Tuple[plt.Figure, plt.Axes]:
        """
        Plot diffraction patterns for multiple wavelengths.
        
        Args:
            wavelengths: List of light wavelengths in meters
            colors: List of colors for each wavelength
            labels: List of labels for each wavelength
            grating_spacing: Distance between slits in meters
            distance_to_screen: Distance from grating to screen in meters
            screen_width: Width of the screen in meters
            num_slits: Number of slits in the grating
            show_maxima: Whether to show diffraction maxima positions
            max_order: Maximum diffraction order to mark
            ax: Optional axes to plot on (will create new if None)
            show_wall_pattern: Whether to show the 2D wall pattern visualization
            
        Returns:
            Figure and axes objects
        """
        # Create new figure if needed
        if ax is None:
            fig, ax = plt.subplots(figsize=(12, 8))
        else:
            fig = ax.figure
        
        # Plot each wavelength
        for i, wavelength in enumerate(wavelengths):
            screen_positions, intensity = DiffractionModel.calculate_intensity_pattern(
                wavelength, grating_spacing, distance_to_screen, screen_width, num_slits=num_slits
            )
            
            ax.plot(screen_positions, intensity, color=colors[i], label=labels[i], linewidth=2)
            
            # Add markers for diffraction maxima if requested
            if show_maxima:
                maxima = DiffractionModel.calculate_maxima_positions(
                    wavelength, grating_spacing, distance_to_screen, screen_width, max_order
                )
                
                for m, pos in maxima:
                    ax.axvline(x=pos, color=colors[i], linestyle='--', alpha=0.3)
                    # Only label maxima for the first wavelength to avoid clutter
                    if i == 0:
                        ax.text(pos, 0.1, f"m={m}", ha='center', color='black')
        
        # Set up the axes
        ax.set_xlim(-screen_width/2, screen_width/2)
        ax.set_ylim(0, 1.05)
        ax.set_xlabel("Position on Screen (m)")
        ax.set_ylabel("Relative Intensity")
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        # Add title
        ax.set_title(f"Diffraction Pattern from {num_slits} Slits\n"
                    f"Grating spacing: {grating_spacing*1e6:.1f} μm, Distance to screen: {distance_to_screen:.1f} m")
        
        # Add a 2D wall pattern visualization if requested
        if show_wall_pattern:
            DiffractionVisualizer.add_wall_pattern_visualization(
                fig, wavelengths, colors, grating_spacing, distance_to_screen, 
                screen_width, num_slits
            )
        
        plt.tight_layout()
        return fig, ax
    
    @staticmethod
    def plot_infinite_slit_pattern(
        wavelengths: List[float],
        colors: List[str],
        labels: List[str],
        grating_spacing: float,
        distance_to_screen: float,
        screen_width: float,
        ax: Optional[plt.Axes] = None,
        show_wall_pattern: bool = True
    ) -> Tuple[plt.Figure, plt.Axes]:
        """
        Plot diffraction patterns for multiple wavelengths with infinite slits.
        
        Args:
            wavelengths: List of light wavelengths in meters
            colors: List of colors for each wavelength
            labels: List of labels for each wavelength
            grating_spacing: Distance between slits in meters
            distance_to_screen: Distance from grating to screen in meters
            screen_width: Width of the screen in meters
            ax: Optional axes to plot on (will create new if None)
            show_wall_pattern: Whether to show the 2D wall pattern visualization
            
        Returns:
            Figure and axes objects
        """
        # Create new figure if needed
        if ax is None:
            fig, ax = plt.subplots(figsize=(12, 8))
        else:
            fig = ax.figure
        
        # Plot each wavelength
        for i, wavelength in enumerate(wavelengths):
            screen_positions, intensity = DiffractionModel.calculate_infinite_slit_pattern(
                wavelength, grating_spacing, distance_to_screen, screen_width
            )
            
            ax.plot(screen_positions, intensity, color=colors[i], label=labels[i], linewidth=2)
        
        # Set up the axes
        ax.set_xlim(-screen_width/2, screen_width/2)
        ax.set_ylim(0, 1.05)
        ax.set_xlabel("Position on Screen (m)")
        ax.set_ylabel("Relative Intensity")
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        # Add title
        ax.set_title(f"Diffraction Pattern from Infinite Slits (Perfect Grating)\n"
                    f"Grating spacing: {grating_spacing*1e6:.1f} μm, Distance to screen: {distance_to_screen:.1f} m")
        
        # Add a 2D wall pattern visualization if requested
        if show_wall_pattern:
            DiffractionVisualizer.add_wall_pattern_visualization_infinite(
                fig, wavelengths, colors, grating_spacing, distance_to_screen, 
                screen_width
            )
        
        plt.tight_layout()
        return fig, ax
    
    @staticmethod
    def add_wall_pattern_visualization(
        fig: plt.Figure,
        wavelengths: List[float],
        colors: List[str],
        grating_spacing: float,
        distance_to_screen: float,
        screen_width: float,
        num_slits: int
    ) -> None:
        """
        Add a 2D visualization of the pattern on the wall.
        
        Args:
            fig: Figure to add the visualization to
            wavelengths: List of wavelengths
            colors: List of colors for each wavelength
            grating_spacing: Grating spacing in meters
            distance_to_screen: Distance to screen in meters
            screen_width: Screen width in meters
            num_slits: Number of slits
        """
        # Add inset axes for the wall pattern
        ax_inset = fig.add_axes([0.35, 0.02, 0.3, 0.25])  # [left, bottom, width, height]
        
        # Create a 2D image representing the pattern on the wall
        wall_height = screen_width  # Make it square for simplicity
        wall_y = np.linspace(-wall_height/2, wall_height/2, 500)
        
        # For each wavelength, calculate the pattern and add to the wall visualization
        wall_positions = np.linspace(-screen_width/2, screen_width/2, 1000)
        wall_pattern = np.zeros((len(wall_y), len(wall_positions), 3))  # RGB array
        
        for i, wavelength in enumerate(wavelengths):
            screen_positions, intensity = DiffractionModel.calculate_intensity_pattern(
                wavelength, grating_spacing, distance_to_screen, screen_width, num_slits=num_slits
            )
            
            # Interpolate intensity to match wall positions
            intensity_interp = interp1d(screen_positions, intensity, kind='linear', 
                                      bounds_error=False, fill_value=0)
            intensity_wall = intensity_interp(wall_positions)
            
            # Create 2D pattern (same in vertical direction)
            pattern_2d = np.tile(intensity_wall, (len(wall_y), 1))
            
            # Get RGB values for the current color
            if colors[i] in COLOR_RGB_MAP:
                rgb = COLOR_RGB_MAP[colors[i]]
            else:
                rgb = [1, 1, 1]  # Default to white
            
            # Add this wavelength's contribution to the wall pattern
            for j in range(3):  # RGB channels
                wall_pattern[:, :, j] += pattern_2d * rgb[j] / len(wavelengths)
        
        # Clip values to valid range [0, 1]
        wall_pattern = np.clip(wall_pattern, 0, 1)
        
        # Display the 2D pattern
        ax_inset.imshow(wall_pattern, extent=[-screen_width/2, screen_width/2, 
                                            -wall_height/2, wall_height/2],
                      aspect='auto', interpolation='bilinear')
        ax_inset.set_title("Pattern on Wall")
        ax_inset.set_xlabel("Horizontal Position (m)")
        ax_inset.set_ylabel("Vertical Position (m)")
    
    @staticmethod
    def add_wall_pattern_visualization_infinite(
        fig: plt.Figure,
        wavelengths: List[float],
        colors: List[str],
        grating_spacing: float,
        distance_to_screen: float,
        screen_width: float
    ) -> None:
        """
        Add a 2D visualization of the pattern on the wall for infinite slits.
        
        Args:
            fig: Figure to add the visualization to
            wavelengths: List of wavelengths
            colors: List of colors for each wavelength
            grating_spacing: Grating spacing in meters
            distance_to_screen: Distance to screen in meters
            screen_width: Screen width in meters
        """
        # Add inset axes for the wall pattern
        ax_inset = fig.add_axes([0.35, 0.02, 0.3, 0.25])  # [left, bottom, width, height]
        
        # Create a 2D image representing the pattern on the wall
        wall_height = screen_width  # Make it square for simplicity
        wall_y = np.linspace(-wall_height/2, wall_height/2, 500)
        
        # For each wavelength, calculate the pattern and add to the wall visualization
        wall_positions = np.linspace(-screen_width/2, screen_width/2, 1000)
        wall_pattern = np.zeros((len(wall_y), len(wall_positions), 3))  # RGB array
        
        for i, wavelength in enumerate(wavelengths):
            screen_positions, intensity = DiffractionModel.calculate_infinite_slit_pattern(
                wavelength, grating_spacing, distance_to_screen, screen_width
            )
            
            # Interpolate intensity to match wall positions
            intensity_interp = interp1d(screen_positions, intensity, kind='linear', 
                                      bounds_error=False, fill_value=0)
            intensity_wall = intensity_interp(wall_positions)
            
            # Create 2D pattern (same in vertical direction)
            pattern_2d = np.tile(intensity_wall, (len(wall_y), 1))
            
            # Get RGB values for the current color
            if colors[i] in COLOR_RGB_MAP:
                rgb = COLOR_RGB_MAP[colors[i]]
            else:
                rgb = [1, 1, 1]  # Default to white
            
            # Add this wavelength's contribution to the wall pattern
            for j in range(3):  # RGB channels
                wall_pattern[:, :, j] += pattern_2d * rgb[j] / len(wavelengths)
        
        # Clip values to valid range [0, 1]
        wall_pattern = np.clip(wall_pattern, 0, 1)
        
        # Display the 2D pattern
        ax_inset.imshow(wall_pattern, extent=[-screen_width/2, screen_width/2, 
                                            -wall_height/2, wall_height/2],
                      aspect='auto', interpolation='bilinear')
        ax_inset.set_title("Pattern on Wall")
        ax_inset.set_xlabel("Horizontal Position (m)")
        ax_inset.set_ylabel("Vertical Position (m)")
    
    @staticmethod
    def plot_total_intensity(
        wavelengths: List[float],
        grating_spacing: float,
        distance_to_screen: float,
        screen_width: float,
        num_slits: int,
        ax: Optional[plt.Axes] = None
    ) -> Tuple[plt.Figure, plt.Axes]:
        """
        Plot the total intensity from all wavelengths combined.
        
        Args:
            wavelengths: List of wavelengths in meters
            grating_spacing: Grating spacing in meters
            distance_to_screen: Distance to screen in meters
            screen_width: Screen width in meters
            num_slits: Number of slits
            ax: Optional axes to plot on (will create new if None)
            
        Returns:
            Figure and axes objects
        """
        # Create new figure if needed
        if ax is None:
            fig, ax = plt.subplots(figsize=(12, 3))
        else:
            fig = ax.figure
        
        # Calculate total intensity
        screen_positions, total_intensity = DiffractionModel.calculate_total_intensity(
            wavelengths, grating_spacing, distance_to_screen, screen_width, num_slits
        )
        
        # Plot total intensity
        ax.plot(screen_positions, total_intensity, 'k-', linewidth=2)
        ax.set_xlim(-screen_width/2, screen_width/2)
        ax.set_ylim(0, 1.05)
        ax.set_xlabel("Position on Screen (m)")
        ax.set_ylabel("Total Intensity")
        ax.grid(True, alpha=0.3)
        ax.set_title("Combined Intensity (Sum of All Selected Wavelengths)")
        
        return fig, ax 