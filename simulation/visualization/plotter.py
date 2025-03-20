#!/usr/bin/env python3
"""
Visualization module for diffraction grating simulation.
Provides utilities for plotting diffraction patterns and spectrum visualizations.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple, Dict, Optional, Union, Any
from scipy.interpolate import interp1d

from simulation.core.physics import DiffractionPhysics, InfiniteSlitDiffractionPhysics
from simulation.core.config import COLOR_RGB_MAP, DEFAULT_NUM_POINTS


class DiffractionVisualizer:
    """Class responsible for visualizing diffraction patterns."""
    
    def plot_diffraction_pattern(
        self,
        physics: DiffractionPhysics,
        wavelengths: List[float],
        colors: List[str],
        labels: List[str],
        num_slits: int,
        show_maxima: bool = True,
        max_order: int = 3,
        ax: Optional[plt.Axes] = None,
        show_wall_pattern: bool = True
    ) -> Tuple[plt.Figure, plt.Axes]:
        """
        Plot diffraction patterns for multiple wavelengths.
        
        Args:
            physics: Physics model to use for calculations
            wavelengths: List of light wavelengths in meters
            colors: List of colors for each wavelength
            labels: List of labels for each wavelength
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
            screen_positions, intensity = physics.calculate_intensity_pattern(
                wavelength, num_slits=num_slits
            )
            
            ax.plot(screen_positions, intensity, color=colors[i], label=labels[i], linewidth=2)
            
            # Add markers for diffraction maxima if requested
            if show_maxima:
                maxima = physics.calculate_maxima_positions(wavelength, max_order)
                
                for m, pos in maxima:
                    ax.axvline(x=pos, color=colors[i], linestyle='--', alpha=0.3)
                    # Only label maxima for the first wavelength to avoid clutter
                    if i == 0:
                        ax.text(pos, 0.1, f"m={m}", ha='center', color='black')
        
        # Set up the axes
        ax.set_xlim(-physics.screen_width/2, physics.screen_width/2)
        ax.set_ylim(0, 1.05)
        ax.set_xlabel("Position on Screen (m)")
        ax.set_ylabel("Relative Intensity")
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        # Add title
        ax.set_title(f"Diffraction Pattern from {num_slits} Slits\n"
                    f"Grating spacing: {physics.grating_spacing*1e6:.1f} μm, "
                    f"Distance to screen: {physics.distance_to_screen:.1f} m")
        
        # Add a 2D wall pattern visualization if requested
        if show_wall_pattern:
            self._add_wall_pattern_visualization(
                fig, physics, wavelengths, colors, num_slits
            )
        
        plt.tight_layout()
        return fig, ax
    
    def plot_infinite_slit_pattern(
        self,
        physics: InfiniteSlitDiffractionPhysics,
        wavelengths: List[float],
        colors: List[str],
        labels: List[str],
        ax: Optional[plt.Axes] = None,
        show_wall_pattern: bool = True
    ) -> Tuple[plt.Figure, plt.Axes]:
        """
        Plot diffraction patterns for infinite slits.
        
        Args:
            physics: Infinite slit physics model to use for calculations
            wavelengths: List of light wavelengths in meters
            colors: List of colors for each wavelength
            labels: List of labels for each wavelength
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
            screen_positions, intensity = physics.calculate_intensity_pattern(wavelength)
            
            ax.plot(screen_positions, intensity, color=colors[i], label=labels[i], linewidth=2)
        
        # Set up the axes
        ax.set_xlim(-physics.screen_width/2, physics.screen_width/2)
        ax.set_ylim(0, 1.05)
        ax.set_xlabel("Position on Screen (m)")
        ax.set_ylabel("Relative Intensity")
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        # Add title
        ax.set_title(f"Diffraction Pattern from Infinite Slits (Perfect Grating)\n"
                   f"Grating spacing: {physics.grating_spacing*1e6:.1f} μm, "
                   f"Distance to screen: {physics.distance_to_screen:.1f} m")
        
        # Add a 2D wall pattern visualization if requested
        if show_wall_pattern:
            self._add_wall_pattern_visualization(
                fig, physics, wavelengths, colors
            )
        
        plt.tight_layout()
        return fig, ax
    
    def plot_custom_spectrum(
        self,
        physics: Union[DiffractionPhysics, InfiniteSlitDiffractionPhysics],
        wavelengths: List[float],
        intensities: np.ndarray,
        colors: List[str],
        labels: List[str],
        num_slits: int,
        use_infinite_slits: bool = False
    ) -> None:
        """
        Plot a diffraction pattern with a custom spectrum.
        
        Args:
            physics: Physics model to use for calculations
            wavelengths: List of wavelengths in meters
            intensities: Relative intensities for each wavelength
            colors: List of colors for each wavelength
            labels: List of labels for each wavelength
            num_slits: Number of slits (only used for finite slits)
            use_infinite_slits: Whether to use the infinite slits model
        """
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec
        
        # Filter wavelengths with significant intensity
        threshold = 0.01  # Threshold to avoid very dim components
        significant_indices = [i for i, intensity in enumerate(intensities) if intensity > threshold]
        
        filtered_wavelengths = [wavelengths[i] for i in significant_indices]
        filtered_colors = [colors[i] for i in significant_indices]
        filtered_labels = [labels[i] for i in significant_indices]
        filtered_intensities = [intensities[i] for i in significant_indices]
        
        # Create figure
        fig = plt.figure(figsize=(14, 10))
        gs = GridSpec(2, 1, height_ratios=[1, 3], figure=fig)
        
        # Plot input spectrum
        ax_spectrum = fig.add_subplot(gs[0])
        ax_spectrum.bar(np.array(wavelengths)*1e9, intensities, width=5, color=colors, alpha=0.7)
        ax_spectrum.set_xlabel("Wavelength (nm)")
        ax_spectrum.set_ylabel("Relative Intensity")
        ax_spectrum.set_title("Input Light Spectrum")
        ax_spectrum.grid(True, alpha=0.3)
        
        # Plot diffraction pattern
        ax_diffraction = fig.add_subplot(gs[1])
        
        if filtered_wavelengths:
            for i, (wl, color, label) in enumerate(zip(filtered_wavelengths, filtered_colors, filtered_labels)):
                # Calculate pattern for this wavelength
                screen_positions, intensity = physics.calculate_intensity_pattern(
                    wl, num_slits=num_slits if not use_infinite_slits else None
                )
                
                # Plot the individual wavelength pattern
                ax_diffraction.plot(screen_positions, intensity, color=color, label=label, linewidth=1, alpha=0.5)
            
            # Calculate and plot combined intensity
            screen_positions, total_intensity = self._calculate_combined_intensity(
                physics, filtered_wavelengths, filtered_intensities, num_slits, use_infinite_slits
            )
            
            ax_diffraction.plot(screen_positions, total_intensity, 'k-', linewidth=2, label='Combined')
            
            # Set up axes
            ax_diffraction.set_xlim(-physics.screen_width/2, physics.screen_width/2)
            ax_diffraction.set_ylim(0, 1.05)
            ax_diffraction.set_xlabel("Position on Screen (m)")
            ax_diffraction.set_ylabel("Relative Intensity")
            ax_diffraction.grid(True, alpha=0.3)
            ax_diffraction.legend()
            
            # Add title
            model_type = "Infinite Slits" if use_infinite_slits else f"{num_slits} Slits"
            ax_diffraction.set_title(f"Diffraction Pattern ({model_type})")
            
            # Add a 2D wall pattern
            self._add_spectrum_wall_pattern(
                fig, physics, filtered_wavelengths, filtered_colors, filtered_intensities, 
                num_slits, use_infinite_slits
            )
        
        # Add overall title
        fig.suptitle(
            f"Diffraction with Custom Spectrum\n"
            f"Grating spacing: {physics.grating_spacing*1e6:.1f} μm, "
            f"Distance: {physics.distance_to_screen:.1f} m",
            fontsize=16
        )
        
        plt.tight_layout()
        plt.subplots_adjust(top=0.9)
        plt.show()
    
    def plot_slit_comparison(
        self,
        finite_physics: DiffractionPhysics,
        wavelength: float,
        color: str,
        slit_numbers: List[int],
        infinite_physics: Optional[InfiniteSlitDiffractionPhysics] = None
    ) -> None:
        """
        Create a plot comparing diffraction patterns for different numbers of slits.
        
        Args:
            finite_physics: Physics model for finite slit calculations
            wavelength: Light wavelength in meters
            color: Color for the plots
            slit_numbers: List of slit numbers to compare
            infinite_physics: Optional physics model for infinite slits
        """
        import matplotlib.pyplot as plt
        
        # Create a figure with subplots
        num_plots = len(slit_numbers) + (1 if infinite_physics else 0)
        fig, axes = plt.subplots(num_plots, 1, figsize=(12, 3*num_plots))
        
        # If there's only one plot, ensure axes is a list
        if num_plots == 1:
            axes = [axes]
        
        # Plot each slit number
        for i, num_slits in enumerate(slit_numbers):
            screen_positions, intensity = finite_physics.calculate_intensity_pattern(
                wavelength, num_slits=num_slits
            )
            
            axes[i].plot(screen_positions, intensity, color=color, linewidth=2)
            axes[i].set_xlim(-finite_physics.screen_width/2, finite_physics.screen_width/2)
            axes[i].set_ylim(0, 1.05)
            axes[i].set_ylabel("Intensity")
            axes[i].grid(True, alpha=0.3)
            axes[i].set_title(f"Number of Slits: {num_slits}")
            
            # Add diffraction maxima markers
            maxima = finite_physics.calculate_maxima_positions(wavelength)
            
            for m, pos in maxima:
                axes[i].axvline(x=pos, color=color, linestyle='--', alpha=0.3)
                axes[i].text(pos, 0.1, f"m={m}", ha='center', color='black')
        
        # Add infinite slits plot if requested
        if infinite_physics:
            screen_positions, intensity = infinite_physics.calculate_intensity_pattern(wavelength)
            
            axes[-1].plot(screen_positions, intensity, color=color, linewidth=2)
            axes[-1].set_xlim(-infinite_physics.screen_width/2, infinite_physics.screen_width/2)
            axes[-1].set_ylim(0, 1.05)
            axes[-1].set_ylabel("Intensity")
            axes[-1].grid(True, alpha=0.3)
            axes[-1].set_title("Infinite Slits (Perfect Grating)")
        
        # Add common x-label
        fig.text(0.5, 0.02, "Position on Screen (m)", ha='center', fontsize=12)
        
        # Add title
        fig.suptitle(
            f"Comparison of Diffraction Patterns with Different Numbers of Slits\n"
            f"Wavelength: {wavelength*1e9:.0f} nm, "
            f"Grating spacing: {finite_physics.grating_spacing*1e6:.1f} μm",
            fontsize=16
        )
        
        plt.tight_layout()
        plt.subplots_adjust(top=0.95, bottom=0.05)
        plt.show()
    
    def _add_wall_pattern_visualization(
        self,
        fig: plt.Figure,
        physics: Union[DiffractionPhysics, InfiniteSlitDiffractionPhysics],
        wavelengths: List[float],
        colors: List[str],
        num_slits: Optional[int] = None
    ) -> None:
        """
        Add a 2D visualization of the pattern on a wall.
        
        Args:
            fig: Figure to add the visualization to
            physics: Physics model to use for calculations
            wavelengths: List of wavelengths in meters
            colors: List of colors for each wavelength
            num_slits: Number of slits (only for finite slits)
        """
        # Create inset axes for the wall pattern
        ax_wall = fig.add_axes([0.15, 0.02, 0.7, 0.15])
        
        # Calculate the wall pattern
        screen_width = physics.screen_width
        wall_height = screen_width / 2
        wall_y = np.linspace(-wall_height, wall_height, 200)
        
        screen_positions = np.linspace(-screen_width/2, screen_width/2, DEFAULT_NUM_POINTS)
        wall_pattern = np.zeros((len(wall_y), len(screen_positions), 3))
        
        for wl, color in zip(wavelengths, colors):
            # Get intensity pattern
            positions, pattern = physics.calculate_intensity_pattern(
                wl, num_slits=num_slits
            )
            
            # Interpolate to common grid
            pattern_interp = interp1d(positions, pattern, kind='linear', 
                                   bounds_error=False, fill_value=0)
            pattern_wall = pattern_interp(screen_positions)
            
            # Create 2D pattern
            pattern_2d = np.tile(pattern_wall, (len(wall_y), 1))
            
            # Add RGB components
            rgb = COLOR_RGB_MAP.get(color, [1, 1, 1])  # Default to white if color not found
            
            for j in range(3):
                wall_pattern[:, :, j] += pattern_2d * rgb[j]
        
        # Clip and normalize
        wall_pattern = np.clip(wall_pattern, 0, 1)
        
        # Display the wall pattern
        ax_wall.imshow(wall_pattern, extent=[-screen_width/2, screen_width/2, -wall_height, wall_height],
                     aspect='auto', interpolation='bilinear')
        ax_wall.set_xlabel("Position (m)")
        ax_wall.set_ylabel("Height (m)")
        ax_wall.set_title("Pattern on Wall")
    
    def _add_spectrum_wall_pattern(
        self,
        fig: plt.Figure,
        physics: Union[DiffractionPhysics, InfiniteSlitDiffractionPhysics],
        wavelengths: List[float],
        colors: List[str],
        intensities: List[float],
        num_slits: int,
        use_infinite_slits: bool
    ) -> None:
        """
        Add a 2D wall pattern for spectrum visualization.
        
        Args:
            fig: Figure to add to
            physics: Physics model to use
            wavelengths: Wavelengths in meters
            colors: Colors for each wavelength
            intensities: Intensity for each wavelength
            num_slits: Number of slits
            use_infinite_slits: Whether to use infinite slits model
        """
        # Create inset axes
        ax_wall = plt.axes([0.35, 0.25, 0.3, 0.15])
        
        # Calculate dimensions
        screen_width = physics.screen_width
        wall_height = screen_width / 2
        screen_positions = np.linspace(-screen_width/2, screen_width/2, 200)
        wall_y = np.linspace(-wall_height, wall_height, 200)
        
        # Initialize wall pattern
        wall_pattern = np.zeros((len(wall_y), len(screen_positions), 3))
        
        # Calculate pattern for each wavelength
        for wl, color, intensity in zip(wavelengths, colors, intensities):
            # Skip very low intensity components
            if intensity < 0.01:
                continue
                
            # Get intensity pattern
            positions, pattern = physics.calculate_intensity_pattern(
                wl, num_slits=None if use_infinite_slits else num_slits
            )
            
            # Interpolate to common grid
            pattern_interp = interp1d(positions, pattern, kind='linear',
                                   bounds_error=False, fill_value=0)
            pattern_wall = pattern_interp(screen_positions)
            
            # Create 2D pattern
            pattern_2d = np.tile(pattern_wall, (len(wall_y), 1))
            
            # Get RGB values
            rgb = COLOR_RGB_MAP.get(color, [1, 1, 1])  # Default to white
            
            # Add to wall pattern with intensity weighting
            for j in range(3):
                wall_pattern[:, :, j] += pattern_2d * rgb[j] * intensity
        
        # Clip and normalize
        wall_pattern = np.clip(wall_pattern, 0, 1)
        
        # Display
        ax_wall.imshow(wall_pattern, extent=[-screen_width/2, screen_width/2, -wall_height, wall_height],
                     aspect='auto', interpolation='bilinear')
        ax_wall.set_title("Pattern on Wall")
    
    def _calculate_combined_intensity(
        self,
        physics: Union[DiffractionPhysics, InfiniteSlitDiffractionPhysics],
        wavelengths: List[float],
        intensities: List[float],
        num_slits: int,
        use_infinite_slits: bool
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate the combined intensity from multiple wavelengths.
        
        Args:
            physics: Physics model to use
            wavelengths: List of wavelengths
            intensities: Intensity of each wavelength
            num_slits: Number of slits
            use_infinite_slits: Whether to use infinite slits model
            
        Returns:
            Tuple of (screen positions, total intensity)
        """
        # Get screen positions from first wavelength
        screen_positions, _ = physics.calculate_intensity_pattern(
            wavelengths[0], num_slits=None if use_infinite_slits else num_slits
        )
        
        # Initialize total intensity
        total_intensity = np.zeros_like(screen_positions)
        
        # Sum all wavelength contributions
        for wl, intensity in zip(wavelengths, intensities):
            # Skip very low intensity components
            if intensity < 0.01:
                continue
                
            # Get pattern for this wavelength
            _, pattern = physics.calculate_intensity_pattern(
                wl, num_slits=None if use_infinite_slits else num_slits
            )
            
            # Add to total, weighted by intensity
            total_intensity += pattern * intensity
        
        # Normalize
        if np.max(total_intensity) > 0:
            total_intensity = total_intensity / np.max(total_intensity)
            
        return screen_positions, total_intensity 