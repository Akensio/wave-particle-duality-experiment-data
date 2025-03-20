#!/usr/bin/env python3
# Simulation module for predefined diffraction grating scenarios

import numpy as np
from typing import List, Dict, Tuple, Optional, Union, Callable

from simulation.core.physics import DiffractionModel
from simulation.core.config import WAVELENGTH_OPTIONS, COLOR_MAP, DEFAULT_NUM_POINTS
from simulation.visualization.plotter import DiffractionVisualizer

class DiffractionSimulator:
    """
    Class for simulating various diffraction grating scenarios.
    This provides an easy way to set up common experiments.
    """
    
    @staticmethod
    def simulate_visible_spectrum(
        grating_spacing: float = 1.67e-6,  # 600 lines/mm
        distance_to_screen: float = 1.0,
        screen_width: float = 1.0,
        num_slits: int = 100,
        show_maxima: bool = True,
        use_infinite_slits: bool = False
    ) -> None:
        """
        Simulate a diffraction pattern with the full visible spectrum.
        
        Args:
            grating_spacing: Distance between slits in meters
            distance_to_screen: Distance from grating to screen in meters
            screen_width: Width of the screen in meters
            num_slits: Number of slits in the grating
            show_maxima: Whether to show diffraction maxima positions
            use_infinite_slits: Whether to use the infinite slits model
        """
        # Get all visible wavelengths
        colors = list(WAVELENGTH_OPTIONS.keys())
        wavelengths = [WAVELENGTH_OPTIONS[color] for color in colors]
        plot_colors = [COLOR_MAP[color] for color in colors]
        labels = [f"{color} ({wavelength*1e9:.0f} nm)" for color, wavelength in zip(colors, wavelengths)]
        
        # Create the plots
        if use_infinite_slits:
            # Use infinite slits model
            DiffractionVisualizer.plot_infinite_slit_pattern(
                wavelengths=wavelengths,
                colors=plot_colors,
                labels=labels,
                grating_spacing=grating_spacing,
                distance_to_screen=distance_to_screen,
                screen_width=screen_width
            )
        else:
            # Use finite slits model
            DiffractionVisualizer.plot_diffraction_pattern(
                wavelengths=wavelengths,
                colors=plot_colors,
                labels=labels,
                grating_spacing=grating_spacing,
                distance_to_screen=distance_to_screen,
                screen_width=screen_width,
                num_slits=num_slits,
                show_maxima=show_maxima
            )
    
    @staticmethod
    def simulate_monochromatic(
        wavelength: float = 532e-9,  # Green laser
        grating_spacing: float = 1.67e-6,  # 600 lines/mm
        distance_to_screen: float = 1.0,
        screen_width: float = 1.0,
        num_slits: int = 100,
        show_maxima: bool = True,
        use_infinite_slits: bool = False
    ) -> None:
        """
        Simulate a diffraction pattern with monochromatic light.
        
        Args:
            wavelength: Light wavelength in meters
            grating_spacing: Distance between slits in meters
            distance_to_screen: Distance from grating to screen in meters
            screen_width: Width of the screen in meters
            num_slits: Number of slits in the grating
            show_maxima: Whether to show diffraction maxima positions
            use_infinite_slits: Whether to use the infinite slits model
        """
        # Format wavelength as a label
        label = f"{wavelength*1e9:.0f} nm"
        
        # Determine color based on wavelength
        if wavelength < 450e-9:
            color = 'violet'
        elif wavelength < 495e-9:
            color = 'b'
        elif wavelength < 570e-9:
            color = 'g'
        elif wavelength < 590e-9:
            color = 'gold'
        elif wavelength < 620e-9:
            color = 'darkorange'
        else:
            color = 'r'
        
        # Create the plots
        if use_infinite_slits:
            # Use infinite slits model
            DiffractionVisualizer.plot_infinite_slit_pattern(
                wavelengths=[wavelength],
                colors=[color],
                labels=[label],
                grating_spacing=grating_spacing,
                distance_to_screen=distance_to_screen,
                screen_width=screen_width
            )
        else:
            # Use finite slits model
            DiffractionVisualizer.plot_diffraction_pattern(
                wavelengths=[wavelength],
                colors=[color],
                labels=[label],
                grating_spacing=grating_spacing,
                distance_to_screen=distance_to_screen,
                screen_width=screen_width,
                num_slits=num_slits,
                show_maxima=show_maxima
            )
    
    @staticmethod
    def simulate_custom_spectrum(
        intensity_function: Callable[[np.ndarray], np.ndarray],
        wavelength_range: Tuple[float, float] = (400e-9, 700e-9),
        num_wavelengths: int = 10,
        grating_spacing: float = 1.67e-6,  # 600 lines/mm
        distance_to_screen: float = 1.0,
        screen_width: float = 1.0,
        num_slits: int = 100,
        use_infinite_slits: bool = False
    ) -> None:
        """
        Simulate a diffraction pattern with a custom spectrum of light.
        
        Args:
            intensity_function: Function that takes wavelengths and returns relative intensities
            wavelength_range: Tuple of (min_wavelength, max_wavelength) in meters
            num_wavelengths: Number of wavelength points to simulate
            grating_spacing: Distance between slits in meters
            distance_to_screen: Distance from grating to screen in meters
            screen_width: Width of the screen in meters
            num_slits: Number of slits in the grating
            use_infinite_slits: Whether to use the infinite slits model
        """
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec
        
        # Generate wavelengths and intensities
        wavelengths = np.linspace(wavelength_range[0], wavelength_range[1], num_wavelengths)
        intensities = intensity_function(wavelengths)
        
        # Normalize intensities
        if np.max(intensities) > 0:
            intensities = intensities / np.max(intensities)
        
        # Create colors based on wavelengths
        colors = []
        for wl in wavelengths:
            if wl < 450e-9:
                colors.append('violet')
            elif wl < 495e-9:
                colors.append('b')
            elif wl < 570e-9:
                colors.append('g')
            elif wl < 590e-9:
                colors.append('gold')
            elif wl < 620e-9:
                colors.append('darkorange')
            else:
                colors.append('r')
        
        # Create labels
        labels = [f"{wl*1e9:.0f} nm" for wl in wavelengths]
        
        # Create figure with spectrum and diffraction pattern
        fig = plt.figure(figsize=(14, 10))
        gs = GridSpec(2, 1, height_ratios=[1, 3], figure=fig)
        
        # Plot input spectrum
        ax_spectrum = fig.add_subplot(gs[0])
        ax_spectrum.bar(wavelengths*1e9, intensities, width=5, color=colors, alpha=0.7)
        ax_spectrum.set_xlabel("Wavelength (nm)")
        ax_spectrum.set_ylabel("Relative Intensity")
        ax_spectrum.set_title("Input Light Spectrum")
        ax_spectrum.grid(True, alpha=0.3)
        
        # Create weighted wavelengths based on intensities
        weighted_wavelengths = []
        weighted_colors = []
        weighted_labels = []
        
        for wl, color, label, intensity in zip(wavelengths, colors, labels, intensities):
            # Only include wavelengths with non-zero intensity
            if intensity > 0.01:  # Threshold to avoid very dim components
                weighted_wavelengths.append(wl)
                weighted_colors.append(color)
                weighted_labels.append(label)
        
        # Create diffraction pattern plot
        ax_diffraction = fig.add_subplot(gs[1])
        
        if use_infinite_slits:
            # Use infinite slits model
            DiffractionVisualizer.plot_infinite_slit_pattern(
                wavelengths=weighted_wavelengths,
                colors=weighted_colors,
                labels=weighted_labels,
                grating_spacing=grating_spacing,
                distance_to_screen=distance_to_screen,
                screen_width=screen_width,
                ax=ax_diffraction,
                show_wall_pattern=False
            )
        else:
            # Use finite slits model
            DiffractionVisualizer.plot_diffraction_pattern(
                wavelengths=weighted_wavelengths,
                colors=weighted_colors,
                labels=weighted_labels,
                grating_spacing=grating_spacing,
                distance_to_screen=distance_to_screen,
                screen_width=screen_width,
                num_slits=num_slits,
                ax=ax_diffraction,
                show_wall_pattern=False
            )
        
        # Add combined intensity plot at the bottom
        ax_total = plt.axes([0.125, 0.07, 0.775, 0.15])
        
        if use_infinite_slits:
            # Calculate total intensity for infinite slits
            screen_positions = np.linspace(-screen_width/2, screen_width/2, DEFAULT_NUM_POINTS)
            total_intensity = np.zeros_like(screen_positions)
            
            for wl, intensity in zip(weighted_wavelengths, intensities):
                if intensity > 0.01:  # Only include non-dim components
                    # Get intensity pattern for this wavelength
                    pos, pattern = DiffractionModel.calculate_infinite_slit_pattern(
                        wl, grating_spacing, distance_to_screen, screen_width
                    )
                    
                    # Interpolate to common grid
                    from scipy.interpolate import interp1d
                    pattern_interp = interp1d(pos, pattern, kind='linear', 
                                           bounds_error=False, fill_value=0)
                    
                    # Add to total intensity, weighted by the spectrum intensity
                    total_intensity += pattern_interp(screen_positions) * intensity
        else:
            # Calculate total intensity for finite slits
            weighted_wavelengths_with_intensities = []
            for wl, intensity in zip(wavelengths, intensities):
                if intensity > 0.01:  # Only include non-dim components
                    # Add wavelength multiple times based on intensity
                    weighted_wavelengths_with_intensities.extend([wl] * int(intensity * 100))
            
            # Calculate total pattern
            screen_positions, total_intensity = DiffractionModel.calculate_total_intensity(
                weighted_wavelengths_with_intensities, grating_spacing, distance_to_screen, 
                screen_width, num_slits
            )
        
        # Normalize
        if np.max(total_intensity) > 0:
            total_intensity = total_intensity / np.max(total_intensity)
            
        # Plot total intensity
        ax_total.plot(screen_positions, total_intensity, 'k-', linewidth=2)
        ax_total.set_xlim(-screen_width/2, screen_width/2)
        ax_total.set_ylim(0, 1.05)
        ax_total.set_xlabel("Position on Screen (m)")
        ax_total.set_ylabel("Total Intensity")
        ax_total.grid(True, alpha=0.3)
        ax_total.set_title("Combined Intensity Pattern")
        
        # Add a 2D wall pattern visualization
        wall_height = screen_width/2
        wall_y = np.linspace(-wall_height, wall_height, 500)
        wall_pattern = np.zeros((len(wall_y), len(screen_positions), 3))
        
        for wl, color, intensity in zip(wavelengths, colors, intensities):
            if intensity > 0.01:  # Only include non-dim components
                # Get intensity pattern
                if use_infinite_slits:
                    pos, pattern = DiffractionModel.calculate_infinite_slit_pattern(
                        wl, grating_spacing, distance_to_screen, screen_width
                    )
                else:
                    pos, pattern = DiffractionModel.calculate_intensity_pattern(
                        wl, grating_spacing, distance_to_screen, screen_width, num_slits=num_slits
                    )
                
                # Interpolate to common grid
                from scipy.interpolate import interp1d
                pattern_interp = interp1d(pos, pattern, kind='linear', 
                                       bounds_error=False, fill_value=0)
                pattern_wall = pattern_interp(screen_positions)
                
                # Create 2D pattern
                pattern_2d = np.tile(pattern_wall, (len(wall_y), 1))
                
                # Get RGB values based on color
                from simulation.core.config import COLOR_RGB_MAP
                if color in COLOR_RGB_MAP:
                    rgb = COLOR_RGB_MAP[color]
                else:
                    rgb = [1, 1, 1]  # Default to white
                
                # Add to wall pattern, weighted by intensity from spectrum
                for j in range(3):
                    wall_pattern[:, :, j] += pattern_2d * rgb[j] * intensity
        
        # Clip and normalize
        wall_pattern = np.clip(wall_pattern, 0, 1)
        
        # Add the wall pattern as an inset
        ax_wall = plt.axes([0.35, 0.25, 0.3, 0.15])
        ax_wall.imshow(wall_pattern, extent=[-screen_width/2, screen_width/2, 
                                           -wall_height, wall_height],
                     aspect='auto', interpolation='bilinear')
        ax_wall.set_title("Pattern on Wall")
        
        # Add overall title
        if use_infinite_slits:
            fig.suptitle(
                f"Diffraction Pattern with Custom Spectrum (Infinite Slits)\n"
                f"Grating spacing: {grating_spacing*1e6:.1f} μm, Distance: {distance_to_screen:.1f} m",
                fontsize=16
            )
        else:
            fig.suptitle(
                f"Diffraction Pattern with Custom Spectrum\n"
                f"Grating spacing: {grating_spacing*1e6:.1f} μm, Number of slits: {num_slits}, "
                f"Distance: {distance_to_screen:.1f} m",
                fontsize=16
            )
        
        plt.tight_layout()
        plt.subplots_adjust(top=0.9)
        plt.show()
    
    @staticmethod
    def simulate_comparing_slit_numbers(
        wavelength: float = 532e-9,  # Green laser
        grating_spacing: float = 1.67e-6,  # 600 lines/mm
        distance_to_screen: float = 1.0,
        screen_width: float = 1.0,
        slit_numbers: List[int] = [2, 5, 10, 50, 100, 500],
        show_infinite: bool = True
    ) -> None:
        """
        Compare diffraction patterns with different numbers of slits.
        
        Args:
            wavelength: Light wavelength in meters
            grating_spacing: Distance between slits in meters
            distance_to_screen: Distance from grating to screen in meters
            screen_width: Width of the screen in meters
            slit_numbers: List of numbers of slits to compare
            show_infinite: Whether to include the infinite slits model
        """
        import matplotlib.pyplot as plt
        
        # Determine color based on wavelength
        if wavelength < 450e-9:
            color = 'violet'
        elif wavelength < 495e-9:
            color = 'b'
        elif wavelength < 570e-9:
            color = 'g'
        elif wavelength < 590e-9:
            color = 'gold'
        elif wavelength < 620e-9:
            color = 'darkorange'
        else:
            color = 'r'
        
        # Create a figure with subplots
        num_plots = len(slit_numbers) + (1 if show_infinite else 0)
        fig, axes = plt.subplots(num_plots, 1, figsize=(12, 3*num_plots))
        
        # If there's only one plot, ensure axes is a list
        if num_plots == 1:
            axes = [axes]
        
        # Plot each slit number
        for i, num_slits in enumerate(slit_numbers):
            screen_positions, intensity = DiffractionModel.calculate_intensity_pattern(
                wavelength, grating_spacing, distance_to_screen, screen_width, num_slits=num_slits
            )
            
            axes[i].plot(screen_positions, intensity, color=color, linewidth=2)
            axes[i].set_xlim(-screen_width/2, screen_width/2)
            axes[i].set_ylim(0, 1.05)
            axes[i].set_ylabel("Intensity")
            axes[i].grid(True, alpha=0.3)
            axes[i].set_title(f"Number of Slits: {num_slits}")
            
            # Add diffraction maxima markers
            maxima = DiffractionModel.calculate_maxima_positions(
                wavelength, grating_spacing, distance_to_screen, screen_width
            )
            
            for m, pos in maxima:
                axes[i].axvline(x=pos, color=color, linestyle='--', alpha=0.3)
                axes[i].text(pos, 0.1, f"m={m}", ha='center', color='black')
        
        # Add infinite slits plot if requested
        if show_infinite:
            screen_positions, intensity = DiffractionModel.calculate_infinite_slit_pattern(
                wavelength, grating_spacing, distance_to_screen, screen_width
            )
            
            axes[-1].plot(screen_positions, intensity, color=color, linewidth=2)
            axes[-1].set_xlim(-screen_width/2, screen_width/2)
            axes[-1].set_ylim(0, 1.05)
            axes[-1].set_ylabel("Intensity")
            axes[-1].grid(True, alpha=0.3)
            axes[-1].set_title("Infinite Slits (Perfect Grating)")
        
        # Add common x-label
        fig.text(0.5, 0.02, "Position on Screen (m)", ha='center', fontsize=12)
        
        # Add title
        fig.suptitle(
            f"Comparison of Diffraction Patterns with Different Numbers of Slits\n"
            f"Wavelength: {wavelength*1e9:.0f} nm, Grating spacing: {grating_spacing*1e6:.1f} μm",
            fontsize=16
        )
        
        plt.tight_layout()
        plt.subplots_adjust(top=0.95, bottom=0.05)
        plt.show() 