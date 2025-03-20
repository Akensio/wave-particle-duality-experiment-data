import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from typing import List, Union, Callable, Tuple

from simulation.core.physics import DiffractionPhysics, InfiniteSlitDiffractionPhysics
from simulation.visualization.wall_pattern import WallPatternRenderer


class SpectrumPlotter:
    """Class for plotting custom spectrum visualizations."""
    
    def __init__(self):
        """Initialize spectrum plotter."""
        self.wall_renderer = WallPatternRenderer()
    
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
        """Plot a diffraction pattern with a custom spectrum."""
        # Filter wavelengths with significant intensity
        threshold = 0.01
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
            self.wall_renderer.add_spectrum_wall_pattern(
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
    
    def _calculate_combined_intensity(
        self,
        physics: Union[DiffractionPhysics, InfiniteSlitDiffractionPhysics],
        wavelengths: List[float],
        intensities: List[float],
        num_slits: int,
        use_infinite_slits: bool
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Calculate the combined intensity from multiple wavelengths."""
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