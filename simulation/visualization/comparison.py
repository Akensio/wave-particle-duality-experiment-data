import matplotlib.pyplot as plt
from typing import List, Optional

from simulation.core.physics import DiffractionPhysics, InfiniteSlitDiffractionPhysics


class ComparisonPlotter:
    """Class for creating comparison plots of diffraction patterns."""
    
    def plot_slit_comparison(
        self,
        finite_physics: DiffractionPhysics,
        wavelength: float,
        color: str,
        slit_numbers: List[int],
        infinite_physics: Optional[InfiniteSlitDiffractionPhysics] = None
    ) -> None:
        """Create a plot comparing diffraction patterns for different numbers of slits."""
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