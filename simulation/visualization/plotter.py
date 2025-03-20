import matplotlib.pyplot as plt
from typing import List, Tuple, Optional

from simulation.core.physics import DiffractionPhysics
from simulation.visualization.wall_pattern import WallPatternRenderer


class DiffractionVisualizer:
    """Class for visualizing diffraction patterns."""
    
    def __init__(self):
        """Initialize the diffraction visualizer."""
        self.wall_renderer = WallPatternRenderer()
    
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
        """Plot diffraction patterns for multiple wavelengths."""
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
            self.wall_renderer.add_wall_pattern(
                fig, physics, wavelengths, colors, num_slits
            )
        
        plt.tight_layout()
        plt.show()
        return fig, ax 