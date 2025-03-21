import matplotlib.pyplot as plt
import numpy as np
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
        show_wall_pattern: bool = True,
    ) -> Tuple[plt.Figure, plt.Axes]:
        """Plot diffraction patterns for multiple wavelengths."""
        # Create new figure if needed
        if ax is None:
            fig, ax = plt.subplots(figsize=(12, 8))
        else:
            fig = ax.figure

        # Plot each wavelength
        for i, wavelength in enumerate(wavelengths):
            angles, intensity = physics.calculate_intensity_pattern(wavelength, num_slits=num_slits)

            ax.plot(angles, intensity, color=colors[i], label=labels[i], linewidth=2)

            # Add markers for diffraction maxima if requested
            if show_maxima:
                maxima = physics.calculate_maxima_positions(wavelength, max_order)

                for m, angle in maxima:
                    ax.axvline(x=angle, color=colors[i], linestyle="--", alpha=0.3)
                    # Only label maxima for the first wavelength to avoid clutter
                    if i == 0:
                        ax.text(angle, 0.1, f"m={m}", ha="center", color="black")

        # Set up the axes
        ax.set_xlim(-np.pi/2, np.pi/2)
        ax.set_ylim(0, 1.05)
        ax.set_xlabel("Angle θ (radians)")
        ax.set_ylabel("Relative Intensity")
        
        # Add a secondary x-axis with angles in degrees
        ax_deg = ax.twiny()
        ax_deg.set_xlim(-90, 90)
        
        ax.grid(True, alpha=0.3)
        ax.legend()

        # Add title
        ax.set_title(
            f"Diffraction Pattern from {num_slits} Slits\n"
            f"Grating spacing: {physics.grating_spacing*1e6:.1f} μm, "
            f"Semi-circular screen at infinite distance"
        )

        # Add a 2D wall pattern visualization if requested
        if show_wall_pattern:
            self.wall_renderer.add_wall_pattern(fig, physics, wavelengths, colors, num_slits)

        plt.tight_layout()
        plt.show()
        return fig, ax
