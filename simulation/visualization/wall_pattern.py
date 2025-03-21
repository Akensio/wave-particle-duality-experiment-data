import numpy as np
import matplotlib.pyplot as plt
from typing import List, Optional
from scipy.interpolate import interp1d

from simulation.core.physics import DiffractionPhysics
from simulation.core.config import COLOR_RGB_MAP, DEFAULT_NUM_POINTS


class WallPatternRenderer:
    """Class for rendering 2D wall pattern visualizations."""

    def add_wall_pattern(
        self,
        fig: plt.Figure,
        physics: DiffractionPhysics,
        wavelengths: List[float],
        colors: List[str],
        num_slits: int,
        ax: Optional[plt.Axes] = None,
    ) -> None:
        """
        Add a 2D visualization of the pattern on a wall.

        Args:
            fig: Figure to add the visualization to
            physics: Physics model to use for calculations
            wavelengths: List of wavelengths in meters
            colors: List of colors for each wavelength
            num_slits: Number of slits
            ax: Optional axes to plot on (will create new if None)
        """
        # Create inset axes for the wall pattern if not provided
        if ax is None:
            ax = fig.add_axes([0.15, 0.02, 0.7, 0.15])

        # Calculate the wall pattern on semi-circular screen
        # Use angular coordinates from -pi/2 to pi/2
        angular_width = np.pi  # Full angular width from -pi/2 to pi/2
        wall_height = angular_width / 4  # Height in the same scale as angular width
        wall_y = np.linspace(-wall_height, wall_height, 50)

        angles = np.linspace(-np.pi/2, np.pi/2, DEFAULT_NUM_POINTS)
        wall_pattern = np.zeros((len(wall_y), len(angles), 3))

        for wl, color in zip(wavelengths, colors):
            # Get intensity pattern
            theta, pattern = physics.calculate_intensity_pattern(wl, num_slits=num_slits)

            # Interpolate to common grid
            pattern_interp = interp1d(theta, pattern, kind="linear", bounds_error=False, fill_value=0)
            pattern_wall = pattern_interp(angles)

            # Create 2D pattern
            pattern_2d = np.tile(pattern_wall, (len(wall_y), 1))

            # Add RGB components
            rgb = COLOR_RGB_MAP.get(color, [1, 1, 1])  # Default to white if color not found

            for j in range(3):
                wall_pattern[:, :, j] += pattern_2d * rgb[j]

        # Clip and normalize
        wall_pattern = np.clip(wall_pattern, 0, 1)

        # Display the wall pattern
        ax.imshow(
            wall_pattern, extent=[-np.pi/2, np.pi/2, -wall_height, wall_height], aspect="auto", interpolation="bilinear"
        )
        ax.set_xlabel("Angle θ (radians)")
        ax.set_ylabel("Height")
        
        # Add a secondary x-axis with angles in degrees
        ax_deg = ax.twiny()
        ax_deg.set_xlim(-90, 90)
        ax_deg.set_xlabel("Angle θ (degrees)")
        
        ax.set_title("Pattern on Semi-Circular Screen")
