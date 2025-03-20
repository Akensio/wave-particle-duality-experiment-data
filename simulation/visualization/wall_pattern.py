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
        ax: Optional[plt.Axes] = None
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
        ax.imshow(wall_pattern, extent=[-screen_width/2, screen_width/2, -wall_height, wall_height],
                 aspect='auto', interpolation='bilinear')
        ax.set_xlabel("Position (m)")
        ax.set_ylabel("Height (m)")
        ax.set_title("Pattern on Wall") 