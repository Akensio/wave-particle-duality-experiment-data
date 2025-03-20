import numpy as np
import matplotlib.pyplot as plt
from typing import List, Union, Optional
from scipy.interpolate import interp1d

from simulation.core.physics import DiffractionPhysics, InfiniteSlitDiffractionPhysics
from simulation.core.config import COLOR_RGB_MAP, DEFAULT_NUM_POINTS


class WallPatternRenderer:
    """Class for rendering 2D wall pattern visualizations."""
    
    def add_wall_pattern(
        self,
        fig: plt.Figure,
        physics: Union[DiffractionPhysics, InfiniteSlitDiffractionPhysics],
        wavelengths: List[float],
        colors: List[str],
        num_slits: Optional[int] = None,
        ax: Optional[plt.Axes] = None
    ) -> None:
        """
        Add a 2D visualization of the pattern on a wall.
        
        Args:
            fig: Figure to add the visualization to
            physics: Physics model to use for calculations
            wavelengths: List of wavelengths in meters
            colors: List of colors for each wavelength
            num_slits: Number of slits (only for finite slits)
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
    
    def add_spectrum_wall_pattern(
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