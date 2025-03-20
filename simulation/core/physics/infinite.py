import numpy as np
from typing import Tuple, Optional
from simulation.core.physics.base import DiffractionPhysics
from simulation.core.config import DEFAULT_NUM_POINTS


class InfiniteSlitDiffractionPhysics(DiffractionPhysics):
    """Physics model for infinite slit diffraction."""
    
    def calculate_intensity_pattern(self, wavelength, num_points=DEFAULT_NUM_POINTS, num_slits=None):
        """Calculate diffraction pattern for infinite slits."""
        screen_positions = np.linspace(-self.screen_width/2, self.screen_width/2, num_points)
        
        # Calculate maxima positions
        maxima = self.calculate_maxima_positions(wavelength, max_order=10)
        
        # Create intensity pattern with sharp peaks
        intensity = np.zeros(num_points)
        for m, pos in maxima:
            # Create a Gaussian peak around that position
            width = 0.01  # Defines sharpness of the peak
            intensity += np.exp(-0.5 * ((screen_positions - pos) / width) ** 2)
        
        # Normalize the intensity if there are any maxima
        if maxima:
            intensity = intensity / np.max(intensity)
            
        return screen_positions, intensity 