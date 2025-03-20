import numpy as np
from typing import Tuple, List
from simulation.core.config import DEFAULT_NUM_POINTS


class DiffractionPhysics:
    """Base diffraction physics model."""
    
    def __init__(self, grating_spacing, distance_to_screen, screen_width):
        self.grating_spacing = grating_spacing
        self.distance_to_screen = distance_to_screen
        self.screen_width = screen_width
    
    def calculate_intensity_pattern(self, wavelength, num_points=DEFAULT_NUM_POINTS, num_slits=5):
        """Calculate diffraction intensity pattern on screen."""
        screen_positions = np.linspace(-self.screen_width/2, self.screen_width/2, num_points)
        angles = np.arctan(screen_positions / self.distance_to_screen)
        
        delta = (2 * np.pi / wavelength) * self.grating_spacing * np.sin(angles)
        
        intensity_factor = np.sin(num_slits * delta / 2)**2 / np.sin(delta / 2)**2
        intensity_factor = np.nan_to_num(intensity_factor, nan=num_slits**2)
        
        intensity = intensity_factor / np.max(intensity_factor)
        
        return screen_positions, intensity
    
    def calculate_maxima_positions(self, wavelength, max_order=3):
        """Calculate diffraction maxima positions."""
        maxima = []
        for m in range(-max_order, max_order + 1):
            if abs(m * wavelength / self.grating_spacing) < 1:
                angle = np.arcsin(m * wavelength / self.grating_spacing)
                pos = self.distance_to_screen * np.tan(angle)
                if abs(pos) <= self.screen_width/2:
                    maxima.append((m, pos))
        
        return maxima
        
    def calculate_total_intensity(self, wavelengths, num_slits, num_points=DEFAULT_NUM_POINTS):
        """Calculate total intensity pattern from multiple wavelengths."""
        if not wavelengths:
            screen_positions = np.linspace(-self.screen_width/2, self.screen_width/2, num_points)
            return screen_positions, np.zeros_like(screen_positions)
        
        screen_positions, _ = self.calculate_intensity_pattern(
            wavelengths[0], num_points=num_points, num_slits=num_slits
        )
        
        total_intensity = np.zeros_like(screen_positions)
        for wavelength in wavelengths:
            _, intensity = self.calculate_intensity_pattern(
                wavelength, num_points=num_points, num_slits=num_slits
            )
            total_intensity += intensity
        
        if np.max(total_intensity) > 0:
            total_intensity = total_intensity / np.max(total_intensity)
            
        return screen_positions, total_intensity 