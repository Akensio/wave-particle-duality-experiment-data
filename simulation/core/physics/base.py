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
        
        # For large numbers of slits, use a more efficient delta function approximation
        if num_slits > 50:
            return self._calculate_large_slit_pattern(wavelength, screen_positions, num_slits)
        
        # Standard calculation for smaller numbers of slits
        angles = np.arctan(screen_positions / self.distance_to_screen)
        delta = (2 * np.pi / wavelength) * self.grating_spacing * np.sin(angles)
        
        
        intensity_factor = np.ones_like(delta)
        
        numerator = np.sin(num_slits * delta / 2)**2
        denominator = np.sin(delta / 2)**2
        intensity_factor = numerator / denominator
        
        # Normalize intensity
        intensity = intensity_factor / np.max(intensity_factor)
        
        return screen_positions, intensity
    
    def _calculate_large_slit_pattern(self, wavelength, screen_positions, num_slits):
        """
        Efficient calculation for large numbers of slits.
        For very large N, the pattern approaches delta functions at the maxima positions.
        """
        # Calculate maxima positions up to a sufficient order
        max_order = 10
        maxima = self.calculate_maxima_positions(wavelength, max_order)
        
        # Start with zero intensity
        intensity = np.zeros_like(screen_positions)
        
        # Width of peaks decreases with increasing number of slits
        # This approximates the narrowing of peaks as slit count increases
        peak_width = min(0.01, 0.1 / np.sqrt(num_slits))  # Width scales with 1/sqrt(N)
        
        # Add sharp Gaussian peaks at each maximum location
        for m, pos in maxima:
            # Create a sharp peak around the maximum position
            intensity += np.exp(-0.5 * ((screen_positions - pos) / peak_width) ** 2)
            
        # Normalize
        if np.max(intensity) > 0:
            intensity = intensity / np.max(intensity)
            
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