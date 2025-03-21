import numpy as np
from typing import Tuple, List, Optional, Callable
from numpy.typing import NDArray
from simulation.core.physics.base import DiffractionPhysics
from simulation.core.config import DEFAULT_NUM_POINTS

class ContinuousSpectrumPhysics(DiffractionPhysics):
    """Physics model for continuous wavelength spectrum diffraction."""
    
    def __init__(self, grating_spacing: float) -> None:
        """Initialize with grating spacing."""
        super().__init__(grating_spacing)
        
    def planck_spectrum(self, wavelengths: NDArray[np.float64], temperature: float) -> NDArray[np.float64]:
        """
        Calculate Planck's blackbody spectrum for a range of wavelengths.
        Returns normalized intensity values.
        
        Args:
            wavelengths: Array of wavelengths in meters
            temperature: Temperature in Kelvin
        
        Returns:
            Normalized intensity values for the given wavelengths
        """
        # Physical constants
        h = 6.626e-34  # Planck's constant (J·s)
        c = 2.998e8    # Speed of light (m/s)
        k_B = 1.381e-23  # Boltzmann constant (J/K)
        
        # Calculate Planck's law formula
        numerator = 2.0 * h * c**2
        denominator = wavelengths**5 * (np.exp((h*c)/(wavelengths*k_B*temperature)) - 1.0)
        intensity = numerator / denominator
        
        # Normalize the intensity
        if np.max(intensity) > 0:
            return intensity / np.max(intensity)
        return intensity
    
    def calculate_continuous_spectrum_pattern(
        self, 
        wavelength_range: Tuple[float, float], 
        num_wavelengths: int, 
        spectrum_function: Callable[[NDArray[np.float64]], NDArray[np.float64]],
        max_order: int = 1,
        num_points: int = DEFAULT_NUM_POINTS
    ) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
        """
        Calculate diffraction pattern for a continuous spectrum with a given spectral distribution.
        
        Args:
            wavelength_range: Tuple of (min_wavelength, max_wavelength) in meters
            num_wavelengths: Number of discrete wavelengths to sample in the range
            spectrum_function: Function that takes wavelengths and returns relative intensities
            max_order: Maximum diffraction order to consider
            num_points: Number of points in the output angle array
            
        Returns:
            Tuple of (angles, intensity) arrays
        """
        # Sample wavelengths from the range
        wavelengths = np.linspace(wavelength_range[0], wavelength_range[1], num_wavelengths)
        
        # Get spectrum intensities for each wavelength
        spectrum_intensities = spectrum_function(wavelengths)
        
        # Setup angle array
        angles = np.linspace(-np.pi/2, np.pi/2, num_points)
        
        # Initialize total intensity
        total_intensity = np.zeros_like(angles)
        
        # For each wavelength, calculate diffraction pattern and add to total
        for wavelength, spectral_intensity in zip(wavelengths, spectrum_intensities):
            # Only consider the specified diffraction order
            for m in range(1, max_order + 1):  # Starting from order 1, ignoring 0
                if m * wavelength / self.grating_spacing < 1:  # Physical constraint check
                    angle = np.arcsin(m * wavelength / self.grating_spacing)
                    
                    # Apply a narrow peak at the angle (similar to delta function approximation)
                    # Width of peak is narrower for higher-order diffraction
                    peak_width = 0.01 / m
                    
                    # Create a sharp Gaussian peak around the maximum position
                    peak_intensity = spectral_intensity * np.exp(-0.5 * ((angles - angle) / peak_width) ** 2)
                    total_intensity += peak_intensity
                    
                    # For negative orders (symmetry)
                    if m > 0:
                        neg_angle = np.arcsin(-m * wavelength / self.grating_spacing)
                        peak_intensity = spectral_intensity * np.exp(-0.5 * ((angles - neg_angle) / peak_width) ** 2)
                        total_intensity += peak_intensity
        
        # Normalize
        if np.max(total_intensity) > 0:
            total_intensity = total_intensity / np.max(total_intensity)
            
        return angles, total_intensity 