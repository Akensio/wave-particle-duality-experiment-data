import numpy as np
from typing import Tuple, List, Optional, Callable
from numpy.typing import NDArray
from simulation.core.physics.base import DiffractionPhysics
from simulation.core.config import DEFAULT_NUM_POINTS
from scipy import constants

class ContinuousSpectrumPhysics(DiffractionPhysics):
    """Physics model for continuous wavelength spectrum diffraction."""
    
    def __init__(self, grating_spacing: float) -> None:
        """Initialize with grating spacing. Grating spacing comes in meters."""
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
        h = constants.Planck  # Planck's constant (J·s)
        c = constants.speed_of_light  # Speed of light (m/s)
        k_B = constants.Boltzmann  # Boltzmann constant (J/K)
        
        # Calculate Planck's law formula
        numerator = 2.0 * h * c**2
        denominator = wavelengths**5 * (np.exp((h*c)/(wavelengths*k_B*temperature)) - 1.0)
        intensity = numerator / denominator

        return intensity
    
    def calculate_order_m_pattern(
        self,
        angles: NDArray[np.float64],
        temperature: float,
        order_m: int
    ) -> NDArray[np.float64]:
        """
        Calculate diffraction pattern using equation 3 from eqs2.md directly.
        
        Args:
            angles: Array of angles in radians
            temperature: Temperature in Kelvin
            order_m: Diffraction order (m)
            
        Returns:
            Intensity array for the given angles and diffraction order
        """
        d = self.grating_spacing

        wavelengths = np.abs(np.sin(angles)) * d / order_m

        intensity = self.planck_spectrum(wavelengths, temperature)

        return intensity 