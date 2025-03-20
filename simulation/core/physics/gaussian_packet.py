import numpy as np
from typing import Tuple, List

from simulation.core.physics.base import DiffractionPhysics
from simulation.core.config import DEFAULT_NUM_POINTS


class GaussianWavePacketPhysics(DiffractionPhysics):
    """Physics model for Gaussian wave packet diffraction."""
    
    def __init__(self, grating_spacing, distance_to_screen, screen_width):
        """Initialize the Gaussian wave packet physics model."""
        super().__init__(grating_spacing, distance_to_screen, screen_width)
    
    def calculate_wave_packet_pattern(
        self,
        center_wavelength: float,
        wavelength_width: float,
        num_samples: int = 50,
        num_points: int = DEFAULT_NUM_POINTS,
        num_slits: int = 3
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate the diffraction pattern for a Gaussian wave packet.
        
        Args:
            center_wavelength: Center wavelength of the Gaussian packet (in meters)
            wavelength_width: Width of the Gaussian wavelength distribution (in meters)
            num_samples: Number of wavelength samples to use for the calculation
            num_points: Number of points to calculate on the screen
            num_slits: Number of slits in the diffraction grating
            
        Returns:
            Tuple of (screen positions, intensity pattern)
        """
        # Sample wavelengths from the Gaussian spectrum (±4 sigma covers >99.99% of the distribution)
        wavelengths = np.linspace(
            center_wavelength - 4*wavelength_width,
            center_wavelength + 4*wavelength_width,
            num_samples
        )
        
        # Ensure all wavelengths are positive
        wavelengths = wavelengths[wavelengths > 0]
        
        if len(wavelengths) == 0:
            screen_positions = np.linspace(-self.screen_width/2, self.screen_width/2, num_points)
            return screen_positions, np.zeros_like(screen_positions)
        
        # Calculate the Gaussian spectrum
        spectrum = np.exp(-0.5 * ((wavelengths - center_wavelength) / wavelength_width) ** 2)
        
        # Normalize the spectrum
        spectrum = spectrum / np.sum(spectrum)
        
        # Get the screen positions
        screen_positions, _ = self.calculate_intensity_pattern(
            wavelengths[0], num_points=num_points, num_slits=num_slits
        )
        
        # Calculate the total intensity pattern
        total_intensity = np.zeros_like(screen_positions)
        
        for i, wavelength in enumerate(wavelengths):
            _, intensity = self.calculate_intensity_pattern(
                wavelength, num_points=num_points, num_slits=num_slits
            )
            total_intensity += intensity * spectrum[i]
        
        # Normalize the total intensity
        if np.max(total_intensity) > 0:
            total_intensity = total_intensity / np.max(total_intensity)
            
        return screen_positions, total_intensity
    
    def calculate_wave_packet_time_evolution(
        self,
        center_wavelength: float,
        wavelength_width: float,
        times: List[float],
        position: float,
        dispersion_factor: float = 0.0
    ) -> Tuple[List[float], np.ndarray]:
        """
        Calculate the time evolution of a Gaussian wave packet at a specific position.
        
        Args:
            center_wavelength: Center wavelength of the Gaussian packet (in meters)
            wavelength_width: Width of the Gaussian wavelength distribution (in meters)
            times: List of time points to evaluate (in seconds)
            position: Position to evaluate the wave packet (in meters)
            dispersion_factor: Factor controlling dispersion (0 for no dispersion)
            
        Returns:
            Tuple of (times, wave packet amplitude at each time)
        """
        # Constants
        c = 3e8  # Speed of light in m/s
        center_k = 2 * np.pi / center_wavelength  # Central wavenumber
        center_omega = c * center_k  # Central angular frequency
        
        # Calculate the spatial width of the packet (uncertainty principle)
        sigma_x = 1 / (2 * np.pi * (1/center_wavelength - 1/(center_wavelength + wavelength_width)))
        
        # Calculate the wave packet amplitude at each time
        amplitudes = []
        for t in times:
            # With dispersion, the packet width increases with time
            sigma_t = sigma_x * (1 + dispersion_factor * t**2)
            
            # Wave packet amplitude (real part of the complex wave function)
            # A * exp(-(x-vt)²/(2σ²)) * cos(k₀x - ω₀t)
            envelope = np.exp(-((position - c*t)**2) / (2 * sigma_t**2))
            carrier = np.cos(center_k * position - center_omega * t)
            
            amplitudes.append(envelope * carrier)
            
        return times, np.array(amplitudes)
    
    def calculate_spectral_distribution(
        self,
        center_wavelength: float,
        wavelength_width: float,
        num_samples: int = 100
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate the spectral distribution of the Gaussian wave packet.
        
        Args:
            center_wavelength: Center wavelength of the Gaussian packet (in meters)
            wavelength_width: Width of the Gaussian wavelength distribution (in meters)
            num_samples: Number of wavelength samples
            
        Returns:
            Tuple of (wavelengths, spectral intensity)
        """
        # Sample wavelengths from the Gaussian spectrum (±4 sigma covers >99.99% of the distribution)
        wavelengths = np.linspace(
            max(0, center_wavelength - 4*wavelength_width),
            center_wavelength + 4*wavelength_width,
            num_samples
        )
        
        # Calculate the Gaussian spectrum
        spectrum = np.exp(-0.5 * ((wavelengths - center_wavelength) / wavelength_width) ** 2)
        
        # Normalize the spectrum
        if np.max(spectrum) > 0:
            spectrum = spectrum / np.max(spectrum)
            
        return wavelengths, spectrum 