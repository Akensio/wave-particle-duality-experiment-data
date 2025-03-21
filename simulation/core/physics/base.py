import numpy as np
from typing import Tuple, List, Optional, Union
from simulation.core.config import DEFAULT_NUM_POINTS


class DiffractionPhysics:
    """Base diffraction physics model with semi-circular screen at infinite distance."""

    def __init__(self, grating_spacing: float, screen_width: Optional[float] = None, distance_to_screen: Optional[float] = None) -> None:
        self.grating_spacing = grating_spacing
        # The following parameters are no longer used but kept for compatibility
        self.distance_to_screen = float('inf')  # Assume infinite distance
        self.screen_width = np.pi  # Full angular range is now pi radians

    def calculate_intensity_pattern(self, wavelength: float, num_points: int = DEFAULT_NUM_POINTS, num_slits: int = 5) -> Tuple[np.ndarray, np.ndarray]:
        """Calculate diffraction intensity pattern on semi-circular screen using angular coordinates."""
        # Directly use angles from -pi/2 to +pi/2 as the x-axis
        angles: np.ndarray = np.linspace(-np.pi/2, np.pi/2, num_points)

        # For large numbers of slits, use a more efficient delta function approximation
        if num_slits > 50:
            return self._calculate_large_slit_pattern(wavelength, angles, num_slits)

        # Standard calculation for smaller numbers of slits
        delta: np.ndarray = (2 * np.pi / wavelength) * self.grating_spacing * np.sin(angles)

        numerator: np.ndarray = np.sin(num_slits * delta / 2) ** 2
        denominator: np.ndarray = np.sin(delta / 2) ** 2
        # Handle division by zero at positions where denominator is very small
        intensity_factor: np.ndarray = np.divide(numerator, denominator, out=np.zeros_like(numerator), where=denominator > 1e-10)

        # Normalize intensity
        intensity: np.ndarray = intensity_factor / np.max(intensity_factor) if np.max(intensity_factor) > 0 else intensity_factor

        return angles, intensity

    def _calculate_large_slit_pattern(self, wavelength: float, angles: np.ndarray, num_slits: int) -> Tuple[np.ndarray, np.ndarray]:
        """
        Efficient calculation for large numbers of slits using angles directly.
        For very large N, the pattern approaches delta functions at the maxima positions.
        """
        # Calculate maxima positions up to a sufficient order
        max_order: int = 10
        maxima: List[Tuple[int, float]] = self.calculate_maxima_positions(wavelength, max_order)

        # Start with zero intensity
        intensity: np.ndarray = np.zeros_like(angles)

        # Width of peaks decreases with increasing number of slits
        # This approximates the narrowing of peaks as slit count increases
        peak_width: float = min(0.01, 0.1 / np.sqrt(num_slits))  # Width scales with 1/sqrt(N)

        # Add sharp Gaussian peaks at each maximum location
        for m, angle in maxima:
            # Create a sharp peak around the maximum position
            intensity += np.exp(-0.5 * ((angles - angle) / peak_width) ** 2)

        # Normalize
        if np.max(intensity) > 0:
            intensity = intensity / np.max(intensity)

        return angles, intensity

    def calculate_maxima_positions(self, wavelength: float, max_order: int = 3) -> List[Tuple[int, float]]:
        """Calculate diffraction maxima positions in angle space."""
        maxima: List[Tuple[int, float]] = []
        for m in range(-max_order, max_order + 1):
            if abs(m * wavelength / self.grating_spacing) < 1:
                angle: float = np.arcsin(m * wavelength / self.grating_spacing)
                if abs(angle) <= np.pi/2:
                    maxima.append((m, angle))

        return maxima

    def calculate_total_intensity(self, wavelengths: List[float], num_slits: int, num_points: int = DEFAULT_NUM_POINTS) -> Tuple[np.ndarray, np.ndarray]:
        """Calculate total intensity pattern from multiple wavelengths."""
        if not wavelengths:
            angles: np.ndarray = np.linspace(-np.pi/2, np.pi/2, num_points)
            return angles, np.zeros_like(angles)

        angles, _ = self.calculate_intensity_pattern(wavelengths[0], num_points=num_points, num_slits=num_slits)

        total_intensity: np.ndarray = np.zeros_like(angles)
        for wavelength in wavelengths:
            _, intensity = self.calculate_intensity_pattern(wavelength, num_points=num_points, num_slits=num_slits)
            total_intensity += intensity

        if np.max(total_intensity) > 0:
            total_intensity = total_intensity / np.max(total_intensity)

        return angles, total_intensity
