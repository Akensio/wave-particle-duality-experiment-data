#!/usr/bin/env python3
"""
Core physics calculations for the diffraction grating simulation.

This module contains the physics models for simulating diffraction patterns
through single and multiple slits, implementing wave optics principles.
"""

import numpy as np
from typing import Tuple, List, Optional, Dict, Union, Protocol
from simulation.core.config import DEFAULT_NUM_POINTS


class DiffractionPhysics:
    """Base class for diffraction physics calculations."""
    
    def __init__(self, grating_spacing: float, distance_to_screen: float, screen_width: float):
        """
        Initialize the diffraction physics model.
        
        Args:
            grating_spacing: Distance between slits in meters
            distance_to_screen: Distance from grating to screen in meters
            screen_width: Width of the screen in meters
        """
        self.grating_spacing = grating_spacing
        self.distance_to_screen = distance_to_screen
        self.screen_width = screen_width
    
    def calculate_intensity_pattern(self, 
                                   wavelength: float, 
                                   num_points: int = DEFAULT_NUM_POINTS,
                                   num_slits: int = 5) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate the intensity pattern on a screen from a diffraction grating.
        
        Args:
            wavelength: Light wavelength in meters
            num_points: Number of points to calculate on the screen
            num_slits: Number of slits in the grating
            
        Returns:
            Tuple of (screen positions, intensity values)
        """
        # Calculate positions on the screen
        screen_positions = np.linspace(-self.screen_width/2, self.screen_width/2, num_points)
        
        # Calculate angles to each position on the screen
        angles = np.arctan(screen_positions / self.distance_to_screen)
        
        # Calculate the phase difference between adjacent slits
        delta = (2 * np.pi / wavelength) * self.grating_spacing * np.sin(angles)
        
        # Calculate the multiple slit interference factor
        intensity_factor = np.sin(num_slits * delta / 2)**2 / np.sin(delta / 2)**2
        # Replace NaN values that occur when sin(delta/2) is zero
        intensity_factor = np.nan_to_num(intensity_factor, nan=num_slits**2)
        
        # Normalize the intensity
        intensity = intensity_factor / np.max(intensity_factor)
        
        return screen_positions, intensity
    
    def calculate_maxima_positions(self,
                                  wavelength: float,
                                  max_order: int = 3) -> List[Tuple[int, float]]:
        """
        Calculate positions of diffraction maxima based on the grating equation.
        
        Args:
            wavelength: Light wavelength in meters
            max_order: Maximum diffraction order to calculate
            
        Returns:
            List of tuples (order, position) for each diffraction maximum that appears on screen
        """
        maxima = []
        for m in range(-max_order, max_order + 1):
            # Check if the maximum is physically possible
            if abs(m * wavelength / self.grating_spacing) < 1:
                angle = np.arcsin(m * wavelength / self.grating_spacing)
                pos = self.distance_to_screen * np.tan(angle)
                # Check if the maximum appears on the screen
                if abs(pos) <= self.screen_width/2:
                    maxima.append((m, pos))
        
        return maxima
    
    def calculate_total_intensity(self,
                                 wavelengths: List[float],
                                 num_slits: int,
                                 num_points: int = DEFAULT_NUM_POINTS) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate the total intensity from multiple wavelengths.
        
        Args:
            wavelengths: List of wavelengths in meters
            num_slits: Number of slits in the grating
            num_points: Number of points to calculate
            
        Returns:
            Tuple of (screen positions, total intensity)
        """
        if not wavelengths:
            # If no wavelengths are provided, return zeros
            screen_positions = np.linspace(-self.screen_width/2, self.screen_width/2, num_points)
            return screen_positions, np.zeros_like(screen_positions)
        
        # Get positions from the first wavelength
        screen_positions, _ = self.calculate_intensity_pattern(
            wavelengths[0], num_points=num_points, num_slits=num_slits
        )
        
        # Sum intensities from all wavelengths
        total_intensity = np.zeros_like(screen_positions)
        for wavelength in wavelengths:
            _, intensity = self.calculate_intensity_pattern(
                wavelength, num_points=num_points, num_slits=num_slits
            )
            total_intensity += intensity
        
        # Normalize the total intensity
        if np.max(total_intensity) > 0:
            total_intensity = total_intensity / np.max(total_intensity)
            
        return screen_positions, total_intensity


class InfiniteSlitDiffractionPhysics(DiffractionPhysics):
    """Physics model for diffraction through an infinite number of slits (perfect grating)."""
    
    def calculate_intensity_pattern(self,
                                   wavelength: float,
                                   num_points: int = DEFAULT_NUM_POINTS,
                                   num_slits: Optional[int] = None) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate the diffraction pattern for an infinite number of slits (perfect grating).
        This produces sharp peaks at the locations given by the grating equation.
        
        Args:
            wavelength: Light wavelength in meters
            num_points: Number of points to calculate
            num_slits: Ignored parameter (for API compatibility)
            
        Returns:
            Tuple of (screen positions, intensity values)
        """
        # Calculate positions on the screen
        screen_positions = np.linspace(-self.screen_width/2, self.screen_width/2, num_points)
        
        # Calculate angles to each position on the screen
        angles = np.arctan(screen_positions / self.distance_to_screen)
        
        # Calculate the expected maxima positions
        maxima = self.calculate_maxima_positions(wavelength, max_order=10)
        
        # Initialize intensity array
        intensity = np.zeros(num_points)
        
        # Add peaks at each maximum location
        for m, pos in maxima:
            # Find the closest index in screen_positions to pos
            idx = np.abs(screen_positions - pos).argmin()
            
            # Create a Gaussian peak around that position
            # Width of the peak is inversely proportional to the number of slits
            # Here we use a very small width to simulate sharp peaks (infinite slits)
            width = 0.01  # This is in meters and defines the sharpness of the peak
            intensity += np.exp(-0.5 * ((screen_positions - pos) / width) ** 2)
        
        # Normalize the intensity if there are any maxima
        if maxima:
            intensity = intensity / np.max(intensity)
            
        return screen_positions, intensity


# Legacy static adapter for backward compatibility
class DiffractionModel:
    """
    Legacy static adapter to maintain backward compatibility.
    
    Provides the original static methods from the old implementation, but delegates to
    the new object-oriented implementations.
    """
    
    @staticmethod
    def calculate_intensity_pattern(
        wavelength: float,
        grating_spacing: float, 
        distance_to_screen: float, 
        screen_width: float, 
        num_points: int = DEFAULT_NUM_POINTS,
        num_slits: int = 5
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Legacy adapter for the old static method."""
        model = DiffractionPhysics(grating_spacing, distance_to_screen, screen_width)
        return model.calculate_intensity_pattern(wavelength, num_points, num_slits)
    
    @staticmethod
    def calculate_maxima_positions(
        wavelength: float,
        grating_spacing: float,
        distance_to_screen: float,
        screen_width: float,
        max_order: int = 3
    ) -> List[Tuple[int, float]]:
        """Legacy adapter for the old static method."""
        model = DiffractionPhysics(grating_spacing, distance_to_screen, screen_width)
        return model.calculate_maxima_positions(wavelength, max_order)
    
    @staticmethod
    def calculate_total_intensity(
        wavelengths: List[float],
        grating_spacing: float,
        distance_to_screen: float,
        screen_width: float,
        num_slits: int,
        num_points: int = DEFAULT_NUM_POINTS
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Legacy adapter for the old static method."""
        model = DiffractionPhysics(grating_spacing, distance_to_screen, screen_width)
        return model.calculate_total_intensity(wavelengths, num_slits, num_points)
    
    @staticmethod
    def calculate_infinite_slit_pattern(
        wavelength: float,
        grating_spacing: float,
        distance_to_screen: float,
        screen_width: float,
        num_points: int = DEFAULT_NUM_POINTS
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Legacy adapter for the old static method."""
        model = InfiniteSlitDiffractionPhysics(grating_spacing, distance_to_screen, screen_width)
        return model.calculate_intensity_pattern(wavelength, num_points) 