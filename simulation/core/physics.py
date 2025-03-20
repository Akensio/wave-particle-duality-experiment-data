#!/usr/bin/env python3
# Core physics calculations for the diffraction grating simulation

import numpy as np
from typing import Tuple, List, Optional, Dict, Union
from simulation.core.config import DEFAULT_NUM_POINTS

class DiffractionModel:
    """Class responsible for physical calculations in the diffraction simulation."""
    
    @staticmethod
    def calculate_intensity_pattern(
        wavelength: float,
        grating_spacing: float, 
        distance_to_screen: float, 
        screen_width: float, 
        num_points: int = DEFAULT_NUM_POINTS,
        num_slits: int = 5
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate the intensity pattern on a screen from a diffraction grating.
        
        Args:
            wavelength: Light wavelength in meters
            grating_spacing: Distance between slits in meters
            distance_to_screen: Distance from grating to screen in meters
            screen_width: Width of the screen in meters
            num_points: Number of points to calculate on the screen
            num_slits: Number of slits in the grating
            
        Returns:
            Tuple of (screen positions, intensity values)
        """
        # Calculate positions on the screen
        screen_positions = np.linspace(-screen_width/2, screen_width/2, num_points)
        
        # Calculate angles to each position on the screen
        angles = np.arctan(screen_positions / distance_to_screen)
        
        # Calculate the phase difference between adjacent slits
        delta = (2 * np.pi / wavelength) * grating_spacing * np.sin(angles)
        
        # Calculate the multiple slit interference factor
        intensity_factor = np.sin(num_slits * delta / 2)**2 / np.sin(delta / 2)**2
        # Replace NaN values that occur when sin(delta/2) is zero
        intensity_factor = np.nan_to_num(intensity_factor, nan=num_slits**2)
        
        # Normalize the intensity
        intensity = intensity_factor / np.max(intensity_factor)
        
        return screen_positions, intensity
    
    @staticmethod
    def calculate_maxima_positions(
        wavelength: float,
        grating_spacing: float,
        distance_to_screen: float,
        screen_width: float,
        max_order: int = 3
    ) -> List[Tuple[int, float]]:
        """
        Calculate positions of diffraction maxima based on the grating equation.
        
        Args:
            wavelength: Light wavelength in meters
            grating_spacing: Distance between slits in meters
            distance_to_screen: Distance from grating to screen in meters
            screen_width: Width of the screen in meters
            max_order: Maximum diffraction order to calculate
            
        Returns:
            List of tuples (order, position) for each diffraction maximum that appears on screen
        """
        maxima = []
        for m in range(-max_order, max_order + 1):
            # Check if the maximum is physically possible
            if abs(m * wavelength / grating_spacing) < 1:
                angle = np.arcsin(m * wavelength / grating_spacing)
                pos = distance_to_screen * np.tan(angle)
                # Check if the maximum appears on the screen
                if abs(pos) <= screen_width/2:
                    maxima.append((m, pos))
        
        return maxima
    
    @staticmethod
    def calculate_total_intensity(
        wavelengths: List[float],
        grating_spacing: float,
        distance_to_screen: float,
        screen_width: float,
        num_slits: int,
        num_points: int = DEFAULT_NUM_POINTS
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate the total intensity from multiple wavelengths.
        
        Args:
            wavelengths: List of wavelengths in meters
            grating_spacing: Distance between slits in meters
            distance_to_screen: Distance from grating to screen in meters
            screen_width: Width of the screen in meters
            num_slits: Number of slits in the grating
            num_points: Number of points to calculate
            
        Returns:
            Tuple of (screen positions, total intensity)
        """
        if not wavelengths:
            # If no wavelengths are provided, return zeros
            screen_positions = np.linspace(-screen_width/2, screen_width/2, num_points)
            return screen_positions, np.zeros_like(screen_positions)
        
        # Get positions from the first wavelength
        screen_positions, _ = DiffractionModel.calculate_intensity_pattern(
            wavelengths[0], grating_spacing, distance_to_screen, screen_width, 
            num_points=num_points, num_slits=num_slits
        )
        
        # Sum intensities from all wavelengths
        total_intensity = np.zeros_like(screen_positions)
        for wavelength in wavelengths:
            _, intensity = DiffractionModel.calculate_intensity_pattern(
                wavelength, grating_spacing, distance_to_screen, screen_width,
                num_points=num_points, num_slits=num_slits
            )
            total_intensity += intensity
        
        # Normalize the total intensity
        if np.max(total_intensity) > 0:
            total_intensity = total_intensity / np.max(total_intensity)
            
        return screen_positions, total_intensity

    @staticmethod
    def calculate_infinite_slit_pattern(
        wavelength: float,
        grating_spacing: float,
        distance_to_screen: float,
        screen_width: float,
        num_points: int = DEFAULT_NUM_POINTS
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate the diffraction pattern for an infinite number of slits (perfect grating).
        This produces sharp peaks at the locations given by the grating equation.
        
        Args:
            wavelength: Light wavelength in meters
            grating_spacing: Distance between slits in meters
            distance_to_screen: Distance from grating to screen in meters
            screen_width: Width of the screen in meters
            num_points: Number of points to calculate on the screen
            
        Returns:
            Tuple of (screen positions, intensity values)
        """
        # Calculate positions on the screen
        screen_positions = np.linspace(-screen_width/2, screen_width/2, num_points)
        
        # Calculate angles to each position on the screen
        angles = np.arctan(screen_positions / distance_to_screen)
        
        # For an infinite number of slits, the pattern consists of sharp peaks
        # at positions where d*sin(θ) = m*λ (the grating equation)
        # We'll approximate this with narrow Gaussian peaks
        
        # Initialize intensity array
        intensity = np.zeros_like(screen_positions)
        
        # Calculate maxima positions up to a high order
        max_order = 20  # Higher order for more peaks
        for m in range(-max_order, max_order + 1):
            # Check if the maximum is physically possible
            if abs(m * wavelength / grating_spacing) < 1:
                sin_theta_max = m * wavelength / grating_spacing
                theta_max = np.arcsin(sin_theta_max)
                pos_max = distance_to_screen * np.tan(theta_max)
                
                # Check if the maximum appears on the screen
                if abs(pos_max) <= screen_width/2:
                    # Create a narrow peak at this position
                    # The width of the peak is inversely proportional to the number of slits
                    # For an infinite number, we use a very small width
                    width = grating_spacing / 1000  # Very narrow peak
                    peak = np.exp(-((screen_positions - pos_max) / width)**2)
                    intensity += peak
        
        # Normalize the intensity
        if np.max(intensity) > 0:
            intensity = intensity / np.max(intensity)
            
        return screen_positions, intensity 