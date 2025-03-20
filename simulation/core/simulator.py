#!/usr/bin/env python3
"""
Simulation module for diffraction grating experiments with various configurations.
"""

import numpy as np
from typing import List, Dict, Tuple, Optional, Union, Callable

from simulation.core.physics import DiffractionPhysics, InfiniteSlitDiffractionPhysics
from simulation.core.config import WAVELENGTH_OPTIONS, COLOR_MAP, DEFAULT_NUM_POINTS
from simulation.visualization.plotter import DiffractionVisualizer


class DiffractionSimulator:
    """
    Class for simulating diffraction grating experiments with various configurations.
    Uses dependency injection to separate physics calculations from visualization.
    """
    
    def __init__(self, 
                grating_spacing: float = 1.67e-6,
                distance_to_screen: float = 1.0, 
                screen_width: float = 1.0,
                num_slits: int = 100):
        """
        Initialize the simulator with default parameters.
        
        Args:
            grating_spacing: Distance between slits in meters
            distance_to_screen: Distance from grating to screen in meters
            screen_width: Width of the screen in meters
            num_slits: Number of slits in the grating
        """
        self.grating_spacing = grating_spacing
        self.distance_to_screen = distance_to_screen
        self.screen_width = screen_width
        self.num_slits = num_slits
        
        # Create physics models
        self.physics = DiffractionPhysics(
            grating_spacing=grating_spacing,
            distance_to_screen=distance_to_screen,
            screen_width=screen_width
        )
        
        self.infinite_physics = InfiniteSlitDiffractionPhysics(
            grating_spacing=grating_spacing,
            distance_to_screen=distance_to_screen,
            screen_width=screen_width
        )
        
        # Create visualizer
        self.visualizer = DiffractionVisualizer()
    
    def update_parameters(self, 
                         grating_spacing: Optional[float] = None,
                         distance_to_screen: Optional[float] = None,
                         screen_width: Optional[float] = None,
                         num_slits: Optional[int] = None) -> None:
        """
        Update the simulator parameters.
        
        Args:
            grating_spacing: Distance between slits in meters
            distance_to_screen: Distance from grating to screen in meters
            screen_width: Width of the screen in meters
            num_slits: Number of slits in the grating
        """
        if grating_spacing is not None:
            self.grating_spacing = grating_spacing
            self.physics.grating_spacing = grating_spacing
            self.infinite_physics.grating_spacing = grating_spacing
            
        if distance_to_screen is not None:
            self.distance_to_screen = distance_to_screen
            self.physics.distance_to_screen = distance_to_screen
            self.infinite_physics.distance_to_screen = distance_to_screen
            
        if screen_width is not None:
            self.screen_width = screen_width
            self.physics.screen_width = screen_width
            self.infinite_physics.screen_width = screen_width
            
        if num_slits is not None:
            self.num_slits = num_slits
    
    def _get_visible_spectrum_data(self) -> Tuple[List[float], List[str], List[str]]:
        """Get wavelengths, colors and labels for the visible spectrum."""
        colors = list(WAVELENGTH_OPTIONS.keys())
        wavelengths = [WAVELENGTH_OPTIONS[color] for color in colors]
        plot_colors = [COLOR_MAP[color] for color in colors]
        labels = [f"{color} ({wavelength*1e9:.0f} nm)" for color, wavelength in zip(colors, wavelengths)]
        
        return wavelengths, plot_colors, labels
    
    def _get_color_for_wavelength(self, wavelength: float) -> str:
        """Determine plot color based on wavelength."""
        if wavelength < 450e-9:
            return 'violet'
        elif wavelength < 495e-9:
            return 'b'
        elif wavelength < 570e-9:
            return 'g'
        elif wavelength < 590e-9:
            return 'gold'
        elif wavelength < 620e-9:
            return 'darkorange'
        else:
            return 'r'
            
    def simulate_visible_spectrum(self, 
                                 show_maxima: bool = True,
                                 use_infinite_slits: bool = False) -> None:
        """
        Simulate a diffraction pattern with the full visible spectrum.
        
        Args:
            show_maxima: Whether to show diffraction maxima positions
            use_infinite_slits: Whether to use the infinite slits model
        """
        wavelengths, plot_colors, labels = self._get_visible_spectrum_data()
        
        if use_infinite_slits:
            self.visualizer.plot_infinite_slit_pattern(
                physics=self.infinite_physics,
                wavelengths=wavelengths,
                colors=plot_colors,
                labels=labels
            )
        else:
            self.visualizer.plot_diffraction_pattern(
                physics=self.physics,
                wavelengths=wavelengths,
                colors=plot_colors,
                labels=labels,
                num_slits=self.num_slits,
                show_maxima=show_maxima
            )
    
    def simulate_monochromatic(self,
                              wavelength: float = 532e-9,
                              show_maxima: bool = True,
                              use_infinite_slits: bool = False) -> None:
        """
        Simulate a diffraction pattern with monochromatic light.
        
        Args:
            wavelength: Light wavelength in meters
            show_maxima: Whether to show diffraction maxima positions
            use_infinite_slits: Whether to use the infinite slits model
        """
        label = f"{wavelength*1e9:.0f} nm"
        color = self._get_color_for_wavelength(wavelength)
        
        if use_infinite_slits:
            self.visualizer.plot_infinite_slit_pattern(
                physics=self.infinite_physics,
                wavelengths=[wavelength],
                colors=[color],
                labels=[label]
            )
        else:
            self.visualizer.plot_diffraction_pattern(
                physics=self.physics,
                wavelengths=[wavelength],
                colors=[color],
                labels=[label],
                num_slits=self.num_slits,
                show_maxima=show_maxima
            )
    
    def simulate_custom_spectrum(self,
                                intensity_function: Callable[[np.ndarray], np.ndarray],
                                wavelength_range: Tuple[float, float] = (400e-9, 700e-9),
                                num_wavelengths: int = 10,
                                use_infinite_slits: bool = False) -> None:
        """
        Simulate a diffraction pattern with a custom spectrum of light.
        
        Args:
            intensity_function: Function that takes wavelengths and returns relative intensities
            wavelength_range: Tuple of (min_wavelength, max_wavelength) in meters
            num_wavelengths: Number of wavelength points to simulate
            use_infinite_slits: Whether to use the infinite slits model
        """
        # Generate wavelengths and intensities
        wavelengths = np.linspace(wavelength_range[0], wavelength_range[1], num_wavelengths)
        intensities = intensity_function(wavelengths)
        
        # Normalize intensities
        if np.max(intensities) > 0:
            intensities = intensities / np.max(intensities)
        
        # Get colors and create weighted wavelengths for simulation
        colors = [self._get_color_for_wavelength(wl) for wl in wavelengths]
        labels = [f"{wl*1e9:.0f} nm" for wl in wavelengths]
        
        physics = self.infinite_physics if use_infinite_slits else self.physics
        self.visualizer.plot_custom_spectrum(
            physics=physics,
            wavelengths=wavelengths,
            intensities=intensities,
            colors=colors,
            labels=labels,
            num_slits=self.num_slits,
            use_infinite_slits=use_infinite_slits
        )
    
    def compare_slit_numbers(self,
                            wavelength: float = 532e-9,
                            slit_numbers: List[int] = [2, 5, 10, 50, 100, 500],
                            show_infinite: bool = True) -> None:
        """
        Compare diffraction patterns with different numbers of slits.
        
        Args:
            wavelength: Light wavelength in meters
            slit_numbers: List of numbers of slits to compare
            show_infinite: Whether to include the infinite slits model
        """
        color = self._get_color_for_wavelength(wavelength)
        
        self.visualizer.plot_slit_comparison(
            finite_physics=self.physics,
            infinite_physics=self.infinite_physics if show_infinite else None,
            wavelength=wavelength,
            color=color,
            slit_numbers=slit_numbers
        )


# Legacy function-based API for backward compatibility
def simulate_visible_spectrum(
    grating_spacing: float = 1.67e-6,
    distance_to_screen: float = 1.0,
    screen_width: float = 1.0,
    num_slits: int = 100,
    show_maxima: bool = True,
    use_infinite_slits: bool = False
) -> None:
    """Legacy function adapter for backward compatibility."""
    simulator = DiffractionSimulator(
        grating_spacing=grating_spacing,
        distance_to_screen=distance_to_screen,
        screen_width=screen_width,
        num_slits=num_slits
    )
    simulator.simulate_visible_spectrum(
        show_maxima=show_maxima,
        use_infinite_slits=use_infinite_slits
    )


def simulate_monochromatic(
    wavelength: float = 532e-9,
    grating_spacing: float = 1.67e-6,
    distance_to_screen: float = 1.0,
    screen_width: float = 1.0,
    num_slits: int = 100,
    show_maxima: bool = True,
    use_infinite_slits: bool = False
) -> None:
    """Legacy function adapter for backward compatibility."""
    simulator = DiffractionSimulator(
        grating_spacing=grating_spacing,
        distance_to_screen=distance_to_screen,
        screen_width=screen_width,
        num_slits=num_slits
    )
    simulator.simulate_monochromatic(
        wavelength=wavelength,
        show_maxima=show_maxima,
        use_infinite_slits=use_infinite_slits
    )


def simulate_custom_spectrum(
    intensity_function: Callable[[np.ndarray], np.ndarray],
    wavelength_range: Tuple[float, float] = (400e-9, 700e-9),
    num_wavelengths: int = 10,
    grating_spacing: float = 1.67e-6,
    distance_to_screen: float = 1.0,
    screen_width: float = 1.0,
    num_slits: int = 100,
    use_infinite_slits: bool = False
) -> None:
    """Legacy function adapter for backward compatibility."""
    simulator = DiffractionSimulator(
        grating_spacing=grating_spacing,
        distance_to_screen=distance_to_screen,
        screen_width=screen_width,
        num_slits=num_slits
    )
    simulator.simulate_custom_spectrum(
        intensity_function=intensity_function,
        wavelength_range=wavelength_range,
        num_wavelengths=num_wavelengths,
        use_infinite_slits=use_infinite_slits
    )


def simulate_comparing_slit_numbers(
    wavelength: float = 532e-9,
    grating_spacing: float = 1.67e-6,
    distance_to_screen: float = 1.0,
    screen_width: float = 1.0,
    slit_numbers: List[int] = [2, 5, 10, 50, 100, 500],
    show_infinite: bool = True
) -> None:
    """Legacy function adapter for backward compatibility."""
    simulator = DiffractionSimulator(
        grating_spacing=grating_spacing,
        distance_to_screen=distance_to_screen,
        screen_width=screen_width
    )
    simulator.compare_slit_numbers(
        wavelength=wavelength,
        slit_numbers=slit_numbers,
        show_infinite=show_infinite
    ) 