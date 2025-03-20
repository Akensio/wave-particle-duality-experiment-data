import numpy as np
from typing import List, Tuple, Optional, Callable

from simulation.core.physics import DiffractionPhysics
from simulation.visualization import DiffractionVisualizer
from simulation.core.config import WAVELENGTH_OPTIONS, COLOR_MAP


class DiffractionSimulator:
    """Core simulator for diffraction grating experiments."""
    
    def __init__(
        self,
        grating_spacing: float = 1.67e-6,
        distance_to_screen: float = 1.0, 
        screen_width: float = 1.0,
        num_slits: int = 100
    ):
        """Initialize the simulator with given parameters."""
        self.grating_spacing = grating_spacing
        self.distance_to_screen = distance_to_screen
        self.screen_width = screen_width
        self.num_slits = num_slits
        
        # Create models
        self.physics = DiffractionPhysics(
            grating_spacing, distance_to_screen, screen_width
        )
        
        # Create visualizers
        self.visualizer = DiffractionVisualizer()
        
    def update_parameters(
        self, 
        grating_spacing: Optional[float] = None,
        distance_to_screen: Optional[float] = None,
        screen_width: Optional[float] = None,
        num_slits: Optional[int] = None
    ) -> None:
        """Update simulator parameters."""
        if grating_spacing is not None:
            self.grating_spacing = grating_spacing
            self.physics.grating_spacing = grating_spacing
            
        if distance_to_screen is not None:
            self.distance_to_screen = distance_to_screen
            self.physics.distance_to_screen = distance_to_screen
            
        if screen_width is not None:
            self.screen_width = screen_width
            self.physics.screen_width = screen_width
            
        if num_slits is not None:
            self.num_slits = num_slits
    
    def _get_color_for_wavelength(self, wavelength: float) -> str:
        """Determine color for a wavelength."""
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
    
    def simulate_visible_spectrum(
        self, 
        show_maxima: bool = True
    ) -> None:
        """Simulate the visible light spectrum diffraction pattern."""
        colors = list(WAVELENGTH_OPTIONS.keys())
        wavelengths = [WAVELENGTH_OPTIONS[color] for color in colors]
        plot_colors = [COLOR_MAP[color] for color in colors]
        labels = [f"{color} ({wl*1e9:.0f} nm)" for color, wl in zip(colors, wavelengths)]
        
        self.visualizer.plot_diffraction_pattern(
            physics=self.physics,
            wavelengths=wavelengths,
            colors=plot_colors,
            labels=labels,
            num_slits=self.num_slits,
            show_maxima=show_maxima
        )
    
    def simulate_monochromatic(
        self,
        wavelength: float = 532e-9,
        show_maxima: bool = True
    ) -> None:
        """Simulate a monochromatic light diffraction pattern."""
        label = f"{wavelength*1e9:.0f} nm"
        color = self._get_color_for_wavelength(wavelength)
        
        self.visualizer.plot_diffraction_pattern(
            physics=self.physics,
            wavelengths=[wavelength],
            colors=[color],
            labels=[label],
            num_slits=self.num_slits,
            show_maxima=show_maxima
        ) 