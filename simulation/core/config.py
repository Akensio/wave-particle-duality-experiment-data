#!/usr/bin/env python3
# Configuration file for the diffraction grating simulation

# Default simulation parameters
DEFAULT_GRATING_SPACING = 1.67e-6  # 1.67 micrometers (600 lines/mm)
DEFAULT_DISTANCE_TO_SCREEN = 1.0  # 1 meter
DEFAULT_NUM_SLITS = 100  # Default number of slits
DEFAULT_SCREEN_WIDTH = 1.0  # 1 meter
DEFAULT_NUM_POINTS = 1000  # Number of points to calculate on the screen

# Predefined wavelengths for visible spectrum (in meters)
WAVELENGTH_OPTIONS = {
    'Red': 700e-9,
    'Orange': 620e-9,
    'Yellow': 580e-9,
    'Green': 530e-9,
    'Blue': 470e-9,
    'Violet': 420e-9
}

# Color mapping for each wavelength
COLOR_MAP = {
    'Red': 'r',
    'Orange': 'darkorange',
    'Yellow': 'gold',
    'Green': 'g',
    'Blue': 'b',
    'Violet': 'violet'
}

# RGB values for each color (for 2D visualization)
COLOR_RGB_MAP = {
    'r': [1, 0, 0], 
    'g': [0, 1, 0], 
    'b': [0, 0, 1], 
    'darkorange': [1, 0.55, 0], 
    'gold': [1, 0.84, 0], 
    'violet': [0.93, 0.51, 0.93]
}

# Default grating presets
GRATING_PRESETS = {
    '600 lines/mm\n(1.67µm)': {'spacing': 1.67e-6, 'slits': 100},
    '300 lines/mm\n(3.33µm)': {'spacing': 3.33e-6, 'slits': 75},
    '1200 lines/mm\n(0.83µm)': {'spacing': 0.83e-6, 'slits': 150}
}

# Default selected wavelengths
DEFAULT_WAVELENGTHS = ['Red', 'Green', 'Blue'] 