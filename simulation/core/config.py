#!/usr/bin/env python3
# Configuration file for the diffraction grating simulation

"""Configuration constants for the diffraction simulation."""

# Default simulation parameters
DEFAULT_GRATING_SPACING = 1.67e-6  # 600 lines/mm in meters
DEFAULT_NUM_SLITS = 5
DEFAULT_NUM_POINTS = 1000

# Wavelengths of visible light in meters
WAVELENGTH_OPTIONS = {"violet": 405e-9, "blue": 470e-9, "green": 532e-9, "yellow": 580e-9, "orange": 605e-9, "red": 650e-9}

# Default wavelengths to display
DEFAULT_WAVELENGTHS = ["red", "green", "blue"]

# Color mapping for wavelengths
COLOR_MAP = {"violet": "violet", "blue": "blue", "green": "green", "yellow": "gold", "orange": "darkorange", "red": "red"}

# RGB mapping for 2D visualizations
COLOR_RGB_MAP = {
    "violet": [0.55, 0, 0.55],
    "blue": [0, 0, 1],
    "green": [0, 0.8, 0],
    "gold": [1, 0.84, 0],
    "yellow": [1, 1, 0],
    "darkorange": [1, 0.55, 0],
    "orange": [1, 0.65, 0],
    "red": [1, 0, 0],
}
