"""
Diffraction Grating Simulation Package

This package provides a modular framework for simulating light diffraction 
through a grating and visualizing the resulting patterns.
"""

from simulation.core.physics import DiffractionModel
from simulation.core.simulator import DiffractionSimulator
from simulation.visualization.plotter import DiffractionVisualizer
from simulation.ui.interactive import run_interactive_simulation
