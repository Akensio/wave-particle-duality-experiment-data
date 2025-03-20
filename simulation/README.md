# Diffraction Grating Simulation

A simple simulation of light passing through a diffraction grating and displaying the resulting pattern.

## Overview

This simulation demonstrates the principles of diffraction and interference of light passing through a diffraction grating. It provides an interactive visualization to explore how different parameters affect the diffraction pattern.

## Features

- **Modular Architecture**: Organized into core physics, visualization, and UI components
- **Interactive Interface**: Adjust parameters in real-time using sliders
- **Visualization Options**:
  - 1D intensity plots
  - 2D wall pattern visualization
  - Diffraction maxima indicators
  - Combined intensity for multiple wavelengths

## Installation

No special installation required beyond standard Python scientific packages:

```bash
pip install numpy matplotlib scipy
```

## Usage

### Running the Simulation

IMPORTANT: Always run the simulation from the ROOT directory of the project, not from inside the simulation folder!

```bash
python simulation/run.py
```

This will launch the interactive simulation where you can adjust parameters and see the diffraction pattern in real-time.

### Interactive Controls

The simulation provides interactive controls:
- Adjust grating spacing using the slider
- Change the number of slits
- Change the distance to the screen
- Adjust screen width
- Select/deselect different wavelengths of light
- Choose from predefined grating presets

## Code Structure

- `core/`: Core physics calculations
  - `physics/`: Diffraction intensity calculations
  - `simulator.py`: Pre-defined simulation scenarios
  - `config.py`: Constants and default values
- `visualization/`: Visualization components
  - `plotter.py`: Functions for plotting diffraction patterns
  - `wall_pattern.py`: 2D wall pattern visualization
- `ui/`: User interface components
  - `interactive.py`: Interactive simulation with sliders

## Physics Background

The simulation is based on the principles of Fraunhofer diffraction. For a diffraction grating with N slits and spacing d, the intensity at angle θ is given by:

$I(\theta) = I_0 \left( \frac{\sin(N \delta /2)}{\sin(\delta/2)} \right)^2$

where $\delta = \frac{2\pi d \sin\theta}{\lambda}$ is the phase difference between adjacent slits, and λ is the wavelength of light.

The positions of intensity maxima are given by the grating equation:

$d \sin\theta = m\lambda$

where m is an integer representing the order of diffraction.

## License

This project is open source and available under the MIT License.