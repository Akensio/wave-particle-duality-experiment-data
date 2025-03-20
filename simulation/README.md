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
- **Gaussian Wave Packet Support**: More realistic simulation using wave packets instead of monochromatic light
  - Spectral analysis of wave packets
  - Time evolution visualization
  - Dispersion effects demonstration

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

### Running the Wave Packet Simulation

To run the Gaussian wave packet simulation:

```bash
python simulation/run_wave_packet.py
```

Command line options:
- `--center-wavelength`: Center wavelength of the wave packet (color name like "red", "green", etc.)
- `--wavelength-width`: Width of the Gaussian distribution in nanometers
- `--grating-spacing`: Spacing between slits in meters
- `--distance`: Distance from grating to screen in meters
- `--screen-width`: Width of the screen in meters
- `--num-slits`: Number of slits in the diffraction grating
- `--analysis`: Show full analysis with multiple plots (spectrum, diffraction pattern, time evolution)

Example:
```bash
python simulation/run_wave_packet.py --center-wavelength green --wavelength-width 10 --num-slits 5 --analysis
```

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
    - `base.py`: Base diffraction physics model
    - `gaussian_packet.py`: Gaussian wave packet model
  - `simulator.py`: Pre-defined simulation scenarios
  - `config.py`: Constants and default values
- `visualization/`: Visualization components
  - `plotter.py`: Functions for plotting diffraction patterns
  - `wall_pattern.py`: 2D wall pattern visualization
  - `wave_packet_plotter.py`: Wave packet visualization tools
- `ui/`: User interface components
  - `interactive.py`: Interactive simulation with sliders
- `docs/`: Documentation
  - `diffraction_equations.md`: Detailed physics equations

## Physics Background

### Single Wavelength Diffraction

The simulation is based on the principles of Fraunhofer diffraction. For a diffraction grating with N slits and spacing d, the intensity at angle θ is given by:

$I(\theta) = I_0 \left( \frac{\sin(N \delta /2)}{\sin(\delta/2)} \right)^2$

where $\delta = \frac{2\pi d \sin\theta}{\lambda}$ is the phase difference between adjacent slits, and λ is the wavelength of light.

The positions of intensity maxima are given by the grating equation:

$d \sin\theta = m\lambda$

where m is an integer representing the order of diffraction.

### Gaussian Wave Packet Diffraction

The simulation also supports Gaussian wave packets, which represent a more realistic model of light pulses. A Gaussian wave packet has a distribution of wavelengths centered around a peak:

$\Psi(x,t) = A \cdot e^{-\frac{(x-vt)^2}{2\sigma^2}} \cdot e^{i(k_0x - \omega_0t)}$

The diffraction pattern for a wave packet is calculated by integrating over all wavelengths in the packet:

$I(\theta) = \int_0^{\infty} S(\lambda) \cdot I_{\lambda}(\theta) \, d\lambda$

Where $S(\lambda)$ is the spectral intensity of the wave packet.

## License

This project is open source and available under the MIT License.