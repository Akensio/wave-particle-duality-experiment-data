# Diffraction Grating Simulation

A modular simulation of light passing through a diffraction grating and displaying the resulting pattern.

## Overview

This simulation demonstrates the principles of diffraction and interference of light passing through a diffraction grating. It provides various visualizations and interactive tools to explore how different parameters affect the diffraction pattern.

## Features

- **Modular Architecture**: Organized into core physics, visualization, and UI components
- **Interactive Mode**: Adjust parameters in real-time using sliders
- **Multiple Simulation Types**:
  - Interactive simulation (finite or infinite slits)
  - Visible spectrum simulation
  - Monochromatic light simulation
  - Custom spectrum simulation
  - Comparison of different slit numbers
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
python simulation/run.py [options]
```

### Interactive Simulation

Run the interactive simulation with:

```bash
python simulation/run.py --mode interactive
```

For a perfect diffraction grating (infinite slits):

```bash
python simulation/run.py --mode interactive-infinite
```

### Fixed Simulations

#### Visible Spectrum

```bash
python simulation/run.py --mode spectrum
```

For a perfect grating with infinite slits:

```bash
python simulation/run.py --mode spectrum-infinite
```

#### Monochromatic Light

```bash
python simulation/run.py --mode monochromatic --wavelength 532nm
```

For a perfect grating with infinite slits:

```bash
python simulation/run.py --mode monochromatic-infinite --wavelength 632.8nm
```

#### Custom Spectrum

```bash
python simulation/run.py --mode custom --spectrum-type gaussian --center 550nm --spectrum-width 30nm
```

For a double-peak spectrum:

```bash
python simulation/run.py --mode custom --spectrum-type double-peak --center 450nm --center2 650nm --spectrum-width 20nm
```

#### Compare Different Slit Numbers

```bash
python simulation/run.py --mode compare --wavelength 632.8nm
```

### Common Parameters

These parameters can be used with any simulation mode:

- `--spacing`: Grating spacing in meters (default: 1.67e-6 for 600 lines/mm)
- `--distance`: Distance to screen in meters (default: 1.0)
- `--width`: Screen width in meters (default: 1.0)
- `--slits`: Number of slits (default: 100)

### Custom Spectrum Parameters

- `--spectrum-type`: Type of spectrum ('gaussian' or 'double-peak')
- `--center`: Center wavelength for the spectrum
- `--spectrum-width`: Width (standard deviation) of the spectrum
- `--center2`: Second peak wavelength (only for double-peak spectrum)

### Wavelength Specification

Wavelengths can be specified in different formats:

- With units: `500nm`, `0.5um`, `0.5µm`
- Direct value in meters: `500e-9`

## Code Structure

- `core/`: Core physics calculations
  - `physics.py`: Diffraction intensity calculations
  - `simulator.py`: Pre-defined simulation scenarios
  - `config.py`: Constants and default values
- `visualization/`: Visualization components
  - `plotter.py`: Functions for plotting diffraction patterns
- `ui/`: User interface components
  - `interactive.py`: Interactive simulation with sliders

## Physics Background

The simulation is based on the principles of Fraunhofer diffraction. For a diffraction grating with N slits and spacing d, the intensity at angle θ is given by:

$I(\theta) = I_0 \left( \frac{\sin(N \delta /2)}{\sin(\delta/2)} \right)^2$

where $\delta = \frac{2\pi d \sin\theta}{\lambda}$ is the phase difference between adjacent slits, and λ is the wavelength of light.

The positions of intensity maxima are given by the grating equation:

$d \sin\theta = m\lambda$

where m is an integer representing the order of diffraction.

## Examples

### Comparing Finite vs Infinite Slits

- Finite slits (e.g., 100): 
  ```bash
  python simulation/run.py --mode monochromatic --slits 100
  ```
- Infinite slits: 
  ```bash
  python simulation/run.py --mode monochromatic-infinite
  ```

### Different Grating Spacings

- 600 lines/mm: 
  ```bash
  python simulation/run.py --mode spectrum --spacing 1.67e-6
  ```
- 1200 lines/mm: 
  ```bash
  python simulation/run.py --mode spectrum --spacing 0.83e-6
  ```

## Infinite Slits Approximation

The "infinite slits" model creates sharp peaks at locations determined by the grating equation (d·sin θ = m·λ). This represents a perfect diffraction grating where the number of slits approaches infinity, resulting in perfect constructive interference only at the exact angles that satisfy the equation.

## License

This project is open source and available under the MIT License.