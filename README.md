# Wave-Particle Duality Experiment Data

This repository contains code and data for experiments related to wave-particle duality.

## Diffraction Grating Simulation

A modular simulation of light passing through a diffraction grating, demonstrating wave properties of light.

### Key Features
- Interactive diffraction pattern visualization
- Support for finite and infinite slits ("perfect" diffraction grating)
- Customizable wavelengths, grating spacing, and other parameters
- Visualization of how different spectra appear after diffraction

### Running the Simulation

**IMPORTANT: Always run the simulation from the ROOT directory of the project:**

```bash
python simulation/run_simulation.py [options]
```

#### Interactive Mode with Infinite Slits

To run the perfect diffraction grating simulation with infinite slits:

```bash
python simulation/run_simulation.py --mode interactive-infinite
```

This will allow you to:
- Control the grating spacing
- Select different wavelengths
- See the resulting diffraction pattern in real-time
- Visualize how the pattern would appear on a wall/screen

#### Other Simulation Modes

There are several other simulation modes available:

```bash
# Standard interactive mode (finite slits)
python simulation/run_simulation.py --mode interactive

# Visible spectrum simulation
python simulation/run_simulation.py --mode spectrum
python simulation/run_simulation.py --mode spectrum-infinite  # With infinite slits

# Single wavelength simulation
python simulation/run_simulation.py --mode monochromatic --wavelength 532nm
python simulation/run_simulation.py --mode monochromatic-infinite --wavelength 532nm  # With infinite slits

# Custom spectrum simulation
python simulation/run_simulation.py --mode custom --spectrum-type gaussian --center 550nm --spectrum-width 30nm
python simulation/run_simulation.py --mode custom-infinite --spectrum-type gaussian --center 550nm --spectrum-width 30nm  # With infinite slits

# Compare different numbers of slits
python simulation/run_simulation.py --mode compare --wavelength 632.8nm
```

For detailed documentation, see the [simulation README](simulation/README.md).

## Requirements

The simulation requires the following Python packages:
- numpy
- matplotlib
- scipy

Install them with:
```bash
pip install numpy matplotlib scipy
```

## Repository Structure

- **experiment/** - Analysis scripts and data for physical experiments
- **simulation/** - Diffraction grating simulation tools
- **common/** - Shared utilities used by both components

## Getting Started

1. Clone the repository
2. Install dependencies using Poetry:
   ```bash
   poetry install
   ```
3. Activate the virtual environment:
   ```bash
   poetry shell
   ```

## Components

### Simulation

The simulation component provides tools for visualizing diffraction patterns created by light passing through a diffraction grating. See the [simulation README](simulation/README.md) for more details.

### Experiment

The experiment component includes analysis scripts for processing experimental data related to black body radiation and electron diffraction. See the [experiment README](experiment/README.md) for more details.

## License

This project is licensed under the terms included in the LICENSE file.