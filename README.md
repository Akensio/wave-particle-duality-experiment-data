# Wave-Particle Duality Experiment

This repository contains tools for exploring and analyzing wave-particle duality phenomena through both simulation and experimental data analysis.

## Repository Structure

- **experiment/** - Analysis scripts and data for physical experiments
- **simulation/** - Diffraction grating simulation tools
- **common/** - Shared utilities used by both components

## Requirements

This project uses Poetry for dependency management. All necessary dependencies are already defined in the `pyproject.toml` file.

Required packages:
- numpy
- matplotlib
- scipy
- pandas

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