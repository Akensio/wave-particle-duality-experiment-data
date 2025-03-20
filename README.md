# Wave-Particle Duality Experiment Data

A repository containing simulation and experiment data related to wave-particle duality phenomena.

## Overview

- **simulation/** - Diffraction grating simulation tools with various modes
- **experiment/** - Analysis scripts and data for physical experiments
- **common/** - Shared utilities

## Getting Started

1. Clone the repository
2. Install Poetry (if not already installed):
   ```bash
   pip install poetry
   ```
3. Install dependencies:
   ```bash
   poetry install
   poetry shell
   ```

## Running the Simulation

Always run from the ROOT directory using either:

```bash
python -m simulation.run [options]
```

or 

```bash
python simulation/run.py [options]
```

For detailed usage and parameters, see [simulation README](simulation/README.md).

## Experiment Data

For information on experiment data and analysis tools, see [experiment README](experiment/README.md).

## License

This project is licensed under the terms included in the LICENSE file.