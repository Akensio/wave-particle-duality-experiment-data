# Wave-Particle Duality Experiments

This directory contains analysis scripts and data for physical experiments related to wave-particle duality phenomena.

## Experiments

### Black Body Radiation

The black body radiation analysis explores the relationship between temperature and wavelength distribution of thermal radiation.

- **Script**: `analyze_black_body_radiation.py`
- **Data**: `data/black_body_radiation_spectrum.xlsx`
- **Documentation**: `equations_black_body_radiation.md`

Key features:
- Calculates and plots the relationship between temperature and peak wavelength
- Validates Wien's displacement law (λ_max * T = constant)
- Determines the refraction index for different wavelengths
- Analyzes resistivity and temperature measurements

### Electron Diffraction

The electron diffraction analysis examines the wave properties of electrons through diffraction patterns.

- **Script**: `analyze_dispersion.py`
- **Data**: `data/dispersion.xlsx`
- **Documentation**: `equations_dispersion.md`

Key features:
- Calculates de Broglie wavelengths for electrons at different accelerating voltages
- Measures diffraction ring radii and analyzes their relationship with wavelength
- Determines interplanar distances in the crystal

## Usage

Run the analysis scripts from the root directory:

```bash
python -m experiment.analyze_black_body_radiation
python -m experiment.analyze_dispersion
```

Outputs and generated plots will be saved to the `out/` directory.

## Data Structure

The experimental data is stored in Excel files with multiple sheets:

### Black Body Radiation:
- **rotation_ratio**: Measurements for determining the small circle angle per large circle revolution
- **starting_angle**: Initial angle measurements
- **lightbulb_with_wires**: Voltage and current measurements for resistance calculations
- **max_intensity_wavelength**: Data for analyzing the wavelength of maximum intensity at different temperatures

### Dispersion:
- **dispersion**: Contains voltage values and diameter measurements for electron diffraction rings

## Physical Principles

The experiments verify two key quantum phenomena:
1. Black body radiation and Wien's displacement law: λ_max * T = constant
2. De Broglie relation for electron waves: λ = h/√(2mEk), demonstrating the wave-particle duality of matter

### Documentation

Detailed theoretical background and equations are available in the following files:

- [Black Body Radiation Equations](equations_black_body_radiation.md) - Contains:
  - Wien's displacement law derivation
  - Resistivity and temperature relationship for tungsten
  - Planck's radiation law formulation
  - Detailed error analysis methodology

- [Electron Diffraction Equations](equations_dispersion.md) - Contains:
  - De Broglie wavelength calculation
  - Diffraction pattern angle relationships
  - Bragg's law for crystalline materials
  - Formulas for interplanar distance determination
