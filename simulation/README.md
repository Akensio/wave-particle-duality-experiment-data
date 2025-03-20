# Diffraction Grating Simulation

This Python program simulates the diffraction pattern created when light passes through a diffraction grating, displaying how different wavelengths of light create distinct patterns on a screen.

## Features

- Calculates and displays the diffraction pattern from multiple slits
- Supports overlapping multiple wavelengths on the same graph for comparison
- Visualizes diffraction patterns as they would appear on a wall/screen
- Shows diffraction maxima positions based on the grating equation
- Includes 2D visualization of how the pattern would look on a screen

## Requirements

This project uses Poetry for dependency management. All necessary dependencies are already defined in the `pyproject.toml` file.

Required packages:
- numpy
- matplotlib
- scipy

## Usage

1. Run the simulation:

```bash
python diffraction_grating_simulation.py
```

2. The program will generate three plots:
   - A comparison of diffraction patterns for the full visible spectrum
   - A comparison of RGB wavelengths
   - A detailed view of a single wavelength

## Physics Background

The simulation is based on the diffraction grating equation:

d sin θ = m λ

Where:
- d = spacing between slits
- θ = angle of diffraction
- m = order of diffraction (0, ±1, ±2, ...)
- λ = wavelength of light

For multiple slits, the intensity pattern follows:

I = I₀ (sin(Nδ/2) / sin(δ/2))²

Where:
- I₀ = maximum intensity
- N = number of slits
- δ = phase difference between adjacent slits

## Customization

You can modify the parameters in the code to experiment with different:
- Wavelengths (by changing the `wavelengths` lists)
- Grating spacings (by changing `grating_spacing`)
- Number of slits (by changing `num_slits`)
- Screen distances (by changing `distance_to_screen`)
- Screen sizes (by changing `screen_width`)

To plot your own custom wavelength combinations, modify the main function or create your own by calling:

```python
plot_diffraction_pattern(
    wavelengths=[your_wavelengths_here],
    colors=[your_colors_here],
    labels=[your_labels_here],
    grating_spacing=your_spacing,
    num_slits=your_number_of_slits
)
```