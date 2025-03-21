# Diffraction Grating Physics: Mathematical Foundation

This document explains the key physical principles and equations used in the wave-particle duality diffraction simulation.

## Multiple Slit Diffraction

The simulation models diffraction patterns created when light passes through multiple slits (a diffraction grating). The physics is based on the wave nature of light, where waves passing through slits interfere constructively and destructively.

### Multiple Slit Interference

For a grating with N slits of spacing d, the intensity pattern on a far screen is given by:

$$I(\theta) = I_0 \left( \frac{\sin(N\delta/2)}{\sin(\delta/2)} \right)^2$$

Where:
- $I_0$ is the maximum intensity
- $\delta = \frac{2\pi d \sin\theta}{\lambda}$ is the phase difference between adjacent slits
- $\theta$ is the angle from the central axis
- $\lambda$ is the wavelength of light

In the code, this is implemented in `DiffractionPhysics.calculate_intensity_pattern()`:

```python
# Calculate the phase difference between adjacent slits
delta = (2 * np.pi / wavelength) * self.grating_spacing * np.sin(angles)

# Calculate the multiple slit interference factor
intensity_factor = np.sin(num_slits * delta / 2)**2 / np.sin(delta / 2)**2
```

### Maxima Positions (Diffraction Orders)

The positions of intensity maxima (bright spots) are determined by the grating equation:

$$d \sin\theta_m = m\lambda$$

Where:
- $m$ is the diffraction order (integer: 0, ±1, ±2, ...)
- $\theta_m$ is the angle of the m-th order maximum

In the code, this is implemented in `DiffractionPhysics.calculate_maxima_positions()`:

```python
# Check if the maximum is physically possible
if abs(m * wavelength / self.grating_spacing) < 1:
    angle = np.arcsin(m * wavelength / self.grating_spacing)
```

### Semi-Circular Screen at Infinite Distance

The simulation models light diffraction onto a semi-circular screen at an effectively infinite distance from the diffraction grating. This better represents the angular nature of diffraction:

- Angles range from -π/2 to +π/2 radians (-90° to +90°)
- The diffraction pattern is plotted directly against the angle θ

### Infinite Slit Approximation

For a perfect grating with infinite slits, the diffraction pattern consists of sharp peaks at the positions given by the grating equation. In the simulation, this is approximated by placing narrow Gaussian peaks at the positions determined by the grating equation:

$$I(\theta) = \sum_{m} \exp\left(-\frac{1}{2}\left(\frac{\theta - \theta_m}{w}\right)^2\right)$$

Where:
- $w$ is a small width parameter to create sharp peaks
- $\theta_m$ are the angles satisfying the grating equation

In the code, this is implemented in the `_calculate_large_slit_pattern()` method:

```python
# Add peaks at each maximum location
for m, angle in maxima:
    # Create a sharp peak around that position
    intensity += np.exp(-0.5 * ((angles - angle) / peak_width) ** 2)
```

### Spectrum Visualization

For visualization with multiple wavelengths, the simulation calculates intensity patterns for each individual wavelength and combines them to form the final pattern:

$$I_{total}(\theta) = \sum_{\lambda} I_{\lambda}(\theta)$$

Where $I_{\lambda}(\theta)$ is the intensity pattern for wavelength $\lambda$.

#### Color Mapping

The simulation maps wavelengths to colors for visualization using the following approximations:

- Violet: 380-450 nm
- Blue: 450-495 nm
- Green: 495-570 nm
- Yellow: 570-590 nm
- Orange: 590-620 nm
- Red: 620-750 nm

These color mappings are used to create the 2D visualization of the diffraction pattern as it would appear on a screen.

## Blackbody Radiation through Diffraction Grating

This section explores how diffraction gratings can be used to analyze blackbody radiation spectra, based on the principles outlined in our simulation.

### Diffraction Grating for Spectroscopy

When using a diffraction grating to analyze blackbody radiation instead of a prism, we can leverage the grating equation to precisely separate different wavelengths. The intensity at a specific angle $\theta$ on a screen placed at a large distance from the grating is determined by the wavelength satisfying the maximum condition for a given diffraction order $m$:

$$\lambda_m = \frac{d \sin\theta}{m}$$

Where:
- $d$ is the grating spacing
- $\theta$ is the angle from the central axis
- $m$ is the diffraction order (integer: 1, 2, 3, ...)

### Blackbody Radiation Intensity

For blackbody radiation at temperature $T$, Planck's law gives the spectral intensity distribution:

$$I(\lambda, T) = \frac{2\pi c^2 h}{\lambda^5}\frac{1}{e^{\frac{hc}{\lambda k T}} - 1}$$

Where:
- $h$ is Planck's constant
- $c$ is the speed of light
- $k$ is Boltzmann's constant
- $T$ is the temperature in Kelvin

This is implemented in the `ContinuousSpectrumPhysics` class:

```python
def planck_spectrum(self, wavelengths: NDArray[np.float64], temperature: float) -> NDArray[np.float64]:
    # Physical constants
    h = constants.Planck  # Planck's constant (J·s)
    c = constants.speed_of_light  # Speed of light (m/s)
    k_B = constants.Boltzmann  # Boltzmann constant (J/K)
    
    # Calculate Planck's law formula
    numerator = 2.0 * h * c**2
    denominator = wavelengths**5 * (np.exp((h*c)/(wavelengths*k_B*temperature)) - 1.0)
    intensity = numerator / denominator

    return intensity
```

### Combined Intensity Pattern for Diffraction Order m

When blackbody radiation passes through a diffraction grating, for a given diffraction order $m$, the intensity at angle $\theta$ is:

$$I_m(\theta) = \frac{2\pi c^2 h m^5}{d^5 \sin^5\theta}\frac{1}{e^{\frac{h c m}{d \sin\theta k T}} - 1}$$

This equation combines the diffraction grating equation with Planck's law, by substituting $\lambda_m = \frac{d \sin\theta}{m}$ into Planck's law.

In the simulation, this is implemented in the `calculate_order_m_pattern` method:

```python
def calculate_order_m_pattern(self, angles: NDArray[np.float64], temperature: float, order_m: int) -> NDArray[np.float64]:
    """Calculate diffraction pattern using equation 3 from eqs2.md directly."""
    d = self.grating_spacing
    wavelengths = np.abs(np.sin(angles)) * d / order_m
    intensity = self.planck_spectrum(wavelengths, temperature)
    return intensity
```

### Total Diffraction Pattern

The full diffraction pattern is the sum of the patterns from each diffraction order:

$$I(\theta) = \sum_{m=1}^{\infty} I_m(\theta)$$

In practice, our simulation uses a finite number of orders since higher orders contribute progressively less to the total pattern.

### Visualization in the Simulation

The `ContinuousSpectrumSimulation` class in `continuous_interactive.py` visualizes how blackbody radiation at different temperatures creates diffraction patterns when passed through a diffraction grating:

1. Users can adjust the temperature to see how the blackbody spectrum changes
2. The diffraction pattern shows how different wavelengths in the blackbody spectrum separate at different angles
3. The simulator shows both individual diffraction orders and their sum

This provides insight into how diffraction gratings can be used as spectroscopic tools for analyzing continuous spectra like blackbody radiation. 

## Wave-Particle Duality

This simulation demonstrates the wave nature of light. According to wave-particle duality in quantum mechanics, light exhibits both wave-like and particle-like properties. In diffraction experiments:

- The wave nature of light is demonstrated by interference patterns
- The particle nature appears when detecting individual photons

The diffraction patterns observed in this simulation are a direct manifestation of the wave nature of light, showing how waves from different slits interfere constructively and destructively based on their path differences.