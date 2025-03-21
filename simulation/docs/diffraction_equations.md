# Diffraction Grating Physics: Mathematical Foundation

This document explains the key physical principles and equations used in the wave-particle duality diffraction simulation.

## Multiple Slit Diffraction

The simulation models diffraction patterns created when light passes through multiple slits (a diffraction grating). The physics is based on the wave nature of light, where waves passing through slits interfere constructively and destructively.

### Multiple Slit Interference

For a grating with N slits of spacing d, the intensity pattern on a screen at distance L is given by:

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

The positions of intensity maxima (bright spots) on the screen are determined by the grating equation:

$$d \sin\theta_m = m\lambda$$

Where:
- $m$ is the diffraction order (integer: 0, ±1, ±2, ...)
- $\theta_m$ is the angle of the m-th order maximum

In the code, this is implemented in `DiffractionPhysics.calculate_maxima_positions()`:

```python
# Check if the maximum is physically possible
if abs(m * wavelength / self.grating_spacing) < 1:
    angle = np.arcsin(m * wavelength / self.grating_spacing)
    pos = self.distance_to_screen * np.tan(angle)
```

### Infinite Slit Approximation

For a perfect grating with infinite slits, the diffraction pattern consists of sharp peaks at the positions given by the grating equation. In the simulation, this is approximated by placing narrow Gaussian peaks at the positions determined by the grating equation:

$$I(\theta) = \sum_{m} \exp\left(-\frac{1}{2}\left(\frac{\theta - \theta_m}{w}\right)^2\right)$$

Where:
- $w$ is a small width parameter to create sharp peaks
- $\theta_m$ are the angles satisfying the grating equation

In the code, this is implemented in `InfiniteSlitDiffractionPhysics.calculate_intensity_pattern()`:

```python
# Add peaks at each maximum location
for m, pos in maxima:
    # Create a sharp peak around that position
    width = 0.01  # This is in meters and defines the sharpness of the peak
    intensity += np.exp(-0.5 * ((screen_positions - pos) / width) ** 2)
```

## Spectrum Visualization

For visualization with multiple wavelengths, the simulation calculates intensity patterns for each individual wavelength and combines them to form the final pattern:

$$I_{total}(\theta) = \sum_{\lambda} I_{\lambda}(\theta)$$

Where $I_{\lambda}(\theta)$ is the intensity pattern for wavelength $\lambda$.

## Color Mapping

The simulation maps wavelengths to colors for visualization using the following approximations:

- Violet: 380-450 nm
- Blue: 450-495 nm
- Green: 495-570 nm
- Yellow: 570-590 nm
- Orange: 590-620 nm
- Red: 620-750 nm

These color mappings are used to create the 2D visualization of the diffraction pattern as it would appear on a screen.

## Wave-Particle Duality

This simulation demonstrates the wave nature of light. According to wave-particle duality in quantum mechanics, light exhibits both wave-like and particle-like properties. In diffraction experiments:

- The wave nature of light is demonstrated by interference patterns
- The particle nature appears when detecting individual photons

The diffraction patterns observed in this simulation are a direct manifestation of the wave nature of light, showing how waves from different slits interfere constructively and destructively based on their path differences. 