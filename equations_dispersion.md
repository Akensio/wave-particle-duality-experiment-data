# Dispersion Analysis Equations

This document details the equations used in the dispersion analysis, including error propagation formulas.

## Wavelength Calculations

The de Broglie wavelength of electrons is calculated as:

$$\lambda = \frac{h}{\sqrt{2 m_e e}} \cdot \frac{1}{\sqrt{V}}$$

Where:
- $h$ is Planck's constant (6.62607015 × 10⁻³⁴ J·s)
- $m_e$ is the electron mass (9.1093837015 × 10⁻³¹ kg)
- $e$ is the elementary charge (1.602176634 × 10⁻¹⁹ C)
- $V$ is the accelerating voltage (in volts)

For convenience, we define $f(V) = \frac{1}{\sqrt{V}}$.

Error propagation for $f(V)$:

$$\Delta f(V) = \left|\frac{df}{dV}\right| \cdot \Delta V = \left|-\frac{1}{2} \cdot V^{-3/2}\right| \cdot \Delta V = \frac{1}{2} \cdot V^{-3/2} \cdot \Delta V = \frac{0.5}{\sqrt{V}} \cdot \frac{\Delta V}{V} = \frac{f(V)}{2} \cdot \frac{\Delta V}{V}$$

Therefore, the error propagation for wavelength can be calculated as:

$$\Delta \lambda = \frac{\lambda}{2} \cdot \frac{\Delta V}{V}$$

## Radius Measurements

For each diffraction ring, two measurements are taken of both the inner and outer edges. The average diameter and error are calculated as follows:

For the first measurement of a ring:

$$\text{diameter\_meas1} = \frac{\text{diameter\_internal} + \text{diameter\_external}}{2}$$

$$\text{diameter\_meas1\_error} = \frac{|\text{diameter\_external} - \text{diameter\_internal}|}{2}$$

For the second measurement of the same ring:

$$\text{diameter\_meas2} = \frac{\text{diameter\_internal\_meas2} + \text{diameter\_external\_meas2}}{2}$$

$$\text{diameter\_meas2\_error} = \frac{|\text{diameter\_external\_meas2} - \text{diameter\_internal\_meas2}|}{2}$$

The final diameter and its error are calculated by combining the two measurements:

$$\text{diameter} = \frac{\text{diameter\_meas1} + \text{diameter\_meas2}}{2}$$

$$\text{diameter\_error} = \frac{\sqrt{(\text{diameter\_meas1\_error})^2 + (\text{diameter\_meas2\_error})^2}}{2}$$

The radius and its error are simply half of the diameter and its error:

$$r = \frac{\text{diameter}}{2}$$

$$\Delta r = \frac{\text{diameter\_error}}{2}$$

## Dispersion Analysis

The relationship between the radius of a diffraction ring and the electron wavelength is linear:

$$r = \text{slope} \cdot \lambda + \text{intercept}$$

The slope from this linear fit is related to the interplanar distance $d$ in the crystal:

$$d = \frac{2L}{\text{slope}}$$

Where $L$ is the radius of the spherical bulb.

Error propagation for interplanar distance:

$$\Delta d = d \cdot \sqrt{\left(\frac{\Delta \text{slope}}{\text{slope}}\right)^2 + \left(\frac{\Delta L}{L}\right)^2}$$

Where:
- $\Delta \text{slope}$ is the error in the slope from the linear fit
- $\Delta L$ is the error in the radius of the bulb

For the conversion to picometers:

$$d_{\text{pm}} = d_{\text{mm}} \cdot 10^9$$

$$\Delta d_{\text{pm}} = \Delta d_{\text{mm}} \cdot 10^9$$ 