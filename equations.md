# Black Body Radiation Analysis Equations

This document details the equations used in the black body radiation spectrum analysis, including error propagation formulas.

## Rotation Ratio Analysis

The rotation ratio is calculated as:

$$\text{ANGLE\_RATIO} = \frac{\text{angle\_change\_per\_large\_circle\_full\_rotation\_in\_rad}}{2\pi}$$

Error propagation:

$$\text{ANGLE\_RATIO\_ERROR} = \frac{\text{angle\_change\_per\_large\_circle\_full\_rotation\_error\_in\_rad}}{2\pi}$$

where the error in the angle change is calculated as the standard error of the mean (SEM):

$$\text{angle\_change\_per\_large\_circle\_full\_rotation\_error\_in\_rad} = \frac{\sigma}{\sqrt{n}}$$

## Resistance Measurements

The resistance is calculated using Ohm's law:

$$R = \frac{V}{I}$$

Error propagation for resistance:

$$\left(\frac{\Delta R}{R}\right)^2 = \left(\frac{\Delta V}{V}\right)^2 + \left(\frac{\Delta I}{I}\right)^2$$

$$\Delta R = R \cdot \sqrt{\left(\frac{\Delta V}{V}\right)^2 + \left(\frac{\Delta I}{I}\right)^2}$$

The resistance of the wires is calculated as:

$$R_{\text{WIRES}} = R_{\text{total}} - R_{\text{LIGHTBULB}}$$

Error propagation for the wire resistance:

$$\Delta R_{\text{WIRES}} = \sqrt{(\Delta R_{\text{total}})^2 + (\Delta R_{\text{LIGHTBULB}})^2}$$

## Resistivity Calculation

The resistivity is calculated using:

$$\rho = \rho_0 \cdot \frac{R - R_w}{R_0}$$

where:
- $\rho_0$ is the resistivity at room temperature (5.65 × 10⁻⁸ Ω·m)
- $R$ is the measured resistance
- $R_w$ is the resistance of the wires
- $R_0$ is the resistance of the lightbulb at room temperature

Error propagation for resistivity:

$$\frac{\Delta \rho}{\rho} = \sqrt{\left(\frac{\Delta R}{R-R_w}\right)^2 + \left(\frac{\Delta R_w}{R-R_w}\right)^2 + \left(\frac{\Delta R_0}{R_0}\right)^2}$$

$$\Delta \rho = \rho \cdot \sqrt{\left(\frac{\Delta R}{R-R_w}\right)^2 + \left(\frac{\Delta R_w}{R-R_w}\right)^2 + \left(\frac{\Delta R_0}{R_0}\right)^2}$$

## Temperature Calculation

The temperature is calculated using the formula:

$$T = 103 + 38.1 \cdot \rho' - 0.095 \cdot (\rho')^2 + 2.48 \times 10^{-4} \cdot (\rho')^3$$

where $\rho'$ is the resistivity in units of 10⁻⁸ Ω·m.

Error propagation for temperature:

$$\Delta T = \left|\frac{dT}{d\rho'}\right| \cdot \Delta \rho'$$

where:

$$\frac{dT}{d\rho'} = 38.1 - 2 \cdot 0.095 \cdot \rho' + 3 \cdot 2.48 \times 10^{-4} \cdot (\rho')^2$$

## Refraction Index Calculation

The refraction index is calculated using:

$$n = \sqrt{\left(1.1547 \cdot \sin\left(\frac{\text{Init} - \text{Angle}}{\text{Ratio}}\right) + 0.5\right)^2 + 0.75}$$

Error propagation for the inner term:

$$\Delta(\text{inner\_term}) = \sqrt{\left(\frac{\Delta \text{Init}}{\text{Ratio}}\right)^2 + \left(\frac{\Delta \text{Angle}}{\text{Ratio}}\right)^2 + \left(\frac{(\text{Init} - \text{Angle}) \cdot \Delta \text{Ratio}}{\text{Ratio}^2}\right)^2}$$

Error propagation for the sin term:

$$\Delta(\text{sin\_term}) = 1.1547 \cdot |\cos(\text{inner\_term})| \cdot \Delta(\text{inner\_term})$$

Error propagation for the refraction index:

$$\Delta n = \left|\frac{\text{sin\_term} \cdot \Delta(\text{sin\_term})}{n}\right|$$

## Wavelength Calculation

The wavelength is calculated using:

$$\lambda = \frac{3000}{\sqrt{A + B \cdot n + C \cdot n^2 + D \cdot n^3 + E \cdot n^4 + F \cdot n^5 + G \cdot n^6 + H \cdot n^7 + I \cdot n^8}}$$

where the coefficients A through I are constants from the calibration table.

Error propagation for wavelength:

$$\Delta \lambda = \left|\frac{d\lambda}{dn}\right| \cdot \Delta n$$

where:

$$\frac{d\lambda}{dn} = -3000 \cdot \frac{\frac{d(\text{denominator})}{dn}}{\text{denominator}^2}$$

and:

$$\frac{d(\text{denominator})}{dn} = \frac{B + 2C \cdot n + 3D \cdot n^2 + 4E \cdot n^3 + 5F \cdot n^4 + 6G \cdot n^5 + 7H \cdot n^6 + 8I \cdot n^7}{2 \cdot \text{denominator}}$$

## Wien's Law Analysis

Wien's displacement law states that:

$$\lambda_{\text{max}} \cdot T = \text{constant}$$

This means:

$$\lambda_{\text{max}} = \text{constant} \cdot \frac{1}{T}$$

Error propagation for inverse temperature:

$$\Delta\left(\frac{1}{T}\right) = \frac{\Delta T}{T^2}$$

The Wien constant is determined through ODR fitting of wavelength vs. 1/T data, and the theoretical value is 2.9 mm·K. 