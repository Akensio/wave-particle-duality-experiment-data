Expansion – Simulation of Blackbody Radiation Spectrum through Diffraction Grating

In the blackbody radiation experiment, we used a prism to scan the spectrum of the body. As is known, a prism has limitations in its ability to separate closely spaced frequencies. In this section, we discuss replacing the prism with a diffraction grating. There are methods to use a diffraction grating in spectroscopy with convex mirrors or lenses and a screen, but here we focus specifically on a direct substitution of the prism with a diffraction grating in the blackbody radiation experimental setup, employing a circular table and a photometer.

For a monochromatic beam of wavelength $( \lambda )$ entering perpendicularly into a grating and passing through $( N )$ slits, spaced at distance $( d )$, and measuring illumination intensity at a large distance from the grating, the intensity as a function of the exit angle $( \theta )$ is given by:

$$ I(\theta) = I_0 \left(\frac{\sin^2\left(\frac{N\pi d}{\lambda}\sin\theta\right)}{\sin^2\left(\frac{\pi d}{\lambda}\sin\theta\right)}\right) $$

Intensity maxima occur at angles $( \theta )$ satisfying $( d\sin\theta = m\lambda )$, where $( m )$ is an integer.

For a large number of slits $( N )$, we have the approximation:

$$ I(\theta) \approx I_0 N^2 \delta\left(\frac{d}{\lambda}\sin\theta - m\right) $$

where $( \delta )$ is the Dirac delta function.

To visualize the resulting intensity pattern, we initially consider the first-order maxima (ignoring $( m=0 )$). For an incoming beam to the grating with intensity $( I_0 )$, the intensity on the screen after the grating at an angle $( \theta )$ for the wavelength that satisfies the maximum condition is:

$$ I(\theta) = I_0 N^2 $$

For any order $( m )$, the intensity of the $( m )$-th order maximum is:

$$ I_m(\theta) = I(\lambda_m) = I_0 N^2 $$

and the total intensity becomes:

$$ I_{\text{total}}(\theta) = \sum_{m} I_m(\theta) $$

For blackbody radiation according to Planck’s law at temperature $( T )$, we obtain:

$$ I(\lambda, T) = \frac{2hc^2}{\lambda^5} \frac{1}{e^{\frac{hc}{\lambda k_B T}} - 1} $$

We note that for higher-order $( m )$, the denominator exponent of $( e^{\frac{hc}{\lambda k_B T}} )$ increases significantly. Thus, for practical purposes, considering only a finite number of terms provides a sufficient approximation.

The condition $( d\sin\theta = m\lambda )$ implies that to observe a maximum for a desired wavelength, we depend on the slit spacing $( d )$. Therefore, we want $( d )$ to be greater than the wavelength $( \lambda )$ we wish to measure.

On the other hand, if $( d )$ is excessively large relative to $( \lambda )$, numerous maxima points will appear. For the same angle $( \theta )$, maxima from multiple wavelengths may overlap. Consequently, a beam containing a continuum of wavelengths would produce an intensity pattern that is spread out, making analysis difficult.

Now, let's discuss the angular separation between two maxima for the same wavelength:

$$\Delta\theta = \frac{\lambda}{d\cos\theta}$$

If $( d )$ is significantly greater than $( \lambda )$ and the order $( m )$ is small, we approximate:

$$\Delta\theta \approx \frac{\lambda}{d}$$

Thus, for short wavelengths, if we have large slit spacing relative to $( \lambda )$, we obtain a smeared intensity pattern that complicates distinguishing and measuring the intensity.

The simulation setup details will be discussed in the experiment procedure section.

