Expansion – Simulation of Blackbody Radiation Spectrum Through Diffraction Grating

In the blackbody radiation experiment, we used a prism to scan the spectrum emitted by the body. As known, a prism has limitations in distinguishing closely spaced frequencies. In this section, we discuss the idea of replacing the prism in the experiment with a diffraction grating. Several methods exist for incorporating diffraction gratings into spectroscopy setups using convex mirrors, lenses, or screens. However, this expansion specifically discusses directly replacing the prism with a diffraction grating in the original experimental setup used previously for the blackbody radiation experiment, consisting of a circular table and a detector measuring light intensity.

For a monochromatic light beam of wavelength $\lambda$, entering perpendicularly into a diffraction grating consisting of $N$ slits separated by distance $d$, measuring intensity at a large distance from the grating, the intensity as a function of exit angle $\theta$ is given by:

$$I(\theta) = I_0\left[\frac{\sin\left(\frac{N\pi d}{\lambda}\sin\theta\right)}{\sin\left(\frac{\pi d}{\lambda}\sin\theta\right)}\right]^2$$

Intensity maxima occur at angles $\theta$ satisfying:

$$d\sin\theta = m\lambda$$

For a large number of slits ($N \gg 1$), the following approximation holds:

$$I(\theta)\approx I_0 N^2\sum_{m=-\infty}^{\infty}\delta\left(\frac{d}{\lambda}\sin\theta - m\right)$$

where $\delta$ is the Dirac delta function. The intuitive meaning of this sum is: "High intensity occurs at each interference maximum, corresponding to angles satisfying $d\sin\theta = m\lambda$, and intensity is essentially zero everywhere else."

For the purpose of our discussion, we may omit the factor $N^2$, since we are not interested in specific intensity levels but rather in qualitative properties of the resulting intensity distribution.

⸻

Now we want to understand what the intensity pattern captured by the sensor will look like. First, let's examine the pattern resulting from first-order maxima (we ignore the order $m=0$). For a beam containing a continuous range of wavelengths entering the grating, having intensity $I(\lambda)$, the intensity measured at angle $\theta$ after the grating corresponds to the wavelength satisfying the maximum condition, i.e., $\lambda = d\sin\theta$:

$$I(\theta)d\theta = I(\lambda)d\lambda$$

(Technically, the calculations here include minor simplifications—we neglect certain details, such as the fact that $I(\lambda)$ is actually intensity density per wavelength, and we omit the derivative $d\sin\theta/d\lambda$. These simplifications effectively cancel out, preserving the integral of intensity over angle, representing total energy.)

For an arbitrary order $m$, the intensity of the $m$-th maximum is given by:

$$I_m(\theta)d\theta = I(\lambda_m)d\lambda_m$$

For blackbody radiation at temperature $T$, according to Planck's law, we have:

$$I(\lambda,T) = \frac{2\pi c^2 h}{\lambda^5}\frac{1}{e^{\frac{hc}{\lambda k T}} - 1}$$

Substituting $\lambda_m = \frac{d \sin\theta}{m}$, explicitly gives:

$$I_m(\theta) = \frac{2\pi c^2 h}{\lambda_m^5}\frac{1}{e^{\frac{hc}{\lambda_m k T}} - 1}
= \frac{2\pi c^2 h}{\left(\frac{d \sin\theta}{m}\right)^5}\frac{1}{e^{\frac{hc}{\left(\frac{d \sin\theta}{m}\right) k T}} - 1}
= \frac{2\pi c^2 h m^5}{d^5 \sin^5\theta}\frac{1}{e^{\frac{h c m}{d \sin\theta k T}} - 1}$$

⸻

The meaning of the condition $d>\lambda$ is that to observe a maximum for a desired wavelength, we depend on the slit spacing $d$. Specifically, we prefer $d$ to be greater than the wavelength $\lambda$ we're measuring.

On the other hand, if the spacing $d$ is too large compared to $\lambda$, multiple maxima from many wavelengths appear at the same angle $\theta$. Thus, an incoming beam containing a continuum of wavelengths will produce a blurred intensity pattern, complicating analysis.

Now, let's examine the angular difference between two maxima for the same wavelength. We have:

$$d\sin\theta = m\lambda$$

For $d\gg\lambda$ and small orders $m$, we obtain the approximation:

$$\Delta\theta\approx\frac{\lambda}{d}$$

Thus, at short wavelengths, if we have gratings with slit spacing $d$ much greater than $\lambda$, the resulting intensity pattern becomes smeared, making it difficult to distinguish maxima and measure intensity clearly.

⸻

