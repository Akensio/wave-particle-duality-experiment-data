import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons
from matplotlib.gridspec import GridSpec
from functools import partial

from simulation.core.physics import ContinuousSpectrumPhysics
from simulation.core.config import DEFAULT_GRATING_SPACING


class ContinuousSpectrumSimulation:
    """Interactive simulation for diffraction with continuous wavelength spectrum."""

    def __init__(self):
        """Initialize the interactive simulation."""
        # Initialize parameters
        self.grating_spacing = DEFAULT_GRATING_SPACING
        self.temperature = 3000  # Starting temperature in Kelvin
        self.num_wavelengths = 200  # Number of wavelength samples
        self.max_order = 1  # Only consider first order diffraction
        self.wavelength_range = (380e-9, 3000e-9)  # Extended range up to 3000nm
        self.spectrum_type = "blackbody"  # Default spectrum type
        
        # Create physics model
        self.physics = ContinuousSpectrumPhysics(self.grating_spacing)
        
        # Create the figure and UI
        self._create_figure()
        self._create_controls()
        self._update_plots()

    def _create_figure(self):
        """Create the figure and axes for the interactive simulation."""
        self.fig = plt.figure(figsize=(12, 14), constrained_layout=False)
        
        # Use GridSpec with spacings
        gs = GridSpec(3, 1, figure=self.fig, height_ratios=[2, 2, 0.6], hspace=0.4)
        
        # Create axes
        self.ax_spectrum = self.fig.add_subplot(gs[0])  # Spectrum plot
        self.ax_diffraction = self.fig.add_subplot(gs[1])  # Diffraction pattern
        
        # Create an empty axis for controls area
        self.ax_controls = self.fig.add_subplot(gs[2])
        self.ax_controls.axis('off')
        
        # Set specific margins
        self.fig.subplots_adjust(left=0.15, right=0.85, top=0.95, bottom=0.05)
        
        # Add title information
        self.title_text = self.ax_spectrum.text(
            0.05, 0.98, "", fontsize=10, ha='left', va='top',
            transform=self.ax_spectrum.transAxes,
            bbox=dict(facecolor='white', alpha=0.7, pad=5)
        )

    def _create_controls(self):
        """Create sliders and radio buttons for interaction."""
        # Create sliders
        slider_width = 0.65
        slider_height = 0.02
        slider_left = 0.25
        
        # Position sliders in the controls area
        ax_grating = plt.axes([slider_left, 0.11, slider_width, slider_height])
        ax_temp = plt.axes([slider_left, 0.07, slider_width, slider_height])
        ax_order = plt.axes([slider_left, 0.03, slider_width, slider_height])

        self.grating_slider = Slider(ax_grating, "Grating Spacing (µm)", 0.5, 10.0, 
                                    valinit=self.grating_spacing * 1e6)
        self.temp_slider = Slider(ax_temp, "Temperature (K)", 1000, 10000, 
                                 valinit=self.temperature, valstep=100)
        self.order_slider = Slider(ax_order, "Max Order", 1, 3, 
                                  valinit=self.max_order, valstep=1)
        
        # Create spectrum type selection
        ax_radio = plt.axes([0.02, 0.01, 0.15, 0.1])
        self.radio = RadioButtons(
            ax_radio, ('Blackbody', 'Uniform'), active=0
        )
        
        # Connect callbacks
        self.grating_slider.on_changed(self._update)
        self.temp_slider.on_changed(self._update)
        self.order_slider.on_changed(self._update)
        self.radio.on_clicked(self._update_spectrum_type)

    def _update_spectrum_type(self, label):
        """Update the spectrum type based on radio button selection."""
        self.spectrum_type = label.lower()
        self._update_plots()

    def _get_spectrum_function(self):
        """Get the appropriate spectrum function based on the selected type."""
        if self.spectrum_type == "blackbody":
            return partial(self.physics.planck_spectrum, temperature=self.temperature)
        elif self.spectrum_type == "uniform":
            return lambda wavelengths: np.ones_like(wavelengths)
        else:
            # Default to uniform spectrum
            return lambda wavelengths: np.ones_like(wavelengths)

    def _update_plots(self):
        """Update all plots based on current parameters."""
        # Clear axes
        self.ax_spectrum.clear()
        self.ax_diffraction.clear()
        
        # Sample wavelengths for the spectrum plot
        wavelengths = np.linspace(self.wavelength_range[0], self.wavelength_range[1], self.num_wavelengths)
        
        # Get spectrum function and calculate spectrum
        spectrum_function = self._get_spectrum_function()
        intensities = spectrum_function(wavelengths)
        
        # Plot the spectrum
        self._plot_spectrum(wavelengths, intensities)
        
        # Calculate and plot the diffraction pattern
        self._plot_diffraction_pattern(spectrum_function)
        
        # Update title information
        self._update_title()
        
        # Refresh the figure
        self.fig.canvas.draw_idle()

    def _plot_spectrum(self, wavelengths, intensities):
        """Plot the wavelength spectrum."""
        # Convert wavelengths to nm for display
        wavelengths_nm = wavelengths * 1e9
        
        # Plot the spectrum
        self.ax_spectrum.plot(wavelengths_nm, intensities, 'k-', linewidth=2)
        
        # Add color background to represent visible spectrum
        self._add_visible_spectrum_background(self.ax_spectrum)
        
        # Setup axis labels and limits
        self.ax_spectrum.set_xlim(self.wavelength_range[0] * 1e9, self.wavelength_range[1] * 1e9)
        self.ax_spectrum.set_ylim(0, 1.05)
        self.ax_spectrum.set_xlabel("Wavelength (nm)", fontsize=10)
        self.ax_spectrum.set_ylabel("Relative Intensity", fontsize=10)
        self.ax_spectrum.grid(True, alpha=0.3)
        
        # Set title based on spectrum type
        if self.spectrum_type == "blackbody":
            self.ax_spectrum.set_title(f"Blackbody Radiation Spectrum (T = {self.temperature} K)", fontsize=12)
        else:
            self.ax_spectrum.set_title("Uniform Intensity Spectrum", fontsize=12)

    def _add_visible_spectrum_background(self, ax):
        """Add color background to represent the visible spectrum."""
        # Define approximate colors for visible spectrum
        colors = ['violet', 'blue', 'cyan', 'green', 'yellow', 'orange', 'red']
        bounds = np.linspace(380, 750, len(colors) + 1)
        
        # Add color bands to background for visible spectrum only
        for i in range(len(colors)):
            ax.axvspan(bounds[i], bounds[i+1], color=colors[i], alpha=0.15)
            
        # Add a light shading for infrared region
        ax.axvspan(750, 3000, color='lightgray', alpha=0.1, label='IR')
        
        # Add a vertical line to mark the boundary between visible and IR
        ax.axvline(x=750, color='gray', linestyle='--', alpha=0.5)
        ax.text(760, 0.95, 'IR', color='gray', fontsize=8)

    def _plot_diffraction_pattern(self, spectrum_function):
        """Calculate and plot the diffraction pattern for the current spectrum."""
        # Calculate the diffraction pattern
        angles, intensity = self.physics.calculate_continuous_spectrum_pattern(
            self.wavelength_range, 
            self.num_wavelengths,
            spectrum_function,
            max_order=self.max_order
        )
        
        # Plot the diffraction pattern
        self.ax_diffraction.plot(angles, intensity, 'k-', linewidth=2)
        
        # Setup axis labels and limits
        self.ax_diffraction.set_xlim(-np.pi/2, np.pi/2)
        self.ax_diffraction.set_ylim(0, 1.05)
        self.ax_diffraction.set_xlabel("Angle θ (radians)", fontsize=10)
        self.ax_diffraction.set_ylabel("Relative Intensity", fontsize=10)
        self.ax_diffraction.grid(True, alpha=0.3)
        
        # Add a secondary x-axis with angles in degrees
        ax_deg = self.ax_diffraction.twiny()
        ax_deg.set_xlim(-90, 90)
        ax_deg.set_xlabel("Angle θ (degrees)", fontsize=10)
        
        # Set title
        self.ax_diffraction.set_title(
            f"First-Order Diffraction Pattern (Order: {self.max_order})", fontsize=12
        )

    def _update_title(self):
        """Update the title text with current parameters."""
        title_str = (
            f"Grating spacing: {self.grating_spacing*1e6:.1f} μm, "
            f"Temperature: {self.temperature} K, "
            f"Max Order: {self.max_order}"
        )
        self.title_text.set_text(title_str)

    def _update(self, val):
        """Update callback for sliders."""
        # Get values from sliders
        self.grating_spacing = self.grating_slider.val * 1e-6  # Convert from µm to m
        self.temperature = self.temp_slider.val
        self.max_order = int(self.order_slider.val)
        
        # Update physics model parameters
        self.physics.grating_spacing = self.grating_spacing
        
        # Update plots
        self._update_plots()


def run_continuous_spectrum_simulation():
    """Run the interactive continuous spectrum simulation."""
    simulation = ContinuousSpectrumSimulation()
    plt.show() 