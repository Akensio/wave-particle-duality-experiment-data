import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, CheckButtons, RadioButtons

from simulation.core.physics import DiffractionPhysics
from simulation.visualization.wall_pattern import WallPatternRenderer
from simulation.core.config import (
    DEFAULT_GRATING_SPACING,
    DEFAULT_NUM_SLITS,
    WAVELENGTH_OPTIONS,
    COLOR_MAP,
    DEFAULT_WAVELENGTHS,
)


class InteractiveSimulation:
    """Interactive diffraction simulation with real-time parameter controls."""

    def __init__(self):
        """Initialize the interactive simulation."""
        # Initialize parameters
        self.grating_spacing = DEFAULT_GRATING_SPACING
        self.num_slits = DEFAULT_NUM_SLITS
        self.selected_wavelengths = DEFAULT_WAVELENGTHS.copy()

        # Create physics model
        self.physics = DiffractionPhysics(self.grating_spacing)

        # Create visualizer components
        self.wall_renderer = WallPatternRenderer()

        # Create the figure and UI
        self._create_figure()
        self._create_controls()
        self._update_plots()

    def _create_figure(self):
        """Create the figure and axes for the interactive simulation."""
        self.fig, (self.ax1, self.ax_total, self.ax2) = plt.subplots(3, 1, figsize=(12, 14), gridspec_kw={"height_ratios": [3, 1, 1]})

        suptitle = "Interactive Diffraction Grating Simulation"
        suptitle += f"\nGrating spacing: {self.grating_spacing*1e6:.1f} μm, Number of slits: {self.num_slits}"

        self.fig.suptitle(suptitle, fontsize=16)
        plt.subplots_adjust(left=0.15, bottom=0.25, right=0.85, top=0.95, hspace=0.4)

    def _create_controls(self):
        """Create sliders, checkboxes, and radio buttons for interaction."""
        # Create sliders
        ax_grating = plt.axes([0.15, 0.17, 0.7, 0.02])
        ax_slits = plt.axes([0.15, 0.13, 0.7, 0.02])

        self.grating_slider = Slider(ax_grating, "Grating Spacing (µm)", 0.5, 10.0, valinit=self.grating_spacing * 1e6)
        self.slits_slider = Slider(ax_slits, "Number of Slits", 2, 50, valinit=self.num_slits, valstep=1)

        # Create wavelength selection checkboxes
        ax_check = plt.axes([0.02, 0.45, 0.1, 0.2])
        check_labels = list(WAVELENGTH_OPTIONS.keys())
        active = [i for i, color in enumerate(check_labels) if color in self.selected_wavelengths]

        self.check = CheckButtons(ax_check, check_labels, [i in active for i in range(len(check_labels))])

        # Connect callbacks
        self.grating_slider.on_changed(self._update)
        self.slits_slider.on_changed(self._update)
        self.check.on_clicked(self._update_wavelengths)

    def _update_plots(self):
        """Update all plots based on current parameters."""
        # Clear axes
        self.ax1.clear()
        self.ax_total.clear()
        self.ax2.clear()

        # Get selected wavelengths data
        wavelengths = [WAVELENGTH_OPTIONS[color] for color in self.selected_wavelengths]
        colors = [COLOR_MAP[color] for color in self.selected_wavelengths]
        labels = [f"{color} ({wl*1e9:.0f} nm)" for color, wl in zip(self.selected_wavelengths, wavelengths)]

        if not wavelengths:
            self.ax1.set_title("Select at least one wavelength")
            return

        # Plot individual wavelengths
        self._plot_wavelengths(wavelengths, colors, labels)

        # Plot total intensity if multiple wavelengths selected
        if len(wavelengths) > 1:
            self._plot_combined_intensity(self.ax_total, wavelengths)

        # Set up axes
        self._setup_plot_axes()

        # Create wall pattern visualization
        self._create_wall_pattern(wavelengths, colors)

        self.fig.canvas.draw_idle()

    def _plot_wavelengths(self, wavelengths, colors, labels):
        """Plot intensity patterns for individual wavelengths."""
        for i, (wavelength, color, label) in enumerate(zip(wavelengths, colors, labels)):
            angles, intensity = self.physics.calculate_intensity_pattern(wavelength, num_slits=self.num_slits)

            self.ax1.plot(angles, intensity, color=color, label=label)

            # Add diffraction maxima markers
            maxima = self.physics.calculate_maxima_positions(wavelength)
            for m, angle in maxima:
                self.ax1.axvline(x=angle, color=color, linestyle="--", alpha=0.3)

    def _plot_combined_intensity(self, ax, wavelengths):
        """Plot the combined intensity from multiple wavelengths."""
        angles, total = self._calculate_combined_intensity(wavelengths)
        ax.plot(angles, total, "k-", linewidth=2)

    def _calculate_combined_intensity(self, wavelengths):
        """Calculate total intensity from multiple wavelengths."""
        # Get positions from first wavelength
        angles, _ = self.physics.calculate_intensity_pattern(wavelengths[0], num_slits=self.num_slits)

        # Sum contributions from all wavelengths
        total = np.zeros_like(angles)
        for wl in wavelengths:
            _, intensity = self.physics.calculate_intensity_pattern(wl, num_slits=self.num_slits)
            total += intensity

        # Normalize
        if np.max(total) > 0:
            total = total / np.max(total)

        return angles, total

    def _setup_plot_axes(self):
        """Configure axes labels, titles and ranges."""
        # Main plot
        self.ax1.set_xlim(-np.pi/2, np.pi/2)
        self.ax1.set_ylim(0, 1.05)
        self.ax1.set_xlabel("Angle θ (radians)")
        self.ax1.set_ylabel("Relative Intensity")
        
        # Add a secondary x-axis with angles in degrees
        ax1_deg = self.ax1.twiny()
        ax1_deg.set_xlim(-90, 90)
        ax1_deg.set_xlabel("Angle θ (degrees)")
        
        self.ax1.grid(True, alpha=0.3)
        self.ax1.legend()

        title = f"Diffraction Pattern ({self.num_slits} Slits)"
        title += f"\nGrating spacing: {self.grating_spacing*1e6:.1f} µm, Semi-circular screen at infinite distance"
        self.ax1.set_title(title)

        # Total intensity plot
        self.ax_total.set_xlim(-np.pi/2, np.pi/2)
        self.ax_total.set_ylim(0, 1.05)
        self.ax_total.set_xlabel("Angle θ (radians)")
        self.ax_total.set_ylabel("Total Intensity")
        self.ax_total.grid(True, alpha=0.3)
        self.ax_total.set_title("Combined Intensity (Sum of All Selected Wavelengths)")

    def _create_wall_pattern(self, wavelengths, colors):
        """Create the 2D wall pattern visualization."""
        self.wall_renderer.add_wall_pattern(self.fig, self.physics, wavelengths, colors, num_slits=self.num_slits, ax=self.ax2)
        self.ax2.set_title("Pattern on Wall")

    def _update(self, val):
        """Update callback for sliders."""
        # Get values from sliders
        self.grating_spacing = self.grating_slider.val * 1e-6  # Convert from µm to m
        self.num_slits = int(self.slits_slider.val)

        # Update physics model parameters
        self.physics.grating_spacing = self.grating_spacing

        # Update title
        suptitle = "Interactive Diffraction Grating Simulation"
        suptitle += f"\nGrating spacing: {self.grating_spacing*1e6:.1f} μm, Number of slits: {self.num_slits}"

        self.fig.suptitle(suptitle, fontsize=16)

        # Update plots
        self._update_plots()

    def _update_wavelengths(self, label):
        """Update callback for wavelength checkboxes."""
        if label in self.selected_wavelengths:
            self.selected_wavelengths.remove(label)
        else:
            self.selected_wavelengths.append(label)

        self._update_plots()


def run_interactive_simulation(use_infinite_slits=False):
    """Run the interactive simulation."""
    # Ignore the use_infinite_slits parameter (kept for compatibility)
    simulation = InteractiveSimulation()
    plt.show()
