import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, CheckButtons, RadioButtons
from matplotlib.gridspec import GridSpec

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
        # Create figure with constrained layout for better spacing
        self.fig = plt.figure(figsize=(12, 18), constrained_layout=False)
        
        # Use GridSpec with more space between plots
        gs = GridSpec(4, 1, figure=self.fig, height_ratios=[3, 1, 1, 0.6], hspace=1.0)
        
        # Create axes with proper spacing
        self.ax1 = self.fig.add_subplot(gs[0])  # Main diffraction pattern
        self.ax_total = self.fig.add_subplot(gs[1])  # Combined intensity
        self.ax2 = self.fig.add_subplot(gs[2])  # Wall pattern
        
        # Create an empty axis for controls area
        self.ax_controls = self.fig.add_subplot(gs[3])
        self.ax_controls.axis('off')  # Hide this axis, just use it for spacing

        # Set specific margins
        self.fig.subplots_adjust(left=0.15, right=0.85, top=0.95, bottom=0.05, hspace=0.6)
        
        # Add subtitle at the top of the first plot
        subtitle = f"Grating spacing: {self.grating_spacing*1e6:.1f} μm, Number of slits: {self.num_slits}"
        self.subtitle_text = self.ax1.text(0.05, 0.98, subtitle, 
                                         fontsize=10, ha='left', va='top',
                                         transform=self.ax1.transAxes,
                                         bbox=dict(facecolor='white', alpha=0.7, pad=5))

    def _create_controls(self):
        """Create sliders, checkboxes, and radio buttons for interaction."""
        # Create sliders with more space - position relative to figure instead of absolute coordinates
        slider_width = 0.65
        slider_height = 0.02
        slider_left = 0.25
        
        # Position sliders in the controls area at the bottom
        ax_grating = plt.axes([slider_left, 0.07, slider_width, slider_height])
        ax_slits = plt.axes([slider_left, 0.03, slider_width, slider_height])

        self.grating_slider = Slider(ax_grating, "Grating Spacing (µm)", 0.5, 10.0, valinit=self.grating_spacing * 1e6)
        self.slits_slider = Slider(ax_slits, "Number of Slits", 2, 50, valinit=self.num_slits, valstep=1)

        # Create wavelength selection checkboxes next to the sliders
        # Reduce width while maintaining height
        ax_check = plt.axes([0.02, 0.01, 0.1, 0.1])
        check_labels = list(WAVELENGTH_OPTIONS.keys())
        active = [i for i, color in enumerate(check_labels) if color in self.selected_wavelengths]
        
        # Create the checkboxes with standard labels
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
            self.ax1.text(0.5, 0.5, "Select at least one wavelength", 
                         ha='center', va='center', fontsize=12, transform=self.ax1.transAxes)
            return

        # Plot individual wavelengths
        self._plot_wavelengths(wavelengths, colors, labels)

        # Plot total intensity if multiple wavelengths selected
        if len(wavelengths) > 1:
            self._plot_combined_intensity(self.ax_total, wavelengths)
        else:
            # If only one wavelength, still setup the axis but don't plot
            self.ax_total.set_xlim(-np.pi/2, np.pi/2)
            self.ax_total.set_ylim(0, 1.05)
            self.ax_total.text(0, 0.5, "Select more than one wavelength to see combined intensity", 
                               ha='center', va='center', fontsize=12)

        # Set up axes
        self._setup_plot_axes()

        # Create wall pattern visualization
        self._create_wall_pattern(wavelengths, colors)
        
        # Update subtitle
        subtitle = f"Grating spacing: {self.grating_spacing*1e6:.1f} μm, Number of slits: {self.num_slits}"
        self.subtitle_text.set_text(subtitle)

        self.fig.canvas.draw_idle()

    def _plot_wavelengths(self, wavelengths, colors, labels):
        """Plot intensity patterns for individual wavelengths."""
        for i, (wavelength, color, label) in enumerate(zip(wavelengths, colors, labels)):
            angles, intensity = self.physics.calculate_intensity_pattern(wavelength, num_slits=self.num_slits)

            self.ax1.plot(angles, intensity, color=color, label=label, linewidth=2)

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
        self.ax1.set_xlabel("Angle θ (radians)", fontsize=10)
        self.ax1.set_ylabel("Relative Intensity", fontsize=10)
        
        # Add a secondary x-axis with angles in degrees
        ax1_deg = self.ax1.twiny()
        ax1_deg.set_xlim(-90, 90)
        ax1_deg.set_xlabel("Angle θ (degrees)", fontsize=10)
        
        self.ax1.grid(True, alpha=0.3)
        
        # Move legend to the upper right corner inside the plot for more stability
        self.ax1.legend(loc='upper right', fontsize=9)

        # Total intensity plot without title
        self.ax_total.set_xlim(-np.pi/2, np.pi/2)
        self.ax_total.set_ylim(0, 1.05)
        self.ax_total.set_xlabel("Angle θ (radians)", fontsize=10)
        self.ax_total.set_ylabel("Total Intensity", fontsize=10)
        self.ax_total.grid(True, alpha=0.3)

        # Add degree labels to total intensity axis
        ax_total_deg = self.ax_total.twiny()
        ax_total_deg.set_xlim(-90, 90)
        ax_total_deg.set_xlabel("Angle θ (degrees)", fontsize=10)

    def _create_wall_pattern(self, wavelengths, colors):
        """Create the 2D wall pattern visualization."""
        self.wall_renderer.add_wall_pattern(self.fig, self.physics, wavelengths, colors, num_slits=self.num_slits, ax=self.ax2)

    def _update(self, val):
        """Update callback for sliders."""
        # Get values from sliders
        self.grating_spacing = self.grating_slider.val * 1e-6  # Convert from µm to m
        self.num_slits = int(self.slits_slider.val)

        # Update physics model parameters
        self.physics.grating_spacing = self.grating_spacing

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
