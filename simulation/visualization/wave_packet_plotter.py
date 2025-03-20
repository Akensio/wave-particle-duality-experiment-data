import matplotlib.pyplot as plt
import numpy as np
from typing import Tuple, List, Optional

from simulation.core.physics import GaussianWavePacketPhysics


class WavePacketVisualizer:
    """Class for visualizing Gaussian wave packets and their diffraction."""
    
    def plot_wave_packet_spectrum(
        self,
        physics: GaussianWavePacketPhysics,
        center_wavelength: float,
        wavelength_width: float,
        ax: Optional[plt.Axes] = None
    ) -> Tuple[plt.Figure, plt.Axes]:
        """Plot the spectral distribution of a Gaussian wave packet."""
        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 6))
        else:
            fig = ax.figure
        
        # Calculate the spectral distribution
        wavelengths, spectrum = physics.calculate_spectral_distribution(
            center_wavelength, wavelength_width
        )
        
        # Convert wavelengths to nanometers for better readability
        wavelengths_nm = wavelengths * 1e9
        
        # Plot the spectrum
        ax.plot(wavelengths_nm, spectrum, 'b-', linewidth=2)
        ax.fill_between(wavelengths_nm, 0, spectrum, color='blue', alpha=0.2)
        
        # Add a marker for the center wavelength
        ax.axvline(x=center_wavelength * 1e9, color='r', linestyle='--', 
                 label=f'Center: {center_wavelength*1e9:.1f} nm')
        
        # Set plot properties
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('Relative Amplitude')
        ax.set_title(f'Gaussian Wave Packet Spectrum\nCenter: {center_wavelength*1e9:.1f} nm, Width: {wavelength_width*1e9:.1f} nm')
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        return fig, ax
    
    def plot_wave_packet_diffraction(
        self,
        physics: GaussianWavePacketPhysics,
        center_wavelength: float,
        wavelength_width: float,
        num_slits: int,
        compare_single: bool = True,
        ax: Optional[plt.Axes] = None
    ) -> Tuple[plt.Figure, plt.Axes]:
        """
        Plot the diffraction pattern of a Gaussian wave packet compared to a single wavelength.
        
        Args:
            physics: GaussianWavePacketPhysics instance
            center_wavelength: Center wavelength of the packet in meters
            wavelength_width: Width of the wavelength distribution in meters
            num_slits: Number of slits in the diffraction grating
            compare_single: Whether to compare with a single wavelength diffraction
            ax: Optional axes to plot on
        
        Returns:
            Tuple of (figure, axes)
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=(12, 8))
        else:
            fig = ax.figure
        
        # Calculate the wave packet diffraction pattern
        positions, packet_intensity = physics.calculate_wave_packet_pattern(
            center_wavelength, wavelength_width, num_slits=num_slits
        )
        
        # Plot the wave packet diffraction pattern
        ax.plot(positions, packet_intensity, 'b-', linewidth=2, 
              label=f'Gaussian packet (width: {wavelength_width*1e9:.1f} nm)')
        
        # Compare with single wavelength if requested
        if compare_single:
            positions, single_intensity = physics.calculate_intensity_pattern(
                center_wavelength, num_slits=num_slits
            )
            ax.plot(positions, single_intensity, 'r--', linewidth=1.5,
                  label=f'Single wavelength ({center_wavelength*1e9:.1f} nm)')
        
        # Set plot properties
        ax.set_xlim(-physics.screen_width/2, physics.screen_width/2)
        ax.set_ylim(0, 1.05)
        ax.set_xlabel("Position on Screen (m)")
        ax.set_ylabel("Relative Intensity")
        ax.set_title(f"Diffraction Pattern from {num_slits} Slits\n"
                    f"Center wavelength: {center_wavelength*1e9:.1f} nm, "
                    f"Grating spacing: {physics.grating_spacing*1e6:.1f} μm")
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        return fig, ax
    
    def plot_wave_packet_evolution(
        self,
        physics: GaussianWavePacketPhysics,
        center_wavelength: float,
        wavelength_width: float,
        position: float,
        time_range: Tuple[float, float],
        num_times: int = 200,
        dispersion_factor: float = 0.0,
        ax: Optional[plt.Axes] = None
    ) -> Tuple[plt.Figure, plt.Axes]:
        """
        Plot the time evolution of a Gaussian wave packet at a specific position.
        
        Args:
            physics: GaussianWavePacketPhysics instance
            center_wavelength: Center wavelength of the packet in meters
            wavelength_width: Width of the wavelength distribution in meters
            position: Position to evaluate the wave packet (in meters)
            time_range: Tuple of (start_time, end_time) in seconds
            num_times: Number of time points to plot
            dispersion_factor: Factor controlling dispersion effects
            ax: Optional axes to plot on
            
        Returns:
            Tuple of (figure, axes)
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=(12, 6))
        else:
            fig = ax.figure
        
        # Generate time points
        times = np.linspace(time_range[0], time_range[1], num_times)
        
        # Calculate the wave packet evolution
        times, amplitudes = physics.calculate_wave_packet_time_evolution(
            center_wavelength, wavelength_width, times, position, dispersion_factor
        )
        
        # Convert times to picoseconds for better readability
        times_ps = np.array(times) * 1e12
        
        # Plot the wave packet evolution
        ax.plot(times_ps, amplitudes, 'b-', linewidth=2)
        
        # Plot an envelope if desired
        envelope = np.abs(amplitudes)
        ax.plot(times_ps, envelope, 'r--', linewidth=1, alpha=0.7, label='Envelope')
        ax.plot(times_ps, -envelope, 'r--', linewidth=1, alpha=0.7)
        
        # Set plot properties
        ax.set_xlabel('Time (ps)')
        ax.set_ylabel('Amplitude')
        ax.set_title(f'Wave Packet Evolution at Position x = {position:.2e} m\n'
                   f'Center wavelength: {center_wavelength*1e9:.1f} nm, '
                   f'Width: {wavelength_width*1e9:.1f} nm')
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        return fig, ax
    
    def plot_wave_packet_full_analysis(
        self,
        physics: GaussianWavePacketPhysics,
        center_wavelength: float,
        wavelength_width: float,
        num_slits: int
    ) -> plt.Figure:
        """
        Create a comprehensive visualization of the wave packet properties.
        
        Args:
            physics: GaussianWavePacketPhysics instance
            center_wavelength: Center wavelength of the packet in meters
            wavelength_width: Width of the wavelength distribution in meters
            num_slits: Number of slits in the diffraction grating
            
        Returns:
            Figure with the comprehensive visualization
        """
        # Create a figure with 2x2 subplots
        fig, axs = plt.subplots(2, 2, figsize=(15, 12))
        
        # Plot the spectral distribution
        self.plot_wave_packet_spectrum(physics, center_wavelength, wavelength_width, ax=axs[0, 0])
        
        # Plot the diffraction pattern
        self.plot_wave_packet_diffraction(
            physics, center_wavelength, wavelength_width, num_slits, compare_single=True, ax=axs[0, 1]
        )
        
        # Plot the wave packet time evolution at the center of the screen
        self.plot_wave_packet_evolution(
            physics, center_wavelength, wavelength_width, 0.0, 
            time_range=(0, 1e-13), dispersion_factor=0.0, ax=axs[1, 0]
        )
        
        # Plot the wave packet time evolution with dispersion
        self.plot_wave_packet_evolution(
            physics, center_wavelength, wavelength_width, 0.0, 
            time_range=(0, 1e-13), dispersion_factor=1e20, ax=axs[1, 1]
        )
        axs[1, 1].set_title(axs[1, 1].get_title() + '\nWith Dispersion')
        
        fig.tight_layout()
        return fig 