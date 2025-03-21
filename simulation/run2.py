#!/usr/bin/env python3
"""Run continuous spectrum diffraction grating simulation."""

import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from simulation.ui.continuous_interactive import run_continuous_spectrum_simulation

def main():
    """Run interactive continuous spectrum diffraction simulation."""
    run_continuous_spectrum_simulation()

if __name__ == "__main__":
    main() 