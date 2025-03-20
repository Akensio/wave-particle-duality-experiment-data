#!/usr/bin/env python3
"""Run diffraction grating simulation."""

import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from simulation.ui.interactive import run_interactive_simulation

def main():
    """Run interactive diffraction simulation."""
    run_interactive_simulation(use_infinite_slits=False)

if __name__ == "__main__":
    main() 