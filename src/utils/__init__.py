"""
Utility functions for spectroscopy calculations and data processing.
"""

from .spectrum_utils import scale_spectrum_to_magnitude, load_star_spectra_from_fits, load_lamp_spectra_from_txt, resample_spectrum_resolution

__all__ = ['scale_spectrum_to_magnitude', 'load_star_spectra_from_fits', 'load_lamp_spectra_from_txt', 'resample_spectrum_resolution']