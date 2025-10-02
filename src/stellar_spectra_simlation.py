import pathlib
import sys
import importlib
from astropy.io import fits
import pandas as pd
import numpy as np
import astropy.units as u
import astropy.constants as const
from src.utils.spectrum_utils import *
import matplotlib.pyplot as plt

# calculate spectral resolution of the device
wl_blue = 460 * u.nm
wl_green = 580 * u.nm
wl_red = 700 * u.nm

wl_range_low = wl_green - wl_blue
wl_range_high = wl_red - wl_green
wl_range = wl_red - wl_blue

# spectrum
sensor_coverage = 28 * u.mm

spot_blue = 182 * u.um
spot_green = 95 * u.um
spot_red = 283 * u.um

delta_low = 9.8 * u.mm
delta_high = 11 * u.mm


bins_low = (delta_low / spot_blue).decompose()
bins_high = (delta_high / spot_red).decompose()

delta_wl_low = wl_range_low / bins_low
delta_wl_high = wl_range_high / bins_high

resolution_low = wl_blue / delta_wl_low
resolution_high = wl_red / delta_wl_high

print(f"Resolution blue:{resolution_low.decompose()}")
print(f"Resolution red:{resolution_high.decompose()}")

#sensor
npixels = 4096
iluminated_pixels = 3100
pixel_pitch = 9 * u.um
dark_current = 0.02
full_well_capacity = 1e6
read_out_noise = 10
integration_time = 60 * 60
quantum_efficiency = 0.6

pixels_per_wavelength = iluminated_pixels / wl_range
wavelength_per_pixel = wl_range / iluminated_pixels


# telescope
telescope_diameter = 0.508
telescope_obstruction = 0.39
telecsope_effective_area = (telescope_diameter/2)**2 * np.pi *(1-telescope_obstruction)

# optical efficency
telescope_m1 = 0.95
telescope_m2 = 0.95

dmd = 0.66
folding_mirror = 0.85
colimator_lens = 0.99
grating = 0.4
focusing_lens = 0.99
window = 0.99

transmission_efficency = telescope_m1 * telescope_m2 * dmd * folding_mirror \
    * colimator_lens * grating * focusing_lens * window
print(f"Total transmission efficency: {transmission_efficency}")


# load star data
data_dir = pathlib.Path("../s4/data/star_spectrum").resolve()
star_data_df, star_spectra_df = load_star_spectra_from_fits(data_dir)
star_spectra_df.columns = ["star 6000", "star 6200", "star 6700"]

print(star_spectra_df.head())

# load lamp spectra
lamp_dir = pathlib.Path("../s4/data/lamp_spectra").resolve()
lamp_spectra_df = load_lamp_spectra_from_txt(lamp_dir)
lamp_spectra_df.columns = ['sodium', 'led1','led2' ]
print(lamp_spectra_df.head())


cols = lamp_spectra_df.columns.tolist() # sodium, led, led
lamp_spectra_df['background'] = 1.5 * lamp_spectra_df['sodium'] \
    + 1 * lamp_spectra_df['led1'] \
    + 1 * lamp_spectra_df['led2']


stellar_magnitude = 11
star_mag= scale_spectrum_to_magnitude(star_spectra_df, target_magnitude=stellar_magnitude)
background_magnitude = 16
background_mag = scale_spectrum_to_magnitude(lamp_spectra_df[['background']], target_magnitude=background_magnitude)


#limit the wavelength range for further computation
lower_wavelength_limit, upper_wavelength_limit = 450, 700
star_spectra_range = star_mag.loc[lower_wavelength_limit:upper_wavelength_limit]
background_spectrum_range = background_mag.loc[lower_wavelength_limit:upper_wavelength_limit]

resolution = 200
delta_lambda = (int)(np.mean([lower_wavelength_limit, upper_wavelength_limit]) / resolution)
n_bins = (int)((upper_wavelength_limit - lower_wavelength_limit) / delta_lambda)


star_spectra_downsampled = bin_with_mean_index(star_spectra_range, n_bins)
background_spectrum_downsampled = bin_with_mean_index(background_spectrum_range, n_bins)
background_spectrum_downsampled.index = star_spectra_downsampled.index

star_spectra_downsampled.plot()

count_rate_star = flux_to_photon_count(star_spectra_downsampled)
count_rate_background = flux_to_photon_count(background_spectrum_downsampled)

count_rate_background.plot()
count_rate_star.plot()

integration_time = 60 * 60
oversampling_factor = 5

binning_factor = np.floor(iluminated_pixels / n_bins / oversampling_factor)

snr, n = calculate_snr(count_rate_star,
                    count_rate_background,
                    integration_time,
                    slit_size= 10,
                    collector_size= telecsope_effective_area,
                    optical_effeciency=transmission_efficency,
                    quantum_efficiency=quantum_efficiency,
                    wavelengths_per_pixel=wavelength_per_pixel,
                    dark_current_rate=dark_current,
                    readout_noise=read_out_noise,
                    full_well_capacity=full_well_capacity,
                    binning_factor=binning_factor)

snr['possible retreived spectrum'] = snr['count_rate'] / np.max(snr['count_rate'])

snr[['snr']].plot()

plt.figure(figsize=(7, 4))
plt.plot(lamp_spectra_df.index, lamp_spectra_df['sodium'], label='sodium')
plt.plot(lamp_spectra_df.index, lamp_spectra_df['led1'], label='led')
plt.title('Spectra of different lamp types')
plt.ylabel('Intensity')
plt.xlabel('wavelength (nm)')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.plot()
plt.savefig('plots/lamp_spectra_df.png', dpi=300, bbox_inches='tight')
plt.close()

star_1 = 'star 6200'
star_2 = 'star 6700'

plt.figure(figsize=(7, 4))
plt.plot(snr.index, snr[star_1]/snr[star_2], label = 'spectrum ratio')
#plt.plot(snr.index, snr[star_1]/np.max(snr[star_1]), label=star_1)
#plt.plot(snr.index, snr[star_2]/np.max(snr[star_2]), label=star_2)
plt.title('Relative difference in the spectra in different phases of the star')
plt.ylabel('Intensity')
plt.xlabel('wavelength (nm)')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.plot()
plt.savefig('plots/star_spectra_ratio.png', dpi=300, bbox_inches='tight')
plt.close()



plot_spectrum(snr, data_column='snr',
              label='',
              title='SNR 60 min observation',
              ylabel='SNR',
              save_path='plots/snr_60_min.png')

plot_spectrum(star_mag, data_column=star_1,
              label='',
              title='Stellar flux',
              ylabel='flux (w/m²/nm)',
              save_path='plots/simulated_star_spectrum.png')

plot_spectrum(background_spectrum_downsampled, data_column='background',
              label='',
              title='Background spectrum',
              ylabel='flux (w/m²/nm/arcsecond²)',
              save_path='plots/background spectrum.png')

plot_spectrum(snr, data_column='possible retreived spectrum',
              label='',
              title='Possible retreived spectrum',
              ylabel='Normalized intensity',
              save_path='plots/retreived_star_spectrum.png')



