import astropy.units as u

wl_blue = 460 * u.nm
wl_green = 580 * u.nm
wl_red = 700 * u.nm

wl_range_low = wl_green - wl_blue
wl_range_high = wl_red - wl_green

# spectrum
sensor_coverage = 28 * u.mm

spot_blue = 182 * u.um
spot_green = 95 * u.um
spot_red = 283 * u.um

delta_low = 9.8 * u.mm
delta_high = 11 * u.mm

bins_low = delta_low / spot_blue
bins_high = delta_high / spot_red

delta_wl_low = wl_range_low / bins_low
delta_wl_high = wl_range_high / bins_high

resolution_low = wl_blue / delta_wl_low
resolution_high = wl_red / delta_wl_high

print(f"Resolution blue:{resolution_low.decompose()}")
print(f"Resolution red:{resolution_high.decompose()}")

#sensor
npixels = 4096
eluminated_pixels = 3100
pixel_pitch = 9 * u.um
dark_current = 0.02