import numpy as np
import pandas as pd
import astropy.units as u
import matplotlib.pyplot as plt


def scale_spectrum_to_magnitude(spectrum_df, target_magnitude, 
                                wavelength_range=(550, 650) * u.nm,
                                zero_magnitude_flux=3.64e-9 * u.erg / u.cm**2 / u.s / u.nm):
    """
    Scale a spectrum DataFrame to match a specific magnitude.
    
    Parameters:
    -----------
    spectrum_df : pd.DataFrame
        DataFrame with wavelength index in nm and spectrum values as columns
    target_magnitude : float
        Target magnitude to scale the spectrum to
    wavelength_range : tuple of astropy.units.Quantity, optional
        Wavelength range to average over for scaling factor computation.
        Default is (550, 650) nm (V-band approximate range)
    zero_magnitude_flux : astropy.units.Quantity, optional  
        Zero-magnitude flux density in CGS units (erg/cm²/s/nm).
        Default is Vega's V-band flux density
        
    Returns:
    --------
    pd.DataFrame
        Scaled spectrum with same structure as input, but in SI units (W/m²/nm)
    """
    # Make a copy of the input DataFrame
    spectrum_df_copy = spectrum_df.copy()
    
    # Ensure spectrum_df index has units (nm)
    if not hasattr(spectrum_df_copy.index, 'unit'):
        wavelength_with_units = spectrum_df_copy.index.values * u.nm
    else:
        wavelength_with_units = spectrum_df_copy.index
    
    # Convert wavelength range to same units as spectrum index
    wl_min = wavelength_range[0].to(u.nm)
    wl_max = wavelength_range[1].to(u.nm)
    
    # Find indices within the wavelength range
    mask = (wavelength_with_units >= wl_min) & (wavelength_with_units <= wl_max)
    
    if not np.any(mask):
        raise ValueError(f"No data points found in wavelength range {wl_min} to {wl_max}")
    
    # Calculate target flux density from magnitude
    # m = -2.5 * log10(F/F0) => F = F0 * 10^(-0.4 * m)
    target_flux_density = zero_magnitude_flux * 10**(-0.4 * target_magnitude)
    
    # Create output DataFrame with same structure
    scaled_df = spectrum_df_copy.copy()
    
    # Scale each spectrum column
    for col in spectrum_df_copy.columns:
        # Get spectrum values in the reference wavelength range
        spectrum_values_in_range = spectrum_df_copy.loc[mask, col].values
        
        # Calculate mean flux in the reference range (assuming input is in arbitrary units)
        mean_flux_in_range = np.mean(spectrum_values_in_range)
        
        if mean_flux_in_range <= 0:
            raise ValueError(f"Mean flux in reference range is non-positive for column {col}")
        
        # Calculate scaling factor to match target flux density
        # Convert target flux from CGS to SI units (W/m²/nm)
        target_flux_si = target_flux_density.to(u.W / u.m**2 / u.nm)
        scaling_factor = target_flux_si.value / mean_flux_in_range
        
        # Apply scaling factor to entire spectrum
        scaled_df[col] = spectrum_df_copy[col] * scaling_factor
    
    # Convert index to SI units (keep in nm for convenience, but could convert to m)
    # For now, keep wavelength in nm as it's more convenient for spectroscopy
    scaled_df.index = wavelength_with_units.to(u.nm).value
    scaled_df.index.name = 'wavelength_nm'
    
    return scaled_df


def load_star_spectra_from_fits(data_dir):
    """
    Load star spectra from FITS files in a directory.
    
    This function reads all FITS files in the specified directory, extracts wavelength
    and flux data, and returns a DataFrame with wavelength as index and spectra as columns.
    
    Parameters:
    -----------
    data_dir : str or PathLike
        Path to directory containing FITS files
        
    Returns:
    --------
    tuple: (star_data_df, star_spectra_df)
        star_data_df : pd.DataFrame
            Raw FITS file information with metadata
        star_spectra_df : pd.DataFrame  
            Processed spectra with wavelength index (nm) and float64 data type
    """
    import pathlib
    from astropy.io import fits
    
    data_dir = pathlib.Path(data_dir).resolve()
    fits_files = sorted(data_dir.glob("*.fits"))
    
    # Read all FITS files and extract HDU information
    rows = []
    for fpath in fits_files:
        try:
            with fits.open(fpath) as hdul:
                for idx, hdu in enumerate(hdul):
                    # try to get data, but keep None if missing
                    arr = hdu.data
                    header = hdu.header
                    rows.append({
                        "file": fpath.name,
                        "path": str(fpath),
                        "hdu_index": idx,
                        "hdu_name": hdu.name,
                        "data_shape": None if arr is None else getattr(arr, "shape", None),
                        "data": None if arr is None else np.array(arr),
                        "header": header,
                    })
        except Exception as e:
            rows.append({
                "file": fpath.name,
                "path": str(fpath),
                "hdu_index": None,
                "hdu_name": None,
                "data_shape": None,
                "data": None,
                "header": str(e),
            })
    
    # Create DataFrame with all FITS information
    star_data_df = pd.DataFrame(rows)
    
    # Build star spectra DataFrame from HDU index 0 entries
    hdu0 = star_data_df[star_data_df['hdu_index'] == 0].copy()
    
    # Find the row that corresponds to the WAVE file (case-insensitive, filename contains 'WAVE')
    wave_candidates = hdu0[hdu0['file'].str.upper().str.contains('WAVE')]
    if wave_candidates.empty:
        # fallback: try to find any HDU with name 'WAVE' or 'WAVE ' in header
        wave_candidates = hdu0[hdu0['hdu_name'].str.upper().fillna('').str.contains('WAVE')]
    
    if wave_candidates.empty:
        raise ValueError("No wavelength data found. Expected file containing 'WAVE' in filename or HDU name.")
    
    # Use the first wave candidate
    wave_row = wave_candidates.iloc[0]
    wave = wave_row['data']
    
    # Collect other spectra (exclude the WAVE file itself)
    spectra_rows = hdu0[hdu0.index != wave_row.name]
    
    # Build a DataFrame aligned on wavelength. We'll trim to the minimum length to avoid mismatches.
    arrays = {}
    min_len = len(wave)
    for _, row in spectra_rows.iterrows():
        arr = row['data']
        if arr is None:
            continue
        # flatten in case of extra dimensions and convert to float64
        arr_flat = np.ravel(arr).astype(np.float64)
        min_len = min(min_len, arr_flat.size)
        arrays[row['file']] = arr_flat
    
    # Trim wave and all arrays to min_len, ensure float64 for wavelength too
    wave_trim = np.ravel(wave)[:min_len].astype(np.float64)
    for k in list(arrays.keys()):
        arrays[k] = arrays[k][:min_len]
    
    # Create DataFrame with explicit float64 data type
    star_spectra_df = pd.DataFrame(arrays, index=wave_trim/10, dtype=np.float64)
    star_spectra_df.index.name = 'wavelength'
    
    return star_data_df, star_spectra_df


def load_lamp_spectra_from_txt(data_dir):
    """
    Load lamp spectra from text files in a directory.
    
    This function reads ASCII format lamp spectrum files with header information
    and wavelength/flux columns, returning a DataFrame with wavelength as index.
    
    Parameters:
    -----------
    data_dir : str or PathLike
        Path to directory containing lamp spectrum text files
        
    Returns:
    --------
    pd.DataFrame
        Lamp spectra with wavelength index (nm) and spectrum columns
    """
    import pathlib
    
    lamp_dir = pathlib.Path(data_dir).resolve()
    lamp_files = sorted(lamp_dir.glob("*.txt"))
    
    # Parse lamp spectrum files (ASCII format with header + W S columns)
    lamp_data = {}
    wavelengths = None
    
    for f in lamp_files:
        try:
            # Try different encodings to handle special characters
            lines = None
            for encoding in ['utf-8', 'latin-1', 'cp1252', 'iso-8859-1']:
                try:
                    with open(f, 'r', encoding=encoding) as file:
                        lines = file.readlines()
                    break
                except UnicodeDecodeError:
                    continue
            
            if lines is None:
                print(f"Could not decode {f.name} with any encoding")
                continue
            
            # Find start of data (after "Begin Processed Spectral Data" or "W	S" header)
            data_start_idx = None
            for i, line in enumerate(lines):
                if 'Begin Processed Spectral Data' in line or (line.strip().startswith('W') and 'S' in line):
                    data_start_idx = i + 1
                    break
            
            if data_start_idx is None:
                # fallback: skip all lines starting with #
                data_start_idx = 0
                for i, line in enumerate(lines):
                    if not line.strip().startswith('#') and len(line.strip()) > 0:
                        data_start_idx = i
                        break
            
            # Parse numerical data from data_start_idx onwards
            wave_vals = []
            flux_vals = []
            for line in lines[data_start_idx:]:
                line = line.strip()
                if not line or line.startswith('#') or line.startswith('W'):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        w = float(parts[0])
                        s = float(parts[1])
                        wave_vals.append(w)
                        flux_vals.append(s)
                    except ValueError:
                        continue
            
            if wave_vals:
                wave_arr = np.array(wave_vals)
                flux_arr = np.array(flux_vals)
                lamp_data[f.stem] = flux_arr  # use filename without extension as column name
                
                # Use first file's wavelength as reference
                if wavelengths is None:
                    wavelengths = wave_arr
                print(f"Loaded {f.name}: {len(wave_vals)} data points")
        
        except Exception as e:
            print(f"Error reading {f.name}: {e}")
    
    # Create DataFrame with wavelength as index and lamp spectra as columns
    if wavelengths is not None and lamp_data:
        # Ensure all spectra have same length by trimming to minimum
        min_len = min(len(wavelengths), min(len(arr) for arr in lamp_data.values()))
        wavelengths_trim = wavelengths[:min_len]
        
        trimmed_data = {}
        for name, arr in lamp_data.items():
            trimmed_data[name] = arr[:min_len]
        
        lamp_spectra_df = pd.DataFrame(trimmed_data, index=wavelengths_trim)
        lamp_spectra_df.index.name = 'wavelength'
        
        print(f"Built lamp_spectra_df with {len(lamp_data)} lamp spectra and {len(wavelengths_trim)} wavelength points")
        return lamp_spectra_df
    else:
        print("No lamp spectra data found")
        return pd.DataFrame()


def resample_spectrum_resolution(spectrum_df, new_wavelength_step, 
                                wavelength_range=None, method='auto'):
    """
    Resample spectrum to a different wavelength resolution.
    
    For decreased resolution (larger steps): averages values within each bin
    For increased resolution (smaller steps): uses nearest neighbor interpolation
    
    Parameters:
    -----------
    spectrum_df : pd.DataFrame
        DataFrame with wavelength index (in nm) and spectrum columns
    new_wavelength_step : float
        New wavelength step size in nm (e.g., 0.1 for 0.1 nm resolution)
    wavelength_range : tuple, optional
        (min_wavelength, max_wavelength) to limit the output range.
        If None, uses the full range of input spectrum
    method : str, optional
        Resampling method: 'auto' (default), 'bin_average', or 'interpolate'
        - 'auto': automatically choose based on resolution change
        - 'bin_average': force binning/averaging (for downsampling)
        - 'interpolate': force nearest neighbor interpolation
        
    Returns:
    --------
    pd.DataFrame
        Resampled spectrum with new wavelength resolution
    """
    # Determine current wavelength step
    current_wavelengths = spectrum_df.index.values
    if hasattr(spectrum_df.index, 'unit'):
        # Remove units if present
        current_wavelengths = spectrum_df.index.value
    
    current_step = np.median(np.diff(current_wavelengths))
    
    # Determine wavelength range
    if wavelength_range is None:
        min_wl = current_wavelengths.min()
        max_wl = current_wavelengths.max()
    else:
        min_wl, max_wl = wavelength_range
        # Ensure range is within original data
        min_wl = max(min_wl, current_wavelengths.min())
        max_wl = min(max_wl, current_wavelengths.max())
    
    # Create new wavelength grid
    new_wavelengths = np.arange(min_wl, max_wl + new_wavelength_step, new_wavelength_step)
    
    # Determine method automatically if not specified
    if method == 'auto':
        if new_wavelength_step > current_step:
            method = 'bin_average'  # Downsampling
        else:
            method = 'interpolate'  # Upsampling
    
    # Create output DataFrame
    resampled_data = {}
    
    if method == 'bin_average':
        # Binning for downsampling
        for col in spectrum_df.columns:
            spectrum_values = spectrum_df[col].values
            binned_values = []
            
            for new_wl in new_wavelengths:
                # Find wavelengths within the bin
                bin_half_width = new_wavelength_step / 2
                mask = (current_wavelengths >= new_wl - bin_half_width) & \
                       (current_wavelengths < new_wl + bin_half_width)
                
                if np.any(mask):
                    # Average values in the bin
                    binned_values.append(np.mean(spectrum_values[mask]))
                else:
                    # If no points in bin, use nearest neighbor
                    nearest_idx = np.argmin(np.abs(current_wavelengths - new_wl))
                    binned_values.append(spectrum_values[nearest_idx])
            
            resampled_data[col] = binned_values
    
    elif method == 'interpolate':
        # Nearest neighbor interpolation for upsampling
        try:
            from scipy.interpolate import interp1d
        except ImportError:
            raise ImportError("scipy is required for interpolation. Install with: pip install scipy")
        
        for col in spectrum_df.columns:
            spectrum_values = spectrum_df[col].values
            
            # Create interpolation function (nearest neighbor)
            interp_func = interp1d(current_wavelengths, spectrum_values, 
                                 kind='nearest', bounds_error=False, 
                                 fill_value='extrapolate')
            
            # Interpolate to new wavelength grid
            resampled_data[col] = interp_func(new_wavelengths)
    
    else:
        raise ValueError(f"Unknown method: {method}. Use 'auto', 'bin_average', or 'interpolate'")
    
    # Create output DataFrame
    resampled_df = pd.DataFrame(resampled_data, index=new_wavelengths)
    resampled_df.index.name = 'wavelength'
    
    print(f"Resampled spectrum from {current_step:.3f} nm to {new_wavelength_step:.3f} nm resolution")
    print(f"Method used: {method}")
    print(f"Original shape: {spectrum_df.shape}, New shape: {resampled_df.shape}")
    
    return resampled_df


def flux_to_photon_count(spectrum_df):
    import astropy.constants as const
    import astropy.units as u
    wavelengths = spectrum_df.index.values * u.nm
    flux_column = spectrum_df.columns[0]
    flux = spectrum_df[flux_column] * u.W / (u.m *2 * u.nm)
    photon_energy = (const.h * const.c / wavelengths).to(u.J)
    count_rate = (flux / photon_energy) #.to(1 / (u.m**2 * u.nm * u.s ))
    result = spectrum_df.copy()
    result['photon_energy'] = photon_energy.value
    result['count_rate'] = count_rate
    return result
    
def bin_with_mean_index(df, n_bins):
    """Bin dataframe using mean of original index values in each bin."""
    
    bins = pd.cut(df.index, bins=n_bins)
    
    # Group and calculate means for both data and index
    binned_data = df.groupby(bins).mean()
    binned_index = pd.Series(df.index).groupby(bins).mean()
    
    # Use the mean index values
    binned_data.index = binned_index.values
    binned_data.index.name = df.index.name
    
    return binned_data

def calculate_snr(target_spectrum,
                  background_spectrum,
                  integration_time=60,
                  slit_size=10,
                  collector_size=1,
                  optical_effeciency=0.1, 
                  quantum_efficiency=0.6, 
                  dark_current_rate=0.1, 
                  readout_noise=10,
                  full_well_capacity=1e6,
                  binning_factor=1, 
                  wavelengths_per_pixel=1
                  ):
    
    total_system_efficency = collector_size * optical_effeciency * quantum_efficiency * wavelengths_per_pixel
    signal_rate = target_spectrum["count_rate"] * total_system_efficency
    background_rate = background_spectrum["count_rate"] * slit_size * total_system_efficency
    signal = signal_rate * integration_time * binning_factor
    background = background_rate * integration_time * binning_factor
    dark_current = dark_current_rate * integration_time * binning_factor
    readout_noise *= binning_factor
    n_frames = np.ceil(np.max(signal + background + dark_current) / full_well_capacity)
    noise = np.sqrt(signal + background + dark_current  + n_frames * readout_noise **2 )
    df = target_spectrum.copy()
    df['snr'] = signal / noise
    return df, n_frames

def plot_spectrum(df, data_column, label='Spectrum', title='Spectrum', 
                  figsize=(7, 4), save_path='spectrum.png', ylabel='data', xlabel='Wavelength (nm)'):
    plt.figure(figsize=figsize)
    plt.plot(df.index, df[data_column], label=label)
    plt.title(title)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.plot()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    