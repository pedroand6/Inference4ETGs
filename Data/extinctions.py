# Import packages and configure the folder to save the dust maps
from astropy.coordinates import SkyCoord
import extinction

import numpy as np
import pandas as pd

from dustmaps.config import config
config['data_dir'] = './dustMaps'

import dustmaps.csfd
dustmaps.csfd.fetch()

# To calculate extinctions
def correct_extinction(dataframe, extinction_maps, extinction_columns):
    '''
    Correct the magnitudes for extinction using the CCM89 Law

    Keyword arguments:
    dataframe         -- dataframe containing the data to be corrected
    extinction_maps   -- SFD Maps
    '''
    corrected_df = dataframe.copy().reset_index(drop=True)

    # Obtaining E(B-V) and Av in a given RA, DEC position
    input_file_coords = SkyCoord(dataframe['ra'], dataframe['dec'], frame='icrs', unit='deg')
    ebv = extinction_maps(input_file_coords)
    av  = 3.1*ebv

    # Calculating the extinction on the S-PLUS bands using the Cardelli, Clayton & Mathis law.
    Lambdas = np.array([3536, 4751, 6258, 7690, 8831]).astype(float)

    extinctions = []
    for i in range(len(av)):
        extinctions.append(extinction.ccm89(Lambdas, av[i], 3.1))

    extinction_df = pd.DataFrame(extinctions, columns=[extinction_columns])

    for i in range(len(extinction_columns)):
        corrected_df[extinction_columns[i]] = extinction_df[extinction_columns[i]]
        
    #corrected_df['EBV'] = ebv

    return corrected_df

# -- Correct data

# Filters
features_SPLUS = ['u_auto', 'g_auto', 'r_auto', 'i_auto', 'z_auto']
extinct_SPLUS = [filt+'_ext' for filt in features_SPLUS]

features_WISE = ['W1_ab', 'W2_ab']
extinct_WISE  = [filt+'_ext' for filt in features_WISE]

features_VHS = ['Ypmag', 'Jpmag', 'Hpmag', 'Kspmag']
extinct_VHS  = [filt+'_ext' for filt in features_VHS]

features_GALEX = ['FUVmag', 'NUVmag']
extinct_GALEX  = [filt+'_ext' for filt in features_GALEX]

features    = features_SPLUS
extinctions = extinct_SPLUS

# Load dust map
csfd_map = dustmaps.csfd.CSFDQuery(map_fname='dustMaps/csfd/csfd_ebv.fits', mask_fname='dustMaps/csfd/mask.fits')

# Apply corrections
dataframe = pd.read_csv('Data/morphgal.csv')
dataframe_corrected = correct_extinction(dataframe, csfd_map, extinctions)

for i in range(len(features)):
    dataframe_corrected[features[i]] = dataframe_corrected[features[i]] - dataframe_corrected[extinctions[i]]

# Save corrected data
dataframe_corrected.to_csv('Data/morphgal_corrected.csv', index=False)