'''
Loads downloaded data into Sunpy maps for each instrument and combine images to make 
synchronic map for each wavelength [195, 171, 280?] and a combined overlaid synchronic maps. 

For each date-time there are 4 synchronic maps available
'''


import scipy
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import sunpy.map
import sunpy.sun
from sunpy.map.maputils import all_coordinates_from_map
from sunpy.coordinates import get_horizons_coord
from reproject import reproject_interp
from reproject.mosaicking import reproject_and_coadd
import glob
import importlib
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS



# Masking function
def mask_outside_disk(inst_map):
    # Find coordinates and radius
    hpc_coords = all_coordinates_from_map(inst_map)
    r = np.sqrt(hpc_coords.Tx ** 2 + hpc_coords.Ty ** 2) / inst_map.rsun_obs

    # Mask everything outside of the solar disk
    mask = ma.masked_greater_equal(r, 1)
    ma.set_fill_value(mask, np.nan)
    where_disk = np.where(mask.mask == 1)

    return where_disk

# Make histogram
def make_hist(norm_log_inst, bins_inst):

    where_mask = np.isfinite(norm_log_inst)
    arr = norm_log_inst[where_mask]
    hist, bins = np.histogram(arr, bins=bins_inst, density=True)
    width = bins[1] - bins[0]

    return hist*width


# Find sample data
path_to_files = 'data/'
filenames_eit = sorted(glob.glob(path_to_files+'eit_l1*'))
filenames_euvil = sorted(glob.glob(path_to_files+'*eu_L.fts'))
filenames_euvir = sorted(glob.glob(path_to_files+'*eu_R.fts'))
filenames_test_set = sorted(glob.glob(path_to_files+'EUVI_mask*'))
# Number of samples
nb_files = len(filenames_eit)
