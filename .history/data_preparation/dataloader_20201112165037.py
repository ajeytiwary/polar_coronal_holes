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


# Find downloaded data
path_to_files = 'data/'
filenames_eit = sorted(glob.glob(path_to_files+'EIT_*'))
filenames_euvil = sorted(glob.glob(path_to_files+'*eua.fts'))
filenames_euvir = sorted(glob.glob(path_to_files+'*eub.fts'))
#filenames_test_set = sorted(glob.glob(path_to_files+'EUVI_mask*'))
# Number of samples
nb_files = len(filenames_eit)

def wavelet_enhancement(img):
    
    # incomplete
    img_enhanced = img
    
    return img_enhanced

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


def make_hist(norm_log_inst, bins_inst):
    
    where_mask = np.isfinite(norm_log_inst)
    arr = norm_log_inst[where_mask]
    hist, bins = np.histogram(arr, bins=bins_inst, density=True)
    width = bins[1] - bins[0]
    
    return hist*width

'''
A quick fix for histigram matching
'''

def match_cumulative_cdf(source, template):
    """
    Return modified source array so that the cumulative density function of
    its values matches the cumulative density function of the template.
    """
    src_values, src_unique_indices, src_counts = np.unique(source.ravel(),
                                                           return_inverse=True,
                                                           return_counts=True)
    tmpl_values, tmpl_counts = np.unique(template.ravel(), return_counts=True)

    # calculate normalized quantiles for each array
    src_quantiles = np.cumsum(src_counts) / source.size
    tmpl_quantiles = np.cumsum(tmpl_counts) / template.size

    interp_a_values = np.interp(src_quantiles, tmpl_quantiles, tmpl_values)
    return interp_a_values[src_unique_indices].reshape(source.shape)

'''
Function to combine three instrument into a synchronic map;
takes list of sunpy map for each instrument as input 
Output is a sunpy map object, a synchronic map combining maps for three instruments 
'''        
def combine_maps(maps_list):


    # Combined maps
    maps=maps_list

    # Combined maps
    shape_out = (180, 360)  # This is set deliberately low to reduce memory consumption
    header = sunpy.map.make_fitswcs_header(shape_out,
                                           SkyCoord(0, 0, unit=u.deg,
                                                    frame="heliographic_stonyhurst",
                                                    obstime=maps[0].date),
                                           scale=[180 / shape_out[0],
                                                  360 / shape_out[1]] * u.deg / u.pix,
                                           wavelength=int(maps[0].meta['wavelnth']) * u.AA,
                                           projection_code="CAR")
    out_wcs = WCS(header)
    coordinates = tuple(map(sunpy.map.all_coordinates_from_map, maps))
    weights = [coord.transform_to("heliocentric").z.value for coord in coordinates]
    weights = [(w / np.nanmax(w)) ** 3 for w in weights]
    for w in weights:
        w[np.isnan(w)] = 0

    array, _ = reproject_and_coadd(maps, out_wcs, shape_out,
                                   input_weights=weights,
                                   reproject_function=reproject_interp,
                                   match_background=True,
                                   background_reference=0)
    outmap = sunpy.map.Map((array, header))


# Nb. of channels
nb_channels = 3
list_channels = [171, 195, 304]

# Find sample data
path_to_files = '/home/ajt/Downloads/Hackweek_data/EUV_data/'




filenames_eit_195 = sorted(glob.glob(path_to_files+'EIT_'+'*195*' +'.fits'))
filenames_euvil_195 = sorted(glob.glob(path_to_files+'EIT_'+'*195*' +'.fits'))
filenames_euvir_195 = sorted(glob.glob(path_to_files+'EIT_'+'*195*' +'.fits'))


filenames_eit_171 = sorted(glob.glob(path_to_files+'EIT_'+'*171*' +'.fits'))
filenames_euvil_171 = sorted(glob.glob(path_to_files+'EIT_'+'*171*' +'.fits'))
filenames_euvir_171 = sorted(glob.glob(path_to_files+'EIT_'+'*171*' +'.fits'))


filenames_eit_304 = sorted(glob.glob(path_to_files+'EIT_'+'*304*' +'.fits'))
filenames_euvil_304 = sorted(glob.glob(path_to_files+'EIT_'+'*304*' +'.fits'))
filenames_euvir_304 = sorted(glob.glob(path_to_files+'EIT_'+'*304*' +'.fits'))

# Nb. of samples to work with
nb_files = max(len(filenames_euvil_171), len(filenames_euvir_171),len(filenames_eit_171))

# Hamada class for homogenization
filename = 'cumulative_hist.npz'
hamada_model = hamada_hist_matching.hamada(filename='cumulative_hist.npz')





np.datetime64('2015-07-04 12:59:59.50')



'''

STOP UNCHECKED AFTER THIS LINE

'''
# Output
output_path = '/home/plasmion/'




eit_maps = sunpy.map.Map(filenames_eit)
nx_eit, ny_eit = eit_maps[0].data.shape

euvil_maps = sunpy.map.Map(filenames_euvil)
nx_euvil, ny_euvil = euvil_maps[0].data.shape

euvir_maps = sunpy.map.Map(filenames_euvir)
nx_euvir, ny_euvir = euvir_maps[0].data.shape


for file_nb in range(nb_files):
    
    # Make map objects for one channel
    eit_maps = sunpy.map.Map(filenames_eit_0[file_nb])
    euvil_maps = sunpy.map.Map(filenames_euvil_0[file_nb])
    euvir_maps = sunpy.map.Map(filenames_euvir_0[file_nb])
    filename_extract = filenames_eit_0[file_nb].split(path_to_files+str(list_channels[0]))
    
    filename_output = output_path + 'composite_' + filename_extract[1]
    print(filename_extract)
    print('here')
    if eit_maps.observatory in ['SOHO']:
        new_coords = get_horizons_coord(eit_maps.observatory.replace(' ', '-'),
                                        eit_maps.date)
        eit_maps.meta['HGLN_OBS'] = new_coords.lon.to('deg').value
        eit_maps.meta['HGLT_OBS'] = new_coords.lat.to('deg').value
        eit_maps.meta['DSUN_OBS'] = new_coords.radius.to('m').value
        eit_maps.meta.pop('hec_x')
        eit_maps.meta.pop('hec_y')
        eit_maps.meta.pop('hec_z')
    
    # Check positioning of instruments in order to cover full Sun
    long_eit = eit_maps.observer_coordinate.lon.to('degree')
    print(long_eit.value)
    long_euvil = euvil_maps.observer_coordinate.lon.to('degree')-long_eit
    print(long_euvil.value)
    long_euvir = euvir_maps.observer_coordinate.lon.to('degree')-long_eit
    print(long_euvir.value)
    if -90 < long_euvil.value < 90 and -90 < long_euvir.value < 90:
        position_condition = 0
    elif 90 < long_euvil.value < 180 and long_euvil.value-180 < long_euvir.value < 180:
        position_condition = 0
    elif 90 < long_euvir.value < 180 and long_euvir.value-180 < long_euvil.value < 180:
        position_condition = 0
    elif -180 < long_euvil.value < -90 and -180 < long_euvir.value < long_euvil.value+180:
        position_condition = 0
    elif -180 < long_euvir.value < -90 and -180 < long_euvil.value < long_euvir.value+180:
        position_condition = 0
    else:
        position_condition = 1
    
    # If positioning is ideal, continue # MODIFY
    if position_condition == 1:
        
        # continue
        for channel_nb in range(nb_channels):
            # Make map objects
            if channel_nb == 1:
                eit_maps = sunpy.map.Map(filenames_eit_1[file_nb])
                euvil_maps = sunpy.map.Map(filenames_euvil_1[file_nb])
                euvir_maps = sunpy.map.Map(filenames_euvir_1[file_nb])
                filename_extract = filenames_eit_1[file_nb].split(path_to_files+str(list_channels[1]))
                filename_output = output_path + 'composite_' + filename_extract[1]
            elif channel_nb == 2:
                eit_maps = sunpy.map.Map(filenames_eit_2[file_nb])
                euvil_maps = sunpy.map.Map(filenames_euvil_2[file_nb])
                euvir_maps = sunpy.map.Map(filenames_euvir_2[file_nb])
                filename_output = output_path+'composite_' + filenames_eit_2[file_nb]
                filename_extract = filenames_eit_2[file_nb].split(path_to_files+str(list_channels[2]))
                filename_output = output_path + 'composite_' + filename_extract[1]
                
            nx_eit, ny_eit = eit_maps.data.shape
            nx_euvil, ny_euvil = euvil_maps.data.shape
            nx_euvir, ny_euvir = euvir_maps.data.shape
            
            # Adjust SoHO if not already done
            if channel_nb > 0 and eit_maps.observatory in ['SOHO']:
                new_coords = get_horizons_coord(eit_maps.observatory.replace(' ', '-'),
                                                eit_maps.date)
                eit_maps.meta['HGLN_OBS'] = new_coords.lon.to('deg').value
                eit_maps.meta['HGLT_OBS'] = new_coords.lat.to('deg').value
                eit_maps.meta['DSUN_OBS'] = new_coords.radius.to('m').value
                eit_maps.meta.pop('hec_x')
                eit_maps.meta.pop('hec_y')
                eit_maps.meta.pop('hec_z')

            # Mask everything outside the solar disk
            where_mask = mask_outside_disk(eit_maps)
            eit_maps.data[where_mask] = np.nan
            where_mask = mask_outside_disk(euvil_maps)
            euvil_maps.data[where_mask] = np.nan
            where_mask = mask_outside_disk(euvir_maps)
            euvir_maps.data[where_mask] = np.nan
            maps_list=[eit_maps,euvil_maps,euvir_maps]
            # Wavelet enchancement of the EIT data for improved contrast
            # Missing Step ######################
            # arr_tmp = eit_maps.data
            # eit_maps.data = wavelet_enhancement(arr_tmp)
            # Homogenization of the EIT data /w respect to EUVI
            outmap = combine_maps(maps_list)
            outmap.plot_settings = maps_list[0].plot_settings
            outmap.nickname = 'EIT + EUVI/A + EUVI/B'
        
            # Output
            outmap.save(filename_output, filetype='fits', overwrite=True)
