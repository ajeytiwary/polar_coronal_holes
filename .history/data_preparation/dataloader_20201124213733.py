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


def eit_correction(eitmap):
    new_coords = get_horizons_coord(eitmap.observatory.replace(' ', '-'),
                                    eitmap.date)
    eitmap.meta['HGLN_OBS'] = new_coords.lon.to('deg').value
    eitmap.meta['HGLT_OBS'] = new_coords.lat.to('deg').value
    eitmap.meta['DSUN_OBS'] = new_coords.radius.to('m').value
    eitmap.meta.pop('hec_x')
    eitmap.meta.pop('hec_y')
    eitmap.meta.pop('hec_z')
    return eitmap

'''
Function to combine three instrument into a synchronic map;
takes list of sunpy map for each instrument as input 
Output is a sunpy map object, a synchronic map combining maps for three instruments 
'''        
def combine_maps(maps_list):


    # Combined maps_list
    maps_list[0]=eit_correction(maps_list[0])
    where_mask = mask_outside_disk(maps_list[0])
    maps_list[0].data[where_mask] = np.nan
    where_mask = mask_outside_disk(maps_list[1])
    maps_list[1].data[where_mask] = np.nan
    where_mask = mask_outside_disk(maps_list[2])
    maps_list[2].data[where_mask] = np.nan
    shape_out = (180, 360)  # This is set deliberately low to reduce memory consumption
    header = sunpy.map.make_fitswcs_header(shape_out,
                                           SkyCoord(0, 0, unit=u.deg,
                                                    frame="heliographic_stonyhurst",
                                                    obstime=maps_list[0].date),
                                           scale=[180 / shape_out[0],
                                                  360 / shape_out[1]] * u.deg / u.pix,
                                           wavelength=int(maps_list[0].meta['wavelnth']) * u.AA,
                                           projection_code="CAR")
    out_wcs = WCS(header)
    coordinates = tuple(map(sunpy.map.all_coordinates_from_map, maps_list))
    weights = [coord.transform_to("heliocentric").z.value for coord in coordinates]
    weights = [(w / np.nanmax(w)) ** 3 for w in weights]
    for w in weights:
        w[np.isnan(w)] = 0

    array, _ = reproject_and_coadd(maps_list, out_wcs, shape_out,
                                   input_weights=weights,
                                   reproject_function=reproject_interp,
                                   match_background=True,
                                   background_reference=0)
    outmaps = sunpy.map.Map((array, header))
    return outmaps



# Nb. of channels
def satellite_position(map1,map2,map3):
    map1 
    if map1.observatory in ['SOHO']:
        long_eit =map1.observer_coordinate.lon.to('degree')

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
    return position_condition    


nb_channels = 3
list_channels = [171, 195, 304]

# Find sample data
path_to_files = '/home/ajt/Downloads/Hackweek_data/EUV_data/'




filenames_eit_195 = sorted(glob.glob(path_to_files+'EIT_'+'*195*' +'.fits'))
filenames_euvil_195 = sorted(glob.glob(path_to_files+'EUVI-A_'+'*195*' +'.fits'))
filenames_euvir_195 = sorted(glob.glob(path_to_files+'EUVI-B_'+'*195*' +'.fits'))


filenames_eit_171 = sorted(glob.glob(path_to_files+'EIT_'+'*171*' +'.fits'))
filenames_euvil_171 = sorted(glob.glob(path_to_files+'EUVI-A_'+'*171*' +'.fits'))
filenames_euvir_171 = sorted(glob.glob(path_to_files+'EUVI-B_'+'*171*' +'.fits'))


filenames_eit_304 = sorted(glob.glob(path_to_files+'EIT_'+'*304*' +'.fits'))
filenames_euvil_304 = sorted(glob.glob(path_to_files+'EUVI-A_'+'*304*' +'.fits'))
filenames_euvir_304 = sorted(glob.glob(path_to_files+'EUVI-B_'+'*304*' +'.fits'))

# Nb. of samples to work with
nb_files = max(len(filenames_euvil_171), len(filenames_euvir_171),len(filenames_eit_171))

# Hamada class for homogenization
filename = 'cumulative_hist.npz'
hamada_model = hamada_hist_matching.hamada(filename='cumulative_hist.npz')




'''
doing time matching now
'''


def get_time(file):
    '''returns time for a file in string format'''
    return Time(file.split('_')[-2][0:-1])

def check_tol(t1,t2,tolerance):
    if t1-t2 < timedelta(minutes=tolerance):
        return True
    else:
        return False

def check_time_neighbour(reference_file, file,  tolerance):
    '''returns True if within tolerance else returns False'''
    ref_t = get_time(reference_file)
    file_t = get_time(file)
    return check_tol(ref_t, file_t, tolerance)

def get_time_neighbour(reference_file, comparison_files, tolerance):
    '''Given a reference file and a list of comparison files,
    returns, if exists, the closest time neighbor'''
    for file in comparison_files:
        if check_time_neighbour(reference_file, file, tolerance):
            return file
            break
        else: 
            continue



observations_195 = []

for file in filenames_eit_195:
    observation_set = []
    observation_set.append(file)
    observation_set.append(get_time_neighbour(file, filenames_euvil_195, 4))
    observation_set.append(get_time_neighbour(file, filenames_euvir_195, 4))

    observations_195.append( observation_set )



observations_171 = []

for file in filenames_eit_171:
    observation_set = []
    observation_set.append(file)
    observation_set.append(get_time_neighbour(file, filenames_euvir_171, 4))
    observation_set.append(get_time_neighbour(file, filenames_euvil_171, 4))

    observations_171.append( observation_set )
    
observations_304 = []

for file in filenames_eit_304:
    observation_set = []
    observation_set.append(file)
    observation_set.append(get_time_neighbour(file, filenames_euvil_304, 4))
    observation_set.append(get_time_neighbour(file, filenames_euvir_304, 4))

    observations_304.append( observation_set )    
 

# Output




for files in observations_195:
    map_lists=sunpy.map.Map(files)
    outmap=combine_maps(map_lists)
    outmap.plot_settings = map_lists[0].plot_settings
    outmap.nickname = 'EIT + EUVI/A + EUVI/B'

    # Output
    outmap.save(path_to_files+files[0].split('_')[-2][0:-1]+'.fits', filetype='fits', overwrite=True) 