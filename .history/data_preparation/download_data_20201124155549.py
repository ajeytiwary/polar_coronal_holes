import sys
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from sunpy.net import Fido
from sunpy.net import attrs as a
import glob
import os
import argparse
import time
from astropy.io.fits import getheader
'''
Download using fido fetch

'''


p = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

p.add_argument('-data_path', type=str, default='coronal_hole_heliohackweek/data/',
               help='Download Path')
p.add_argument('-startdate', type=str, default='2012/05/10 11:06',
               help='Start date and time in dd/mm/yyyy hhmm')
p.add_argument('-enddate', type=str, default='2012/06/10 12:06',
               help='end date and time in dd/mm/yyyy hhmm')
p.add_argument('-wavelengths', type=int, default=[170,304],
                help='wavelengths min and max')
p.add_argument('-cadence', type=int, default=8,
                help='sampling time for data')

args = p.parse_args()


print(r'---------------------- Fido download for EIT,EUVI-A, EUVI-B ----------------------')


eit = (a.Instrument("eit") &
          a.Time(args.startdate, args.enddate))

euvi = (a.vso.Source('STEREO_A') &
        a.Instrument("EUVI") &
        a.Time(args.startdate, args.enddate))
euvib = (a.vso.Source('STEREO_B') &
        a.Instrument("EUVI") &
        a.Time(args.startdate, args.enddate))



samp=a.Sample(args.cadence*u.minute)
wave=a.Wavelength(args.wavelengths[0]*u.AA, args.wavelengths[1]*u.AA)

print(args.startdate, args.enddate,wave,samp)
fethched_results = Fido.search(wave, eit | euvi | euvib, samp)
print(r'---------------------- Fido fetched files for EIT,EUVI-A, EUVI-B ----------------------')
print(fethched_results)

files = Fido.fetch(fethched_results,path=args.data_path)

while not files.errors:
    time.sleep(10)
    files=Fido.fetch(files)
    
'''
Renaming files now
'''
euvi_l= sorted(glob.glob(args.data_path+'*eua.fts'))
euvi_r = sorted(glob.glob(args.data_path+'*eub.fts'))
eit=sorted(glob.glob(args.data_path+'efz*'))

for file in eit:
    header_ = getheader(file)
    os.rename(file,args.data_path+'EIT_'+str(header_.get('date_obs'))+'_'+str(header_.get('wavelnth'))+'.fits')

for file in euvi_l:
    header_ = getheader(file)
    os.rename(file,args.data_path+'EUVI-A_'+str(header_.get('date-obs'))+'_'+str(header_.get('wavelnth'))+'.fits')
for file in euvi_r:
    header_ = getheader(file)
    os.rename(file,args.data_path+'EUVI-B_'+str(header_.get('date-obs'))+'_'+str(header_.get('wavelnth'))+'.fits')


print('All files downloaded successfully and files named appropriately renamed')