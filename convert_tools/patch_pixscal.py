'''
Created on Nov 13, 2020

@author: Colley JM

'''
from astropy.io import fits
import os.path
from astropy import wcs
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import pixel_to_skycoord

def add_kw_pixscal(file_fits):
    hdul = fits.open(file_fits)
    my_wcs = wcs.WCS(hdul[1].header)
    # get center
    center1 = int(hdul[1].header['NAXIS1']/2)
    center2 = int(hdul[1].header['NAXIS2']/2)
    # pointing of     
    sky_1 = pixel_to_skycoord(center1, center2, my_wcs, 0)
    sky_2 = pixel_to_skycoord(center1+1, center2, my_wcs, 0) 
    sky_3 = pixel_to_skycoord(center1, center2+1, my_wcs, 0)
    print(sky_1, sky_2)
    sep = sky_1.separation(sky_2)
    print(sep)
    print(sep.arcsec)
    print(sky_1.separation(sky_3).arcsec)
    hdul[1].header['PIXSCAL'] = (sep.arcsec, 'Pixel scale  (arcsec/pixel)') 
    hdul.writeto("patch.fits",overwrite=True)
    hdul.close()
     
     

if __name__ == '__main__':
    f_ima1 = "/home/user/Work/data/EUC_EXT_DPDEXTDETRENDEDFRAME_LSST-208257-R01-S01_20200929T091612.420Z_LSST-SWF1-Check.fits"
    add_kw_pixscal(f_ima1)