'''
Created on Mar 26, 2021

@author: user
'''
import matplotlib.pyplot as plt
import check_tools.catalog as ctc
from astropy.table import Table
from astropy.wcs import wcs
import numpy as np
import os
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
from astropy.table import Column
import astropy.io.fits as fits

S_path_butler = "/home/user/Work/data/pilot"


def convert_pos_pixel_2_radec(f_cat, f_frame):
    '''
    
    :param f_cat: LSST star/galaxie catalog 
    :param f_frame: LSST image with WCS
    
    return a astropy table of f_cat with 2 new columns "ra" and "dec"
    '''
    # retrieve Table of LSST catalog
    cat_lsst = ctc.get_cat(f_cat)
    print(cat_lsst)

    # define WCS converter, 
    
#  code example from
    hdul = fits.open(f_frame)
    header = hdul[1].header

#          
    my_wcs = wcs.WCS(header)
    print(my_wcs.wcs.name)

#     x = cat_lsst['XWIN_IMAGE']+0.5
#     y = cat_lsst['YWIN_IMAGE']+0.5
    x = cat_lsst['XWIN_IMAGE'] 
    y = cat_lsst['YWIN_IMAGE']

    # convert XWIN_IMAGE, YWIN_IMAGE in ra dec, in fact for star and galaxy
    # Need to SELECT converter based on WCS type in header
    
    # Core without distortions  
    ra, dec = my_wcs.wcs_pix2world(x, y, 1)
    my_wcs.printwcs()
    print("calc_footprint")
    print(my_wcs.calc_footprint(undistort=True))
    print(my_wcs.calc_footprint(undistort=False))    
    print("wcs.print_contents:\n")
    # my_wcs.wcs.print_contents()
    print("has distortion ? ", my_wcs.has_distortion)
    print("low_level_wcs:\n", my_wcs.low_level_wcs)

    # ra, dec from TPV with distortions
#     my_wcs_2 = wcs.WCS(header)
#     ra_dist, dec_dist = my_wcs_2.all_pix2world(x, y, 0)  # <- third argument is another method to transform with and without distortions

    # We select only star during match catalog
#  code example from  
# my_wcs.wcs_pix2world(0, 0, 0), 

    # add column ra dec to cat_lsst
#    cat_lsst.add_column(ra_from wcs..., name='ra')
#    cat_lsst.add_column(dec_from wcs..., name='dec')
   
    cat_lsst.add_column(Column(ra), name='RA')
    cat_lsst.add_column(Column(dec), name='DEC')
#     cat_lsst.add_column(Column(ra_dist), name='RA_DIST')
#     cat_lsst.add_column(Column(dec_dist), name='DEC_DIST')
   
    cat_lsst.write('test.fits', overwrite=True)
    
    return cat_lsst


def plot_histo(array, label_x="", title="", range=None):
    plt.figure()
    plt.title(title)
    if range:
        plt.hist(array, bins=57, range=range)
    else:
        plt.hist(array, bins=57)
    plt.xlabel(label_x)
    plt.grid()

    
class RefCatalogLSST(object):

    def __init__(self, path_butler):
        self.path_butler = path_butler
        self.filter = "u"
        self.l_filter = ["u", "g", "r", "i", "z", "y"]
        self.fact_sigma = 3
    
    def _set_cat_coord(self):
        self.cat_coord = SkyCoord(ra=self.cat_ref['RA'] * u.degree,
                                  dec=self.cat_ref['DEC'] * u.degree)

    def set_filter(self, p_filter):
        self.filter = p_filter
    
    @property
    def size(self):
        return len(self.cat_ref)
        
    def read_cat_ref(self, name_cat=None):
        if not name_cat:
            name_cat = "cat_temp.txt"
        self.name_cat = name_cat
        self.cat_ref = Table.read(os.path.join(self.path_butler, name_cat), format="ascii.commented_header", delimiter=',')
        print(self.cat_ref.info())
        
    def remove_high_magnitude(self, threasold=27):
        '''remove star with magnitude > threasold
        
        :param threasold:
        '''
        
        idx_high = np.where(self.cat_ref[self.filter] < threasold)
        print(f"Keep {len(idx_high[0])} on {self.size} stars with magnitude < {threasold} in filter {self.filter}")
        self.cat_ref = self.cat_ref[idx_high]
                
    def remove_out_cat(self, cat):
        '''remove star in ref catalog (for full focal plane) out of FOV of current image
        
        :param cat:
        '''
        min_p = np.amin(self.cat_ref["RA"])
        max_p = np.amax(self.cat_ref["RA"])
        print(f"RA ref: [{min_p:.7}, {max_p:.7}]")
        min_p = np.amin(self.cat_ref["DEC"])
        max_p = np.amax(self.cat_ref["DEC"])
        print(f"DEC ref: [{min_p:.7}, {max_p:.7}]")
        min_p = np.amin(cat["RA"])
        max_p = np.amax(cat["RA"])
        moy_ra = ((max_p - min_p) / 2 + min_p) * u.degree
        print(f"RA cat: [{min_p:.7}, {max_p:.7}]")
        min_p = np.amin(cat["DEC"])
        max_p = np.amax(cat["DEC"])
        print(f"DEC cat: [{min_p:.7}, {max_p:.7}]")
        moy_dec = ((max_p - min_p) / 2 + min_p) * u.degree
        pos = SkyCoord(ra=moy_ra, dec=moy_dec)
        a_radius = Angle(10, unit='arcmin')
        self._set_cat_coord()
        d2d = pos.separation(self.cat_coord)        
        catalogmsk = d2d < a_radius
        idx_near = np.where(catalogmsk)
        print(f"Keep {len(idx_near[0])} on {self.size} stars around {moy_ra:.6},{moy_dec:.6} with radius {a_radius.arcmin} arcsmin")        
        self.cat_ref = self.cat_ref[idx_near]
        print(self.size)
        self._set_cat_coord()            
        
    def search_around(self, radec, radius=1.0):
        '''search nearest star in cercle error position with lowest magnitude 
        
        :param radec:
        :param radius: in arcsec
        '''
        a_radius = Angle(radius * u.arcsec)
        print(a_radius.deg)
        if not hasattr(self, "cat_coord"):
            self._set_cat_coord()
        # idx_in = self.cat_coord.search_around_sky(radec, radius * u.arcsec)
        d2d = radec.separation(self.cat_coord)        
        catalogmsk = d2d < a_radius
        idx_near = np.where(catalogmsk)
        print("idx_near: ", idx_near)
        if len(idx_near[0]) == 0:
            # no match
            print("========> No matched !!!!!!!!!!!!!!")
#             fig, ax = plt.subplots()            
#             ax.plot(self.cat_ref['RA'], self.cat_ref['DEC'], "b.")
#             print(radec.ra.value, radec.dec.value)
#             ax.plot(radec.ra.value, radec.dec.value, "yx")
#             circle1 = plt.Circle((radec.ra.value, radec.dec.value),  a_radius.deg, color='y', fill=False)
#             ax.add_patch(circle1)
#             ax.grid()
            # plt.show()            
            return None, None      
        print(f"find {len(idx_near)} objects around with radius {a_radius}")
        self.idx_near = idx_near        
        # print(radec.ra.degree)
        idx_rel_min = np.argmin(self.cat_ref[idx_near][self.filter])        
        idx_min = idx_near[0][idx_rel_min]
        return idx_min, d2d[idx_min]
    
    def match_cat_ext(self, cat_ext):
        '''
        
        :param cat_ext:
        '''
        cpt_nm = 0
        l_sep = []
        l_dra = []
        l_ddec = []
        l_mag = []
        self.l_mag_no_match = []
        for star in cat_ext:
#             if len(l_sep) == 100:
#                 break
            if star["CLASS_STAR"] == 0.0:
                print(star)
                ra = star['RA']
                dec = star['DEC']
                # error ellipse => cercle => *3 ~ 3 sigma
                err_r = star['ERRAWIN_IMAGE']
                if err_r < star['ERRBWIN_IMAGE']:
                    err_r = star['ERRBWIN_IMAGE']
                # 0.2 arsec/pixel LSST
                err_r *= 0.2 * self.fact_sigma
                print("radius error [arcsec]:", err_r)
                pos = SkyCoord(ra=ra * u.degree, dec=dec * u.degree)                
                idx, sep = self.search_around(pos, err_r)
                if idx == None:
                    cpt_nm += 1
                    self.l_mag_no_match.append(star['MAG_PSF'])
                    continue
                print(idx, sep)
                l_sep.append(sep.arcsec)
                ddra = Angle(ra - self.cat_ref["RA"][idx], unit='degree').arcsec
                print("ddra:", ddra)
                l_dra.append(ddra)
                l_ddec.append(Angle(dec - self.cat_ref["DEC"][idx], unit='degree').arcsec)
                l_mag.append(star['MAG_PSF'])
        # create numpy file
        a_err = np.zeros((len(l_sep), 4), dtype=np.float64)
        print(a_err.shape)
        print(np.array(l_sep))
        print(l_sep, l_dra, l_ddec)
        a_err[:, 0] = np.array(l_sep)
        a_err[:, 1] = np.array(l_dra)
        a_err[:, 2] = np.array(l_ddec)
        a_err[:, 3] = np.array(l_mag)
        self.a_err = a_err
        np.save("astrometry_err", a_err)
        print(f"No match star {cpt_nm}")
    
    def load_astrometry_err(self, f_astro):
        self.a_err = np.load(f_astro)
    
    #
    # PLOT
    #
    def plot_histo_mag(self):
        plt.figure()
        plt.title(f'Histogram star magnitude in filter {self.filter}')
        plt.hist(self.cat_ref[self.filter], bins=30)
        
    def plot_histo_mag_all_filter(self):
        plt.figure()
        plt.title(f'Histogram star magnitude in all filters')
        trans = 0.2
        for pfil in self.l_filter:
            plt.hist(self.cat_ref[pfil], bins=30, alpha=trans)
        plt.legend(self.l_filter)
        plt.grid()
        
    def plot_astrometry_err(self):
        err_med = np.median(self.a_err[:, 0]) * 1000
        # title_err = r"Angle error for match source (<${}\sigma$ ERR_WIN)".format(self.fact_sigma)
        title_err = f"Angle error for match source (circle radius = {self.fact_sigma}*ERRWIN_IMAGE)"
        title_err += "\n" + f"Median error value {err_med:5.3} mas"
        title_err += "\n" + f"Total match {self.a_err.shape[0]} stars"
        plot_histo(self.a_err[:, 0] * 1000, "mas", title_err , (0, 200))
        plot_histo(self.a_err[:, 1] * 1000, "mas", "RA error for match source", (-1000, 1000))
        plot_histo(self.a_err[:, 2] * 1000, "mas", "DEC error for match source", (-1000, 1000))
        title_no_match = "Magnitude estimated for NO match stars"
        title_no_match += "\n" + f"Total NO match {len(self.l_mag_no_match)} stars"
        plot_histo(self.l_mag_no_match, "mag", title_no_match)
        plt.figure()
        plt.title("Histogram of star magnitudes from the LSST catalog")
        trans = 0.2
        plt.hist(self.a_err[:, 3], bins=30, alpha=trans)
        plt.hist(self.l_mag_no_match, bins=30, alpha=trans)
        msg_1 = f"Total matched stars : {self.a_err.shape[0]}"
        msg_2 = f"Total NO matched stars : {len(self.l_mag_no_match)}"
        plt.legend([msg_1, msg_2 ])
        plt.xlabel("mag")
        plt.grid()


class CheckAstrometryCatalog(object):
    
    def __init__(self, path_butler):
        self.path_butler = path_butler
        self.filter = "u"
        self.ref = RefCatalogLSST(path_butler)

    def set_catalog(self, f_cat, f_wcs):
        # return cat with ra, dec obtain with WCS
        self.cat = convert_pos_pixel_2_radec(f_cat, f_wcs)
        print(self.cat)
        self.cat.sort('FLUX_PSF')
        self.cat.reverse()
        
    def set_filter(self, filter):
        self.filter = filter
        self.ref.set_filter(filter)
    
    def plot_star_cat_on_image(self):
        pass
        
    def check(self):
        self.ref.read_cat_ref()
        self.ref.remove_high_magnitude()
        self.ref.remove_out_cat(self.cat)
        self.ref.match_cat_ext(self.cat)
        self.ref.plot_astrometry_err()

         
def test_RefCatalogLSST_remove():
    cat = RefCatalogLSST(S_path_butler)
    cat.read_cat_ref()
    cat.plot_histo_mag_all_filter()
    cat.plot_histo_mag()
    print(cat.size)
    print(cat.cat_ref)
    cat.remove_high_magnitude()
    print(cat.size)
    print(cat.cat_ref)
    cat.plot_histo_mag()


def test_RefCatalogLSST_search():
    cat = RefCatalogLSST(S_path_butler)
    cat.read_cat_ref()
    cat.remove_high_magnitude()
    ra = 228.651726523 * u.degree
    dec = 30.7653586443 * u.degree
    pos = SkyCoord(ra=ra, dec=dec)
    print(pos)
    print(cat.cat_ref)
    print(cat.search_around(pos, 60))


def test_test_RefCatalogLSST_match_cat_ext():
    cat = RefCatalogLSST(S_path_butler)
    cat.read_cat_ref()
    cat.remove_high_magnitude()
    ext_cat = cat.cat_ref.copy()
    l_selec = [10, 1000, 3000]
    ext_cat = ext_cat[l_selec]
    xerr = np.ones(3, dtype=np.float64)
    print(xerr)
    zeros = np.zeros(3)
    ext_cat.add_column(Column(xerr), name="ERRAWIN_IMAGE")
    ext_cat.add_column(Column(xerr), name="ERRBWIN_IMAGE")
    ext_cat.add_column(Column(zeros), name="CLASS_STAR")
    offset = Angle(1, unit='arcsec') / 10
    ext_cat['RA'] += offset.degree
    ext_cat['DEC'] += 2 * offset.degree
    cat.match_cat_ext(ext_cat)
    cat.plot_astrometry_err()
    print(ext_cat)
    

def test_distorsion():
    path = "/home/user/Work/data/pilot/"
    frame = path + "EUC_EXT_DPDEXTDETRENDEDFRAME_LSST-140219-R32-S12_20210312T120709.895Z_SC8_PF_RUBIN_79171_R1_ST1a.fits"
    cat = path + "EUC_EXT_DPDEXTSOURCECATALOG_LSST-140219-R32-S12_20210312T120759.471Z_SC8_PF_RUBIN_79171_R1_ST1a.fits"
#     path = "/home/user/Desktop/"
#     frame=path+"EUC_EXT_DPDEXTDETRENDEDFRAME_DECAM-90001080-01_20190111T165610.3Z_00.00.fits.fz"
#     frame=path+"tpv.fits"
    cat_rd = convert_pos_pixel_2_radec(cat, frame)
    dra = Angle(cat_rd["RA"] - cat_rd["RA_DIST"], unit="degree").arcsec
    print("max diff :", np.amax(dra), np.amin(dra))
    plot_histo(dra)
    

def check_catalog_stack(butler, f_cat, f_wcs, filter):
    astro = CheckAstrometryCatalog(butler)
    astro.set_catalog(f_cat, f_wcs)
    astro.set_filter(filter)
    astro.check()

 
def test_CheckFocalPlane():
    check = CheckFocalPlane()
    check.prepare_cat_ref()
    check.set_catalog(f_cat, f_wcs)

    
if __name__ == '__main__':
    path = "/home/user/Work/data/pilot/"
    frame = path + "EUC_EXT_DPDEXTDETRENDEDFRAME_LSST-140219-R32-S12_20210319T141211.357Z_SC8_PP_RUBIN_v2.fits"
    cat = path + "EUC_EXT_DPDEXTSOURCECATALOG_LSST-140219-R32-S12_20210319T141246.101Z_SC8_PP_RUBIN_v2.fits"
    path = "/home/user/Work/data/pilot/2127698/"
    frame = path + "EUC_EXT_DPDEXTDETRENDEDFRAME_LSST-2127698-R42-S22_20210409T132234.800Z_SC8_PP_RUBIN_v3.fits.fz"
    cat = path + "EUC_EXT_DPDEXTSOURCECATALOG_LSST-2127698-R42-S22_20210409T132305.271Z_SC8_PP_RUBIN_v3.fits"
    
#     path = "/home/user/Work/data/pilot/2127698/"
#     frame = path + "EUC_EXT_DPDEXTDETRENDEDFRAME_LSST-2127698-R42-S22_20210312T021052.926Z_SC8_PF_RUBIN_79171_R1_ST1a.fits"
#     cat = path + "EUC_EXT_DPDEXTSOURCECATALOG_LSST-2127698-R42-S22_20210312T021130.385Z_SC8_PF_RUBIN_79171_R1_ST1a.fits"
    
    # test_CheckFocalPlane()
#    test_RefCatalogLSST_remove()
    # test_RefCatalogLSST_search()
    # test_test_RefCatalogLSST_match_cat_ext()
    check_catalog_stack(path, cat, frame, "i")
    # test_distorsion()
    plt.show()
    pass
