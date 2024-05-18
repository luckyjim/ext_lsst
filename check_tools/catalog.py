'''
Created on Sep 29, 2020

@author: Colley Jm
'''
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import patches
import matplotlib.colors as colors

import astropy.io.fits as fits
import astropy.visualization as avis
from astropy.table import Table
import astropy.coordinates as coord
import copy

f_cat1 = "/home/user/Work/data/EUC_EXT_DPDEXTSOURCECATALOG_LSST-208257-R01-S01_20200929T091551.321Z_LSST-SWF1-Check.fits"
f_cat2 = "/home/user/Work/data/EUC_EXT_DPDEXTSOURCECATALOG_LSST-208257-R01-S01_20200929T091636.459Z_LSST-SWF1-Check.fits"
f_ima1 = "/home/user/Work/data/EUC_EXT_DPDEXTDETRENDEDFRAME_LSST-208257-R01-S01_20200929T091612.420Z_LSST-SWF1-Check.fits"


def get_cat(f_cat):
    with fits.open(f_cat) as hdul:
        hdul.info()
        cat = hdul[2].data
    return Table(cat)


def get_image_raw(f_ima):
    with fits.open(f_ima) as hdul:
        hdul.info()
        ima = hdul[1].data
    return ima


class MatchCatalogInCCD(object):

    def __init__(self):
        self.dkw = {'x': "XWIN_IMAGE",
                       'y':  "YWIN_IMAGE",
                       'ella' : "ERRAWIN_IMAGE",
                       'ellb' : "ERRBWIN_IMAGE",
                       'flux' : "FLUX_PSF"
                       }
        
    def set_dkw(self, dict_keyword):
        self.dkw = dict_keyword
    
    def set_file_cat_ref_match(self, fcat_ref, fcat_match):
        self.cat_ref = get_cat(fcat_ref)
        self.cat_match = get_cat(fcat_match)
        
    def set_cat_ref_match(self, cat_ref, cat_match):
        self.cat_ref = cat_ref
        self.cat_match = cat_match
    
    def filter_ref_flux(self, flux_min=None, nb_src=None):
        self.cat_ref.sort(self.dkw["flux"], reverse=True)
        if flux_min:
            idx = np.where(self.cat_ref[self.dkw["flux"]] > flux_min)[0]
            self.cat_ref = self.cat_ref[idx]
        elif nb_src:
            self.cat_ref = self.cat_ref[:nb_src]
            
    @property
    def len_cat_ref(self):
        return len(self.cat_ref)
    
    @property
    def len_cat_match(self):
        return len(self.cat_match)

    def do_match(self):
        v_zero = np.zeros(self.len_cat_ref)
        pos_ref = coord.SkyCoord(x=self.cat_ref[self.dkw['x']], y=self.cat_ref[self.dkw['y']], z=v_zero, representation_type='cartesian')
        v2_zero = np.zeros(self.len_cat_match)
        pos_mat = coord.SkyCoord(x=self.cat_match[self.dkw['x']], y=self.cat_match[self.dkw['y']], z=v2_zero, representation_type='cartesian')
        l_idx, sep, dist = coord.matching.match_coordinates_3d(pos_mat, pos_ref)
        self.idx_ref_associated = []
        self.idx_matched = []
        lrank = range(len(l_idx))
        err_match = (self.cat_match[self.dkw['ella']] + self.cat_match[self.dkw['ellb']]) / 2
        err_ref = (self.cat_ref[self.dkw['ella']] + self.cat_ref[self.dkw['ellb']]) / 2        
        for idx_ref, rk in zip(l_idx, lrank):                       
            err = (err_match[rk] + err_ref[idx_ref])/6
            if dist[rk] < err:            
                print("Position: ", pos_ref[idx_ref], pos_mat[rk])
                print('Flux: ', self.cat_ref[self.dkw['flux']][idx_ref])
                self.idx_ref_associated.append(idx_ref)
                self.idx_matched.append(rk)
        print(f"Try to match ref cat with {self.len_cat_ref} sources with catalogue with {self.len_cat_match} sources")
        print(f"{len(self.idx_matched)} sources matched with ref catalog.")
        if len( self.idx_ref_associated) <  self.len_cat_ref:
            self.status_match = "under"
        elif len( self.idx_ref_associated) == self.len_cat_ref:
            self.status_match = "same"
        else:
            self.status_match = "over"
    
    def compare_flux(self):
        if self.status_match == "same" or self.status_match == "under":
            plt.figure()
            flux_ref = self.cat_ref[self.dkw['flux']][self.idx_ref_associated]
            err_rel = 100*np.abs((flux_ref-self.cat_match[self.dkw['flux']][self.idx_matched])/flux_ref)
            plt.subplot(122)
            plt.title('relative error of estimated flux on matched sources')
            plt.hist(err_rel)
            plt.xlabel('%')
            plt.ylabel('source number')
            plt.grid(True)
            plt.subplot(121)
            plt.title('histogram flux sources')
            plt.hist(flux_ref, orientation='horizontal')
            plt.xlabel('source number')
            plt.ylabel('flux [?]')
            plt.grid(True)
            
    def compare_flux_raw(self):
        flux_ref = self.cat_match[self.dkw['flux']][self.idx_matched]
        err_rel = 100*np.abs((flux_ref-self.cat_ref[self.dkw['flux']][self.idx_ref_associated])/flux_ref)
        plt.figure()
        plt.subplot(122)
        plt.title('relative error of estimated flux on matched sources')
        plt.hist(err_rel,bins=30)
        plt.xlabel('%')
        plt.ylabel('source number')
        plt.grid(True)
        plt.subplot(121)
        plt.title('histogram flux sources')
        plt.hist(flux_ref, orientation='horizontal', log=True)
        plt.xlabel('source number')
        plt.ylabel('flux [?]')
        plt.grid(True)
            
    def plot_cat(self, name='ref', f_image=None):
        if name == "ref":
            cat = self.cat_ref
        else:
            cat = self.cat_match
        plot_catalog_in_CCD(cat, name, f_image)

    def plot_matched_cat(self, f_image=None):
        plot_two_catalog_in_CCD(self.cat_ref, self.cat_match[self.idx_matched], "Match catalog result", f_image, 10)
        

def test_astromatch():
    np.random.seed(10)
    pos = np.random.uniform(0, 100, 30).reshape((10, 3))
    pos[:, 2] = 0
    print(pos)
    # Z is required
    cat1 = coord.SkyCoord(x=pos[:, 0], y=pos[:, 1], z=pos[:, 2], representation_type='cartesian')
    print(cat1)
    noise = np.random.uniform(0, 1, 20).reshape((10, 2))
    pos2 = pos.copy()
    pos2 = np.flip(pos2, 0)
    pos2[:, 0] = pos2[:, 0] + noise[:, 0]
    pos2[:, 1] = pos2[:, 1] + noise[:, 1]
    print(pos2)
    cat2 = coord.SkyCoord(x=pos2[:, 0], y=pos2[:, 1], z=pos2[:, 2], representation_type='cartesian')
    print(cat2)
    res = coord.matching.match_coordinates_3d(cat1, cat2)
    print(res)
    pos2[7, :] = 0
    cat2 = coord.SkyCoord(x=pos2[4:, 0], y=pos2[4:, 1], z=pos2[4:, 2], representation_type='cartesian')
    print(cat2)
    l_idx, sep, dist = coord.matching.match_coordinates_3d(cat2, cat1)
    print(res)
    lrank = range(len(l_idx))
    for idx, rk in zip(l_idx, lrank):
        print(f'\nMatch candidate with dist {dist[rk]}')
        print(cat1[idx], cat2[rk])
        if dist[rk] > 1:
            print('REJECTED')
        else:
            print('OK !!!')
    
    
def plot_catalog_in_CCD_old(cat, title='', f_image=None):
    plt.figure()
    plt.title(title)
    # get current axis
    if f_image:
        ima = get_image_raw(f_image)
        # plt.imshow(ima, norm=colors.LogNorm(vmin=ima.min(), vmax=ima.max()), cmap='PuBu_r')
        m_norm = avis.ImageNormalize(ima, avis.ZScaleInterval(), stretch=avis.SinhStretch())
        plt.imshow(ima, norm=m_norm, cmap='gray')    
        plt.colorbar(extend='max')
    ax = plt.gca()
    for obj in cat:
        x_pos = obj["XWIN_IMAGE"]
        y_pos = obj["YWIN_IMAGE"]   
        diam_x = 2 * obj["ERRAWIN_IMAGE"] 
        diam_y = 2 * obj["ERRBWIN_IMAGE"]
        angle = obj["ERRTHETAWIN_IMAGE"]
        ell = patches.Ellipse((x_pos, y_pos), diam_x, diam_y, angle , edgecolor='yellow', alpha=0.8, lw=3)
        ax.add_patch(ell)
        ax.scatter(x_pos, y_pos, c='blue', s=3)        


def plot_catalog_in_CCD(cat, title='', f_image=None):
    plt.figure()
    plt.title(title)
    # get current axis
    if f_image:
        ima = get_image_raw(f_image)
        # plt.imshow(ima, norm=colors.LogNorm(vmin=ima.min(), vmax=ima.max()), cmap='PuBu_r')
        m_norm = avis.ImageNormalize(ima, avis.ZScaleInterval(), stretch=avis.SinhStretch())
        plt.imshow(ima, norm=m_norm, cmap='gray')    
        plt.colorbar(extend='max')
    ax = plt.gca()
    for obj in cat:
        x_pos = obj["XWIN_IMAGE"]
        y_pos = obj["YWIN_IMAGE"]   
        diam_x = 2 * obj["ERRAWIN_IMAGE"] 
        diam_y = 2 * obj["ERRBWIN_IMAGE"]
        angle = obj["ERRTHETAWIN_IMAGE"]
        ell = patches.Ellipse((x_pos, y_pos), diam_x, diam_y, angle , edgecolor='yellow', alpha=0.05, lw=3)
        ax.add_patch(ell)
        ax.scatter(x_pos, y_pos, c='blue', s=3)        

def plot_CCD_catalog_in_CCD(cat, title='', f_image=None):
    plt.figure()   
    plt.subplot(121)
    # get current axis
    color_src = ['yellow', 'red']
    plt.title(title)
    if f_image:
        ima = get_image_raw(f_image)
        # plt.imshow(ima, norm=colors.LogNorm(vmin=ima.min(), vmax=ima.max()), cmap='PuBu_r')
        m_norm = avis.ImageNormalize(ima, avis.ZScaleInterval(), stretch=avis.SinhStretch())
        plt.imshow(ima, norm=m_norm, cmap='gray')    
        plt.colorbar(extend='max')
    plt.subplot(122)
    plt.title(title + "+ stack LSST catalog")
    ax = plt.gca()
    if f_image:        
        # plt.imshow(ima, norm=colors.LogNorm(vmin=ima.min(), vmax=ima.max()), cmap='PuBu_r')
        m_norm = avis.ImageNormalize(ima, avis.ZScaleInterval(), stretch=avis.SinhStretch())
        plt.imshow(ima, norm=m_norm, cmap='gray')    
        plt.colorbar(extend='max')
    idx_src_sel = 0 
    for obj in cat:
        flux = obj["FLUX_PSF"]
        s_flux = f'FLUX_PSF: {flux}'
        if  flux > obj["FLUXERR_PSF"]*6:
            idx_src_sel += 1            
            x_pos = obj["XWIN_IMAGE"]
            y_pos = obj["YWIN_IMAGE"] 
            diam_x = 2 * obj["ERRAWIN_IMAGE"] 
            diam_y = 2 * obj["ERRBWIN_IMAGE"]
            angle = obj["ERRTHETAWIN_IMAGE"]
            idx_col = int(obj["CLASS_STAR"]+0.1)
            ell = patches.Ellipse((x_pos, y_pos), diam_x, diam_y, angle , fill=False, edgecolor=color_src[idx_col], fc=None, lw=2)
            ax.add_patch(ell)
            #ax.scatter(x_pos, y_pos, c='blue', s=3)
    print(f'Select {idx_src_sel}/{len(cat) }')
        
                
def plot_two_catalog_in_CCD_old(cat1, cat2, title='', f_image=None, factor_ell=1):
    plt.figure()
    plt.title(title)
    # get current axis
    ax = plt.gca()
    l_cat = [cat1, cat2]
    l_col = ['blue', 'yellow']
    l_ls = ['--', '-.']
    print(f_image)
    if f_image:
        ima = get_image_raw(f_image)
        # plt.imshow(ima, norm=colors.LogNorm(vmin=ima.min(), vmax=ima.max()), cmap='PuBu_r')
        m_norm = avis.ImageNormalize(ima, avis.ZScaleInterval(), stretch=avis.SinhStretch())
        plt.imshow(ima, norm=m_norm, cmap='gray')    
        plt.colorbar(extend='max')
        print('plot image')
    for cat, col, m_ls in zip(l_cat, l_col, l_ls):
        for obj in cat:
            x_pos = obj["XWIN_IMAGE"]
            y_pos = obj["YWIN_IMAGE"]   
            diam_x = 2 * obj["ERRAWIN_IMAGE"] * factor_ell
            diam_y = 2 * obj["ERRBWIN_IMAGE"] * factor_ell
            angle = obj["ERRTHETAWIN_IMAGE"]
            ell = patches.Ellipse((x_pos, y_pos), diam_x, diam_y, angle , edgecolor=col, fill=False, lw=2, ls=m_ls)
            ax.add_patch(ell)
            ax.scatter(x_pos, y_pos, c=col, s=3)

        
def load_cat_stage_1_and_plot(f_cat):    
    plot_catalog_in_CCD(get_cat(f_cat), f_cat)


def load_2_cat_stage_1_and_plot(f_cat_1, f_cat_2, f_ima_1=None):    
    cat_1 = get_cat(f_cat_1)
    cat_2 = get_cat(f_cat_2)
    plot_two_catalog_in_CCD(cat_1, cat_2, f_image=f_ima_1)
    

def plot_image(f_ima_1):
    ima = get_image_raw(f_ima_1)
    plt.figure()
    print(ima.shape)
    ax = plt.gca()
    plt.imshow(ima, norm=colors.LogNorm(vmin=ima.min(), vmax=ima.max()), cmap='PuBu_r')    
    plt.colorbar(extend='max')

    
def plot_image_scale(f_ima_1):
    ima = get_image_raw(f_ima_1)
    plt.figure()
    print(ima.shape)
    ax = plt.gca()
    m_norm = avis.ImageNormalize(ima, avis.ZScaleInterval(), stretch=avis.SinhStretch())
    plt.imshow(ima, norm=m_norm, cmap='PuBu_r')    
    plt.colorbar(extend='max')


def compare_cat_flux(fcat1, fcat2, fima):
    nb_src = 30
    cat_1 = get_cat(fcat1)
    cat_2 = get_cat(fcat2)
    # sort by flux
    print(type(cat_1))
    cat_1.sort('FLUX_PSF', reverse=True)
    cat_2.sort('FLUX_PSF', reverse=True)
    cat_1 = cat_1[:nb_src]
    print("cat_1\n")
    cat_1.pprint_all()
    cat_2 = cat_2[:nb_src]
    print("cat_2\n")
    cat_2.pprint_all()
    plot_two_catalog_in_CCD(cat_1, cat_2, fima, factor_ell=10)


def test_MatchCatalogInCCD_len():
    o_match = MatchCatalogInCCD()
    o_match.set_file_cat_ref_match(f_cat1, f_cat2)
    print(o_match.len_cat_ref)
    o_match.filter_ref_flux(nb_src=10)
    print(o_match.len_cat_ref)


def test_MatchCatalogInCCD_filter():
    o_match = MatchCatalogInCCD()
    o_match.set_file_cat_ref_match(f_cat1, f_cat2)
    print(o_match.len_cat_ref)
    o_match.filter_ref_flux(10000)
    print(o_match.len_cat_ref)
    o_match.plot_cat(f_image=f_ima1)


def test_MatchCatalogInCCD_match():
    o_match = MatchCatalogInCCD()
    o_match.set_file_cat_ref_match(f_cat1, f_cat2)
    print(o_match.len_cat_ref)
#    o_match.filter_ref_flux(nb_src=50)
    o_match.filter_ref_flux(flux_min=600_000)
    o_match.do_match()
    o_match.plot_matched_cat(f_ima1)
    o_match.compare_flux()
    o_match.compare_flux_raw()

    
if __name__ == '__main__':
# plot_image_scale(f_ima1)
#     plot_image(f_ima_1)
#     load_cat_stage_1_and_plot(f_cat_1)
#     load_2_cat_stage_1_and_plot(f_cat1, f_cat2, f_ima1)
#     compare_cat_flux(f_cat1, f_cat2, f_ima1)
    # test_astromatch()
    #test_MatchCatalogInCCD_filter()
#     test_MatchCatalogInCCD_len()
    test_MatchCatalogInCCD_match()
    plt.show()
