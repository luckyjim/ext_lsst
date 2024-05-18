import os
import numpy as np
from astropy.io import fits
import lsst.daf.persistence
from sip2tpv import sip2tpv
import pytz
import datetime
from collections import namedtuple
from contextlib import contextmanager
import tempfile
from astropy.io.fits.column import FITS2NUMPY
import astropy.table as at


PsfInfo = namedtuple("PsfInfo", "data header shape")
G_verbose = 0

G_time_0 = datetime.datetime.now()

def chrono():
    return datetime.datetime.now() - G_time_0

# Hacking to allow spaces in string fits data...
#    with patch(fits.fitsrec, _rstrip_inplace):
#        hdulist.writeto(...)
#    ...
@contextmanager
def patch(module, funcname, newfunc=lambda x: x):
    oldfunc = getattr(module, funcname)
    setattr(module, funcname, newfunc)
    yield
    setattr(module, funcname, oldfunc)


class FrameConverterInterface(object):
    def get_wcs(self):
        pass

    def get_raw_image(self):
        pass

    def get_background(self):
        pass

    def get_weight(self):
        pass

    def get_mask(self):
        pass

    def get_psf(self):
        pass

    def get_mag_zp(self):
        pass

    def get_err_mag_zp(self):
        pass

    def get_detector(self):
        pass

    def get_filter(self):
        pass

    def get_observation_time(self):
        pass


class LSSTFrameConverter(FrameConverterInterface):
    def __init__(self, butler_base_path, input="input", output="output"):
        self.path_in = os.path.join(butler_base_path, input)
        self.path_out = os.path.join(butler_base_path, output)
        self._butler_in = None
        self._butler_out = None
        self._calexp = None
        self._catalog = None
        self._calexpBackground = None
        self._dataid = None
        self._eimage = None

    @property
    def butler_in(self):
        if self._butler_in is None:
            #print('start open butler in: ', chrono())            
            self._butler_in = lsst.daf.persistence.Butler(self.path_in)
            #print('end open butler in : ', chrono())
        return self._butler_in

    @property
    def butler_out(self):
        if self._butler_out is None:
            #print('start open butler out: ', chrono())
            self._butler_out = lsst.daf.persistence.Butler(self.path_out)
            #print('start open butler out: ', chrono())
        return self._butler_out

    @property
    def calexp(self):
        if self._calexp is None:
            self._calexp = self.butler_out.get("calexp", self.dataid)
        return self._calexp

    @property
    def catalog(self):
        if self._catalog is None:
            self._catalog = self.butler_out.get("src", self.dataid)
        return self._catalog

    @property
    def dataid(self):
        if self._dataid is None:
            raise ValueError("dataid is undefined")
        return self._dataid

    @dataid.setter
    def dataid(self, value):
        dd = dict(value)
        if self._dataid != dd:
            self._dataid = dd
            self._calexp = None
            self._eimage = None
            self._calexpBackground = None

    @property
    def calexpBackground(self):
        if self._calexpBackground is None:
            self._calexpBackground = self.butler_out.get(
                "calexpBackground", self.dataid
            )
        return self._calexpBackground

    @property
    def metadata(self):
        return self.eimage.getMetadata()

    @property
    def eimage(self):
        if self._eimage is None:
            self._eimage = self.butler_in.get("eimage", self.dataid)
        return self._eimage

    def get_wcs(self):
        wcs = self.calexp.getWcs()
        wcs_dict = wcs.getFitsMetadata().toDict()
        wcs_tpv = sip2tpv(wcs_dict)
        return wcs_tpv

    def get_raw_image(self):
        return self.eimage.image.convertF().array

    def get_background(self):
        background = self.calexpBackground.getImage().array
        calexp = self.butler_out.get("calexp", self.dataid)
        im_calexp = calexp.getImage().array
        bkg_raw_unit = background - background.mean() 
        bkg_raw_unit += (self.get_raw_image() - im_calexp).mean()        
        return bkg_raw_unit

    def get_weight(self):
        variance = self.calexp.getVariance()
        weight = 1 / variance.array
        return weight

    def get_mask(self):
        mask = self.calexp.getMask()
        return mask.array

    def get_detector(self):
        return self.metadata.get("CHIPID").replace("_", ":")

    def get_filter(self):
        return self.metadata.get("FILTER")

    def get_observation_time(self):
        return self.metadata.get("DATE")

    def get_mag_zp(self):
        flux0, errflux0 = self.calexp.getCalib().getFluxMag0()
        mag_zp = 2.5 * np.log10(flux0)
        return mag_zp

    def get_err_mag_zp(self):
        flux0, errflux0 = self.calexp.getCalib().getFluxMag0()
        err_mag_zp = 2.5 * errflux0 / flux0 / np.log(10)
        return err_mag_zp

    def get_psf(self):
        psf = self.calexp.getPsf()
        psf_info = extract_psf_info(psf)
        return psf_info


class FITSCreator(object):
    def __init__(self):
        self._hdulist = None

    @property
    def hdulist(self):
        if self._hdulist is None:
            self.prepare_hdulist()
        return self._hdulist

    def prepare_hdulist(self):
        raise

    def write(self, filename, overwrite=False):
        self.hdulist.writeto(filename, overwrite=overwrite)


class DummyFITSCreator(FITSCreator):
    def __init__(self):
        super().__init__()

    def prepare_hdulist(self):
        hdulist = fits.HDUList()
        # Primary HDU is empty
        prim = fits.PrimaryHDU()
        hdulist.append(prim)
        # Finalizei
        self._hdulist = hdulist


class CatalogFITSCreator(FITSCreator):
    def __init__(self, butler_proxy=None, header_image=None):
        if butler_proxy:
            assert isinstance(butler_proxy, LSSTFrameConverter)
            self.butler_proxy = butler_proxy
        else:
            self.butler_proxy = None
        self._hdulist = None
        self.header_image = header_image

    def create_empty_table_cat_euclid(self, nb_row):
        """
        <Column name="NUMBER" unit="" format="J"/>
        <Column name="XWIN_IMAGE" unit="pixel" format="D"/>
        <Column name="YWIN_IMAGE" unit="pixel" format="D"/>
        <Column name="ERRAWIN_IMAGE" unit="pixel" format="E"/>
        <Column name="ERRBWIN_IMAGE" unit="pixel" format="E"/>
        <Column name="ERRTHETAWIN_IMAGE" unit="deg" format="E"/>
        <Column name="FLAGS" unit="" format="I"/>
        <Column name="FLAGS_WEIGHT" unit="" format="I"/>
        <Column name="FLAGS_WIN" unit="" format="I"/>
        <Column name="IMAFLAGS_ISO" unit="" format="J"/>
        <Column name="FLUX_RADIUS" unit="pixel" format="E"/>
        <Column name="FLUX_PSF" unit="count" format="E"/>
        <Column name="FLUXERR_PSF" unit="count" format="E"/>
        <Column name="FLUX_AUTO" unit="count" format="E"/>
        <Column name="FLUXERR_AUTO" unit="count" format="E"/>
        <Column name="SPREAD_MODEL" unit="" format="E"/>
        <Column name="SPREADERR_MODEL" unit="" format="E"/>
        <Column name="CLASS_STAR" unit="" format="E"/>
        
        
        # mapping from TFORM data type to numpy data type (code)
        # L: Logical (Boolean)
        # B: Unsigned Byte
        # I: 16-bit Integer
        # J: 32-bit Integer
        # K: 64-bit Integer
        # E: Single-precision Floating Point
        # D: Double-precision Floating Point
        # C: Single-precision Complex
        # M: Double-precision Complex
        # A: Character
        FITS2NUMPY = {'L': 'i1', 'B': 'u1', 'I': 'i2', 'J': 'i4', 'K': 'i8', 
                      'E': 'f4', 'D': 'f8', 'C': 'c8', 'M': 'c16', 'A': 'a'}
        """
        cat_euclid = at.Table()
        # replace type J by K, LSST stack used int 64 bits for id
        m_column = at.Column(name="NUMBER", dtype=FITS2NUMPY["K"])
        cat_euclid.add_column(m_column)
        m_column = at.Column(
            name="XWIN_IMAGE",
            unit="pixel",
            dtype=FITS2NUMPY["D"],
            format="%2.10",
        )
        cat_euclid.add_column(m_column)
        m_column = at.Column(
            name="YWIN_IMAGE", unit="pixel", dtype=FITS2NUMPY["D"]
        )
        cat_euclid.add_column(m_column)
        m_column = at.Column(
            name="ERRAWIN_IMAGE", unit="pixel", dtype=FITS2NUMPY["E"]
        )
        cat_euclid.add_column(m_column)
        m_column = at.Column(
            name="ERRBWIN_IMAGE", unit="pixel", dtype=FITS2NUMPY["E"]
        )
        cat_euclid.add_column(m_column)
        m_column = at.Column(
            name="ERRTHETAWIN_IMAGE", unit="deg", dtype=FITS2NUMPY["E"]
        )
        cat_euclid.add_column(m_column)
        m_column = at.Column(name="FLAGS", dtype=FITS2NUMPY["I"])
        cat_euclid.add_column(m_column)
        m_column = at.Column(name="FLAGS_WEIGHT", dtype=FITS2NUMPY["I"])
        cat_euclid.add_column(m_column)
        m_column = at.Column(name="FLAGS_WIN", dtype=FITS2NUMPY["I"])
        cat_euclid.add_column(m_column)
        m_column = at.Column(name="IMAFLAGS_ISO", dtype=FITS2NUMPY["J"])
        cat_euclid.add_column(m_column)
        m_column = at.Column(
            name="FLUX_RADIUS", unit="pixel", dtype=FITS2NUMPY["E"]
        )
        cat_euclid.add_column(m_column)
        m_column = at.Column(
            name="FLUX_PSF", unit="count", dtype=FITS2NUMPY["E"]
        )
        cat_euclid.add_column(m_column)
        m_column = at.Column(
            name="FLUXERR_PSF", unit="count", dtype=FITS2NUMPY["E"]
        )
        cat_euclid.add_column(m_column)
        m_column = at.Column(
            name="FLUX_AUTO", unit="count", dtype=FITS2NUMPY["E"]
        )
        cat_euclid.add_column(m_column)
        m_column = at.Column(
            name="FLUXERR_AUTO", unit="count", dtype=FITS2NUMPY["E"]
        )        
        cat_euclid.add_column(m_column)
        m_column = at.Column(
            name="MAG_PSF", unit="mag", dtype=FITS2NUMPY["E"]
        )        
        cat_euclid.add_column(m_column)
        m_column = at.Column(
            name="MAGERR_PSF", unit="mag", dtype=FITS2NUMPY["E"]
        )        
        cat_euclid.add_column(m_column)
        m_column = at.Column(name="SPREAD_MODEL", dtype=FITS2NUMPY["E"])
        cat_euclid.add_column(m_column)
        m_column = at.Column(name="SPREADERR_MODEL", dtype=FITS2NUMPY["E"])
        cat_euclid.add_column(m_column)
        m_column = at.Column(name="CLASS_STAR", dtype=FITS2NUMPY["E"])
        cat_euclid.add_column(m_column)
#        self._cat_euclid = cat_euclid
#         data = np.zeros(1, dtype=cat_euclid.dtype)[0]
#         for idx in range(nb_row):
#             self._cat_euclid.add_row(data)
        self._cat_euclid = at.Table(data = np.zeros(nb_row,dtype=cat_euclid.dtype)) 

    def set_stack_cat_from_fits(self, cat_lsst_fits):
        # open LSST_stack catalog fits
        hdul = fits.open(cat_lsst_fits)
        self.stack_cat = hdul[1].data
        hdul.close()

    def set_stack_cat(self, stack_cat):
        self.stack_cat = stack_cat

    def set_stack_cat_from_butler(self):
        self.stack_cat = self.butler_proxy.catalog
                

    def convert_to_sextractor(self):
        # create table cat_euclid
        nb_elet = len(self.stack_cat)
        self.create_empty_table_cat_euclid(nb_elet)
        # copy column
        self._cat_euclid["NUMBER"] = self.stack_cat["id"]
        # name="base_SdssShape_x",
        # doc="elliptical Gaussian adaptive moments", units="pixel")
        # name="base_SdssShape_y",
        # doc="elliptical Gaussian adaptive moments", units="pixel")
        # Add 1 to compliant with FITS convention
        self._cat_euclid["XWIN_IMAGE"] = self.stack_cat["base_SdssShape_x"] + 1.0
        self._cat_euclid["YWIN_IMAGE"] = self.stack_cat["base_SdssShape_y"] + 1.0
        # compute ERRAWIN_IMAGE, ERRBWIN_IMAGE, ERRTHETAWIN_IMAGE
        # https://sextractor.readthedocs.io/en/latest/PositionWin.html#poserr-win-def
        #
        # name="base_SdssShape_xx",
        # doc="elliptical Gaussian adaptive moments", units="pixel^2")
        # name="base_SdssShape_yy",
        # doc="elliptical Gaussian adaptive moments", units="pixel^2")
        # name="base_SdssShape_xy",
        # doc="elliptical Gaussian adaptive moments", units="pixel^2")
        var_x = self.stack_cat["base_SdssShape_xx"]
        var_y = self.stack_cat["base_SdssShape_yy"]
        cov_xy = self.stack_cat["base_SdssShape_xy"]
        s_2 = (var_x + var_y) / 2
        dvar = var_x - var_y
        sqrt_vc = np.sqrt((dvar / 2) ** 2 + cov_xy ** 2)
        self._cat_euclid["ERRAWIN_IMAGE"] = np.sqrt(s_2 + sqrt_vc)
        self._cat_euclid["ERRBWIN_IMAGE"] = np.sqrt(s_2 - sqrt_vc)
        self._cat_euclid["ERRTHETAWIN_IMAGE"] = np.rad2deg(
            np.arctan2(2 * cov_xy, dvar) / 2
        )
        # flags
        # name="base_SdssShape_flag", doc="General Failure Flag")
        try:
            self._cat_euclid["FLAGS"] = self.stack_cat["base_SdssShape_flag"]
        except:
            print("Flags by default")
            self._cat_euclid["FLAGS"] = 0
        # (name="base_SdssShape_flag_unweighted",
        #  doc="Weighted moments converged to an invalid value; "
        #      "using unweighted moments")
        try:
            self._cat_euclid["FLAGS_WEIGHT"] = self.stack_cat[
                "base_SdssShape_flag_unweighted"
            ]
        except:
            # fits catalog from stack doesn't have flag !?
            print("Flags by default")
            self._cat_euclid["FLAGS_WEIGHT"] = 0
        self._cat_euclid["FLAGS_WIN"] = 0
        self._cat_euclid["IMAFLAGS_ISO"] = 0
        # Flux
        #
        # SEXTRACTOR:
        #         "FLUX_AUTO", unit="count"
        #         "FLUXERR_AUTO", unit="count"
        # FLUX_AUTO provides an estimate of the total flux by integrating
        # pixel values within
        # an adaptively scaled aperture. SExtractor's automatic aperture
        # photometry routine
        # derives from Kron's first moment algorithm
        #
        # (name="ext_photometryKron_KronFlux_flux",
        #  doc="flux from Kron Flux algorithm", units="count")
        self._cat_euclid["FLUX_AUTO"] = self.stack_cat[
            "ext_photometryKron_KronFlux_flux"
        ]
        # name="ext_photometryKron_KronFlux_fluxSigma",
        # doc="1-sigma flux uncertainty", units="count")
        self._cat_euclid["FLUXERR_AUTO"] = self.stack_cat[
            "ext_photometryKron_KronFlux_fluxSigma"
        ]
        #
        #         "FLUX_RADIUS", unit="pixel"
        #  no reference
        #     http://astroa.physics.metu.edu.tr/MANUALS/sextractor/Guide2source_extractor.pdf
        #     https://sextractor.readthedocs.io/en/latest/Param.html
        # but definition in astromatic forum
        # http://www.astromatic.net/forum/showthread.php?tid=318
        # FLUX_RADIUS estimates the radius of the circle centered
        # on the barycenter that encloses
        # about half of the total flux. For a Gaussian profile, this is
        # equal to 1/2 FWHM. But with most images on astronomical images
        # it will be slightly higher.
        #
        #  to my knowledge no equivalent in catalog LSST
        self._cat_euclid["FLUX_RADIUS"] = 0
        # flux psf
        #         "FLUX_PSF", unit="count"
        #         "FLUXERR_PSF", unit="count"
        #
        #  (Field['D'](name="base_PsfFlux_flux",
        #              doc="flux derived from linear least-squares "
        #                  "fit of PSF model", units="count")
        self._cat_euclid["FLUX_PSF"] = self.stack_cat["base_PsfFlux_flux"]
        #  (Field['D'](name="base_PsfFlux_fluxSigma",
        #              doc="1-sigma flux uncertainty", units="count")
        self._cat_euclid["FLUXERR_PSF"] = self.stack_cat[
            "base_PsfFlux_fluxSigma"
        ]
        
        # MAG_PSF
        self._cat_euclid["MAG_PSF"] = self.butler_proxy.get_mag_zp() - 2.5*np.log10(self._cat_euclid["FLUX_PSF"])
        self._cat_euclid["MAGERR_PSF"] = (2.5/ np.log(10))*self._cat_euclid["FLUXERR_PSF"]/self._cat_euclid["FLUX_PSF"]
        #
        # classification object
        #
        # name="base_ClassificationExtendedness_value",
        # doc="Set to 1 for extended sources, 0 for point sources.")
        self._cat_euclid["SPREAD_MODEL"] = self.stack_cat[
            "base_ClassificationExtendedness_value"
        ]
        self._cat_euclid["CLASS_STAR"] = self.stack_cat[
            "base_ClassificationExtendedness_value"
        ]
        self._cat_euclid["SPREADERR_MODEL"] = 1
        
        
    def check_cat_convert(self):
        """ remove NAN in euclid catalog 
        
        remove object with NAN position or error position
        remove object with NAN flux for psf and auto estimator
        replace NAN flux *psf* or auto by auto or *psf*
        replace NAN classicator by 0.5
        """
        # check position
        mask_x = np.isnan(self._cat_euclid["XWIN_IMAGE"])
        mask_y = np.isnan(self._cat_euclid["YWIN_IMAGE"])
        n_total = mask_x.size
        mask_xy = np.logical_or(mask_x, mask_y)
        n_nan_pos = np.where(mask_xy == True)[0].size
        print(f"Number of sources  : {n_total}")
        print(f"   NAN from position : {n_nan_pos}")        
        self._cat_euclid = self._cat_euclid[np.logical_not(mask_xy)]
        # check error position
        mask_erra = np.isnan(self._cat_euclid["ERRAWIN_IMAGE"])
        mask_errb = np.isnan(self._cat_euclid["ERRBWIN_IMAGE"])
        mask_errt = np.isnan(self._cat_euclid["ERRTHETAWIN_IMAGE"])
        mask_err = np.logical_or(np.logical_or(mask_erra, mask_errb), mask_errt)
        n_nan_err = np.where(mask_err == True)[0].size        
        print(f"   NAN from error position : {n_nan_err}")        
        self._cat_euclid = self._cat_euclid[np.logical_not(mask_err)]
        # check flux PSF
        n_nan_flux = 0
        l_row_nok = []
        for idx , row_src in enumerate(self._cat_euclid):
            nan_auto = np.isnan(row_src["FLUX_AUTO"]) or np.isnan(row_src["FLUXERR_AUTO"])
            nan_psf = np.isnan(row_src["FLUX_PSF"]) or np.isnan(row_src["FLUXERR_PSF"])     
            if nan_auto and nan_psf:
                n_nan_flux += 1
                l_row_nok.append(idx)
            elif nan_auto:
                # replace NAN auto by psf
                row_src["FLUX_AUTO"] = row_src["FLUX_PSF"]
                row_src["FLUXERR_AUTO"] = row_src["FLUXERR_PSF"]
            elif nan_psf:
                # replace NAN psf by auto
                row_src["FLUX_PSF"] = row_src["FLUX_AUTO"]
                row_src["FLUXERR_PSF"] = row_src["FLUXERR_AUTO"]
            # replace NAN classificator by 0.5
            if np.isnan(row_src["CLASS_STAR"]):
                row_src["CLASS_STAR"] = 0.5
                row_src["SPREAD_MODEL"] = 0.5        
        # remove object with NAN flux
        self._cat_euclid.remove_rows(l_row_nok)
        print(f"   NAN from flux : {n_nan_flux}")
        
        
        
    def create_LDAC_IMHEAD(self):
        self.ldac_imhead = at.Table(
            dtype=[("Field Header Card", "S80", (len(self.header_image),))]
        )
        data = np.zeros(1, dtype=self.ldac_imhead.dtype)[0]
        self.ldac_imhead.add_row(data)
        l_kw_str = repr(self.header_image).split("\n")
        for idx, kw_str in enumerate(l_kw_str):
            self.ldac_imhead[0][0][idx] = kw_str
            if G_verbose >= 2: 
                print(kw_str)
        if G_verbose >= 2:
            print(self.ldac_imhead)

    def prepare_hdulist(self):
        # empty primary
        prihdr = fits.Header()
        prihdu = fits.PrimaryHDU(header=prihdr)
        # Extension 1: empty extension (?)
        self.create_LDAC_IMHEAD()
        ldac_imhead = fits.BinTableHDU(
            data=self.ldac_imhead, name="LDAC_IMHEAD"
        )
        # Extension 2: catalog
        self.convert_to_sextractor()
        self.check_cat_convert()
        ldac_object = fits.BinTableHDU(
            data=self._cat_euclid, name="LDAC_OBJECTS"
        )
        self._hdulist = fits.HDUList([prihdu, ldac_imhead, ldac_object])


class DetrendedFrameFITSCreator(FITSCreator):
    def __init__(self, frame_converter=None):
        self._hdulist = None
        self.butler_proxy = frame_converter

    @property
    def butler_proxy(self):
        return self._butler_proxy

    @butler_proxy.setter
    def butler_proxy(self, value):
        if not isinstance(value, FrameConverterInterface):
            raise TypeError(
                "expected type or subtype of FrameConverterInterface, "
                f"got {type(value)}"
            )
        self._butler_proxy = value

    def prepare_hdulist(self):
        hdulist = fits.HDUList()
        bg = self._butler_proxy.get_background()
        
        # Primary HDU is empty
        prim = fits.PrimaryHDU()
        hdulist.append(prim)

        # Extension 1: raw image
        raw_image = self._butler_proxy.get_raw_image()
        image_hdr = fits.Header(
            cards={
                "EXTNAME": "IMAGE",
                "PROCTYPE": "RAW",
                "PRODTYPE": "image",
                "OBS-TIME": self._butler_proxy.get_observation_time(),
                "FILTER": f'LSST_{self._butler_proxy.get_filter()}',
                "DETECTOR": self._butler_proxy.get_detector(),
                "ZEROPNT": self._butler_proxy.get_mag_zp(),
                "MAGZP": self._butler_proxy.get_mag_zp(),
                "ERRMAGZP": self._butler_proxy.get_err_mag_zp(),
                #TODO: replace hardcoding value by the right butler query
                "EXPTIME": 30,
                "INSTRUME": "LSSTCAM",
                "SKYBRITE": bg.mean(),
                "SKYSIGMA": bg.std(),
                "CCDSECA" : '[1:4004,1:4096]',
                "GAIN_A" : 1.0
            }
        )
        image_hdr.update(self._butler_proxy.get_wcs())
        imghdu = fits.ImageHDU(data=raw_image, header=image_hdr, uint=False)
        self.raw_image_hdr = imghdu.header
        hdulist.append(imghdu)

        # Extension 2: Mask
        mask = self._butler_proxy.get_mask()
        mask_hdr = fits.Header(cards=dict(EXTNAME="MASK"))
        mask_hdr.update(self._butler_proxy.get_wcs())
        hdulist.append(fits.ImageHDU(data=mask, header=mask_hdr))

        # Extension 3: RMS
        rms = np.sqrt(raw_image)
        rms_hdr = fits.Header(cards=dict(EXTNAME="RMS"))
        hdulist.append(fits.ImageHDU(data=rms, header=rms_hdr))

        # Extension 4: Background
        bg = self._butler_proxy.get_background()
        bg_hdr = fits.Header(cards=dict(EXTNAME="BACKGROUND"))
        hdulist.append(fits.ImageHDU(data=bg, header=bg_hdr))

        # Finalize
        self._hdulist = hdulist


class PsfFITSCreator(FITSCreator):
    def __init__(self, frame_converter=None):
        self._hdulist = None
        if frame_converter is None:
            self._frame_converter = None
        else:
            self.frame_converter = frame_converter

    @property
    def frame_converter(self):
        return self._frame_converter

    @frame_converter.setter
    def frame_converter(self, value):
        if not isinstance(value, FrameConverterInterface):
            raise ValueError(
                f"expected type {type(FrameConverterInterface)}, "
                f"got {type(value)}"
            )
        self._frame_converter = value

    def prepare_hdulist(self):
        hdulist = fits.HDUList()

        # Primary HDU is empty
        prim = fits.PrimaryHDU()
        hdulist.append(prim)

        # Extension
        psf = self.frame_converter.get_psf()
        column = fits.Column(
            name="PSF_MASK",
            format=f"{psf.data.size}E",
            dim=f"{psf.shape}",
            array=psf.data,
            #             comment="Tabulated PSF data",
        )
        hdu = fits.BinTableHDU.from_columns(
            [column], name="PSF_DATA", header=psf.header
        )

        hdulist.append(hdu)

        # Finalize
        self._hdulist = hdulist


class CreateProductsLSSTwithDummy(object):
    def __init__(
        self,
        path,
        visit,
        raft,
        sensor,
        output_dir=None,
        overwrite=False,
        tag="",
    ):
        self.path = path
        self.visit = visit
        self.sensor = sensor
        self.raft = raft
        self.output_dir = output_dir
        self.overwrite = overwrite
        self.tag = tag

    def create_all_fits(self):
        #print('start create_all_fits: ', chrono())
        self.make_detrended_frame_fits()
        #print('make_detrended_frame_fits done: ', chrono())
        self.make_psf_fits()
        #print('make_psf_fits done: ', chrono())        
        self.make_cat_fits()
        #print('make_cat_fits done: ', chrono())

    def _add_output_dir(self, filename):
        if self.output_dir is not None:
            filename = os.path.join(self.output_dir, filename)
        return filename

    def get_date(self):
        now = datetime.datetime.now(pytz.UTC)
        return now.strftime("%Y%m%dT%H%M%S.%f")[:-3] + "Z"

    def get_name_product(self, pref_prod):
        date = self.get_date()
        filename = (
            f"EUC_EXT_{pref_prod}_LSST-{self.visit}-"
            f"{self.detector}_{date}{self.tag}"
        )
        return self._add_output_dir(filename)

    def make_psf_fits(self):
        ofits = DummyFITSCreator()
        ofits.write(
            self.get_name_product("DPDEXTPSFMODEL") + "_EMPTY.fits",
            self.overwrite,
        )

    def make_cat_fits(self):
        ofits = DummyFITSCreator()
        ofits.write(
            self.get_name_product("DPDEXTSOURCECATALOG") + "_EMPTY.fits",
            self.overwrite,
        )

    def make_detrended_frame_fits(self):
        self.butler_proxy = LSSTFrameConverter(self.path)
        detrendedfitscreator = DetrendedFrameFITSCreator(self.butler_proxy)
        detrendedfitscreator.butler_proxy.dataid = dict(
            visit=self.visit, raft=self.raft, sensor=self.sensor
        )
        detector = (
            "R"
            + self.raft.replace(",", "")
            + "-S"
            + self.sensor.replace(",", "")
        )
        self.detector = detector
        # EUC_EXT_DPDEXTDETRENDEDFRAME_LSST-100482-R11-S11_20290906T032512.fits
        # date = butler_proxy.get_observation_time().replace(':','')\
        #                    .replace('-','').split('.')[0]
        filename = self.get_name_product("DPDEXTDETRENDEDFRAME") + ".fits"
        detrendedfitscreator.write(filename, self.overwrite)
        self.raw_image_hdr = detrendedfitscreator.raw_image_hdr
        # discussion in course for compression, not add
        os.system(f"fpack -r -D -Y {filename}")


def extract_psf_info(psf):
    """From a LSST psf object, get information to fill a PSFEx file.
    """
    # Write to a temporary fits file
    tmpfile = tempfile.NamedTemporaryFile(
        mode="wb", suffix=".fits", delete=False
    )
    tmpfile.close()
    psf.writeFits(tmpfile.name)
    # Read back the fits file
    hdl = fits.open(tmpfile.name)
    # Create data, shape and header
    data2 = hdl[2].data
    data3 = hdl[3].data
    data_array = data3._comp
    shape = tuple(*hdl[3].data._size)
    cards = [
        fits.Card(keyword=kw, value=val, comment=com)
        for kw, val, com in (
            [
                ("LOADED", 0, "Number of loaded sources"),
                ("ACCEPTED", 0, "Number of accepted sources"),
                ("CHI2", 1.0, "Final Chi2"),
                (
                    "POLNAXIS",
                    int(data2._context_size),
                    "Number of context parameters",
                ),
                (
                    "POLGRP1",
                    data3.group[0, 0],
                    "Polynom group for this context parameter",
                ),
                ("POLNAME1", "X_IMAGE", "Name of this context parameter"),
                (
                    "POLZERO1",
                    data3._context_first[0, 0],
                    "Offset value for this context parameter",
                ),
                (
                    "POLSCAL1",
                    data3._context_second[0, 0],
                    "Scale value for this context parameter",
                ),
                (
                    "POLGRP2",
                    data3.group[0, 1],
                    "Polynom group for this context parameter",
                ),
                ("POLNAME2", "Y_IMAGE", "Name of this context parameter"),
                (
                    "POLZERO2",
                    data3._context_first[0, 1],
                    "Offset value for this context parameter",
                ),
                (
                    "POLSCAL2",
                    data3._context_second[0, 1],
                    "Scale value for this context parameter",
                ),
                ("POLNGRP", len(data3.degree), ""),
            ]
            + [(f"POLDEG{i+1}", d, "") for i, d in enumerate(data3.degree)]
            + [
                ("PSF_FWHM", 3.5, "PSF FWHM"),  # value in pixels
                (
                    "PSF_SAMP",
                    float(data2._pixstep),
                    "Sampling step of the PSF data",
                ),
                ("PSFNAXIS", len(shape), "Dimensionnality of the PSF data"),
                ("PSFAXIS1", shape[0], "Number of element along this axis"),
                ("PSFAXIS2", shape[1], "Number of element along this axis"),
                ("PSFAXIS3", shape[2], "Number of elements along this axis"),
            ]
        )
    ]

    # delete the temporary file
    hdl.close()
    os.remove(tmpfile.name)

    return PsfInfo(data=data_array, header=fits.Header(cards), shape=shape)


class CreateProductsLSST(CreateProductsLSSTwithDummy):
    def __init__(
        self,
        path,
        visit,
        raft,
        sensor,
        output_dir=None,
        overwrite=False,
        tag="",
    ):
        super().__init__(path, visit, raft, sensor, output_dir, overwrite, tag)

    def make_psf_fits(self):
        ofits = PsfFITSCreator(self.butler_proxy)
        ofits.write(
            self.get_name_product("DPDEXTPSFMODEL") + ".fits", self.overwrite
        )

    def make_cat_fits(self):
        ofits = CatalogFITSCreator(self.butler_proxy, self.raw_image_hdr)
        ofits.set_stack_cat_from_butler()
        ofits.convert_to_sextractor()
        # TODO: find a way to remove this ugly hack (to write spaces and not
        #       zeros at the end of string)
        with patch(fits.fitsrec, "_rstrip_inplace"):
            ofits.write(
                self.get_name_product("DPDEXTSOURCECATALOG") + ".fits",
                self.overwrite,
            )


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
        description="Create Detrended Frame for a visit/raft/sensor image."
    )
    parser.add_argument("path", help="path to butler (with input and output).")
    parser.add_argument("--visit", type=int, help="visit id")
    parser.add_argument("--raft", help="raft i,j")
    parser.add_argument("--sensor", help="sensor i,j")
    parser.add_argument("-o", "--output-dir", help="output directory")
    parser.add_argument("-t", "--tag", help="tag name product", default="")
    parser.add_argument(
        "--force", action="store_true", help="proceed without asking"
    )

    args = parser.parse_args()
    print("Extact FITS from butler for Euclid")
    if False:
        print(f"path: {args.path}")
        print(f"visit: {args.visit}")
        print(f"raft: {args.raft}")
        print(f"sensor: {args.sensor}")
        print(f"output-dir: {args.output_dir}")
    o_prod = CreateProductsLSST(
        args.path,
        args.visit,
        args.raft,
        args.sensor,
        output_dir=args.output_dir,
        overwrite=args.force,
        tag=args.tag,
    )
    #print('CreateProductsLSST done: ', chrono())
    o_prod.create_all_fits()
    print('create_all_fits done: ', chrono())


# example:
# python lsst2euc.py /sps/euclid/Users/colley/lsst/butler/stage1_4T6_100482
#         --visit 100482 --raft 1,2 --sensor 1,1 --force
