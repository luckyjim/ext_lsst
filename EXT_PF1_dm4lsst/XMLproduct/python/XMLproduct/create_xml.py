#
# Copyright (C) 2012-2020 Euclid Science Ground Segment
#
# This library is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 3.0 of the License, or (at your option)
# any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
#

"""
File: python/XMLproduct/create_xml.py

Created on: 06/29/20
Author: Colley
"""
import pyxb
import ST_DataModelBindings.pro.ext_stub as ext_dict
import ST_DataModelBindings.bas.img_stub as img_dict
import ST_DataModelBindings.bas.dqc.ext_stub as extdqc_dict
import ST_DataModelBindings.bas.raw.dtd_stub as dtdDict
import ST_DataModelBindings.bas.ppr.par_stub as par_dict
import ST_DataModelBindings.dpd.ext.singleepochframe_stub as ext_singleepochframe
import ST_DataModelBindings.bas.cot_stub as cotDict
import ST_DataModelBindings.bas.imp.raw.stc_stub as stcDict
try:
    import HeaderProvider.GenericHeaderProvider as HeaderProvider
except:
    import ST_DM_HeaderProvider.GenericHeaderProvider as HeaderProvider

import XMLproduct.dmUtils as dmUtils
import XMLproduct.ExtDmUtils as extUtils

from ST_DM_DmUtils import DmUtils as dm_utils

from astropy.io import fits
import os.path

import pytz
import datetime
import copy
import numpy as np
from astropy.coordinates import Angle 
from astropy import units as u

from astropy import wcs
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import pixel_to_skycoord


class BBox(object):
    def __init__(self):
        self.ra_center = 0
        self.dec_center = 0
        self.ra_half = 0
        self.dec_half = 0
    
    def from_footprint(self, polygon, center):
        #print(polygon)
        a_poly = Angle(polygon, unit='deg')
        a_center = Angle(center, unit='deg')
        poly_center = a_poly - a_center        
        self.ra_center = a_center[0]
        self.dec_center = a_center[1]
        self.ra_half = np.abs(poly_center[:,0]).max()
        self.dec_half = np.abs(poly_center[:,1]).max()
    
    def __str__(self):
        ret = f'center ra, dec: {self.ra_center}, {self.dec_center}'
        ret += f'\nhalf   ra, dec: {self.ra_half}, {self.dec_half}'
        return ret                                    

    

class DataModelBase(object):
    
    def save_xml_product(self, xml_file_name):
        with open(xml_file_name, "w") as f:
            f.write(self.product.toDOM().toprettyxml(encoding="utf-8").decode("utf-8"))


    def read_xml_product(self, xml_file_name):
        # Read the xml file as a string
        with open(xml_file_name, "r") as f:
            xml_string = f.read()
    
        # Create a new NIR product instance using the NIR data product dictionary
        product = self.root_prod.CreateFromDocument(xml_string)    
        return product


    def create_fill_header(self, sci_grp_name):
        """
        Create header structure and fill it with default values 
        
        <xs:simpleType name="scientificGroupName">
        
        """
        hdr = HeaderProvider.create_generic_header(sci_grp_name)        
        hdr.ProductId = "EUC_EXT_DPDEXTSINGLEEPOCHFRAME_LSST"        
        hdr.SoftwareName = "EXT_LSST_PF1"
        hdr.SoftwareRelease = "new_dm"
        hdr.EuclidPipelineSoftwareRelease = "new_dm"
        hdr.ProdSDC = "SDC-FR"
        hdr.DataSetRelease = "TEST"
        hdr.Purpose = "REPROCESSING"
        hdr.PlanId = "TBD"
        hdr.PPOId = "TBD"
        hdr.PipelineDefinitionId = "1.0"
        hdr.PipelineRun = "defined_later"
        hdr.ExitStatusCode = "0"
        hdr.ManualValidationStatus = "VALID"
        hdr.ExpirationDate = "2020-06-04T13:48:28.165Z"
        hdr.ToBePublished = "false"
        hdr.Published = "false"
        hdr.Curator = "CNRS/IN2P3/APC"        
        return hdr
        
        

class SingleEpochFrameLSSTDefault(DataModelBase):
    """
    create a single Epoch frame XML for LSST with default values
    """
    
    def __init__(self):
        # Create the appropriate data product binding
        self.product = ext_singleepochframe.DpdExtSingleEpochFrame()
        dp = self.product
    
        # Add the generic header to the data product
        dp.Header = self.create_fill_header("DpdExtSingleEpochFrame")

        # Add the data element to the data product
        dp.Data = dmUtils.createGroundInstrumentImage(
            ext_dict.extSingleEpochFrame, "LSST", "LSSTCAM", "LSST_u")
                
        # Add the optional ImageArea
        dp.Data.ImgArea = pyxb.BIND()
        dp.Data.ImgArea.Name = "tbd"
        
        # Add the ProcessParams
        dp.Data.ProcessParams = par_dict.extScienceFrameParameters()
    
        # Add the QualityParmas
        dp.Data.QualityParams = ext_dict.extCalibratedFrameDqp()
    
        # Add the ObservationDateTime
        dp.Data.ObservationDateTime = img_dict.groundObservationDateTime()
    
        # Add the DataStorage
        dp.Data.DataStorage = dm_utils.create_fits_storage(
            ext_dict.extDetrendedFrameFitsFile, "data.fits", "ext.detrendedFrame",
            "0.1")
    
        # Add the PsdModelStorage
        dp.Data.PsfModelStorage = dm_utils.create_fits_storage(
            ext_dict.extPsfModelFitsFile, "psfModel.fits", "ext.psfModel",
            "0.1")
            
        # Add the CatalogStorage
        dp.Data.CatalogDataStorage = dm_utils.create_fits_storage(
            ext_dict.extSourceCatalogFitsFile, "sourceCatalog.fits", "ext.sourceCatalog",
            "0.1")    
    
        # Add the ObsCollectionDataStorage
        # - =  not in the ext_dict
        #dp.Data.ObsCollectionDataStorage = dm_utils.create_fits_storage(
        #    ext_dict.extObsCollectionDataStorageFitsFile, "obsCollectionDataStorage.fits", "ext.obsCollection",
        #    "0.1")    
    
        dp.Data.FWHM = pyxb.BIND()
    
        # Add the QualityFlags element to the data productSingleEpochFrameStruct
        #dp.QualityFlags = __create_calibrated_frame_quality_flags()
        # Add the QualityFlags element to the data product
        dp.QualityFlags = extUtils.create_quality_flags(
            extdqc_dict.sqfDpdExtCalibratedFrame)
    
        self.fill_lsst_constant()
    
    
    def fill_lsst_constant(self):
        data = self.product.Data        
        data.AxisLengths = "2048 4096"
        data.DataLength = 8388608
        data.Instrument.Longitude = -70.749417
        data.Instrument.Latitude = -30.244639
        data.Instrument.Elevation = 2663
        data.Instrument.Timezone = '-4'


class SingleEpochFrameLSST(SingleEpochFrameLSSTDefault):
    """
    create a single Epoch frame XML for LSST with product values
    """
    def __init__(self, str_filter):
        super().__init__()
        if str_filter == "":
            self.filter_file = False
            return
        mfilter = str_filter.split(',')
        ra = float(mfilter[0])
        dec = float(mfilter[1])
        self.radius = float(mfilter[2])
        self.center_filter = SkyCoord(ra=ra, dec=dec, frame="icrs", unit="deg")
        self.filter_file = True
        
        
    def fill_xml_with_product(self, f_detrended, f_psf, f_cat,
                              out_path=None, tag="noDef",
                              release="TBD"):
        if out_path:    
            self.out_path = out_path
        else:
            # output in same directory than input
            self.out_path = os.path.dirname(f_detrended)
        # print('Save xml in direction:  ', self.out_path)
        self.f_detrended = os.path.basename(f_detrended)
        self.f_psf = f_psf
        self.f_cat = f_cat
        hdul = fits.open(f_detrended)
        self.hdul = hdul
        bbox_center = [self.hdul[1].header['CRVAL1'], self.hdul[1].header['CRVAL2']]                    
        self.utc_obs = hdul[1].header['OBS-TIME']
        self.filter = hdul[1].header['FILTER']
        self.detector = hdul[1].header['DETECTOR']
        my_wcs = wcs.WCS(hdul[1].header)
        if self.filter_file:            
            center_ccd = SkyCoord(ra=bbox_center[0], dec=bbox_center[1], frame="icrs", unit="deg")            
            #print(center_ccd, self.center_filter)			
            angle = self.center_filter.separation(center_ccd).deg
            #print('dif angle ', angle)
            if angle > self.radius:
                print(f'Filter: Exclude image {angle}>{self.radius}')
                hdul.close()
                return False            
        self.width = hdul[1].header['NAXIS1']
        self.height = hdul[1].header['NAXIS2']
        self.AxisLengths = f"{self.width} {self.height}"
        self.polygon = [ 
            my_wcs.wcs_pix2world(0, 0, 0),
            my_wcs.wcs_pix2world(self.width, 0, 0),
            my_wcs.wcs_pix2world(self.width, self.height, 0),
            my_wcs.wcs_pix2world(0, self.height, 0) ]
  	    #print(self.polygon)
        #BBox
        #bbox_center = my_wcs.wcs_pix2world(self.width/2, self.height/2, 0)
        
        self.bbox = BBox()
        self.bbox.from_footprint(self.polygon, bbox_center)
        #print(self.bbox)
        
        # file_detrended format expected
        # EUC_EXT_DPDEXTASTROMSOLUTION_LSST-100482-R02-S01_20190111T165610.fits
        # visit is 100482
        self.id_ccd = self.f_detrended.split('_')[3]        
        # update default value
        self.product.Header.PipelineRun = tag
        self.product.Header.DataSetRelease = release
        self.detrented_frame_xml()
        self.psf_xml()
        self.photometry_xml()
        self.astrometry_xml()
        self.catalog_xml()
        return True
                
        
    def format_date(self, date):
        # date creation is in header
        # now = datetime.datetime.now(pytz.UTC)
        # return now.strftime('%Y%m%dT%H%M%S.%f')[:-3] + 'Z'
        return date.replace(':', '').replace('-', '').replace(' ', 'T').split('+')[0][:-3] + 'Z'


    def get_name(self, prefixe):
        date = str(self.product.Header.CreationDate)        
        return prefixe + self.id_ccd + '_' + self.format_date(date)
    
    
    def detrented_frame_xml(self):        
        self.product.Header.ProductId = self.get_name("EUC_EXT_DPDEXTSINGLEEPOCHFRAME_")
        data = self.product.Data
        data.ObservationDateTime.UTCObservationDateTime = self.utc_obs
        data.Filter.Name = self.filter
        data.Detector = self.detector
        data.AxisLengths = self.AxisLengths
        data.DataLength = str(self.width * self.height)
        self._create_footprint()
        data.DataStorage.DataContainer.FileName = self.f_detrended


    def psf_xml(self):            
        data = self.product.Data        
        data.PsfModelStorage.DataContainer.FileName = os.path.basename(self.f_psf)
#         qual = self.product.QualityFlags
#         qual.FullWidthHalfMaximum = 10000.0
#         data.FullWidthHalfMaximum.flagged = True
        try: 
            with fits.open(self.f_psf) as f_psf:
                try:
                    fwhm = float(f_psf[1].header['PSF_FWHM'])
                    if fwhm > 0.0:                        
                        data.FWHM = fwhm
                        #TODO: manage quality flag 
                        #data.QualityParams.FullWidthHalfMaximum.flagged = False
                    else:
                        print(f"Can't used PSF_FWHM={fwhm} <= 0 !!!")
                except:
                    print(f"ERROR to set PSF_FWHM value with fits file '{self.f_psf}'")
                    pass
        except:
            print(f"ERROR: file '{self.f_psf}' doesn't exist or isn't a FITS")        
            
    def photometry_xml(self):
        data = self.product.Data
        data.Zeropoint.Value = self.hdul[1].header['MAGZP']
        data.Zeropoint.Error = self.hdul[1].header['ERRMAGZP']
    
    def catalog_xml(self):
        data = self.product.Data
        data.CatalogDataStorage.DataContainer.FileName = os.path.basename(self.f_cat)
    
    def astrometry_xml(self):
        data = self.product.Data.WCS
        data.CRVAL1 = self.hdul[1].header['CRVAL1']
        data.CRVAL2 = self.hdul[1].header['CRVAL2']
        data.CRPIX1 = self.hdul[1].header['CRPIX1']
        data.CRPIX2 = self.hdul[1].header['CRPIX2']
        data.CD1_1 = self.hdul[1].header['CD1_1']
        data.CD1_2 = self.hdul[1].header['CD1_2']
        data.CD2_1 = self.hdul[1].header['CD2_1']
        data.CD2_2 = self.hdul[1].header['CD2_2']
        data.NonLinearCoeffs.NonLinearTPVAstromCoeffs.PV1_0 = self.hdul[1].header['PV1_0']  
        data.NonLinearCoeffs.NonLinearTPVAstromCoeffs.PV1_1 = self.hdul[1].header['PV1_1']
        data.NonLinearCoeffs.NonLinearTPVAstromCoeffs.PV1_2 = self.hdul[1].header['PV1_2']
        data.NonLinearCoeffs.NonLinearTPVAstromCoeffs.PV1_3 = self.hdul[1].header['PV1_3'] 
        data.NonLinearCoeffs.NonLinearTPVAstromCoeffs.PV1_4 = self.hdul[1].header['PV1_4']
        data.NonLinearCoeffs.NonLinearTPVAstromCoeffs.PV1_5 = self.hdul[1].header['PV1_5']
        data.NonLinearCoeffs.NonLinearTPVAstromCoeffs.PV1_6 = self.hdul[1].header['PV1_6']
        data.NonLinearCoeffs.NonLinearTPVAstromCoeffs.PV1_7 = self.hdul[1].header['PV1_7']
        data.NonLinearCoeffs.NonLinearTPVAstromCoeffs.PV1_8 = self.hdul[1].header['PV1_8']
        data.NonLinearCoeffs.NonLinearTPVAstromCoeffs.PV1_9 = self.hdul[1].header['PV1_9']
        data.NonLinearCoeffs.NonLinearTPVAstromCoeffs.PV1_10 = self.hdul[1].header['PV1_10']
        data.NonLinearCoeffs.NonLinearTPVAstromCoeffs.PV2_0 = self.hdul[1].header['PV2_0']
        data.NonLinearCoeffs.NonLinearTPVAstromCoeffs.PV2_1 = self.hdul[1].header['PV2_1']
        data.NonLinearCoeffs.NonLinearTPVAstromCoeffs.PV2_2 = self.hdul[1].header['PV2_2']
        data.NonLinearCoeffs.NonLinearTPVAstromCoeffs.PV2_3 = self.hdul[1].header['PV2_3']
        data.NonLinearCoeffs.NonLinearTPVAstromCoeffs.PV2_4 = self.hdul[1].header['PV2_4']
        data.NonLinearCoeffs.NonLinearTPVAstromCoeffs.PV2_5 = self.hdul[1].header['PV2_5']
        data.NonLinearCoeffs.NonLinearTPVAstromCoeffs.PV2_6 = self.hdul[1].header['PV2_6']
        data.NonLinearCoeffs.NonLinearTPVAstromCoeffs.PV2_7 = self.hdul[1].header['PV2_7']
        data.NonLinearCoeffs.NonLinearTPVAstromCoeffs.PV2_8 = self.hdul[1].header['PV2_8']
        data.NonLinearCoeffs.NonLinearTPVAstromCoeffs.PV2_9 = self.hdul[1].header['PV2_9']
        data.NonLinearCoeffs.NonLinearTPVAstromCoeffs.PV2_10 = self.hdul[1].header['PV2_10']
    
    
    def _create_footprint(self):
        #print(self.polygon)
        Polygon = stcDict.polygonType()        
        # add 4 default values for the corners
        vertex = stcDict.vertexType()
        vertex.Position = dtdDict.double2Type()                     
        vertex.Position.C1 = self.polygon[0][0]
        vertex.Position.C2 = self.polygon[0][1]
        Polygon.append(copy.deepcopy(vertex))
        vertex.Position.C1 = self.polygon[1][0]
        vertex.Position.C2 = self.polygon[1][1]        
        Polygon.append(copy.deepcopy(vertex))        
        vertex.Position.C1 = self.polygon[2][0]
        vertex.Position.C2 = self.polygon[2][1]
        Polygon.append(copy.deepcopy(vertex))        
        vertex.Position.C1 = self.polygon[3][0]
        vertex.Position.C2 = self.polygon[3][1]
        Polygon.append(vertex)       
        self.product.Data.ImgSpatialFootprint.Polygon = Polygon
        # BBOX
        bbox = stcDict.bBox()
        bbox.DecCenter = self.bbox.dec_center.value
        bbox.RaCenter = self.bbox.ra_center.value
        bbox.DecHalfWidth = self.bbox.dec_half.value
        bbox.RaHalfWidth = self.bbox.ra_half.value
        self.product.Data.ImgSpatialFootprint.BBox = bbox

    def save(self):
        # save
        self.hdul.close()
        fname = self.product.Header.ProductId + '.xml'
        fname = os.path.join(self.out_path, fname)
        self.f_xml = fname
        self.save_xml_product(fname)   

        
