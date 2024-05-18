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

Some utility methods to work with the Ext-stage2 Data Model.

Created on: 20/12/16
Author: Javier Gracia Carpio (jgracia@mpe.mpg.de)
and MIchael Wetzstein (drwstein@mpe.mpg.de)
"""

import pyxb

import ST_DataModelBindings.pro.ext_stub as ext_dict
import ST_DataModelBindings.bas.img_stub as img_dict
import ST_DataModelBindings.bas.dqc.ext_stub as extdqc_dict
import ST_DataModelBindings.bas.ppr.par_stub as par_dict
import ST_DataModelBindings.dpd.ext.calibratedframe_stub as ext_calibratedframe
import ST_DataModelBindings.dpd.ext.singleepochframe_stub as ext_singleepochframe
import ST_DataModelBindings.dpd.ext.stackedframe_stub as ext_stackedframe
try:
    import HeaderProvider.GenericHeaderProvider as HeaderProvider
except:
    import ST_DM_HeaderProvider.GenericHeaderProvider as HeaderProvider
#

#from dxml import dmUtils #ext dm utils

import ST_DM_DmUtils.DmUtils as dm_utils
import XMLproduct.dmUtils as dmUtils

def init():
    """Initializes the module.

    It adds some extra functionality to some Ext-stage2 data product binding
    classes.

    """
    # Add the extra functionality when the module is imported
    __add_binding_extra_functionality()


def create_single_epoch_frame(telescope, instrument, filter_name):
    """Creates an EXT-stage2 single epoch frame binding.

    Parameters
    ----------
    telescope: str
        The telescope name.
    instrument: str
        The instrument name.
    filter_name: str
        The filter name.

    Returns
    -------
    object
        The EXT-stage2 calibrated frame binding.

    """
    # Create the appropriate data product binding
    dp = ext_singleepochframe.DpdExtSingleEpochFrame()

    # Add the generic header to the data product
    dp.Header = HeaderProvider.create_generic_header("EXTDES")

    # Add the data element to the data product
    dp.Data = dmUtils.createGroundInstrumentImage(
        ext_dict.extSingleEpochFrame, telescope, instrument, filter_name)
    
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
    # -> not in the ext_dict
    #dp.Data.ObsCollectionDataStorage = dm_utils.create_fits_storage(
    #    ext_dict.extObsCollectionDataStorageFitsFile, "obsCollectionDataStorage.fits", "ext.obsCollection",
    #    "0.1")    
    dp.Data.FWHM = pyxb.BIND()

    # Add the QualityFlags element to the data product
    #dp.QualityFlags = __create_calibrated_frame_quality_flags()
    # Add the QualityFlags element to the data product
    dp.QualityFlags = create_quality_flags(
        extdqc_dict.sqfDpdExtCalibratedFrame)

    return dp


def create_calibrated_frame(telescope, instrument, filter_name):
    """Creates an EXT-stage2 calibrated frame binding.

    Parameters
    ----------
    telescope: str
        The telescope name.
    instrument: str
        The instrument name.
    filter_name: str
        The filter name.

    Returns
    -------
    object
        The EXT-stage2 calibrated frame binding.

    """
    # Create the appropriate data product binding
    dp = ext_calibratedframe.DpdExtCalibratedFrame()

    # Add the generic header to the data product
    dp.Header = HeaderProvider.create_generic_header("EXTDES")

    # Add the data element to the data product
    dp.Data = dm_utils.create_single_image(
        ext_dict.extCalibratedFrame, telescope, instrument, filter_name)
    
    # Add the optional ImageArea
    dp.Data.ImgArea = pyxb.BIND()
    dp.Data.ImgArea.Name = "  "
    
    # Add the ProcessParams
    dp.Data.ProcessParams = par_dict.extCalibratedFrameParameters()

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

   # Add the QualityFlags element to the data product
    #dp.QualityFlags = __create_calibrated_frame_quality_flags()
    # Add the QualityFlags element to the data product
    dp.QualityFlags = create_quality_flags(
        extdqc_dict.sqfDpdExtCalibratedFrame)


    return dp


def create_quality_flags(binding_class, flags_names=None):
    """Creates a static quality flags binding.

    Parameters
    ----------
    binding_class: class
        The static quality flags binding class.
    flags_names: list, optional
        The names of the quality flags that should be added. Default is None,
        which means all flags.

    Returns
    -------
    object
        The static quality flags binding.

    """
    # Add all the quality flags if they are not given
    if flags_names is None:
        flags_names = get_binding_class_element_names(binding_class)

    # Create the quality flags binding
    quality_flags = binding_class.Factory()

    # Fill the quality flags with default values
    for flag_name in flags_names:
        # Create the quality flag
        setattr(quality_flags, flag_name, pyxb.BIND(False, flagged=False))
        flag = getattr(quality_flags, flag_name)

        # Set the values of the optional and fixed attributes
        flag.xpath = "XXXX" if flag.xpath is None else flag.xpath

        if hasattr(flag, "minLimit"):
            flag.minLimit = -999 if flag.minLimit is None else flag.minLimit

        if hasattr(flag, "maxLimit"):
            flag.maxLimit = -999 if flag.maxLimit is None else flag.maxLimit

    return quality_flags


def get_binding_class_element_names(binding_class):
    """Returns the names of the elements and attributes inside a given binding
    class.

    Parameters
    ----------
    binding_class: class
        The binding class.

    Returns
    -------
    list
        A python list with the sorted names of the elements and attributes
        inside the binding class.

    """
    # Get the binding class name and super classes names
    class_names = [binding_class.__name__]
    super_class = binding_class.__bases__[0]

    while "euclid" in super_class._Name():
        class_names.append(super_class.__name__)
        super_class = super_class.__bases__[0]

    # Loop over the binding class attributes and extract the element names
    element_names = []

    for attribute_name in dir(binding_class):
        for class_name in class_names:
            if attribute_name.startswith("_" + class_name):
                element_names.append(attribute_name.split("__")[1])

    # Remove possible duplicates
    element_names = list(set(element_names))

    # Sort the element names
    element_names.sort()

    return element_names






def __create_calibrated_frame_quality_flags():
    """Creates an EXT calibrated frame quality flags binding.

    Returns
    -------
    object
        The EXT calibrated frame quality flags binding.

    """
    # Create the EXT calibrated frame quality flags binding
    quality_flags = extdqc_dict.sqfDpdExtCalibratedFrame()

    # Fill it with some default information
    #quality_flags.RootMeanSquareMinExceeded = pyxb.BIND(
    #    False, xpath="//Data/QualityParams/RootMeanSquare", minLimit=10, flagged=False)    
    #quality_flags.RootMeanSquareMaxExceeded = pyxb.BIND(
    #    False, xpath="//Data/QualityParams/RootMeanSquare", maxLimit=50, flagged=False)
    quality_flags.RootMeanSquareMinExceeded = pyxb.BIND(
        False, minLimit=10, flagged=False)    
    quality_flags.RootMeanSquareMaxExceeded = pyxb.BIND(
        False, maxLimit=50, flagged=False)
    quality_flags.NumberOfReferencePairingsMinExceeded = pyxb.BIND(
        False, minLimit=15, flagged=False)
    quality_flags.NumberOfReferencePairingsMaxExceeded = pyxb.BIND(
        False, maxLimit=1000, flagged=False)
    quality_flags.SigmaDeltaRightAscensionMaxExceeded = pyxb.BIND(
        False, maxLimit=50, flagged=False)
    quality_flags.SigmaDeltaDeclinationMaxExceeded = pyxb.BIND(
        False, maxLimit=50, flagged=False)
    quality_flags.MeanDeltaRightAscensionMaxExceeded = pyxb.BIND(
        False, maxLimit=50, flagged=False)
    quality_flags.MeanDeltaDeclinationMaxExceeded = pyxb.BIND(
        False, maxLimit=50, flagged=False)
    quality_flags.MeanOffsetMaxExceeded = pyxb.BIND(
        False, maxLimit=0.05, flagged=False)
    quality_flags.StandardDeviationMaxExceeded = pyxb.BIND(
        False, maxLimit=0.05, flagged=False)
    quality_flags.ResidualsSlopeXMinExceeded = pyxb.BIND(
        False, minLimit=-0.01, flagged=False)
    quality_flags.ResidualsSlopeXMaxExceeded = pyxb.BIND(
        False, maxLimit=0.01, flagged=False)
    quality_flags.ResidualsSlopeYMinExceeded = pyxb.BIND(
        False, minLimit=-0.01, flagged=False)
    quality_flags.ResidualsSlopeYMaxExceeded = pyxb.BIND(
        False, maxLimit=0.01, flagged=False)
    quality_flags.ResidualsSlopeRMinExceeded = pyxb.BIND(
        False, minLimit=-0.01, flagged=False)
    quality_flags.ResidualsSlopeRMaxExceeded = pyxb.BIND(
        False, maxLimit=0.01, flagged=False)
    quality_flags.NumberDensityUsedCalibratorsMinExceeded = pyxb.BIND(
        False, minLimit=200, flagged=False)
    quality_flags.NumberDensityUsedValidatorsMinExceeded = pyxb.BIND(
        False, minLimit=200, flagged=False)
    quality_flags.FullWidthHalfMaximumMinExceeded = pyxb.BIND(
        False, minLimit=0.3, flagged=False)
    quality_flags.FullWidthHalfMaximumMaxExceeded = pyxb.BIND(
        False, maxLimit=2.0, flagged=False)

    return quality_flags


def create_stacked_frame(telescope, instrument, filter_name):
    """Creates an EXT-stage2 stacked frame binding.

    Parameters
    ----------
    telescope: str
        The telescope name.
    instrument: str
        The instrument name.
    filter_name: str
        The filter name.

    Returns
    -------
    object
        The EXT-stage2 stacked frame binding.

    """
    # Create the appropriate data product binding
    dp = ext_stackedframe.DpdExtStackedFrame()

    # Add the generic header to the data product
    dp.Header = HeaderProvider.create_generic_header("EXTKIDS")

    # Add the data element to the data product
    dp.Data = dm_utils.create_stacked_image(
        ext_dict.extStackedFrame, telescope, instrument, filter_name)

    # Add the ProcessParams
    dp.Data.ProcessParams = par_dict.extStackedFrameParameters()

    # Add the QualityParmas
    dp.Data.QualityParams = ext_dict.extStackedFrameDqp()

    # Add the DataStorage
    dp.Data.DataStorage = dm_utils.create_fits_storage(
        ext_dict.extStackedFrameFitsFile, "data.fits", "ext.stackedFrame",
        "0.1")

    # Add the PsfModelStorage
    dp.Data.PsfModelStorage = dm_utils.create_fits_storage(
        ext_dict.extPsfModelFitsFile, "psfModel.fits", "ext.psfModel", "0.1")

    return dp


def __add_ext_extra_functionality(dp_binding_class):
    """Adds some extra functionality to a EXT-stage2 data product binding
    class.

    Parameters
    ----------
    dp_binding_class: class
        The EXT-stage2 data product binding class.

    """
    # Check that it contains a Data attribute
    if not hasattr(dp_binding_class, "Data"):
        print("The provided binding should contain a Data element.\n"
              "No extra functionality added.")
        return

    # Create a temporal binding instance
    dp_instance = dp_binding_class()

    # Create the Data element
    dp_instance.Data = pyxb.BIND()
    data_element = dp_instance.Data

    # Add the methods to the data product class
    dp_binding_class.is_stack = __is_stack

    if hasattr(data_element, "DataStorage"):
        dp_binding_class.set_data = __set_data
        dp_binding_class.get_data = __get_data
    if hasattr(data_element, "PsfModelStorage"):
        dp_binding_class.set_psf_model = __set_psf_model
        dp_binding_class.get_psf_model = __get_psf_model


def __is_stack(self):
    """Checks if the product is of stack type.

    Returns
    -------
    bool
        True if the product is a stack.

    """
    return isinstance(self, ext_stackedframe.dpdExtStackedFrame)


def __set_data(self, file_name):
    """Sets the data fits file name.

    Parameters
    ----------
    file_name: str
        The fits file name.

    """
    self.Data.DataStorage.DataContainer.FileName = file_name


def __get_data(self):
    """Returns the data fits file name.

    Returns
    -------
    str
        The data fits file name.

    """
    return self.Data.DataStorage.DataContainer.FileName


def __set_psf_model(self, file_name):
    """Sets the PSF model fits file name.

    Parameters
    ----------
    file_name: str
        The fits file name.

    """
    self.Data.PsfModelStorage.DataContainer.FileName = file_name


def __get_psf_model(self):
    """Returns the PSF model fits file name.

    Returns
    -------
    str
        The PSF model fits file name.

    """
    return self.Data.PsfModelStorage.DataContainer.FileName


def __add_binding_extra_functionality():
    """Adds some extra functionality to some Ext-stage2 data product binding
    classes.

    """
    dp_list = [ext_calibratedframe.dpdExtCalibratedFrame,
               ext_stackedframe.dpdExtStackedFrame]

    for dp in dp_list:
        dm_utils.add_extra_functionality(dp)
        __add_ext_extra_functionality(dp)
