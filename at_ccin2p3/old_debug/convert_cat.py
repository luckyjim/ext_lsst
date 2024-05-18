#! /usr/bin/env python3

from SIM_EXT import catalog, starcatalog, galcatalog

import sys
import os

"""
. /cvmfs/euclid.in2p3.fr/CentOS7/EDEN-2.0/etc/profile.d/euclid.sh
 
export PYTHONPATH=$PYTHONPATH:$PWD/SIM_EXT/SIM_EXT/python
python3 convert_cat.py 

"""

#input_cat_file = "/sps/euclid/Users/colley/lsst/simuT2/input/cat/EUC_SIM_TUSTARCAT-100007_20180823T124312.991Z_SC456-EXT-C7a_T2.fits"
#input_cat_file = "/sps/euclid/Users/colley/lsst/wddata/EUC_SIM_TUSTARCAT-100175_20180823T124734.334Z_SC456-EXT-C7a_T2.fits"
input_cat_file = sys.argv[1]
output_cat_file = sys.argv[2]
os.system('rm -rf '+ output_cat_file)
input_cat = starcatalog.read_starcatalog(input_cat_file)
#(ou input_cat = galcatalog.read_galcatalog(input_cat_file))
catalog.write_catalog(input_cat, output_cat_file, format="lsst")
#catalog.write_lsststack_catalog(input_cat, output_cat_file, format="lsst")
