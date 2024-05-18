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
File: python/XMLproduct/script_createXML.py

Created on: 05/02/19
Author: Colley/Rosset
"""

from __future__ import division, print_function
import sys
if sys.version_info[0] < 3:
    from future_builtins import *

import argparse
import ElementsKernel.Logging as log
import glob
import os
import re
import XMLproduct.create_xml as xmlp


from os import listdir
from os.path import isfile, join

RE_TAG = re.compile(r'(?P<name>.*_\d{8}T\d{6}\.\d{3}Z)(?P<tag>.*)(?P<ext>\.fits)')


def replace_tag(fitsname, new_tag):
    m = RE_TAG.match(fitsname)
    if not m:
        return fitsname
    return ''.join([m.group('name'), new_tag, m.group('ext')])


def check_directory(value):    
    if not os.path.exists(value):
        raise argparse.ArgumentTypeError("%s isn't a valid directory" % value)
    return value


def defineSpecificProgramOptions():
    """
    @brief Allows to define the (command line and configuration file) options
    specific to this program

    @details
        See the Elements documentation for more details.
    @return
        An  ArgumentParser.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--path-fits', help='path to FITS files.',
            type=check_directory, required=True)
    parser.add_argument('-t', '--tag', help='tag name to fill PipelineRun XML element in Header.',
            default="PipelineRun_xx")
    parser.add_argument('-r', '--rel', help='tag name to fill DataSetRelease XML element in Header.',
            default="DataSetRelease_xx")
    parser.add_argument('-f', '--fil', help='filter position (RA,DEC,RADIUS) ',
            default="")
    #
    # !!! Write your program options here !!!
    # e.g. parser.add_argument('--string-value', type=str, help='A string option')
    #

    return parser


def extract_id_ccd(name):
    tag_lsst = name.split('_')[3].split('-')
    return tag_lsst[2] + '-' + tag_lsst[3]


def mainMethod(args):
    """
    @brief The "main" method.
    @details
        This method is the entry point to the program. In this sense, it is
        similar to a main (and it is why it is called mainMethod()).
    """

    logger = log.getLogger('script_createXML')
    logger.info('# Entering script_createXML mainMethod()')
    mypath = args.path_fits
    # ========================================
    # Delete previously created XML files
    # ========================================
    cmd = f"rm -f {mypath}/*.xml"
    print('rm old xml: ', cmd)
    os.system(cmd)

    # set list of FITS files
    only_fits = [f for f in listdir(mypath) if (isfile(join(mypath, f)) and '.fits' in f)]

    # Update tag product in fits file names if needed
    # JMC vori discusion sur gitlab APC4EuclidEXT
    flag_update_tag_fits = False
    if flag_update_tag_fits:
        fits_list = []
        for fitsname in only_fits:
            newfitsname = replace_tag(fitsname, args.tag)
            if newfitsname != fitsname:
                os.rename(os.path.join(mypath, fitsname),
                          os.path.join(mypath, newfitsname))
            fits_list.append(newfitsname)
    
        # Use the updated fits names...
        only_fits = fits_list

    l_frame = [f for f in only_fits if 'DPDEXTDETRENDEDFRAME' in f]
    l_cat = [f for f in only_fits if 'DPDEXTSOURCECATALOG' in f]
    l_psf = [f for f in only_fits if 'DPDEXTPSF' in f]
    
    d_frame = {extract_id_ccd(f):f for f in l_frame}
    d_cat = {extract_id_ccd(f):f for f in l_cat}
    d_psf = {extract_id_ccd(f):f for f in l_psf}
    
#     print("image:\n",d_frame)
#     print("cat:\n",d_cat)
#     print("psf:\n",d_psf)
    cpt_ccd_ok = 0
    o_xml = xmlp.SingleEpochFrameLSST(args.fil)	
    for ccd in d_frame:        
        # check if all FITS are present
        b_all_exist = (ccd in d_psf) and (ccd in d_cat)
        if  not b_all_exist:
            logger.warning(f'Skip {ccd} missed psf or catallog FITS')
            continue
        logger.info(f'Processing {ccd}')
        cpt_ccd_ok += 1
        f_frame = os.path.join(mypath, d_frame[ccd])
        f_psf = os.path.join(mypath, d_psf[ccd])
        f_cat = os.path.join(mypath, d_cat[ccd])
        if o_xml.fill_xml_with_product(f_frame, f_psf, f_cat, tag=args.tag, release=args.rel):
            # filter None or False
            o_xml.save()
    finale_mes = f'xml processing {cpt_ccd_ok} CCD on {len(d_frame)}.'
    if cpt_ccd_ok != len(d_frame):
        logger.warning(finale_mes)
    else:
        logger.info(finale_mes)
