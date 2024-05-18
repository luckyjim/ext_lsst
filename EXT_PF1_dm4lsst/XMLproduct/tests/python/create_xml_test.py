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
File: tests/python/create_xml_test.py

Created on: 06/29/20
Author: user
"""

import py.test
from XMLproduct.create_xml import *
import pyxb
from contextlib import contextmanager
import os.path


@contextmanager
def pretty_pyxb_exceptions():

    '''Context to raise pretty PyXB exceptions instead of cryptic ones.

    To be used with "with" like:

    with pretty_pyxb_exceptions():

        product.toDOM().pretty_xml()

    '''
    try:
        yield

    except pyxb.IncompleteElementContentError as e:
        print(e.details())
        raise e
   
    except pyxb.SimpleTypeValueError as e:
        print('\n+++++++++++++SimpleTypeValueError+++++++++++++++++++++')
        print(e.details())
        print('++++++++++++++++++++++++++++++++++')
        raise e
    except pyxb.SimpleContentAbsentError as e:
        print('\n+++++++++++++SimpleContentAbsentError+++++++++++++++++++++')
        print(f'Instance {e.instance} with simple content '
        f'was not provided with a value.')
        print('++++++++++++++++++++++++++++++++++')
        raise e
    except pyxb.UnprocessedKeywordContentError as e:
        print('\n+++++++++++++UnprocessedKeywordContentError+++++++++++++++++++++')
        print(e.details())
        print('++++++++++++++++++++++++++++++++++')
        raise e
    except pyxb.MissingAttributeError as e:
        print('\n+++++++++++++MissingAttributeError+++++++++++++++++++++')
        print(e.details())
        print(f'Instance {e.instance} with simple content '
        f'was not provided with a value.')
        print('++++++++++++++++++++++++++++++++++')
        raise e
    except pyxb.UnboundElementError as e:
        print('\n+++++++++++++UnboundElementError+++++++++++++++++++++')
        # print(e.details())
        print(f'Instance {e.instance} with simple content '
        f'was not provided with a value.')
        print('++++++++++++++++++++++++++++++++++')
        raise e
    except pyxb.ValidationError as e:
        print('\n+++++++++++++MissingAttributeError+++++++++++++++++++++')
        print(e.details())
        print(f'Instance {e.instance} with simple content ')
        print('++++++++++++++++++++++++++++++++++')
        raise e
        
    except pyxb.SimpleFacetValueError as e:
        print(f'\nViolated facet is: {e.facet}')
        raise e

    
class Testcreate_xml(object):
    """
    @class Testcreate_xml
    @brief Unit Test class
    """

    def test_bbox(self):
        poly = [[9.897055390092202, -18.75946680049288],
              [9.895983465004214, -18.90883519001286],
              [10.212465108028743, -18.90979249229269],
              [10.213092160709051, -18.760367919018012]]
        center = [10.054539777969973, -18.834585939240693]
        mbb = BBox()
        mbb.from_footprint(poly, center)
        print(mbb)
        
        
    def test_SingleEpochFrameLSSTDefault_create(self):
        n_target = 'lsst_default.xml'
        with pretty_pyxb_exceptions():            
            o_xml = SingleEpochFrameLSSTDefault()
            o_xml.save_xml_product(n_target)            
        assert os.path.exists(n_target)
    
    def test_SingleEpochFrameLSSTDefault_fill01(self):
        f_detrended = "/home/user/Work/Projects/EUC_EXT_DPDEXTDETRENDEDFRAME_LSST-1961147-R31-S21_20190721T002950.890Z_LSST-SWF1-2007.fits"
        f_fake = "fake.fits"
        with pretty_pyxb_exceptions():            
            o_xml = SingleEpochFrameLSST()
            o_xml.fill_xml_with_product(f_detrended, f_fake, f_fake, tag='_JMC')
            o_xml.save()          
        assert os.path.exists(o_xml.f_xml)
    
