#! /usr/bin/env python
'''
Created on May 3, 2019

@author: colley
'''

"""
## source /cvmfs/euclid-dev.in2p3.fr/CentOS7/EDEN-2.1/bin/activate
"""

import sys
import process_image_lsst as pil
import file_conf as fc


file_par = sys.argv[1]
print(file_par )
my_pars= fc.FileConfLSSTProcessing(file_par)
if not my_pars.is_ok():
    print(' ')
    print(my_pars)
    raise

oProcess = pil.ProcessStage1LSST()
oProcess.init_with_file_param(my_pars)
oProcess.init_env_var()
ret = oProcess.euclid_processing()
exit(ret)