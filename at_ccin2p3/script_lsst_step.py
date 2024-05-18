#! /usr/bin/env python
"""
##source /cvmfs/sw.lsst.eu/linux-x86_64/lsst_distrib/w_2018_31/loadLSST.bash
##setup lsst_distrib
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

if my_pars.data_in_sps_or_eas == 'sps':
    oProcess = pil.ProcessStage1LSST()
else:
    oProcess = pil.ProcessStage1LsstEasIn()
oProcess.init_with_file_param(my_pars)
oProcess.init_env_var()
ret = oProcess.lsst_processing()
exit(ret)
