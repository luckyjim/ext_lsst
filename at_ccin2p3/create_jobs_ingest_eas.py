#! /usr/bin/env python3

import process_image_lsst as pli
import eas
import argparse
import os.path 
import file_conf as fc

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
    			description='Create a SGE script array job to ingest SEF in EAS with a slow number of jobs')
    parser.add_argument('-f', '--file-param',
    			help='file parameters described images and star catalogue for LSST processing',
    			type=argparse.FileType('r'), required=True)
    args = parser.parse_args()
    # add attribut .name with file type
    my_pars = fc.FileConfLSSTProcessing(args.file_param.name)
    if not my_pars.check:
        print(' ') 
        print("========> STOP py.py")
        exit()
    if not my_pars.is_ok():
        print(' ')    
        print("========> STOP py.py")
        exit()
    print('Finale file parameters:\n========')
    print(my_pars)
    print('========')
    # test exclude option
    # print(args)
    # print(args.butler)
    # deal nb_proc
    eas_obj = eas.EasToolsToIngest(my_pars, True)    
    mbut = pli.ManageSetButlers(my_pars.root_butlers)
    mbut.load_balancing_eas(my_pars.nb_proc_eas, eas_obj)
    if mbut.nb_sef_job > 0:
        opli = pli.ProcessStage1LSST()
        opli.init_with_file_param(my_pars)
        opli.create_script_sge_for_eas(mbut.nb_sef_job)
    else:
        print('Load balancing FAILED.')
