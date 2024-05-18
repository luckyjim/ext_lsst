#! /usr/bin/env python3

import process_image_lsst as pli
import argparse
import os.path 
import file_conf as fc

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
    			description='Create a SGE script array job to process LSST image with 2 steps. Load balancing is 1 job by focal plan.')
    parser.add_argument('-f', '--file-param',
    			help='file parameters described images and star catalogue for LSST processing',
    			type=argparse.FileType('r'), required=True)
    parser.add_argument('-p', '--pattern',
                        help='pattern to select image name in file_list_images, by example id visit -p=-100106,-100175', default="*")
    parser.add_argument('-l', '--limit-file',
                        help='debugging option, reduce number image by proc at limit-file value', type=int, default=10000)
    args = parser.parse_args()
    # add attribut .name with file type
    my_pars = fc.FileConfLSSTProcessing(args.file_param.name)
    if not my_pars.check:
        print(' ') 
        print("========> STOP create_arrayjob_lsst.py")
        exit()
    if not my_pars.is_ok():
        print(' ')    
        print("========> STOP create_arrayjob_lsst.py")
        exit()
    print('Finale file parameters:\n========')
    print(my_pars)
    print('========')
    # test exclude option
    # print(args)
    # print(args.butler)
    # deal nb_proc
    opli = pli.ProcessStage1LSST()
    opli.limit_max_file = args.limit_file    
    if args.pattern == '*':
        opli.pattern = None
    else:
        opli.pattern = args.pattern.split(',')
        print(opli.pattern)        
    opli.init_with_file_param(my_pars)
    opli.select_and_loadbalancing()
    opli.create_script_sge_for_lsst()
    opli.create_script_sge_for_euclid()
