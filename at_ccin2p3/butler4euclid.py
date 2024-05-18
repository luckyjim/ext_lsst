"""

Create bulter (LSST repository) and insert star catalog in lsst format

Colley JM APC/IN2P3/CNRS

"""

import sys
import shlex, subprocess
import os.path 
import argparse
from stack_for_stage1 import LsstStackForStage1


class CreateButlerFromEuclidSimu(object):
    '''
    Create butler directory from template and convert, ingest reference stars catalog
    '''
    
    def __init__(self, stack4s1):
        '''
        
        :param stack4s1: 
        '''
        assert issubclass(type(stack4s1), LsstStackForStage1)
        self.stack4s1 = stack4s1

    def init_file(self, pathname_butler, file_cat):
        '''Create butler and ingest ref star catalag
        
        :param pathname_butler: 
        :param file_cat: path to ref star catalog in LSST format
        '''
        self.pathname_butler = os.path.abspath(pathname_butler)
        self.ref_cat = file_cat
        self.cat_temp = self.pathname_butler + '/cat_temp.txt'
        self.verbose = 2
        self.status = -1
    
    def convert_in_lsst_format_from_flux(self): 
        '''Convert ref star catalog in LSST format (from euclid format since 2021)
        
        Since march 2021 True Univers star catalog doesn't contain LSST magnitude columns       
        
        :param name_cat:
        '''
        self.stack4s1.convert_star_cat(self.ref_cat, self.cat_temp)
    
    def remove_temp_cat(self): 
        '''remove ref star catalog        
        '''
        os.system("rm -rf " + self.cat_temp)
     
    def insert_cat(self):
        '''Insert ref star catalog in butler        
        '''
        if not os.path.exists(self.cat_temp):
            print("No temp catalog !")
            return        
        command = self.stack4s1.ingest_star_cat_cmd(self.cat_temp)
        args = shlex.split(command)
        out_proc = subprocess.run(args, stdout=subprocess.PIPE)
        self.status = out_proc.returncode
        if self.verbose >= 1: 
            print("import butler : ", out_proc.returncode == 0)
        if self.verbose >= 2:
            output = out_proc.stdout                        
            print(output.decode("utf-8"))
    
    def copy_butler_template(self):
        '''Copy template directory of butler        
        '''
        path_butler = self.stack4s1.path_to_butler_template()
        os.system(f"cp -r {path_butler} {self.pathname_butler}")
    
    def create_insert_ref_cat(self):
        self.copy_butler_template()
        os.chdir(self.pathname_butler)
        self.remove_temp_cat()
        self.convert_in_lsst_format_from_flux()
        self.insert_cat()            
        # self.remove_temp_cat()
        