#! /usr/bin/env python3
'''
Created on Apr 16, 2019
'''

import configparser
import os.path


def stringbool_2_bool(str_bool):
    s_ok = str_bool.lower() in ['true', 'false']
    return s_ok, str_bool.lower() == 'true'


class FileConfLSSTProcessing(object):
    '''
    '''

    def __init__(self, p_nfile="NA"):
        '''
        '''        
        # 0: name attrib,  1: description, 2: example, 3:default value
        #  defautl value "$TBD" means mandatory must to be defined     
        self.def_pars = [ ["file_list_images", "File content name of all raw image to process. TEMP : In EAS mode write a dummy file", "raw_image_list.txt", "$TBD"],
                          ["file_list_star_cat", "file content name of all stars ref catalog for each visit", "cat_list.txt", "TBD.txt"],
                          ["root_butlers", "full path of root directory contents butler directory for each visit", "/path/to/run/directory", "$TBD"],
                          ["root_script", "full path of EXT_LSST_Testing git package used to process data", "/path/to/run/src", "$TBD"],
                          ["singularity_lsst_stack", "path/name of singularity LSST stack image, with name format w_YYYY_NN.sif", "~/w_2018_31.sif", "$TBD"],
                          ["pipeline_run", "value set in PipelineRun XML element, start by _", "_RUN_SC8_2", "_TBD"],
                          ["data_set_release", "DataSetRelease stage 1 out", "SC8_SWF1", "TBD"],
                          ["data_in_sps_or_eas", "Indicate origin of raw images, from sps with string 'sps' or DataSetRelease of raw images in EAS", "SIMU_XXX_R1", "sps"],                                                
                          ["filter_pos", "Parameters to select image a region of sky: <ra,dec,radius> in degree ", "227.8,76.4,2.0", ""],
                          ["skip_proc", "Skip LSST stack processing. Type: [true, false]", "false", "false"],
                          ["skip_xml", "Skip xml creating. Type: [true, false]", "false", "false"],
                          ["skip_ingest", "Skip eas ingestion. Type: [true, false]", "false", "false"],
                          ["del_fits", "delete fits files at the begining of processing. Type: [true, false]", "true", "true"],
                          ["user_eas", "User name for the EAS account. Type: string", "colley", "$TBD"],
                          ["file_pwd_eas", "Full path with name of file contents password EAS. Type : string path", "/home/colley/eas_pwd.txt", "$TBD"],
                          ["file_curl_eas", "Full path with name of file contents server, login and password EAS. Type : string path", "/home/colley/eas_pwd.txt", "$TBD"],
                          ["eas_project", "Name of Euclid data base. Type : [TEST, EUCLID]", "TEST", "TEST"],
                          ["eas_ingest_script", "path and script to use to ingest data in EAS. Type : string path", "/my/ingest_client.py", "/sps/euclid/OU-EXT/lsst/eas/ingest_client_sc8.py"],
                          ["eas_retrieve_script", "path and script to use to retrieve data from EAS. Type : string path", "/my/dataProductRetrieval.py", "/sps/euclid/OU-EXT/lsst/eas/dataProductRetrieval_SC8.py"],
                          ["sdc_eas", "Name of SDC for EAS storage. Type : string ", "SDC-NL", "SDC-FR-PROD"],
                          ["nb_proc_ingest", "Number of process in parallel to ingest data in EAS. Type : integer ", 10, 10],
                          ["use_old_files", "During ingest EAS step, I can use FITS already in EAS with true. Type : [true, false]", "false", "false"]
                          ]
    
    
        if p_nfile == "NA":
            return
        
        self.key_pars = [param[0] for param in self.def_pars]
        # loop on params in file
        config = configparser.ConfigParser()
        config.read(p_nfile)
        self.nfile = os.path.abspath(p_nfile)
        self.nfile_short = os.path.basename(p_nfile)
        
        for key in config['DEFAULT']:
            if key in self.key_pars:
                # print(f"param '{key}' is in list")
                setattr(self, key, config['DEFAULT'][key])
            else:
                print(f"====>  Unknow parameter '{key}', rejeted !!!")
        self.check = True
        
        # check presence
        for pars in self.def_pars:
            key = pars[0]
            if pars[3] == "$TBD":
                # mandatory key
                if not hasattr(self, key):
                    print(f"You must defined '{key}' parameters, check NOK")
                    self.check = False
            else:
                # optional
                if not hasattr(self, key):
                    setattr(self, key, pars[3])   
        # manage boolean type
        syntax, val_bool = stringbool_2_bool(self.skip_proc)
        if not syntax:
            print(self.skip_proc)
            self.check = False
        else:
            self.skip_proc = val_bool
        # use_old_files
        syntax, val_bool = stringbool_2_bool(self.skip_xml)
        if not syntax:
            self.check = False
        else:
            self.skip_xml = val_bool

        syntax, val_bool = stringbool_2_bool(self.skip_ingest)
        if not syntax:
            self.check = False
        else:
            self.skip_ingest = val_bool

        syntax, val_bool = stringbool_2_bool(self.del_fits)
        if not syntax:
            self.check = False
        else:
            self.del_fits = val_bool

        syntax, val_bool = stringbool_2_bool(self.use_old_files)
        if not syntax:
            self.check = False
        else:
            self.use_old_files = val_bool
            
        # convert to integer
        self.nb_proc_eas = int(self.nb_proc_ingest)
    
    def _check_file(self, key):
        m_file = getattr(self, key)
        ret = os.path.exists(m_file)
        if not ret:
            print(f'[{key}] ERROR "{m_file}" not exist.')
        return ret
    
    def is_ok(self):
        """
        ret_final = ret_final and self._check_file("")
        """
        if not self.check:
            return False
        ret_final = self._check_file("file_pwd_eas")        
        ret_final = ret_final and self._check_file("file_list_images")
        # problem: in /sps case file_list_star_cat mandotary but not in EAS case .... so
        # ret_final = ret_final and self._check_file("file_list_star_cat")
        ret_final = ret_final and self._check_file("root_script")
        ret_final = ret_final and self._check_file("eas_ingest_script")
        ret_final = ret_final and self._check_file("singularity_lsst_stack")                     
        return ret_final
    
    def __str__(self):
        mstr = ""
        for key in self.key_pars:
            mstr += f"\n{key} = {getattr(self, key)}"
        return mstr

    def example(self):
        expl = "#======================================================"
        expl += "\n#  File parameters example "
        expl += "\n#======================================================"        
        expl += "\n\n[DEFAULT]"
        for param in self.def_pars:
            expl += f"\n# {param[1]}"            
            # print(expl)
            if param[3] == "$TBD":
                expl += f"\n# Mandatory, example value : {param[2]}"
                expl += f"\n{param[0]} = TBD"
            else:
                expl += f"\n# Optional, default value is {param[3]}"
                expl += f"\n# {param[0]} = {param[2]}"
            expl += "\n"
        return expl
            
            
if __name__ == "__main__":
    obj = FileConfLSSTProcessing('NA')
    print(obj.example())
    
