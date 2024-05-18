'''
@author: colley
'''
import sys
import os
import shlex, subprocess
from astropy.io import fits
import file_conf as fc
import glob
import threading
import time
import psutil
import numpy as np

class EasTools(object):
    
    def __init__(self):
        self.server_old = "https://eas-dps-cus.test.euclid.astro.rug.nl"
        #self.server = "https://eas-dps-cus-ops.esac.esa.int"
        self.server = "https://eas-dps-cps-ops.esac.esa.int"                       
        self.server_cc = "https://cceuclid1.in2p3.fr:443" 
        self.verbose = False
        
    def set_params(self, params):
        self.params = params
        
            
    def _subprocess_cmd(self):
        cmd = self._cmd
        if self.verbose:
            print(cmd)
        args = shlex.split(cmd)        
        self._proc = subprocess.run(args, stdout=subprocess.PIPE)
        try:
            std_out = self._proc.stdout.decode("utf-8")
        except:
            std_out = ""
        return std_out
    
    
    def cmd_curl_download(self, file_name, dest_path="", use_nl=False):
        '''
        return True if process return 0
        
        :param file_name:
        :param dest_path:
        '''
        if use_nl:
            server = self.server
        else:
            server = self.server_cc
        path_file = os.path.join(dest_path, file_name)
        self._cmd = f'curl -k --netrc-file {self.params.file_curl_eas}'
        self._cmd += f' "{server}/{file_name}" -o {path_file}'
        ret = self._subprocess_cmd()
        return self._proc.returncode == 0
              
                    
    def cmd_curl(self, prod, query, fields):
        '''
        return "" or list of string
        
        :param prod:
        :param query:
        :param fields:
        '''
        pjt = self.params.eas_project
        self._cmd = f'curl "{self.server}/COORS?class_name={prod}&{query}&project={pjt}&fields={fields}"'
        ret = self._subprocess_cmd()
        if ret != "":
            ret = ret.split("\n")[1:]
            ret.remove('')
        return ret
            
                
    
    def cmd_retrieve(self, prod, query):
        '''
        
        :param prod:
        :param query:
        '''
        '''
        python dataProductRetrieval_SC8.py --username=EC_JMCOLLEY --password=/pbs/home/c/colley/mypwd_EAS.txt --project=TEST 
        --data_product=DpdTrueUniverseOutput  --query="Header.DataSetRelease=LSST_SC8_703_R1&Data.EuclidPointingId=353398&Header.PipelineDefinitionId=SIM_EXT-LSST"
        '''
        script = self.params.eas_retrieve_script
        user = self.params.user_eas
        pwd = self.params.file_pwd_eas
        pjt = self.params.eas_project
        self._cmd = f"python {script} --username={user} --password={pwd} --project={pjt}"
        self._cmd += f'--data_product={prod} --query="{query}"'
        return self._subprocess_cmd()
    

class EasToolsToDownload(EasTools):
    
    def __init__(self):
        super().__init__()
    
    def query_list_raw_image(self, data_set_release):
        self.list_raw_image = ""
        query = f"Header.DataSetRelease={data_set_release}&Header.ManualValidationStatus.ManualValidationStatus!=INVALID"
        fields = "Data.DataStorage.DataContainer.FileName.DataFileName"
        std_out = self.cmd_curl("DpdExtLssSingleVisitImage", query, fields)
        if std_out == "":
            print(f"'{self._cmd}' is nok")
            return None
        self.list_raw_image = std_out        
        return self.list_raw_image
    
    def query_definition_id(self, raw_file):
        query = f"Data.DataStorage.DataContainer.FileName.DataFileName={raw_file}"
        fields = "Header.PipelineDefinitionId"
        std_out = self.cmd_curl("DpdExtLssSingleVisitImage", query, fields)
        if std_out == "":
            print(f"'{self._cmd}' is nok")
            return None                
        return std_out[0]

    def query_name_stars_catalog(self, release, id_visit, def_id):
        query = f"Header.DataSetRelease={release}&"
        query += f"Data.EuclidPointingId={id_visit}&"
        query += f"Header.PipelineDefinitionId={def_id}"
        fields = "Data.StarCatalogFitsFile.DataContainer.FileName.DataFileName"
        std_out = self.cmd_curl("DpdTrueUniverseOutput", query, fields)
        if std_out == "":
            print(f"'{self._cmd}' is nok")
            return None
        return std_out[0]


class EasToolsToIngest(EasTools):
    
    def __init__(self, params, verbose=False):
        super().__init__()
        self.set_params(params)        
        self.verbose = verbose
    
    def get_list_sef_in_eas(self):
        release = self.params.data_set_release
        project = self.params.eas_project
        cmd = f'curl "{self.server}/COORS?class_name=DpdExtSingleEpochFrame'
        cmd += f'&Header.DataSetRelease={release}&project={project}'
        cmd += f'&fields=Header.ProductId.ObjectId"'
        if self.verbose:
            print(cmd)
        args = shlex.split(cmd)
        # out_proc = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        out_proc = subprocess.run(args, stdout=subprocess.PIPE)
        try:
            xml_files = out_proc.stdout.decode("utf-8").split('\n')[1:]
        except:
            xml_files = []
        return xml_files
    
    def sef_in_eas(self, xml_file):
        # remove .xml
        xml_file_no_ext = xml_file.split('.xml')[0]
        release = self.params.data_set_release
        project = self.params.eas_project
        cmd = f'curl "{self.server}/COORS?class_name=DpdExtSingleEpochFrame'
        cmd += f'&Header.DataSetRelease={release}&project={project}&Header.ProductId.ObjectId'
        cmd += f'={xml_file_no_ext}&fields=Header.ProductId.ObjectId"'
        if self.verbose:
            print(cmd)
        args = shlex.split(cmd)
        out_proc = subprocess.run(args, stdout=subprocess.PIPE)
        out_dec = out_proc.stdout.decode("utf-8")
        print(out_dec)
        if xml_file_no_ext in out_dec:
            return True
        else:
            return False    
    
    
