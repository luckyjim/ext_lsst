'''
Goal of the module: Handle different LSST stack version with same interface

Solution: design pattern factory

@author: colley
'''

import os.path
import ref_catalog_euclid2lsst as rce


def factory_for_stage1_stack(image_singu):
    '''
    Crude implementation of factory with mane of the singularity LSST stack image.
    
    :param image_singu: path/name of singularity LSST stack image
    '''
    name = os.path.basename(image_singu)
    if "2018" in name:
        return LsstStackForStage1_w2018()
    else:
        return LsstStackForStage1_w2021()
    

class LsstStackForStage1(object):
    
    def is_ref_cats(self, butler):
        '''
        
        :param butler:
        :return: True if reference cat is already in butler
        '''
        path_cat = os.path.join(butler, 'input/ref_cats')
        return os.path.exists(path_cat)

    def path_to_butler_template(self):
        pass

    def ingest_image_cmd(self, raw_file):
        pass
    
    def process_image_cmd(self):
        pass
    

class LsstStackForStage1_w2018(LsstStackForStage1):
    '''
    process image with processEImage, in old format without segment 
    '''
        
    def path_to_butler_template(self):
        return "/sps/euclid/OU-EXT/lsst/butler/butler_template_2018"
    
    def convert_star_cat(self, input_cat, output_cat_file):
        return rce.ref_cat_conversion_2018(input_cat, output_cat_file)
        
    def ingest_star_cat_cmd(self, star_cat):
        cmd = f'ingestReferenceCatalog.py input {star_cat}'
        cmd += ' --configfile IngestIndexedReferenceTask_DC1.py'
        return cmd

    def ingest_image_cmd(self, raw_file):
        cmd = f'ingestSimImages.py input {raw_file}'
        return cmd 
    
    def process_image_cmd(self, visit, raft, sensor):
        cmd = 'processEimage.py input --output output'
        cmd += f' --id visit={visit} raft="{raft}" sensor="{sensor}"'
        cmd += ' --configfile processConfig.py'
        return cmd   

    
class LsstStackForStage1_w2021(LsstStackForStage1):
    """
    process image with processCdd, in format with 16 segments
    """
    
    def path_to_butler_template(self):
        return "/sps/euclid/OU-EXT/lsst/butler/butler_template_2021"
            
    def convert_star_cat(self, input_cat, output_cat_file):
        '''
        3 changes  compared to the 2018 version
          * y => Y
          * isresolved => resolved
          * value "resolved" for star is now O (?!), before 1 with "isresolved" 
          
        :param input_cat:
        :param output_cat_file:
        '''
        return rce.ref_cat_conversion_2021(input_cat, output_cat_file)
        
    def ingest_star_cat_cmd(self, star_cat):
        cmd = f'ingestReferenceCatalog.py input {star_cat}'
        cmd += ' --configfile IngestIndexedReferenceTask.py'
        cmd += ' --clobber-config'
        return cmd

    def ingest_image_cmd(self, raw_file):
        cmd = f'ingestImages.py input {raw_file} --mode link --ignore-ingested'
        return cmd
    
    def process_image_cmd(self, visit, raft, sensor):
        '''
        
        :param visit:
        :param raftName
        :param detectorName:
        '''        
        cmd = 'processCcd.py input --rerun processCcdOutputs'
        #cmd += f' --id visit={visit} raft={raft} sensor={sensor}'
        #cmd += f" --id visit='{visit}' ccd='95'"
        cmd += ' --id visit='{visit}' raftName='{raft}' detectorName='{sensor}' --configfile processConfig.py'
        cmd += ' --longlog'
        cmd += ' --show data'
        return cmd
 
