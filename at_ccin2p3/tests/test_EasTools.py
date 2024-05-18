'''
Created on Jun 3, 2021

@author: colley
'''

from process_image_lsst import *
from file_conf import FileConfLSSTProcessing


pars_test = FileConfLSSTProcessing('pars4test.txt')

eas_in = EasToolsToDownload()
eas_in.verbose = True
eas_in.set_params(pars_test)


def test_get_list_sef_in_eas():
    eas = EasToolsToIngest(pars_test, True)
    l_xml = eas.get_list_sef_in_eas()
    print(len(l_xml))
    print(l_xml[0])                        
    print(l_xml[1])
    print(l_xml[2])
    print(l_xml[3])    
    #print(l_xml)
    
def test_sef_in_eas():
    eas = EasToolsToIngest(pars_test, True)
    ret = eas.sef_in_eas("EUC_EXT_DPDEXTSINGLEEPOCHFRAME_LSST-353363-R10-S00_20210527T102042.252Z")
    print(ret)
    ret = eas.sef_in_eas("EUC_EXT_DPDEXTSINGLEEPOCHFRAME_LSST-353363-R10-S00_20210527T102042.252")
    print(ret)
    
    
def proto_list_raw_image():
    eas_in = EasToolsToDownload()
    eas_in.verbose = True
    eas_in.set_params(pars_test)
    l_raw =eas_in.query_list_raw_image("LSST_SC8_703_R1") 
    if l_raw:
        print(l_raw[:100])
        print(len(l_raw))
        

def proto_def_id():
    raw_file = "EUC_SIM_LSST-R02S21-129227-SCIENCE-r_07CD30AD09E2-0006047_20210520T114045.578297Z_SC8_MAIN_EXT_703_1.fits"
    ret = eas_in.query_definition_id(raw_file)
    if ret:
        print(ret)

def proto_name_cat():
    def_id = "SIM_EXT-LSST"
    ret = eas_in.query_name_stars_catalog("LSST_SC8_703_R1", "129227", def_id)
    if ret:
        print(ret)
    

def proto_download():
    cat = "EUC_SIM_TUSTARCAT-129227_07CD30ACE1DB-0012004_20210520T083938.188591Z_SC8_MAIN_EXT_703_1.fits"
    eas_in.cmd_curl_download(cat)

if __name__ == '__main__':
    #test_get_list_sef_in_eas()
    #proto_list_raw_image()
    #proto_def_id()
    #proto_name_cat()
    proto_download()
    #test_sef_in_eas()