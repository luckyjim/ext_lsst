import process_image_lsst as pli

    

def test_load_balancing_image():
    os.environ["GE_TASK_ID"] = "1"
    os.environ["GE_TASK_LAST"] = "4"
    os.environ["LSST_ID_CAT"] = "100175"
    os.environ["PATH_LIST_IMAGE"] = "/sps/euclid/Users/colley/lsst/processT2/raw_image_list_head.txt"
    oProcess = pli.ProcessStage1_forLSST()    
    oProcess.select_and_loadbalancing()    
    print(oProcess.index_file)
    oProcess.rank = 1
    oProcess.select_and_loadbalancing()    
    print(oProcess.index_file)
    oProcess.rank = 2    
    oProcess.select_and_loadbalancing()    
    print(oProcess.index_file)    
    oProcess.rank = 3    
    oProcess.select_and_loadbalancing()    
    print(oProcess.index_file)
    
    
    
def test_my_list_file(): 
    """ add docstring    
    """   
    os.environ["PATH_LIST_IMAGE"] = "/sps/euclid/Users/colley/lsst/processT2/raw_image_list_head.txt"
    oProcess = pli.ProcessStage1LSST(lfile)    
    oProcess.select_and_loadbalancing()
    print(oProcess.my_list)
    oProcess.lsst_processing()



def test_script_sge():
    oProcess = pli.ProcessStage1LSST()
    oProcess.my_list = ["a","b","c",'d']
    print(oProcess.create_script_sge_for_lsst())
    

def test_m_select_raw_image():
    pro = pli.ProcessStage1LSST()
    #pro.file_name_image = "/sps/euclid/Users/colley/lsst/stage1/processT6/list_test.txt"
    pro.file_name_image = "/sps/euclid/Users/colley/lsst/stage1/processT6/list_raw_t6.txt"
    pro.select_raw_image()
    print(pro.l_select_images)
    return pro
    

def test_m_extract_all_visit_ids():
    pro = test_m_select_raw_image()
    pro.extract_all_visit_ids()
    print(pro.dict_visit_imag)
    return pro
    

def test_m_create_file_image_by_visit():
    pro = test_m_extract_all_visit_ids()
    pro.path_butler = "/sps/euclid/Users/colley/lsst/stage1/processT6/test"
    pro.create_file_image_by_visit()
    return pro


def test_m_create_link_butler():
    pro = test_m_create_file_image_by_visit()
    pro.create_link_butler()
    

def test_m_search_max_image_by_visit():
    pro = test_m_extract_all_visit_ids()
    pro.search_max_image_by_visit()
    print(pro.max_file)
    
    
if __name__ == '__main__':
    #test_script_sge()
    #test_m_select_raw_image()
    #test_m_extract_all_visit_ids()
    #test_m_create_file_image_by_visit()
    #test_m_create_link_butler()
    test_m_search_max_image_by_visit()
    