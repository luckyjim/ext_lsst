'''
Created on Jun 4, 2021

@author: colley
'''

from process_image_lsst import *
from file_conf import FileConfLSSTProcessing


pars_test = FileConfLSSTProcessing('pars4test.txt')


def test_get_all_xml():
    mbut = ManageSetButlers(pars_test.root_butlers)
    l_sef = mbut.get_all_xml()
    print(len(l_sef))
    print(l_sef[:10])
    
    
def test_load_balancing_eas():
    eas_obj = EasTools(pars_test, True)
    mbut = ManageSetButlers(pars_test.root_butlers)
    mbut.load_balancing_eas(5, eas_obj)
    
    
def test_div_list():
    ml = [1,2,3,4,5,6,7,8,9,10]
    div_list(ml, 3)
    
    
    
    
if __name__ == '__main__':
    #test_get_all_xml()
    test_load_balancing_eas()
    #test_div_list()