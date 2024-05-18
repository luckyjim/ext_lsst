'''
Created on Mar 18, 2021
'''
import pprint
import file_conf as fc


def test_params_init():
    obj = fc.FileConfLSSTProcessing('exa.txt')
    if obj.check:
        pprint.pprint(obj.__dict__)


def test_params_is_ok():
    obj = fc.FileConfLSSTProcessing('exa.txt')
    print(obj)    
    if obj.is_ok():
        print(obj)
    
    
if __name__ == '__main__':
    #test_params()
    test_params_is_ok()
