#!/usr/bin/env python3
from SimUtils import file_utils as fu
import sys

TEST_PATH = './'

# user and password
AUTH_STRING = 'ZXVjbGlkMDEwOkV1Y2wxZDJTaGFyMw=='
WEBDAV_URL = 'https://shared.euclid.pic.es'

# Move test data to local folder
wd = fu.WebdavRepo('webdav', WEBDAV_URL, auth_string=AUTH_STRING)
fs = fu.FileSystemRepo('file_system', TEST_PATH)
fm = fu.FileMover(wd, fs)

        
ORIGIN_DIR = '/SC456/EXT/LSST/T2/output/data'
DEST_DIR = './'
lfile = wd.get_file_info(ORIGIN_DIR)
with open("raw_image_list.txt", "a") as mfile:
    for ifile in lfile:
        namef = ifile[0]
        mfile.write(namef+'\n')    
    
#fm.update_repo_content(ORIGIN_DIR, DEST_DIR, recursive=True)
