# -*- coding: utf-8 -*-

"""
Module description
"""
# TODO:
#   (+) 
#---------
# NOTES:
#   -
from ismn.readers import IsmnFileCollection

def test_missing_csv():
    path = r"C:\Temp\delete_me\problems\hiwater"
    ds = IsmnFileCollection(path)
    print(ds.files)

def test_unknown_format():
    path = r"C:\Temp\delete_me\problems\risma_neu\Data_separate_files_19630813_20201015_5712_tJsl_20201015.zip"
    ds = IsmnFileCollection(path)
    print(ds.files)

def test_duplicate_depth():
    path = r"C:\Temp\delete_me\problems\sasmas"
    ds = IsmnFileCollection(path)
    print(ds.files)

if __name__ == '__main__':
    test_duplicate_depth()