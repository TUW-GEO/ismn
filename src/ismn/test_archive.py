# -*- coding: utf-8 -*-

from ismn.archive import IsmnRoot
from tempfile import mkdtemp

def test_archive(path):
    archive = IsmnRoot(path)
    cont = archive.scan()
    csvs = archive.find_files('FMI/SAA111')

def test_zip_archive(path):
    archive = IsmnRoot(path)
    cont = archive.scan()
    fmi_csvs = archive.find_files('FMI/SAA111')
    all_csvs = archive.find_files(None, '*.csv')

    tempd = mkdtemp()

    archive.extract_file(fmi_csvs[0], tempd)
    archive.extract_dir('FMI/SAA111', tempd)

if __name__ == '__main__':
    test_archive(r"C:\Temp\delete_me\ismn\testdata_ceop")
    test_zip_archive(r"C:\Temp\delete_me\ismn\testdata_ceop.zip")