# -*- coding: utf-8 -*-
from ismn.readers import IsmnDataFile

def usecase_file(archive):
    filepath = "FMI/SAA111/FMI_FMI_SAA111_sm_0.050000_0.050000_5TE_20101001_20201005.stm"
    nodat = IsmnDataFile(archive, filepath, load_data=False)
    dat = IsmnDataFile(archive, filepath, load_data=True)
    nodat.load_data()
    print(nodat.data)
    print(dat.data)

def usecase_zip(archive):
    filepath = "FMI/SAA111/FMI_FMI_SAA111_sm_0.050000_0.050000_5TE_20101001_20201005.stm"
    nodat = IsmnDataFile(archive, filepath, load_data=False)
    dat = IsmnDataFile(archive, filepath, load_data=True)
    nodat.load_data()
    print(nodat.data)
    print(dat.data)

if __name__ == '__main__':
    usecase_zip( r"C:\Temp\delete_me\ismn\testdata_ceop.zip")
    usecase_file( r"C:\Temp\delete_me\ismn\testdata_ceop")