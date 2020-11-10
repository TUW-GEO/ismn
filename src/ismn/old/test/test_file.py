# -*- coding: utf-8 -*-
from ismn.readers import DataFile, StaticMetaFile

def usecase_static_meta_file_zip(archive):
    filepath = 'FMI/SAA111/FMI_FMI_SAA111_static_variables.csv'
    f = StaticMetaFile(archive, filepath)
    meta = f.read_metadata()

def usecase_file(archive):
    filepath = "FMI/SAA111/FMI_FMI_SAA111_sm_0.050000_0.050000_5TE_20101001_20201005.stm"
    nodat = DataFile(archive, filepath, load_data=False)
    dat = DataFile(archive, filepath, load_data=True)
    nodat.load_data()
    print(nodat.data)
    print(dat.data)

    dat = DataFile(archive, filepath, load_data=True, static_meta=nodat.static_meta)
    dat.load_data()
    print(dat.data)
    print(dat.metadata)

def usecase_zip(archive):
    filepath = "FMI/SAA111/FMI_FMI_SAA111_sm_0.050000_0.050000_5TE_20101001_20201005.stm"
    nodat = DataFile(archive, filepath, load_data=False)
    dat = DataFile(archive, filepath, load_data=True)
    nodat.load_data()
    meta = dat.get_formatted_metadata()
    print(nodat.data)
    print(dat.data)
    assert dat.metadata == nodat.metadata
    print(dat.metadata)

    dat = DataFile(archive, filepath, load_data=True, static_meta=nodat.static_meta)
    dat.load_data()
    print(dat.data)
    print(dat.metadata)



if __name__ == '__main__':
    usecase_file( r"C:\Temp\delete_me\ismn\testdata_ceop")
    usecase_static_meta_file_zip(r"C:\Temp\delete_me\ismn\testdata_ceop.zip")
    usecase_zip( r"C:\Temp\delete_me\ismn\testdata_ceop.zip")

