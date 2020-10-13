# -*- coding: utf-8 -*-

import os
import numpy as np
import zipfile
from collections import OrderedDict
from typing import Sequence


    # @lru_cache(1)
    def scan_zip_archive(self, station_subdirs=True):
        """ Go through an ismn archive (extracted or zipped) and group
         station folders in archive by network folders
         """
        cont = {}
        if zipfile.is_zipfile(ismn_data_path):
            with zipfile.ZipFile(ismn_data_path) as zi:
                file_list = zi.namelist()
                for f in file_list:
                    relpath = os.path.split(f)[0]
                    if relpath == '': continue
                    net, stat = os.path.split(relpath)
                    if station_subdirs:
                        stat = relpath
                    if net not in cont.keys():
                        cont[net] = np.array([])
                    if not np.isin([stat], cont[net])[0]:
                        cont[net] = np.append(cont[net], stat)
        else:
            for f in os.scandir(ismn_data_path):
                if f.is_dir():
                    if f.path == ismn_data_path:continue
                    net = os.path.relpath(f.path, ismn_data_path)
                    if net == '': continue
                    for stat in os.scandir(f.path):
                        if stat.is_dir():
                            if net not in cont.keys():
                                cont[net] = np.array([])
                            if station_subdirs:
                                cont[net] = np.append(cont[net], os.path.join(net, stat.name))
                            else:
                                cont[net] = np.append(cont[net], stat.name)

        return OrderedDict(sorted(cont.items()))



def extract_from_archive(zip, subdir, out_path):
    """
    Extract all files in a zip file under the passed subdir

    Parameters
    ----------
    zip : str
        Path to ismn zipfile
    subdir : str
        No leading /!
        Subdir to extract, or a single file in the zip archive to extract
        Or a template for files to extract.
    out_path : str
        Path where the extracted file(s) is/are stored

    Returns
    -------
    extracted : Sequence
        Path(s) to the extracted file(s)
    """
    with zipfile.ZipFile(zip) as zi:
        filelist = np.array(zi.namelist())

        if subdir in filelist: # single file was passed
            zi.extract(member=subdir, path=out_path)

            extracted = [os.path.realpath(os.path.join(out_path, subdir))]

        else: # subdir was passed
            filterlist = list(filter(lambda x: x.startswith(subdir), filelist)).copy()
            zi.extractall(members=list(filterlist), path=out_path)

            extracted = [os.path.realpath(os.path.join(out_path, f))
                    for f in filterlist]

    return extracted

if __name__ == '__main__':
    tempdir = "/tmp/tmp1tkf_s8e"
    path = "/media/wolfgang/Windows/temp/ismndata/testdata_ceop.zip"
    subdir = 'FMI/SOD130/FMI_FMI_SOD130_sm_0.100000_0.100000_CS655-B_20101001_20201005.stm'
    ds = extract_from_archive(path, subdir, tempdir)