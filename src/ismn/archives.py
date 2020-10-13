# -*- coding: utf-8 -*-


import os
import numpy as np
import zipfile
from collections import OrderedDict
import glob
import fnmatch
from pathlib import Path

def _zip(func):
    def wrapper(cls, *args, **kwargs):
        if not cls.zip:
            raise IOError("Zip archive expected")
        return func(cls, *args, **kwargs)
    return wrapper

def _dir(func):
    def wrapper(cls, *args, **kwargs):
        if cls.zip:
            raise IOError("Unzipped archive expected")
        return func(cls, *args, **kwargs)
    return wrapper


class ISMNArchive():

    def __init__(self, path):
        """

        Parameters
        ----------
        path : str
            Path to the downloaded zip file or the extracted zip
            directory.
        """
        self.path = Path(path)

        if not self.path.exists():
            raise IOError(f'Archive does not exist: {self.path}')

        if zipfile.is_zipfile(path):
            self.zip = zipfile.ZipFile(self.path, mode='r')
        else:
            self.zip = None

        self.__cont = None

    def __contains__(self, filepath):
        """ Check if files exists in archive """
        if self.zip:
            return filepath in self.zip.filelist()
        else:
            path = self.path / filepath
            return path.exists()

    def _check_subpath(self, subpath):
        """ Check if subpath is a valid path """
        if self.zip:
            assert subpath[0] != '/', "No leading / in subpath!"
            assert "\\" not in subpath , "Backslash not allowed in zip"
        else:
            assert (self.path / Path(subpath)).exists(), \
            "Subpath does not exist in archive"

    @property
    def cont(self):
        """ Get cont if object, or scan to create cont. """
        if self.__cont is None:
            self.__cont = self.scan()
        return self.__cont

    @_zip
    def _scan_zip(self, station_subdirs:bool=True) -> OrderedDict:
        """
        Go through the archive and group station folders in archive
        by network folders
        """
        cont = {}

        file_list = self.zip.namelist()

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

        self.__cont = OrderedDict(sorted(cont.items()))

        return self.cont

    @_dir
    def _scan_dir(self, station_subdirs:bool=True) -> OrderedDict:
        """
        Go through the archive and group station folders in archive
        by network folders
        """
        cont = {}

        for f in os.scandir(self.path):
            if f.is_dir():
                if f.path == self.path:continue
                net = os.path.relpath(f.path, self.path)
                if net == '': continue
                for stat in os.scandir(f.path):
                    if stat.is_dir():
                        if net not in cont.keys():
                            cont[net] = np.array([])
                        if station_subdirs:
                            cont[net] = np.append(cont[net], os.path.join(net, stat.name))
                        else:
                            cont[net] = np.append(cont[net], stat.name)

        self.__cont = OrderedDict(sorted(cont.items()))

        return self.cont

    @_dir
    def _find_files_dir(self, subpath:str=None, fn_templ:str='*.csv') -> list:
        """
        Find files in the archive or a subdirectory of the archive
        that match to the passed filename template.
        """
        if subpath is None:
            subpath = '**'

        self._check_subpath(subpath)

        filenames = glob.glob(str(self.path / subpath / fn_templ))

        return filenames

    @_zip
    def _find_files_zip(self, subpath:str=None, fn_templ:str='*.csv') -> list:
        """
        Find files in zip archive that match the passed template and subdir.
        """
        if subpath is None:
            subpath = '**'

        self._check_subpath(subpath)

        all_files = np.array(self.zip.namelist())


        filterlist = \
            list(filter(lambda f: fnmatch.fnmatch(f, f"{subpath}/{fn_templ}"),
                        all_files)).copy()

        return filterlist

    def scan(self, station_subdirs=True) -> OrderedDict:
        """
        Go through archive (zip or dir) and group station folders

        Parameters
        ----------
        station_subdirs : bool, optional (default: True)
            Include the station dir as a subdir of network dir.
            If Flase is selected, station dir is included directly.

        Returns
        -------
        cont : OrderedDict
            Archive content, station dirs grouped by network dirs
        """

        if self.zip:
            return self._scan_zip(station_subdirs)
        else:
            return self._scan_dir(station_subdirs)

    def find_files(self, subpath=None, fn_templ='*.csv'):
        """
        List files in archive or a subdirectory of the archive that match
        the passed filename pattern.

        Parameters
        ----------
        subpath: str, optional (default: None)
            Use linux slashes /, no leading / to define a subpath
            in the zip file. If None is selected, the whole archive is searched.
        fn_templ : str, optional (default: '*.csv')
            Filename template for files that are searched in the passed dir.

        Returns
        -------
        files : list
            Found files that match the passed template.
        """

        if self.zip:
            return self._find_files_zip(subpath, fn_templ)
        else:
            return self._find_files_dir(subpath, fn_templ)

    @_zip
    def extract_file(self, path_in_archive, out_path):
        """
        Extract a file from the zip archive

        Parameters
        ----------
        path_in_archive : str
            Use linux slashes /, no leading / to define a subpath
            Relative path in the archive (network/station/filename)
        out_path : str
            Directory where the extracted file is stored.

        Returns
        -------
        extraced_file : str
            Path to the the extracted file.
        """
        out_path = Path(out_path)

        self._check_subpath(path_in_archive)

        filelist = np.array(self.zip.namelist())

        if path_in_archive in filelist: # single file was passed
            self.zip.extract(member=path_in_archive, path=out_path)

        extracted = out_path / path_in_archive

        return extracted

    @_zip
    def extract_dir(self, subdir_in_archive, out_path):
        """
        Extract all files in subdir from the zip archive

        Parameters
        ----------
        subdir_in_archive : str
            Relative path in the archive (network/station or network)
        out_path : str
            Directory where the extracted files are stored.

        Returns
        -------
        extraced_files : str
            Path to the the extracted files.
        """
        out_path = Path(out_path)

        self._check_subpath(subdir_in_archive)

        filelist = np.array(self.zip.namelist())

        filterlist = list(filter(lambda x: x.startswith(subdir_in_archive),
                                 filelist)).copy()

        self.zip.extractall(members=filterlist, path=out_path)

        return [out_path / f for f in filterlist]

def test_archive():
    path = "/media/wolfgang/Windows/temp/ismndata/testdata_ceop"
    archive = ISMNArchive(path)
    cont = archive.scan()
    csvs = archive.find_files('FMI/SAA111')


def test_zip_archive():
    from tempfile import mkdtemp

    path = "/media/wolfgang/Windows/temp/ismndata/testdata_ceop.zip"
    archive = ISMNArchive(path)
    cont = archive.scan()
    fmi_csvs = archive.find_files('FMI/SAA111')
    all_csvs = archive.find_files(None, '*.csv')

    tempd = mkdtemp()

    archive.extract_file(fmi_csvs[0], tempd)
    archive.extract_dir('FMI/SAA111', tempd)





if __name__ == '__main__':
    test_archive()
    test_zip_archive()
