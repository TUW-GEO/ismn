# -*- coding: utf-8 -*-
# The MIT License (MIT)
#
# Copyright (c) 2021 TU Wien
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import os
import numpy as np
import zipfile
from collections import OrderedDict
import glob
import fnmatch
import warnings
from pathlib import Path, PurePosixPath
from typing import Union, List


def zip(func):
    def wrapper(cls, *args, **kwargs):
        if not cls.zip:
            raise IOError("Zip archive expected, use @dir functions instead.")
        return func(cls, *args, **kwargs)

    return wrapper


def dir(func):
    def wrapper(cls, *args, **kwargs):
        if cls.zip:
            raise IOError("Unzipped archive expected, use @zip functions instead.")
        return func(cls, *args, **kwargs)

    return wrapper


class IsmnRoot:
    """
    Connection to the zip resp. extracted zip archive downloaded from the
    ismn website. This class only handles file access / requests made by the
    readers, lists files in path and can extract files to temp folders for
    safe reading.

    Attributes
    ----------
    path : Path
        Data directory
    """

    def __init__(self, path):
        """
        Parameters
        ----------
        path : str or Path
            Path to the downloaded zip file or the extracted zip
            directory.
        """
        self.path = Path(path)

        if not self.path.exists():
            raise IOError(f"Archive does not exist: {self.path}")

        self.__cont = None
        self.__isopen = False

        self.open()

    @property
    def isopen(self) -> bool:
        # if data is a zipfile, this indicates if the zip is opened
        return self.__isopen

    @isopen.setter
    def isopen(self, isopen):
        self.__isopen = isopen

    @property
    def root_dir(self) -> Path:
        # the parent directory where the data is stored
        if self.zip:
            return self.path.parent
        else:
            return self.path

    def __repr__(self):
        """Simplified representation of object as string"""
        __type = type(self)
        __zip = "Zip" if self.zip else "Unzipped"
        return f"{__type.__module__}.{__type.__qualname__} {__zip} at {self.path}"

    def __contains__(self, filepath) -> bool:
        """Check if files exists in archive"""
        if self.zip:
            filepath = PurePosixPath(filepath)
            return str(filepath) in self.zip.namelist()
        else:
            path = self.path / filepath
            return path.exists()

    def clean_subpath(self, subpath) -> Union[Path, PurePosixPath]:
        """Check if subpath is a valid path and adapt to archive format and os"""
        subpath = Path(subpath)
        if subpath.parts[0] in ["/", "\\"]:
            warnings.warn("Remove leading (back)slash in passed subpath.")
            subpath = Path(*subpath.parts[1:])

        if self.zip:
            subpath = PurePosixPath(subpath)
        else:
            assert (
                self.path / Path(subpath)
            ).exists(), "Subpath does not exist in archive"

        return subpath

    @property
    def cont(self):
        """Get cont of object, or scan to create cont."""
        if self.__cont is None:
            self.__cont = self.scan()
        return self.__cont

    @zip
    def __scan_zip(self, station_subdirs: bool = True) -> OrderedDict:
        """
        Go through the archive and group station folders in archive
        by network folders. if scan_files is selected, group files by
        station folders.
        """
        cont = {}

        file_list = self.zip.namelist()

        for f in file_list:
            relpath = os.path.split(f)[0]
            if relpath == "":
                continue
            net, stat = os.path.split(relpath)
            if net == "":
                continue
            if station_subdirs:
                stat = relpath
            if net not in cont.keys():
                cont[net] = np.array([])
            if not np.isin([stat], cont[net])[0]:
                cont[net] = np.append(cont[net], stat)

        self.__cont = OrderedDict(sorted(cont.items()))

        return self.cont

    @dir
    def __scan_dir(self, station_subdirs: bool = True) -> OrderedDict:
        """
        Go through the archive and group station folders in archive
        by network folders
        """
        cont = {}

        for f in os.scandir(self.path):
            if f.is_dir():
                if f.path == self.path:
                    continue
                net = os.path.relpath(f.path, self.path)
                if net == "":
                    continue
                for stat in os.scandir(f.path):
                    if stat.is_dir():
                        if net not in cont.keys():
                            cont[net] = np.array([])
                        if station_subdirs:
                            cont[net] = np.append(cont[net], Path(net, stat.name))
                        else:
                            cont[net] = np.append(cont[net], stat.name)

        self.__cont = OrderedDict(sorted(cont.items()))

        return self.cont

    @dir
    def __find_files_dir(self, subpath: str = None, fn_templ: str = "*.csv") -> list:
        """
        Find files in the archive or a subdirectory of the archive
        that match to the passed filename template.
        """
        if subpath is None:
            subpath = "**"

        subpath = self.clean_subpath(subpath)

        filenames = glob.glob(str(self.path / subpath / fn_templ))
        filenames = [Path(os.path.relpath(f, self.path)) for f in filenames]
        return filenames

    @zip
    def __find_files_zip(self, subpath: str = None, fn_templ: str = "*.csv") -> list:
        """
        Find files in zip archive that match the passed template and subdir.
        """
        if subpath is None:
            subpath = "**"

        subpath = self.clean_subpath(subpath)

        all_files = np.array(self.zip.namelist())

        filterlist = list(
            filter(
                lambda f: fnmatch.fnmatch(f, f"{subpath}/{fn_templ}"),
                all_files,
            )
        ).copy()

        return filterlist

    def scan(self, station_subdirs=True) -> OrderedDict:
        """
        Go through archive (zip or dir) and group station folders

        Parameters
        ----------
        station_subdirs : bool, optional (default: True)
            Include the station dir as a subdir of network dir.
            If False is selected, station dir is included directly.

        Returns
        -------
        cont : OrderedDict
            Archive content, station dirs grouped by network dirs
        """

        if self.zip:
            return self.__scan_zip(station_subdirs)
        else:
            return self.__scan_dir(station_subdirs)

    def find_files(self, subpath=None, fn_templ="*.csv"):
        """
        List files in archive or a subdirectory of the archive that match
        the passed filename pattern.

        Parameters
        ----------
        subpath: str, optional (default: None)
            Use linux slashes '/' (and no leading '/') to define a subpath
            in the zip file. If None is selected, the whole archive is searched.
        fn_templ: str, optional (default: '*.csv')
            Filename template for files that are searched in the passed dir.

        Returns
        -------
        files : list[str]
            Found files that match the passed template.
        """

        if self.zip:
            return self.__find_files_zip(subpath, fn_templ)
        else:
            return self.__find_files_dir(subpath, fn_templ)

    @zip
    def extract_file(self, file_in_archive, out_path):
        """
        Extract a file from the zip archive

        Parameters
        ----------
        file_in_archive : Path or str
            Use linux slashes '/' (no leading '/') to define a subpath
            Relative path in the archive (network/station/filename)
        out_path : Path or str
            Directory where the extracted file is stored.

        Returns
        -------
        extracted_file : Path
            Path to the the extracted file.
        """
        out_path = Path(out_path)
        file_in_archive = PurePosixPath(file_in_archive)

        file_in_archive = self.clean_subpath(file_in_archive)

        ls = np.array(self.zip.namelist())

        ext = None
        if str(file_in_archive) in ls:  # single file was passed
            ext = self.zip.extract(member=str(file_in_archive), path=out_path)

        return Path(ext)

    @zip
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

        subdir_in_archive = PurePosixPath(subdir_in_archive)
        subdir_in_archive = self.clean_subpath(subdir_in_archive)

        ls = np.array(self.zip.namelist())

        filterlist = list(
            filter(lambda x: x.startswith(str(subdir_in_archive)), ls)
        ).copy()

        self.zip.extractall(members=filterlist, path=out_path)

        return [out_path / f for f in filterlist]

    def open(self):
        # open connection to data archive
        if zipfile.is_zipfile(self.path):
            self.zip = zipfile.ZipFile(self.path, mode="r")
            self.name = self.path.with_suffix("").name
        else:
            self.zip = None
            self.name = self.path.name

        self.isopen = True

    def close(self):
        # close connection to data archive
        if self.zip:  # if not zip, connection is always open
            self.zip.close()
            self.zip = None
            self.isopen = False

    def __enter__(self):
        return self

    def __exit__(self, value_type, value, traceback):
        self.close()
