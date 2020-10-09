# The MIT License (MIT)
#
# Copyright (c) 2020 TU Wien
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

import zipfile
import tempfile
import shutil
import os
import numpy as np

def extract_subdir(zip, subdir, out_path):
    """
    Extract all files in a zip file under the passed subdir

    Parameters
    ----------
    zip : str
        Path to ismn zipfile
    subdir : str
        No leading /!
        Subdir to extract
    """
    with zipfile.ZipFile(zip) as zi:
        filelist = zi.namelist()
        for file in filelist:
            if os.path.split(os.path.commonpath([file, subdir]))\
                    == os.path.split(subdir):
                zi.extract(file, out_path)


def unzip_file(path_in_zipfile, path_to_zipfile):
    """
    function creates temporary folder and extracts file from zip archive

    Parameters
    ----------
    path_in_zipfile : string
        relative path of zipped file (inside archive)
    path_to_zipfile : string
        absolute path of zip-archive
    Returns
        -------
        tmp_directory : string
            temporary directory created to store extracted files
        tmp_filename : string
            temporary path to extracted file
    """
    tmp_directory = tempfile.mkdtemp()
    with zipfile.ZipFile(path_to_zipfile) as zi:
        zi.extract(path_in_zipfile, tmp_directory)
    return tmp_directory, os.path.join(tmp_directory, path_in_zipfile)




def zip_folder(zip_path):
    """
    function returns the directory of a zip file.

    Parameters
    ----------
    zip_path : string
        path of zip-file of downloaded ISMN data
    """
    if os.path.isdir(zip_path):
        return zip_path
    return os.path.dirname(zip_path)

def walk(zip_path, fun, tmp_directory):
    """
    function emulates the functionality of os.walk for a zipped archive.

    Parameters
    ----------
    zip_path : string
        path to archive
    fun : function
        function to be executed
    tmp_directory : string
        path to temporary directory

    Returns
    -------
    metadata_archive : list of tuples
        list contains the metadata info one sensor per row
    """
    metadata_archive = []
    with zipfile.ZipFile(zip_path) as zi:
        file_list = zi.namelist()
        path_files_list = [(os.path.split(x)[0], x) for x in file_list]
        sub_folders = [x[0] for x in path_files_list]
        unique_folders = np.unique(sub_folders)
        for sub in unique_folders:
            sub_filepaths = [item[1] for item in path_files_list if item[0] == sub]
            for filename in sub_filepaths:
                zi.extract(filename, tmp_directory)
            sub_filenames = [os.path.split(filename)[1] for filename in sub_filepaths]
            metadata = fun(os.path.join(tmp_directory, sub), sub_filenames,
                           os.path.join(zip_folder(zip_path), sub))
            if len(metadata):
                metadata_archive.extend(metadata)

    return metadata_archive


def take_walk(zip_path, fun):
    """
    function creates a temporary directory, calls the walk function and deletes the temporary directory.

    Parameters
    ----------
    zip_path : string
        absolute path of zip-file where ISMN data os stored

    Returns
    -------
    metadata_catalog : list of tuples
        list contains the metadata info one sensor per row
    """
    tmp_directory = tempfile.mkdtemp()
    metadata_catalog = walk(zip_path, fun, tmp_directory)
    shutil.rmtree(tmp_directory)
    return metadata_catalog






