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
from tempfile import TemporaryDirectory
import os
import glob
import ismn.readers as readers
import numpy as np
import logging
from ismn.readers import extract_from_archive
from collections import OrderedDict
from functools import lru_cache
from tempfile import gettempdir
from shutil import rmtree

# @lru_cache(1)
def scan_archive(ismn_data_path, station_subdirs=True):
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

class MetaCollector(object):
    """ Read metadata for files downloaded from the ISMN """
    # default values, if there is no csv file available or it crashes for e.g saturation

    # attributes template for metadata that should be taken from data files
    # name, default value, dtype
    _file_metadata_template = np.array([
        ('network', '', object),
        ('station', '', object),
        ('variable', '', object),
        ('depth_from', np.nan, np.float),
        ('depth_to', np.nan, np.float),
        ('sensor', '', object),
        ('longitude', np.nan, np.float),
        ('latitude', np.nan, np.float),
        ('elevation', np.nan, np.float),
        ('root_path', '', object),
        ('sub_path', '', object),
        ('filename', '', object),
    ])

    # attributes template for metadata that should be taken from csv files
    # name, default value, dtype
    _csv_metadata_template = np.array([
        ('lc_2000', np.nan, object),
        ('lc_2005', np.nan, object),
        ('lc_2010', np.nan, object),
        ('lc_insitu', '', object),
        ('climate_KG', '', object),
        ('climate_insitu', '', object),
        ('saturation', np.nan, object),
        ('clay_fraction', np.nan, object),
        ('sand_fraction', np.nan, object),
        ('silt_fraction', np.nan, object),
        ('organic_carbon', np.nan, object),
        ])

    def __init__(self, data_path, meta_path, temp_root=gettempdir(), verbose=False):
        # unzip_strategy = station : extract per station before reading

        if not os.path.exists(data_path):
            raise FileNotFoundError(f"{data_path} does not exist")

        self.data_path = data_path

        if zipfile.is_zipfile(self.data_path):
            self.from_zip = True
        else:
            self.from_zip = False

        self.meta_path = meta_path

        # actual tempdir created (and removed) in existing tempdir later
        self.temp_root = temp_root

        self.verbose = verbose # log to print

    def _reset_dirs(self):
        """
        Initialise dirs for meta extraction. Delete any existing teo dirs.
        """
        if os.path.exists(self.meta_path):
            rmtree(self.meta_path)

        os.makedirs(self.meta_path)
        os.makedirs(self.temp_root, exist_ok=True)

    def get_path_info(self, filepath):
        """ Path info to strore in metadata"""
        if self.from_zip:
            sub_path = os.path.relpath(os.path.dirname(filepath), self.temp_root)
            return (os.path.dirname(self.data_path),
                    os.path.join(*sub_path.split(os.path.sep)[1:]),
                    os.path.basename(filepath))
        else:
            sub_path = os.path.relpath(os.path.dirname(filepath), self.data_path)
            return (self.data_path,
                    sub_path,
                    os.path.basename(filepath))

    def get_station_meta(self, station_path):
        files = os.listdir(station_path)
        # read additional metadata from csv file
        filename_csv = glob.glob('{}/*.csv'.format(station_path))

        station_meta = []

        if len(filename_csv) > 0:
            csv_meta_keys, csv_meta_entries = \
                self.meta_from_csv(os.path.join(station_path, filename_csv[0]))

            for filename in files:
                if filename.endswith('.stm'):
                    try:
                        file_meta_keys, file_meta_entries = \
                            self.meta_from_datafile(os.path.join(station_path, filename))
                    except (readers.ReaderException, IOError) as e:
                        continue

                    for var_meta_entry in file_meta_entries:
                        meta = var_meta_entry + csv_meta_entries
                        station_meta.append(tuple(meta))
        else:
            if any(filename.endswith('.stm') for filename in files):
                logging.info('No csv file available ({})'.format(station_path))

        return station_meta

    def meta_from_csv(self, filepath):
        """
        Read csv specific meta data,
        i.e., landcover info, climate info, soil info

        Parameters
        ----------
        filepath : str
            Name of the csv file to read.

        Returns
        -------
        metadata_keys : list
            List of variable names in metadata entries
        metadata_entries : list
            list contains the metadata info one sensor per row
        """
        try:
            csv_meta = readers.get_metadata_from_csv(filepath, as_dict=True)
        except:
            logging.info(f"Error occured while reading metadata from csv file "
                         f"({filepath})")
            csv_meta = dict(zip(self._csv_metadata_template[:,0],
                                self._csv_metadata_template[:,1]))

        return list(csv_meta.keys()), list(csv_meta.values())

    def meta_from_datafile(self, filepath):
        """
        function walks the data_path directory and looks for network
        folders and ISMN datafiles. It collects metadata for every
        file found and returns a numpy.ndarray of metadata

        Parameters
        ----------
        filepath : str
            Path to data file to read metadata from.

        Returns
        -------
        metadata_keys : list
            List of variable names in metadata entries
        metadata_entries : list
            list contains the metadata info one sensor per row
        """

        metadata_entries = []

        metadata = readers.get_metadata(filepath)
        metadata_keys = None

        for i, variable in enumerate(metadata['variable']):
            root_path, sub_path, filename = self.get_path_info(filepath)
            variable_meta = np.array([
                             ('network', metadata['network']),
                             ('station', metadata['station']),
                             ('variable', variable),
                             ('depth_from', metadata['depth_from'][i]),
                             ('depth_to', metadata['depth_to'][i]),
                             ('sensor', metadata['sensor']),
                             ('longitude', metadata['longitude']),
                             ('latitude', metadata['latitude']),
                             ('elevation', metadata['elevation']),
                             ('root_path', root_path),
                             ('sub_path', sub_path),
                             ('filename', filename),
                            ])

            if metadata_keys is None:
                metadata_keys = variable_meta[:,0].tolist()

            metadata_entries.append(variable_meta[:,1].tolist())

        return metadata_keys, metadata_entries

    def collect_from_archive(self):
        """
        Collect all metadata from a downloaded ISMN archive. Either the zip
        directly (slower) or the directory containing the extracted zip files.

        Returns
        -------
        metadata_catalog : np.array
            Structured array that contains all the metadata.
        """
        self._reset_dirs()

        logging.basicConfig(filename=os.path.join(self.meta_path, 'metadata.log'),
                            level=logging.DEBUG)

        dtype = list(map(tuple, self._file_metadata_template[:, (0,2)])) + \
                list(map(tuple, self._csv_metadata_template[:, (0,2)]))

        metadata_catalog = []

        archive = scan_archive(self.data_path)
        for network_dir, station_dirs in archive.items():
            print(network_dir)
            for station_dir in station_dirs:
                print(station_dir)
                if self.from_zip:
                    with TemporaryDirectory(prefix='ismn', dir=self.temp_root) as tempdir:
                        extract_from_archive(self.data_path, station_dir, tempdir)
                        station_path = os.path.normpath(os.path.join(tempdir, station_dir))
                        station_meta = self.get_station_meta(station_path)
                        metadata_catalog += station_meta
                else:
                    station_path = os.path.join(self.data_path, station_dir)
                    station_meta = self.get_station_meta(station_path)
                    metadata_catalog += station_meta

        return np.array(metadata_catalog, dtype=dtype)

if __name__ == '__main__':
    data = r"C:\Temp\delete_me\hawaii\Data_separate_files_20100101_20201008_5712_q2h0_20201008.zip"
    col = MetaCollector(data, meta_path=r"C:\Temp\delete_me\hawaii\meta")
    meta = col.collect_from_archive()
    print(meta)

    #col.get_station_meta(os.path.join('FMI', 'SOD140'))
