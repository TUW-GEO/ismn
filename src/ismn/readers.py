# The MIT License (MIT)
#
# Copyright (c) 2019 TU Wien
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
import io
import logging
import glob

import numpy as np
import pandas as pd

from collections import OrderedDict
from ismn.archive import IsmnArchive
from ismn.components import *

from tempfile import gettempdir, TemporaryDirectory
from pathlib import Path, PurePosixPath
from typing import Union
import warnings

import zipfile as zf

class IsmnFileCollection(object):

    """
    The IsmnFileCollection class reads and organized the metadata
    information of ISMN files.

    Parameters
    ----------
    path : str
        Root path of ISMN files.
    load_data : bool, optional
        If True data will be loaded during metadata reading.

    Attributes
    ----------
    data_path : str
        Root path of ISMN files.
    files : list
        List of ISMN filenames.

    Methods
    -------
    get_networks()
        Get networks from ISMN file collection.
    get_stations(network=None)
        Get stations from ISMN file collection.
    get_sensors(self, network=None, station=None)
        Get sensors from ISMN file collection.
    """

    def __init__(self, data_path, load_data=False):

        self.data_path = data_path
        self.files = {}
        self._build_filelist(load_data)


    def _build_filelist(self, load_data):
        """
        Scan ismn archive and build file object list.
        Reuse static metadata if possible
        """

        archive = scan_archive(self.data_path)

        for net_dir, stat_dirs in archive.items():
            for stat_dir in stat_dirs:
                root = os.path.join(self.data_path, stat_dir)
                static_meta = None

                for filename in glob.glob(os.path.join(root, '*.stm')):
                    f = IsmnFile(filename, load_data, static_meta)
                    static_meta = f.static_meta

                    if f['network'] not in self.files.keys():
                        self.files[f['network']] = {}
                    if f['station'] not in self.files[f['network']]:
                        self.files[f['network']][f['station']] = []

                    self.files[f['network']][f['station']].append(f)


    def get_networks(self):
        """
        Get networks from ISMN file collection.

        Returns
        -------
        networks : dict of empty networks
            Dict of networks.
        """
        networks = {}

        for n in self.files.keys():
            networks[n] = Network(n)

        return networks

    def get_stations(self, network=None):
        """
        Get stations from ISMN file collection.

        Parameters
        ----------
        network : str, optional
            Network name (default: None).

        Returns
        -------
        stations : dict
            Dict of stations.
        """
        stations = {}

        for net, stats in self.files.items():
            if network in [None, net]:
                for stat, files in stats.items():
                    for f in files:
                        if f['station'] not in stations:
                            stations[f['station']] = Station(f['station'],
                                                             f['longitude'],
                                                             f['latitude'],
                                                             f['elevation'])

        return stations

    def get_sensors(self, network=None, station=None):
        """
        Get sensors from ISMN file collection.

        Parameters
        ----------
        network : str, optional
            Network name (default: None).
        station : str, optional
            Station name (default: None).

        Returns
        -------
        sensors : dict
            Dict of sensors.
        """
        sensors = {}

        for net, stats in self.files.items():
            if network in [None, net]:
                for stat, files in stats.items():
                    if stat in [None, station]:
                        for f in files:
                            name = '{}_{}_{}_{}'.format(f['sensor'],
                                                        f['variable'],
                                                        f['depth'].start,
                                                        f['depth'].end)

                            if name not in sensors:
                                snr = Sensor(name, f['variable'],
                                             f['sensor'], f['depth'])
                                sensors[name] = snr

        return sensors

    def __repr__(self):
        """
        Print summary of ISMN file collection.
        """
        logger.info('Number of networks: {}'.format(len(self.get_networks())))
        logger.info('Number of stations: {}'.format(len(self.get_stations())))
        logger.info('Number of sensors: {}'.format(len(self.get_sensors())))
        logger.info('Number of files: {}'.format(len(self.files)))

class IsmnFile(object):
    """
    General base class for data and static metadata files in ismn archive.

    Parameters
    ----------
    archive: IsmnArchive or str
        Archive that contains the file to read
    file_path : Path or str
        Path to the file in the archive.
    temp_root : Path or str, optional (default : gettempdir())
        Root directory where a separate subdir for temporary files
        will be created (and deleted).
    """

    def __init__(self, archive, file_path, temp_root=gettempdir()):

        if not isinstance(archive, IsmnArchive):
            archive = IsmnArchive(archive)

        self.archive = archive
        self.file_path = self.archive._clean_subpath(file_path)

        if self.file_path not in self.archive:
            raise IOError(f'Archive does not contain file: {self.file_path}')

        if not os.path.exists(temp_root):
            os.makedirs(temp_root, exist_ok=True)

        self.temp_root = temp_root

class StaticMetaFile(IsmnFile):
    """
    Represents a csv file containing site specific static variables.
    These attributes shall be assigned to all sensors at that site.

    Parameters
    ----------
    archive: IsmnArchive
        Archive that contains the file to read
    file_path : Path or str
        Path to the file in the archive. No leading slash!
    temp_root : Path or str, optional (default : gettempdir())
        Root directory where a separate subdir for temporary files
        will be created (and deleted).
    """

    def __init__(self, archive, file_path, temp_root=gettempdir()):

        super(StaticMetaFile, self).__init__(archive, file_path, temp_root)

        if self.file_path.suffix.lower() != '.csv':
            raise IOError(f'CSV file expected for StaticMetaFile object')

    def _read_field(self, fieldname:str) -> np.array:
        """
        Extract a field from the loaded csv metadata
        """
        if fieldname in self.data.index:
            dt = list()

            for i, j in zip(np.atleast_1d(self.data.loc[fieldname]['depth_from[m]']),
                            np.atleast_1d(self.data.loc[fieldname]['depth_to[m]'])):
                dt.append(('{}m_{}m'.format(i, j), np.float))

            return np.array([tuple(np.atleast_1d(self.data.loc[fieldname]['value']))],
                            dtype=np.dtype(dt))
        else:
            return None

    def _read_csv(self, csvfile:Path) -> pd.DataFrame:
        """ Load static metadata data frame from csv """
        try:
            data = pd.read_csv(csvfile, delimiter=";")
            data.set_index('quantity_name', inplace=True)
        except:
            # set columns manually
            logging.info('no header: {}'.format(csvfile))
            data = pd.read_csv(csvfile, delimiter=";", header=None)
            cols = list(data.columns.values)
            cols[0:7] = ['quantity_name', 'unit', 'depth_from[m]', 'depth_to[m]',
                         'value', 'description', 'quantity_source_name']
            data.columns = cols
            data.set_index('quantity_name', inplace=True)

        return data

    def read_metadata(self):
        """
        Read csv file containing static variables into data frame.

        Returns
        -------
        data : OrderedDict
            Data read from csv file.
        """
        if self.archive.zip:
            with TemporaryDirectory(prefix='ismn', dir=self.temp_root) as tempdir:
                extracted = self.archive.extract_file(self.file_path, tempdir)
                self.data = self._read_csv(extracted)
        else:
            self.data = self._read_csv(self.archive.path / self.file_path)

        # read landcover classifications
        lc = self.data.loc[['land cover classification']][['value', 'quantity_source_name']]

        lc_dict = {'CCI_landcover_2000': np.nan,
                   'CCI_landcover_2005': np.nan,
                   'CCI_landcover_2010': np.nan,
                   'insitu': ''}

        for key in lc_dict.keys():
            if key in lc['quantity_source_name'].values:
                if key != 'insitu':
                    lc_dict[key] = np.int(lc.loc[lc['quantity_source_name']
                                                 == key]['value'].values[0])
                else:
                    lc_dict[key] = lc.loc[lc['quantity_source_name']
                                          == key]['value'].values[0]
                    logging.info(f'insitu land cover classification available: {self.file_path}')

        # read climate classifications
        cl = self.data.loc[['climate classification']][['value', 'quantity_source_name']]
        cl_dict = {'koeppen_geiger_2007': '', 'insitu': ''}
        for key in cl_dict.keys():
            if key in cl['quantity_source_name'].values:
                cl_dict[key] = cl.loc[cl['quantity_source_name'] == key]['value'].values[0]
                if key == 'insitu':
                    logging.info(f'insitu climate classification available: {self.file_path}')

        saturation = self._read_field('saturation')
        clay_fraction = self._read_field('clay fraction')
        sand_fraction = self._read_field('sand fraction')
        silt_fraction = self._read_field('silt fraction')
        organic_carbon = self._read_field('organic carbon')

        return OrderedDict([
            ('lc_2000', lc_dict['CCI_landcover_2000']),
            ('lc_2005', lc_dict['CCI_landcover_2005']),
            ('lc_2010', lc_dict['CCI_landcover_2010']),
            ('lc_insitu', lc_dict['insitu']),
            ('climate_KG', cl_dict['koeppen_geiger_2007']),
            ('climate_insitu', cl_dict['insitu']),
            ('saturation', saturation),
            ('clay_fraction', clay_fraction),
            ('sand_fraction', sand_fraction),
            ('silt_fraction', silt_fraction),
            ('organic_carbon', organic_carbon),
        ])

class IsmnDataFile(IsmnFile):

    """
    IsmnFile class represents a single ISMN data file.
    This represents only .stm data files not metadata csv files.

    Parameters
    ----------
    archive : IsmnArchive or str
        Archive to the downloaded data.
    file_path : str
        Path in the archive to the ismn file. No leading slash!
    load_data : bool, optional
        If True data will be loaded during metadata reading.
    todo: update.

    Attributes
    ----------
    filename : str
        Filename.
    file_type : str
        File type information (e.g. ceop).
    metadata : dict
        Metadata information.
    data : numpy.ndarray
        Data stored in file.

    Methods
    -------
    load_data()
        Load data from file.
    read_data()
        Read data in file.
    _read_metadata()
        Read metadata from file name and first line of file.
    _get_metadata_ceop_sep()
        Get metadata in the file format called CEOP in separate files.
    _get_metadata_header_values()
        Get metadata file in the format called Header Values.
    _get_metadata_from_file(delim='_')
        Read first line of file and split filename.
        Information is used to collect metadata information for all
        ISMN formats.
    _read_format_ceop_sep()
        Read data in the file format called CEOP in separate files.
    _read_format_header_values()
        Read data file in the format called Header Values.
    _read_csv(names=None, usecols=None, skiprows=0)
        Read data.
    """

    def __init__(self, archive, file_path, load_data=False, static_meta=None,
                 temp_root=gettempdir()):

        super(IsmnDataFile, self).__init__(archive, file_path, temp_root)

        self.file_type = 'undefined'
        self.metadata = {}
        self.data = None

        self.static_meta = static_meta

        self._read_metadata()

        if load_data:
            self.load_data()

    def __getitem__(self, item):
        return self.metadata[item]

    def load_data(self):
        """
        Load data from file.
        """
        if self.data is None:
            if self.file_type == 'ceop':
                # self._read_format_ceop()
                raise NotImplementedError
            elif self.file_type == 'ceop_sep':
                self._read_format_ceop_sep()
            elif self.file_type == 'header_values':
                self._read_format_header_values()
            else:
                logger.warning(f"Unknown file type: {self.file_path}")

    def read_data(self):
        """
        Read data in file.

        Returns
        -------
        data : pandas.DataFrame
            File content.
        """
        self.load_data()
        return self.data

    def _read_metadata(self):
        """
        Read metadata from file name and first line of file.
        """
        header_elements, filename_elements = self._get_metadata_from_file()

        if len(filename_elements) == 5 and len(header_elements) == 16:
            self.file_type = 'ceop'
            raise RuntimeError('CEOP format not supported')
        elif len(header_elements) == 15 and len(filename_elements) >= 9:
            self.metadata = self._get_metadata_ceop_sep()
            self.file_type = 'ceop_sep'
        elif len(header_elements) < 14 and len(filename_elements) >= 9:
            self.metadata = self._get_metadata_header_values()
            self.file_type = 'header_values'
        else:
            logger.warning(f"Unknown file type: {self.file_path} in {self.archive}")

        if self.static_meta is None:
            self.static_meta = self._get_static_metadata_from_csv()

        self.metadata.update(self.static_meta)

    def _get_static_metadata_from_csv(self):
        """
        Read static metadata from csv file in the same directory as the ismn
        data file.

        Returns
        -------
        static_meta : OrderedDict
            Dictionary of static metadata
        """
        csv = self.archive.find_files(self.file_path.parent, '*.csv')
        assert len(csv) == 1, f"Expected 1 csv file for station, found {len(csv)}"
        static_meta = StaticMetaFile(self.archive, csv[0]).read_metadata()

        return static_meta

    def _get_metadata_ceop_sep(self):
        """
        Get metadata in the file format called CEOP in separate files.

        Returns
        -------
        metadata : dict
            Metadata information.
        """
        header_elements, filename_elements = self._get_metadata_from_file()

        if len(filename_elements) > 9:
            sensor = '_'.join(filename_elements[6:len(filename_elements) - 2])
        else:
            sensor = filename_elements[6]

        if filename_elements[3] in variable_lookup:
            variable = variable_lookup[filename_elements[3]]
        else:
            variable = filename_elements[3]

        metadata = {'network': filename_elements[1],
                    'station': filename_elements[2],
                    'variable': variable,
                    'depth': Depth(float(filename_elements[4]),
                                   float(filename_elements[5])),
                    'sensor': sensor,
                    'latitude': float(header_elements[7]),
                    'longitude': float(header_elements[8]),
                    'elevation': float(header_elements[9])}

        return metadata

    def _get_metadata_header_values(self):
        """
        Get metadata file in the format called Header Values.

        Returns
        -------
        metadata : dict
            Metadata information.
        """
        header_elements, filename_elements = self._get_metadata_from_file()

        if len(filename_elements) > 9:
            sensor = '_'.join(filename_elements[6:len(filename_elements) - 2])
        else:
            sensor = filename_elements[6]

        if filename_elements[3] in variable_lookup:
            variable = variable_lookup[filename_elements[3]]
        else:
            variable = filename_elements[3]

        metadata = {'network': header_elements[1],
                    'station': header_elements[2],
                    'latitude': float(header_elements[3]),
                    'longitude': float(header_elements[4]),
                    'elevation': float(header_elements[5]),
                    'depth': Depth(float(header_elements[6]),
                                   float(header_elements[7])),
                    'variable': variable,
                    'sensor': sensor}

        return metadata

    def _get_metadata_from_file(self, delim='_'):
        """
        Read first line of file and split filename.
        Information is used to collect metadata information for all
        ISMN formats.

        Parameters
        ----------
        delim : str, optional
            File basename delimiter.

        Returns
        -------
        header_elements : list[str]
            First line of file split into list
        file_basename_elements : list[str]
            File basename without path split by 'delim'
        """
        if self.archive.zip:
            with TemporaryDirectory(prefix='ismn', dir=self.temp_root) as tempdir:
                filename = self.archive.extract_file(self.file_path, tempdir)

                with filename.open(mode='r', newline=None) as f:
                    header = f.readline()
        else:
            filename = self.archive.path / self.file_path

            with filename.open(mode='r', newline=None) as f:
                header = f.readline()

        header_elements = header.split()
        path, basename = os.path.split(filename)
        file_basename_elements = basename.split(delim)

        return header_elements, file_basename_elements

    def _read_format_ceop_sep(self):
        """
        Read data in the file format called CEOP in separate files.
        """
        names = ['date', 'time', self.metadata['variable'],
                 self.metadata['variable'] + '_flag',
                 self.metadata['variable'] + '_orig_flag']
        usecols = [0, 1, 12, 13, 14]

        self.data = self._read_csv(names, usecols)

    def _read_format_header_values(self):
        """
        Read data file in the format called Header Values.
        """
        names = ['date', 'time', self.metadata['variable'],
                 self.metadata['variable'] + '_flag',
                 self.metadata['variable'] + '_orig_flag']

        self.data = self._read_csv(names, skiprows=1)

    def _read_csv(self, names=None, usecols=None, skiprows=0):
        """
        Read data.

        Parameters
        ----------
        names : list, optional
            List of column names to use.
        usecols : list, optional
            Return a subset of the columns.

        Returns
        -------
        data : pandas.DataFrame
            Time series.
        """
        readf = lambda f: pd.read_csv(f, skiprows=skiprows, usecols=usecols,
                                      names=names, delim_whitespace=True,
                                      parse_dates=[[0, 1]])
        if self.archive.zip:
            with TemporaryDirectory(prefix='ismn', dir=self.temp_root) as tempdir:
                filename = self.archive.extract_file(self.file_path, tempdir)
                data = readf(filename)
        else:
            data = readf(self.archive.path / self.file_path)

        data.set_index('date_time', inplace=True)

        return data


def usecase_file_zip():
    path = r"C:\Temp\delete_me\ismn\testdata_ceop.zip"
    archive = IsmnArchive(path)
    f = IsmnFile(archive, file_path='FMI\SAA111\FMI_FMI_SAA111_sm_0.050000_0.050000_5TE_20101001_20201005.stm',
                     load_data=False)
    d2 = f.load_data()

def usecase_coll():
    path = r"C:\Temp\delete_me\ismn\testdata_ceop"
    coll = IsmnFileCollection(path, load_data=False)
    nets = coll.get_networks()
    stats = coll.get_stations(None)
    sens = coll.get_sensors(nets[0], stats[0])


if __name__ == '__main__':
    usecase_file_zip()


    coll = IsmnFileCollection(r"C:\Temp\delete_me\ismn\testdata_ceop")
    coll.get_networks()
    coll.get_stations('FMI')
    coll.get_sensors('FMI', 'SOD021')
