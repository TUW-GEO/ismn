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
from ismn.readers import get_format
import pandas as pd



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
            csv_meta = get_metadata_from_csv(filepath, as_dict=True)
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

        metadata = get_metadata(filepath)
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

def get_metadata_from_csv(filename, as_dict=False):
    """
    reads ISMN metadata from csv file

    Parameters
    ----------
    filename: str, path to csv file

    Returns
    -------
    landcover_2000: int, cci landcover classification for station (year 2000)
    landcover_2005: int, cci landcover classification for station (year 2005)
    landcover_2010: int, cci landcover classification for station (year 2010)
    landcover_insitu: str, in situ landcover classification
    climate: str, Koeppen Geiger climate classification for station
    climate_insitu: str, in situ climate classification for station
    saturation: nd.array, saturation for all available depths
    clay_fraction: nd.array, clay fraction for all available depths (in % weight)
    sand_fraction: nd.array, sand fraction for all available depths (in % weight)
    silt_fraction: nd.array, silt fraction for all available depths (in % weight)
    organic_carbon: nd.array, organic carbon for all available depths (in % weight)
    """
    def read_field(fieldname):
        if fieldname in data.index:
            dt = list()
            for i, j in zip(np.atleast_1d(data.loc[fieldname]['depth_from[m]']),
                            np.atleast_1d(data.loc[fieldname]['depth_to[m]'])):
                dt.append(('{}m_{}m'.format(i, j), np.float))
            return np.array([tuple(np.atleast_1d(data.loc[fieldname]['value']))], dtype=np.dtype(dt))
        else:
            return np.nan

    # some stations don't come with correct format in csv file (missing header)
    try:
        data = pd.read_csv(filename, delimiter=";")
        data.set_index('quantity_name', inplace=True)
    except:
        # set columns manually
        logging.info('no header: {}'.format(filename))
        data = pd.read_csv(filename, delimiter=";", header=None)
        cols = list(data.columns.values)
        cols[0:7] = ['quantity_name', 'unit', 'depth_from[m]', 'depth_to[m]',
                     'value', 'description', 'quantity_source_name']
        data.columns = cols
        data.set_index('quantity_name', inplace=True)

    # read landcover classifications
    lc = data.loc[['land cover classification']][['value', 'quantity_source_name']]
    lc_dict = {'CCI_landcover_2000': np.nan, 'CCI_landcover_2005': np.nan,
               'CCI_landcover_2010': np.nan, 'insitu': ''}
    for key in lc_dict.keys():
        if key in lc['quantity_source_name'].values:
            if key != 'insitu':
                lc_dict[key] = np.int(lc.loc[lc['quantity_source_name'] == key]['value'].values[0])
            else:
                lc_dict[key] = lc.loc[lc['quantity_source_name'] == key]['value'].values[0]
                logging.info('insitu land cover classification available: {}'.format(filename))

    # read climate classifications
    cl = data.loc[['climate classification']][['value', 'quantity_source_name']]
    cl_dict = {'koeppen_geiger_2007': '', 'insitu': ''}
    for key in cl_dict.keys():
        if key in cl['quantity_source_name'].values:
            cl_dict[key] = cl.loc[cl['quantity_source_name'] == key]['value'].values[0]
            if key == 'insitu':
                logging.info('insitu climate classification available: {}'.format(filename))

    saturation = read_field('saturation')
    clay_fraction = read_field('clay fraction')
    sand_fraction = read_field('sand fraction')
    silt_fraction = read_field('silt fraction')
    organic_carbon = read_field('organic carbon')

    if as_dict:
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
    else:
        return lc_dict['CCI_landcover_2000'], \
               lc_dict['CCI_landcover_2005'], \
               lc_dict['CCI_landcover_2010'], \
               lc_dict['insitu'], \
               cl_dict['koeppen_geiger_2007'], \
               cl_dict['insitu'], \
               saturation, \
               clay_fraction, \
               sand_fraction, \
               silt_fraction, \
               organic_carbon


def get_metadata(filename):
    """
    reads ISMN metadata from any format

    Parameters
    ----------
    filename: string

    Returns
    -------
    metadata: dict
    """
    dicton = globals()
    func = dicton['get_metadata_' + get_format(filename)]
    return func(filename)

if __name__ == '__main__':
    data = r"C:\Temp\delete_me\ismn\testdata_ceop"
    meta_path = os.path.join(data, 'python_metadata')
    col = MetaCollector(data, meta_path)
    meta = col.collect_from_archive()
    #col.get_station_meta(os.path.join('FMI', 'SOD140'))
