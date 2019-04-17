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
# import glob
import logging

import pandas as pd

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(levelname)s - %(asctime)s: %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)


class DataCollection(object):

    def __init__(self):
        self.networks = {}

    def add_network(self, name):
        """
        Add network to data collection.

        Parameters
        ----------
        name : str
            Network name.
        """
        if name not in self.networks:
            self.networks[name] = Network(name)

    def add_station(self, network_name, station_name, lon, lat,
                    elev=None, static_variables=None, landcover=None,
                    climate_class=None):

        if network_name not in self.networks:
            self.add_network(network_name)

        if station_name not in self.networks[network_name]:
            self.networks[network_name].add_station(
                station_name, lon, lat, elev, static_variables,
                landcover, climate_class)

    def add_sensor(self, network_name, station_name):
        pass


class Network(object):

    def __init__(self, name):
        self.name = name
        self.stations = {}

    def add_station(self, name, lon, lat, elev=None, static_variables=None,
                    landcover=None, climate_class=None):

        if name not in self.stations:
            self.stations[name] = Station(
                name, lon, lat, elev, static_variables, landcover,
                climate_class)


class Station(object):

    def __init__(self, name, lon, lat, elev=None, static_variables=None,
                 landcover=None, climate_class=None):

        self.lon = lon
        self.lat = lat
        self.elev = elev
        self.static_variables = static_variables
        self.landcover = landcover
        self.climate_class = climate_class
        self.sensors = None


class Sensor(object):

    def __init__(self, variable, data, instrument,
                 depth_from=None, depth_to=None):

        self.variable = variable
        self.data = data
        self.instrument = instrument
        self.depth_from = depth_from
        self.depth_to = depth_to
        self.filename = None


def init_data(path):
    """
    Initialize ISMN data from file path.

    Parameters
    ----------
    path : str
        Path to downloaded ISMN data.
    """
    log_filename = os.path.join(path, 'python_metadata', 'metadata.log')
    fh = logging.FileHandler(log_filename)
    logger.addHandler(fh)

    ismn_files = []

    for root, sub_folders, basenames in os.walk(path):
        sub_folders.sort()
        basenames.sort()
        for basename in basenames:
            if basename.endswith('.stm'):
                filename = os.path.join(root, basename)
                logger.debug('Reading file {}'.format(filename))
                ismn_files.append(IsmnFile(filename))

    import pdb
    pdb.set_trace()
    pass


variable_lookup = {'sm': 'soil moisture',
                   'ts': 'soil temperature',
                   'su': 'soil suction',
                   'p': 'precipitation',
                   'ta': 'air temperature',
                   'fc': 'field capacity',
                   'wp': 'permanent wilting point',
                   'paw': 'plant available water',
                   'ppaw': 'potential plant available water',
                   'sat': 'saturation',
                   'si_h': 'silt fraction',
                   'sd': 'snow depth',
                   'sa_h': 'sand fraction',
                   'cl_h': 'clay fraction',
                   'oc_h': 'organic carbon',
                   'sweq': 'snow water equivalent',
                   'tsf': 'surface temperature',
                   'tsfq': 'surface temperature quality flag original'}


class IsmnFile(object):

    """
    IsmnFile class represents a single ISMN file.

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
    """

    def __init__(self, filename, load_data=False):
        self.filename = filename
        self.file_type = 'undefined'
        self.metadata = {}
        self.data = None
        self._read_metadata()

        if load_data:
            self.load_data()

    def load_data(self):
        """
        Load data from file.
        """
        if self.data is None:
            if self.file_type == 'ceop':
                self._read_format_ceop()
            elif self.file_type == 'ceop_sep':
                self._read_format_ceop_sep()
            elif self.file_type == 'header_values':
                self._read_format_header_values()
            else:
                logger.warning("Unknown file type: {}".format(self.filename))

    def _read_metadata(self):
        """
        Read metadata from file name and first line of file.
        """
        header_elements, filename_elements = self._get_metadata_from_file()

        if len(filename_elements) == 5 and len(header_elements) == 16:
            self.metadata = self._get_metadata_ceop()
            self.file_type = 'ceop'
        elif len(header_elements) == 15 and len(filename_elements) >= 9:
            self.metadata = self._get_metadata_ceop_sep()
            self.file_type = 'ceop_sep'
        elif len(header_elements) < 14 and len(filename_elements) >= 9:
            self.metadata = self._get_metadata_header_values()
            self.file_type = 'header_values'
        else:
            logger.warning("Unknown file type: {}".format(self.filename))

    def _get_metadata_ceop(self):
        """
        Get metadata in the file format called CEOP Reference Data Format.

        Returns
        -------
        metadata : dict
            Metadata information.
        """
        header_elements, filename_elements = self._get_metadata_from_file()

        metadata = {'network': filename_elements[1],
                    'station': header_elements[6],
                    'variable': ['ts', 'sm'],
                    'sensor': 'n.s',
                    'depth_from': ['multiple'],
                    'depth_to': ['multiple'],
                    'latitude': float(header_elements[7]),
                    'longitude': float(header_elements[8]),
                    'elevation': float(header_elements[9])}
        return metadata

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
            variable = [variable_lookup[filename_elements[3]]]
        else:
            variable = [filename_elements[3]]

        metadata = {'network': filename_elements[1],
                    'station': filename_elements[2],
                    'variable': variable,
                    'depth_from': [float(filename_elements[4])],
                    'depth_to': [float(filename_elements[5])],
                    'sensor': sensor,
                    'latitude': float(header_elements[7]),
                    'longitude': float(header_elements[8]),
                    'elevation': float(header_elements[9])}

        return metadata

    def _get_metadata_header_values(self):
        """
        Get metadata file in the format called ?.

        Returns
        -------
        metadata : dict
            Metadata information.
        """
        header_elements, filename_elements = self._get_info_from_file()

        if len(filename_elements) > 9:
            sensor = '_'.join(filename_elements[6:len(filename_elements) - 2])
        else:
            sensor = filename_elements[6]

        if filename_elements[3] in variable_lookup:
            variable = [variable_lookup[filename_elements[3]]]
        else:
            variable = [filename_elements[3]]

        metadata = {'network': header_elements[1],
                    'station': header_elements[2],
                    'latitude': float(header_elements[3]),
                    'longitude': float(header_elements[4]),
                    'elevation': float(header_elements[5]),
                    'depth_from': [float(header_elements[6])],
                    'depth_to': [float(header_elements[7])],
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
        header_elements : list
            First line of file split into list
        file_basename_elements : list
            File basename without path split by 'delim'
        """
        with io.open(self.filename, mode='r', newline=None) as f:
            header = f.readline()

        header_elements = header.split()
        path, basename = os.path.split(self.filename)
        file_basename_elements = basename.split(delim)

        return header_elements, file_basename_elements

    def _read_format_ceop(self):
        """
        """
        names = ['date', 'time', 'depth_from',
                 self.metadata['variable'][0],
                 self.metadata['variable'][0] + '_flag',
                 self.metadata['variable'][1],
                 self.metadata['variable'][1] + '_flag'],
        usecols = [0, 1, 11, 12, 13, 14, 15]

        # data needs to be converted, multi-index?
        self.data = self._read_csv(names, usecols)

    def _read_format_ceop_sep(self):
        """
        """
        names = ['date', 'time', self.metadata['variable'][0],
                 self.metadata['variable'][0] + '_flag',
                 self.metadata['variable'][0] + '_orig_flag']
        usecols = [0, 1, 12, 13, 14]

        self.data = self._read_csv(names, usecols)

    def _read_format_header_values(self):
        """
        """
        names = ['date', 'time', self.metadata['variable'][0],
                 self.metadata['variable'][0] + '_flag',
                 self.metadata['variable'][0] + '_orig_flag']

        self.data = self._read_csv(names)

    def _read_csv(self, names=None, usecols=None):
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
        data = pd.read_csv(self.filename, skiprows=1, usecols=usecols,
                           names=names, delim_whitespace=True,
                           parse_dates=[[0, 1]])

        return data.set_index('date_time', inplace=True)


def main():
    """
    Main routine.
    """
    path = os.path.join('/data2', 'shahn', 'datapool', 'ismn')
    ismn_data = init_data(path)


if __name__ == '__main__':
    main()
