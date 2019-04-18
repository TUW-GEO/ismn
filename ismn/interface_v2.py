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

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter('%(levelname)s - %(asctime)s: %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)

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


class NetworkCollection(object):

    """
    A hierarchically description of an IsmnFileCollection.

    Attributes
    ----------
    file_collection : IsmnFileCollection
        File collection.
    networks : dict of Network
        Networks.

    """

    def __init__(self, file_collection):

        self.file_collection = file_collection
        self.networks = {}

        for f_id, f in self.file_collection.files.items():
            nw_name = f.metadata['network']
            st_name = f.metadata['station']
            se_name = f.metadata['sensor']

            if nw_name not in self.networks:
                self.add_network(nw_name)

            if st_name not in self.networks:
                self.networks[nw_name].add_station(
                    st_name, f.metadata['longitude'], f.metadata['latitude'],
                    f.metadata['elevation'])

            if se_name not in self.networks[nw_name].stations:
                self.networks[nw_name].stations[st_name].add_sensor(
                    se_name, f.metadata['variable'],
                    f.metadata['depth_from'], f.metadata['depth_to'])

    def add_network(self, name):
        """
        Add network to collection.

        Parameters
        ----------
        name : str
            Network name.
        """
        if name not in self.networks:
            self.networks[name] = Network(name)

    def get_sensors(self, network=None, station=None, variable=None,
                    depth_from=None, depth_to=None):
        """
        Get all sensors for specific variable and/or depth.

        Parameters
        ----------
        network : str, optional
            Network name (default: None).
        station : str, optional
            Station name (default: None).
        variable : str, optional
            Variable name (default: None).
        depth_from : float, optional
            Start sensing depth (default: None).
        depth_to : float, optional
            End sensing depth (default: None).

        Returns
        -------
        sensors : list of Sensor
            List of found sensors.
        """
        sensors = []

        for n in self.networks.values():
            if network not in [None, n.name]:
                continue
            for st in n.stations.values():
                if station not in [None, st.name]:
                    continue
                for se in st.sensors.values():
                    if variable not in [None, se.variable]:
                        continue
                    sensors.append(se)

        return sensors

    def get_nearest_station(self, lon, lat, max_dist=np.inf):
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

    def remove_station(self, name):
        del self.stations[name]


class Station(object):

    def __init__(self, name, lon, lat, elev=None, static_variables=None,
                 landcover=None, climate_class=None):

        self.name = name
        self.lon = lon
        self.lat = lat
        self.elev = elev
        self.static_variables = static_variables
        self.landcover = landcover
        self.climate_class = climate_class
        self.sensors = {}

    def add_sensor(self, instrument, variable, depth_from, depth_to):

        sensor_name = '{}_{}_{}_{}'.format(instrument, variable,
                                           depth_from, depth_to)

        if sensor_name not in self.sensors:
            self.sensors[sensor_name] = Sensor(instrument, variable,
                                               depth_from, depth_to)

    def remove_sensor(self, name):
        del self.sensors[name]


class Sensor(object):

    def __init__(self, instrument, variable, depth_from, depth_to):

        self.name = '{}_{}_{}_{}'.format(instrument, variable,
                                         depth_from, depth_to)
        self.instrument = instrument
        self.variable = variable
        self.depth_from = depth_from
        self.depth_to = depth_to


class IsmnFileCollection(object):

    """
    The IsmnFileCollection class reads and organized the metadata
    information of ISMN files.

    Attributes
    ----------
    path : str
        Root path of ISMN files.
    files : list
        List of ISMN filenames.

    Methods
    -------
    get_networks()

    get_stations()

    get_sensors()

    summary()

    """

    def __init__(self, path):

        self.path = path
        self.files = {}

        i = 0
        for root, sub_folders, basenames in os.walk(self.path):
            sub_folders.sort()
            basenames.sort()
            for basename in basenames:
                if basename.endswith('.stm'):
                    filename = os.path.join(root, basename)
                    logger.debug('Reading file {}'.format(filename))
                    self.files[i] = IsmnFile(filename)
                    i = i + 1

    def get_networks(self):
        """
        Get networks from ISMN file collection.

        Returns
        -------
        networks : List of Networks
            List of networks.
        """
        networks = []

        for f in self.files.values():
            networks.append(Network(f.metadata['network']))

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
        stations : list
            List of stations.
        """
        stations = []

        for f in self.files.values():
            if network in [None, f.metadata['network']]:
                st = Station(f.metadata['station'], f.metadata['longitude'],
                             f.metadata['latitude'], f.metadata['elevation'])
                stations.append(st)

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
        sensors : list
            List of sensors.
        """
        sensors = []

        for f in self.files.values():
            if ((network in [None, f.metadata['network']]) &
                    (station in [None, f.metadata['station']])):
                snr = Sensor(f.metadata['variable'], f.metadata['sensor'],
                             f.metadata['depth_from'], f.metadata['depth_to'])
                sensors.append(snr)

        return sensors

    def summary(self):
        """
        Print summary of ISMN file collection.
        """
        logger.info('Number of networks: {}'.format(len(self.get_networks())))
        logger.info('Number of stations: {}'.format(len(self.get_stations())))
        logger.info('Number of sensors: {}'.format(len(self.get_sensors())))
        logger.info('Number of files: {}'.format(len(self.files)))


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

    Methods
    -------
    load_data()

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
            variable = variable_lookup[filename_elements[3]]
        else:
            variable = filename_elements[3]

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

    fc = IsmnFileCollection(path)
    # print(fc.get_networks())
    # print(fc.get_stations())
    # print(fc.get_sensors())
    fc.summary()
    # st = fc.get_stations('SMOSMANIA')
    sen = fc.get_sensors(station='Narbonne')
    print(len(sen))

    nwc = NetworkCollection(fc)
    # print(nwc.networks['SMOSMANIA'].stations)
    # sen = nwc.get_sensors()
    # print(len(sen))
    # sen = nwc.get_sensors(station='Narbonne')
    # for s in sen:
    #     print(s.variable)
    print(len(nwc.get_sensors(station='Narbonne', variable='soil moisture')))
    print(len(nwc.get_sensors(station='Narbonne', variable='soil temperature')))


def main():
    """
    Main routine.
    """
    path = os.path.join('/data2', 'shahn', 'datapool', 'ismn')
    init_data(path)


if __name__ == '__main__':
    main()
