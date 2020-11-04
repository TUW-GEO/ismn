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

from pygeogrids.grids import BasicGrid

logger = logging.getLogger(__name__)

ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter('%(levelname)s - %(asctime)s: %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)

variable_lookup = {'sm': 'soil_moisture',
                   'ts': 'soil_temperature',
                   'su': 'soil_suction',
                   'p': 'precipitation',
                   'ta': 'air_temperature',
                   'fc': 'field_capacity',
                   'wp': 'permanent_wilting_point',
                   'paw': 'plant_available_water',
                   'ppaw': 'potential_plant_available_water',
                   'sat': 'saturation',
                   'si_h': 'silt_fraction',
                   'sd': 'snow_depth',
                   'sa_h': 'sand_fraction',
                   'cl_h': 'clay_fraction',
                   'oc_h': 'organic_carbon',
                   'sweq': 'snow_water_equivalent',
                   'tsf': 'surface_temperature',
                   'tsfq': 'surface_temperature_quality_flag_original'}


class NetworkCollection(object):

    """
    A hierarchical description of Networks, Stations and Sensors
    derived from an IsmnFileCollection.

    Parameters
    ----------
    file_collection : IsmnFileCollection
        Collection of ISMN files described as IsmnFileCollection.

    Attributes
    ----------
    file_collection : IsmnFileCollection
        File collection.
    networks : dict of Network
        Networks.
    grid : pygeogrids.BasicGrid
        Latitude/longitude coordinates of stations.
    grid_lut : numpy.ndarray
        Look-up table between grid point indices and file index.

    Methods
    -------
    add_network(name)
        Add network to collection.
    iter_sensors(network=None, station=None, variable=None, depth=None)
        Get all sensors for a specific network and/or station and/or
        variable and/or depth.
    get_nearest_station(lon, lat, max_dist=np.inf)
        Get nearest station for given longitude/latitude coordinates.
    """

    def __init__(self, file_collection):

        self.file_collection = file_collection
        self.networks = {}

        lon = []
        lat = []
        idx = []

        for f_id, f in self.file_collection.files.items():
            nw_name = f.metadata['network']
            st_name = f.metadata['station']
            se_name = f.metadata['sensor']

            if nw_name not in self.networks:
                self.add_network(nw_name)

            if st_name not in self.networks[nw_name].stations:
                self.networks[nw_name].add_station(
                    st_name, f.metadata['longitude'], f.metadata['latitude'],
                    f.metadata['elevation'])
                idx.append(f_id)
                lon.append(f.metadata['longitude'])
                lat.append(f.metadata['latitude'])

            if se_name not in self.networks[nw_name].stations[st_name].sensors:
                self.networks[nw_name].stations[st_name].add_sensor(
                    se_name, f.metadata['variable'], f.metadata['depth'], f)

        self.grid = BasicGrid(lon, lat)
        self.grid_lut = np.array(idx)

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
        else:
            logger.warning('Network already exists: {}'.format(name))

    def get_sensors(self, network=None, station=None, variable=None,
                    depth=None):
        """
        Yield all sensors for a specific network and/or station and/or
        variable and/or depth.

        Parameters
        ----------
        network : str, optional
            Network name (default: None).
        station : str, optional
            Station name (default: None).
        variable : str, optional
            Variable name (default: None).
        depth : Depth, optional
            Sensing depth (default: None).

        Yield
        -----
        sensor : Sensor
            Sensor.
        """
        if depth is not None:
            d = Depth(depth[0], depth[1])
        else:
            d = Depth(-np.inf, np.inf)

        for n in self.networks.values():
            if network not in [None, n.name]:
                continue
            for st in n.stations.values():
                if station not in [None, st.name]:
                    continue
                for se in st.sensors.values():
                    if (variable not in [None, se.variable] and
                            se.depth.enclose(d)):
                        continue

                    yield se

    def get_nearest_station(self, lon, lat, max_dist=np.inf):
        """
        Get nearest station for given longitude/latitude coordinates.

        Parameters
        ----------
        lon : float
            Longitude coordinate.
        lat : float
            Latitude coordinate.
        max_dist : float, optional
            Maximum search distance (default: numpy.inf).

        Returns
        -------
        station : Station
            Station.
        dist : float
            Distance in meter.
        """
        gpi, dist = self.grid.find_nearest_gpi(lon, lat, max_dist)

        if dist != np.inf:
            idx = self.grid_lut[gpi]
            f = self.file_collection.files[idx]
            station = self.networks[f.metadata['network']].stations[
                f.metadata['station']]
        else:
            station = None

        return station, dist


class Network(object):
    """
    A network is described by a distinct name and can be composed of
    multiple stations.

    Parameters
    ----------
    name : str
        Network name.

    Attributes
    ----------
    name : str
        Network name.
    stations : dict of Station
        Stations belonging to the network.

    Methods
    -------
    add_station(name, lon, lat, elev, static_variables=None)
        Add station to network.
    remove_station(name)
        Remove station from network.
    iter_stations(variable=None, depth=None)
        Get all stations having at least one sensor observing
        a specific variable and/or sensing depth.
    n_stations()
        Number of stations.
    """

    def __init__(self, name):
        self.name = name
        self.stations = {}

    def add_station(self, name, lon, lat, elev, static_variables=None):
        """
        Add station to network.

        Parameters
        ----------
        name : str
            Station name.
        lon : float
            Longitude coordinate.
        lat : float
            Latitude coordinate.
        elev : float
            Elevation.
        static_variables : list, optional
            Static variable information (default: None).
        """
        if name not in self.stations:
            self.stations[name] = Station(name, lon, lat, elev,
                                          static_variables)
        else:
            logger.warning('Station already exists: {}'.format(name))

    def remove_station(self, name):
        """
        Remove station from network.

        Parameters
        ----------
        name : str
            Station name.
        """
        if name in self.stations:
            del self.stations[name]
        else:
            logger.warning('Station not found {}'.format(name))

    def get_stations(self, variable=None, depth=None):
        """
        Get all stations having at least one sensor observing
        a specific variable and/or sensing depth.

        Parameters
        ----------
        variable : str, optional
            Observed variable.
        depth : list, optional
            Sensing depth.

        Yields
        ------
        station : Station
            Station.
        """
        if depth is not None:
            d = Depth(depth[0], depth[1])
        else:
            d = Depth(-np.inf, np.inf)

        for station in self.stations.values():
            flag = False
            for sensor in station.sensors.values():
                if (variable in [None, sensor.variable] and
                        sensor.depth.enclose(d)):
                    flag = True
            if flag:
                yield station

    def n_stations(self):
        """
        Number of stations.

        Returns
        -------
        n : int
            Number of stations in the network.
        """
        return len(self.stations)

    def __repr__(self):
        """
        Provide basic network information.

        Returns
        -------
        info : str
            Basic network information.
        """
        info = 'Network {} has {} stations.'.format(self.name,
                                                    self.n_stations())
        return info


class Station(object):

    """
    A station is described by a distinct name and location.
    Multiple sensors at various depths can be part of a station.

    Parameters
    ----------
    name : str
        Station name.
    lon : float
        Longitude coordinate.
    lat : float
        Latitude coordinate.
    elev : float
        Elevation information.
    static_variables : list, optional
        Static variable information (default: None).

    Attributes
    ----------
    name : str
        Station name.
    lon : float
        Longitude coordinate.
    lat : float
        Latitude coordinate.
    elev : float
        Elevation information.
    static_variables : list, optional
        Static variable information (default: None).

    Methods
    -------
    add_sensor(instrument, variable, depth, filehandler)
        Add sensor to station.
    remove_sensor(name)
        Remove sensor from station.
    iter_sensors(variable=None, depth=None)
        Get all sensors for variable and/or depth information.
    n_sensors()
        Number of sensors.
    """

    def __init__(self, name, lon, lat, elev, static_variables=None):
        self.name = name
        self.lon = lon
        self.lat = lat
        self.elev = elev
        self.static_variables = static_variables
        self.sensors = {}

    def add_sensor(self, instrument, variable, depth, filehandler):
        """
        Add sensor to station.

        Parameters
        ----------
        instrument : str
            Instrument name.
        variable : str
            Observed variable.
        depth : Depth
            Sensing depth.
        filehandler : IsmnFile
            File handler.
        """
        name = '{}_{}_{:1.6f}_{:1.6f}'.format(
            instrument, variable, depth.start, depth.end)

        if name not in self.sensors:
            self.sensors[name] = Sensor(name, instrument, variable,
                                        depth, filehandler)
        else:
            logger.warning('Sensor already exists: {}'.format(name))

    def remove_sensor(self, name):
        """
        Remove sensor from station.

        Parameters
        ----------
        name : str
            Sensor name.
        """
        if name in self.sensors:
            del self.sensors[name]
        else:
            logger.warning('Sensor not found: {}'.format(name))

    def get_sensors(self, variable=None, depth=None):
        """
        Get all sensors for variable and/or depth information.

        Parameters
        ----------
        variable : str, optional
            Observed variable.
        depth : list, optional
            Sensing depth.

        Yields
        ------
        sensors : Sensor
            Sensor.
        """
        if depth is not None:
            d = Depth(depth[0], depth[1])
        else:
            d = Depth(-np.inf, np.inf)

        for sensor in self.sensors.values():
            if variable in [None, sensor.variable] and sensor.depth.enclose(d):
                yield sensor

    def n_sensors(self):
        """
        Number of sensors.

        Returns
        -------
        n : int
            Number of sensors at the station.
        """
        return len(self.sensors)


class Sensor(object):

    """
    A Sensor with ground observations.

    Parameters
    ----------
    name : str
        Name of the sensor.
    instrument : str
        Instrument name.
    variable : str
        Observed variable.
    depth : Depth
        Sensing depth.
    filehandler : IsmnFile, optional
        File handler (default: None).

    Attributes
    ----------
    name : str
        Name of the sensor.
    instrument : str
        Instrument name.
    variable : str
        Observed variable.
    depth : Depth
        Sensing depth.
    filehandler : IsmnFile, optional
        File handler (default: None).
    """

    def __init__(self, name, instrument, variable, depth, filehandler=None):

        self.name = name
        self.instrument = instrument
        self.variable = variable
        self.depth = depth
        self.filehandler = filehandler


class Depth(object):

    """
    A class representing a depth (0=surface).

    Parameters
    ----------
    start : float
        Depth start.
    end : float
        Depth end.

    Attributes
    ----------
    start : float
        Depth start.
    end : float
        Depth end.

    Methods
    -------
    __eq__(other)
        Test if two Depth are equal.
    enclose(other)
        Test if other Depth encloses given Depth.
    """

    def __init__(self, start, end):
        self.start = start
        self.end = end
        self.extent = self.end - self.start

        if self.start == self.end:
            self.is_profile = False
        else:
            self.is_profile = True

    def __eq__(self, other):
        """
        Test if two Depth are equal.

        Parameters
        ----------
        other : Depth
            Depth.

        Returns
        -------
        flag : bool
            True if both depths are equal, False otherwise.
        """
        if self.start == other.start and self.end == other.end:
            flag = True
        else:
            flag = False

        return flag

    def enclose(self, other):
        """
        Test if other Depth encloses given Depth.

        Parameters
        ----------
        other : Depth
            Depth.

        Returns
        -------
        flag : bool
            True if other depth surrounds given depth, False otherwise.
        """
        if (other.start <= self.start) and (other.end >= self.end):
            flag = True
        else:
            flag = False

        return flag


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
    path : str
        Root path of ISMN files.
    files : list
        List of ISMN filenames.

    Methods
    -------
    get_networks()
        Get networks from ISMN file collection.
    iter_stations(network=None)
        Get stations from ISMN file collection.
    iter_sensors(self, network=None, station=None)
        Get sensors from ISMN file collection.
    """

    def __init__(self, path, load_data=False):

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
                    self.files[i] = IsmnFile(filename, load_data)
                    i = i + 1

    def get_networks(self):
        """
        Get networks from ISMN file collection.

        Returns
        -------
        networks : Dict of Networks
            Dict of networks.
        """
        networks = {}

        for f in self.files.values():
            if f.metadata['network'] not in networks:
                networks[f.metadata['network']] = Network(
                    f.metadata['network'])

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

        for f in self.files.values():
            if network in [None, f.metadata['network']]:
                if f.metadata['station'] not in stations:
                    st = Station(f.metadata['station'], f.metadata['longitude'],
                                 f.metadata['latitude'], f.metadata['elevation'])
                    stations[f.metadata['station']] = st

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

        for f in self.files.values():
            if ((network in [None, f.metadata['network']]) and
                    (station in [None, f.metadata['station']])):

                name = '{}_{}_{}_{}'.format(f.metadata['sensor'],
                                            f.metadata['variable'],
                                            f.metadata['depth'].start,
                                            f.metadata['depth'].end)

                if name not in sensors:
                    snr = Sensor(name, f.metadata['variable'],
                                 f.metadata['sensor'], f.metadata['depth'])
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
    IsmnFile class represents a single ISMN file.

    Parameters
    ----------
    filename : str
        Filename.
    load_data : bool, optional
        If True data will be loaded during metadata reading.

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
    read_metadata()
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

    def __init__(self, filename, load_data=False):

        if not os.path.isfile(filename):
            raise IOError('File does not exist: {}'.format(filename))

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
                # self._read_format_ceop()
                raise NotImplementedError
            elif self.file_type == 'ceop_sep':
                self._read_format_ceop_sep()
            elif self.file_type == 'header_values':
                self._read_format_header_values()
            else:
                logger.warning("Unknown file type: {}".format(self.filename))

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
            logger.warning("Unknown file type: {}".format(self.filename))

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
        data = pd.read_csv(self.filename, skiprows=skiprows, usecols=usecols,
                           names=names, delim_whitespace=True,
                           parse_dates=[[0, 1]])

        data.set_index('date_time', inplace=True)

        return data


def create_network_collection(path, load_data=False):
    """
    Create a network collection for given file path.

    Parameters
    ----------
    path : str
        Path to data.
    load_data : bool, optional
        Load data while reading metadata (default: False)
        Attention: Can be slow for a large number of files.

    Returns
    -------
    nwc : NetworkCollection
        Metadata and data (if loaded) from in situ data.
    """
    fc = IsmnFileCollection(path, load_data)
    nwc = NetworkCollection(fc)

    return nwc

if __name__ == '__main__':
    coll = IsmnFileCollection(r"C:\Temp\delete_me\ismn\testdata_ceop")
    coll.get_networks()
    coll.get_stations('FMI')
    coll.get_sensors('FMI', 'SOD021')
