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

from collections import Sequence

from ismn.tables import *
from pygeogrids import BasicGrid

import warnings
import logging

logger = logging.getLogger(__name__)

ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter('%(levelname)s - %(asctime)s: %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)




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
        self.stations = {} # todo. using station dicts means that duplicate station names are not possible

    @property
    def coords(self) -> (list, list):
        lons, lats = [], []
        for name, stat in self.stations.items():
            lons.append(stat.lon)
            lats.append(stat.lat)

        return lons, lats

    @property
    def grid(self):
        return BasicGrid(*self.coords)

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
                                          static_variables=static_variables)
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

    def iter_stations(self, variable=None, depth=None, filter_dict=None):
        """
        Get all stations having at least one sensor observing
        a specific variable and/or sensing depth.

        Parameters
        ----------
        variable : str, optional
            Observed variable.
        depth : Depth, optional
            Sensing depth.

        Yields
        ------
        station : Station
            Station.
        """
        if depth is None:
            depth = Depth(-np.inf, np.inf)

        for station in self.stations.values():
            flag = False
            for sensor in station.sensors.values():
                if (variable in [None, sensor.variable]) and depth.encloses(sensor.depth):
                    flag = True
                if flag and filter_dict:
                    f = sensor.filehandler
                    if f is None:
                        warnings.warn("Filehandler is None, can't filter by metadata")
                    flag = sensor.filehandler.check_metadata(
                        variable, depth.start, depth.end, filter_static_vars=filter_dict)
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

    def get_variables(self):
        """
        Get variables measured by all sensors at station.

        Returns
        -------
        variables : list
            List of variables that are observed.
        """
        variables = []
        for sensor in self.sensors.values():
            if sensor.variable not in variables:
                variables.append(sensor.variable)

        return variables

    def get_depths(self, variable=None):
        """
        Get depths of sensors measuring at station.

        Parameters
        ----------
        variable : str, optional (default: None)
            Only consider sensors measuring this variable.

        Returns
        -------
        variables : list
            List of variables that are observed.
        """
        depths = []
        for sensor in self.sensors.values():
            if sensor.eval(variable=variable):
                depths.append(sensor.depth)

        return depths

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

    def iter_sensors(self, variable=None, depth=None):
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

    def read_data(self):
        """
        Load data from filehandler for sensor.
        """
        if self.filehandler is None:
            warnings.warn(f"No filehandler found for sensor {self.name}")
        else:
            return self.filehandler.read_data()

    def eval(self, variable=None, depth=None, filter_meta_dict=None):
        """
        Evaluate whether the sensor complies with the passed metadata requirements.

        Parameters
        ----------
        variable : str, optional (default: None)
            Check if the variable name matches
        depth : Depth, optional (default: None)
            Check if the passed depth encloses the sensor depth.
        filter_meta_dict : dict, optional (default: None)
            Additional metadata keys and values for which the file list is filtered
            e.g. {'lc_2010': 10} to filter for a landcover class.

        Returns
        -------
        flag : bool
            Indicates success.
        """

        if depth is None:
            depth = Depth(-np.inf, np.inf)

        flag = False

        if (variable in [None, self.variable]) and depth.encloses(self.depth):
            flag = True

        if flag and filter_meta_dict:
            if self.filehandler is None:
                warnings.warn("Filehandler is None, can't filter by metadata")
            else:
                # checks also if the metadata in file matches
                flag = self.filehandler.check_metadata(
                    variable, depth.start, depth.end, filter_static_vars=filter_meta_dict)

        return flag

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

        if self.extent < 0:
            raise ValueError("End can not be smaller than start")

        if self.start == self.end:
            self.is_profile = False
        else:
            self.is_profile = True

    def __str__(self):
        # todo: could add a unit to depth?
        return f"{self.start}_{self.end}[m]"

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
        if (self.start == other.start) and (self.end == other.end):
            flag = True
        else:
            flag = False

        return flag

    def encloses(self, other):
        """
        Test if this Depth encloses other Depth.

        Parameters
        ----------
        other : Depth
            Depth.

        Returns
        -------
        flag : bool
            True if other depth surrounds given depth, False otherwise.
        """
        if (self.start <= other.start) and (self.end >= other.end):
            flag = True
        else:
            flag = False

        return flag

    def enclosed(self, other):
        """
        Test if other Depth encloses this Depth.

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

class StaticVariable():
    """
    Represents a static variable for a ismn time series
    """
    def __init__(self, name, values:Sequence, depths:Sequence=None):

        self.name = name
        self.values = np.atleast_1d(values)

        if depths is not None:
            depths = np.atleast_1d(depths)
            assert len(depths) == len(values), "One value per depth expected"

        self.depths = depths

    def as_struct(self):
        #return variable in old format (struct array)
        vals, dtypes = [], []

        for v, d in zip(self.values, self.depths):
            dtypes.append((d, type(v)))
            vals.append(v)

        return np.array(vals, dtype=dtypes)

    def value_for_depth(self, depth_from, depth_to=None):
        """
        Get the respective value for the passed depth range.
        If no depth_to is passed, only depth from is used.
        One best matching value is found via: TODO:
        """
        pass


if __name__ == '__main__':
    net = Network('test')
    for name, num in [('one', 1), ('two', 2), ('three', 3)]:
        net.add_station(name, num, num, num)
        net.stations[name].add_sensor('sensor_name', 'instrument', Depth(0,1), None)
    net.stations['one'].sensors['sensor_name_instrument_0.000000_1.000000'].eval()
    col = NetworkCollection([net])