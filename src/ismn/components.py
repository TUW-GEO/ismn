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

from pygeogrids import BasicGrid

import numpy as np
import warnings
import logging

import pandas as pd

logger = logging.getLogger(__name__)

ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter('%(levelname)s - %(asctime)s: %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)

class IsmnComponent: pass

class Network(IsmnComponent):
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
    coords : list, list
        Station lats and lons
    grid : BasicGrid
        Staton locations as grid

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
        # todo: using station dicts means that duplicate station names are not possible

    @property
    def coords(self) -> (list, list):
        # get lists of lats and lons for all stations in the network
        lons, lats = [], []
        for name, stat in self.stations.items():
            lons.append(stat.lon)
            lats.append(stat.lat)

        return lons, lats

    @property
    def grid(self):
        # get grid for all stations in network
        return BasicGrid(*self.coords)

    def add_station(self, name, lon, lat, elev, static_variables=None):
        """
        Add a station to the network.

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


class Station(IsmnComponent):

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
    static_variables : MetaData
        Station static variables

    Methods
    -------
    get_variables()
        Get variables measured by all sensors at station.
    get_depths(variable)
        Get depths of sensors measuring at station.
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

    def __repr__(self):
        """
        Provide basic station information.

        Returns
        -------
        info : str
            Basic station information.
        """
        info = 'Station {} has {} sensors.'.format(self.name,
                                                   len(self.sensors.keys()))
        return info

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
        depths : list
            List of depths of all sensors that measure the passed variable.
        """
        depths = []
        for sensor in self.sensors.values():
            if sensor.eval(variable=variable):
                depths.append(sensor.depth)

        return depths

    def add_sensor(self, instrument, variable, depth, filehandler, name=None,
                   keep_loaded_data=False):
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
        name : str or int, optional (default: None)
            A name or id for the sensor. If None is passed, one is generated.
        """
        if name is None:
            name = f"{instrument}_{variable}_{depth.start:1.6f}_{depth.end:1.6f}"

        if name not in self.sensors:
            self.sensors[name] = Sensor(instrument, variable, depth, name,
                                        filehandler, keep_loaded_data)
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

    def iter_sensors(self, variable=None, depth=None, **kwargs):
        """
        Get all sensors for variable and/or depth information.

        Parameters
        ----------
        variable : str, optional
            Observed variabl e.
        depth: Depth, optinal (default: None)
            Yield only sensors that are enlosed by this depth.
        kwargs:
            All other kwargs as used in the sensor.eval() function to
            evaluate the sensor.

        Yields
        ------
        sensors : Sensor
            Sensor that measure the passsed var within the passed depth.
        """
        for sensor in self.sensors.values():
            if sensor.eval(variable, depth, **kwargs):
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


class Sensor(IsmnComponent):

    """
    A Sensor with ground observations.

    Parameters
    ----------
    instrument : str
        Instrument name.
    variable : str
        Observed variable.
    depth : Depth
        Sensing depth.
    name : str or int, optional (default: None)
        Id or Name of the sensor. If None is passed, a name is generated.
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

    def __init__(self, instrument, variable, depth, name=None, filehandler=None,
                 keep_loaded_data=False):

        self.instrument = instrument
        self.variable = variable
        self.depth = depth
        self.filehandler = filehandler
        self.keep_loaded_data = keep_loaded_data
        self.data = None
        self.name = name if name is not None else self.__repr__()

    def __repr__(self):
        return f"{self.instrument}_{self.variable}_{self.depth.start:1.6f}_{self.depth.end:1.6f}"

    def read_data(self) -> pd.DataFrame:
        """
        Load data from filehandler for sensor.
        """
        if self.filehandler is None:
            warnings.warn(f"No filehandler found for sensor {self.name}")
        else:
            if self.data is None:
                data = self.filehandler.read_data()

                if self.keep_loaded_data:
                    self.data = data

                return data
            else:
                return self.data

    def eval(self, variable=None, depth=None, filter_meta_dict=None,
             check_only_sensor_depth_from=False):
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
        check_only_sensor_depth_from : bool, optional (default: False)
            Ignores the sensors depth_to value and only checks if depth_from of
            the sensor is in the passed depth (e.g. for cosmic ray probes).

        Returns
        -------
        flag : bool
            Indicates success.
        """

        if depth is None:
            depth = Depth(-np.inf, np.inf)

        flag = False
        
        if check_only_sensor_depth_from:
            d = Depth(self.depth.start, self.depth.start)
        else:
            d = self.depth 

        if (variable in [None, self.variable]) and depth.encloses(d):
            flag = True

        if flag and filter_meta_dict:
            if self.filehandler is None:
                warnings.warn("Filehandler is None, can't filter by metadata")
            else:
                # checks also if the metadata in file matches
                flag = self.filehandler.check_metadata(
                    variable, depth.start, depth.end, filter_static_vars=filter_meta_dict)

        return flag

class Depth():

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

        # todo: could add unit to depth?

        self.start = float(start)
        self.end = float(end)

        self.extent = self.end - self.start

        # todo: allow negative depths? i.e above surface, for Temperature?
        # if self.extent < 0:
        #     raise ValueError("End can not be smaller than start")

        if self.start == self.end:
            self.is_profile = False
        else:
            self.is_profile = True

    def __str__(self):
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

    def __iter__(self):
        for d in [self.start, self.end]:
            yield d

    def perc_overlap(self, other):
        """
        Estimate how much 2 depths correspond.
        1 means that the are the same, 0 means that they have an infinitely
        small correspondence (e.g. a single layer within a range, or 2 adjacent
        depths). -1 means that they dont overlap.

        Parameters
        ----------
        other : Depth
            Second depth

        Returns
        -------
        p : float
            Normalised overlap range
            <0 = no overlap, 0 = adjacent, >0 = overlap, 1 = equal
        """
        if self == other: # same depths
            return 1
        else:
            r = max([self.end, other.end]) - min([self.start, other.start])
            # Overlapping range normalised to the overall depth range r
            p_f = abs(self.start - other.start) / r
            p_t = abs(self.end - other.end) / r
            p = 1 - p_f - p_t

            if p < 0:
                p = -1

        return p

    def overlap(self, other, return_perc=False):
        """
        Check if two depths overlap, (if the start of one depth is the same as
        the end of the other, they overlap,
        e.g. Depth(0, 0.1) and Depth(0.1, 0.2) do overlap.

        Parameters
        ----------
        other : Depth
            Other Depth
        return_perc : bool, optional (default: False)
            Returns how much the depths overlap. See func: perc_overlap()

        Returns
        -------
        overlap : bool
            True if Depths overlap
        perc_overlap: float
            Normalised overlap.
        """
        other_start_encl = self.encloses(Depth(other.start, other.start))
        other_end_encl = self.encloses(Depth(other.end, other.end))

        this_start_encl = other.encloses(Depth(self.start, self.start))
        this_end_encl = other.encloses(Depth(self.end, self.end))

        overlap = any([other_start_encl, other_end_encl,
                       this_start_encl, this_end_encl])

        if return_perc:
            return overlap, self.perc_overlap(other)
        else:
            return overlap

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

if __name__ == '__main__':
    d1 = Depth(0, 0.3)
    d2 = Depth(-0.01, -0.1)

    overlaps, perc = d1.overlap(d2, return_perc=True)
