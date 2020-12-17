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

from ismn.meta import MetaData, Depth
from collections import OrderedDict

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

    def __init__(self,
                 name):
        """
        Parameters
        ----------
        name : str
            Network name.
        """
        self.name = name
        self.stations = OrderedDict([])
        # todo: using station dicts means that duplicate station names are not possible

    def __repr__(self):
        """
        Provide basic network information.

        Returns
        -------
        info : str
            Basic network information.
        """
        # {self.__class__.__name__}(
        return f"Stations in '{self.name}': {list(self.stations.keys())}"

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

    def add_station(self, name, lon, lat, elev):
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
        """
        if name not in self.stations:
            self.stations[name] = Station(name,
                                          lon,
                                          lat,
                                          elev)
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
        variable : str, optional (default: None)
            Observed variable.
        depth : Depth, optional (default: None)
            Sensing depth.
        filter_dict : dict, optional (default: None)
            Metadata to use for filtering

        Yields
        ------
        station : Station
            Station.
        """
        if depth is None:
            depth = Depth(-np.inf, np.inf)

        for station in self.stations.values():

            if (variable is not None) and (variable not in station.get_variables()):
                continue # shortcut if station does not measure var

            flag = False
            for sensor in station.sensors.values():
                if (variable in [None, sensor.variable]) and \
                        depth.encloses(sensor.depth):
                    flag = True
                if flag and filter_dict:
                    f = sensor.filehandler
                    if f is None:
                        warnings.warn("No Filehandler, can't filter by metadata")

                    flag = sensor.filehandler.check_metadata(
                        variable, depth.start, depth.end,
                        filter_static_vars=filter_dict)
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


class Station(IsmnComponent):

    """
    A station is described by a distinct name and location.
    Multiple sensors at various depths can be part of a station.

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

    def __init__(self, name, lon, lat, elev):
        """
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
        """

        self.name = name
        self.lon = lon
        self.lat = lat
        self.elev = elev
        self.sensors = OrderedDict([])

    def __repr__(self):
        """
        Provide basic station information.

        Returns
        -------
        info : str
            Basic station information.
        """
        return f"Sensors at '{self.name}': {[s.name for s in self.sensors.values()]}"

    @property
    def metadata(self):
        """
        Collect the metadata from all sensors at station.
        """
        sens_meta = [s.metadata for s in self.sensors.values()]
        station_meta = MetaData().merge(sens_meta, inplace=False)
        return station_meta

    def from_file_collection(self, file_collection, keep_loaded_data=False):
        n = len(file_collection.index.size)
        d = file_collection.to_dict('list')

        assert (len(np.array(d['station'])) == 1) and \
               (len(np.array(d['network']))), \
            "Unique network and station names are expected"

        for i in range(n):
            depth = Depth(d['sensor_depth_from'][i],
                          d['sensor_depth_to'][i])

            self.add_sensor(d['instrument'][i],
                            d['variable'][i],
                            depth=depth,
                            filehandler=d['filehandler'][i],
                            name=None, # auto-generate name
                            keep_loaded_data=keep_loaded_data)

    def get_variables(self):
        """
        Get variables measured by all sensors at station.

        Returns
        -------
        variables : list
            List of variables that are observed.
        """
        #return list(np.unique(np.array(self.metadata['variable'].values())))
        return list(np.unique([s.variable for s in self.sensors.values()]))

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

    def add_sensor(self, instrument, variable, depth, filehandler,
                   name=None, keep_loaded_data=False):
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

    def __init__(self, instrument, variable, depth, name=None,
                 filehandler=None, keep_loaded_data=False):
        """
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
            File handler that allows access to observation data and
            sensor metadata (default: None).
        """

        self.instrument = instrument
        self.variable = variable
        self.depth = depth
        self.filehandler = filehandler
        self.keep_loaded_data = keep_loaded_data
        self.data = None
        self.name = name if name is not None else self.__repr__()

    def __repr__(self):
        return f"{self.instrument}_{self.variable}_" \
            f"{self.depth.start:1.6f}_{self.depth.end:1.6f}"

    @property
    def metadata(self):
        return MetaData() if self.filehandler is None else self.filehandler.metadata

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
                    variable, allowed_depth=depth,
                    filter_meta_dict=filter_meta_dict)

        return flag


if __name__ == '__main__':
    from ismn.filecollection import IsmnFileCollection

    collection = IsmnFileCollection.from_metadata_csv(
        "/home/wolfgang/code/ismn/tests/test_data/Data_seperate_files_20170810_20180809",
        "/home/wolfgang/code/ismn/tests/test_data/Data_seperate_files_20170810_20180809/python_metadata/Data_seperate_files_20170810_20180809.csv")

    s = Station('test', 1, 1, 1)
    fh = collection.files.iloc[0]['filehandler']
    d1 = Depth(0, 0.3)
    s.add_sensor('inst', 'var', d1, filehandler=fh)
    s.metadata