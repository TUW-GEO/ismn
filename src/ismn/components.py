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
from typing import Union

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

class NetworkCollection(IsmnComponent):
    """
    A NetworkCollection holds multiple networks and provides functionality
    to perform access to components from multiple networks.
    A grid is added that contains all stations to perform spatial searches.
    """
    def __init__(self, networks):

        """
        Create network collection from previously loaded IsmnFileCollection.

        Parameters
        ----------
        networks : list[Network or str]
            List of Networks that build the collection
        """

        self.networks = OrderedDict([])
        lons = []
        lats = []
        for net in networks:
            self.networks[net.name] = net
            net_lons, net_lats = net.coords
            lons += net_lons
            lats += net_lats

        self.grid = BasicGrid(lons, lats)

    def __repr__(self):
        strs = []
        for net in self.networks.values():
            strs.append(f"{net.name}: {list(net.stations.keys())}")
        return ', \n'.join(strs)

    def __getitem__(self, item:Union[int,str]):
        if isinstance(item, int):
            item = list(self.networks.keys())[item]
        return self.networks[item]

    def iter_networks(self):
        # Iterate through all networks
        for network in self.networks.values():
            yield network

    def station4idx(self, idx):
        """
        Get the station object for the passed gpi.

        Parameters
        ----------
        idx : int
            Point index in self.grid, resp. line index in file list.

        Returns
        -------
        station : Station
            Station at gpi.
        """
        if idx not in self.grid.activegpis:
            raise ValueError("Index does not exist in loaded grid")

        lon, lat = self.grid.gpi2lonlat(idx)
        for net in self.iter_networks():
            for stat in net.iter_stations():
                if (stat.lon == lon) and (stat.lat == lat):
                    return stat

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
        gpi, dist = self.grid.find_nearest_gpi(lon, lat, max_dist=max_dist)
        station = self.station4idx(gpi)

        return station, dist

    # def get_sensors(self, network=None, station=None, variable=None,
    #                 depth=None, filter_static_vars=None):
    #     """
    #     Yield all sensors for a specific network and/or station and/or
    #     variable and/or depth.
    #
    #     Parameters
    #     ----------
    #     network : str, optional
    #         Network name (default: None).
    #     station : str, optional
    #         Station name (default: None).
    #     variable : str, optional
    #         Variable name (default: None).
    #     depth : Depth, optional (default: None)
    #         Sensing depth.
    #     filter_static_vars: dict, optional (default: None)
    #         Additional metadata keys and values for which the file list is filtered
    #         e.g. {'lc_2010': 10} to filter for a landcover class.
    #         if there are multiple conditions, ALL have to be fulfilled.
    #         e.g. {'lc_2010': 10', 'climate_KG': 'Dfc'})
    #
    #     Yield
    #     -----
    #     sensor : Sensor
    #         Sensor.
    #     """
    #     if depth is None:
    #         depth = Depth(-np.inf, np.inf)
    #
    #     for n in self.networks.values():
    #         if network not in [None, n.name]:
    #             continue
    #         for st in n.stations.values():
    #             if station not in [None, st.name]:
    #                 continue
    #             for se in st.sensors.values():
    #                 if se.eval(variable, depth, filter_static_vars):
    #                     yield se

class Network(IsmnComponent):
    """
    A network is described by a distinct name and can be composed of
    multiple stations.

    Attributes
    ----------
    name : str
        Network name.
    stations : List[Station]
        Stations belonging to the network. If a string is passed, an empty
        station is added.
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
                 name,
                 stations=None):
        """
        Parameters
        ----------
        name : str
            Network name.
        """
        self.name = name
        self.stations = OrderedDict([])

        if stations is not None:
            for s in stations:
                self.stations[s.name] = s

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

    def add_sensor(self, instrument, variable, depth, filehandler=None,
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
        filehandler : IsmnFile, optional (default: None)
            File handler.
        name: str or int, optional (default: None)
            A name or id for the sensor. If None is passed, one is generated.
        keep_loaded_data : bool, optional (default: False)
            Keep data for a file in memory once it is loaded. This makes subsequent
            calls of data faster (if e.g. a station is accessed multiple times)
            but can fill up memory if multiple networks are loaded.
        """
        if name is None:
            name = f"{instrument}_{variable}_{depth.start:1.6f}_{depth.end:1.6f}"

        if name not in self.sensors:
            self.sensors[name] = Sensor(instrument, variable, depth, name,
                                        filehandler, keep_loaded_data)
        else:
            logger.warning(f'Sensor already exists: {name}')

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
        keep_loaded_data : bool, optional (default: False)
            Keep data for a file in memory once it is loaded. This makes subsequent
            calls of data faster (if e.g. a station is accessed multiple times)
            but can fill up memory if multiple networks are loaded.
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

    networks = []

    for n_id in range(10):
        stations = []
        for s_id in range(100):
            station = Station(f'station{s_id}', *np.random.rand(3)*100)
            for sens_id in range(10):
                station.add_sensor(f'Instrument{sens_id}', 'var', Depth(1,2))
            stations.append(station)
        networks.append(Network(name=f"Network-{n_id}", stations=stations))

    coll = NetworkCollection(networks)

    coll.get_nearest_station(10,10)

