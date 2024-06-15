# The MIT License (MIT)
#
# Copyright (c) 2021 TU Wien
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

import os.path
import sys
from pygeogrids import CellGrid
from typing import Union

import numpy as np
import warnings
from collections import OrderedDict
import pandas as pd
from tqdm import tqdm

from ismn.meta import MetaData, Depth
from ismn.const import deprecated, CITATIONS, ismnlog
from ismn.const import xarray_available, xr

import json


class IsmnComponent:
    def _eval_xarray_installed(self):
        if not xarray_available:
            raise ImportError(
                "Please install `xarray` and `dask` with `conda install xarray "
                "dask` to use the conversion to xarray feature.")


class Sensor(IsmnComponent):
    """
    A Sensor with insitu observations.

    Attributes
    ----------
    instrument : str
        Instrument name.
    variable : str
        Observed variable.
    depth : Depth
        Sensing depth.
    name : str
        Name of the sensor.
    filehandler : IsmnFile
        File handler object to read data.
    keep_loaded_data : bool
        Keep data in memory after loading.
    data : pandas.DataFrame
        Container for data in memory (if it is being kept)
    """

    def __init__(
        self,
        instrument,
        variable,
        depth,
        name=None,
        filehandler=None,
        keep_loaded_data=False,
    ):
        """
        Initialise Sensor object.

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
        filehandler : DataFile, optional (default: None)
            File handler that allows access to observation data and
            sensor metadata via :func:`ismn.filehandlers.DataFile.read_data`.
        keep_loaded_data : bool, optional (default: False)
            Keep data for a file in memory once it is loaded. This makes subsequent
            calls of data faster (if e.g. a station is accessed multiple times)
            but can fill up memory if multiple networks are loaded.
        """
        super().__init__()

        self.instrument = instrument
        self.variable = variable
        self.depth = depth
        self.filehandler = filehandler
        self.keep_loaded_data = keep_loaded_data
        self._data = None
        self.name = name if name is not None else self.__repr__()


    def __repr__(self):
        return (f"{self.instrument}_{self.variable}_"
                f"{self.depth.start:1.6f}_{self.depth.end:1.6f}")

    @property
    def metadata(self) -> MetaData:
        return MetaData(
        ) if self.filehandler is None else self.filehandler.metadata

    @property
    def data(self):
        return self.read_data()

    def to_xarray(self) -> xr.Dataset:
        """
        Convert the Sensors data to an xarray.DataSet object
        with a location and time dimension in a single chunk.
        The dataset will contain static (dimension `sensor`) and
        dynamic (dimensions `sensor` & `time`) variables.

        Returns
        -------
        dat: xarray.DataSet
            Sensor data as xarray dataset
        """
        self._eval_xarray_installed()

        dat = self.data
        metadata = self.metadata
        
        if dat is None or metadata is None:
            return None

        data_vars = {v: (("sensor", "date_time"), np.array([dat[v].values]))
                     for v in dat.columns}
        data_vars['depth_from'] = (("sensor",), [self.depth.start])
        data_vars['depth_to'] = (("sensor",), [self.depth.end])

        for var in metadata:
            data_vars[var.name] = (("sensor",), [var.val])

        ds = xr.Dataset(
            data_vars=data_vars,
            coords={'date_time': dat.index.values},
        )

        ds['depth_from'].attrs['units'] = 'm'
        ds['depth_to'].attrs['units'] = 'm'

        return ds

    def get_coverage(self, only_good=True, start=None, end=None,
                     freq='1h'):
        """
        Estimate the temporal coverage of this sensor, i.e. the percentage
        of valid observations in the sensor time series.

        Returns
        -------
        only_good: bool, optional (default: True)
            Only consider values where the ISMN quality flag is 'G'
            as valid observations
        start: str or datetime, optional (default: None)
            Beginning of the period in which measurements are expected.
            If None, the start of the time series is used.
        end: str or datetime, optional (default: None)
            End of the period in which measurements are expected.
            If None, the start of the time series is used.
        freq: str, optional (default: '1h')
            Frequency at which the sensor is expected to take measurements.
            Most sensors in ISMN provide hourly measurements (default).
            If a different frequency is used, it must be on that
            :func:`pd.date_range` can interpret.

        Returns
        -------
        perc_coverage : float
            Data coverage of the sensor at the chosen expected measurement
            frequency within the chosen period. 0=No data, 100=no data gaps
        """
        data = self.read_data()
        if start is None:
            start = pd.Timestamp(data.index.values[0]).to_pydatetime()
        else:
            start = pd.to_datetime(start)
        if end is None:
            end = pd.Timestamp(data.index.values[-1]).to_pydatetime()
        else:
            end = pd.to_datetime(end)

        if only_good:
            data = data[data[f"{self.variable}_flag"] == 'G'].loc[:, self.variable]

        cov = (len(data.values) / len(pd.date_range(start, end, freq=freq))) * 100

        return cov

    def read_data(self):
        """
        Load data from filehandler for this Sensor by calling
        :func:`ismn.filehandlers.DataFile.read_data`.

        Returns
        -------
        data : pandas.DataFrame
            Insitu time series for this sensor, loaded from file or memory
            (if it was loaded and kept before).
        """
        if self.filehandler is None:
            ismnlog.warning(f"No filehandler found for sensor {self.name}")
        else:
            if self._data is None:
                data = self.filehandler.read_data()

                if self.keep_loaded_data:
                    self._data = data

                return data
            else:
                return self._data

    def eval(
        self,
        variable=None,
        depth=None,
        filter_meta_dict=None,
        check_only_sensor_depth_from=False,
    ):
        """
        Evaluate whether the sensor complies with the passed metadata
        requirements.

        Parameters
        ----------
        variable : str or list[str], optional (default: None)
            Check if the variable name matches, e.g. soil_moisture.
            One or multiple of :const:`ismn.const.VARIABLE_LUT`
        depth : Depth or list or tuple, optional (default: None)
            Check if the passed depth encloses the sensor depth.
            A list/tuple must contain 2 values where the first is the depth start
            and the second is the end. Start must be closer to 0 than end (or equal).
            A negative depth range is above the surface.
        filter_meta_dict : dict, optional (default: None)
            Additional metadata keys and values for which the file list is filtered
            e.g. {'lc_2010': [10, 130]} or
                 {'climate_KG': 'Dwa', 'lc_2010': [10, 130] }
            to filter for a multiple landcover classes and a climate class.
        check_only_sensor_depth_from : bool, optional (default: False)
            Ignores the sensors depth_to value and only checks if depth_from of
            the sensor is in the passed depth (e.g. for cosmic ray probes).

        Returns
        -------
        flag : bool
            Indicates whether metadata for this Sensor matches with the passed
            requirements.
        """
        if isinstance(depth, (list, tuple)):
            depth = Depth(depth[0], depth[1])

        if depth is None:
            depth = Depth(-np.inf, np.inf)

        flag = False

        if check_only_sensor_depth_from:
            d = Depth(self.depth.start, self.depth.start)
        else:
            d = self.depth

        if variable is not None:
            variable = np.atleast_1d(variable)

            if any([v in [None, self.variable] for v in variable]):
                flag = True
        else:
            flag = True

        if not depth.encloses(d):
            flag = False

        if flag and filter_meta_dict:
            if self.filehandler is None:
                warnings.warn("No filehandle found, can't filter by metadata.")
            else:
                # checks also if the metadata in file matches
                flag = self.filehandler.check_metadata(
                    variable,
                    allowed_depth=depth,
                    filter_meta_dict=filter_meta_dict,
                    check_only_sensor_depth_from=check_only_sensor_depth_from,
                )

        return flag


class Station(IsmnComponent):
    """
    A station is described by a distinct name and location.
    Multiple sensors at various depths can be part of a station.

    Attributes
    ----------
    name : str
        Station name.
    lon : float
        Longitude coordinate of station and all sensors at station.
    lat : float
        Latitude coordinate of station and all sensors at station.
    elev : float
        Elevation information of station.
    sensors : collections.OrderedDict
        Collection of Sensors and their names.
    """

    def __init__(self, name, lon, lat, elev):
        """
        Initialise Station object.

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
        super().__init__()

        self.name = name
        self.lon = lon
        self.lat = lat
        self.elev = elev
        self.sensors = OrderedDict([])

    def __repr__(self):
        # Provide basic station information.
        return f"Station '{self.name}' with Sensors: {[s.name for s in self.sensors.values()]}"

    @property
    def metadata(self) -> MetaData:
        """
        Collect the metadata from all sensors at station.
        """
        sens_meta = [s.metadata for s in self.sensors.values()]
        station_meta = MetaData().merge(sens_meta, inplace=False)
        return station_meta

    @property
    def n_sensors(self) -> int:
        """
        Number of Sensors at this Station.
        """
        return len(self.sensors)

    def __getitem__(self, item: Union[int, str]) -> Sensor:
        if isinstance(item, int):
            return self.sensors[list(self.sensors.keys())[item]]
        else:
            return self.sensors[item]

    def to_xarray(self, **filter_kwargs) -> xr.Dataset:
        """
        Collect all sensor data at this station into a xarray.DataSet object
        with a location and time dimension in a single chunk.
        The dataset will contain static (dimension `sensor`) and
        dynamic (dimensions `sensor` & `time`) variables.

        Parameters
        ----------
        filter_kwargs: optional
            Filter sensors at the station to include in the dataset
            (variable, depth, etc.).
            For a description of possible filter kwargs, see
            :func:`ismn.components.Sensor.eval`

        Returns
        -------
        dat: xarray.DataSet
            Sensor data as xarray dataset
        """
        self._eval_xarray_installed()

        station = []
        for sensor in self.iter_sensors(**filter_kwargs):
            s = sensor.to_xarray()
            if s is not None:
                station.append(s)

        if len(station) == 0:
            return None
        
        station = xr.concat(station, dim='sensor')

        if 'depth_from' in station.attrs:
            station.attrs.pop('depth_from')
        if 'depth_to' in station.attrs:
            station.attrs.pop('depth_to')

        station.attrs['station_name'] = self.name
        station.attrs['lat'] = self.lat
        station.attrs['lon'] = self.lon
        station.attrs['n_sensors'] = len(station['sensor'].values)

        return station


    def get_variables(self):
        """
        Get variables measured by all sensors at station.

        Returns
        -------
        variables : list
            List of variables that are observed.
        """
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

    def get_min_max_obs_timestamp(self,
                                  variable="soil moisture",
                                  min_depth=None,
                                  max_depth=None):
        """
        Goes through the sensors associated with this station
        and checks the metadata to get and approximate time coverage of the station.
        This is just an overview. If holes have to be detected the
        complete file must be read.

        Parameters
        ----------
        variable: str, optional (default: 'soil_moisture')
            name of the variable, only sensors measuring that variable are used.
        min_depth : float, optional (default: None)
            depth_from of variable has to be >= min_depth in order to be
            included.
        max_depth : float, optional (default: None)
            depth_to of variable has to be <= max_depth in order to be
            included.

        Returns
        -------
        start_date: datetime.datetime
            Earliest date observed by any sensor at the station after filtering
            for the passed requirements.
        end_date: datetime.datetime
            Latest date observed by any sensor at the station after filtering
            for the passed requirements.
        """
        depth = Depth(
            -np.inf if min_depth is None else min_depth,
            np.inf if max_depth is None else max_depth,
        )

        min_from, max_to = None, None

        for sensor in self.iter_sensors(variable=variable, depth=depth):
            time_from = sensor.metadata["timerange_from"].val
            time_to = sensor.metadata["timerange_to"].val
            if (min_from is None) or (time_from < min_from):
                min_from = time_from
            if (max_to is None) or (time_to > max_to):
                max_to = time_to

        min_from = min_from.to_pydatetime() if min_from is not None else None
        max_to = max_to.to_pydatetime() if max_to is not None else None

        return min_from, max_to

    def add_sensor(
        self,
        instrument,
        variable,
        depth,
        filehandler=None,
        name=None,
        keep_loaded_data=False,
    ):
        """
        Add a new Sensor to this Station.

        Parameters
        ----------
        instrument : str
            Instrument name. e.g. ThetaProbe-ML2X
        variable : str
            Observed variable. e.g. soil_moisture
        depth : Depth
            Sensing depth. e.g. Depth(0, 0.1)
        filehandler : DataFile, optional (default: None)
            File handler object that allows access to observation data and
            sensor metadata via its read_data() function (default: None).
        name: str or int, optional (default: None)
            A name or id for the sensor. If None is passed, one is generated
            automatically from other properties.
        keep_loaded_data : bool, optional (default: False)
            Keep data for the sensor in memory once it is loaded.
            This makes subsequent reading of the same data faster
            but can fill up memory if stations / networks are loaded.
        """
        if name is None:
            name = f"{instrument}_{variable}_{depth.start:1.6f}_{depth.end:1.6f}"

        if name not in self.sensors:
            self.sensors[name] = Sensor(
                instrument,
                variable,
                depth,
                name,
                filehandler,
                keep_loaded_data,
            )
        else:
            ismnlog.warning(f"Sensor already exists: {name}")

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
            ismnlog.warning(f"Sensor not found: {name}")

    def iter_sensors(self, **filter_kwargs):
        """
        Iterates over all sensors in this station and yields those that
        comply with the passed filter settings (or all).

        Parameters
        ----------
        Keyword arguments are used to check all sensors at all stations,
        only stations that have at least one matching sensor are returned.

        For a description of possible filter kwargs, see
        :func:`ismn.components.Sensor.eval`

        Yields
        ------
        sensors : Sensor
            (Filtered) Sensors at the Station.
        """
        for sensor in self.sensors.values():
            if sensor.eval(**filter_kwargs):
                yield sensor

    @deprecated  # use iter_sensors instead
    def get_sensors(self, variable, depth_from, depth_to):
        """
        get the sensors at which the variable was measured at the
        given depth

        Parameters
        ----------
        variable : str
            variable abbreviation
        depth_from : float
            shallower depth of layer the variable was measured at
        depth_to : float
            deeper depth of layer the variable was measured at
        Returns
        -------
        sensors : numpy.ndarray
            array of sensors found for the given combination of variable and depths
        """
        return np.array([
            s for s in self.iter_sensors(
                variable=variable, depth=Depth(depth_from, depth_to))
        ])


class Network(IsmnComponent):
    """
    A network is described by a distinct name and can be composed of
    multiple ISMN stations.

    Attributes
    ----------
    name : str
        Network name.
    stations : OrderedDict[name, Station]
        Stations belonging to the network. If a string is passed, an empty
        Station is added.
    """

    def __init__(self, name, stations=None):
        """
        Initialise Network object.

        Parameters
        ----------
        name : str
            Network name.
        stations : list[Station], optional (default: None)
            Initial list of Station object to fill the network with.
            Additional Stations can be added later.
        """
        super().__init__()

        # todo: using station dicts means that duplicate station names are not possible
        self.name = name
        self.stations = OrderedDict([])

        if stations is not None:
            for s in stations:
                self.stations[s.name] = s

    def __repr__(self):
        # Provide basic Network information.
        return f"Network '{self.name}' with Stations: {list(self.stations.keys())}"

    def __getitem__(self, item: Union[int, str]):
        # shortcut to access networks directly
        if isinstance(item, int):
            item = list(self.stations.keys())[item]
        return self.stations[item]

    @property
    def coords(self) -> (list, list):
        """
        Get lists of lats and lons for all stations in the Network
        """
        lons, lats = [], []
        for name, stat in self.stations.items():
            lons.append(stat.lon)
            lats.append(stat.lat)

        return lons, lats

    @property
    def grid(self) -> CellGrid:
        """
        Get grid for all Stations in Network
        """
        # CellGrid is intentional
        lons, lats = self.coords
        return CellGrid(lons, lats, cells=np.full(len(lons), 0))

    @property
    def n_stations(self) -> int:
        """
        Number of Stations in this Network.
        """
        return len(self.stations)

    def to_xarray(self, **filter_kwargs) -> xr.Dataset:
        """
        Collect all sensor data at this station into a xarray.DataSet object
        with a location and time dimension in a single chunk.
        The dataset will contain static (dimension `sensor`) and
        dynamic (dimensions `sensor` & `time`) variables.

        Parameters
        ----------
        filter_kwargs: optional
            Filter sensors in the network to include in the dataset
            (variable, depth, etc.).
            For a description of possible filter kwargs, see
            :func:`ismn.components.Sensor.eval`

        Returns
        -------
        dat: xarray.DataSet
            Sensor data as xarray Dataset. Sensor data are stored as
            Dask arrays!
        """
        self._eval_xarray_installed()

        net = []
        for station in tqdm(self.iter_stations(), total=self.n_stations):
            s = station.to_xarray(**filter_kwargs)
            # store sensor data as dask arrays to save memory
            s = s.chunk(dict(date_time=None, sensor=1))
            if s is not None:
                net.append(s)

        if len(net) == 0:
            return None

        n_stations = len(net)

        net = xr.concat(net, dim="sensor")

        net.attrs.pop('station_name')
        net.attrs.pop('lat')
        net.attrs.pop('lon')
        net.attrs['n_sensors'] = len(net['sensor'].values)
        net.attrs['n_stations'] = n_stations
        net.attrs['network'] = self.name

        return net

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
            self.stations[name] = Station(name, lon, lat, elev)
        else:
            ismnlog.warning(f"Station already exists: {name}")

    def remove_station(self, name):
        """
        Remove Station from Network.

        Parameters
        ----------
        name : str
            Station name.
        """
        if name in self.stations:
            del self.stations[name]
        else:
            ismnlog.warning(f"Station not found {name}")

    def iter_stations(self, **filter_kwargs):
        """
        Get all stations having at least one sensor observing
        a specific variable and/or sensing depth.

        Parameters
        ----------
        Parameters are used to check all sensors at all stations, only stations
        that have at least one matching sensor are returned.

        For a description of possible filter kwargs, see
        :func:`ismn.components.Sensor.eval`

        Yields
        ------
        station : Station
            Stations that contain at least one sensor that matches to the passed
            conditions.
        """
        for station in self.stations.values():
            for sensor in station.sensors.values():
                if sensor.eval(**filter_kwargs):
                    yield station
                    break

    def iter_sensors(self, **filter_kwargs):
        """
        Get all sensors in all stations in the network that comply with the
        passed filtering parameters.

        Parameters
        ----------
        Keyword arguments are used to evaluate the sensors, see
        :func:`ismn.components.Sensor.eval`

        Yields
        ------
        station : Station
            Station that contains Sensor.
        sensor : Sensor
            Sensor at Station that matches to the passed filtering conditions.
        """
        for station in self.stations.values():
            for sensor in station.sensors.values():
                if sensor.eval(**filter_kwargs):
                    yield station, sensor

    def get_citations(self):
        """
        Return reference(s) for this network. Users of ISMN should cite the
        networks they are using in a publication. This information can also
        be found on the ISMN website.

        Returns
        -------
        references : list
            A list of references / citations / acknowledgements for this Network.
        """
        try:
            refs = CITATIONS[self.name]
        except KeyError:
            refs = [f"No reference(s) for network {self.name} available."]

        return refs


class NetworkCollection(IsmnComponent):
    """
    A NetworkCollection holds multiple networks and provides functionality
    to perform access to components from multiple networks.
    A grid is added that contains all stations to perform spatial searches.

    Attributes
    ----------
    networks : OrderedDict
        Collection of network names and Networks
    grid : CellGrid
        Grid that contains one point for each station in all networks.
    """

    def __init__(self, networks):
        """
        Create network collection from previously created Networks.

        Parameters
        ----------
        networks : list[Network]
            List of Networks that build the collection.
        """
        super().__init__()

        self.networks = OrderedDict([])

        lons = []
        lats = []
        for net in networks:
            self.networks[net.name] = net
            net_lons, net_lats = net.coords
            lons += net_lons
            lats += net_lats

        if (len(lons) > 0) and (len(lats) > 0):
            # Should be CellGrid
            self.grid: CellGrid = CellGrid(lons, lats, cells=np.full(len(lons), 0))
        else:
            self.grid = None

    def __repr__(self, indent: str = ""):
        return ",\n".join([
            f"{indent}{net.name}: {list(net.stations.keys())}"
            for net in self.networks.values()
        ])

    def __getitem__(
            self, item: Union[int, str,
                              list]) -> Union["NetworkCollection", Network]:
        # shortcut to access networks directly
        if isinstance(item, (int, str)):
            if isinstance(item, int):
                item = list(self.networks.keys())[item]
            net: Network = self.networks[item]
            return net
        else:
            keys = list(self.networks.keys())
            sub: NetworkCollection = NetworkCollection(networks=[
                self.networks[n] if isinstance(n, str) else self
                .networks[keys[n]] for n in item
            ])
            return sub

    def iter_networks(self) -> Network:
        """
        Iterate through all networks in the Collection.
        """
        for nw in self.networks.values():
            yield nw

    def iter_stations(self, **filter_kwargs) -> (Network, Station):
        """
        Iterate through Networks in the Collection and get (all/filtered)
        Stations.
        """
        for nw in self.networks.values():
            for stat in nw.iter_stations(**filter_kwargs):
                yield nw, stat

    def iter_sensors(self, **filter_kwargs) -> (Network, Station, Sensor):
        """
        Iterate through Networks in the Collection and get (all/filtered)
        Stations and Sensors at each Station.
        """
        for nw in self.networks.values():
            for stat, sen in nw.iter_sensors(**filter_kwargs):
                yield nw, stat, sen

    def station4gpi(self, gpi):
        """
        Get the Station for the passed gpi in the grid.

        Parameters
        ----------
        gpi : int or list[int]
            Point index or multiple indices in self.grid.

        Returns
        -------
        station : Station or list[Station]
            Station(s) at gpi(s).
        """
        idxs = np.atleast_1d(gpi)
        in_grid = np.isin(idxs, self.grid.activegpis)

        if not all(in_grid):
            raise ValueError(
                f"Index not found in loaded grid: {idxs[~in_grid]}")

        lon, lat = self.grid.gpi2lonlat(idxs)

        stations = []
        for net, stat in self.iter_stations():
            if (stat.lon == lon) and (stat.lat == lat):
                stations.append(stat)
                if len(stations) == len(idxs):
                    break  # stop when all indices are found

        return stations[0] if len(stations) == 1 else stations

    def get_nearest_station(self, lon, lat, max_dist=np.inf):
        """
        Get nearest station for given longitude/latitude coordinates.

        Parameters
        ----------
        lon : float or list[float]
            Longitude coordinate(s).
        lat : float or list[float]
            Latitude coordinate(s).
        max_dist : float, optional (default: np.Inf)
            Maximum search distance.

        Returns
        -------
        station : Station or list[Station]
            The nearest Station(s) to the passed coordinates.
        dist : float or list[float]
            Distance in meter between the passed coordinates and the
            actual location of the station.
        """
        gpi, dist = self.grid.find_nearest_gpi(lon, lat, max_dist=max_dist)
        station = self.station4gpi(gpi)

        return station, dist

    def export_citations(self, out_file=None):
        """
        Returns the references for all networks in the collection.
        Optionally, they are also written to file.
        Information on how to correctly cite ISMN networks can be found
        on the ISMN website.

        Parameters
        ----------
        out_file : str, optional (default: None)
            If a path is passed here, a new file will be generated with all
            references for the  current collection.

        Returns
        -------
        references: OrderedDict
            Network names as keys and network references as values
        """
        refs = OrderedDict([
            (net.name, net.get_citations()) for net in self.iter_networks()
        ])

        if out_file is not None:
            with open(out_file, mode="w") as out_file:
                for name, reflist in refs.items():
                    out_file.write(f"References for Network {name}:\n")
                    out_file.write(
                        "-----------------------------------------\n")
                    for ref in reflist:
                        out_file.write(f"{ref}\n")
                        out_file.write("\n")
                    out_file.write("\n")

        return refs

    def export_geojson(self, path, network=True, station=True, sensor=False,
                       depth=True, timerange=True, extra_props=None,
                       filter_kwargs=None):
        """
        Filter sensors in collection and create geojson file containing all
        features.

        Parameters
        ----------
        path: str
            Path to geojson file
        network: bool, optional (default: True)
            If True, network names are included in geojson file
        station: bool, optional (default: True)
            If True, station names are included in geojson file
        sensor: bool, optional (default: False)
            If True, sensor names are included in geojson file
        depth: bool, optional (default: True)
            If True, depth_from and depth_to are included in geojson file
        timerange: bool, optional (default: True)
            If True, timerange_from and timerange_to are included in geojson
        extra_props: list[str], optional (default: None)
            List of extra properties from sensor metadata to include in
            geojson file
            By default only depth_from and depth_to are included
            e.g. ['variable', 'frm_class'] etc.
        filter_kwargs: dict, optional (default: None)
            Keyword arguments to filter sensors in collection before extracting
            metadata.
            see :func:`ismn.components.Sensor.eval`
        """
        extra_props = extra_props or []
        geoinfo = {
            "type": "FeatureCollection",
            "features": [],
        }

        filter_kwargs = filter_kwargs or dict()

        for nw, stat, sens in self.iter_sensors(**filter_kwargs):
            feature = {
                "type": "Feature",
                "geometry": {
                    "type": "Point",
                    "coordinates": [
                        stat.lon,
                        stat.lat
                    ],
                },
                "properties": {
                    "datasetProperties": []
                }
            }
            if network:
                feature["properties"]["datasetProperties"] += [
                    {
                        "propertyName": "network",
                        "propertyValue": nw.name
                    }
                ]
            if station:
                feature["properties"]["datasetProperties"] += [
                    {
                        "propertyName": "station",
                        "propertyValue": stat.name
                    }
                ]
            if sensor:
                feature["properties"]["datasetProperties"] += [
                    {
                        "propertyName": "sensor",
                        "propertyValue": sens.name
                    }
                ]
            if depth:
                feature["properties"]["datasetProperties"] += [
                    {
                        "propertyName": "depth_from",
                        "propertyValue": str(sens.depth[0])
                    },
                    {
                        "propertyName": "depth_to",
                        "propertyValue": str(sens.depth[1])
                    }
                ]
            if timerange:
                feature["properties"]["datasetProperties"] += [
                    {
                        "propertyName": "timerange_from",
                        "propertyValue": str(sens.metadata["timerange_from"].val)
                    },
                    {
                        "propertyName": "timerange_to",
                        "propertyValue": str(sens.metadata["timerange_to"].val)
                    }
                ]
            for prop in extra_props:
                if prop not in sens.metadata:
                    raise KeyError(f"No sensor property '{prop}' found. "
                                   f"Choose one of {sens.metadata.keys()}.")
                feature["properties"]["datasetProperties"] += [
                    {
                        "propertyName": prop,
                        "propertyValue": str(sens.metadata[prop].val),
                    }
                ]

            geoinfo["features"].append(feature)

        with open(path, 'w') as f:
            json.dump(geoinfo, f, ensure_ascii=False)
