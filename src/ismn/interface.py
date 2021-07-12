# -*- coding: utf-8 -*-
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

from pathlib import Path
from tempfile import gettempdir
import platform
import os
import sys
from typing import Optional, Union

import pandas as pd

from ismn.filecollection import IsmnFileCollection, _load_metadata_df
from ismn.components import *
from ismn.const import *
from ismn.base import IsmnRoot

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    if platform.system() == "Darwin":
        import matplotlib

        matplotlib.use("TkAgg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    plotlibs = True
except ImportError:
    plotlibs = False


class ISMN_Interface:
    """
    Class provides interface to ISMN data downloaded from the ISMN website
    https://ismn.earth.
    Upon initialization it collects metadata from all files in
    path_to_data and saves metadata information in a csv file into the folder
    python_metadata in meta_path (or data_path if no meta_path is defined).
    First initialization can take some time if all ISMN
    data is present in data_path and will start multiple processes.

    Parameters
    ----------
    data_path : str or Path
        Path to ISMN data to read, either to a zip archive or to the extracted
        directory that contains the network folders.
        Download data from https://ismn.earth after registration.
    meta_path : str or Path
        Path where the metadata csv file(s) is / are stored. The actual filename
        is defined by the name of data_path and will be generated automatically
        if it does not yet exist.
    network : str or list, optional (default: None)
        Name(s) of network(s) to load. Other data in the data_path will be ignored.
        By default or if None is passed, all networks are activated.
        If an empty list is passed no networks are activated.
    parallel: bool, optional (default: False)
        Activate parallel processes to speed up metadata generation.
        All available CPUs will be used.
    keep_loaded_data : bool, optional (default: False)
        Keep data for a file in memory once it is loaded. This makes subsequent
        calls of data faster (if e.g. a station is accessed multiple times)
        but can fill up memory if multiple networks are loaded.

    Raises
    ------
    ISMNError
        if given, network was not found :attr:`.ISMN_Interface.data_path`

    Attributes
    ----------
    climate : collections.OrderedDict
        All Climate classes and their descriptions.
    collection : NetworkCollection
        Contains all loaded networks with stations and sensors.
    keep_loaded_data : bool
        Switch to keep data in memory after loading (not recommended).
    metadata : pandas.DataFrame
        Metadata for active networks, with idx that could also be passed
        to :func:`ismn.interface.read_metadata`
    landcover : collections.OrderedDict
        All Landcover classes and their descriptions.
    parallel : bool
        Switch to activate parallel processing where possible.
    root : IsmnRoot
        ISMN data folder or .zip access


    Properties
    ----------
    networks : OrderedDict
        Access Networks container from collection directly.
    grid : pygeogrids.grid.BasicGrid
        Grid from collection that contains all station lats and lons
    """

    def __init__(
        self,
        data_path,
        meta_path=None,
        network=None,
        parallel=False,
        keep_loaded_data=False,
        temp_root=gettempdir(),
    ):

        self.climate, self.landcover = KOEPPENGEIGER, LANDCOVER
        self.parallel = parallel

        self.root = IsmnRoot(data_path)

        self.keep_loaded_data = keep_loaded_data

        self.activate_network(network=network, meta_path=meta_path, temp_root=temp_root)

    def activate_network(
        self,
        network: Union[list, str] = None,
        meta_path: str = None,
        temp_root: str = gettempdir(),
    ):
        """
        Load (file) collection for specific networks.
        """

        meta_csv_filename = f"{self.root.name}.csv"

        if meta_path is None:
            meta_path = Path(self.root.root_dir) / "python_metadata"
        else:
            meta_path = Path(meta_path)

        meta_csv_file = meta_path / meta_csv_filename

        if network is not None:
            network = np.atleast_1d(network)

        if not os.path.isfile(meta_csv_file):
            self.__file_collection = IsmnFileCollection.build_from_scratch(
                self.root,
                parallel=self.parallel,
                log_path=meta_path,
                temp_root=temp_root,
            )
            self.__file_collection.to_metadata_csv(meta_csv_file)

        self.__file_collection = IsmnFileCollection.from_metadata_csv(
            self.root, meta_csv_file, network=network
        )

        metadata = self.__file_collection.metadata_df

        metadata.index = range(len(metadata.index))
        self.metadata = metadata

        networks = self.__collect_networks(network)
        self.collection = NetworkCollection(networks)

    def __collect_networks(self, network_names: Optional[list] = None) -> list:
        """
        Build Networks and fill them with Stations and Sensors and apply
        according filehandlers from filelist for data reading.
        """
        networks = OrderedDict([])

        for f in self.__file_collection.iter_filehandlers(network_names):

            nw_name, st_name, instrument = (
                f.metadata["network"].val,
                f.metadata["station"].val,
                f.metadata["instrument"].val,
            )

            if nw_name not in networks:
                networks[nw_name] = Network(nw_name)

            if st_name not in networks[nw_name].stations:
                networks[nw_name].add_station(
                    st_name,
                    f.metadata["longitude"].val,
                    f.metadata["latitude"].val,
                    f.metadata["elevation"].val,
                )

            # the senor name is the index in the list
            networks[nw_name].stations[st_name].add_sensor(
                instrument,
                f.metadata["variable"].val,
                f.metadata["variable"].depth,
                filehandler=f,  # todo: remove station meta from sensor
                name=None,
                keep_loaded_data=self.keep_loaded_data,
            )

        return list(networks.values())

    def __getitem__(self, item):
        return self.collection[item]

    def __repr__(self):
        return (
            f"{self.root}\n"
            f"with Networks[Stations]:\n"
            f"------------------------\n"
            f"{self.collection.__repr__('  ')}"
        )

    @property
    def networks(self):
        return self.collection.networks

    @property
    def grid(self):
        return self.collection.grid

    @deprecated
    def list_networks(self) -> np.ndarray:
        # get network names from list of active files
        return np.array(list(self.networks.keys()))

    @deprecated
    def list_stations(self, network: str = None) -> np.ndarray:
        # get station names for one of the active networks
        if network is not None:
            if network not in self.networks:
                raise ISMNError(
                    f"Network {network} not found in currently loaded networks."
                )
            return np.array(list(self.networks[network].stations.keys()))
        else:
            stations = []
            for network in self.networks.values():
                stations += list(network.stations.keys())
            return np.array(stations)

    @deprecated
    def list_sensors(self, network: str = None, station: str = None) -> np.ndarray:
        # List sensors names for a specific sensor in an active network
        sensors = np.array([])
        for net in self.networks.values():
            if network in [None, net.name]:
                for stat in net.stations.values():
                    if station in [None, stat.name]:
                        sensors = np.append(
                            sensors,
                            # get the objects instead, use .values()?
                            np.array(list(stat.sensors.keys())),
                        )

        return sensors

    def network_for_station(self, stationname, name_only=True):
        """
        Find networks that contain a station of the passed name.

        Parameters
        ----------
        stationname : str
            Station name to search in the active networks.
        name_only : bool, optional (default: True)
            Returns only the name of the network and not the Network.

        Returns
        -------
        network_names : str or Network or None
            Network that contains a station of that name, or None if no such
            network exists.
            Prints are warning and uses the FIRST found network name if there
            are multiple stations with the same name in different networks.

        """
        network_with_station = []

        for network in self.networks.values():
            if stationname in network.stations.keys():
                network_with_station.append(network)

        if len(network_with_station) > 1:
            warnings.warn("stationname occurs in multiple networks")

        if len(network_with_station) == 0:
            return None
        else:
            nw = network_with_station[0]
            if name_only:
                warnings.warn(
                    "Future Versions of the package will always return the Network object (same as name_only=False now). "
                    "You can use Network.name to get the name of a network.",
                    category=DeprecationWarning,
                )
                return nw.name
            else:
                return nw

    def stations_that_measure(self, variable, **filter_kwargs):
        """
        Goes through all stations and returns those that measure the specified
        variable

        Parameters
        ----------
        variable : str
            variable name, one of:
                * soil_moisture
                * soil_temperature
                * soil_suction
                * precipitation
                * air_temperature
                * field_capacity
                * permanent_wilting_point
                * plant_available_water
                * potential_plant_available_water
                * saturation
                * silt_fraction
                * snow_depth
                * sand_fraction
                * clay_fraction
                * organic_carbon
                * snow_water_equivalent
                * surface_temperature
                * surface_temperature_quality_flag_original
        filter_kwargs :
            Parameters are used to check all sensors at all stations, only stations
            that have at least one matching sensor are returned.
            For a description of possible filter kwargs, see
            :func:`ismn.components.Sensor.eval`

        Yields
        -------
        ISMN_station : Station
        """
        for network in self.networks.values():
            for station in network.iter_stations(variable=variable, **filter_kwargs):
                yield station

    def get_dataset_ids(
        self,
        variable,
        min_depth=0,
        max_depth=0.1,
        filter_meta_dict=None,
        check_only_sensor_depth_from=False,
        groupby=None,
    ):
        """
        Yield all sensors for a specific network and/or station and/or
        variable and/or depth. The id is defined by the position of the filehandler
        in the filelist.

        Parameters
        ----------
        variable : str or list[str] or None
            Variable(s) to filer out, None to allow all variables.
        min_depth : float, optional (default: 0)
            Min depth of sensors to search
        max_depth : float, optional (default: 0.1)
            Max depth of sensors to search
        filter_meta_dict: dict, optional (default: None)
            Additional metadata keys and values for which the file list is filtered
            e.g. {'lc_2010': 10} to filter for a landcover class.
            if there are multiple conditions, ALL have to be fulfilled.
            e.g. {'lc_2010': 10', 'climate_KG': 'Dfc'})
        check_only_sensor_depth_from : bool, optional (default: False)
            Ignores the sensors depth_to value and only checks if depth_from of
            the sensor is in the passed depth (e.g. for cosmic ray probes).
        groupby : str, optional (default: None)
            A metadata field name that is used to group sensors, e.g. network
        """
        if groupby is None:
            ids = []
        else:
            ids = {}

        depth = Depth(min_depth, max_depth)

        for id, filehandler in enumerate(self.__file_collection.iter_filehandlers()):
            eval = filehandler.check_metadata(
                variable=variable,
                allowed_depth=depth,
                filter_meta_dict=filter_meta_dict,
                check_only_sensor_depth_from=check_only_sensor_depth_from,
            )

            if eval:
                if groupby is not None:
                    groupval = filehandler.metadata[groupby].val
                    if groupval not in ids.keys():
                        ids[groupval] = []
                    ids[groupval].append(id)
                else:
                    ids.append(id)

        return ids

    def read_metadata(self, idx, format="pandas"):
        """
        Read only metadata by id as pd.DataFrame.

        Parameters
        ----------
        idx : int or list
            id of sensor to read, best one of those returned
            from :func:`ismn.interface.get_dataset_ids` or one in
            :attr:`.ISMN_Interface.metadata`.
        format : str, optional (default: 'pandas')
            This only affects the return value when a SINGLE idx is passed.
            If multiple indices or None is passed, a DataFrame is returned.
                - pandas : return metadata as dataframe (Default)
                - dict : return metadata as dict (only for single idx)
                - obj : return metadata as MetaData object (only for single idx)

        Returns
        -------
        metadata : pd.DataFrame or dict or MetaData
            Metadata for the passed index.
        """
        idx = np.atleast_1d(idx)

        if len(idx) == 1:
            filehandler = self.__file_collection.get_filehandler(idx[0])
            if format.lower() == "pandas":
                return filehandler.metadata.to_pd()
            elif format.lower() == "dict":
                return filehandler.metadata.to_dict()
            elif format.lower() == "obj":
                return filehandler.metadata
            else:
                raise NotImplementedError(f"{format} is not a supported format.")
        else:
            if format.lower() != "pandas":
                warnings.warn(
                    "Multiple indices passed (or None), return format will be 'pandas'"
                )

            dfs = []
            for i in idx:
                filehandler = self.__file_collection.get_filehandler(i)
                if len(idx) == 1:
                    return filehandler.metadata.to_pd()
                else:
                    df = filehandler.metadata.to_pd(transpose=True, dropna=False)
                    df.index = [i]
                    dfs.append(df)

            return pd.concat(dfs, axis=0).dropna(axis=1, how="all")

    def read_ts(self, idx, return_meta=False):
        """
        Read a time series directly by the filehandler id.

        Parameters
        ----------
        idx : int
            id of filehandler to read, best one of those returned
            by :func:`ismn.interface.ISMN_Interface.get_dataset_ids`
        return_meta : bool, optional (default: False)
            Also return the metadata for this sensor (as a second return value)

        Returns
        -------
        timeseries : pd.DataFrame
            Observation time series
        metadata : pd.DataFrame, optional
            All available metadata for that sensor. Only returned when
            `return_meta=False`
        """

        filehandler = self.__file_collection.get_filehandler(idx)
        if return_meta:
            return filehandler.read_data(), filehandler.metadata.to_pd()
        else:
            return filehandler.read_data()

    def read(self, *args, **kwargs):
        # alias of :func:`ismn.interface.ISMN_Interface.read_ts`
        return self.read_ts(*args, **kwargs)

    def find_nearest_station(self, lon, lat, return_distance=False, max_dist=np.inf):
        """
        Finds the nearest station to passed coordinates available in downloaded
        data.

        Parameters
        ----------
        lon : float
            Longitude of point
        lat : float
            Latitude of point
        return_distance : bool, optional (default: False)
            if True also distance is returned
        max_dist : float, optional (default: np.inf)
            Maximum distance allowed. If no station is within this distance
            None is returned.

        Returns
        -------
        station : ismn.components.Station
            Nearest station object that was found in within the selected distance
        distance : float, optional
            distance to station in meters, measured in cartesian coordinates and not on
            a great circle. Should be OK for small distances
        """
        # what happens if there is no point within max dist if that works?
        gpi, d = self.collection.grid.find_nearest_gpi(lon, lat, max_dist=max_dist)

        if len(np.atleast_1d(gpi)) == 0:
            stat = None
            d = None
        else:
            stat = self.collection.station4gpi(gpi)

        if return_distance:
            return stat, d
        else:
            return stat

    def plot_station_locations(
        self,
        variable=None,
        min_depth=-np.inf,
        max_depth=np.inf,
        stats_text=True,
        check_only_sensor_depth_from=False,
        markersize=1,
        filename=None,
        ax=None,
    ):
        """
        Plots available stations on a world map in robinson projection.

        Parameters
        ----------
        variable : str or list[str], optional (default: None)
            Show only stations that measure this/these variable(s), e.g. soil_moisture
            If None is passed, no filtering for variable is performed.
        min_depth : float, optional (default: -np.inf)
            Minimum depth, only stations that have a valid sensor measuring the
            passed variable (if one is selected) in this depth range are included.
        max_depth : float, optional (default: -np.inf)
            See description of min_depth. This is the bottom threshold for the
            allowed depth.
        stats_text : bool, optianal (default: False)
            Include text of net/stat/sens counts in plot.
        check_only_sensor_depth_from : bool, optional (default: False)
            Ignores the sensors depth_to value and only checks if depth_from of
            the sensor is in the passed depth_range (e.g. for cosmic ray probes).
        markersize : int, optional (default: 1)
            Size of the marker, might depend on the amount of stations you plot.
        filename : str or Path, optional (default: None)
            Filename where image is stored. If None is passed, no file is created.
        ax : plt.axes
            Axes object that can be used by cartopy (projection assigned).

        Returns
        -------
        fig: matplotlib.Figure
            created figure instance. If axes was given this will be None.
        ax: matplitlib.Axes
            used axes instance, can be added to another figure for example.
        count : dict
            Number of valid sensors and stations that contain at least one valid
            sensor and networks that contain at least one valid station.
        """

        if filename and ax:
            raise ValueError(
                "Either pass a filename OR pass ax to use for plot, not both."
            )

        if not plotlibs:
            warnings.warn(
                "Could not import all plotting libs, plotting functions not available."
                "Please install cartopy and matplotlib."
            )
            return

        data_crs = ccrs.PlateCarree()

        if ax is None:
            fig, ax = plt.subplots(1, 1)
            ax = plt.axes(projection=ccrs.Robinson())
        else:
            fig = None

        ax.coastlines(linewidth=0.5)
        # show global map
        ax.set_global()
        ax.add_feature(cfeature.BORDERS, linewidth=0.5, edgecolor="gray")
        if not (sys.version_info[0] == 3 and sys.version_info[1] == 4):
            ax.add_feature(cfeature.STATES, linewidth=0.5, edgecolor="gray")
            colormap = plt.get_cmap("tab20")
        else:
            colormap = plt.get_cmap("Set1")

        all_networks = list(self.networks.keys())
        colorsteps = np.arange(0, 1, 1 / float(len(all_networks)))

        rect = []
        act_networks = []
        act_stations = []

        iterator = self.collection.iter_sensors(
            variable=variable,
            depth=Depth(min_depth, max_depth),
            filter_meta_dict=None,
            check_only_sensor_depth_from=check_only_sensor_depth_from,
        )

        n_sens = 0
        for nw, stat, sens in iterator:
            netcolor = colormap(colorsteps[all_networks.index(nw.name)])
            if nw.name not in act_networks:
                act_networks.append(nw.name)
                rect.append(Rectangle((0, 0), 1, 1, fc=netcolor))

            if stat.name not in act_stations:
                act_stations.append(stat.name)
                ax.plot(
                    stat.lon,
                    stat.lat,
                    color=netcolor,
                    markersize=markersize,
                    marker="s",
                    transform=data_crs,
                )
            n_sens += 1

        nrows = 8.0 if len(act_networks) > 8 else len(act_networks)

        try:
            ncols = int(len(act_networks) / nrows)
        except ZeroDivisionError:
            ncols = 0
        if ncols == 0:
            ncols = 1

        handles, labels = ax.get_legend_handles_labels()
        lgd = ax.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, -0.1))

        ax.legend(
            rect,
            act_networks,
            loc="upper center",
            ncol=ncols,
            bbox_to_anchor=(0.5, -0.05),
            fontsize=4,
        )

        postfix_depth = (
            "when only considering depth_from of the sensor"
            if check_only_sensor_depth_from
            else ""
        )
        depth_text = f"between {min_depth} and {max_depth} m \n {postfix_depth}"
        feedback = (
            f"{n_sens} valid sensors in {len(act_stations)} stations "
            f"in {len(act_networks)} networks (of {len(all_networks)} potential networks) \n"
            f"for {f'variable {variable}' if variable is not None else 'all variables'} "
            f"{depth_text}"
        )

        if stats_text:
            text = ax.text(
                0.5,
                1.05,
                feedback,
                transform=ax.transAxes,
                fontsize="xx-small",
                horizontalalignment="center",
            )
        else:
            text = None

        if fig:
            fig.set_size_inches([6, 3.5 + 0.25 * nrows])

        if filename is not None:
            fig.savefig(
                filename,
                bbox_extra_artists=(lgd, text) if stats_text else (lgd),
                dpi=300,
            )
        else:
            counts = (len(act_networks), len(act_stations), n_sens)
            return fig, ax, counts

    def get_min_max_obs_timestamps(
        self,
        variable="soil_moisture",
        min_depth=-np.inf,
        max_depth=np.inf,
        filter_meta_dict=None,
    ):
        """
        Filter the active file list and return the min/max time stamp from ALL
        time series that match the passed criteria.
        This time period does NOT apply to all time series in the collection
        but is the OVERALL earliest and latest timestamp found.

        Parameters
        ----------
        variable : str, optional (default: 'soil_moisture')
            One of those in :const:`ismn.const.VARIABLE_LUT` or returned by
            :func:`ismn.interface.ISMN_Interface.get_variables`:
            'soil_moisture', 'soil_temperature', 'soil_suction',
            'precipitation', 'air_temperature', 'field_capacity',
            'permanent_wilting_point', 'plant_available_water',
            'potential_plant_available_water', 'saturation', 'silt_fraction',
            'snow_depth', 'sand_fraction', 'clay_fraction', 'organic_carbon',
            'snow_water_equivalent', 'surface_temperature',
            'surface_temperature_quality_flag_original'
        min_depth : float, optional (default: -np.inf)
            Only sensors in this depth are considered.
        max_depth : float, optional (default: np.inf)
            Only sensors in this depth are considered.
        filter_meta_dict: dict, optional (default: None)
            Additional metadata keys and values for which the file list is filtered
            e.g. {'lc_2010': 10} to filter for a landcover class.
            if there are multiple conditions, ALL have to be fulfilled.
            e.g. {'lc_2010': 10', 'climate_KG': 'Dfc'})

        Returns
        -------
        start_date: datetime.datetime
            Earliest time stamp found in all sensors that fulfill the passed
            requirements.
        end_date: datetime.datetime
            Latest time stamp found in all sensors that fulfill the passed
            requirements.
        """
        t_min = None
        t_max = None

        for net, stat, sens in self.collection.iter_sensors(
            variable=variable,
            depth=Depth(min_depth, max_depth),
            filter_meta_dict=filter_meta_dict,
        ):

            time_from = pd.Timestamp(sens.metadata["timerange_from"].val)
            time_to = pd.Timestamp(sens.metadata["timerange_to"].val)

            if t_min is None:
                t_min = time_from
            if t_max is None:
                t_max = time_to
            if time_from < t_min:
                t_min = time_from
            if time_to > t_max:
                t_max = time_to

        t_min = t_min.to_pydatetime() if t_min is not None else None
        t_max = t_max.to_pydatetime() if t_max is not None else None

        return t_min, t_max

    def get_static_var_vals(
        self,
        variable="soil_moisture",
        min_depth=-np.inf,
        max_depth=np.inf,
        static_var_name="lc_2010",
    ) -> dict:
        """
        Get unique meta values for the selected static variable in the active
        networks.

        Parameters
        ----------
        variable : str, optional (default: 'soil_moisture')
            One of those in :const:`ismn.const.VARIABLE_LUT` or returned by
            :func:`ismn.interface.ISMN_Interface.get_variables`:
            'soil_moisture', 'soil_temperature', 'soil_suction',
            'precipitation', 'air_temperature', 'field_capacity',
            'permanent_wilting_point', 'plant_available_water',
            'potential_plant_available_water', 'saturation', 'silt_fraction',
            'snow_depth', 'sand_fraction', 'clay_fraction', 'organic_carbon',
            'snow_water_equivalent', 'surface_temperature',
            'surface_temperature_quality_flag_original'
        min_depth : float, optional (default: -np.inf)
            Only sensors in this depth are considered.
        max_depth : float, optional (default: np.inf)
            Only sensors in this depth are considered.
        static_var_name : str, optional (default: 'lc_2010')
            One of:
            'lc_2000', 'lc_2005', 'lc_2010', 'lc_insitu', 'climate_KG',
            'climate_insitu'

        Returns
        -------
        vals : dict
            Unique values found in static meta and their meanings.
        """

        if static_var_name not in CSV_META_TEMPLATE_SURF_VAR.keys():
            raise ValueError(
                f"{static_var_name} is not in the list of supported variables."
                f"Choose one of {list(CSV_META_TEMPLATE_SURF_VAR.keys())}"
            )

        vals = []
        for net in self.networks.values():
            for sta in net.stations.values():
                for sen in sta.sensors.values():
                    if sen.eval(variable=variable, depth=Depth(min_depth, max_depth)):
                        vals.append(sen.filehandler.metadata[static_var_name].val)

        val_dict = {}
        for val in np.unique(np.array(vals)):
            if val in self.climate.keys():
                val_dict[val] = self.climate[val]
            elif val in self.landcover.values():
                for k, v in self.landcover.items():
                    if v == val:
                        val_dict[v] = k

        return val_dict

    def get_landcover_types(
        self,
        variable: str = "soil_moisture",
        min_depth: float = 0,
        max_depth: float = 10,
        landcover: str = "lc_2010",
    ) -> dict:
        """
        See :func:`ismn.interface.ISMN_Interface.get_static_var_vals`
        """
        return self.get_static_var_vals(variable, min_depth, max_depth, landcover)

    def get_climate_types(
        self,
        variable: str = "soil_moisture",
        min_depth: float = 0,
        max_depth: float = 10,
        climate: str = "climate_KG",
    ) -> dict:
        """
        See :func:`ismn.interface.ISMN_Interface.get_static_var_vals`
        """
        return self.get_static_var_vals(variable, min_depth, max_depth, climate)

    def get_variables(self) -> np.ndarray:
        """
        get a list of variables available in the data
        """
        all_vars = np.array([])
        for _, station in self.collection.iter_stations():
            stat_vars = station.get_variables()
            if not all(np.isin(stat_vars, all_vars)):
                all_vars = np.union1d(stat_vars, all_vars)

        return all_vars

    def print_landcover_dict(self) -> None:
        """
        print all classes provided by the CCI Landcover Classification
        """
        print("CCI Landcover Classification")
        print("----------------------------")
        for key in self.landcover.keys():
            print("{:4}: {}".format(key, self.landcover[key]))

    def print_climate_dict(self) -> None:
        """
        print all classes provided by the Koeppen-Geiger climate Classification
        """
        print("KOEPPEN GEIGER Climate Classification")
        print("-------------------------------------")
        for key in self.climate.keys():
            print("{:4}: {}".format(key, self.climate[key]))

    def close_files(self):
        # close all open filehandlers
        self.__file_collection.close()
