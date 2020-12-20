# -*- coding: utf-8 -*-

from pathlib import Path
from tempfile import gettempdir
import platform
import os
import sys

from ismn.filecollection import IsmnFileCollection
from ismn.components import *
from ismn.tables import *
from ismn.base import IsmnRoot

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    if platform.system() == 'Darwin':
          import matplotlib
          matplotlib.use("TkAgg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    plotlibs = True
except ImportError:
    plotlibs = False


class ISMN_Interface():
    """
    Class provides interface to ISMN data downloaded from the ISMN website
    upon initialization it collects metadata from all files in
    path_to_data and saves metadata information in a csv file into the folder
    python_metadata in meta_path (or data_path if no meta_path is defined).
    First initialization can take some time if all ISMN
    data is present in data_path and will start multiple processes.

    Parameters
    ----------
    data_path : str or Path
        Path to ISMN data to read, either to a zip archive or to the extracted
        directory.
    meta_path : str or Path
        Path where the metadata csv file(s) is / are stored. The actual filename
        is defined by the name of data_path and will be generated automatically.
    keep_loaded_data : bool, optional (default: False)
        Keep data for a file in memory once it is loaded. This makes subsequent
        calls of data faster (if e.g. a station is accessed multiple times)
        but can fill up memory if multiple networks are loaded.
    network : str or list, optional (default: None)
        Name(s) of network(s) to load. Other data in the data_path will be ignored.
        By default all networks are activated.

    Raises
    ------
    ISMNError
        if given network was not found in path_to_data

    Attributes
    ----------
    metadata : numpy.array
        metadata array for all stations contained in the path given during initialization
    grid : pygeogrids.grid.BasicGrid
        Grid object used for finding nearest insitu station for given lon lat
    Methods
    -------
    find_nearest_station(lon,lat)
        find nearest station for given coordinates
    """

    def __init__(self, data_path, meta_path=None, network=None,
                 keep_loaded_data=False, temp_root=gettempdir()):

        self.climate, self.landcover = KOEPPENGEIGER, LANDCOVER

        self.root = IsmnRoot(data_path)

        meta_csv_filename = f'{self.root.name}.csv'

        if meta_path is None:
            meta_path = Path(self.root.root_dir) / 'python_metadata'
        else:
            meta_path = Path(meta_path)

        meta_csv_file = meta_path / meta_csv_filename

        if os.path.isfile(meta_csv_file):
            self.__file_collection = IsmnFileCollection.from_metadata_csv(
                self.root, meta_csv_file)
        else:
            self.__file_collection = IsmnFileCollection.build_from_scratch(
                self.root, parallel=True, log_path=meta_path, temp_root=temp_root)
            self.__file_collection.to_metadata_csv(meta_csv_file)

        self.keep_loaded_data = keep_loaded_data

        networks = self.__collect_networks(network)
        self.collection = NetworkCollection(networks)

    def __collect_networks(self, networks:list=None) -> list:
        """
        Build networks and fill them with stations and sensors and apply
        filehandlers for data reading.
        """
        if networks is not None:
            filelist = self.__file_collection.filter_networks('network', networks)
        else:
            filelist = self.__file_collection.files

        networks = OrderedDict([])
        # points = []  # idx, lon, lat

        for idx, row in filelist.iterrows(): # todo: slow iterrows??
            f = row['filehandler']

            nw_name, st_name, instrument = f.metadata['network'].val, \
                                           f.metadata['station'].val, \
                                           f.metadata['instrument'].val

            if nw_name not in networks:
                networks[nw_name] = Network(nw_name)

            if st_name not in networks[nw_name].stations:
                networks[nw_name].add_station(st_name,
                                              f.metadata['longitude'].val,
                                              f.metadata['latitude'].val,
                                              f.metadata['elevation'].val)

            # the senor name is the index in the list
            if idx not in networks[nw_name].stations[st_name].sensors:
                networks[nw_name].stations[st_name]. \
                    add_sensor(instrument,
                               f.metadata['variable'].val,
                               f.metadata['variable'].depth,
                               filehandler=f, # todo: remove station meta from sensor
                               name=idx,
                               keep_loaded_data=self.keep_loaded_data)
                #points.append((idx, f.metadata['longitude'].val, f.metadata['latitude'].val))

        #points = np.array(points)

        #grid = BasicGrid(points[:, 1], points[:, 2], gpis=points[:,0])

        return list(networks.values()) # , grid

    def __repr__(self):
        indent = 2
        return f"root:\n" \
               f"{' '* indent}{self.root}\n" + \
               f"-" * (len(str(self.root)) + indent) + "\n" \
               f"networks:\n{self.collection.__repr__(indent)}"


    @property
    def networks(self):
        return self.collection.networks

    @property
    def grid(self):
        return self.networks.grid

    def load_all(self):
        """
        Load data for all file handlers in the current network collection.
        This may take some time and fill up your memory if multiple networks are
        loaded at once.
        """
        if not self.keep_loaded_data:
            raise IOError("Can only load all data when storing to memory is allowed. "
                          "Pass keep_loaded_data=True when creating the NetworkCollection.")

        for net in self.networks.iter_networks():
            for stat in net.iter_stations():
                print(stat.name)
                for sens in stat.iter_sensors():
                    assert sens.keep_loaded_data == True
                    sens.read_data()

    def list_networks(self) -> np.array:
        # get network names from list of active files
        return np.array(list(self.networks.keys()))

    def list_stations(self, network=None) -> np.array:
        # get station names for one of the active networks
        if network is not None:
            if network not in self.networks:
                raise ISMNError(f'Network {network} not found in currently loaded networks.')
            return np.array(list(self.networks[network].stations.keys()))
        else:
            stations = []
            for network in self.networks.values():
                stations += list(network.stations.keys())
            return np.array(stations)

    def list_sensors(self, network:str=None, station:str=None) -> np.array:
        # List sensors names for a specific sensor in an active network
        sensors = np.array([])
        for net in self.networks.values():
            if network in [None, net.name]:
                for stat in net.stations.values():
                    if station in [None, stat.name]:
                        sensors = np.append(sensors,
                                            # get the objects instead, use .values()?
                                            np.array(list(stat.sensors.keys())))

        return sensors

    def network_for_station(self, stationname):
        """
        Find networks that contain a station of the passed name.

        Parameters
        ----------
        stationname : str
            Station name to search in the active networks.

        Returns
        -------
        network_names : str or None
            Network that contains a station of that name, or None if no such
            network exists.
            Prints are warning and uses the FIRST found network name if there
            are multiple stations with the same name in different networks.
        """
        network_with_station = []

        for network in self.networks.values():
            if stationname in network.stations.keys():
                network_with_station.append(network.name)

        if len(network_with_station) > 1:
           warnings.warn("stationname occurs in multiple networks")

        if len(network_with_station) == 0:
            return None
        else:
            return network_with_station[0]

    def stations_that_measure(self, variable, **eval_kwargs):
        """
        Goes through all stations and returns those that measure the specified
        variable

        Parameters
        ----------
        variable : string
            variable name
            one of
                * 'soil moisture',
                * 'soil temperature',
                * 'soil suction',
                * 'precipitation',
                * 'air temperature',
                * 'field capacity',
                * 'permanent wilting point',
                * 'plant available water',
                * 'potential plant available water',
                * 'saturation',
                * 'silt fraction',
                * 'snow depth',
                * 'sand fraction',
                * 'clay fraction',
                * 'organic carbon',
                * 'snow water equivalent',
                * 'surface temperature',
                * 'surface temperature quality flag original'

        Yields
        -------
        ISMN_station : Station
        """
        for network in self.networks.values(): # type: Network
            for station in network.iter_stations(variable, **eval_kwargs):
                yield station

    def get_dataset_ids(self, variable, min_depth=0, max_depth=0.1,
                        filter_meta_dict=None, check_only_sensor_depth_from=False):
        """
        Yield all sensors for a specific network and/or station and/or
        variable and/or depth.

        Parameters
        ----------
        variable : str
            Variable to filer out
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
        """
        ids = []

        depth = Depth(min_depth, max_depth)
        for _, _, sensor in self.collection.iter_sensors(
                variable=variable, depth=depth,
                filter_meta_dict=filter_meta_dict,
                check_only_sensor_depth_from=check_only_sensor_depth_from):
            ids.append(sensor.name)

        return ids

    def read_ts(self, idx):# todo: load data?
        # todo: rename data column from variable to e.g. soil_moisture?
        """
        Read a time series directly by the id

        Parameters
        ----------
        idx : int
            id into self.metadata, best one of those returned
            from get_dataset_ids()
        Returns
        -------
        timeseries : pandas.DataFrame
            of the read data
        """
        station = self.collection.station4idx(idx)
        for se in station.iter_sensors():
            if se.name == idx:
                return se.read_data()

    def read(self, *args, **kwargs):
        # calls read_ts
        return self.read_ts(*args, **kwargs)

    def find_nearest_station(self, lon, lat, return_distance=False,
                             max_dist=np.inf):
        """
        finds the nearest station available in downloaded data
        Parameters
        ----------
        lon : float
            Longitude of point
        lat : float
            Latitude of point
        return_distance : boolean, optional
            if True also distance is returned
        max_dist : float, optional (default: np.inf)
            Maximum distance allowed. # todo: fix in pygeogrids
            
        Returns
        -------
        station : ISMN_station
            ISMN_station object
        distance : float, optional
            distance to station in meters, measured in cartesian coordinates and not on
            a great circle. Should be OK for small distances
        """
        # todo: fix bug in pygeogrids that leads to always np.inf as max dist
        # what happens if there is no point within max dist if that works?
        idx, d = self.collection.grid.find_nearest_gpi(lon, lat, max_dist=max_dist)

        if idx is None: # todo: not sure what this looks like when pygeogrids is fixed.
            stat = None
        else:
            stat = self.collection.station4idx(idx)

        if return_distance:
            return stat, d
        else:
            return stat

    def plot_station_locations(self, variable=None, min_depth=-np.inf,
                               max_depth=np.inf, stats_text=True,
                               check_only_sensor_depth_from=False,
                               markersize=1, filename=None):
        # TODO: optionally indicate sensor count for each station in map directly (symbol, number)?
        # TODO: fix similar colors for different networks, e.g. using sybols?
        """
        Plots available stations on a world map in robinson projection.

        Parameters
        ----------
        variable : str, optional (default: None)
            Show only stations that measure this variable, e.g. soil_moisture
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

        Returns
        -------
        fig: matplotlib.Figure
            created figure instance. If axes was given this will be None.
        ax: matplitlib.Axes
            used axes instance.
        count : dict
            Number of valid sensors and stations that contain at least one valid
            sensor and networks that contain at least one valid station.
        """

        if not plotlibs:
            warnings.warn("Could not import all plotting libs, plotting functions not available.")
            return

        data_crs = ccrs.PlateCarree()

        fig, ax = plt.subplots(1, 1)
        ax = plt.axes(projection=ccrs.Robinson())
        ax.coastlines(linewidth=0.5)
        # show global map
        ax.set_global()
        ax.add_feature(cfeature.BORDERS, linewidth=0.5, edgecolor='gray')
        if not (sys.version_info[0] == 3 and sys.version_info[1] == 4):
            ax.add_feature(cfeature.STATES, linewidth=0.5, edgecolor='gray')
            colormap = plt.get_cmap('tab20')
        else:
            colormap = plt.get_cmap('Set1')
        uniq_networks = list(self.networks.keys())
        colorsteps = np.arange(0, 1, 1 / float(len(uniq_networks)))
        rect = []

        uniq_networks = []
        counts = {'networks': 0, 'stations': 0, 'sensors': 0}

        for j, (nw_name, nw) in enumerate(self.networks.items()):
            netcolor = colormap(colorsteps[j])
            station_count = 0
            for station in nw.stations.values():
                sensor_count = 0 # count number of valid sensors at station, indicate number in map?
                for sensor in station.sensors.values():
                    # could add filtering for other metadata here,
                    # e.g. for landcover using the filter_meta_dict .. slow
                    if sensor.eval(variable, depth=Depth(min_depth, max_depth),
                                   filter_meta_dict=None,
                                   check_only_sensor_depth_from=check_only_sensor_depth_from):
                        sensor_count += 1
                if sensor_count > 0:
                    station_count += 1
                    counts['sensors'] += sensor_count
                    ax.plot(station.lon, station.lat, color=netcolor, markersize=markersize,
                            marker='s', transform=data_crs)
            if station_count > 0:
                counts['networks'] += 1
                counts['stations'] += station_count

                uniq_networks.append(nw_name)
                rect.append(Rectangle((0, 0), 1, 1, fc=netcolor))

        nrows = 8. if len(uniq_networks) > 8 else len(uniq_networks)

        ncols = int(counts['networks'] / nrows)
        if ncols == 0:
            ncols = 1

        handles, labels = ax.get_legend_handles_labels()
        lgd = ax.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, -0.1))

        plt.legend(rect, uniq_networks, loc='upper center', ncol=ncols,
                   bbox_to_anchor=(0.5, -0.05), fontsize=4)

        postfix_depth = "when only considering depth_from of the sensor" if check_only_sensor_depth_from else ''
        depth_text =  f"between {min_depth} and {max_depth} m \n {postfix_depth}"
        feedback = f"{counts['sensors']} valid sensors in {counts['stations']} stations " \
                   f"in {counts['networks']} networks (of {len(list(self.networks.keys()))} potential networks) \n" \
                   f"for {f'variable {variable}' if variable is not None else 'all variables'} " \
                   f"{depth_text}"

        if stats_text:
            text = ax.text(0.5, 1.05, feedback, transform=ax.transAxes, fontsize='xx-small',
                           horizontalalignment='center')
        else:
            text = None

        fig.set_size_inches([6, 3.5 + 0.25 * nrows])

        if filename is not None:
            fig.savefig(filename, bbox_extra_artists=(lgd, text) if stats_text else (lgd),
                        dpi=300)
        else:
            return fig, ax, counts

    def get_min_max_obs_timestamps(self, variable="soil_moisture",
                                   min_depth=-np.inf, max_depth=np.inf,
                                   filter_meta_dict=None):
        """
        Filter the active file list and return the min/max time stamp from ALL
        time series that match the passed criteria.
        This time period does NOT apply to all time series in the collection
        but is the OVERALL earliest and latest timestamp found.

        Parameters
        ----------
        variable : str, optional (default: 'soil_moisture')
            One of:
            'soil_moisture', 'soil_temperature', 'soil_suction',
            'precipitation', 'air_temperature', 'field_capacity',
            'permanent_wilting_point', 'plant_available_water',
            'potential_plant_available_water', 'saturation', 'silt_fraction',
            'snow_depth', 'sand_fraction', 'clay_fraction', 'organic_carbon',
            'snow_water_equivalent', 'surface_temperature',
            'surface_temperature_quality_flag_original'
        min_depth : float, optional (default: 0)
            Only sensors in this depth are considered.
        max_depth : float, optional (default: 10)
            Only sensors in this depth are considered.
        filter_meta_dict: dict, optional (default: None)
            Additional metadata keys and values for which the file list is filtered
            e.g. {'lc_2010': 10} to filter for a landcover class.
            if there are multiple conditions, ALL have to be fulfilled.
            e.g. {'lc_2010': 10', 'climate_KG': 'Dfc'})

        Returns
        -------
        start_date: datetime
            Earliest time stamp found in all sensors that fulfill the passed
            requirements.
        end_date: datetime
            Latest time stamp found in all sensors that fulfill the passed
            requirements.
        """
        t_min = None
        t_max = None

        for net, stat, sens in self.collection.iter_sensors(
            variable=variable, depth=Depth(min_depth, max_depth),
                filter_meta_dict=filter_meta_dict):

            time_from = pd.Timestamp(sens.metadata['timerange_from'].val)
            time_to = pd.Timestamp(sens.metadata['timerange_to'].val)

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

    def get_static_var_vals(self, variable='soil_moisture', min_depth=0,
                            max_depth=10, static_var_name='lc_2010'):
        """
        Get unique meta values for the selected static variable in the active
        networks.

        Parameters
        ----------
        variable : str, optional (default: 'soil_moisture')
            One of:
            'soil_moisture', 'soil_temperature', 'soil_suction',
            'precipitation', 'air_temperature', 'field_capacity',
            'permanent_wilting_point', 'plant_available_water',
            'potential_plant_available_water', 'saturation', 'silt_fraction',
            'snow_depth', 'sand_fraction', 'clay_fraction', 'organic_carbon',
            'snow_water_equivalent', 'surface_temperature',
            'surface_temperature_quality_flag_original'
        min_depth : float, optional (default: 0)
            Only sensors in this depth are considered.
        max_depth : float, optional (default: 10)
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
            raise ValueError(f"{static_var_name} is not in the list of supported variables."
                             f"Choose one of {list(CSV_META_TEMPLATE_SURF_VAR.keys())}")

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

    def get_landcover_types(self, variable:str='soil_moisture', min_depth:float=0,
                          max_depth:float=10, landcover:str='lc_2010') -> dict:
        """
        See description of get_static_var_vals
        """
        return self.get_static_var_vals(variable, min_depth, max_depth, landcover)

    def get_climate_types(self, variable:str='soil_moisture', min_depth:float=0,
                          max_depth:float=10, climate:str='climate_KG') -> dict:
        """
        See description of get_static_var_vals
        """
        return self.get_static_var_vals(variable, min_depth, max_depth, climate)

    def get_variables(self):
        """
         get a list of variables available for the data

         Returns
         -------
         variables : numpy.array
             array of variables available for the data
        """
        all_vars = np.array([])
        for _, station in self.collection.iter_stations():
            stat_vars = station.get_variables()
            if not all(np.isin(stat_vars, all_vars)):
                all_vars = np.union1d(stat_vars, all_vars)

        return all_vars


    def print_landcover_dict(self):
        """
        print all classes provided by the CCI Landcover Classification
        :return: None
        """
        print('CCI Landcover Classification')
        print('----------------------------')
        for key in self.landcover.keys():
            print('{:4}: {}'.format(key, self.landcover[key]))

    def print_climate_dict(self):
        """
        print all classes provided by the Koeppen-Geiger climate Classification
        :return: None
        """
        print('KOEPPEN GEIGER Climate Classification')
        print('-------------------------------------')
        for key in self.climate.keys():
            print('{:4}: {}'.format(key, self.climate[key]))



if __name__ == '__main__':
    path = r"D:\data-read\ISMN\global_20191024"
    ds = ISMN_Interface(path)
    # mmin, mmax = ds.get_min_max_obs_timestamps('soil_moisture')
    # ids = ds.get_dataset_ids('soil_moisture', 0, 0.05, filter_meta_dict={'lc_2010': 130})
    # # ds.plot_station_locations('soil_moisture', 0., 0.1, filename="/home/wolfgang/data-write/temp/plot.png")
    # netname = ds.network_for_station('Villevielle')
    # ts = ds.read_ts(1)
    # print(ts)
