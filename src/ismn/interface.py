# -*- coding: utf-8 -*-

import platform
import warnings

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
    warnings.warn("Could not import all plotting libs, plotting functions not available.")
    plotlibs = False

from pathlib import Path
from ismn.collection import *
from tempfile import gettempdir
from ismn.components import *

import sys
import os
import numpy as np

from ismn import tables

class ISMN_Interface():
    """
    class provides interface to ISMN data downloaded from the ISMN website
    upon initialization it collects metadata from all files in
    path_to_data and saves metadata information
    in numpy file in folder path_to_data/python_metadata/
    First initialization can take a minute or so if all ISMN
    data is present in path_to_data

    Parameters
    ----------
    data_path : str or Path
        Path to ISMN data to read, either to a zip archive or to the extracted
        directory.
    meta_path : str or Path
        Path to a metadata.pkl file that contains the filelist of a previously
        loaded IsmnFileCollection. If no meta_path is given, the metadata.pkl
        file is searched in the data root dir in the python_metadata folder.
        If it is not found, a new collection will be created and stored.
        NOTE: Paths to data in the metadata must match to the here passed data
        path.
    load_all : bool, optional (default: False)
        Load each ISMN time series directly upon initialisation. Can conusme much
        memory but leads to faster reading afterwards.
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
                 load_all=False, temp_root=gettempdir()):

        self.climate, self.landcover = tables.KOEPPENGEIGER, tables.LANDCOVER

        __root = IsmnRoot(data_path)

        if meta_path is not None:
            self.meta_path = Path(meta_path)
        else:
            self.meta_path = __root.root_dir / 'python_metadata'

        if not os.path.exists(self.meta_path):
            os.makedirs(self.meta_path)

        pklpath = self.meta_path / f'{__root.name}.pkl'

        if os.path.isfile(pklpath):
            logger.info(f'Loading existing filelist: {pklpath}')
            self.__file_collection = IsmnFileCollection.from_pkl(
                pklpath, temp_root=temp_root)
        else:
            logger.info(f"Create new filelist into: {pklpath}")
            self.__file_collection = IsmnFileCollection(__root.path, temp_root=temp_root)
            self.__file_collection.store(pklpath, drop_filehandler=False)
            # todo: write the logfile

        self.activate_network(network, load_all=load_all)

    def activate_network(self, network=None, load_all=False):
        """
        Takes either all networks from the full collection or a selection
        if network names to create a network collection for them.

        Parameters
        ----------
        network : list or str, optional (default: None)
            Name of the network(s) to activate
        load_all : bool, optional (default: False)
            Load data for all networks immediately and keep it in memory
            for faster access.
        """

        if network is None:
            self.act_coll = NetworkCollection(self.__file_collection,
                                              load_all=load_all)
        else:
            active_files = IsmnFileCollection.from_filelist(
                self.__file_collection.filter_col_val('network', network),
                temp_root=self.__file_collection.temp_root)
            self.act_coll = NetworkCollection(active_files, load_all=load_all)

    @property
    def grid(self):
        return self.act_coll.grid

    @property
    def networks(self):
        return self.act_coll.networks

    @property
    def files(self):
        return self.act_coll.files

    def __repr__(self):
        """
        Print summary of ISMN file collection.
        """
        logger.info('Number of networks: {}'.format(len(self.list_networks())))
        logger.info('Number of stations: {}'.format(len(self.list_stations())))
        logger.info('Number of sesors: {}'.format(len(self.list_sensors())))

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
        network_names : list[str]
            List of networks that contain the station.
            Usually one, but can be more.
        """
        network_with_station = []

        for network in self.networks.values():
            if stationname in network.stations.keys():
                network_with_station.append(network.name)

        # if len(network_with_station) == 0:
        #     raise ISMNError("stationname was not found")
        # if len(network_with_station) > 1:
        #     raise ISMNError("stationname occurs in multiple networks")

        # return self.networks[network_with_station]
        return network_with_station

    def stations_that_measure(self, variable):
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

        Returns
        -------
        ISMN_station : Station
        """
        for network in self.networks: # type: Network
            for station in network.iter_stations(variable, depth=None):
                yield station

    def get_dataset_ids(self, variable, min_depth=0, max_depth=0.1,
                        filter_static_vars=None):
        """
        Filter the filelist for active networks for variables, depths and
        metadata values to get the row numbers in the filelist that comply.

        Parameters
        ----------
        variable : str
            Variable to filer out
        min_depth : float, optional (default: 0)
            Min depth of sensors to search
        max_depth : float, optional (default: 0.1)
            Max depth of sensors to search
        filter_static_vars: dict, optional (default: None)
            Additional metadata keys and values for which the file list is filtered
            e.g. {'lc_2010': 10} to filter for a landcover class.
            if there are multiple conditions, ALL have to be fulfilled.
            e.g. {'lc_2010': 10', 'climate_KG': 'Dfc'})
        """
        return self.act_coll.get_dataset_ids(variable, min_depth=min_depth,
            max_depth=max_depth, filter_static_vars=filter_static_vars)

    def read_ts(self, idx):# todo: load data?
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
        self.act_coll.files[idx]['filehandler'].read_data()

    def read(self, *args, **kwargs):
        # calls read_ts
        return self.read_ts(*args, **kwargs)

    def find_nearest_station(self, lon, lat, return_distance=False, max_dist=np.inf):
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
        idx, d = self.act_coll.grid.find_nearest_gpi(lon, lat, max_dist=max_dist)

        if idx is None: # todo: not sure what this looks like when pygeogrids is fixed.
            stat = None
        else:
            net_name, stat_name = self.act_coll.files.loc[idx, ['network', 'station']]
            stat = self.networks[net_name].stations[stat_name]

        if return_distance:
            return stat, d
        else:
            return stat

    def plot_station_locations(self, variable=None, min_depth=-np.inf,
                               max_depth=np.inf, stats_text=True,
                               check_only_sensor_depth_from=False, filename=None):
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
        uniq_networks = self.list_networks()
        colorsteps = np.arange(0, 1, 1 / float(uniq_networks.size))
        rect = []

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
                    ax.plot(station.lon, station.lat, color=netcolor, markersize=1,
                            marker='s', transform=data_crs)
            if station_count > 0:
                counts['networks'] += 1
                counts['stations'] += station_count

            rect.append(Rectangle((0, 0), 1, 1, fc=netcolor))

        nrows = 8. if len(uniq_networks) > 8 else len(uniq_networks)

        ncols = int(counts['networks'] / nrows)
        if ncols == 0:
            ncols = 1

        handles, labels = ax.get_legend_handles_labels()
        lgd = ax.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, -0.1))

        plt.legend(rect, uniq_networks.tolist(), loc='upper center', ncol=ncols,
                   bbox_to_anchor=(0.5, -0.05), fontsize=4)

        postfix_depth = "when only considering depth_from of the sensor" if check_only_sensor_depth_from else ''
        depth_text =  f"between {min_depth} and {max_depth} m \n {postfix_depth}"
        feedback = f"{counts['sensors']} valid sensors in {counts['stations']} stations " \
                   f"in {counts['networks']} networks (of {len(uniq_networks)} potential networks) \n" \
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

        return fig, ax, counts


    def get_min_max_obs_timestamps(self, variable="soil moisture",
                                   min_depth=None, max_depth=None):
        # todo: implement?
        raise NotImplementedError

    def get_static_var_vals(self, variable='soil_moisture', min_depth=0, max_depth=10,
                            static_var_name='lc_2010'):
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
            'lc_2000', 'lc_2005', 'lc_2010', 'lc_insitu', 'climate_KG', 'climate_insitu'

        Returns
        -------
        vals : dict
            Unique values found in static meta and their meanings.
        """

        if static_var_name not in tables.CSV_META_TEMPLATE_SURF_VAR.keys():
            raise ValueError(f"{static_var_name} is not in the list of supported variables."
                             f"Choose one of {list(CSV_META_TEMPLATE_SURF_VAR.keys())}")

        vals = []
        for net in self.networks.values():
            for sta in net.stations.values():
                for sen in sta.sensors.values():
                    if sen.eval(variable=variable, depth=Depth(min_depth, max_depth)):
                        vals.append(sen.filehandler.metadata[static_var_name].val)

        val_dict = None
        for d in [self.landcover, self.climate]:
            if static_var_name in d.keys():
                val_dict = {k : d[k] for k in d.keys() if d[k] in np.unique(np.array(vals))}

        return val_dict

    def get_landcover_types(self, variable:str='soil_moisture', min_depth:float=0,
                          max_depth:float=10, landcover:str='lc_2010') -> dict:
        """
        See description of get_static_var_vals
        """
        return self.get_static_var_vals(variable, min_depth, max_depth, landcover)

    def get_climate_types(self, variable:str='soil_moisture', min_depth:float=0,
                          max_depth:float=10, climate:str='climate') -> dict:
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
        return np.unique(self.files['variable'].values)


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
    ds = ISMN_Interface(r"C:\Temp\delete_me\ismn\testdata_ceop.zip")
    fig, ax, count = ds.plot_station_locations(variable='soil_moisture', min_depth=0, max_depth=0.05,
                                               stats_text=False, filename=r"C:\Temp\delete_me\map.png",
                                               check_only_sensor_depth_from=False)
    # ds.find_nearest_station(1,1)
    # ds.network_for_station('SilverSword')
    # ids = ds.get_dataset_ids('soil_moisture')

