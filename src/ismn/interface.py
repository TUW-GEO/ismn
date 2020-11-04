# -*- coding: utf-8 -*-


import platform

from ismn.network_collection import NetworkCollection

if platform.system() == 'Darwin':
      import matplotlib
      matplotlib.use("TkAgg")

from pathlib import Path
from ismn.file_collection import *
from tempfile import gettempdir
from ismn.components import *

import os
import numpy as np
import configparser


def load_config():
    config = configparser.ConfigParser()
    config.optionxform = str
    config.read_file(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                       'classifications.ini')))
    landcover = dict(config.items('LANDCOVER'))
    landcover = dict([(int(v), k) for k, v in landcover.items()])
    climate = dict(config.items('KOEPPENGEIGER'))

    return landcover, climate

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

        self.climate, self.landcover = load_config()

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

    def plot_station_locations(self):

        pass

    def get_min_max_obs_timestamps(self, variable="soil moisture",
                                   min_depth=None, max_depth=None):
        pass

    def get_landcover_types(self, variable='soil moisture', min_depth=0, max_depth=10,
                            landcover='landcover_2010'):
        pass

    def get_climate_types(self, variable='soil moisture', min_depth=0, max_depth=10,
                          climate='climate'):
        pass

    def get_variables(self):
        pass


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
    ds = ISMN_Interface(r"H:\code\qa4sm-testdata\input_data\ISMN\ISMN_TESTDATA_HAWAII")
    ds.find_nearest_station(1,1)
    ds.network_for_station('SilverSword')
    ids = ds.get_dataset_ids('soil_moisture')

