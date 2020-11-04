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

from ismn.base import IsmnRoot
from ismn.components import *
from ismn.tables import *
from ismn.files import DataFile

import pandas as pd
from tempfile import gettempdir, TemporaryDirectory
from pathlib import Path, PurePosixPath

from pygeogrids.grids import BasicGrid

class NetworkCollection(object):

    # todo: this should be in the the ISMNInterface
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

    def load_metadata(self, path):
        meta_arrays = np.load(path, allow_pickle=True)


    def store_metadata(self, path):
        meta_arrays = []
        dtype = None
        for net, stats, in self.files.items():
            for stat, files in stats.items():
                for file in files:
                    if dtype is None:
                        dtype = file.get_formatted_metadata('struct').dtype
                    meta_arrays.append(file.get_formatted_metadata('list'))

        meta_arrays = np.array(meta_arrays, dtype)
        np.save(path, meta_arrays)

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

class IsmnFileCollection(object):

    # todo: limit the metadata to the sensor depth.

    """
    The IsmnFileCollection class loads all files in a directory and creates
    a metadata/filelist, which can be used afterwards to recreate the collection.
    And provides functions to access components (networks, stations, sensors).

    Parameters
    ----------
    path : str or Path
        Root path of ISMN files. i.e. path to the downloaded zip file
        or the extracted zip directory (faster).
    load_data : bool, optional
        If True data will be loaded during metadata reading.

    Attributes
    ----------
    data_path : str
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

    def __init__(self, data_path, load_data=False, temp_root=gettempdir()):

        if not isinstance(data_path, IsmnRoot):
            self.root = IsmnRoot(data_path)
        else:
            self.root = data_path

        self.name = self.root.name

        if not os.path.exists(temp_root):
            os.makedirs(temp_root, exist_ok=True)

        self.temp_root = temp_root

        self.files = {}

        self.filelist = self._build_filelist(load_data)

    def _build_filelist(self, load_data=False) -> pd.DataFrame:
        for net_dir, stat_dirs in self.root.cont.items():
            for stat_dir in stat_dirs:

                data_files = self.root.find_files(stat_dir, '*.stm')
                static_meta = None

                for file_path in data_files:

                    try:
                        f = DataFile(self.root, file_path, load_data,
                                     static_meta=static_meta,
                                     temp_root=self.temp_root)
                    except IOError as e:
                        logger.error(f'Error loading ismn file: {e}')
                        continue


        
    def _old_build_filelist(self, load_data=False):
        """
        Scan ismn archive and build file object list.
        Reuse static metadata if possible
        """
        # todo: can this be stored directly? Maybe if load_data is False?

        for net_dir, stat_dirs in self.root.cont.items():
            for stat_dir in stat_dirs:
                print(stat_dir)
                data_files = self.root.find_files(stat_dir, '*.stm')
                static_meta = None

                for file_path in data_files:

                    try:
                        f = DataFile(self.root, file_path, load_data,
                                     static_meta=static_meta,
                                     temp_root=self.temp_root)
                    except IOError as e:
                        logger.error(f'Error loading ismn file: {e}')
                        continue

                    static_meta = f.static_meta

                    if f['network'] not in self.files.keys():
                        self.files[f['network']] = {}
                    if f['station'] not in self.files[f['network']]:
                        self.files[f['network']][f['station']] = []

                    self.files[f['network']][f['station']].append(f)

    def get_networks(self):
        """
        Get networks from ISMN file collection.

        Returns
        -------
        networks : dict of empty networks
            Dict of networks.
        """
        networks = {}

        # todo: repr says that network has 0 stations
        for n in self.files.keys():
            networks[n] = Network(n)

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

        for net, stats in self.files.items():
            if network in [None, net]:
                for stat, files in stats.items():
                    for f in files:
                        if f['station'] not in stations:
                            stations[f['station']] = Station(f['station'],
                                                             f['longitude'],
                                                             f['latitude'],
                                                             f['elevation'])

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

        for net, stats in self.files.items():
            if network in [None, net]:
                for stat, files in stats.items():
                    if stat in [None, station]:
                        for f in files:
                            name = f"{f['sensor']}_{f['variable']}_{f['depth'].start}_{f['depth'].end}"

                            if name not in sensors:
                                snr = Sensor(name, f['variable'],
                                             f['sensor'], f['depth'])
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

if __name__ == '__main__':
    path = r"C:\Temp\delete_me\ismn\testdata_ceop"
    coll = IsmnFileCollection(path, load_data=False)



