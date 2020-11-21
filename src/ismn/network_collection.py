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

from ismn.components import *
from ismn.base import IsmnRoot
from ismn.file_collection import IsmnFileCollection
from ismn.tables import *

from tempfile import gettempdir
from pathlib import Path
import os
import pandas as pd

class NetworkCollection(object):
    """
    The NetworkCollection builds ISMN Networks and fills them with stations
    and sensors using the passed file collection.
    A grid is added that contains all stations to perform spatial searches.
    """
    def __init__(self, data_path, networks=None, meta_path=None,
                 keep_loaded_data=False, temp_root=gettempdir()):

        """
        Create network collection from previously loaded IsmnFileCollection.

        Parameters
        ----------
        filelist : pd.DataFrame
            DataFrame that contains the files to read.
            As in IsmnFileCollection.files or as returned by the filter functions.
        networks : list or str, optional (default: None)
            List of network names to activate.
        meta_path : str or Path
            Path where the metadata csv file(s) is / are stored. The actual filename
            is defined by the name of data_path and will be generated automatically.
        keep_loaded_data : bool, optional (default: False)
            Keep data for a file in memory once it is loaded. This makes subsequent
            calls of data faster (if e.g. a station is accessed multiple times)
            but can fill up memory if multiple networks are loaded.
        temp_root : str or Path, optional (default: os temp dir)
            Root path where temporary files are stored.
        """
        
        root = IsmnRoot(data_path)

        meta_csv_filename = f'{root.name}.csv'

        if meta_path is None:
            meta_csv_file = root.root_dir / 'python_metadata' / meta_csv_filename
        else:
            meta_csv_file = Path(meta_path) / meta_csv_filename

        if os.path.isfile(meta_csv_file):
            self.file_collection = IsmnFileCollection.from_metadata_csv(root, meta_csv_file)
        else:
            self.file_collection = IsmnFileCollection.from_scratch(root, temp_root)
            self.file_collection.to_metadata_csv(meta_csv_file)

        self.keep_loaded_data = keep_loaded_data

        self.networks, self.grid = self._collect_networks(networks)

    def iter_networks(self):
        # Iterate through all networks
        for network in self.networks.values():
            yield network

    def load_all(self):
        """
        Load data for all file handlers in the current network collection.
        This may take some time and fill up your memory if multiple networks are
        loaded at once.
        """
        if not self.keep_loaded_data:
            raise ValueError("Can only load all data when storing to memory is allowed. "
                             "Pass keep_loaded=True when creating the NetworkCollection.")

        for net in self.iter_networks():
            for stat in net.iter_stations():
                print(stat.name)
                for sens in stat.iter_sensors():
                    assert sens.keep_loaded_data == True
                    sens.read_data()

    def _collect_networks(self, networks:list=None) -> (dict, BasicGrid):
        # build networks and fill them with stations and sensors and apply
        # filehandlers for data reading.
        if networks is not None:
            filelist = self.file_collection.filter_col_val('network', networks)
        else:
            filelist = self.file_collection.files

        networks = {}
        points = []  # idx, lon, lat

        for idx, row in filelist.iterrows():
            f = row['filehandler']

            nw_name, st_name, instrument = f['network'].val, f['station'].val, f['instrument'].val

            if nw_name not in networks:
                networks[nw_name] = Network(nw_name)

            if st_name not in networks[nw_name].stations:
                networks[nw_name].add_station(st_name,
                                              f['longitude'].val,
                                              f['latitude'].val,
                                              f['elevation'].val)

            # the senor name is the index in the list
            if idx not in networks[nw_name].stations[st_name].sensors:
                networks[nw_name].stations[st_name]. \
                    add_sensor(instrument, f['variable'].val, f['variable'].depth,
                               f, name=idx, keep_loaded_data=self.keep_loaded_data)
                points.append((idx, f['longitude'].val, f['latitude'].val))

        points = np.array(points)

        grid = BasicGrid(points[:, 1], points[:, 2], gpis=points[:,0])

        return networks, grid

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
        ids = []

        for se in self.get_sensors(variable=variable,
                                   depth=Depth(min_depth, max_depth),
                                   filter_static_vars=filter_static_vars):
            ids.append(se.name)

        return ids

    def get_sensors(self, network=None, station=None, variable=None,
                    depth=None, filter_static_vars=None):
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
        depth : Depth, optional (default: None)
            Sensing depth.
        filter_static_vars: dict, optional (default: None)
            Additional metadata keys and values for which the file list is filtered
            e.g. {'lc_2010': 10} to filter for a landcover class.
            if there are multiple conditions, ALL have to be fulfilled.
            e.g. {'lc_2010': 10', 'climate_KG': 'Dfc'})

        Yield
        -----
        sensor : Sensor
            Sensor.
        """
        if depth is None:
            depth = Depth(-np.inf, np.inf)

        for n in self.networks.values():
            if network not in [None, n.name]:
                continue
            for st in n.stations.values():
                if station not in [None, st.name]:
                    continue
                for se in st.sensors.values():
                    if se.eval(variable, depth, filter_static_vars):
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
        gpi, dist = self.grid.find_nearest_gpi(lon, lat, max_dist=max_dist)
        station = self.station4idx(gpi)

        return station, dist

if __name__ == '__main__':
    networks = NetworkCollection(r"C:\Temp\delete_me\ismn\testdata_ceop.zip",
                                 meta_path=r"C:\Temp\delete_me\ismn\temp",
                                 keep_loaded_data=True,
                                 networks=None)


