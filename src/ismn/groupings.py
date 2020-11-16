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
from ismn.components import Network, Depth
from ismn.tables import *
from ismn.filehandlers import DataFile

import ismn
pkg_version = ismn.__version__

import pandas as pd
from tempfile import gettempdir
from pathlib import Path
import warnings

class MetadataError(IOError):
    pass

class ISMNError(Exception):
    pass

def build_filelist(data_path, temp_root=gettempdir()):
    """
    Build the file list
    Iterate over all networks and station folders in the root directory
    or archive. For each ismn station the according static metadata is
    loaded and stored in the file handler together with the specific meta
    data for each file.
    In the file list, for faster filtering, some metadata information is
    stored in separate columns, i.e. network/station names, variables and
    depths.

    Parameters
    ----------
    data_path : str or Path
        Root path where the data is stored.
    temp_root : str or Path
        Root path where temporary files are stored.
    """

    root = IsmnRoot(data_path)

    # todo: basic info, not necessary?
    filelist = {'network': [], 'station': [], 'sensor': [],
                'variable': [], 'depth_from' : [], 'depth_to': [],
                'root_path': [], 'file_path': [],  # file info, to reload file
                'filehandler': [], # filehandler object
                }

    for net_dir, stat_dirs in root.cont.items():
        for stat_dir in stat_dirs:
            print(net_dir, stat_dir)

            data_files = root.find_files(stat_dir, '*.stm')
            static_meta = None

            for file_path in data_files:
                try:
                    f = DataFile(root, file_path, False,
                                 static_meta=static_meta,
                                 temp_root=temp_root)
                     #f.close() # to store pickle cut connection to archive
                except IOError as e:
                    logger.error(f'Error loading ismn file: {e}')
                    continue

                filelist['network'].append(f.metadata['network'].val)
                filelist['station'].append(f.metadata['station'].val)
                filelist['sensor'].append(f.metadata['sensor'].val)
                filelist['variable'].append(f.metadata['variable'].val)
                filelist['depth_from'].append(f.metadata['sensor'].depth.start)
                filelist['depth_to'].append(f.metadata['sensor'].depth.end)
                filelist['root_path'].append(f.root.path)
                filelist['file_path'].append(f.file_path)
                filelist['filehandler'].append(f)

    files = pd.DataFrame.from_dict(filelist)
    files.index = files.index.astype('int')
    files.ismn_pkg_version = pkg_version

    return files

class IsmnFileCollection(object):

    """
    The IsmnFileCollection class contains a pandas data frame with all ismn files
    in the given data directory. The file list can be loaded from a previously
    stored item, or built by iterating over all files in the data root.
    This class also contains function to filter the file list for faster access
    to files for a certain network/variable etc.

    Parameters
    ----------
    data_path : str or Path
        Root path of ISMN files or path to metadata pkl file.
        i.e. path to the downloaded zip file or the extracted zip directory (faster)
        or a file list that contains these infos already.
    filelist : pd.DataFarme, optional (Default: None)
        A pre-built filelist to use, if None is give one from data_path is created.
    meta_path : str or Path
    """

    def __init__(self, data_path, filelist=None, temp_root=gettempdir()):

        self.root = None

        if not os.path.exists(temp_root):
            os.makedirs(temp_root, exist_ok=True)

        self.temp_root = temp_root

        data_path = Path(data_path)
        self.root = IsmnRoot(data_path)

        self.data_loaded = False

        if filelist is not None:
            assert data_path == self._rootpath_from_filelist(filelist), \
                   "Root paths dont match"
            self.files = filelist
        else:
            self.files = build_filelist(data_path, temp_root)

    @classmethod
    def from_filelist(cls, filelist, temp_root=gettempdir()):
        """
        Create FileCollection from a pre-built file list.
        
        Parameters
        ----------
        filelist : pd.DataFrame
            File list that go into the file collection
        """
        data_path = cls._rootpath_from_filelist(filelist)

        return cls(data_path, filelist=filelist, temp_root=temp_root)

    @classmethod
    def from_pkl(cls, filepath, temp_root=gettempdir()):
        """
        Load a previously created and stored filelist from pkl.

        Parameters
        ----------
        filepath : str or Path
            Path to a pkl file containing a file list data frame.
        """

        filelist = pd.read_pickle(filepath)

        # if self.root.path != Path(data_path):
        #     raise MetadataError("Root path in filelist list does not match to passed root path."
        #                         "Try re-creating ISMN metadata.")

        if hasattr(filelist, 'ismn_pkg_version'):
            if filelist.ismn_pkg_version.lower() == 'unknown':
                warnings.warn("Metadata was built with an unknown package version. "
                             f"You should remove {filepath} and rebuild the metadata.")

            if filelist.ismn_pkg_version != pkg_version:
                warnings.warn("Metadata was build with a different package version. "
                              f"You should remove {filepath} and rebuild the metadata.")
                
        return cls.from_filelist(filelist, temp_root=temp_root)
    
    @staticmethod
    def _rootpath_from_filelist(filelist):
        root = np.unique(filelist['root_path'].values)
        assert len(root) == 1, "Found multiple root dirs"
        return Path(root[0])
        
    def filter_col_val(self, col, vals, return_index=False):
        """
        Filter the file list for certain values in the columns, except the
        filehandler column and the depth range.

        Parameters
        ----------
        col : str
            Column based on which the filtering is performed.
        vals : Any
            Value(s) that are allowed in col.
        return_index : bool, optional (default: False)
            Return only the index, no the filtered data frame.

        Returns
        -------
        filtered_filelist or filtered_index : pd.DataFrame or np.array
            Filtered file list or indices of included elements
        """

        if col not in self.files.columns:
            raise ValueError(f"Column {col} is not in file list.")

        if col in ['filehandler']:
            raise ValueError(f"Cannot filter based on column {col}.")

        mask = np.isin(self.files[col], np.atleast_1d(vals))
        idx = self.files.loc[mask].index.values

        if return_index:
            return idx
        else:
            return self.files.loc[idx, :]

    def filter_depth(self, min_depth=-np.inf, max_depth=np.inf, return_index=False,
                     only_consider_depth_from=False):
        """
        Filter filelist by depth_from and depth_to columns.

        Parameters
        ----------
        min_depth : float, optional (default: -np.inf)
            Return files below this depth, i.e. with depth_from greater than this
            are kept.
        max_depth : float, optional (default: np.inf)
            Return files above this depth, i.e. with depth_to smaller than this
            are kept. If only_consider_depth_from is selected, files with
            depth_to smaller than this are kept.
        return_index : bool, optional (default: False)
            Return only the index, no the filtered data frame.
        only_consider_depth_from : bool, optional (default: False)
            Only check whether the depth_from of the file is within the passed
            depth range.

        Returns
        -------
        filtered_filelist or filtered_index : pd.DataFrame or np.array
            Filtered file list or indices of included elements
        """
        mask_from = np.greater_equal(self.files['depth_from'], min_depth)

        if only_consider_depth_from:
            mask_to = np.less_equal(self.files['depth_from'], max_depth)
        else:
            mask_to = np.less_equal(self.files['depth_to'], max_depth)

        idx = self.files.loc[(mask_from & mask_to)].index.values

        if return_index:
            return idx
        else:
            return self.files.loc[idx, :]

    def filter_metadata(self, filter_dict:dict, filelist=None, return_index=False):
        """
        Filter file list by comparing file metadata to passed metadata dict.

        Parameters
        ----------
        filter_dict: dict
            Additional metadata keys and values for which the file list is filtered
            e.g. {'lc_2010': 10} to filter for a landcover class.
        filelist : pd.DataFrame, optional (default: None)
            The filelist to filter, if None is passed, self.files is used as
            for the other filter functions.
        return_index : bool, optional (default: False)
            Return only the index, no the filtered data frame.

        Returns
        -------
        filtered_filelist or filtered_index : pd.DataFrame or np.array
            Filtered file list or indices of included elements
        """
        # FIler based on the metadata in the filehander

        if filelist is None:
            filelist = self.files

        filehandlers = filelist['filehandler'].values

        mask = []
        for filehandler in filehandlers:
            flags = []
            for meta_key, meta_vals in filter_dict.items():
                meta_vals = np.atleast_1d(meta_vals)

                if meta_key not in filehandler.metadata.keys():
                    raise ValueError(f"{meta_key} is not a valid metadata variable")

                # check if the metadata val is one of the passed allowed values
                flag = filehandler.metadata[meta_key].val in meta_vals
                flags.append(flag)
            mask.append(all(flags))

        idx = filelist.loc[mask].index.values

        if return_index:
            return idx
        else:
            return filelist.loc[idx, :]

    def store(self, filename, drop_filehandler=False):
        """
        Store file list data frame to pkl, only possible if no data was loaded.

        Parameters
        ----------
        filename : str or Path
            Path to the file that will be created.
        drop_filehandler : bool, optional (default: False)
            Drow the file handlers from the dataframe before storing.
        """

        # todo: remove data when loaded? allow loaded data
        if self.data_loaded:
            raise IOError("Cannot store filelist when data was loaded.")

        self.close()

        if pkg_version == 'unknown':
            warnings.warn('Storing metadata from unknown package version')

        if drop_filehandler:
            self.files.drop(columns='filehandler').to_pickle(filename)
        else:
        # todo: add functio to add filehandler to list with paths only?
            self.files.to_pickle(filename)

    def close(self):
         # close root and all filehandlers
        self.root.close()
        for f in self.files['filehandler'].values:
            f.close()

class NetworkCollection(object):
    """
    The NetworkCollection builds ISMN Networks and fills them with stations
    and sensors using the passed file collection.
    A grid is added that contains all stations to perform spatial searches.
    """
    def __init__(self, file_collection:IsmnFileCollection, networks=None,
                 load_all=False):

        """
        Create network collection from previously loaded IsmnFileCollection.

        Parameters
        ----------
        filelist : pd.DataFrame
            DataFrame that contains the files to read.
            As in IsmnFileCollection.files or as returned by the filter functions.
        networks : list or str, optional (default: None)
            List of networks to activate.
        load_all : bool, optional (default: False)
            Load data into memory for all networks directly. This makes subsequent
            reading faster.
        """

        self.full_collection = file_collection
        if networks is None:
            self.active_collection = self.full_collection
        else:
            active_list = self.full_collection.filter_col_val('network', networks)
            self.active_collection = IsmnFileCollection.from_filelist(active_list)

        self.networks, self.grid = self._collect_networks(load_data=load_all)

    @property
    def files(self):
        # get files list
        return self.active_collection.files

    def _collect_networks(self, load_data=False) -> (dict, BasicGrid):
        # build networks and fill them with stations and sensors and apply filehandlers
        # for data reading.

        filelist = self.active_collection.files
        networks = {}
        points = []  # idx, lon, lat

        for idx, row in filelist.iterrows():
            f = row['filehandler']

            if load_data:
                f.load_data()

            nw_name, st_name, se_name = f['network'].val, f['station'].val, f['sensor'].val

            if nw_name not in networks:
                networks[nw_name] = Network(nw_name)

            if st_name not in networks[nw_name].stations:
                networks[nw_name].add_station(st_name,
                                              f['longitude'].val,
                                              f['latitude'].val,
                                              f['elevation'].val)
                points.append((idx, f['longitude'].val, f['latitude'].val))

            if se_name not in networks[nw_name].stations[st_name].sensors:
                networks[nw_name].stations[st_name]. \
                    add_sensor(se_name, f['variable'].val, f['variable'].depth, f)

        points = np.array(points)

        # todo: could use subset and create grid for full filelist?
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

        net, stat = self.full_collection.files.loc[idx, ['network', 'station']]
        return self.networks[net].stations[stat]

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

        idx_1 = self.active_collection.filter_col_val('variable', variable,
                                                      return_index=True)
        idx_2 = self.active_collection.filter_depth(min_depth=min_depth,
                                                    max_depth=max_depth,
                                                    return_index=True)
        idx = np.intersect1d(idx_1, idx_2)

        if filter_static_vars is not None:
            # use the prefiltered list to reduce number of loops
            prefiltered_list = self.active_collection.files.loc[idx]
            idx = self.active_collection.filter_metadata(
                filter_static_vars, filelist=prefiltered_list, return_index=True)

        return idx

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
        depth : Depth, optional (default: None)
            Sensing depth.

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
                    if (variable not in [None, se.variable] and
                            depth.encloses(se.depth)):
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
        gpi, dist = self.grid.find_nearest_gpi(lon, lat, max_dist=max_dist)
        station = self.station4idx(gpi)

        return station, dist
    #
    #
    # # todo: into other class?
    # # @filelist
    # def load_all(self, network=None):
    #     # todo: use the filehandler (bad if code changes) or paths and recreate handler
    #     if network is None:
    #         filelist = self.files
    #     else:
    #         filelist = self.files.loc[self.files['network'] == network]
    #
    #     for f in filelist:
    #         print(f'Loading: {f}')
    #         f.load_data()
    #
    #     self.data_loaded = True
    #
    # # @filelist
    # def unload_all(self):
    #     for f in self.files['filehandler']:
    #         f.data = None
    #     self.data_loaded = False

if __name__ == '__main__':
    filelist= build_filelist(r"C:\Temp\delete_me\ismn\testdata_ceop")
    coll = IsmnFileCollection.from_filelist(filelist)

    coll.store(r"C:\Temp\delete_me\meta.pkl")

    ocoll = IsmnFileCollection.from_pkl(r"C:\Temp\delete_me\meta.pkl")

    # coll2 = IsmnFileCollection.from_file(r"C:\Temp\delete_me\ismn_metadata.pkl")
    # print(str(coll2.files))