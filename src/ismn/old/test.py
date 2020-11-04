# -*- coding: utf-8 -*-

"""
Module description
"""
# TODO:
#   (+) 
#---------
# NOTES:
#   -

import zipfile
import os
import numpy as np

filepath= r"C:\Temp\delete_me\ismn\Data_separate_files_20101001_20201005_5712_GMIw_20201005.zip"


class NetworkCollection():

    def __init__(self, networks: list = None):
        """
        A network collection contains multiple networks and a grid with all
        stations for spatial searches.

        Parameters
        ----------
        networks: list[Network], optional (default: None)
            A list of Network objects that create the (initial) collection.
        load_data : bool, optional (default: True)
            Load and keep ismn data into memory as soon as possible. Accessing
            already loaded data is faster than re-loading.
        """
        self.networks = {}

        if networks is not None:
            for net in networks:
                self.networks[net.name] = net

        self.grid = None
        self._update_grid()

    def iter_networks(self):
        # iterate over current networks
        for name, net in self.networks.items():
            yield net

    def _update_grid(self):
        # build grid for current networks
        lons, lats = [], []
        for net in self.iter_networks():
            net_lons, net_lats = net.coords
            lons += net_lons
            lats += net_lats

        self.grid = BasicGrid(lons, lats)

    @classmethod
    def from_filelist(self, filelist, load_data=False):
        """
        Create network collection from previously loaded IsmnFileCollection.

        Parameters
        ----------
        filelist : pd.DataFrame
            DataFrame that contains the files to read.
            As in IsmnFileCollection.files or as returned by the filter functions.
        load_data : bool, optional (default: False)
            Load data into memory for all networks

        Returns
        -------

        """
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

        self.networks = networks
        # self.grid_lut = points[:,0]
        self.grid = BasicGrid(points[:, 1], points[:, 2])  # todo: could use subset and create grid for full filelist?

    def add_network(self, network):
        """
        Add another network to the collection.

        Parameters
        ----------
        network : Network
            Network to add, must not exist in collection
        """
        if network.name in self.networks.keys():  # todo: or warn?
            raise ValueError("Network {network} already exists in collection.")

        self.networks[network.name] = network

        self._update_grid()

@filelist
def filter_col_network(self, network):
    """
    Filter the file list for files from passed network(s).

    Parameters
    ----------
    network : list or str
        Name of the network to filter file list for

    Returns
    -------
    filtered_list : pd.DataFrame
        The filtered original data frame.
    """

    network = np.atleast_1d(network)
    net_in_coll = np.unique(self.files['network'].values)
    net_not_found = network[~np.in1d(network, net_in_coll)]
    if any(net_not_found):
        ISMNError(f"Network(s) not found: {net_not_found.tolist()}")
    filtered_filelist = self.files.loc[np.isin(self.files['network'], network)]

    return filtered_filelist

def walk_zipfile():
    root_path, dirs, files = [], [], []
    with zipfile.ZipFile(filepath) as zip:
        file_list = zip.namelist()
        path_files_list = [(os.path.split(x)[0], x) for x in file_list]
        sub_folders = [x[0] for x in path_files_list]
        unique_folders = np.unique(sub_folders)

        for sub in unique_folders:
            filenames = [item[1] for item in path_files_list if item[0] == sub]
            for filename in filenames:
                root_path.append(os.path.join())
                dirs.append(sub)
                files.append8


    def collect_from_folder(self):
        """
        function walks the data_path directory and looks for network
        folders and ISMN datafiles. It collects metadata for every
        file found and returns a numpy.ndarray of metadata

        Parameters
        ----------
        data_path : string
            root directory on filesystem where the ISMN data was unzipped to

        Returns
        -------
        metadata : numpy.ndarray
            structured numpy array which contains the metadata for one file per row
        """

        logging.basicConfig(filename=os.path.join(self.meta_path, 'metadata.log'),
                            level=logging.DEBUG)

        dtype = list(map(tuple, self._file_metadata_template[:, (0,2)])) + \
                list(map(tuple, self._csv_metadata_template[:, (0,2)]))

        metadata_catalog = []

        for root_path, dirs, files in os.walk(self.data_path):

            files.sort()

            # read additional metadata from csv file
            filename_csv = glob.glob('{}/*.csv'.format(root_path))
            # csv_meta_keys, csv_meta_entries = \
            #     self._csv_metadata_template[:,0], self._csv_metadata_template[:,1]

            if len(filename_csv) > 0:
                csv_meta_keys, csv_meta_entries = \
                    self.meta_from_csv(filename_csv[0])

                for filename in files:
                    if filename.endswith('.stm'):
                        subpath = os.path.relpath(root_path, self.data_path)
                        try:
                            file_meta_keys, file_meta_entries = \
                                self.meta_from_datafile(self.data_path,
                                                        subpath,
                                                        filename)
                        except (readers.ReaderException, IOError) as e:
                            continue

                        for var_meta_entry in file_meta_entries:
                            meta = var_meta_entry + csv_meta_entries
                            metadata_catalog.append(tuple(meta))
            else:
                if any(filename.endswith('.stm') for filename in files):
                    logging.info('No csv file available ({})'.format(root_path))
                else:
                    continue

        return np.array(metadata_catalog, dtype=dtype)

    def _walk_zipped_files(self, root_path) -> (str, Generator):
        # yield all files in data zip
        with zipfile.ZipFile(root_path) as zi:
            file_list = zi.namelist()
            path_files_list = [(os.path.split(x)[0], x) for x in file_list]
            sub_folders = [x[0] for x in path_files_list]
            unique_folders = np.unique(sub_folders)
            for sub in unique_folders:
                sub_filepaths = [item[1] for item in path_files_list if item[0] == sub]
                for filename in sub_filepaths:
                    yield root_path, filename


def list_network_dirs(ismn_data_path):
    """
    Get list of all network folders in an ismn archive (zip or extracted)
    """
    if zipfile.is_zipfile(ismn_data_path):
        with zipfile.ZipFile(ismn_data_path) as zi:
            file_list = zi.namelist()
            path_files_list = [os.path.split(os.path.split(x)[0])[0] for x in file_list]
            sub_folders = [x[0] for x in path_files_list if x[0] != '']
            network_folders = np.unique(sub_folders)
        return network_folders
    else:
        network_folders = []
        for f in os.scandir(ismn_data_path):
            if f.is_dir():
                if f.path == ismn_data_path:
                    continue
                else:
                    network_folders.append(os.path.relpath(f.path, ismn_data_path))
        return np.array(network_folders)
def list_station_dirs(ismn_data_path, network_dir_name, only_station=False):
    """
    Get list of all stations in a network folder of an ismn archive (zip or extracted)

    Parameters
    ----------
    ismn_data_path : str
        Path to ismn archive (zip or extracted)
    network_dir_name : str
        Name of ntwork dir in data_path to get station folder for
    only_station : bool, optional (default: False)
        If selected, station folder names are retured without network folder.

    Returns
    -------
    station_dirs : np.array
        Array of station dirs in network dir
    """

    if zipfile.is_zipfile(ismn_data_path):
        station_folders = []
        with zipfile.ZipFile(ismn_data_path) as zi:
            file_list = zi.namelist()
            for f in file_list:
                s = os.path.split(f)[0]
                if s != '':
                    if only_station:
                        station_folders.append(os.path.relpath(s, network_dir_name))
                    else:
                        station_folders.append(s)

        return np.unique(np.array(station_folders))

    else:
        station_folders = []
        for f in os.scandir(os.path.join(ismn_data_path, network_dir_name)):
            if f.is_dir():
                relpath = os.path.relpath(f.path, ismn_data_path)
                if only_station:
                    station_folders.append(os.path.basename(os.path.normpath(relpath)))
                else:
                    station_folders.append(relpath)

        return np.array(station_folders)


if meta_path is None:
    if self.from_zip:
        self.meta_path = os.path.join(os.path.dirname(self.data_path),
                                      'python_metadata')
    else:
        self.meta_path = os.path.join(self.data_path, 'python_metadata')
else:
    if self.from_zip:
        metadirname = f"python_metadata_" \
                      f"{os.path.basename(os.path.dirname(self.data_path))}"
    else:
        metadirname = f"python_metadata_" \
                      f"{os.path.basename(self.data_path)}"

    self.meta_path = os.path.join(meta_path, metadirname)