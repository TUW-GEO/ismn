
import numpy as np
from pygeogrids import BasicGrid

from ismn.file_collection import IsmnFileCollection
from ismn.components import Network, Depth

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


    # todo: into other class?
    # @filelist
    def load_all(self, network=None):
        # todo: use the filehandler (bad if code changes) or paths and recreate handler
        if network is None:
            filelist = self.files
        else:
            filelist = self.files.loc[self.files['network'] == network]

        for f in filelist:
            print(f'Loading: {f}')
            f.load_data()

        self.data_loaded = True

    # @filelist
    def unload_all(self):
        for f in self.files['filehandler']:
            f.data = None
        self.data_loaded = False

if __name__ == '__main__':
    filecoll = IsmnFileCollection(r"C:\Temp\delete_me\ismn\testdata_ceop")
    netcoll = NetworkCollection(filecoll, networks='FMI')
    netcoll.get_dataset_ids('soil_moisture', 0, 0.05, {'network': 'FMI', 'lc_2010': 70})
