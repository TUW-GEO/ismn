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

import sys
import platform

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
            self.file_collection = IsmnFileCollection.from_scratch(
                root, parallel=True, temp_root=temp_root)
            self.file_collection.to_metadata_csv(meta_csv_file)

        self.keep_loaded_data = keep_loaded_data

        self.networks, self.grid = self._collect_networks(networks)

    @property
    def files(self):
        return self.file_collection.files

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


if __name__ == '__main__':
    networks = NetworkCollection(r"D:\data-read\ISMN\global_20191024",
                                 keep_loaded_data=False,
                                 networks=None)

    networks.plot_station_locations(variable='precipitation', min_depth=-np.inf, max_depth=np.inf,
                                               stats_text=True, filename=r"C:\Temp\delete_me\map_rainf.png",
                                               check_only_sensor_depth_from=False)


