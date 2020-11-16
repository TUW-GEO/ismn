# -*- coding: utf-8 -*-

"""
Test file list functions, ie reading data for a filehander in the list,
and filtering the filelist.
"""

import os
import unittest

from ismn.groupings import NetworkCollection, IsmnFileCollection
from ismn.components import Depth

from tests.test_filecollection import cleanup

testdata_root = os.path.join(os.path.dirname(__file__), 'test_data')


class Test_NetworkCollectionCeopSep(unittest.TestCase):

    def setUp(self) -> None:

        testdata_path_unzipped = os.path.join(testdata_root,
            'Data_seperate_files_20170810_20180809')
        # clean existing metadata
        cleanup(os.path.join(testdata_path_unzipped, 'python_metadata'))

        self.netcol = NetworkCollection(IsmnFileCollection(testdata_path_unzipped),
                                        load_all=True)

        gpis, lons, lats = self.netcol.grid.get_grid_points()
        assert gpis.size == lons.size == lats.size == self.netcol.files.index.size


    def test_station4idx(self):
        station = self.netcol.station4idx(0)
        assert station.name == 'ARM-1'
        assert len(station.sensors) == 1

        station = self.netcol.station4idx(1)
        assert station.name == 'Barrow-ARM'
        assert len(station.sensors) == 1

    def test_get_dataset_ids(self):
        ids = self.netcol.get_dataset_ids('soil_moisture', 0., 0.19) # should get 1
        assert len(ids) == 1

        ids = self.netcol.get_dataset_ids('soil_moisture', 0., 1.) # should get 2
        assert len(ids) == 2

        ids = self.netcol.get_dataset_ids('soil_moisture', 0., 1.,
            filter_static_vars={'lc_2010': 210}) # should get 1
        assert len(ids) == 1

        ids = self.netcol.get_dataset_ids('nonexisting') # should get 0
        assert len(ids) == 0

    def test_get_sensors(self):
        i = 0
        for se in self.netcol.get_sensors('COSMOS'):
            data = se.read_data()
            # check if the networks is COSMOS or station in [ARM, Barrow-ARM]
            assert not data.empty
            # check something for that one station
            i += 1
        assert i == 2

        i = 0
        for se in self.netcol.get_sensors(station='Barrow-ARM'):
            data = se.read_data()
            assert not data.empty
            # check something for that one station
            i += 1
        assert i == 1

        i = 0
        for se in self.netcol.get_sensors(depth=Depth(0,1)):
            data = se.read_data()
            assert not data.empty
            i +=1
        assert i == 2

        for se in self.netcol.get_sensors(variable='nonexisting'):
            raise ValueError("Found sensor, although none should exist")

    def test_get_nearest_station(self):
        station, dist = self.netcol.get_nearest_station(-156.62870,71.32980)
        assert dist == 0
        assert station.lon == -156.62870
        assert station.lat == 71.32980
        gpi, dist = self.netcol.grid.find_nearest_gpi(-156,70)
        assert dist != 0
        assert self.netcol.files.loc[int(gpi)]['filehandler'].metadata['longitude'].val == \
               station.lon
        assert self.netcol.files.loc[int(gpi)]['filehandler'].metadata['latitude'].val == \
               station.lat

        station, dist = self.netcol.get_nearest_station(0,0, max_dist=100)
        # todo: when fixed in pygeogrids this should return nothing...
        #https://github.com/TUW-GEO/pygeogrids/issues/64
        #assert station == None


if __name__ == '__main__':
    unittest.main()

