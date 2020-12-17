# -*- coding: utf-8 -*-

"""
Test file list functions, ie reading data for a filehander in the list,
and filtering the filelist.
"""

import os
import unittest
import pytest
import logging

from ismn.network_collection import NetworkCollection
from ismn.components import Depth

from tests.test_filecollection import cleanup

testdata_root = os.path.join(os.path.dirname(__file__), 'test_data')

class Test_NetworkCollectionCeopSepUnzipped(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(Test_NetworkCollectionCeopSepUnzipped, cls).setUpClass()

        testdata_path_unzipped = os.path.join(testdata_root,
            'Data_seperate_files_20170810_20180809')
        # clean existing metadata
        cleanup(os.path.join(testdata_path_unzipped, 'python_metadata'))

        # build metadata once
        NetworkCollection(testdata_path_unzipped, keep_loaded_data=False)

        cls.testdata_path_unzipped = testdata_path_unzipped

    def setUp(self) -> None:
        # load metadata for each test
        self.netcol = NetworkCollection(self.testdata_path_unzipped,
                                        keep_loaded_data=False)

        gpis, lons, lats = self.netcol.grid.get_grid_points()
        n_stats = 0
        for net in self.netcol.iter_networks():
            n_stats += net.n_stations()
        assert gpis.size == lons.size == lats.size == n_stats

    def tearDown(self) -> None:
        logging.shutdown()

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
        should_lon, should_lat = -156.62870, 71.32980

        station, dist = self.netcol.get_nearest_station(should_lon, should_lat)
        assert dist == 0
        assert station.lon == should_lon
        assert station.lat == should_lat
        gpi, dist = self.netcol.grid.find_nearest_gpi(int(should_lon),int(should_lat))
        assert dist != 0
        for net in self.netcol.iter_networks():
            if station.name in net.stations.keys():
                assert net.stations[station.name].lon == should_lon
                assert net.stations[station.name].lat == should_lat

        #station, dist = self.netcol.get_nearest_station(0, 0, max_dist=100)
        # todo: when fixed in pygeogrids this should return nothing...
        #https://github.com/TUW-GEO/pygeogrids/issues/64
        #assert station == None

    def test_eval_sensor_metadata(self):
        # test based on metadata
        station = self.netcol.networks['COSMOS'].stations['Barrow-ARM']
        assert station.sensors[1].eval('soil_moisture', Depth(0,1),
                filter_meta_dict={'lc_2010': 210, 'climate_KG': 'ET'})
        assert not station.sensors[1].eval('soil_moisture', Depth(0,1),
                    filter_meta_dict={'lc_2010': 999})

class Test_NetworkCollectionHeaderValuesUnzipped(Test_NetworkCollectionCeopSepUnzipped):

    @classmethod
    def setUpClass(cls):
        super(Test_NetworkCollectionHeaderValuesUnzipped, cls).setUpClass()

        testdata_path_unzipped = os.path.join(testdata_root,
            'Data_seperate_files_header_20170810_20180809')
        # clean existing metadata
        cleanup(os.path.join(testdata_path_unzipped, 'python_metadata'))

        # build metadata once
        NetworkCollection(testdata_path_unzipped, keep_loaded_data=False)

        cls.testdata_path_unzipped = testdata_path_unzipped

    def setUp(self) -> None:
        self.netcol = NetworkCollection(self.testdata_path_unzipped,
                                        keep_loaded_data=False)

        gpis, lons, lats = self.netcol.grid.get_grid_points()
        n_stats = 0
        for net in self.netcol.iter_networks():
            n_stats += net.n_stations()
        assert gpis.size == lons.size == lats.size == n_stats

@pytest.mark.zip
class Test_NetworkCollectionCeopSepZipped(Test_NetworkCollectionCeopSepUnzipped):

    @classmethod
    def setUpClass(cls):
        super(Test_NetworkCollectionCeopSepZipped, cls).setUpClass()

        testdata_path = os.path.join(testdata_root, 'zip_archives', 'ceop')
        testdata_zip_path = os.path.join(testdata_path,
            'Data_seperate_files_20170810_20180809.zip')

        # clean up existing metadata
        metadata_path = os.path.join(testdata_path, 'python_metadata')
        cleanup(metadata_path)

        # build metadata once
        NetworkCollection(testdata_zip_path, meta_path=metadata_path,
                          keep_loaded_data=False)
        cls.testdata_zip_path = testdata_zip_path

    def setUp(self) -> None:
        self.netcol = NetworkCollection(self.testdata_zip_path,
                                        keep_loaded_data=True)

        gpis, lons, lats = self.netcol.grid.get_grid_points()
        n_stats = 0
        for net in self.netcol.iter_networks():
            n_stats += net.n_stations()
        assert gpis.size == lons.size == lats.size == n_stats

@pytest.mark.zip
class Test_NetworkCollectionHeaderValuesZipped(Test_NetworkCollectionCeopSepUnzipped):

    @classmethod
    def setUpClass(cls):
        super(Test_NetworkCollectionHeaderValuesZipped, cls).setUpClass()

        testdata_path = os.path.join(testdata_root, 'zip_archives', 'header')
        testdata_zip_path = os.path.join(testdata_path,
            'Data_seperate_files_header_20170810_20180809.zip')

        # clean up existing metadata
        metadata_path = os.path.join(testdata_path, 'python_metadata')
        cleanup(metadata_path)

        NetworkCollection(testdata_zip_path, meta_path=metadata_path,
                          keep_loaded_data=False)

        cls.testdata_zip_path = testdata_zip_path

    def setUp(self) -> None:

        self.netcol = NetworkCollection(self.testdata_zip_path,
                                        keep_loaded_data=True)

        gpis, lons, lats = self.netcol.grid.get_grid_points()
        n_stats = 0
        for net in self.netcol.iter_networks():
            n_stats += net.n_stations()
        assert gpis.size == lons.size == lats.size == n_stats

if __name__ == '__main__':
    unittest.main()

