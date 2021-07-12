# -*- coding: utf-8 -*-
import unittest
import os
from tempfile import TemporaryDirectory
from datetime import datetime

import numpy as np
import pytest
import logging
from collections import OrderedDict

from tests.test_filecollection import cleanup
from ismn.interface import ISMN_Interface
from ismn.meta import Depth

testdata_root = os.path.join(os.path.dirname(__file__), "test_data")

def test_metadata_dataframe():
    # make sure that metadata.index represents same values as get_dataset_ids
    testdata = os.path.join(testdata_root, "Data_seperate_files_20170810_20180809")
    metadata_path = os.path.join(testdata, "python_metadata")
    cleanup(metadata_path)
    ds_one = ISMN_Interface(testdata, meta_path=metadata_path, network='FR_Aqui')

    assert np.all(ds_one.metadata.index.values == ds_one.get_dataset_ids(None, -np.inf, np.inf))
    ids = ds_one.get_dataset_ids('soil_moisture')
    assert ids == ds_one.metadata.index.values
    assert ds_one.metadata.loc[ids[0], 'variable']['val'] == 'soil_moisture'
    assert ds_one.metadata.loc[ids[0], 'network']['val'] == 'FR_Aqui'
    ds_one.close_files()

class Test_ISMN_Interface_CeopUnzipped(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        super(Test_ISMN_Interface_CeopUnzipped, cls).setUpClass()

        testdata = os.path.join(testdata_root, "Data_seperate_files_20170810_20180809")
        metadata_path = os.path.join(testdata, "python_metadata")

        cleanup(metadata_path)
        ds = ISMN_Interface(testdata, network=[], parallel=True)
        assert ds.networks == OrderedDict()
        cls.testdata = testdata

    def setUp(self) -> None:
        self.ds = ISMN_Interface(self.testdata, network=["COSMOS"])

    def tearDown(self) -> None:
        self.ds.close_files()
        logging.shutdown()

    def test_list(self):
        with pytest.deprecated_call():
            assert len(self.ds.list_networks()) == 1
            assert len(self.ds.list_stations()) == len(self.ds.list_stations("COSMOS")) == 2
            assert len(self.ds.list_sensors()) == 2
            assert len(self.ds.list_sensors(station="Barrow-ARM")) == 1

    def test_network_for_station(self):
        assert self.ds.network_for_station("Barrow-ARM") == "COSMOS"
        assert self.ds.network_for_station("ARM-1") == "COSMOS"

    def test_stations_that_measure(self):
        for s in self.ds.stations_that_measure("soil_moisture"):
            assert s.name in ["ARM-1", "Barrow-ARM"]

        for s in self.ds.stations_that_measure("nonexisting"):
            raise AssertionError("Found var that doesnt exist")

    def test_get_dataset_ids(self):
        ids = self.ds.get_dataset_ids("soil_moisture", max_depth=100, groupby="network")
        assert list(ids.keys()) == ["COSMOS"]
        assert ids["COSMOS"] == [0, 1]

        ids = self.ds.get_dataset_ids("soil_moisture", max_depth=0.19)
        assert ids == [0]

        ids = self.ds.get_dataset_ids(
            ["soil_moisture"],
            max_depth=99,
            filter_meta_dict={
                "lc_2010": 210,
                "network": "COSMOS",
                "station": "Barrow-ARM",
            },
        )
        assert ids == [1]

        ids = self.ds.get_dataset_ids("novar")
        assert len(ids) == 0

        ids = self.ds.get_dataset_ids(["soil_moisture", "shouldhavenoeffect"], 0.0, 0.19)  # should get 1
        assert len(ids) == 1

        ids = self.ds.get_dataset_ids("soil_moisture", 0.0, 1.0)  # should get 2
        assert len(ids) == 2

        ids = self.ds.get_dataset_ids(
            "soil_moisture", 0.0, 1.0, filter_meta_dict={"lc_2010": 210}
        )  # should get 1
        assert len(ids) == 1

        ids = self.ds.get_dataset_ids("nonexisting")  # should get 0
        assert len(ids) == 0

    def test_read_ts(self):
        data1 = self.ds.read(0)
        assert not data1.empty

        data2, meta = self.ds.read_ts(1, return_meta=True)
        assert not data2.empty

    def test_read_metadata(self):
        data2, meta = self.ds.read_ts(1, return_meta=True)
        assert all(meta == self.ds.read_metadata(1, format="pandas"))
        assert self.ds.read_metadata(1, format="dict") is not None
        assert self.ds.read_metadata([1], format="obj") is not None

        assert not self.ds.metadata.empty
        assert self.ds.metadata.loc[1]['station']['val'] \
               == self.ds.read_metadata([0,1]).loc[1, ('station', 'val')]

    def test_find_nearest_station(self):
        should_lon, should_lat = -156.62870, 71.32980

        station = self.ds.find_nearest_station(should_lon, should_lat)

        assert station.lon == should_lon
        assert station.lat == should_lat

    def test_plot_station_locations(self):
        with TemporaryDirectory() as out_dir:
            outpath = os.path.join(out_dir, "plot.png")
            self.ds.plot_station_locations(
                ["soil_moisture", 'precipitation'], markersize=5, filename=outpath
            )
            assert len(os.listdir(out_dir)) == 1

    def test_get_min_max_obs_timestamps(self):
        tmin, tmax = self.ds.get_min_max_obs_timestamps("soil_moisture", max_depth=0.19)
        assert tmin == datetime(2017, 8, 10, 0)
        assert tmax == datetime(2018, 8, 9, 23)

    def test_get_min_max_obs_timestamps_for_station(self):
        station = self.ds.collection.networks["COSMOS"].stations["ARM-1"]
        tmin, tmax = station.get_min_max_obs_timestamp("soil_moisture", 0, 0.19)
        assert tmin == datetime(2017, 8, 10, 0)
        assert tmax == datetime(2018, 8, 9, 23)

    def test_get_static_var_val(self):
        vals = self.ds.get_static_var_vals("soil_moisture", max_depth=0.19)
        assert vals == {130: "Grassland"}

        vals = self.ds.get_landcover_types("soil_moisture", max_depth=100)
        assert len(vals) == 2
        assert vals[130] == "Grassland"
        assert vals[210] == "Water"
        self.ds.print_landcover_dict()

        vals = self.ds.get_climate_types(
            "soil_moisture", max_depth=100, climate="climate_KG"
        )
        assert len(vals) == 2
        assert vals["ET"] == "Polar Tundra"
        assert vals["Cfa"] == "Temperate Without Dry Season, Hot Summer"
        self.ds.print_climate_dict()

    def test_get_var(self):
        vars = self.ds.get_variables()
        assert vars == ["soil_moisture"]

    def test_get_sensors(self):
        i = 0
        for nw, station in self.ds.collection.iter_stations(
            filter_meta_dict={"network": "COSMOS"}
        ):
            for se in station.iter_sensors():
                data = se.read_data()
                # check if the networks is COSMOS or station in [ARM, Barrow-ARM]
                assert not data.empty
                # check something for that one station
                i += 1
        assert i == 2

        i = 0
        for se in self.ds.networks["COSMOS"].stations["Barrow-ARM"].iter_sensors():
            data = se.read_data()
            assert not data.empty
            # check something for that one station
            i += 1
        assert i == 1

        i = 0
        for net, stat, sens in self.ds.collection.iter_sensors(
            depth=Depth(0, 1),
            filter_meta_dict={"station": ["Barrow-ARM", "ARM-1"]},
        ):
            data = sens.read_data()
            assert not data.empty
            i += 1
        assert i == 2

        for nw, station in self.ds.collection.iter_stations():
            for se in station.iter_sensors(variable="nonexisting"):
                raise ValueError("Found sensor, although none should exist")

    def test_get_nearest_station(self):
        should_lon, should_lat = -156.62870, 71.32980

        station, dist = self.ds.collection.get_nearest_station(should_lon, should_lat)
        assert dist == 0
        assert station.lon == should_lon
        assert station.lat == should_lat
        gpi, dist = self.ds.collection.grid.find_nearest_gpi(
            int(should_lon), int(should_lat)
        )
        assert dist != 0
        for net in self.ds.collection.iter_networks():
            if station.name in net.stations.keys():
                assert net.stations[station.name].lon == should_lon
                assert net.stations[station.name].lat == should_lat

        station, dist = self.ds.find_nearest_station(
            0, 0, return_distance=True, max_dist=100
        )
        assert station == dist == None

    def test_citation(self):
        with TemporaryDirectory() as out_dir:
            out_file = os.path.join(out_dir, 'citation.txt')
            refs = self.ds.collection.export_citations(out_file=out_file)
            assert all([net in refs.keys() for net in list(self.ds.collection.networks.keys())])
            assert os.path.exists(out_file)
            with open(out_file, mode='r') as f:
                lines = f.readlines()
                assert len(lines) > 0

class Test_ISMN_Interface_HeaderValuesUnzipped(Test_ISMN_Interface_CeopUnzipped):
    @classmethod
    def setUpClass(cls):
        super(Test_ISMN_Interface_HeaderValuesUnzipped, cls).setUpClass()

        testdata_path_unzipped = os.path.join(
            testdata_root, "Data_seperate_files_header_20170810_20180809"
        )
        # clean existing metadata

        metadata_path = os.path.join(testdata_path_unzipped, "python_metadata")

        cleanup(metadata_path)

        ISMN_Interface(testdata_path_unzipped)

        cls.testdata = testdata_path_unzipped

    def setUp(self) -> None:

        self.ds = ISMN_Interface(self.testdata)


@pytest.mark.zip
class Test_ISMN_Interface_CeopZipped(Test_ISMN_Interface_CeopUnzipped):
    @classmethod
    def setUpClass(cls):
        super(Test_ISMN_Interface_CeopZipped, cls).setUpClass()

        testdata_path = os.path.join(testdata_root, "zip_archives", "ceop")
        testdata_zip_path = os.path.join(
            testdata_path, "Data_seperate_files_20170810_20180809.zip"
        )
        # clean up existing metadata
        metadata_path = os.path.join(testdata_path, "python_metadata")
        cleanup(metadata_path)

        ISMN_Interface(testdata_zip_path)

        cls.testdata_zip_path = testdata_zip_path

    def setUp(self) -> None:
        self.ds = ISMN_Interface(self.testdata_zip_path)


@pytest.mark.zip
class Test_ISMN_Interface_HeaderValuesZipped(Test_ISMN_Interface_CeopUnzipped):
    @classmethod
    def setUpClass(cls):
        super(Test_ISMN_Interface_HeaderValuesZipped, cls).setUpClass()

        testdata_path = os.path.join(testdata_root, "zip_archives", "header")
        testdata_zip_path = os.path.join(
            testdata_path, "Data_seperate_files_header_20170810_20180809.zip"
        )
        # clean up existing metadata
        metadata_path = os.path.join(testdata_path, "python_metadata")
        cleanup(metadata_path)

        ISMN_Interface(testdata_zip_path)

        cls.testdata_zip_path = testdata_zip_path

    def setUp(self) -> None:

        self.ds = ISMN_Interface(self.testdata_zip_path)


if __name__ == "__main__":
    unittest.main()
