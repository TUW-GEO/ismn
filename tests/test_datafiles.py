# -*- coding: utf-8 -*-

"""
Module that tests filehandler classes for reading ismn data.
"""

import os
import unittest

from ismn.files import DataFile
from ismn.meta import MetaData

from pathlib import Path
from datetime import datetime

testdata_path = Path(os.path.join(os.path.dirname(__file__), 'test_data'))


class DataFileCeopSepTestUnzipped(unittest.TestCase):
    # from dir, no load_data
    def setUp(self) -> None:
        self.station_name = "Barrow-ARM"
        self.network_name = "COSMOS"
        self.instrument = "Cosmic-ray-Probe"
        self.variable = 'soil_moisture'
        self.depth_from = 0.0
        self.depth_to = 0.21

        self.longitude = -156.62870
        self.latitude = 71.32980
        self.elevation = 4.

        root = testdata_path / "Data_seperate_files_20170810_20180809"
        filepath = Path(self.network_name,  self.station_name,
            f"{self.network_name}_{self.network_name}_{self.station_name}_sm_{self.depth_from:.6f}_{self.depth_to:.6f}_{self.instrument}_20170810_20180809.stm")

        self.file = DataFile(root, filepath, load_data=True)

        self.data_should_201708113 = {
            "datetime": "2017-08-11 13:00:00",
            self.variable: 0.183,
            f"{self.variable}_flag": "G",
            f"{self.variable}_orig_flag": "M"}

    def test_metadata(self):
        """ test reading the loaded metadata for a file (incl. static meta) """
        lc = self.file.metadata['lc_2010'].val
        assert lc == 210

        flag = self.file.check_metadata(self.variable, self.depth_from, self.depth_to,
                                        filter_static_vars={'lc_2010': lc})
        assert flag == True

        flag = self.file.check_metadata('nonexistingvar', 0, 0.1,
                                        filter_static_vars={'lc_2010': lc})
        assert flag == False

        flag = self.file.check_metadata(self.variable, 98., 99.)
        assert flag == False

        assert self.file.metadata['station'].val == self.station_name
        assert self.file.metadata['network'].val == self.network_name
        assert self.file.metadata['sensor'].depth.start == self.depth_from
        assert self.file.metadata['sensor'].depth.end == self.depth_to
        assert self.file.metadata['sensor'].val == self.instrument
        assert self.file.metadata['longitude'].val == self.longitude
        assert self.file.metadata['latitude'].val == self.latitude
        assert self.file.metadata['variable'].val == self.variable
        assert self.file.metadata['sand_fraction'].val == 34.
        assert self.file.metadata['sand_fraction'].depth.start == 0.
        assert self.file.metadata['sand_fraction'].depth.end == 0.3


    def test_data(self):
        """ test reading the actual data for a file """
        # todo: why is sm column called "variable"?
        timestamp = datetime(2017,8,11,13)

        data_is = self.file.data.loc[timestamp]

        assert data_is[self.variable] == self.data_should_201708113[self.variable]
        assert data_is[f"{self.variable}_flag"] == self.data_should_201708113[f"{self.variable}_flag"]
        assert data_is[f"{self.variable}_orig_flag"] == self.data_should_201708113[f"{self.variable}_orig_flag"]


    def test_metadata_for_depth(self):
        """ Check finding best matching metadata for file """
        print(self.file.root.path)
        bestmeta = self.file.read_metadata(best_meta_for_sensor=True)
        allmeta = self.file.read_metadata(best_meta_for_sensor=False)

        assert len([k for k in allmeta.keys() if k == 'saturation']) > 1
        assert len([k for k in bestmeta.keys() if k == 'saturation']) == 1

        assert bestmeta['saturation'].depth.start == bestmeta['sand_fraction'].depth.start == 0.
        assert bestmeta['saturation'].depth.end == bestmeta['sand_fraction'].depth.end ==  0.3

        assert allmeta['saturation'][0].depth == bestmeta['saturation'].depth

        assert allmeta['saturation'][1].depth.start == 0.3
        assert allmeta['saturation'][1].depth.end == 1.0

        additional_metadata = MetaData.from_dict({'somevar' : 'test'},)


        addmeta = self.file.read_metadata(static_meta=additional_metadata,
                                          best_meta_for_sensor=False)

        assert addmeta['somevar'].val == 'test'

        for var in addmeta.metadata:
            if var.name == 'somevar': continue # todo: add variable drop() to metadata?
            assert var in allmeta


class DataFileCeopSepTestZipped(DataFileCeopSepTestUnzipped):
    # same as test from unzipped archive, but uses zip directly. Runs all tests.
    def setUp(self) -> None:
        super(DataFileCeopSepTestZipped, self).setUp()

        root = testdata_path / "zip_archives" / "ceop" / "Data_seperate_files_20170810_20180809.zip"
        filepath = Path(self.network_name,  self.station_name,
        f"{self.network_name}_{self.network_name}_{self.station_name}_sm_{self.depth_from:.6f}_{self.depth_to:.6f}_{self.instrument}_20170810_20180809.stm")

        self.file = DataFile(root, filepath, load_data=True)


class DataFileHeaderValuesTestUnzipped(DataFileCeopSepTestUnzipped):
    # same as test from unzipped archive, but uses zip directly. Runs all tests.
    def setUp(self) -> None:
        super(DataFileHeaderValuesTestUnzipped, self).setUp()

        root = testdata_path / "Data_seperate_files_header_20170810_20180809"
        filepath = Path(self.network_name,  self.station_name,
        f"{self.network_name}_{self.network_name}_{self.station_name}_sm_{self.depth_from:.6f}_{self.depth_to:.6f}_{self.instrument}_20170810_20180809.stm")

        self.file = DataFile(root, filepath, load_data=True)


class DataFileHeaderValuesTestZipped(DataFileCeopSepTestUnzipped):
    # same as test from unzipped archive, but uses zip directly. Runs all tests.
    def setUp(self) -> None:
        super(DataFileHeaderValuesTestZipped, self).setUp()

        root = testdata_path / "zip_archives" / "header" / "Data_seperate_files_header_20170810_20180809.zip"
        filepath = Path(self.network_name,  self.station_name,
        f"{self.network_name}_{self.network_name}_{self.station_name}_sm_{self.depth_from:.6f}_{self.depth_to:.6f}_{self.instrument}_20170810_20180809.stm")

        self.file = DataFile(root, filepath, load_data=True)


# todo: test from zip, load_data
if __name__ == '__main__':
    unittest.main()