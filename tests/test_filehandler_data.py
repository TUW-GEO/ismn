# -*- coding: utf-8 -*-

"""
Module that tests filehandler classes for reading ismn data.
"""

import os
import unittest

from ismn.filehandlers import DataFile
from ismn.meta import MetaData, Depth

from pathlib import Path
from datetime import datetime

testdata_path = Path(os.path.join(os.path.dirname(__file__), "test_data"))


class Test_DataFileCeopSepUnzipped(unittest.TestCase):
    # from dir, no _load_data
    def setUp(self) -> None:
        self.station_name = "Barrow-ARM"
        self.network_name = "COSMOS"
        self.instrument = "Cosmic-ray-Probe"
        self.variable = "soil_moisture"
        self.depth_from = 0.0
        self.depth_to = 0.21

        self.longitude = -156.62870
        self.latitude = 71.32980
        self.elevation = 4.0

        root = testdata_path / "Data_seperate_files_20170810_20180809"
        filepath = Path(
            self.network_name,
            self.station_name,
            f"{self.network_name}_{self.network_name}_{self.station_name}_sm_{self.depth_from:.6f}_{self.depth_to:.6f}_{self.instrument}_20170810_20180809.stm",
        )

        self.file = DataFile(root, filepath)

        self.data_should_201708113 = {
            "datetime": "2017-08-11 13:00:00",
            self.variable: 0.183,
            f"{self.variable}_flag": "G",
            f"{self.variable}_orig_flag": "M",
        }

    def test_metadata(self):
        """test reading the loaded metadata for a file (incl. static meta)"""

        assert self.file.check_metadata()
        assert not self.file.check_metadata(filter_meta_dict={"network": "WRONGNAME"})

        flag = self.file.check_metadata(
            self.variable,
            allowed_depth=Depth(0, 0.21),
            filter_meta_dict={
                "longitude": self.longitude,
                "variable": self.variable,
            },
        )
        assert flag == True

        flag = self.file.check_metadata("nonexistingvar", Depth(0, 0.1))
        assert flag == False

        flag = self.file.check_metadata(self.variable, Depth(98.0, 99.0))
        assert flag == False

        assert self.file.metadata["station"].val == self.station_name
        assert self.file.metadata["network"].val == self.network_name
        assert self.file.metadata["instrument"].depth.start == self.depth_from
        assert self.file.metadata["instrument"].depth.end == self.depth_to
        assert self.file.metadata["instrument"].val == self.instrument
        assert self.file.metadata["longitude"].val == self.longitude
        assert self.file.metadata["latitude"].val == self.latitude
        assert self.file.metadata["variable"].val == self.variable

        self.file.check_metadata(
            filter_meta_dict={
                "timerange_from": datetime(2017, 8, 10, 0),
                "timerange_to": datetime(2018, 8, 9, 8),
            }
        )

    def test_data(self):
        """test reading the actual data for a file"""
        # todo: why is sm column called "variable"?
        timestamp = datetime(2017, 8, 11, 13)

        data = self.file.read_data()
        data_is = data.loc[timestamp]

        assert data_is[self.variable] == self.data_should_201708113[self.variable]
        assert (
            data_is[f"{self.variable}_flag"]
            == self.data_should_201708113[f"{self.variable}_flag"]
        )
        assert (
            data_is[f"{self.variable}_orig_flag"]
            == self.data_should_201708113[f"{self.variable}_orig_flag"]
        )

    def test_metadata_for_depth(self):
        """Check finding best matching metadata for file"""
        bestmeta = self.file.read_metadata(best_meta_for_sensor=True)
        allmeta = self.file.read_metadata(best_meta_for_sensor=False)

        assert bestmeta == allmeta  # no vars with multiple depths

        assert "saturation" not in allmeta


class Test_DataFileCeopSepZipped(Test_DataFileCeopSepUnzipped):
    # same as test from unzipped archive, but uses zip directly. Runs all tests.
    def setUp(self) -> None:
        super(Test_DataFileCeopSepZipped, self).setUp()

        root = (
            testdata_path
            / "zip_archives"
            / "ceop"
            / "Data_seperate_files_20170810_20180809.zip"
        )
        filepath = Path(
            self.network_name,
            self.station_name,
            f"{self.network_name}_{self.network_name}_{self.station_name}_sm_{self.depth_from:.6f}_{self.depth_to:.6f}_{self.instrument}_20170810_20180809.stm",
        )

        self.file = DataFile(root, filepath)


class Test_DataFileHeaderValuesUnzipped(Test_DataFileCeopSepUnzipped):
    # same as for ceop sep format, but uses header values data
    def setUp(self) -> None:
        super(Test_DataFileHeaderValuesUnzipped, self).setUp()

        root = testdata_path / "Data_seperate_files_header_20170810_20180809"
        filepath = Path(
            self.network_name,
            self.station_name,
            f"{self.network_name}_{self.network_name}_{self.station_name}_sm_{self.depth_from:.6f}_{self.depth_to:.6f}_{self.instrument}_20170810_20180809.stm",
        )

        self.file = DataFile(root, filepath)


class Test_DataFileHeaderValuesZipped(Test_DataFileCeopSepUnzipped):
    # same as for ceop sep format, but uses header values data
    def setUp(self) -> None:
        super(Test_DataFileHeaderValuesZipped, self).setUp()

        root = (
            testdata_path
            / "zip_archives"
            / "header"
            / "Data_seperate_files_header_20170810_20180809.zip"
        )
        filepath = Path(
            self.network_name,
            self.station_name,
            f"{self.network_name}_{self.network_name}_{self.station_name}_sm_{self.depth_from:.6f}_{self.depth_to:.6f}_{self.instrument}_20170810_20180809.stm",
        )

        self.file = DataFile(root, filepath)


# todo: test from zip, _load_data
if __name__ == "__main__":
    unittest.main()
