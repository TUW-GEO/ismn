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

"""
Created on Jul 31, 2013
@author: Christoph Paulik
Updated on Dec 14, 2018
@author: Philip Buttinger philip.buttinger@geo.tuwien.ac.at
"""

import os
import unittest
from ismn.filehandlers import DataFile
from pandas import Timestamp
import pandas as pd
import pytest


class TestReaders(unittest.TestCase):
    """
    Old readers, kept to test backward compatibility
    """

    def setUp(self) -> None:
        filename_format_header_values_root = os.path.join(
            os.path.dirname(__file__), "test_data", "format_header_values"
        )
        filename_format_header_values_filepath = os.path.join(
            "SMOSMANIA",
            "SMOSMANIA_SMOSMANIA_Narbonne_sm_0.050000_0.050000_ThetaProbe-ML2X_20070101_20070131.stm",
        )

        self.filehandler_header_values = DataFile(
            filename_format_header_values_root,
            filename_format_header_values_filepath,
        )

        filename_format_ceop_sep_root = os.path.join(
            os.path.dirname(__file__), "test_data", "format_ceop_sep"
        )
        filename_format_ceop_sep_filepath = os.path.join(
            "SMOSMANIA",
            "SMOSMANIA_SMOSMANIA_Narbonne_sm_0.050000_0.050000_ThetaProbe-ML2X_20070101_20070131.stm",
        )

        self.filehandler_ceop_sep = DataFile(
            filename_format_ceop_sep_root, filename_format_ceop_sep_filepath
        )

        filename_malformed_root = os.path.join(
            os.path.dirname(__file__), "test_data", "malformed"
        )
        filename_malformed_filepath = os.path.join("mal_formed_file.txt")

        with pytest.raises(IOError):
            DataFile(filename_malformed_root, filename_malformed_filepath)

        self.metadata_ref = {
            "network": "SMOSMANIA",
            "station": "Narbonne",
            "latitude": 43.15,
            "longitude": 2.9567,
            "elevation": 112.0,
            "variable": "soil_moisture",
            "timerange_from": Timestamp(2007, 1, 1, 1),
            "timerange_to": Timestamp(2007, 1, 31, 23),
            "instrument": "ThetaProbe-ML2X",
        }

        self.metadata_depth_from = 0.05
        self.metadata_depth_to = 0.05

        self.metadata_ref_ceop = dict(self.metadata_ref)
        self.metadata_ref_ceop["depth_from"] = ["multiple"]
        self.metadata_ref_ceop["depth_to"] = ["multiple"]
        self.metadata_ref_ceop["variable"] = ["ts", "sm"]
        self.metadata_ref_ceop["sensor"] = "n.s"

    def test_get_info_from_file(self):

        (
            header_elements,
            _,
            _,
            filename_elements,
        ) = self.filehandler_ceop_sep.get_elements_from_file()

        assert sorted(header_elements) == sorted(
            [
                "2007/01/01",
                "01:00",
                "2007/01/01",
                "01:00",
                "SMOSMANIA",
                "SMOSMANIA",
                "Narbonne",
                "43.15000",
                "2.95670",
                "112.00",
                "0.05",
                "0.05",
                "0.2140",
                "U",
                "M",
            ]
        )
        assert sorted(filename_elements) == sorted(
            [
                "SMOSMANIA",
                "SMOSMANIA",
                "Narbonne",
                "sm",
                "0.050000",
                "0.050000",
                "ThetaProbe-ML2X",
                "20070101",
                "20070131.stm",
            ]
        )

    def test_get_metadata_header_values(self):

        (
            metadata,
            depth,
        ) = self.filehandler_header_values.get_metadata_header_values()

        for key in metadata.keys():
            assert metadata[key].val == self.metadata_ref[key]
        assert metadata["variable"].depth.start == self.metadata_depth_from
        assert metadata["variable"].depth.end == self.metadata_depth_to

    def test_reader_format_header_values(self):
        filehandler = self.filehandler_header_values

        assert filehandler.metadata["network"].val == "SMOSMANIA"
        assert filehandler.metadata["station"].val == "Narbonne"
        assert filehandler.metadata["latitude"].val == 43.15
        assert filehandler.metadata["longitude"].val == 2.9567
        assert filehandler.metadata["elevation"].val == 112.0
        assert filehandler.metadata["variable"].val == "soil_moisture"
        assert filehandler.metadata["instrument"].depth.start == 0.05
        assert filehandler.metadata["instrument"].depth.end == 0.05
        assert filehandler.metadata["instrument"].val == "ThetaProbe-ML2X"

        data = filehandler.read_data()
        assert type(data) == pd.DataFrame
        assert data.index[7] == pd.Timestamp("2007-1-1 8:0:0")
        assert sorted(data.columns) == sorted(
            ["soil_moisture", "soil_moisture_flag", "soil_moisture_orig_flag"]
        )
        assert data["soil_moisture"].values[8] == 0.2135
        assert data["soil_moisture_flag"].values[8] == "U"
        assert data["soil_moisture_orig_flag"].values[8] == "M"

    def test_get_metadata_ceop_sep(self):

        filehandler = self.filehandler_ceop_sep

        metadata, depth = filehandler.get_metadata_ceop_sep()
        for key in metadata.keys():
            assert metadata[key].val == self.metadata_ref[key]
        assert metadata["variable"].depth.start == self.metadata_depth_from
        assert metadata["variable"].depth.end == self.metadata_depth_to

    def test_reader_format_ceop_sep(self):
        filehandler = self.filehandler_ceop_sep

        assert filehandler.metadata["network"].val == "SMOSMANIA"
        assert filehandler.metadata["station"].val == "Narbonne"
        assert filehandler.metadata["latitude"].val == 43.15
        assert filehandler.metadata["longitude"].val == 2.9567
        assert filehandler.metadata["elevation"].val == 112.0
        assert filehandler.metadata["variable"].val == "soil_moisture"
        assert filehandler.metadata["instrument"].depth.start == 0.05
        assert filehandler.metadata["instrument"].depth.end == 0.05
        assert filehandler.metadata["instrument"].val == "ThetaProbe-ML2X"

        data = filehandler.read_data()
        assert type(data) == pd.DataFrame
        assert data.index[7] == pd.Timestamp("2007-1-1 8:0:0")
        assert sorted(data.columns) == sorted(
            ["soil_moisture", "soil_moisture_flag", "soil_moisture_orig_flag"]
        )
        assert data["soil_moisture"].values[8] == 0.2135
        assert data["soil_moisture_flag"].values[8] == "U"
        assert data["soil_moisture_orig_flag"].values[347] == "M"

    def test_reader_get_format(self):
        assert self.filehandler_ceop_sep.file_type == "ceop_sep"
        assert self.filehandler_header_values.file_type == "header_values"

    def test_get_min_max_from_file(self):
        assert self.filehandler_ceop_sep.metadata["timerange_from"].val == Timestamp(
            2007, 1, 1, 1
        )
        assert self.filehandler_ceop_sep.metadata["timerange_to"].val == Timestamp(
            2007, 1, 31, 23
        )

        assert self.filehandler_header_values.metadata[
            "timerange_from"
        ].val == Timestamp(2007, 1, 1, 1)
        assert self.filehandler_ceop_sep.metadata["timerange_to"].val == Timestamp(
            2007, 1, 31, 23
        )


if __name__ == "__main__":
    unittest.main()
