# -*- coding: utf-8 -*-

"""
Module that tests filehandler classes for reading ismn data.
"""

import os
import unittest


from ismn.components import Depth
from ismn.filehandlers import StaticMetaFile
from pathlib import Path

from ismn.const import CSV_META_TEMPLATE

testdata_path = Path(os.path.join(os.path.dirname(__file__), "test_data"))

# todo: test handling when some static vars are missing in csv file?


class Test_StaticMetaUnzipped(unittest.TestCase):
    def setUp(self) -> None:
        root = testdata_path / "Data_seperate_files_20170810_20180809"
        filepath = Path("COSMOS", "ARM-1", "COSMOS_COSMOS_ARM-1_static_variables.csv")
        self.file = StaticMetaFile(root, filepath)

    def test_read_metadata(self):
        meta = self.file.read_metadata()
        for k in CSV_META_TEMPLATE.keys():
            assert k in meta.keys()

        assert meta["lc_2010"].val == meta["lc_2005"].val == meta["lc_2000"].val == 130
        assert meta["climate_KG"].val == "Cfa"

        assert self.file.check_metadata(
            filter_meta_dict={
                "lc_2005": 130,
                "climate_KG": "Cfa",
                "sand_fraction": 36.0,
            }
        )

        sats = meta["saturation"]
        assert len(sats) == 2
        assert sats[0].val == 0.46
        assert sats[1].depth == Depth(0.3, 1.0)

        sands = meta["sand_fraction"]
        assert len(sands) == 2
        assert sands[1].val == 29.0
        assert sands[0].depth == Depth(0, 0.3)

        assert len(meta["silt_fraction"]) == len(meta["clay_fraction"]) == 2
        assert len(meta["organic_carbon"]) == 3
        assert meta["organic_carbon"][0].val == 1.27
        assert meta["organic_carbon"][0].depth == Depth(0, 0.3)
        assert meta["organic_carbon"][1].depth == Depth(0.3, 1)

        # todo: why is there a value for this depth even?
        assert meta["organic_carbon"][2].depth == Depth(-99.9, -99.9)
        assert meta["organic_carbon"][2].val == 0.59

    def test_metadata_for_depth(self):
        """Check finding best matching metadata for file"""
        allmeta = self.file.read_metadata()

        meta_for_depth = allmeta.best_meta_for_depth(Depth(0.0, 0.21))

        assert allmeta != meta_for_depth

        assert len([k for k in allmeta.keys() if k == "saturation"]) > 1
        assert len([k for k in meta_for_depth.keys() if k == "saturation"]) == 1

        assert (
            meta_for_depth["saturation"].depth.start
            == meta_for_depth["sand_fraction"].depth.start
            == 0.0
        )
        assert (
            meta_for_depth["saturation"].depth.end
            == meta_for_depth["sand_fraction"].depth.end
            == 0.3
        )

        assert allmeta["saturation"][0].depth == meta_for_depth["saturation"].depth

        assert allmeta["saturation"][1].depth.start == 0.3
        assert allmeta["saturation"][1].depth.end == 1.0


class Test_StaticMetaZipped(Test_StaticMetaUnzipped):
    # same as test from unzipped archive, but uses zip directly. Runs all tests.
    def setUp(self) -> None:
        super(Test_StaticMetaUnzipped, self).setUp()

        root = (
            testdata_path
            / "zip_archives"
            / "ceop"
            / "Data_seperate_files_20170810_20180809.zip"
        )
        filepath = Path("COSMOS", "ARM-1", "COSMOS_COSMOS_ARM-1_static_variables.csv")

        self.file = StaticMetaFile(root, filepath)


if __name__ == "__main__":
    unittest.main()
