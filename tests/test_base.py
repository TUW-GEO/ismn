# -*- coding: utf-8 -*-
import os
from ismn.base import IsmnRoot
from pathlib import Path

testdata_path = os.path.join(os.path.dirname(__file__), "test_data")

testdata = Path(testdata_path)


def test_root_dir():
    root = IsmnRoot(testdata / "Data_seperate_files_20170810_20180809")
    assert root.root_dir.parent == Path(testdata_path)
    assert list(root.cont.keys()) == ["COSMOS", "FR_Aqui"]
    assert len(root.cont["COSMOS"]) == 2
    assert root.isopen

    assert "COSMOS/Barrow-ARM/COSMOS_COSMOS_Barrow-ARM_static_variables.csv" in root
    assert len(root.find_files("COSMOS/Barrow-ARM")) == 1

    root.close()

    assert root.isopen  # dir is always open


def test_root_zip():
    root = IsmnRoot(
        testdata / "zip_archives" / "ceop" / "Data_seperate_files_20170810_20180809.zip"
    )
    assert root.root_dir.parents[1] == Path(testdata_path)
    assert list(root.cont.keys()) == ["COSMOS"]
    assert len(root.cont["COSMOS"]) == 2
    assert root.isopen

    assert "COSMOS/Barrow-ARM/COSMOS_COSMOS_Barrow-ARM_static_variables.csv" in root
    assert len(root.find_files("COSMOS/Barrow-ARM")) == 1

    root.close()

    assert not root.isopen


if __name__ == "__main__":
    test_root_zip()
    test_root_dir()
