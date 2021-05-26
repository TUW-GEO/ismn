# -*- coding: utf-8 -*-

"""
Test file list functions, ie reading data for a filehander in the list,
and filtering the filelist.
"""

import os
import unittest
import shutil
import pytest
from tempfile import TemporaryDirectory

from pathlib import Path
from ismn.filecollection import IsmnFileCollection

testdata_root = os.path.join(os.path.dirname(__file__), "test_data")


def cleanup(metadata_path):
    # clean existing metadata
    if os.path.isdir(metadata_path):
        shutil.rmtree(metadata_path)


class Test_FileCollectionCeopSepUnzipped(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """get_some_resource() is slow, to avoid calling it for each test use setUpClass()
        and store the result as class variable
        """
        super(Test_FileCollectionCeopSepUnzipped, cls).setUpClass()
        testdata_path_unzipped = os.path.join(
            testdata_root, "Data_seperate_files_20170810_20180809"
        )
        cleanup(os.path.join(testdata_path_unzipped, "python_metadata"))
        print("Setup from scratch")
        cls.coll = IsmnFileCollection.build_from_scratch(
            Path(testdata_path_unzipped), parallel=True
        )

    def setUp(self) -> None:
        pass

    def tearDown(self) -> None:
        self.coll.close()

    def test_filelist(self):
        # cecks content of file collection

        if (
            os.path.split(self.coll.root.path)[1]
            == "Data_seperate_files_20170810_20180809"
        ):
            # check that FR_Aqui is in data AND network name is network folder
            assert list(self.coll.filelist.keys()) == ["COSMOS", "FR_Aqui"]
            # check that files are associated to right network name
            files = [f for f in self.coll.iter_filehandlers(networks="FR_Aqui")]
            assert files != []
        else:
            assert list(self.coll.filelist.keys()) == ["COSMOS"]

        files = [f for f in self.coll.iter_filehandlers()]

        assert Path(files[1].file_path).parts == (
            "COSMOS",
            "Barrow-ARM",
            "COSMOS_COSMOS_Barrow-ARM_sm_0.000000_0.210000_Cosmic-ray-Probe_20170810_20180809.stm",
        )

        assert files[1].metadata["instrument"].val == "Cosmic-ray-Probe"
        assert files[1].metadata["variable"].val == "soil_moisture"
        assert files[1].metadata["variable"].depth.start == 0.0
        assert files[1].metadata["variable"].depth.end == 0.21

        # check some values that are in file list AND in metadata of filehandler
        assert files[1].metadata["station"].val == "Barrow-ARM"
        assert (
            files[1].metadata["instrument"].depth.start
            == files[1].metadata["variable"].depth.start
        )
        assert (
            files[1].metadata["instrument"].depth.end
            == files[1].metadata["variable"].depth.end
        )

        # read data using a filehandler
        data = files[1].read_data()

        assert files[1].metadata["timerange_from"].val == data.index[0]
        assert files[1].metadata["timerange_to"].val == data.index[-1]

        assert files[1].metadata["variable"].val in data.columns
        assert len(data.index) == 7059

    def test_from_csv(self):
        # test alternative method to build collection, should get same result
        print("Setup from csv")

        with TemporaryDirectory() as temp:
            self.coll.to_metadata_csv(os.path.join(temp, "meta.csv"))

            other = IsmnFileCollection.from_metadata_csv(
                self.coll.root.path, os.path.join(temp, "meta.csv")
            )

        for thisfile, otherfile in zip(
            self.coll.iter_filehandlers(), other.iter_filehandlers()
        ):
            assert thisfile.file_path == otherfile.file_path, "Paths dont match"
            assert thisfile.root.path == otherfile.root.path, "Paths dont match"
            assert thisfile.metadata == otherfile.metadata, "Meta dont match"


class Test_FileCollectionHeaderValuesUnzipped(Test_FileCollectionCeopSepUnzipped):
    # same tests as for ceop sep format,
    @classmethod
    def setUpClass(cls):
        super(Test_FileCollectionCeopSepUnzipped, cls).setUpClass()
        testdata_path_unzipped = os.path.join(
            testdata_root, "Data_seperate_files_header_20170810_20180809"
        )

        # clean existing metadata
        metadata_path = os.path.join(testdata_path_unzipped, "python_metadata")
        cleanup(metadata_path)

        cls.coll = IsmnFileCollection.build_from_scratch(testdata_path_unzipped)


@pytest.mark.zip
class Test_FileCollectionCeopSepZipped(Test_FileCollectionCeopSepUnzipped):
    # same tests as for ceop sep format,

    @classmethod
    def setUpClass(cls):
        testdata_path = os.path.join(testdata_root, "zip_archives", "ceop")
        testdata_zip_path = os.path.join(
            testdata_path, "Data_seperate_files_20170810_20180809.zip"
        )

        # clean existing metadata
        metadata_path = os.path.join(testdata_path, "python_metadata")
        cleanup(metadata_path)

        cls.coll = IsmnFileCollection.build_from_scratch(testdata_zip_path)


@pytest.mark.zip
class Test_FileCollectionHeaderValuesZipped(Test_FileCollectionCeopSepUnzipped):
    # same tests as for ceop sep format,
    @classmethod
    def setUpClass(cls):
        testdata_path = os.path.join(testdata_root, "zip_archives", "header")
        testdata_zip_path = os.path.join(
            testdata_path, "Data_seperate_files_header_20170810_20180809.zip"
        )
        # clean existing metadata
        metadata_path = os.path.join(testdata_path, "python_metadata")
        cleanup(metadata_path)

        cls.coll = IsmnFileCollection.build_from_scratch(testdata_zip_path)


if __name__ == "__main__":
    unittest.main()
