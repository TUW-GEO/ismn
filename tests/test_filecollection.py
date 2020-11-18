# -*- coding: utf-8 -*-

"""
Test file list functions, ie reading data for a filehander in the list,
and filtering the filelist.
"""

import os
import unittest
import shutil
from pathlib import Path
import pytest

import numpy as np
from ismn.file_collection import IsmnFileCollection

testdata_root = os.path.join(os.path.dirname(__file__), 'test_data')

def cleanup(metadata_path):
    # clean existing metadata
    if os.path.isdir(metadata_path):
        shutil.rmtree(metadata_path)

class Test_FileCollectionCeopSepUnzipped(unittest.TestCase):

    def setUp(self) -> None:
        testdata_path_unzipped = os.path.join(testdata_root,
            'Data_seperate_files_20170810_20180809')

        cleanup(os.path.join(testdata_path_unzipped, 'python_metadata'))

        self.coll = IsmnFileCollection(testdata_path_unzipped)

    def test_filelist(self):
        # cecks content of file collection

        cols_should = ['network', 'station', 'instrument', 'variable',
                       'depth_from', 'depth_to', 'root_path', 'file_path',
                       'timerange_from', 'timerange_to', 'filehandler']

        assert all([c in cols_should for c in self.coll.files.columns])

        assert self.coll.files.iloc[1]['root_path'] == self.coll.root.path
        assert self.coll.files.iloc[1]['file_path'].parts == \
               ("COSMOS", "Barrow-ARM",
                    "COSMOS_COSMOS_Barrow-ARM_sm_0.000000_0.210000_Cosmic-ray-Probe_20170810_20180809.stm")
        assert self.coll.files.iloc[1]['instrument'] == 'Cosmic-ray-Probe'
        assert self.coll.files.iloc[1]['variable'] == 'soil_moisture'
        assert self.coll.files.iloc[1]['depth_from'] == 0.0
        assert self.coll.files.iloc[1]['depth_to'] == 0.21

        # check some values that are in file list AND in metadata of filehandler
        assert self.coll.files.iloc[1].filehandler.metadata['station'].val == 'Barrow-ARM'
        assert self.coll.files.iloc[1].filehandler.metadata['instrument'].depth.start ==\
               self.coll.files.iloc[1]['depth_from']
        assert self.coll.files.iloc[1].filehandler.metadata['instrument'].depth.end ==\
               self.coll.files.iloc[1]['depth_to']

        # read data using a filehandler
        data = self.coll.files.iloc[1]['filehandler'].read_data()

        assert self.coll.files.iloc[1]['timerange_from'] == data.index[0]
        assert self.coll.files.iloc[1]['timerange_to'] == data.index[-1]

        assert self.coll.files.iloc[1]['variable']  in data.columns
        assert len(data.index) == 7059


    def test_from_filelist(self):
        # test alternative method to build collection, should get same result
        other = IsmnFileCollection.from_filelist(self.coll.files)

        for coll in self.coll.files.columns:
            if coll == 'filehandler': continue
            assert np.all(other.files[coll] == self.coll.files[coll]), \
            f"{coll} of filelists do not match"

    def test_filter_column(self):
        # test filtering the collection data frame directly (fast)
        filtered_station = self.coll.filter_col_val(col='station', vals=['Barrow-ARM'])
        assert np.all(filtered_station['station'] == 'Barrow-ARM')

        filtered_all = self.coll.filter_col_val(col='network', vals=['COSMOS'])
        # as only cosmos in testdata, get same list again.
        for col in self.coll.files.columns:
            if col == 'filehandler': continue
            assert np.all(filtered_all[col] == self.coll.files[col]), \
                f"{col} of filelists do not match"

    def test_filter_depth(self):
        # filter filelist for files in a certain depth range (fast)
        fil = self.coll.filter_depth(0, 1, return_index=False) # should get both
        assert all((fil.index == [0,1]) & (fil['network'] == ['COSMOS'] * 2))
        idx = self.coll.filter_depth(0, 0.19, return_index=True) # should get first
        assert (len(idx) == 1) & (idx[0] == 0)
        idx = self.coll.filter_depth(0, 0.01, return_index=True) # should get None
        assert len(idx) == 0

    def filter_metadata(self):
        # filter filelist by metadata in each file handler (slow)

        # filtering for network, station, etc can be done with both functions
        # but using the filter_col function is faster as it avoids iterating over
        # all filehandlers
        filtered = self.coll.filter_metadata({'station': 'Barrow-ARM'},
                                             return_index=False)
        assert np.all(filtered.index.values == \
            self.coll.filter_col_val('station', 'Barrow-ARM', return_index=True))

        # filtering for elements that are NOT in the filelist can only be done
        # with this function:
        filtered = self.coll.filter_metadata({'lc_2010': 210, 'climate_KG': 'ET'})
        assert len(filtered.index) == 1
        for col in [c for c in filtered.columns if c not in ['filehandler']]:
            assert filtered.loc[1, col] == self.coll.files.loc[1, col], \
                f"{col} of filelists do not match"

        idx = self.coll.filter_metadata({'lc_2010': 9999}, filelist=filtered,
                                        return_index=True)
        assert len(idx) == 0
        idx = self.coll.filter_metadata({'lc_2010': 210}, filelist=filtered,
                                        return_index=True)
        assert len(idx) == 1

        try:
            self.coll.filter_metadata({'wrong_name': 123})
            raise AssertionError("Expected error not raised.")
        except ValueError: # should raise ValueError
            pass

class Test_FileCollectionHeaderValuesUnzipped(Test_FileCollectionCeopSepUnzipped):
    # same tests as for ceop sep format,
    def setUp(self) -> None:
        testdata_path_unzipped = os.path.join(testdata_root,
            'Data_seperate_files_header_20170810_20180809')

        # clean existing metadata
        metadata_path = os.path.join(testdata_path_unzipped, 'python_metadata')
        cleanup(metadata_path)

        self.coll = IsmnFileCollection(testdata_path_unzipped)

@pytest.mark.zip
class Test_FileCollectionCeopSepZipped(Test_FileCollectionCeopSepUnzipped):
    # same tests as for ceop sep format,
    def setUp(self) -> None:
        testdata_path = os.path.join(testdata_root, 'zip_archives', 'ceop')
        testdata_zip_path = os.path.join(testdata_path,
            'Data_seperate_files_20170810_20180809.zip')

        # clean existing metadata
        metadata_path = os.path.join(testdata_path, 'python_metadata')
        cleanup(metadata_path)

        self.coll = IsmnFileCollection(testdata_zip_path)

@pytest.mark.zip
class Test_FileCollectionHeaderValuesZipped(Test_FileCollectionCeopSepUnzipped):
    # same tests as for ceop sep format,
    def setUp(self) -> None:
        testdata_path = os.path.join(testdata_root, 'zip_archives', 'header')
        testdata_zip_path = os.path.join(testdata_path,
                'Data_seperate_files_header_20170810_20180809.zip')

        # clean existing metadata
        metadata_path = os.path.join(testdata_path, 'python_metadata')
        cleanup(metadata_path)

        self.coll = IsmnFileCollection(testdata_zip_path)

if __name__ == '__main__':
    unittest.main()