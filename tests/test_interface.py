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

'''
tests for the ismn interface
Created on Thu Feb 26 12:36:30 2015

@author: Christoph Paulik

Updated on Dec 14, 2018

@author: Philip Buttinger philip.buttinger@geo.tuwien.ac.at
'''

import matplotlib
matplotlib.use('Agg')
from ismn import interface
import os
import sys
import datetime
import pytest
import numpy.testing as nptest
import ismn.metadata_collector as metadata_collector
import numpy as np


def test_min_max_obstime_getting():
    """
    test of the getting minimum and maxiumum observation time
    of a station
    """

    path_header_values = os.path.join(os.path.dirname(__file__),
                                      'test_data', 'format_header_values', 'SMOSMANIA')
    hv_interface = interface.ISMN_Interface(path_header_values)

    station = hv_interface.get_station('Narbonne')
    startd, endd = station.get_min_max_obs_timestamp()
    assert startd == datetime.datetime(2007, 1, 1, 1)
    assert endd == datetime.datetime(2007, 1, 31, 23)

    path_ceop_sep = os.path.join(os.path.dirname(__file__),
                                 'test_data', 'format_ceop_sep', 'SMOSMANIA')
    ceop_sep_interface = interface.ISMN_Interface(path_ceop_sep)

    station = ceop_sep_interface.get_station('Narbonne')
    startd, endd = station.get_min_max_obs_timestamp()
    assert startd == datetime.datetime(2007, 1, 1, 1)
    assert endd == datetime.datetime(2007, 1, 31, 23)


def test_min_max_obstime_networks():
    """
    test of the getting minimum and maxiumum observation time
    of several networks
    """

    path_header_values = os.path.join(os.path.dirname(__file__),
                                      'test_data', 'multinetwork', 'header_values')
    hv_interface = interface.ISMN_Interface(path_header_values)
    data = hv_interface.get_min_max_obs_timestamps(min_depth=0, max_depth=0.1)
    assert data.loc['MAQU']['end date'][
        0] == datetime.datetime(2010, 7, 31, 23)
    assert data.loc['MAQU']['end date'][
        1] == datetime.datetime(2010, 7, 31, 23)
    assert data.loc['MAQU']['start date'][
        1] == datetime.datetime(2008, 7, 1, 0)
    assert data.loc['SCAN']['start date'][
        1] == datetime.datetime(2007, 1, 1, 0)
    assert data.loc['SOILSCAPE']['start date'][
        1] == datetime.datetime(2012, 12, 14, 19)


def test_interface_network_init():
    """
    test limitation of interface to certain networks
    """

    path_header_values = os.path.join(os.path.dirname(__file__),
                                      'test_data', 'multinetwork', 'header_values')
    hv_interface = interface.ISMN_Interface(
        path_header_values, network=['SCAN'])
    assert hv_interface.list_networks().size == 1
    assert hv_interface.list_networks()[0] == 'SCAN'
    hv_interface = interface.ISMN_Interface(
        path_header_values, network=['SCAN', 'MAQU'])
    assert hv_interface.list_networks().size == 2


def test_find_nearest_station():
    """
    Test nearest neighbor search
    """
    path_header_values = os.path.join(os.path.dirname(__file__),
                                      'test_data', 'multinetwork', 'header_values')
    hv_interface = interface.ISMN_Interface(
        path_header_values, network=['SCAN'])
    station, distance = hv_interface.find_nearest_station(-90, 35, True)
    assert station.station == "AAMU-jtg"
    assert station.network == "SCAN"
    nptest.assert_almost_equal(distance, 316228.53147802927)


@pytest.mark.skipif(sys.version_info[0] == 3 and sys.version_info[1] == 4,
                    reason="Cartopy for python 3.4 does not support plotting of state boundaries.")
@pytest.mark.mpl_image_compare(tolerance=7)
def test_interface_plotting():
    """
    test plotting of networks
    """
    path_header_values = os.path.join(os.path.dirname(__file__),
                                      'test_data', 'multinetwork', 'header_values')
    hv_interface = interface.ISMN_Interface(
        path_header_values, network=['SCAN'])
    fig, axes = hv_interface.plot_station_locations()
    return fig


def test_station_order():
    """
    Test the station order returned by the metadata collector
    """
    path_header_values = os.path.join(os.path.dirname(__file__),
                                      'test_data', 'multinetwork', 'header_values')

    metadata = metadata_collector.collect_from_folder(path_header_values)

    filenames = []
    for m in metadata:
        filenames.append(m['filename'])

    sorted_filenames = sorted(filenames)

    assert sorted_filenames == filenames


def test_list_landcover_types():
    """
    Test available landcover classifications for dataset
    """
    path_to_ismn_data = os.path.join(os.path.dirname(__file__), 'test_data',
                                     'Data_seperate_files_header_20170810_20180809')
    ISMN_reader = interface.ISMN_Interface(path_to_ismn_data)
    lc = ISMN_reader.get_landcover_types()
    assert list(lc) == [130, 210]
    lc = ISMN_reader.get_landcover_types(landcover='landcover_insitu')
    assert lc == ['']

    path_to_ismn_data = os.path.join(os.path.dirname(__file__), 'test_data',
                                     'Data_seperate_files_20170810_20180809')
    ISMN_reader = interface.ISMN_Interface(path_to_ismn_data)
    lc = ISMN_reader.get_landcover_types()
    assert list(lc) == [130, 210]
    lc = ISMN_reader.get_landcover_types(landcover='landcover_insitu')
    assert lc == ['']


def test_list_climate_types():
    """
    Test available climate classifications for dataset
    """
    path_to_ismn_data = os.path.join(os.path.dirname(__file__), 'test_data',
                                     'Data_seperate_files_header_20170810_20180809')
    ISMN_reader = interface.ISMN_Interface(path_to_ismn_data)
    cl = ISMN_reader.get_climate_types()
    assert sorted(list(cl)) == sorted(['Cfa', 'ET'])
    cl = ISMN_reader.get_climate_types(climate='climate_insitu')
    assert list(cl) == []

    path_to_ismn_data = os.path.join(os.path.dirname(__file__), 'test_data',
                                     'Data_seperate_files_20170810_20180809')
    ISMN_reader = interface.ISMN_Interface(path_to_ismn_data)
    cl = ISMN_reader.get_climate_types()
    assert sorted(list(cl)) == sorted(['Cfa', 'ET'])
    cl = ISMN_reader.get_climate_types(climate='climate_insitu')
    assert list(cl) == []


def test_get_dataset_ids():
    """
    Test returned indeces from filtering
    """
    path_to_ismn_data = os.path.join(os.path.dirname(__file__), 'test_data',
                                     'Data_seperate_files_header_20170810_20180809')
    ISMN_reader = interface.ISMN_Interface(path_to_ismn_data)
    ids1 = ISMN_reader.get_dataset_ids(variable='soil moisture', min_depth=0, max_depth=1, landcover_2010=130)
    assert np.array_equal(np.array([0]), ids1)
    ids2 = ISMN_reader.get_dataset_ids(variable='soil moisture', min_depth=0, max_depth=1, climate='ET')
    assert np.array_equal(np.array([1]), ids2)
    ids3 = ISMN_reader.get_dataset_ids(variable='soil moisture', min_depth=0, max_depth=1)
    assert np.array_equal(np.array([0, 1]), ids3)
    ids4 = ISMN_reader.get_dataset_ids(variable='soil moisture', min_depth=0, max_depth=1, landcover_insitu='')
    assert np.array_equal(np.array([0, 1]), ids4)

    path_to_ismn_data = os.path.join(os.path.dirname(__file__), 'test_data',
                                     'Data_seperate_files_20170810_20180809')
    ISMN_reader = interface.ISMN_Interface(path_to_ismn_data)
    ids1 = ISMN_reader.get_dataset_ids(variable='soil moisture', min_depth=0, max_depth=1, landcover_2010=130)
    assert np.array_equal(np.array([0]), ids1)
    ids2 = ISMN_reader.get_dataset_ids(variable='soil moisture', min_depth=0, max_depth=1, climate='ET')
    assert np.array_equal(np.array([1]), ids2)
    ids3 = ISMN_reader.get_dataset_ids(variable='soil moisture', min_depth=0, max_depth=1)
    assert np.array_equal(np.array([0, 1]), ids3)
    ids4 = ISMN_reader.get_dataset_ids(variable='soil moisture', min_depth=0, max_depth=1, landcover_insitu='')
    assert np.array_equal(np.array([0, 1]), ids4)


