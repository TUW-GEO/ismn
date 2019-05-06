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
This module tests the ISMN interface.
"""

import os
import pytest
import unittest

import numpy as np

from ismn.class_interface import IsmnFile, IsmnFileCollection
from ismn.class_interface import Network, Station, Sensor, Depth
from ismn.class_interface import create_network_collection

rpath = os.path.join(os.path.dirname(__file__), 'test_data')


class NetworkTest(unittest.TestCase):

    def setUp(self):
        """
        Setup test data.
        """
        self.network = Network('Network1')

    def test_attributes(self):
        """
        Test network attributes.
        """
        assert self.network.name == 'Network1'

    def test_add_station(self):
        """
        Test adding stations.
        """
        name = 'station1'
        self.network.add_station(name, 0, 0, 0)

        assert self.network.stations[name].name == name
        assert self.network.stations[name].lon == 0
        assert self.network.stations[name].lat == 0
        assert self.network.stations[name].elev == 0

    def test_remove_station(self):
        """
        Test removing stations.
        """
        self.network.add_station('station1', 0, 0, 0)
        self.network.add_station('station2', 0, 0, 0)

        assert self.network.n_stations() == 2

        self.network.remove_station('station2')

        assert self.network.n_stations() == 1


class StationTest(unittest.TestCase):

    def setUp(self):
        """
        Setup test data.
        """
        self.station = Station('station1', 0, 0, 0)

    def test_add_sensor(self):
        """
        Test adding sensors.
        """
        name = 'sensor1'
        d = Depth(0, 0.05)
        variable = 'sm'
        se_name = '{}_{}_{:1.6f}_{:1.6f}'.format(
            name, variable, d.start, d.end)

        self.station.add_sensor(name, variable, d, None)

        assert self.station.sensors[se_name].variable == 'sm'

    def test_remove_sensor(self):
        """
        Test removing sensors.
        """
        name = 'sensor1'
        d = Depth(0, 0.05)
        variable = 'sm'
        se_name = '{}_{}_{:1.6f}_{:1.6f}'.format(
            name, variable, d.start, d.end)

        d2 = Depth(0, 0.10)

        self.station.add_sensor('sensor1', variable, d, None)
        self.station.add_sensor('sensor2', variable, d2, None)

        assert self.station.n_sensors() == 2

        self.station.remove_sensor(se_name)

        assert self.station.n_sensors() == 1


class SensorTest(unittest.TestCase):

    def setUp(self):
        """
        Setup test data.
        """
        instrument = 'sensor1'
        d = Depth(0, 0.05)
        variable = 'sm'
        name = '{}_{}_{:1.6f}_{:1.6f}'.format(
            instrument, variable, d.start, d.end)

        self.sensor = Sensor(name, instrument, variable, d)

    def test_sensor_attributes(self):
        """
        Test sensor attributes.
        """
        self.sensor.instrument == 'sensor1'
        self.sensor.variable = 'sm'
        self.sensor.depth == Depth(0, 0.05)


class DepthTest(unittest.TestCase):

    def setUp(self):
        """
        Setup test data.
        """
        self.d = Depth(0, 0.05)

    def test_attributes(self):
        """
        Test depth attributes.
        """
        assert self.d.start == 0
        assert self.d.end == 0.05

    def test_is_profile(self):
        """
        Test if depth represents a profile.
        """
        assert self.d.is_profile

    def test_equal(self):
        """
        Test depth equality.
        """
        other = Depth(0, 0.05)
        assert self.d == other

    def test_enclose(self):
        """
        Test if other depth encloses depth.
        """
        other = Depth(0, 0.05)
        assert self.d.enclose(other)

        other = Depth(0, 0.02)
        assert not self.d.enclose(other)


class NetworkCollectionTest(unittest.TestCase):

    def setUp(self):
        """
        Setup test data.
        """
        path = os.path.join(rpath, 'Data_seperate_files_20170810_20180809')
        self.nwc = create_network_collection(path)

    def test_attributes(self):
        """
        Test attributes.
        """
        assert len(self.nwc.networks) == 1


class IsmnFileTest(unittest.TestCase):

    def setUp(self):
        """
        Setup test data.
        """
        filename_ceop_sep = os.path.join(
            rpath, 'format_ceop_sep', 'SMOSMANIA',
            'SMOSMANIA_SMOSMANIA_Narbonne_sm_0.050000_0.050000_ThetaProbe-ML2X_20070101_20070131.stm')

        filename_header_values = os.path.join(
            rpath, 'format_header_values', 'SMOSMANIA',
            'SMOSMANIA_SMOSMANIA_Narbonne_sm_0.050000_0.050000_ThetaProbe-ML2X_20070101_20070131.stm')

        self.if_ceop_sep = IsmnFile(filename_ceop_sep)
        self.if_header_values = IsmnFile(filename_header_values)

    def test_metadata(self):
        """
        Test metadata.
        """
        assert self.if_ceop_sep.metadata == self.if_header_values.metadata

    def test_data(self):
        """
        Test data.
        """
        data_ceop_sep = self.if_ceop_sep.read_data()
        data_header_values = self.if_header_values.read_data()

        variables = ['soil_moisture', 'soil_moisture_flag']

        for variable in variables:
            np.testing.assert_array_equal(data_ceop_sep[variable].values,
                                          data_header_values[variable].values)


path_c = os.path.join(rpath, 'Data_seperate_files_20170810_20180809')
path_h = os.path.join(rpath, 'Data_seperate_files_header_20170810_20180809')
paths = [path_c, path_h]


@pytest.mark.parametrize("path", paths)
def test_IsmnFileCollection(path):
    """
    Test IsmnFileCollection class.
    """
    fc = IsmnFileCollection(path)

    nw = fc.get_networks()
    assert len(nw) == 1

    st = fc.get_stations()
    assert len(st) == 2

    sen = fc.get_sensors()
    assert len(sen) == 2

    assert len(fc.files) == 2

    data = fc.files[0].read_data()
    np.testing.assert_equal(data['soil_moisture'][0:4].values,
                            [0.141, 0.139, 0.139, 0.140])


if __name__ == '__main__':
    unittest.main()
