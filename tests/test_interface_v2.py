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
import unittest

from ismn.interface_v2 import IsmnFile, IsmnFileCollection
from ismn.interface_v2 import Network, Station, Sensor, Depth


class NetworkTest(unittest.TestCase):

    def setUp(self):
        """
        Setup tests.
        """
        self.network = Network('test')

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
        Setup tests.
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
        Setup tests.
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


# class IsmnFileCollectionTest(unittest.TestCase):

#     def setUp(self):
#         """
#         Setup tests.
#         """
#         self.path = os.path.join(os.path.dirname(__file__),
#                                  'test_data', 'format_ceop_sep')

#     def test_metadata(self):
#         """
#         Test
#         """
#         fc = IsmnFileCollection(self.path)
#         fc.summary()


# class IsmnFileTest(unittest.TestCase):

#     def setUp(self):
#         """
#         Setup tests.
#         """
#         self.filename = os.path.join(
#             os.path.dirname(
#                 __file__), 'test_data', 'format_ceop_sep', 'SMOSMANIA',
#             'SMOSMANIA_SMOSMANIA_Narbonne_sm_0.050000_0.050000_ThetaProbe-ML2X_20070101_20070131.stm')
#     def test_metadata(self):
#         """
#         Test
#         """
#         IsmnFile(self.filename)
if __name__ == '__main__':
    unittest.main()
