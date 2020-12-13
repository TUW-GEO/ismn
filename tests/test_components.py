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
This module tests the ISMN components.
"""

import os
import unittest

from ismn.components import Network, Station, Sensor, Depth
from pygeogrids.grids import BasicGrid
from ismn.tables import DepthError

rpath = os.path.join(os.path.dirname(__file__), 'test_data')

from ismn.filehandlers import DataFile


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

        assert self.network.grid == BasicGrid([0], [0])

    def test_remove_station(self):
        """
        Test removing stations.
        """
        self.network.add_station('station1', 0, 0, 0)
        self.network.add_station('station2', 0, 0, 0)

        assert self.network.n_stations() == 2

        self.network.remove_station('station2')

        assert self.network.n_stations() == 1

    def test_iter_stations(self):
        self.network.add_station('station1', 0, 0, 0)
        self.network.stations['station1'].add_sensor('sens1', 'var1', Depth(0.5, 1), None)
        self.network.stations['station1'].add_sensor('sens2', 'var1', Depth(1 ,2), None)

        for s in self.network.iter_stations('var1', Depth(0 ,1)):
            # sensor 1 applies to conditions, therefore station is found.
            assert s.name == 'station1'

        self.network.add_station('station2', 1, 1, 1)
        self.network.stations['station2'].add_sensor('sens1', 'var1', Depth(0 ,1), None)
        self.network.stations['station2'].add_sensor('sens2', 'var2', Depth(1 ,2), None)

        for s in self.network.iter_stations('var1', Depth(0 ,0.5)):
            raise ValueError("Found station but shouldn't ...") # this should never be reached

class StationTest(unittest.TestCase):

    def setUp(self):
        """
        Setup test data.
        """
        self.station = Station('station1', 0, 0, 0)

        name = 'sensor1'
        d = Depth(0, 0.05)
        variable = 'sm'
        se_name = '{}_{}_{:1.6f}_{:1.6f}'.format(
            name, variable, d.start, d.end)

        self.station.add_sensor(name, variable, d, None)

        assert self.station.sensors[se_name].variable == 'sm'

        name = 'sensor2'
        d = Depth(0, 0.1)
        variable = 'sm2'

        self.station.add_sensor(name, variable, d, None)

        assert self.station.n_sensors() == 2

    def test_add_sensor(self):
        """
        Test adding sensors.
        """
        name = 'sensor3'
        d = Depth(0, 0.1)
        variable = 'sm'

        self.station.add_sensor(name, variable, d, None)
        assert self.station.n_sensors() == 3

    def test_get_variables(self):
        """
        Test deriving the variables of all sensors at setion
        """
        assert self.station.get_variables() == ['sm', 'sm2']

    def test_get_depths(self):
        """
        Test deriving the depths of all sensors measuring a variable at station
        """
        assert self.station.get_depths('sm') == [Depth(0,0.05)]

    def test_iter_sensors(self):
        """
        Test if sensor1 one is found when iteration over all sensors in 0-0.05 m.
        """
        for sen in self.station.iter_sensors('sm', Depth(0, 0.05),
                                             check_only_sensor_depth_from=False):
            assert sen.name == 'sensor1_sm_0.000000_0.050000'

    def test_remove_sensor(self):
        """
        Test removing sensors.
        """
        assert self.station.n_sensors() == 2

        self.station.remove_sensor('sensor2_sm2_0.000000_0.100000')

        assert self.station.n_sensors() == 1

class SensorTest(unittest.TestCase):

    # todo: test reading sensor metadata?

    def setUp(self):
        """
        Setup test data.
        """
        instrument = 'Cosmic-ray-Probe'
        d = Depth(0, 0.21)
        variable = 'soil_moisture'
        name = '{}_{}_{:1.6f}_{:1.6f}'.format(
            instrument, variable, d.start, d.end)

        root = os.path.join(
            rpath, "Data_seperate_files_20170810_20180809")
        subpath = os.path.join("COSMOS", "Barrow-ARM",
            "COSMOS_COSMOS_Barrow-ARM_sm_0.000000_0.210000_Cosmic-ray-Probe_20170810_20180809.stm")

        self.sensor = Sensor(instrument, variable, d,
                             filehandler=DataFile(root, subpath))

    def test_sensor_attributes(self):
        """
        Test sensor attributes.
        """
        assert self.sensor.instrument == 'Cosmic-ray-Probe'
        assert self.sensor.variable == 'soil_moisture'
        assert self.sensor.depth == Depth(0, 0.21)

    def test_eval(self):
        assert self.sensor.eval('soil_moisture', Depth(0,5))
        assert self.sensor.eval('soil_moisture', Depth(0,0.21),
                                check_only_sensor_depth_from=True)
        assert self.sensor.eval('soil_moisture', Depth(0,0.05),
                                check_only_sensor_depth_from=True)
        # fails because of varname:
        assert not self.sensor.eval('wrongname', Depth(0,0.21))
        # fails because of depth
        assert not self.sensor.eval('soil_moisture', Depth(0,0.05),
                                    check_only_sensor_depth_from=False)

        # test based on metadata
        assert self.sensor.eval('soil_moisture', Depth(0,1),
                                filter_meta_dict={'lc_2010': 210, 'climate_KG': 'ET'})
        assert not self.sensor.eval('soil_moisture', Depth(0,1),
                                    filter_meta_dict={'lc_2010': 999})

    def test_read_data(self):
        """ Test reading the actual data """
        data = self.sensor.read_data()
        assert data.index.size == 7059
        assert data.columns.size == 3


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

    def test_invalid(self):
        try:
            Depth(0.5,0.1)
            raise AssertionError
        except DepthError:
            pass
        try:
            Depth(-0.5, -0.1)
            raise AssertionError
        except DepthError:
            pass

    def test_enclose(self):
        """
        Test if other depth encloses depth.
        """
        other = Depth(0, 0.05)
        assert self.d.encloses(other)
        assert other.enclosed(self.d)

        other = Depth(0, 0.1)
        assert not self.d.encloses(other)
        assert not other.enclosed(self.d)

        other = Depth(-0.1,-0.2)
        assert not self.d.encloses(other)
        assert not self.d.enclosed(other)


def test_perc_overlap(self):
        other = Depth(0.05, 0.1)
        assert self.d.overlap(other) == True
        assert self.d.perc_overlap(other) == 0.

        other = Depth(0.03, 0.05)
        assert self.d.overlap(other) == True
        assert self.d.perc_overlap(other) == round(0.02/0.05, 7)

        other = Depth(0,0.05)
        assert self.d.overlap(other) == True
        assert self.d.perc_overlap(other) == 1.

        other = Depth(-0.01, -0.05)
        assert self.d.overlap(other) == False
        assert self.d.perc_overlap(other) == -1

        other = Depth(0.01, -0.01)
        assert self.d.overlap(other) == True
        assert self.d.perc_overlap(other) == round(0.01/0.06, 7)

if __name__ == '__main__':
    unittest.main()
