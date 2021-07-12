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

import pytest

from ismn.components import NetworkCollection, Network, Station, Sensor, Depth
from pygeogrids.grids import BasicGrid

rpath = os.path.join(os.path.dirname(__file__), "test_data")

from ismn.filehandlers import DataFile


class NetCollTest(unittest.TestCase):
    def setUp(self) -> None:
        """
        Setup test data.
        """
        net1 = Network("Net1")
        net1.add_station("station_1_1", 0, 0, 0)
        net1.stations["station_1_1"].add_sensor(
            "sens_1_1_1", "var1", Depth(0.5, 1), None
        )
        net1.stations["station_1_1"].add_sensor("sens_1_1_2", "var1", Depth(1, 2), None)

        net2 = Network("Net2")
        net2.add_station("station_2_1", 1, 1, 1)
        net2.stations["station_2_1"].add_sensor(
            "sens_2_1_1", "var1", Depth(0.5, 1), None
        )
        net2.stations["station_2_1"].add_sensor("sens_2_1_2", "var1", Depth(1, 2), None)

        self.netcol = NetworkCollection([net1, net2])

    def test_grid(self):
        assert self.netcol.grid.gpi2lonlat(0) == (0, 0)
        assert self.netcol.grid.gpi2lonlat(1) == (1, 1)

    def test_station4idx(self):
        assert self.netcol.station4gpi(0).name == "station_1_1"
        assert self.netcol.station4gpi(1).name == "station_2_1"

    def test_get_nearest_station(self):
        assert self.netcol.get_nearest_station(0.1, 0.1)[0].name == "station_1_1"
        assert self.netcol.get_nearest_station(1, 1)[0].name == "station_2_1"

    def test_references(self):
        refs = self.netcol.export_citations(out_file=None)
        assert len(refs.keys()) == 2

class NetworkTest(unittest.TestCase):
    def setUp(self):
        """
        Setup test data.
        """
        self.network = Network("Network1")

    def test_attributes(self):
        """
        Test network attributes.
        """
        assert self.network.name == "Network1"

    def test_add_station(self):
        """
        Test adding stations.
        """
        name = "station1"
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
        self.network.add_station("station1", 0, 0, 0)
        self.network.add_station("station2", 0, 0, 0)

        assert self.network.n_stations == 2

        self.network.remove_station("station2")

        assert self.network.n_stations == 1

    def test_iter_stations(self):
        self.network.add_station("station1", 0, 0, 0)
        self.network.stations["station1"].add_sensor(
            "sens1", "var1", Depth(0.5, 1), None
        )
        self.network.stations["station1"].add_sensor("sens2", "var1", Depth(1, 2), None)

        for s in self.network.iter_stations(variable="var1", depth=Depth(0, 1)):
            # sensor 1 applies to conditions, therefore station is found.
            assert s.name == "station1"

        self.network.add_station("station2", 1, 1, 1)
        self.network.stations["station2"].add_sensor("sens1", "var1", Depth(0, 1), None)
        self.network.stations["station2"].add_sensor("sens2", "var2", Depth(1, 2), None)

        for s in self.network.iter_stations(variable="var1", depth=Depth(0, 0.5)):
            raise ValueError(
                "Found station but shouldn't ..."
            )  # this should never be reached

    def test_get_citation(self):
        refs = self.network.get_citations()
        assert isinstance(refs, list) and len(refs[0]) > 1


class StationTest(unittest.TestCase):
    def setUp(self):
        """
        Setup test data.
        """
        self.station = Station("station1", 0, 0, 0)

        name = "sensor1"
        d = Depth(0, 0.05)
        variable = "sm"
        se_name = "{}_{}_{:1.6f}_{:1.6f}".format(name, variable, d.start, d.end)

        self.station.add_sensor(name, variable, d, None)

        assert self.station.sensors[se_name].variable == "sm"

        name = "sensor2"
        d = Depth(0, 0.1)
        variable = "sm2"

        self.station.add_sensor(name, variable, d, None)

        assert self.station.n_sensors == 2

    def test_add_sensor(self):
        """
        Test adding sensors.
        """
        name = "sensor3"
        d = Depth(0, 0.1)
        variable = "sm"

        self.station.add_sensor(name, variable, d, None)
        assert self.station.n_sensors == 3

    def test_get_variables(self):
        """
        Test deriving the variables of all sensors at setion
        """
        assert self.station.get_variables() == ["sm", "sm2"]

    def test_get_depths(self):
        """
        Test deriving the depths of all sensors measuring a variable at station
        """
        assert self.station.get_depths("sm") == [Depth(0, 0.05)]

    def test_iter_sensors(self):
        """
        Test if sensor1 one is found when iteration over all sensors in 0-0.05 m.
        """
        for sen in self.station.iter_sensors(
            variable="sm",
            depth=Depth(0, 0.05),
            check_only_sensor_depth_from=False,
        ):
            assert sen.name == "sensor1_sm_0.000000_0.050000"

    def test_remove_sensor(self):
        """
        Test removing sensors.
        """
        assert self.station.n_sensors == 2

        self.station.remove_sensor("sensor2_sm2_0.000000_0.100000")

        assert self.station.n_sensors == 1

    def test_get_sensors(self):
        with pytest.deprecated_call():
            variable, depth = "soil_moisture", Depth(0, 0.05)
            sensors = self.station.get_sensors(
                variable, depth_from=depth.start, depth_to=depth.end
            )
            assert len(sensors) == len(
                [s for s in self.station.iter_sensors(variable=variable, depth=depth)]
            )


class SensorTest(unittest.TestCase):

    # todo: test reading sensor metadata?

    def setUp(self):
        """
        Setup test data.
        """
        instrument = "Cosmic-ray-Probe"
        d = Depth(0, 0.21)
        variable = "soil_moisture"

        root = os.path.join(rpath, "Data_seperate_files_20170810_20180809")
        subpath = os.path.join(
            "COSMOS",
            "Barrow-ARM",
            "COSMOS_COSMOS_Barrow-ARM_sm_0.000000_0.210000_Cosmic-ray-Probe_20170810_20180809.stm",
        )

        self.sensor = Sensor(
            instrument, variable, d, filehandler=DataFile(root, subpath)
        )

        name = "{}_{}_{:1.6f}_{:1.6f}".format(instrument, variable, d.start, d.end)

        assert self.sensor.name == name

    def test_sensor_attributes(self):
        """
        Test sensor attributes.
        """
        assert self.sensor.instrument == "Cosmic-ray-Probe"
        assert self.sensor.variable == "soil_moisture"
        assert self.sensor.depth == Depth(0, 0.21)

    def test_eval(self):
        assert self.sensor.eval("soil_moisture", Depth(0, 5))
        assert self.sensor.eval(
            "soil_moisture", (0, 0.21), check_only_sensor_depth_from=True
        )
        assert self.sensor.eval(
            "soil_moisture", Depth(0, 0.05), check_only_sensor_depth_from=True
        )
        # fails because of varname:
        assert not self.sensor.eval("wrongname", Depth(0, 0.21))
        # fails because of depth
        assert not self.sensor.eval(
            "soil_moisture", [0, 0.05], check_only_sensor_depth_from=False
        )

    def test_read_data(self):
        """Test reading the actual data"""
        data = self.sensor.read_data()
        assert data.index.size == 7059
        assert data.columns.size == 3


if __name__ == "__main__":
    unittest.main()
