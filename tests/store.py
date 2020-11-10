


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

    def test_get_sensors(self):
        """
        Test accessing sensor metadata and data.
        """
        for sensor in self.nwc.iter_sensors():
            assert sensor.instrument == 'Cosmic-ray-Probe'

    def test_get_nearest_station(self):
        """
        Test nearest station method.
        """
        station, dist = self.nwc.get_nearest_station(-97, 36)
        np.testing.assert_almost_equal(station.lon, -97.4878, 4)
        np.testing.assert_almost_equal(station.lat, 36.6054, 4)
        np.testing.assert_almost_equal(dist, 80201.85, 2)



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

    st = fc.iter_stations()
    assert len(st) == 2

    sen = fc.iter_sensors()
    assert len(sen) == 2

    assert len(fc.files) == 2

    data = fc.files[0].read_data()
    np.testing.assert_equal(data['soil_moisture'][0:4].values,
                            [0.141, 0.139, 0.139, 0.140])


if __name__ == '__main__':
    unittest.main()
