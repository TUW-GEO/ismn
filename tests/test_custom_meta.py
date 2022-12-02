import os
from ismn.custom import CustomSensorMetadataCsv, CustomStationMetadataCsv
from ismn.meta import Depth
from ismn.interface import ISMN_Interface
import tempfile
import numpy as np

testdata_root = os.path.join(os.path.dirname(__file__), "test_data")

def test_build_custom_metadata_station():
    """
    Include custom metadata for the FR_Aqui - fraye station when building
    metadata.
    """

    with tempfile.TemporaryDirectory() as tmpdir:
        testdata = os.path.join(testdata_root, "Data_seperate_files_20170810_20180809")
        metadata_path = os.path.join(tmpdir, "python_metadata")
        csv_path = os.path.join(testdata_root, "custom_metadata",
                                "custom_stationmeta.csv")
        custom_meta_reader = [CustomStationMetadataCsv(csv_path)]

        ds = ISMN_Interface(
            testdata, meta_path=metadata_path, network=['FR_Aqui', 'COSMOS'],
            custom_meta_reader=custom_meta_reader
        )

        assert np.isnan(ds['COSMOS'][0][0].metadata['myvar1'].val)
        assert np.isnan(ds['COSMOS'][0][0].metadata['myvar2'].val)

        assert ds['FR_Aqui']['fraye'][0].metadata['myvar1'].val == 'lorem'
        assert ds['FR_Aqui']['fraye'][0].metadata['myvar1'].depth is None
        assert ds['FR_Aqui']['fraye'][0].metadata['myvar2'].val == 1.3
        assert ds['FR_Aqui']['fraye'][0].metadata['myvar3'].val == '2022-01-03'


def test_build_custom_metadata_sensor():
    """
    Include custom metadata for the FR_Aqui - fraye station when building
    metadata.
    """

    with tempfile.TemporaryDirectory() as tmpdir:
        testdata = os.path.join(testdata_root, "Data_seperate_files_20170810_20180809")
        metadata_path = os.path.join(tmpdir, "python_metadata")
        csv_path = os.path.join(testdata_root, "custom_metadata",
                                "custom_sensormeta.csv")
        custom_meta_reader = [CustomSensorMetadataCsv(csv_path)]

        ds = ISMN_Interface(
            testdata, meta_path=metadata_path, network=['FR_Aqui', 'COSMOS'],
            custom_meta_reader=custom_meta_reader
        )

        assert np.isnan(ds['COSMOS'][0][0].metadata['myvar1'].val)
        assert np.isnan(ds['COSMOS'][0][0].metadata['myvar2'].val)

        assert ds['FR_Aqui']['fraye'][0].metadata['myvar1'].val == 'lorem'
        assert ds['FR_Aqui']['fraye'][0].metadata['myvar1'].depth == Depth(0, 1)
        assert ds['FR_Aqui']['fraye'][0].metadata['myvar2'].val == 1.1
        assert ds['FR_Aqui']['fraye'][0].metadata['myvar2'].depth is None
        assert ds['FR_Aqui']['fraye'][0].metadata['myvar3'].val == '2022-01-01'
        assert ds['FR_Aqui']['fraye'][0].metadata['myvar3'].depth[1] == 1.0


if __name__ == '__main__':
    test_build_custom_metadata_sensor()
