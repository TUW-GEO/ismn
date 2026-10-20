# -*- coding: utf-8 -*-
"""
Tests for ISMN_Interface.find_sensors_by_tolerance_station_selection, the
station-level depth filter used by qa4sm when ISMN sensors are merged.

Test data: UDC_SMOS/Karolinenfeld holds 3 soil moisture sensors, all at
0.05 - 0.05 m, with (slightly) different measurement periods.
"""

import os
import unittest

import numpy as np
import pandas as pd

from ismn.interface import ISMN_Interface

testdata_root = os.path.join(os.path.dirname(__file__), "test_data")


def reshape_meta(metadata):
    """
    Same reshaping that qa4sm passes into the function under test:
    dict of {key: [(val, depth_from, depth_to)]} -> {key: val}, plus the
    depth of the instrument as separate entries.
    """
    reshaped = {}
    for key, value in metadata.items():
        meta_value = value[0][0]
        if isinstance(meta_value, pd.Timestamp):
            meta_value = meta_value.to_numpy()
        reshaped[key] = meta_value

    reshaped["instrument_depthfrom"] = metadata["instrument"][0][1]
    reshaped["instrument_depthto"] = metadata["instrument"][0][2]

    return reshaped


class Test_FindSensorsByToleranceStationSelection(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.testdata = os.path.join(
            testdata_root,
            "Data_separate_files_header_20080301_20100331_9550_1awu_20260817",
        )

    def setUp(self) -> None:
        self.ds = ISMN_Interface(self.testdata)
        ids = self.ds.get_dataset_ids("soil_moisture", -np.inf, np.inf,
                                      groupby="station")
        self.sensor_ids = ids["Karolinenfeld"]
        assert self.sensor_ids == [0, 1, 2]

    def tearDown(self) -> None:
        self.ds.close_files()

    def test_merge_all_sensors_of_station(self):
        """
        All 3 sensors sit in the (tolerance extended) 0 - 0.05 m layer, so the
        station is kept and the sensors are merged onto the lowest id.
        """
        primary, meta = self.ds.find_sensors_by_tolerance_station_selection(
            self.sensor_ids,
            depth_top=0.0,
            depth_bottom=0.1,
            top_tol=0.05,
            bottom_tol=0.05,
            reshape_meta=reshape_meta,
        )

        assert primary == 0
        assert meta["other_ids"] == [1, 2]
        assert meta["station"] == "Karolinenfeld"
        assert meta["instrument"] == "3 Averaged Sensors"
        assert meta["frm_class"] == "None"
        assert meta["frm_snr"] is None
        # depth span over all sensors of the station
        assert meta["instrument_depthfrom"] == 0.05
        assert meta["instrument_depthto"] == 0.05
        # merged series only covers the overlap of the sensor periods
        assert meta["timerange_from"] == np.datetime64("2008-05-18T00:00:00")
        assert meta["timerange_to"] == np.datetime64("2008-06-18T23:00:00")

    def test_station_discarded_outside_tolerance(self):
        """
        No sensor of the station reaches the 0.5 - 1.0 m layer -> station
        is dropped, or with keep_failed returned as a placeholder on its
        lowest sensor id.
        """
        kwargs = dict(depth_top=0.5, depth_bottom=1.0, top_tol=0.0,
                      bottom_tol=0.0, reshape_meta=reshape_meta)

        assert self.ds.find_sensors_by_tolerance_station_selection(
            self.sensor_ids, **kwargs) is None

        primary, meta = self.ds.find_sensors_by_tolerance_station_selection(
            self.sensor_ids, keep_failed=True, **kwargs)

        assert primary == 0
        assert meta["other_ids"] == []
        assert meta["station"] == "Karolinenfeld"
        assert meta["instrument"] == "Depth averaging failed"
        assert meta["frm_class"] == "None"
        assert meta["frm_snr"] is None
        assert meta["selection_failed"] is True
        assert meta["selection_reason"] == "depth"
        assert meta["instrument_depthfrom"] == 0.05
        assert meta["instrument_depthto"] == 0.05

    def test_keep_failed_no_overlap(self):
        """
        Sensors that pass the depth selection but never measure at the same
        time are discarded, or with keep_failed returned as a placeholder.
        """
        read_metadata = self.ds.read_metadata

        def read_metadata_shifted(idx, format="pandas"):
            # move sensor 2 to a period after the other two sensors
            meta = read_metadata(idx, format=format)
            if idx == 2:
                meta["timerange_from"] = [(pd.Timestamp("2009-01-01"),
                                           None, None)]
                meta["timerange_to"] = [(pd.Timestamp("2009-12-31"),
                                         None, None)]
            return meta

        self.ds.read_metadata = read_metadata_shifted
        kwargs = dict(depth_top=0.0, depth_bottom=0.05, top_tol=0.05,
                      bottom_tol=0.05, reshape_meta=reshape_meta)

        assert self.ds.find_sensors_by_tolerance_station_selection(
            self.sensor_ids, **kwargs) is None

        primary, meta = self.ds.find_sensors_by_tolerance_station_selection(
            self.sensor_ids, keep_failed=True, **kwargs)

        assert primary == 0
        assert meta["other_ids"] == []
        assert meta["instrument"] == "Depth averaging failed"
        assert meta["selection_failed"] is True
        assert meta["selection_reason"] == "no_overlap"

    def test_keep_failed_no_effect_on_kept_station(self):
        """
        A station that passes the selection is returned unchanged with
        keep_failed and is not flagged.
        """
        kwargs = dict(depth_top=0.0, depth_bottom=0.05, top_tol=0.05,
                      bottom_tol=0.05, reshape_meta=reshape_meta)
        primary, meta = self.ds.find_sensors_by_tolerance_station_selection(
            self.sensor_ids, **kwargs)
        primary_kf, meta_kf = \
            self.ds.find_sensors_by_tolerance_station_selection(
                self.sensor_ids, keep_failed=True, **kwargs)

        assert primary_kf == primary
        assert meta_kf["other_ids"] == meta["other_ids"] == [1, 2]
        assert meta_kf["instrument"] == "3 Averaged Sensors"
        assert "selection_failed" not in meta_kf

    def test_single_point_measurement_discarded(self):
        """
        A station with a single point-measurement sensor (depth_from ==
        depth_to) is discarded, or with keep_failed returned as a placeholder.
        """
        kwargs = dict(depth_top=0.05, depth_bottom=0.05, top_tol=0.0,
                      bottom_tol=0.0, reshape_meta=reshape_meta)

        assert self.ds.find_sensors_by_tolerance_station_selection(
            [2], **kwargs) is None

        primary, meta = self.ds.find_sensors_by_tolerance_station_selection(
            [2], keep_failed=True, **kwargs)

        assert primary == 2
        assert meta["other_ids"] == []
        assert meta["instrument"] == "Depth averaging failed"
        assert meta["selection_failed"] is True
        assert meta["selection_reason"] == "point_measurement"

    def test_single_sensor_metadata_untouched(self):
        """
        A station with a single matching sensor is kept as is, no merging
        metadata is written.
        """
        read_metadata = self.ds.read_metadata

        def read_metadata_layer(idx, format="pandas"):
            # turn the 0.05 m point measurement of sensor 2 into a
            # 0.0 - 0.05 m layer measurement
            meta = read_metadata(idx, format=format)
            if idx == 2:
                meta["variable"] = [("soil_moisture", 0.0, 0.05)]
            return meta

        self.ds.read_metadata = read_metadata_layer

        primary, meta = self.ds.find_sensors_by_tolerance_station_selection(
            [2],
            depth_top=0.0,
            depth_bottom=0.05,
            top_tol=0.0,
            bottom_tol=0.0,
            reshape_meta=reshape_meta,
        )

        assert primary == 2
        assert meta["other_ids"] == []
        assert meta["instrument"] == "EC-ET_1_5"
        assert meta["instrument_depthfrom"] == 0.0
        assert meta["instrument_depthto"] == 0.05
        # original period of that sensor, no intersection applied
        assert meta["timerange_from"] == np.datetime64("2008-05-17T23:00:00")
        assert meta["timerange_to"] == np.datetime64("2008-06-18T23:00:00")


if __name__ == "__main__":
    unittest.main()
