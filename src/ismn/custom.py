"""
Module that handles custom, additional information that can be assigned to the
ismn data by the user.
Sometimes it is convenient to have additional information on at a sensor,
station, or the surroundings, which is not directly provided by the ISMN,
assigned to the ISMN metadata.
This module contains a base class and implementations for certain metadata
formats that the ISMN_Interface class can then use to add additional values to
python_metadata during metadata collection.
"""

from abc import abstractmethod
from typing import Union
import numpy as np
from ismn.meta import MetaData, MetaVar, Depth
import pandas as pd


class CustomMetaReader:
    """
    Template class for a reader to assign additional metadata to ismn sensors.
    The `read_metadata` function must be implemented and return the metadata
    to add to a sensors either as MetaData object (which allows assigning
    depth information to metadata) or a dictionary (which will be converted
    later on to MetaData without depth information assigned)

    Metadata readers take the existing metadata from a sensor, and based on the
    information there they can extract other metadata.
    Can return Metadata objects or dicts (which are converted in ISMN package
    to metadata)

    Objects based on `CustomMetaReaders` can be passed to
    :class:`ismn.interface.Ismn_Interface`
    """

    @abstractmethod
    def read_metadata(self, meta) -> Union[MetaData, dict]:
        """
        Read metadata from additional sources (that are not provided directly
        by the ISMN). Uses available information for an ismn sensor for
        selecting the correct data (usually lat / lon of a sensor).

        Parameters
        ----------
        meta: MetaData
            Existing Metadata for a sensor, as collected from csv and .stm
            files. Contains for each sensor at least:
                Shared by all sensors at a station:
                    longitude, latitude, elevation, network, station,
                    lc_2010, lc_insitu, climate_KG, climate_insitu
                Sensor specific:
                    instrument (with depth_from and depth_to)
                    variable, clay_fraction, sand_fraction, organic_carbon,
                    silt_fraction (and depths of dataset layer they were
                    extracted from)

        Returns
        -------
        ancillary_meta: MetaData or dict
            Metadata collected by this reader. Dict also works but will be
            converted to MetaData without depths assigned later on.
            Metadata is then assigned to the sensor
        """
        ...


class CustomStationMetadataCsv(CustomMetaReader):
    """
    Allows passing (static) metadata for ISMN stations as a csv file.
    E.g. if the station specific variables provided by the ISMN are not enough.
    In this case that the metadata must be stored in a csv file with the
    following structure:

        network;station;<var1>;<var1>_depth_from;<var1>_depth_to;<var2>;...

    - where network and station refer to existing names in the metadata.
    - where <var1> etc. are the names of the custom metadata variables that are
    transferred into the python metadata
    - where <var1>_depth_from and <var1>_depth_to etc are the depths that
    are assigned to the metadata (if columns exist)
    """

    def __init__(self, station_meta_csv, fill_values=None, **kwargs):
        """
         Parameters
         ----------
         station_meta_csv: str
             Path to the csv file with the above described content
        fill_values: dict, optional (default: None)
             Values to use for a certain custom metadata variable, if no
             match is found.
         kwargs:
             Additional kwargs as passed to :func:`pandas.read_csv`
             To use a different separator than the default semicolon, use `sep`
        """
        if "sep" in kwargs:
            sep = kwargs.pop("sep")
        else:
            sep = ";"

        self.fill_values = dict() if fill_values is None else fill_values
        self.df = pd.read_csv(station_meta_csv, sep=sep, **kwargs)

    def _empty_var(self, varnames) -> list:
        """
        For all passed variable names, create an empty MetaVar if a fill value
        for the name is set in `self.fill_value`.
        """
        vars = []
        for var in varnames:
            if var in self.fill_values.keys() and not (
                var.endswith("_depth_from") or var.endswith("_depth_to")
            ):
                vars.append(MetaVar(var, self.fill_values[var]))
        return vars

    @staticmethod
    def _row2var(row: dict) -> list:
        """
        Extract name, value, depth from row.
        """
        vars = []

        for k, v in row.items():
            if k.endswith("_depth_from") or k.endswith("_depth_to"):
                continue

            if f"{k}_depth_from" in row:
                depth_from = row[f"{k}_depth_from"]
            else:
                depth_from = None
            if f"{k}_depth_to" in row:
                depth_to = row[f"{k}_depth_to"]
            else:
                depth_to = None

            if (depth_from is None) and (depth_to is None):
                depth = None
            else:
                if depth_from is None:
                    depth_from = -np.inf
                if depth_to is None:
                    depth_to = np.inf
                depth = Depth(depth_from, depth_to)

            vars.append(MetaVar(k, v, depth))

        return vars

    def read_metadata(self, meta: MetaData):
        """
        Match passed metadata entries to the csv file to find common stations
        for which the csv metadata is then added. The network and station
        names must match between csv file and previously collected metadata.

        Parameters
        ----------
        meta: MetaData
            Metadata to which the values from the csv file are added when
            the station and sensor name matches.

        Returns
        -------
        meta: dict
            Additional depth-independent metadata at the location

        """

        cond = (self.df["network"] == meta["network"].val) & (
            self.df["station"] == meta["station"].val
        )

        df = self.df[cond].set_index(["network", "station"])

        # drop potential duplicates, keep first
        df = df[~df.index.duplicated(keep="first")]

        vars = []

        if df.empty and (self.fill_values is not None):
            vars += self._empty_var(df.columns.values)
        else:
            for row in df.to_dict("records"):
                vars += self._row2var(row)

        return MetaData(vars)


class CustomSensorMetadataCsv(CustomStationMetadataCsv):
    """
    Allows passing metadata for ISMN sensors as a csv file. E.g. if the
    sensor specific variables provided by the ISMN are not enough.
    In this case that the metadata must be stored in a csv file with the
    following structure:

        network;station;instrument;variable;depth_from;depth_to;<var1>;<var1>_depth_from;<var1>_depth_to;<var2> ...

    where <var1> etc. are the names of the custom metadata variables that are
    transferred into the python metadata
    where <var1>_depth_from etc. are the
    """

    def read_metadata(self, meta: MetaData):
        """
        Match passed metadata entries to the csv file to find common sensors
        for which the csv metadata is then added.

        Parameters
        ----------
        meta: MetaData
            Metadata that the csv values are added to for sensors where
            the network, station, instrument, and instrument depths match.

        Returns
        -------
        meta: Metadata
            Additional depth-dependent metadata at the location
        """
        cond = (
            (self.df["network"] == meta["network"].val)
            & (self.df["station"] == meta["station"].val)
            & (self.df["instrument"] == meta["instrument"].val)
            & (self.df["variable"] == meta["variable"].val)
            & (self.df["depth_from"] == meta["instrument"].depth[0])
            & (self.df["depth_to"] == meta["instrument"].depth[1])
        )

        df = self.df[cond].set_index(
            ["network", "station", "instrument", "variable", "depth_from", "depth_to"]
        )

        # drop potential duplicates, keep first
        df = df[~df.index.duplicated(keep="first")]

        vars = []

        if df.empty and (self.fill_values is not None):
            vars += self._empty_var(df.columns.values)
        else:
            for row in df.to_dict("records"):
                vars += self._row2var(row)

        return MetaData(vars)
