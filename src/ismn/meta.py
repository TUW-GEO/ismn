# -*- coding: utf-8 -*-
# The MIT License (MIT)
#
# Copyright (c) 2021 TU Wien
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

from typing import Optional, List, Any, Union
import pandas as pd
from ismn.const import *


class Depth:
    """
    A class representing a depth range.
    For depth range start and end:
        0: surface
        >0: Below surface
        <0: Above surface

    Attributes
    ----------
    start: float
        Depth start. Upper boundary of a layer.
    end: float
        Depth end. Lower boundary of a layer.
    extend: float
        Layer range in metres.
    across0: bool
        Depth range across surface layer
    """

    def __init__(self, start, end):
        """
        Parameters
        ----------
        start : float
            Depth start. Upper boundary of a layer.
        end : float
            Depth end. Lower boundary of a layer.
        """

        self.start = float(start)
        self.end = float(end)

        self.extent = self.end - self.start

        if self.across0:
            if self.start > 0:
                raise DepthError("Start must be negative for Depths across 0")
        else:
            if abs(start) > abs(end):
                raise DepthError(
                    "Depth end can not be further from 0" " than depth start"
                )

    @property
    def is_profile(self) -> bool:
        return False if self.start == self.end else True

    @property
    def across0(self) -> bool:
        return True if (self.start * self.end) < 0 else False

    def __repr__(self):
        return f"{self.__class__.__name__}([{self.start}, {self.end}])"

    def __getitem__(self, item: int):
        return [self.start, self.end][item]

    def __str__(self):
        return f"{self.start} to {self.end} [m]"

    def __eq__(self, other):
        """
        Test if two Depths are equal.

        Parameters
        ----------
        other : Depth
            Different depth that is compared.

        Returns
        -------
        flag : bool
            True if both depths are equal, False otherwise.
        """
        if (self.start == other.start) and (self.end == other.end):
            flag = True
        else:
            flag = False

        return flag

    def __iter__(self):
        for d in [self.start, self.end]:
            yield d

    def __temp_pos_depths(self, other=None) -> ("Depth", Union["Depth", None]):
        # Create temporary depths that are shifted to positive

        shift = min(
            [self.end, other.end] + [self.start, other.start]
            if other is not None
            else []
        )

        if np.isfinite(shift) and (shift < 0):  # no neg depths
            # move both to pos range if necessary
            if ((self.start < 0) or (self.end < 0)) and (not self.across0):
                temp_d1 = Depth(self.end - shift, self.start - shift)
            else:
                temp_d1 = Depth(self.start - shift, self.end - shift)

            if other is not None:
                if ((other.start < 0) or (other.end < 0)) and (not other.across0):
                    temp_d2 = Depth(other.end - shift, other.start - shift)
                else:
                    temp_d2 = Depth(other.start - shift, other.end - shift)
            else:
                temp_d2 = None
        else:
            temp_d1 = Depth(self.start, self.end)
            temp_d2 = Depth(other.start, other.end) if other is not None else None

        return temp_d1, temp_d2

    def perc_overlap(self, other):
        """
        Estimate how much 2 depths correspond.
        - 1 means that the are the same
        - 0 means that they have an infinitely small correspondence
            (e.g. a single layer within a range, or 2 adjacent depths).
        - -1 means that they don't overlap.

        Parameters
        ----------
        other : Depth
            Second depth, overlap with this depth is calculated.

        Returns
        -------
        p : float
            Normalised overlap range
            <0 = no overlap, 0 = adjacent, >0 = overlap, 1 = equal
        """
        if self == other:  # same depths
            return 1
        else:
            # shift depths to pos ranges, flip start/end so that formulas work.
            temp_d1, temp_d2 = self.__temp_pos_depths(other)

            r = max([temp_d1.end, temp_d2.end]) - min([temp_d1.start, temp_d2.start])

            # Overlapping range normalised to the overall depth range r
            p_f = abs(temp_d1.start - temp_d2.start) / r
            p_t = abs(temp_d1.end - temp_d2.end) / r
            p = 1 - p_f - p_t

            p = round(p, 7)

            if p < 0:
                p = -1

        return p

    def overlap(self, other, return_perc=False):
        """
        Check if two depths overlap, (if the start of one depth is the same as
        the end of the other, they would also overlap),
        e.g. Depth(0, 0.1) and Depth(0.1, 0.2) do overlap.

        Parameters
        ----------
        other : Depth
            Other Depth
        return_perc : bool, optional (Default: False)
            Returns how much the depths overlap.
            See func: :func:`ismn.meta.Depth.perc_overlap`

        Returns
        -------
        overlap : bool
            True if Depths overlap
        perc_overlap: float, optional
            Normalised overlap.
        """
        other_start_encl = self.encloses(Depth(other.start, other.start))
        other_end_encl = self.encloses(Depth(other.end, other.end))

        this_start_encl = other.encloses(Depth(self.start, self.start))
        this_end_encl = other.encloses(Depth(self.end, self.end))

        overlap = any(
            [other_start_encl, other_end_encl, this_start_encl, this_end_encl]
        )

        if return_perc:
            return overlap, self.perc_overlap(other)
        else:
            return overlap

    def encloses(self, other):
        """
        Test if this Depth encloses other Depth.
        Reverse of :func:`ismn.meta.Depth.enclosed`.

        Parameters
        ----------
        other : Depth
            Check if other is enclosed by self.

        Returns
        -------
        flag : bool
            True if other depth is surrounded by given depth, False otherwise.
        """
        temp_d1, temp_d2 = self.__temp_pos_depths(other)
        if (temp_d1.start <= temp_d2.start) and (temp_d1.end >= temp_d2.end):
            flag = True
        else:
            flag = False

        return flag

    def enclosed(self, other):
        """
        Test if other Depth encloses this Depth.
        Reverse of :func:`ismn.meta.Depth.encloses`.

        Parameters
        ----------
        other : Depth
            Check if self is enclosed by other.

        Returns
        -------
        flag : bool
            True if other depth surrounds given depth, False otherwise.
        """
        temp_d1, temp_d2 = self.__temp_pos_depths(other)

        if (temp_d2.start <= temp_d1.start) and (temp_d2.end >= temp_d1.end):
            flag = True
        else:
            flag = False

        return flag


class MetaVar:
    """
    MetaVar is a simple combination of a name, a value
    and a depth range (optional).
    """

    def __init__(self, name: str, val: Any, depth: Depth = None):
        """
        A named value that can be representative of a depth range.

        Parameters
        ----------
        name : str
            Name of the variable
        val : Any
            Value of the variable, a number or string
        depth : Depth, optional (default: None)
            Depth range assigned to the value
        """
        self.name = name
        self.val = val
        self.depth = depth

    def __repr__(self):
        return (
            f"{self.__class__.__name__}([{self.name}, {self.val}, "
            f"{None.__repr__() if not self.depth else self.depth.__repr__()}])"
        )

    def __getitem__(self, item: int):
        return [self.name, self.val, self.depth][item]

    def __str__(self):
        d = str(self.depth) if self.depth else "no depth"
        return f"{self.name} ({d}): {self.val}"

    def __iter__(self):
        yield self.name
        yield self.val
        if self.depth is None:
            yield None
            yield None
        else:
            yield self.depth.start
            yield self.depth.end

    def __eq__(self, other: "MetaVar"):
        try:
            assert self.name == other.name
            assert (self.val == other.val) | np.all(pd.isna([self.val, other.val]))
            assert self.depth == other.depth
            return True
        except (AssertionError, AttributeError, TypeError):
            return False

    @property
    def empty(self) -> bool:
        # Check if Var has a valid value
        return pd.isnull(self.val)  # np.nan or None

    @classmethod
    def from_tuple(cls, args: tuple):
        """
        Create Metadata for a list of arguments.

        Parameters
        ----------
        args : tuple
            2 or 4 elements.
            2: name and value
            4: name, value, depth_from, depth_to
        """
        if len(args) == 2:
            return cls(*args)
        elif len(args) == 4:
            if pd.isnull(args[2]) or pd.isnull(args[3]):
                return cls(*args[:2])
            else:
                return cls(args[0], args[1], depth=Depth(args[2], args[3]))
        else:
            raise ValueError("Expected tuple of length 2 or 4.")


class MetaData:
    """
    MetaData contains multiple MetaVars as a list (there can be multiple
    vars with the same name, e.g. for different depths)
    """

    def __init__(self, vars: List[MetaVar] = None):
        """
        Parameters
        ----------
        vars : list[MetaVar]
            List of MetaVars that build the MetaData
        """
        if vars is None:
            self.metadata = []
            # not a dict, because multiple meta vars with same name possible
        else:
            self.metadata = vars

    def __iter__(self):
        for var in self.metadata:
            yield var

    def __repr__(self):
        return (
            f"{self.__class__.__name__}([\n"
            + ",\n".join(["  " + var.__repr__() for var in self.metadata])
            + "\n])"
        )

    def __getitem__(
        self, item: Union[str, int, list]
    ) -> Union["MetaData", MetaVar, None]:
        # get all variables with the selected name or at the selected index
        if not isinstance(item, list):
            if isinstance(item, int):
                return self.metadata[item]
            else:
                items = [v for v in self.metadata if v.name == item]
                if len(items) == 0:
                    return None
                elif len(items) == 1:
                    return items[0]
                else:
                    return MetaData(items)
        else:
            items = [v for v in self.metadata if v.name in item]
            return MetaData(items)

    def __contains__(self, item: Union[MetaVar, str]):
        if isinstance(item, MetaVar):
            for var in self.metadata:
                if var == item:
                    return True
        else:
            if item in self.keys():
                return True
        return False

    def __len__(self):
        return len(self.metadata)

    def __eq__(self, other):
        if len(self) != len(other):
            return False
        for v in self.metadata:
            if v not in other:
                return False
        for v in other.metadata:
            if v not in self:
                return False

        return True

    def keys(self) -> list:
        # get only variable names
        keys = []
        for var in self.metadata:
            keys.append(var.name)
        return keys

    def values(self) -> list:
        values = []
        for var in self.metadata:
            values.append(var.val)
        return values

    def to_dict(self):
        """
        Convert metadata to dictionary.

        Returns
        -------
        meta : dict
            Variable name as key, value and depth as values
        """
        d = {}
        for var in self.metadata:
            dat = tuple(var)
            name = dat[0]
            if name not in d.keys():
                d[name] = []
            d[name].append(dat[1:])

        return d

    def to_pd(self, transpose=False, dropna=True):
        """
        Convert metadata to a pandas DataFrame.

        Parameters
        ----------
        transpose : bool, optional (default: False)
            Organise variables in columns instead of rows.
        dropna : bool, optional (default: True)
            Drop NaNs, e.g. when no depth is assigned or a variable is empty
            the corresponding rows/cols will be excluded from the returned
            data frame.

        Returns
        -------
        df : pd.DataFrame
            Metadata collection as a data frame.
        """

        args = ["val", "depth_from", "depth_to"]

        var_names, values = [], []
        for var_name in np.unique(self.keys()):
            var = self[var_name]
            if isinstance(var, MetaVar):
                values.append(tuple(var)[1:])
                var_names.append(var_name)
            else:
                for v in var:
                    values.append(tuple(v)[1:])
                    var_names.append(var_name)

        values = list(sum(values, ()))

        index = pd.MultiIndex.from_product([var_names, args], names=["variable", "key"])

        df = pd.DataFrame(index=index, data=values).fillna(np.nan)
        df = df.rename(columns={0: "data"})

        if dropna:
            df.dropna(inplace=True)

        return df.loc[:, "data"] if not transpose else df.T

    def merge(self, other, inplace=False, exclude_empty=True):
        """
        Merge two or more metadata sets, i.e. take all variables from other(s)
        that are not in this metadata, and add them.

        Parameters
        ----------
        other: MetaData or list[MetaData]
            Other MetaData Collection or a list of MetaData, e.g. from multiple
            sensors.
        inplace: bool, optional (default: False)
            Replace self.metadata with the merged meteadata, if False then
            the merged metadata is returned
        exclude_empty : bool, optional (default: True)
            Variables where the value is NaN are ignored during merging.

        Returns
        -------
        merged_meta: MetaData or None
            The merged metadata (if inplace is False)
        """

        if isinstance(other, MetaData):
            other = [other]
        merged_meta = MetaData()

        for m in [self, *other]:
            for v in m.metadata:
                if (not v.empty if exclude_empty else True) and (v not in merged_meta):
                    merged_meta.add(v.name, v.val, v.depth)

        if inplace:
            self.metadata = merged_meta
        else:
            return merged_meta

    def add(self, name, val, depth=None):
        """
        Create a new MetaVar and add it to this collection.

        Parameters
        ----------
        name : str
            Name of the variable
        val : Any
            Value of the variable
        depth : Depth, optional (default: None)
            A depth that is assigned to the variable.
        """
        self.metadata.append(MetaVar(name, val, depth))

    def replace(self, name, val, depth=None):
        """
        Replace the value of a MetaVar in the initialized class

        Parameters
        ----------
        name: str
            Name of the MetaVar.
        val: Any
            New value of the MetaVar.
        depth : Depth, optional (default: None)
            New Depth of the variable.
        """
        Var = self.__getitem__(name)
        if not Var is None:
            self.metadata.remove(Var)
            self.metadata.append(MetaVar(name, val, depth))
        else:
            raise MetadataError("There is no MetaVar with name '{}'".format(name))

    def best_meta_for_depth(self, depth):
        """
        For meta variables that have a depth assigned, find the ones that match
        best (see :func:`ismn.depth.Depth.perc_overlap`) to the passed depth.

        Parameters
        ----------
        depth : Depth
            Reference depth, e.g. the depth of a sensor.

        Returns
        -------
        best_vars : MetaData
            A dict of variable names and a single variable for each name that
            was found to match best to the passed depth.
            Any variables that have a depth assigned which does not overlap
            with the passed depth are excluded here!
        """
        best_vars = []
        for varname in np.unique(self.keys()):
            var = self[varname]
            if isinstance(var, MetaData):
                best_p = -np.inf
                best_var = var[0]
                for v in var:
                    p = depth.perc_overlap(v.depth)
                    if p > best_p:
                        best_p = p
                        best_var = v
                if depth.overlap(best_var.depth):  # only add if best var overlaps
                    best_vars.append(best_var)
            else:
                if var.depth is None:  # if var has no depth, use it
                    best_vars.append(var)
                elif depth.overlap(var.depth):  # need overlap
                    best_vars.append(var)
                else:
                    pass  # ignore var only if there is a depth but no overlap

        return MetaData(best_vars)
