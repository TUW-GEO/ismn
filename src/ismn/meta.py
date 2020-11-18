# -*- coding: utf-8 -*-
from ismn.components import Depth
import numpy as np
from typing import Optional, List, Any, Union

class MetaVar():
    """
    Meta Variable is a simple combination of a name, a value
    and a depth (optional)
    """
    def __init__(self,
                 name : str,
                 val : Any,
                 depth: Depth = None):
        """
        A named value that can be representative of a depth.

        Parameters
        ----------
        name : str
            Name of the variable
        val : Any
            Value of the variable
        depth : Depth, optional (default: None)
            Depth of the value
        """
        self.name = name
        self.val = val
        self.depth = depth

    def __str__(self):
        return f"{self.name} ({str(self.depth) if self.depth else 'no depth'}):" \
               f" {self.val}"

    def __iter__(self):
        yield self.name; yield self.val
        if self.depth is None:
            yield None; yield None
        else:
            yield self.depth.start; yield self.depth.end

    def __eq__(self, other):
        try:
            assert self.name == other.name
            assert self.val == other.val
            assert self.depth == other.depth
            return True
        except AssertionError:
            return False

class MetaData():
    """
    MetaData contains multiple MetaVars as a list (there can be multiple
    vars with the same name, e.g. for different depths)
    """
    def __init__(self,
                 vars: List[MetaVar] = None):
        """
        Parameters
        ----------
        vars : List[MetaVar]
            List of MetaVars that build the MetaData
        """
        if vars is None:
            self.metadata = []
            # not a dict, because multiple meta vars with same name possible
        else:
            self.metadata = vars

    def __iter__(self):
        for var in self.metadata:
            yield tuple(var)

    def __getitem__(self, item:Union[str,int]):
        # get all variables with the selected name
        if isinstance(item, int):
            return self.metadata[item]
        else:
            items = [v for v in self.metadata if v.name == item]
            if len(items) == 0:
                return None
            elif len(items) == 1:
                return items[0]
            else:
                return items

    def __str__(self):
        names, depths = [], []
        for var in self.metadata:
            names.append(var.name)
            if var.depth is not None:
                depth = str(var.depth)
            else:
                depth = "no depth"
            depths.append(depth)
        return ", ".join([f"{name} ({depth})" for name, depth in zip(names, depths)])

    def __contains__(self, item:Union[MetaVar, str]):
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

    def keys(self) -> list :
        # get only variable names
        keys = []
        for var in self.metadata:
            keys.append(var.name)
        return keys

    @classmethod
    def from_dict(cls, data:dict) -> 'MetaData':
        # Build Metadata from dict {name: (*args) ... }
        vars = []
        for k, v in data.items():
            v = np.atleast_1d(v)
            if len(v) == 3:
                if v[1] is None and v[2] is None:
                    vars.append(MetaVar(k, v[0]))
                else:
                    vars.append(MetaVar(k, v[0], Depth(float(v[1]), float(v[2]))))
            else:
                vars.append(MetaVar(k, *v))
        return cls(vars)
    
    def to_dict(self):
        d = {}
        for var in self:
            t = tuple(var)
            if t[0] in d.keys():
                raise ValueError(f'Cannot convert to dict with duplicate {t[0]}')
            d[t[0]] = t[1:]
        return d

    def merge(self, other:'MetaData', inplace=False) -> Optional['MetaData']:
        # Merge two metadata sets
        if inplace:
            self.metadata += other.metadata
        else:
            vars = self.metadata + other.metadata
            return MetaData(vars)

    def add(self, name, val, depth: Depth = None):
        """
        Create a new MetaVar and add it to this collection.

        Parameters
        ----------
        name : str
            Name of the variable
        val : Any
            Value of the variable
        depth : Depth, optional (default: None)
            A depth that is asssigned to the variable.
        """
        self.metadata.append(MetaVar(name, val, depth))

    def best_meta_for_depth(self, depth):
        """
        For meta variables that have a depth assigned, find the ones that match
        best (see func: perc_overlap()) to the passed depth.

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
            if isinstance(var, list):
                best_p = -np.inf
                best_var = var[0]
                for v in var:
                    p = depth.perc_overlap(v.depth)
                    if p > best_p:
                        best_p = p
                        best_var = v
                if depth.overlap(best_var.depth): # only add if best var overlaps
                    best_vars.append(best_var)
            else:
                if var.depth is None: # if var has no depth, use it
                    best_vars.append(var)
                elif depth.overlap(var.depth): # need overlap
                    best_vars.append(var)
                else:
                    pass # ignore var only if there is a depth but no overlap

        return MetaData(best_vars)

    def get_formatted(self, format='list'):
        """
        Return metadata for file in different formats such as list, dataframe,
        dict and structured array.

        Parameters
        ----------
        format : str, optional (default: 'list')
            Name of the format, one of 'struct', 'dict', 'pandas', 'list'

        Returns
        -------
        formatted_meta : Any
            Metadata in the selected format.
        """
        raise NotImplementedError
            #     meta = deepcopy(self.metadata)
            #
            #     depth = meta.pop('depth')
            #
            #     meta['depth_from'], meta['depth_to'] = depth.start, depth.end
            #
            #     meta['archive'] = self.root.path
            #     meta['filepath'] = self.file_path
            #
            #     if format.lower() == 'struct':
            #         return np.array([tuple(meta.values())],
            #                         [(k, object) for k in meta.keys()])
            #     elif format.lower() == 'dict':
            #         return meta
            #     elif format.lower() == 'pandas':
            #         return pd.Series(meta)
            #     elif format.lower() == 'list':
            #         return tuple(list(meta.values()))
            #     else:
            #         raise NotImplementedError(f"Format {format} is not (yet) implemented")

if __name__ == '__main__':
    var1 = MetaVar('station', 'bla1')
    var2 = MetaVar('sand_fraction', 9000, Depth(0, 0.1))
    var3 = MetaVar('sand_fraction', 1, Depth(0.05, 0.1))
    var4 = MetaVar('sand_fraction', 1, Depth(0.1, 0.3))
    var5 = MetaVar('sand_fraction', 1, Depth(0.5, 1.))

    meta = MetaData([var1, var2, var3, var4, var5])
    d = meta.best_meta_for_depth(Depth(0,0.5))
    dd = d.to_dict()

