# -*- coding: utf-8 -*-
from ismn.components import Depth

class MetaVar():
    def __init__(self, name, val, depth:Depth=None):
        self.name = name
        self.val = val
        self.depth = depth

class MetaData():
    def __init__(self, vars:list=None):
        if vars is None: # not a dict, because multiple meta vars with same name possible
            self.metadata = []
        else:
            self.metadata = vars

    def __getitem__(self, item):
        items = [v for v in self.metadata if v.name == item]
        if len(items) == 0:
            return None
        elif len(items) == 1:
            return items[0]
        else:
            return items

    def __repr__(self):
        names, depths = [], []
        for var in self.metadata:
            names.append(var.name)
            if var.depth is not None:
                depth = str(var.depth)
            else:
                depth = "no depth"
            depths.append(depth)
        return "\n".join([f"{name} ({depth})" for name, depth in zip(names, depths)])

    def keys(self) -> list :
        keys = []
        for var in self.metadata:
            keys.append(var.name)
        return keys

    @classmethod
    def from_dict(cls, data:dict):
        vars = []
        for k, v in data.items():
            vars.append(MetaVar(k, v))
        return cls(vars)
    
    def merge(self, other:'MetaData', inplace=False):
        if inplace:
            self.metadata += other.metadata
        else:
            vars = self.metadata + other.metadata
            return MetaData(vars)

    def add(self, name, val, depth:Depth=None):
        self.metadata.append(MetaVar(name, val, depth))

    def get_meta_for_depth(self, depth_from, depth_to=None):
        # go through vars and read those that have no depth or match best to passed d.
        # also get metadata withut depth assigned
        # sensor depth is assigned to sensor"
        # todo: implement method to return metadataonly for specific depth range
        # find a way to handle to select a best matching depth.
        NotImplementedError
        
    
if __name__ == '__main__':
    var1 = MetaVar('station', 'bla1')
    var2 = MetaVar('sand_fraction', 9000, Depth(0, 0.1))
    var3 = MetaVar('sand_fraction', 1, Depth(0.05, 0.1))
    var4 = MetaVar('sand_fraction', 1, Depth(0.1, 0.3))
    var5 = MetaVar('sand_fraction', 1, Depth(0.5, 1.))

    meta = MetaData([var1, var2, var3, var4, var5])

