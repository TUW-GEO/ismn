# -*- coding: utf-8 -*-

import pandas as pd
from datetime import datetime
import numpy as np
from ismn.components import Depth
from ismn.meta import MetaVar, MetaData
from ismn.filehandlers import DataFile



metavars = ['saturation', 'lc_2010', 'network', 'timerange_from', 'timerange_to',
            'somevar1', 'somevar2', 'somevar3', 'somevar4', 'somevar5', 'somevar6']
args = ['val', 'depth_from', 'depth_to']

values1 =[1, 0, 0.1, 'lc', None, None, 'netname', None, None, datetime(2000,1,1), None, None,
          datetime(2000,1,1), None, None, 1,1,1, 1,1,1, 1,1,1, 1,1,1, 1,1,1, 1,1,1]

index = pd.MultiIndex.from_product([metavars, args], names=['name', 'meta_args'])


df = pd.DataFrame(index=index, data=values1).fillna(np.nan).T

dfs = [df, df.drop(columns=['network'])] + [df]*3000

data = pd.concat(dfs, axis=0, ignore_index=True)

row = data.iloc[0]

metavars = []

for idx, row in data.iterrows():
    print(idx)
    for metavar_name in row.index.get_level_values('name'):
        var = row.loc[metavar_name]
        depth_from, depth_to = var['depth_from'], var['depth_to']

        if np.all(np.isnan(np.array([depth_from, depth_to]))):
            depth = None
        else:
            depth = Depth(depth_from, depth_to)

        metavar = MetaVar(metavar_name, var['val'], depth)
        metavars.append(metavar)


    metadata = MetaData(metavars)
    filehandler = DataFile(r"C:\Temp\delete_me\ismn\scan",
             "SCAN/BeasleyLake/SCAN_SCAN_BeasleyLake_p_0.000000_0.000000_Pulse-Count_19780101_20191211.stm",
             load_metadata=False)
    filehandler.metadata = metadata