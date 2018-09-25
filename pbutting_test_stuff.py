
import os
import matplotlib.pyplot as plt
import random
import ismn.interface as ismn
import pandas as pd
import numpy as np
import datetime

# TODO: check if data for soil fractions is useful or not

# path_to_ismn_data = os.path.join('/home/pbutting/shares/exchange/Staff/pbutting/Python_Projects/ismn/tests/test_data/'
#                                  'Data_seperate_files_20170810_20180809')

# path_to_ismn_data = os.path.join('/home/pbutting/shares/exchange/Staff/pbutting/Python_Projects/ismn/tests/test_data/'
#                                  'Data_seperate_files_header_20170810_20180809')

# path_to_ismn_data = os.path.join('/home/pbutting/shares/exchange/Staff/pbutting/Python_Projects/ismn/tests/test_data/'
#                                  'Data_seperate_files_20150820_20180820_2366439')

# path_to_ismn_data = os.path.join('/home/pbutting/shares/radar/Datapool_processed/ESA_CCI_SM/ISMN/4J89')

# use this networks
path_to_ismn_data = os.path.join('/home/pbutting/shares/radar/Datapool_processed/ESA_CCI_SM/ISMN/v20180830/')


# ISMN_reader = ismn.ISMN_Interface(path_to_ismn_data)
# #SCAN_network = ismn.ISMN_Interface(path_to_ismn_data, network=['SCAN'])
# df_station_info = pd.DataFrame(ISMN_reader.metadata)
#
# networks = ISMN_reader.list_networks()
#
# stations = ISMN_reader.list_stations()
# station = random.choice(stations)
# station_obj = ISMN_reader.get_station(station)
#
# ids1 = ISMN_reader.get_dataset_ids(variable='soil moisture', min_depth=0, max_depth=1, landcover_2010=10)
# ts_1 = ISMN_reader.read_ts(ids1[0])
# ids2 = ISMN_reader.get_dataset_ids(variable='soil moisture', min_depth=0, max_depth=1, landcover_2010=130, climate='Csa')
# # simply all data
# ids3 = ISMN_reader.get_dataset_ids(variable='air temperature', min_depth=0, max_depth=10)
# ids4 = ISMN_reader.get_dataset_ids(variable='soil moisture', min_depth=0, max_depth=10)
# ts_4 = ISMN_reader.read_ts(ids4[0])
#
# # TODO: a way to search for certain values in the soil fraction parameters --> need function for that?
# clay = ISMN_reader.metadata['clay_fraction']
# ids5 = np.where(['0.2m_0.2m' in t.dtype.fields and t['0.2m_0.2m'] < 10 if type(t) is np.ndarray else False for t in clay])[0]
#
# lc = ISMN_reader.get_landcover_types(landcover='landcover_2000')
# clim = ISMN_reader.get_climate_types(climate='climate')
# clim_insitu = ISMN_reader.get_climate_types(climate='climate_insitu')
#
# # TODO: changed the colormap from Set1 to tab20 to get more colors for plotting the networks
# # fig, ax = ISMN_reader.plot_station_locations()
#
# # ISMN_reader.metadata = ISMN_reader.metadata[ids1]
# # ISMN_reader.plot_station_locations()
#
# print(ISMN_reader.get_variables())
#
# plt.show()


dataset = ismn.readers.read_format_ceop('/home/pbutting/shares/exchange/Staff/pbutting/Python_Projects/ismn/tests/test_data/format_ceop/SMOSMANIA/SMOSMANIA_SMOSMANIA_NBN_20100304_20130801.stm')
assert dataset.network == 'SMOSMANIA'
assert dataset.station == 'Narbonne'
assert dataset.latitude == 43.15
assert dataset.longitude == 2.9567
assert dataset.elevation == 112.0
assert sorted(dataset.variable) == sorted(['sm', 'ts'])
assert sorted(dataset.depth_from) == sorted([0.05, 0.1, 0.2, 0.3])
assert sorted(dataset.depth_to) == sorted([0.05, 0.1, 0.2, 0.3])
assert dataset.sensor == 'n.s'
assert type(dataset.data) == pd.DataFrame
assert dataset.data.index[7] == (
    0.05, 0.05, datetime(2010, 10, 21, 9, 0, 0))
assert sorted(dataset.data.columns) == sorted(
    ['sm', 'sm_flag', 'ts', 'ts_flag'])
assert dataset.data['sm'].values[8] == 0.2227
assert dataset.data['sm_flag'].values[8] == 'U'
assert np.isnan(dataset.data.ix[0.3, 0.3]['ts'].values[6])
assert dataset.data.ix[0.3, 0.3]['ts_flag'].values[6] == 'M'