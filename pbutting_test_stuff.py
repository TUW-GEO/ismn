
import os
import matplotlib.pyplot as plt
import random
import ismn.interface as ismn
import pandas as pd
import numpy as np

# TODO: check if data for soil fractions is useful or not

# path_to_ismn_data = os.path.join('/home/pbutting/shares/exchange/Staff/pbutting/Python_Projects/ismn/tests/test_data/'
#                                  'Data_seperate_files_20170810_20180809')

# path_to_ismn_data = os.path.join('/home/pbutting/shares/exchange/Staff/pbutting/Python_Projects/ismn/tests/test_data/'
#                                  'Data_seperate_files_header_20170810_20180809')

# path_to_ismn_data = os.path.join('/home/pbutting/shares/exchange/Staff/pbutting/Python_Projects/ismn/tests/test_data/'
#                                  'Data_seperate_files_20150820_20180820_2366439')

path_to_ismn_data = os.path.join('/home/pbutting/shares/radar/Datapool_processed/ESA_CCI_SM/ISMN/4J89')


ISMN_reader = ismn.ISMN_Interface(path_to_ismn_data)
SCAN_network = ismn.ISMN_Interface(path_to_ismn_data, network=['SCAN'])
df_station_info = pd.DataFrame(ISMN_reader.metadata)

networks = ISMN_reader.list_networks()

stations = ISMN_reader.list_stations()
station = random.choice(stations)
station_obj = ISMN_reader.get_station(station)

ids1 = ISMN_reader.get_dataset_ids(variable='soil moisture', min_depth=0, max_depth=1, landcover='Water')
ids2 = ISMN_reader.get_dataset_ids(variable='soil moisture', min_depth=0, max_depth=1, landcover='Grassland', climate='Temperate - Without dry season - Hot Summer')
# simply all data
ids3 = ISMN_reader.get_dataset_ids(variable='air temperature', min_depth=0, max_depth=10)
ids4 = ISMN_reader.get_dataset_ids(variable='soil moisture', min_depth=0, max_depth=10)

ts = ISMN_reader.read_ts(ids4[0])

# TODO: a way to search for certain values in the soil fraction parameters --> need function for that?
clay = ISMN_reader.metadata['clay_fraction']
ids5 = np.where(['0.2m_0.2m' in t.dtype.fields and t['0.2m_0.2m'] < 10 if type(t) is np.ndarray else False for t in clay])[0]

lc = ISMN_reader.get_landcover_types()
clim = ISMN_reader.get_climate_types()

# TODO: changed the colormap from Set1 to tab20 to get more colors for plotting the networks
fig, ax = ISMN_reader.plot_station_locations()

# ISMN_reader.metadata = ISMN_reader.metadata[ids1]
# ISMN_reader.plot_station_locations()

print(ISMN_reader.get_variables())

plt.show()