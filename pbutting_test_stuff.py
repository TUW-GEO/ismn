
import os
import matplotlib.pyplot as plt
import random
import ismn.interface as ismn
import pandas as pd

# path_to_ismn_data = os.path.join('/home/pbutting/shares/exchange/Staff/pbutting/Python_Projects/ismn/tests/test_data/'
#                                  'Data_seperate_files_20170810_20180809')

# path_to_ismn_data = os.path.join('/home/pbutting/shares/exchange/Staff/pbutting/Python_Projects/ismn/tests/test_data/'
#                                  'Data_seperate_files_header_20170810_20180809')

path_to_ismn_data = os.path.join('/home/pbutting/shares/radar/Datapool_processed/ESA_CCI_SM/ISMN/4J89')

ISMN_reader = ismn.ISMN_Interface(path_to_ismn_data)
networks = ISMN_reader.list_networks()

stations = ISMN_reader.list_stations()
station = random.choice(stations)
station_obj = ISMN_reader.get_station(station)

ids1 = ISMN_reader.get_dataset_ids(variable='soil moisture', min_depth=0, max_depth=1, landcover='Water')
ids2 = ISMN_reader.get_dataset_ids(variable='soil moisture', min_depth=0, max_depth=1, landcover='Grassland', climate='Temperate - Without dry season - Hot Summer')
# simply all data
ids3 = ISMN_reader.get_dataset_ids(variable='air temperature', min_depth=0, max_depth=10)
ids4 = ISMN_reader.get_dataset_ids(variable='soil moisture', min_depth=0, max_depth=10)

ts = ISMN_reader.read_ts(ids3[0])

df_station_info = pd.DataFrame(ISMN_reader.metadata[ids2])

lc = ISMN_reader.get_landcover_types()
clim = ISMN_reader.get_climate_types()

#TODO: changed the colormap from Set1 to tab20 to get more colors for plotting the networks
fig, ax = ISMN_reader.plot_station_locations()

# ISMN_reader.metadata = ISMN_reader.metadata[ids1]
# ISMN_reader.plot_station_locations()

print(ISMN_reader.get_variables())

plt.show()