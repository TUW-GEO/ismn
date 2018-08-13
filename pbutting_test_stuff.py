
import os
import matplotlib.pyplot as plt
import random
import ismn.interface as ismn

path_to_ismn_data = os.path.join('/home/pbutting/shares/exchange/Staff/pbutting/Python_Projects/ismn/tests/test_data/'
                                 'Data_seperate_files_20170810_20180809')

# path_to_ismn_data = os.path.join('/home/pbutting/shares/exchange/Staff/pbutting/Python_Projects/ismn/tests/test_data/'
#                                  'Data_seperate_files_header_20170810_20180809')

ISMN_reader = ismn.ISMN_Interface(path_to_ismn_data)
networks = ISMN_reader.list_networks()

stations = ISMN_reader.list_stations()
station = random.choice(stations)
station_obj = ISMN_reader.get_station(station)


ids1 = ISMN_reader.get_dataset_ids(variable='soil moisture', min_depth=0, max_depth=1, landcover='Water')
ids2 = ISMN_reader.get_dataset_ids(variable='soil moisture', min_depth=0, max_depth=1, landcover='Water', climate='Polar - Tundra')

ts = ISMN_reader.read_ts(ids1[0])

lc = ISMN_reader.list_landcover_types()
clim = ISMN_reader.list_climate_types()

fig, ax = ISMN_reader.plot_station_locations()
plt.show()