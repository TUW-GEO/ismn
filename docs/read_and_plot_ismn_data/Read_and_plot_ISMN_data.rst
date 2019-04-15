
.. code:: python

    from ismn.interface import ISMN_Interface
    import matplotlib.pyplot as plt
    import random
    %matplotlib inline

.. code:: python

    path_to_ismn_data = "/data/CGLS/SQE_2016/cgls-validationreports/ISMN/raw/data_20160101_20161231/"

.. code:: python

    #initialize interface, this can take up to a few minutes the first
    #time, since all metadata has to be collected
    
    ISMN_reader = ISMN_Interface(path_to_ismn_data)
    
    #plot available station on a map
    fig, ax = ISMN_reader.plot_station_locations()
    plt.show()




.. image:: read_and_plot_ismn_data/output_2_0.png


Next we explore the available networks and stations and select a random
network and station to plot data from

.. code:: python

    networks = ISMN_reader.list_networks()
    print "Available Networks:"
    print networks


.. parsed-literal::

    Available Networks:
    ['BIEBRZA_S-1' 'COSMOS' 'FMI' 'LAB-net' 'PBO_H2O' 'REMEDHUS' 'RSMN' 'SCAN'
     'SNOTEL' 'SOILSCAPE' 'TERENO' 'USCRN' 'WEGENERNET']


.. code:: python

    network = random.choice(networks)
    stations = ISMN_reader.list_stations(network = network)
    print "Available Stations in Network %s"%network
    print stations



.. parsed-literal::

    Available Stations in Network SNOTEL
    ['AGUA_CANYON' 'ANCHOR_RIVER_DIVIDE' 'ANNIE_SPRINGS' 'ARAPAHO_RIDGE'
     'ATIGUN_PASS' 'ATLANTA_SUMMIT' 'BAKER_BUTTE_SMT' 'BALDY' 'BANNER_SUMMIT'
     'BEAR_CREEK' 'BEAR_RIVER_RS' 'BEAVER_DAMS' 'BEAVER_DIVIDE' 'BEAVER_PASS'
     'BEN_LOMOND_PEAK' 'BEN_LOMOND_TRAIL' 'BERRY_CREEK' 'BERTHOUD_SUMMIT'
     'BIG_BEND' 'BIG_CREEK_SUM' 'BIG_FLAT' 'BIG_GOOSE' 'BIG_MEADOW'
     'BIG_SANDY_OPENING' 'BILLIE_CREEK_DIVIDE' 'BIRD_CREEK' 'BLACKHALL_MTN'
     'BLACKS_FORK_JCT' 'BLACKTAIL_MTN' 'BLACK_BEAR' 'BLACK_FLAT-U.M._CK'
     'BLACK_PINE' 'BLUE_LAKES' 'BOGUS_BASIN' 'BONE_SPRINGS_DIV' 'BOURNE'
     'BOX_CREEK' 'BRIAN_HEAD' 'BRIGHTON' 'BRISTLECONE_TRAIL' 'BROWN_DUCK'
     'BROWN_TOP' 'BRUMLEY' 'BUCKBOARD_FLAT' 'BUCKINGHORSE' 'BUCKSKIN_JOE'
     'BUCKSKIN_LOWER' 'BUCK_FLAT' 'BUCK_PASTURE' 'BUG_LAKE' 'BURNSIDE_LAKE'
     'BURNT_MOUNTAIN' 'BURRO_MOUNTAIN' 'BURTS_MILLER_RANCH' 'BUTTE'
     'Baker_Butte' 'Bar_M' 'Bevans_Cabin' 'Black_Mesa' 'CAMP_JACKSON'
     'CARSON_PASS' 'CASCADE_#2' 'CASCADE_MOUNTAIN' 'CASTLE_CREEK'
     'CASTLE_VALLEY' 'CAVE_MOUNTAIN' 'CAYUSE_PASS' 'CHALK_CREEK_#1'
     'CHALK_CREEK_#2' 'CHEMULT_ALTERNATE' 'CHEPETA' 'CHOCOLATE_GULCH'
     'CINNABAR_PARK' 'CLACKAMAS_LAKE' 'CLAYTON_SPRINGS' 'CLEAR_CREEK_#1'
     'CLEAR_CREEK_#2' 'CLEAR_LAKE' 'CLOVER_MEADOW' 'COCHETOPA_PASS' 'COLDFOOT'
     'CORRAL_CANYON' 'CRAB_CREEK' 'CRATER_MEADOWS' 'CROW_CREEK' 'CSS_LAB'
     'CULEBRA_#2' 'CURRANT_CREEK' 'Chalender' 'Columbia_Basin'
     'Copper_Mountain' 'Corduroy_Flat' 'Corral' 'DANIELS-STRAWBERRY'
     'DIAMOND_PEAK' 'DILLS_CAMP' 'DISASTER_PEAK' 'DIVIDE' 'DONKEY_RESERVOIR'
     'DORSEY_BASIN' 'DRAW_CREEK' 'DRY_BREAD_POND' 'DRY_FORK' 'DRY_LAKE'
     'Defiance_Mines' 'Dollarhide_Summit' 'Dry_Creek' 'EAGLE_SUMMIT'
     'EAST_RIM_DIVIDE' 'EAST_WILLOW_CREEK' 'EBBETTS_PASS' 'ECHO_PEAK'
     'EF_BLACKS_FORK_GS' 'ELK_PEAK' 'EXIT_GLACIER' 'Elk_Cabin' 'FALLEN_LEAF'
     'FARMINGTON' 'FARMINGTON_LOWER' 'FARNSWORTH_LAKE' 'FAWN_CREEK'
     'FISH_LAKE_UTAH' 'FIVE_POINTS_LAKE' 'FLATTOP_MTN.' 'FORESTDALE_CREEK'
     'FORT_VALLEY' 'FRANKLIN_BASIN' 'FREMONT_PASS' 'Fry_Canyon' 'GALENA_SUMMIT'
     'GARDEN_CITY_SUMMIT' 'GARDNER_PEAK' 'GBRC_HQ' 'GBRC_MEADOWS'
     'GEORGE_CREEK' 'GIVEOUT' 'GOBBLERS_KNOB' 'GOLCONDA' 'GOLD_AXE_CAMP'
     'GOOSEBERRY_R.S.' 'GOOSEBERRY_R.S._UP' 'GRAND_TARGHEE' 'GRANITE_CRK'
     'GRANITE_PEAK' 'GREEN_MOUNTAIN' 'GRIZZLY_PEAK' 'GROUSE_CAMP' 'GUTZ_PEAK'
     'Gallegos_Peak' 'Granite_Creek' 'Gunsight_Pass' 'HAGANS_MEADOW'
     'HAPPY_JACK' 'HARDSCRABBLE' 'HARRIS_FLAT' 'HARTS_PASS' 'HAYDEN_FORK'
     'HEAVENLY_VALLEY' 'HEWINTA' 'HICKERSON_PARK' 'HIGH_RIDGE' 'HILTS_CREEK'
     'HOLE-IN-MOUNTAIN' 'HOLE-IN-ROCK' 'HOLLAND_MEADOWS' 'HOOSIER_PASS'
     'HORSE_MEADOW' 'HORSE_RIDGE' 'HYNDMAN' 'Hawley_Lake' 'Hobble_Creek'
     'Hopewell' 'Huntington_Horse' 'IMNAVIAT_CREEK' 'INDEPENDENCE_CAMP'
     'INDEPENDENCE_CREEK' 'INDEPENDENCE_LAKE' 'INDIAN_CANYON' 'INDIAN_ROCK'
     'JACKSON_PEAK' 'JACKS_PEAK' 'JACKWHACKER_GULCH' 'JACK_CREEK_UPPER'
     'JONES_CORRAL' 'Jakes_Creek' 'KALAMAZOO' 'KELLEY_R.S.' 'KELLY_STATION'
     'KENAI_MOOSE_PENS' 'KILFOIL_CREEK' 'KIMBERLY_MINE' 'KINGS_CABIN'
     'KLONDIKE_NARROWS' 'KOLOB' 'LAKEFORK_#1' 'LAKEFORK_#3' 'LAKEVIEW_RIDGE'
     'LAMANCE_CREEK' 'LAMOILLE_#3' 'LAPRELE_CREEK' 'LARSEN_CREEK'
     'LASAL_MOUNTAIN' 'LASAL_MOUNTAIN-LOWER' 'LAUREL_DRAW' 'LEAVITT_LAKE'
     'LEAVITT_MEADOWS' 'LEE_CANYON' 'LEWIS_LAKE_DIVIDE' 'LEWIS_PEAK'
     'LICK_CREEK' 'LIGHTNING_RIDGE' 'LILY_LAKE' 'LILY_POND' 'LITTLE_BEAR'
     'LITTLE_CHENA_RIDGE' 'LITTLE_GOOSE' 'LITTLE_GRASSY' 'LITTLE_SNAKE_RIVER'
     'LIZARD_HEAD_PASS' 'LOBDELL_LAKE' 'LONE_CONE' 'LONG_DRAW_RESV' 'LONG_FLAT'
     'LONG_VALLEY' 'LONG_VALLEY_JCT' 'LOOKOUT' 'LOOKOUT_PEAK' 'LOST_CREEK_RESV'
     'LOST_DOG' 'LOST_HORSE' 'LOUIS_MEADOW' 'LYNX_PASS' 'Lakefork_Basin'
     'Little_Valley' 'Lonesome_Beaver' 'MADISON_BUTTE' 'MAGIC_MOUNTAIN'
     'MAMMOTH-COTTONWOOD' 'MANY_GLACIER' 'MARLETTE_LAKE' 'MEDANO_PASS'
     'MERCHANT_VALLEY' 'MF_Nooksack' 'MICA_CREEK' 'MICHIGAN_CREEK'
     'MIDDLE_FORK_CAMP' 'MIDWAY_VALLEY' 'MILL-D_NORTH' 'MILLER_WOODS'
     'MINING_FORK' 'MONAHAN_FLAT' 'MONITOR_PASS' 'MONTE_CRISTO'
     'MONUMENT_CREEK' 'MOORE_CREEK_BRIDGE' 'MORMON_MTN_SUMMIT' 'MOSBY_MTN.'
     'MOSCOW_MOUNTAIN' 'MOSES_MTN' 'MOSQUITO_RIDGE' 'MOSS_SPRINGS'
     'MOUNT_LOCKHART' 'MT._HOWARD' 'MT._RYAN' 'MT_Baldy' 'MT_ROSE_SKI_AREA'
     'MUD_FLAT' 'MUNSON_RIDGE' 'MYRTLE_CREEK' 'Marten_Ridge' 'McNeil_River_SGS'
     'Med_Bow' 'Merritt_Mountain' 'Midas' 'Mormon_Mountain' 'Mt_Pennell'
     'NAVAJO_WHISKEY_CK' 'NEVADA_RIDGE' 'NUKA_GLACIER' 'OAK_CREEK' 'PALO'
     'PARADISE' 'PARK_CONE' 'PARK_CREEK_RIDGE' 'PARK_RESERVOIR'
     'PARLEYS_SUMMIT' 'PARRISH_CREEK' 'PAYSON_R.S.' 'PHANTOM_VALLEY'
     'PICKLE_KEG' 'PIERCE_R.S.' 'PINE_CREEK' 'POCKET_CREEK' 'POISON_FLAT'
     'POLE_CREEK_R.S.' 'PORPHYRY_CREEK' 'PORT_GRAHAM' 'PRUDHOE_BAY'
     'Panguitch_Lake_RS' 'Pole_Canyon' 'QUARTZ_MOUNTAIN' 'QUARTZ_PEAK'
     'Quemazon' 'RAGGED_MOUNTAIN' 'RAINBOW_CANYON' 'RAINY_PASS'
     'RED_PINE_RIDGE' 'RED_RIVER_PASS_#2' 'REYNOLDS_CREEK'
     'ROCKY_BASIN-SETTLEME' 'ROCKY_POINT' 'ROCK_CREEK' 'ROCK_SPRINGS'
     'ROUGH_AND_TUMBLE' 'RUBICON_#2' 'Redden_Mine_Lwr' 'Rees_Flat'
     'Rio_Santa_Barbara' 'SAGE_CREEK_BASIN' 'SALMON_MEADOWS' 'SALT_CREEK_FALLS'
     'SALT_RIVER_SUMMIT' 'SASSE_RIDGE' 'SAVAGE_PASS' 'SCHNEIDER_MEADOWS'
     'SCHOFIELD_PASS' 'SEELEY_CREEK' 'SENTINEL_BUTTE' 'SEVENTYSIX_CREEK'
     'SHANGHI_SUMMIT' 'SHARKSTOOTH' 'SHEEP_MTN.' 'SHUREE' 'SIERRA_BLANCA'
     'SILVER_CREEK' 'SILVIES' 'SLEEPING_WOMAN' 'SLUMGULLION' 'SMILEY_MOUNTAIN'
     'SMITH_and_MOREHOUSE' 'SNAKE_RIVER_STATION' 'SNOWBIRD' 'SNOW_MOUNTAIN'
     'SOLDIER_PARK' 'SOMSEN_RANCH' 'SONORA_PASS' 'SOURDOUGH_GULCH' 'SOUTH_MTN.'
     'SPIRIT_LK' 'SPRATT_CREEK' 'SPUR_PARK' 'SQUAW_SPRINGS' 'SQUAW_VALLEY_G.C.'
     'STEEL_CREEK_PARK' 'STRAWBERRY_DIVIDE' 'SUCKER_CREEK' 'SUMMIT_CREEK'
     'SUMMIT_LK' 'SUMMIT_MEADOW' 'SUMMIT_RANCH' 'SUSITNA_VALLEY_HIGH'
     'SWEDE_PEAK' 'Santa_Fe' 'Sawtooth' 'Senorita_Divide_#2' 'Sherwin'
     'Silver_Creek_Nv' 'Snowstorm_Mtn' 'Stag_Mountain' 'State_Line'
     'Sunflower_Flat' 'Suu_Ranch' 'TAHOE_CITY_CROSS' 'TAOS_POWDERHORN'
     'TAYLOR_BUTTE' 'TAYLOR_CANYON' 'TEMPLE_FORK' 'THAYNES_CANYON' 'TIMBERLINE'
     'TIMPANOGOS_DIVIDE' 'TIPTON' 'TOE_JAM' 'TOGWOTEE_PASS' 'TOKOSITNA_VALLEY'
     'TONY_GROVE_LAKE' 'TONY_GROVE_RS' 'TOUCHET' 'TOWNSEND_CREEK' 'TRIAL_LAKE'
     'TROUGH' 'TROUT_CREEK' 'TRUCKEE_#2' 'Takka_Wiiya' 'Tent_Mtn_Lower'
     'Thistle_Flat' 'Thumb_Divide' 'Tres_Ritos' 'UPPER_NOME_CREEK'
     'UPPER_RIO_GRANDE' 'UPPER_SAN_JUAN' 'UPPER_TAYLOR' 'UPPER_TSAINA_RIVER'
     'USU_DOC_DANIEL' 'Upper_Joes_Valley' 'VACARRO_SPRING' 'VAN_WYCK'
     'VERNON_CREEK' 'VIRGINIA_LAKES_RIDGE' 'Vacas_Locas' 'WARD_CREEK_#3'
     'WARD_MOUNTAIN' 'WATERHOLE' 'WEBSTER_FLAT' 'WESNER_SPRINGS'
     'WEST_YELLOWSTONE' 'WHEELER_PEAK' 'WHISKEY_CK' 'WHITE_HORSE_LAKE'
     'WHITE_MILL' 'WHITE_RIVER_#1' 'WIDTSOE_#3' 'WILDHORSE_DIVIDE' 'WILD_BASIN'
     'WILSON_CREEK' 'WINDY_PEAK' 'WOLF_CREEK_SUMMIT' 'White_River_Nv'
     'Wrigley_Creek' 'Yankee_Reservoir' 'ZIRKEL']


.. code:: python

    station = random.choice(stations)
    station_obj = ISMN_reader.get_station(station)
    print "Available Variables at Station %s"%station
    #get the variables that this station measures
    variables = station_obj.get_variables()
    print variables



.. parsed-literal::

    Available Variables at Station Hopewell
    ['air temperature' 'snow depth' 'snow water equivalent' 'soil moisture'
     'soil temperature']


.. code:: python

    #to make sure the selected variable is not measured
    #by different sensors at the same depths
    #we also select the first depth and the first sensor
    #even if there is only one
    depths_from,depths_to = station_obj.get_depths(variables[0])
    
    sensors = station_obj.get_sensors(variables[0],depths_from[0],depths_to[0])
    
    #read the data of the variable, depth, sensor combination
    time_series = station_obj.read_variable(variables[0],depth_from=depths_from[0],depth_to=depths_to[0],sensor=sensors[0])
    
    #print information about the selected time series
    print "Selected time series is:"
    print time_series



.. parsed-literal::

    Selected time series is:
    SNOTEL Hopewell -2.00 m - -2.00 m air temperature measured with n.s. 


.. code:: python

    #plot the data
    time_series.plot()
    plt.legend()
    plt.show()




.. image:: read_and_plot_ismn_data/output_8_0.png


.. code:: python

    #we also want to see soil moisture
    sm_depht_from,sm_depht_to = station_obj.get_depths('soil moisture')
    print sm_depht_from,sm_depht_to



.. parsed-literal::

    [ 0.2   0.51  0.05] [ 0.2   0.51  0.05]


.. code:: python

    #read sm data measured in first layer 0.2-0.2m
    sm = station_obj.read_variable('soil moisture',depth_from=0.2,depth_to=0.2)
    sm.plot()
    plt.show()




.. image:: read_and_plot_ismn_data/output_10_0.png


.. code:: python

    # the data attribute is a pandas.DataFrame
    time_series.data




.. raw:: html

    <div>
    <style>
        .dataframe thead tr:only-child th {
            text-align: right;
        }
    
        .dataframe thead th {
            text-align: left;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>air temperature</th>
          <th>air temperature_flag</th>
          <th>air temperature_orig_flag</th>
        </tr>
        <tr>
          <th>date_time</th>
          <th></th>
          <th></th>
          <th></th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>2016-01-01 00:00:00</th>
          <td>-10.5</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-01-01 01:00:00</th>
          <td>-11.4</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-01-01 02:00:00</th>
          <td>-12.1</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-01-01 03:00:00</th>
          <td>-12.0</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-01-01 04:00:00</th>
          <td>-12.4</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-01-01 05:00:00</th>
          <td>-12.3</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-01-01 06:00:00</th>
          <td>-13.5</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-01-01 07:00:00</th>
          <td>-14.1</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-01-01 08:00:00</th>
          <td>-13.8</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-01-01 09:00:00</th>
          <td>-14.1</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-01-01 10:00:00</th>
          <td>-14.3</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-01-01 11:00:00</th>
          <td>-14.2</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-01-01 12:00:00</th>
          <td>-14.4</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-01-01 13:00:00</th>
          <td>-13.8</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-01-01 14:00:00</th>
          <td>-12.9</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-01-01 15:00:00</th>
          <td>-12.8</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-01-01 16:00:00</th>
          <td>-9.8</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-01-01 17:00:00</th>
          <td>-7.2</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-01-01 18:00:00</th>
          <td>-5.4</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-01-01 19:00:00</th>
          <td>-4.1</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-01-01 20:00:00</th>
          <td>-4.0</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-01-01 21:00:00</th>
          <td>-2.6</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-01-01 22:00:00</th>
          <td>-2.9</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-01-01 23:00:00</th>
          <td>-5.8</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-01-02 00:00:00</th>
          <td>-8.5</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-01-02 01:00:00</th>
          <td>-9.4</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-01-02 02:00:00</th>
          <td>-8.8</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-01-02 03:00:00</th>
          <td>-9.9</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-01-02 04:00:00</th>
          <td>-9.2</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-01-02 05:00:00</th>
          <td>-10.3</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>2016-12-30 18:00:00</th>
          <td>2.3</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-12-30 19:00:00</th>
          <td>1.5</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-12-30 20:00:00</th>
          <td>1.1</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-12-30 21:00:00</th>
          <td>0.7</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-12-30 22:00:00</th>
          <td>1.4</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-12-30 23:00:00</th>
          <td>-0.7</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-12-31 00:00:00</th>
          <td>-2.4</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-12-31 01:00:00</th>
          <td>-3.4</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-12-31 02:00:00</th>
          <td>-3.5</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-12-31 03:00:00</th>
          <td>-4.0</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-12-31 04:00:00</th>
          <td>-4.1</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-12-31 05:00:00</th>
          <td>-5.9</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-12-31 06:00:00</th>
          <td>-5.3</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-12-31 07:00:00</th>
          <td>-4.1</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-12-31 08:00:00</th>
          <td>-5.7</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-12-31 09:00:00</th>
          <td>-4.3</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-12-31 10:00:00</th>
          <td>-3.2</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-12-31 11:00:00</th>
          <td>-3.2</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-12-31 12:00:00</th>
          <td>-3.4</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-12-31 13:00:00</th>
          <td>-3.1</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-12-31 14:00:00</th>
          <td>-2.8</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-12-31 15:00:00</th>
          <td>-2.6</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-12-31 16:00:00</th>
          <td>-1.5</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-12-31 17:00:00</th>
          <td>-0.3</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-12-31 18:00:00</th>
          <td>-0.8</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-12-31 19:00:00</th>
          <td>-0.1</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-12-31 20:00:00</th>
          <td>0.1</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-12-31 21:00:00</th>
          <td>-0.4</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-12-31 22:00:00</th>
          <td>-1.3</td>
          <td>G</td>
          <td>M</td>
        </tr>
        <tr>
          <th>2016-12-31 23:00:00</th>
          <td>-1.7</td>
          <td>G</td>
          <td>M</td>
        </tr>
      </tbody>
    </table>
    <p>8780 rows × 3 columns</p>
    </div>


Selection of ISMN stations by landcover or climate classification:

.. code:: python

    # Return all available landcover classifications (ESA CCI landcover 2000) for the variable soil moisture.
    # To use ESA CCI landcover data from the year 2005 or 2010 set landcover parameter to 'landcover_2005' and
    # 'landcover_2010', respectively.
    lc_2000 = ISMN_reader.get_landcover_types(variable='soil moisture', landcover='landcover_2000')

    # return all available landcover classifications (ESA CCI landcover 2005) for the variable soil moisture
    # (depths from 0 to 0.1m)
    lc_2005 = ISMN_reader.get_landcover_types(variable='soil moisture', landcover='landcover_2005'
                                              min_depth=0, max_depth=0.1)
    # return all available landcover classifications (ESA CCI landcover 2010) for the variable soil moisture
    # (depths from 0.1 to 0.5m)
    lc_2010 = ISMN_reader.get_landcover_types(variable='soil moisture', landcover='landcover_2010'
                                              min_depth=0.1, max_depth=0.5)
    # return all available landcover classifications (in situ) for the variable soil moisture
    lc_insitu = ISMN_reader.get_landcover_types(variable='soil moisture', landcover='landcover_insitu')

    # return all available climate classifications (Koeppen-Geiger 2007) for the variable soil moisture
    clim = ISMN_reader.get_climate_types(variable='soil moisture', climate='climate')
    # return all available climate classifications (in situ) for the variable soil moisture
    clim_insitu = ISMN_reader.get_climate_types(variable='soil moisture', climate='climate_insitu')


    # print all landcover classes covered by the ESA CCI landcover classification
    ISMN_reader.print_landcover_dict()
    # print all climate classes covered by the Koeppen-Geiger classification
    ISMN_reader.print_climate_dict()


    # Select ISMN stations where soil moisture at depths from 0 to 0.1m is available and the landcover
    # classification is equal to 130 (Grassland). In this example the ESA CCI landcover classification
    # for the year 2010 (landcover_2010) is used.
    ids1 = ISMN_reader.get_dataset_ids(variable='soil moisture', min_depth=0, max_depth=0.1, landcover_2010=130)
    # read time series from first element in the returned list
    ts_1 = ISMN_reader.read_ts(ids1[0])

    # Select ISMN stations where soil moisture at depths from 0 to 0.1m is available, the landcover
    # class (year 2005) is equal to 130 (Grassland) and the climate class is equal to Csa (Temperate
    # Dry Hot Summer)
    ids2 = ISMN_reader.get_dataset_ids(variable='soil moisture', min_depth=0, max_depth=1,
                                       landcover_2005=130, climate='Csa')
    # read time series from first element in the returned list
    ts_2 = ISMN_reader.read_ts(ids2[0])


