# -*- coding: utf-8 -*-

variable_lookup = {'sm': 'soil_moisture',
                   'ts': 'soil_temperature',
                   'su': 'soil_suction',
                   'p': 'precipitation',
                   'ta': 'air_temperature',
                   'fc': 'field_capacity',
                   'wp': 'permanent_wilting_point',
                   'paw': 'plant_available_water',
                   'ppaw': 'potential_plant_available_water',
                   'sat': 'saturation',
                   'si_h': 'silt_fraction',
                   'sd': 'snow_depth',
                   'sa_h': 'sand_fraction',
                   'cl_h': 'clay_fraction',
                   'oc_h': 'organic_carbon',
                   'sweq': 'snow_water_equivalent',
                   'tsf': 'surface_temperature',
                   'tsfq': 'surface_temperature_quality_flag_original'}

# static meta (csv)
import numpy as np
from collections import OrderedDict

csv_cols = ['quantity_name', 'unit', 'depth_from[m]', 'depth_to[m]',
            'value', 'description', 'quantity_source_name']

csv_metadata_template = OrderedDict([('lc_2000', np.nan),
                                     ('lc_2005', np.nan),
                                     ('lc_2010', np.nan),
                                     ('lc_insitu', ''),
                                     ('climate_KG', ''),
                                     ('climate_insitu', ''),
                                     ('saturation', np.nan),
                                     ('clay_fraction', np.nan),
                                     ('sand_fraction', np.nan),
                                     ('silt_fraction', np.nan),
                                     ('organic_carbon', np.nan),
                                   ])

