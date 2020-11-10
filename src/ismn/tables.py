# -*- coding: utf-8 -*-

import numpy as np
from collections import OrderedDict

# ==============================================================================
# Variable short names

# todo: long names dont match with old data, e.g soil_moisture vs soil moisture
VARIABLE_LUT = OrderedDict([
    ('sm', 'soil_moisture'),
    ('ts', 'soil_temperature'),
    ('su', 'soil_suction'),
    ('p', 'precipitation'),
    ('ta', 'air_temperature'),
    ('fc', 'field_capacity'),
    ('wp', 'permanent_wilting_point'),
    ('paw', 'plant_available_water'),
    ('ppaw', 'potential_plant_available_water'),
    ('sat', 'saturation'),
    ('si_h', 'silt_fraction'),
    ('sd', 'snow_depth'),
    ('sa_h', 'sand_fraction'),
    ('cl_h', 'clay_fraction'),
    ('oc_h', 'organic_carbon'),
    ('sweq', 'snow_water_equivalent'),
    ('tsf', 'surface_temperature'),
    ('tsfq', 'surface_temperature_quality_flag_original'),
    ])

# ==============================================================================
# static meta data template (csv)

CSV_COLS = ['quantity_name', 'unit', 'depth_from[m]', 'depth_to[m]',
            'value', 'description', 'quantity_source_name']

CSV_META_TEMPLATE_SURF_VAR = OrderedDict([
    ('lc_2000', np.nan),
    ('lc_2005', np.nan),
    ('lc_2010', np.nan),
    ('lc_insitu', ''),
    ('climate_KG', ''),
    ('climate_insitu', ''),
    ])
CSV_META_TEMPLATE_GROUND_VAR = OrderedDict([
    ('saturation', np.nan),
    ('clay_fraction', np.nan),
    ('sand_fraction', np.nan),
    ('silt_fraction', np.nan),
    ('organic_carbon', np.nan),
    ])

CSV_META_TEMPLATE = OrderedDict(**CSV_META_TEMPLATE_SURF_VAR,
                                **CSV_META_TEMPLATE_GROUND_VAR)
# ==============================================================================
# ESA CCI landcover classification:
# This landcover information is available for the years 2000, 2005 and 2010.
# To filter the data by landcover classes in the get_dataset_ids function,
# set the parameter (landcover_2000, landcover_2005 or landcover_2010) to an
# integer value in the list below (e.g. landcover_2000=10).

LANDCOVER = OrderedDict([
    ('Cropland, rainfed', 10),
    ('Cropland, rainfed / Herbaceous cover', 11),
    ('Cropland, rainfed / Tree or shrub cover', 12),
    ('Cropland, irrigated or post-flooding', 20),
    ('Mosaic cropland (>50%) / natural vegetation (tree, shrub, herbaceous', 30),
    ('Mosaic natural vegetation (tree, shrub, herbaceous cover) (>50%) / cropland (<50%)', 40),
    ('Tree cover, broadleaved, evergreen, Closed to open (>15%)', 50),
    ('Tree cover, broadleaved, deciduous, Closed to open (>15%)', 60),
    ('Tree cover, broadleaved, deciduous, Closed (>40%)', 61),
    ('Tree cover, broadleaved, deciduous, Open (15-40%)', 62),
    ('Tree cover, needleleaved, evergreen, closed to open (>15%)', 70),
    ('Tree cover, needleleaved, evergreen, closed (>40%)', 71),
    ('Tree cover, needleleaved, evergreen, open (15-40%)', 72),
    ('Tree cover, needleleaved, deciduous, closed to open (>15%)', 80),
    ('Tree cover, needleleaved, deciduous, closed (>40%)', 81),
    ('Tree cover, needleleaved, deciduous, open (15-40%)', 82),
    ('Tree cover, mixed leaf type (broadleaved and needleleaved)', 90),
    ('Mosaic tree and shrub (>50%) / herbaceous cover (<50%)', 100),
    ('Mosaic herbaceous cover (>50%) / tree and shrub (<50%)', 110),
    ('Shrubland', 120),
    ('Shrubland / Evergreen Shrubland', 121),
    ('Shrubland / Deciduous Shrubland', 122),
    ('Grassland', 130),
    ('Lichens and mosses', 140),
    ('Sparse vegetation (tree, shrub, herbaceous cover) (<15%)', 150),
    ('Sparse vegetation (tree, shrub, herbaceous cover) (<15%) / Sparse shrub (<15%)', 152),
    ('Sparse vegetation (tree, shrub, herbaceous cover) (<15%) / Sparse herbaceous cover (<15%)', 153),
    ('Tree cover, flooded, fresh or brakish water', 160),
    ('Tree cover, flooded, saline water', 170),
    ('Shrub or herbaceous cover, flooded, fresh/saline/brakish water', 180),
    ('Urban areas', 190),
    ('Bare areas', 200),
    ('Consolidated bare areas', 201),
    ('Unconsolidated bare areas', 202),
    ('Water', 210),
    ('Permanent snow and ice', 220),
    ])

# ==============================================================================
# Koeppen-Geiger 2007 climate classification:
# To filter the data by Koeppen-Geiger climate classes (2007) in the get_dataset_ids function,
# set the parameter 'climate' to a valid value from the list below (e.g.
# climate='Af').

KOEPPENGEIGER = OrderedDict([
    ('Af', 'Tropical Rainforest'),
    ('Am', 'Tropical Monsoon'),
    ('As', 'Tropical Savanna Dry'),
    ('Aw', 'Tropical Savanna Wet'),
    ('BWk', 'Arid Desert Cold'),
    ('BWh', 'Arid Desert Hot'),
    ('BWn', 'Arid Desert With Frequent Fog'),
    ('BSk', 'Arid Steppe Cold'),
    ('BSh', 'Arid Steppe Hot'),
    ('BSn', 'Arid Steppe With Frequent Fog'),
    ('Csa', 'Temperate Dry Hot Summer'),
    ('Csb', 'Temperate Dry Warm Summer'),
    ('Csc', 'Temperate Dry Cold Summer'),
    ('Cwa', 'Temperate Dry Winter, Hot Summer'),
    ('Cwb', 'Temperate Dry Winter, Warm Summer'),
    ('Cwc', 'Temperate Dry Winter, Cold Summer'),
    ('Cfa', 'Temperate Without Dry Season, Hot Summer'),
    ('Cfb', 'Temperate Without Dry Season, Warm Summer'),
    ('Cfc', 'Temperate Without Dry Season, Cold Summer'),
    ('Dsa', 'Cold Dry Summer, Hot Summer'),
    ('Dsb', 'Cold Dry Summer, Warm Summer'),
    ('Dsc', 'Cold Dry Summer, Cold Summer'),
    ('Dsd', 'Cold Dry Summer, Very Cold Winter'),
    ('Dwa', 'Cold Dry Winter, Hot Summer'),
    ('Dwb', 'Cold Dry Winter, Warm Summer'),
    ('Dwc', 'Cold Dry Winter, Cold Summer'),
    ('Dwd', 'Cold Dry Winter, Very Cold Winter'),
    ('Dfa', 'Cold Dry Without Dry Season, Hot Summer'),
    ('Dfb', 'Cold Dry Without Dry Season, Warm Summer'),
    ('Dfc', 'Cold Dry Without Dry Season, Cold Summer'),
    ('Dfd', 'Cold Dry Without Dry Season, Very Cold Winter'),
    ('ET', 'Polar Tundra'),
    ('EF', 'Polar Eternal Winter'),
    ('W', 'Water'),
    ])

# ==============================================================================
