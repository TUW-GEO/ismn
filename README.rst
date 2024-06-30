====
ismn
====

|ci| |cov| |pip| |doc|

.. |ci| image:: https://github.com/TUW-GEO/ismn/actions/workflows/build.yml/badge.svg?branch=master
   :target: https://github.com/TUW-GEO/ismn/actions

.. |cov| image:: https://coveralls.io/repos/TUW-GEO/ismn/badge.png?branch=master
  :target: https://coveralls.io/r/TUW-GEO/ismn?branch=master

.. |pip| image:: https://badge.fury.io/py/ismn.svg
    :target: http://badge.fury.io/py/ismn

.. |doc| image:: https://readthedocs.org/projects/ismn/badge/?version=latest
   :target: http://ismn.readthedocs.org/

Readers for data from the International Soil Moisture Network (ISMN) https://ismn.earth.

This package is installable through pip (verified for ``python>=3.9``)

.. code::

    pip install ismn

This should also install all dependencies on Window, Linux or MacOS.

Quickstart
----------
Initialise an `ISMN_Interface` by passing the path to your downloaded data.
The interface shows you available ISMN networks, stations and sensors.
You can load sensor time series as pandas DataFrames as well as metadata
on the stations surroundings, soil conditions and probes
(depths, sensor type, etc.).

.. code-block:: python

    >> from ismn.interface import ISMN_Interface
    """ .zip archives are downloaded from https://ismn.earth """
    >> ds = ISMN_Interface('Data_separate_files_header_20090101_20201231_9289_Cwpc_20221201.zip')
    """ Read time series from your previously downloaded ISMN archive as pandas DataFrames """
    >> ds["REMEDHUS"]["Canizal"][0].data

        Out[0]:
                                 soil_moisture soil_moisture_flag soil_moisture_orig_flag
        date_time
        2009-01-01 00:00:00          0.372                  G                       M
        2009-01-01 01:00:00          0.372                  G                       M
        ...                          ...                   ...                     ...
        2020-12-31 22:00:00          0.285                  G                       M
        2020-12-31 23:00:00          0.285                  G                       M

    """ Each ISMN sensor comes with additional information on soil/landcover/climate etc. """
    >> ds["REMEDHUS"]["Canizal"][0].metadata.to_pd()

        Out[0]:
        variable        key
        climate_KG      val                           BSk
        instrument      val           Stevens-Hydra-Probe
                        depth_from                    0.0
                        depth_to                     0.05
        ...             ...                           ...
        latitude        val                      41.19603
        lc_2010         val                            20
        longitude       val                      -5.35997
        network         val                      REMEDHUS
        station         val                       Canizal

Many more features to e.g. visualise, select or transform data are available.
See the `full documentation <https://ismn.readthedocs.io/en/latest/>`_.

Documentation
-------------
The full documentation is available at https://ismn.readthedocs.io/en/latest and includes
a tutorial on reading ISMN data in python after downloading it from
https://ismn.earth

The following **tutorials** are also available as ipython notebooks in ``docs/examples``:

 #. `ISMN reader basic functionality <https://ismn.readthedocs.io/en/latest/examples/interface.html>`_
 #. `Adding custom metadata readers <https://ismn.readthedocs.io/en/latest/examples/custom_meta.html>`_

Data used in the tutorials is *not* provided in this package. Please create an account at `ismn.earth <https://ismn.earth/en/>`_
to download the required files.

For a detailed description of the ISMN, technical data aspects (properties, coverage, etc.) and correct usage (applications), see

    W. Dorigo et al. **The International Soil Moisture Network: serving Earth system science for over a decade**,
    Hydrol. Earth Syst. Sci., 25, 5749–5804, https://doi.org/10.5194/hess-25-5749-2021, 2021.

Optional dependencies
---------------------

The `cartopy <https://github.com/SciTools/cartopy>`_ and `matplotlib <https://github.com/matplotlib/matplotlib>`_ packages
are only needed when creating data visualisations. They can be installed separately via

.. code::

    conda install -c conda-forge matplotlib cartopy

or ``pip install ismn[plot]`` for most operating systems if you already have `geos <https://libgeos.org/>`_ installed.

If you want to convert ISMN data into xarray objects, please install `xarray <https://github.com/pydata/xarray>`_ and
`dask <https://github.com/dask/dask>`_ via

.. code::

    conda install -c conda-forge xarray dask

or ``pip install ismn[xr]`` for most operating systems.

Citation
--------

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.855308.svg
   :target: https://doi.org/10.5281/zenodo.855308

If you use the software in a publication then please cite it using the Zenodo DOI.
Be aware that this badge links to the latest package version.

Please select your specific version at https://doi.org/10.5281/zenodo.855308 to get the DOI of that version.
You should normally always use the DOI for the specific version of your record in citations.
This is to ensure that other researchers can access the exact research artefact you used for reproducibility.

You can find additional information regarding DOI versioning at http://help.zenodo.org/#versioning

Description
===========

ISMN data can be downloaded for free after creating an account on the `ISMN Website
<http://ismn.geo.tuwien.ac.at/>`_

ISMN data can be downloaded in two different formats:

* Variables stored in separate files (CEOP formatted)
* Variables stored in separate files (Header+values) (default format)

Both formats are supported by this package.

If you downloaded ISMN data in one of the supported formats in the past it can
be that station names are not recognized correctly because they contained the
'_' character which is supposed to be the separator. If you experience problems
because of this please download new data from the ISMN since this issue should
be fixed.

Variables and Units
-------------------
The following variables are available in the ISMN. Note that not every station
measures all of the variables. You can use this package to read only data for
locations where one or multiple of the variables were measured.

.. list-table:: Temporally dynamic variables and their units in ISMN
   :widths: 25 15
   :header-rows: 1

   * - Variable
     - Units
   * - Soil Moisture
     - m\ :sup:`3`\ /m\ :sup:`3`\
   * - Soil Suction
     - kPa
   * - Soil Temperature
     - °C
   * - Air Temperature
     - °C
   * - Surface Temperature
     - °C
   * - Precipitation
     - mm
   * - Snow Depth
     - mm
   * - Snow Water Equivalent
     - mm

----

.. list-table:: Temporally static variables and their units in ISMN
   :widths: 35 35
   :header-rows: 1

   * - Variable
     - Units
   * - Climate classification
     - None
   * - Land cover classification
     - None
   * - Soil classification
     - None
   * - Bulk density
     - g/cm³
   * - Sand fraction
     - % weight
   * - Silt fraction
     - % weight
   * - Clay fraction
     - % weight
   * - Organic carbon
     - % weight
   * - Saturation
     - % vol
   * - Field capacity
     - % vol
   * - Potential plant available water
     - % vol
   * - Permanent wilting point
     - % vol

Landcover Classification
------------------------
The ISMN data comes with information about landcover classification from the
ESA CCI land cover project (years 2000, 2005 and 2010) as well as from in-situ
measurements. To use ESA CCI land cover variables for filtering the data in the get_dataset_ids
function, set the keyword parameters (landcover_2000, landcover_2005 or landcover_2010)
to the corresponding integer values (e.g. 10) in the list below. To get a list of
possible values for filtering by in-situ values (keyword parameter: "landcover_insitu"),
call the get_landcover_types method of your ISMN_Interface object and set landcover='landcover_insitu'.

.. list-table:: ISMN Landcover classes and meanings
   :widths: 5 50
   :header-rows: 1

   * - Value
     - Meaning
   * - 10
     - Cropland, rainfed
   * - 11
     - Cropland, rainfed / Herbaceous cover
   * - 12
     - Cropland, rainfed / Tree or shrub cover
   * - 20
     - Cropland, irrigated or post-flooding
   * - 30
     - Mosaic cropland (>50%) / natural vegetation (tree, shrub, herbaceous)
   * - 40
     - Mosaic natural vegetation (>50%) / cropland (<50%)
   * - 50
     - Tree cover, broadleaved, evergreen, Closed to open (>15%)
   * - 60
     - Tree cover, broadleaved, deciduous, Closed to open (>15%)
   * - 61
     - Tree cover, broadleaved, deciduous, Closed (>40%)
   * - 62
     - Tree cover, broadleaved, deciduous, Open (15-40%)
   * - 70
     - Tree cover, needleleaved, evergreen, Closed to open (>15%)
   * - 71
     - Tree cover, needleleaved, evergreen, Closed (>40%)
   * - 72
     - Tree cover, needleleaved, evergreen, Open (15-40%)
   * - 80
     - Tree cover, needleleaved, deciduous, Closed to open (>15%)
   * - 81
     - Tree cover, needleleaved, deciduous, Closed (>40%)
   * - 82
     - Tree cover, needleleaved, deciduous, Open (15-40%)
   * - 90
     - Tree cover, mixed leaf type (broadleaved and needleleaved)
   * - 100
     - Mosaic tree and shrub (>50%) / herbaceous cover (<50%)
   * - 110
     - Mosaic herbaceous cover (>50%) / tree and shrub (<50%)
   * - 120
     - Shrubland
   * - 121
     - Shrubland / Evergreen Shrubland
   * - 122
     - Shrubland / Deciduous Shrubland
   * - 130
     - Grassland
   * - 140
     - Lichens and mosses
   * - 150
     - Sparse vegetation (tree, shrub, herbaceous cover) (<15%)
   * - 152
     - Sparse vegetation (<15%) / Sparse shrub (<15%)
   * - 153
     - Sparse vegetation (<15%) / Sparse herbaceous cover (<15%)
   * - 160
     - Tree cover, flooded, fresh or brackish water
   * - 170
     - Tree cover, flooded, saline water
   * - 180
     - Shrub or herbaceous cover, flooded, fresh/saline/brackish water
   * - 190
     - Urban areas
   * - 200
     - Bare areas
   * - 201
     - Consolidated bare areas
   * - 202
     - Unconsolidated bare areas
   * - 210
     - Water
   * - 220
     - Permanent snow and ice

Climate Classification
----------------------
The ISMN data comes with information about climate classification from the Koeppen-Geiger
Climate Classification (2007) as well as in-situ measurements. To use
Koeppen-Geiger variable for filtering the data in the get_dataset_ids function, set the
keyword parameter "climate" to the corresponding keys (e.g. 'Af') in the list below. To get a list of
possible values for filtering by in-situ values (keyword parameter: "climate_insitu"), call the
get_climate_types method of your ISMN_Interface object and set climate='climate_insitu'.

.. list-table:: Climate Classes and Meanings
   :widths: 5 50
   :header-rows: 1

   * - Class
     - Meaning
   * - Af
     - Tropical Rainforest
   * - Am
     - Tropical Monsoon
   * - As
     - Tropical Savanna Dry
   * - Aw
     - Tropical Savanna Wet
   * - BWk
     - Arid Desert Cold
   * - BWh
     - Arid Desert Hot
   * - BWn
     - Arid Desert With Frequent Fog
   * - BSk
     - Arid Steppe Cold
   * - BSh
     - Arid Steppe Hot
   * - BSn
     - Arid Steppe With Frequent Fog
   * - Csa
     - Temperate Dry Hot Summer
   * - Csb
     - Temperate Dry Warm Summer
   * - Csc
     - Temperate Dry Cold Summer
   * - Cwa
     - Temperate Dry Winter, Hot Summer
   * - Cwb
     - Temperate Dry Winter, Warm Summer
   * - Cwc
     - Temperate Dry Winter, Cold Summer
   * - Cfa
     - Temperate Without Dry Season, Hot Summer
   * - Cfb
     - Temperate Without Dry Season, Warm Summer
   * - Cfc
     - Temperate Without Dry Season, Cold Summer
   * - Dsa
     - Cold Dry Summer, Hot Summer
   * - Dsb
     - Cold Dry Summer, Warm Summer
   * - Dsc
     - Cold Dry Summer, Cold Summer
   * - Dsd
     - Cold Dry Summer, Very Cold Winter
   * - Dwa
     - Cold Dry Winter, Hot Summer
   * - Dwb
     - Cold Dry Winter, Warm Summer
   * - Dwc
     - Cold Dry Winter, Cold Summer
   * - Dwd
     - Cold Dry Winter, Very Cold Winter
   * - Dfa
     - Cold Dry Without Dry Season, Hot Summer
   * - Dfb
     - Cold Dry Without Dry Season, Warm Summer
   * - Dfc
     - Cold Dry Without Dry Season, Cold Summer
   * - Dfd
     - Cold Dry Without Dry Season, Very Cold Winter
   * - ET
     - Polar Tundra
   * - EF
     - Polar Eternal Winter
   * - W
     - Water


Contribute
==========

We are happy if you want to contribute. Please raise an issue explaining what
is missing or if you find a bug. We will also gladly accept pull requests
against our master branch for new features or bug fixes.


Guidelines
----------

If you want to contribute please follow these steps:

- Fork the ismn repository to your account
- Clone the repository
- make a new feature branch from the ismn master branch
- Add your feature
- Please include tests for your contributions in one of the test directories.
  We use pytest so a simple function called test_my_feature is enough
- submit a pull request to our master branch

Code Formatting
---------------
To apply pep8 conform styling to any changed files [we use `yapf`](https://github.com/google/yapf). The correct
settings are already set in `setup.cfg`. Therefore the following command
should be enough:

    yapf file.py --in-place

Release new version
-------------------

To release a new version of this package, make sure all tests are passing on the
master branch and the CHANGELOG.rst is up-to-date, with changes for the new version
at the top.

Then draft a new release at https://github.com/TUW-GEO/ismn/releases.
Create a version tag following the ``v{MAJOR}.{MINOR}.{PATCH}`` pattern.
This will trigger a new build on GitHub and should push the packages to pypi after
all tests have passed.

If this does not work (tests pass but upload fails) you can download the
``whl`` and ``dist`` packages for each workflow run from
https://github.com/TUW-GEO/ismn/actions (Artifacts) and push them manually to
https://pypi.org/project/ismn/ (you need to be a package maintainer on pypi for that).

In any case, ``pip install ismn`` should download the newest version afterwards.
