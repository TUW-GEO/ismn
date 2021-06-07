====
ismn
====


.. image:: https://github.com/TUW-GEO/ismn/workflows/Automated%20Tests/badge.svg?branch=master&event=push
   :target: https://github.com/TUW-GEO/ismn/actions

.. image:: https://coveralls.io/repos/TUW-GEO/ismn/badge.png?branch=master
  :target: https://coveralls.io/r/TUW-GEO/ismn?branch=master

.. image:: https://badge.fury.io/py/ismn.svg
    :target: http://badge.fury.io/py/ismn

.. image:: https://readthedocs.org/projects/ismn/badge/?version=latest
   :target: http://ismn.readthedocs.org/

Readers for the data from the International Soil Moisture Database (ISMN).

Documentation
-------------
The full documentation is available at https://ismn.readthedocs.io and includes
a tutorial on reading ISMN data in python after downloading it from
https://ismn.earth

Citation
========

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.855308.svg
   :target: https://doi.org/10.5281/zenodo.855308

If you use the software in a publication then please cite it using the Zenodo DOI.
Be aware that this badge links to the latest package version.

Please select your specific version at https://doi.org/10.5281/zenodo.855308 to get the DOI of that version.
You should normally always use the DOI for the specific version of your record in citations.
This is to ensure that other researchers can access the exact research artefact you used for reproducibility.

You can find additional information regarding DOI versioning at http://help.zenodo.org/#versioning

Installation
============

This package should be installable through pip:

.. code::

    pip install ismn

Optional dependencies
---------------------

The ``cartopy`` and ``matplotlib`` packages are only needed when creating data visualisations.
They can be installed separately with:

.. code::

    conda install -c conda-forge matplotlib
    conda install -c conda-forge cartopy

Example installation script
---------------------------

The following script will install miniconda and setup the environment on a UNIX
like system. Miniconda will be installed into ``$HOME/miniconda``.

.. code::

   wget https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh
   bash miniconda.sh -b -p $HOME/miniconda
   export PATH="$HOME/miniconda/bin:$PATH"
   git clone git@github.com:TUW-GEO/ismn.git ismn
   cd ismn
   conda env create -f environment.yml
   source activate ismn

This script adds ``$HOME/miniconda/bin`` temporarily to the ``PATH`` to do this
permanently add ``export PATH="$HOME/miniconda/bin:$PATH"`` to your ``.bashrc``
or ``.zshrc``

The second to last line in the example activates the ``ismn`` environment.

After that you should be able to run:

.. code::

    pytest

to run the test suite.

Description
===========

ISMN data can be downloaded for free after creating an account on the `ISMN Website
<http://ismn.geo.tuwien.ac.at/>`_

ISMN data can be downloaded in two different formats:

* Variables stored in separate files (CEOP formatted)

	this format is supported 100% and should work with all examples

* Variables stored in separate files (Header+values)

	this format is supported 100% and should work with all examples

If you downloaded ISMN data in one of the supported formats in the past it can
be that station names are not recognized correctly because they contained the
'_' character which is supposed to be the separator. If you experience problems
because of this please download new data from the ISMN since this issue should
be fixed.

Landcover Classification
------------------------
The ISMN data comes with information about landcover classification from the
ESA CCI land cover project (years 2000, 2005 and 2010) as well as from in-situ
measurements. To use ESA CCI land cover variables for filtering the data in the get_dataset_ids
function, set the keyword parameters (landcover_2000, landcover_2005 or landcover_2010)
to the corresponding integer values (e.g. 10) in the list below. To get a list of
possible values for filtering by in-situ values (keyword parameter: "landcover_insitu"),
call the get_landcover_types method of your ISMN_Interface object and set landcover='landcover_insitu'.

* 10: Cropland, rainfed
* 11: Cropland, rainfed / Herbaceous cover
* 12: Cropland, rainfed / Tree or shrub cover,
* 20: Cropland, irrigated or post-flooding,
* 30: Mosaic cropland (>50%) / natural vegetation (tree, shrub, herbaceous,
* 40: Mosaic natural vegetation (tree, shrub, herbaceous cover) (>50%) / cropland (<50%),
* 50: Tree cover, broadleaved, evergreen, Closed to open (>15%),
* 60: Tree cover, broadleaved, deciduous, Closed to open (>15%),
* 61: Tree cover, broadleaved, deciduous, Closed (>40%),
* 62: Tree cover, broadleaved, deciduous, Open (15-40%),
* 70: Tree cover, needleleaved, evergreen, closed to open (>15%),
* 71: Tree cover, needleleaved, evergreen, closed (>40%),
* 72: Tree cover, needleleaved, evergreen, open (15-40%),
* 80: Tree cover, needleleaved, deciduous, closed to open (>15%),
* 81: Tree cover, needleleaved, deciduous, closed (>40%),
* 82: Tree cover, needleleaved, deciduous, open (15-40%),
* 90: Tree cover, mixed leaf type (broadleaved and needleleaved),
* 100: Mosaic tree and shrub (>50%) / herbaceous cover (<50%),
* 110: Mosaic herbaceous cover (>50%) / tree and shrub (<50%),
* 120: Shrubland,
* 121: Shrubland / Evergreen Shrubland,
* 122: Shrubland / Deciduous Shrubland,
* 130: Grassland,
* 140: Lichens and mosses,
* 150: Sparse vegetation (tree, shrub, herbaceous cover) (<15%),
* 152: Sparse vegetation (tree, shrub, herbaceous cover) (<15%) / Sparse shrub (<15%),
* 153: Sparse vegetation (tree, shrub, herbaceous cover) (<15%) / Sparse herbaceous cover (<15%),
* 160: Tree cover, flooded, fresh or brakish water,
* 170: Tree cover, flooded, saline water,
* 180: Shrub or herbaceous cover, flooded, fresh/saline/brakish water,
* 190: Urban areas,
* 200: Bare areas,
* 201: Consolidated bare areas,
* 202: Unconsolidated bare areas,
* 210: Water,
* 220: Permanent snow and ice,

Climate Classification
----------------------
The ISMN data comes with information about climate classification from the Koeppen-Geiger
Climate Classification (2007) as well as in-situ measurements. To use
Koeppen-Geiger variable for filtering the data in the get_dataset_ids function, set the
keyword parameter "climate" to the corresponding keys (e.g. 'Af') in the list below. To get a list of
possible values for filtering by in-situ values (keyword parameter: "climate_insitu"), call the
get_climate_types method of your ISMN_Interface object and set climate='climate_insitu'.

* Af: Tropical Rainforest
* Am: Tropical Monsoon
* As: Tropical Savanna Dry
* Aw: Tropical Savanna Wet
* BWk: Arid Desert Cold
* BWh: Arid Desert Hot
* BWn: Arid Desert With Frequent Fog
* BSk: Arid Steppe Cold
* BSh: Arid Steppe Hot
* BSn: Arid Steppe With Frequent Fog
* Csa: Temperate Dry Hot Summer
* Csb: Temperate Dry Warm Summer
* Csc: Temperate Dry Cold Summer
* Cwa: Temperate Dry Winter, Hot Summer
* Cwb: Temperate Dry Winter, Warm Summer
* Cwc: Temperate Dry Winter, Cold Summer
* Cfa: Temperate Without Dry Season, Hot Summer
* Cfb: Temperate Without Dry Season, Warm Summer
* Cfc: Temperate Without Dry Season, Cold Summer
* Dsa: Cold Dry Summer, Hot Summer
* Dsb: Cold Dry Summer, Warm Summer
* Dsc: Cold Dry Summer, Cold Summer
* Dsd: Cold Dry Summer, Very Cold Winter
* Dwa: Cold Dry Winter, Hot Summer
* Dwb: Cold Dry Winter, Warm Summer
* Dwc: Cold Dry Winter, Cold Summer
* Dwd: Cold Dry Winter, Very Cold Winter
* Dfa: Cold Dry Without Dry Season, Hot Summer
* Dfb: Cold Dry Without Dry Season, Warm Summer
* Dfc: Cold Dry Without Dry Season, Cold Summer
* Dfd: Cold Dry Without Dry Season, Very Cold Winter
* ET: Polar Tundra
* EF: Polar Eternal Winter
* W: Water


Contribute
==========

We are happy if you want to contribute. Please raise an issue explaining what
is missing or if you find a bug. We will also gladly accept pull requests
against our master branch for new features or bug fixes.

Development setup
-----------------

For Development we also recommend a ``conda`` environment. You can create one
including test dependencies and debugger by running
``conda env create -f environment.yml``. This will create a new
``ismn`` environment which you can activate by using
``conda activate ismn``.

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
