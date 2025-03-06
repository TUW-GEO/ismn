=========
Changelog
=========

Unreleased changes in master branch
===================================
-

Version 1.5.2
=============
- Fixed dependency list for `ismn[plot]` install (PyQt5 and cartopy>=0.22)
- Use `repurpose.process.parallel_process` now

Version 1.5.1
=============
- Fixed a bug where missing optional requirements (xarray, dask) lead to an error
- CI uses dependencies from setup.cfg now
- Switched to `treading` backend for multiprocessing

Version 1.5.0
=============
- Fixed a bug where metadata collection failed when a station static meta csv file was completely missing
- ISMN grid is now a CellGrid
- Added a method to sensors, stations and networks components for conversion to xarray Datasets

Version 1.4.1
=============
- Fixed a bug where some ISMN files could not be parsed correctly due to empty spaces in the sensor name (`Issue #81 <https://github.com/TUW-GEO/ismn/issues/81>`_)
- Parallel metadata collection now uses the repurpose package wrapper around joblib
- Logging was improved for metadata collection

Version 1.4.0
=============
- Added `ismn collect_metadata` and `ismn export_geojson` CLI programs.
- Added method to NetworkCollection to export metadata as (geo)json.
- Added more options when plotting the station overview map.
- Network citation list updated.
- Faster ISMN_Interface generation for large datasets by skipping some checks

Version 1.3.4
=============
- Fixed bug with parsing sensor files for RISMA network in "header & files" format

Version 1.3.3
=============
- Continue metadata collection for station on any type of error (PR #64)

Version 1.3.2
=============
- Fix bug where station names in metadata can be different between Header and CEOP format.
- Custom Sensor Metadata reader now also checks the measured variable.

Version 1.3.1
=============
- Added functionality to provide fill values with predefined custom metadata readers.
- Documentation and constants updated.
- ``ISMN_Interface.read_ts`` raised an error when np.int64 was passed.
- Fixed for `pygeogrids>=0.4.2` where dist=np.inf can be returned if

Version 1.3.0
=============
- Add module to assign custom metadata readers to ISMN_Interface
- Notebook added that describes using a custom metadata reader
- RTD build uses a separate, smaller environment.yml now (and mamba)
- ISMN_Interface now has a method to create an instance of itself for a selection of  (`ISMN_Interface.subset_from_ids`)

Version 1.2.0
=============
- Citation for `COSMOS-UK` network added
- Metadata for empty variables now also generated
- Added citation for `SMN-SDR` network (`Issue #44 <https://github.com/TUW-GEO/ismn/issues/44>`_)
- Update station plots with keyword argument for font size

Version 1.1.0
=============
- Fixed bug in the metadata dataframe when not using all available networks.
- Citation export functions added to `Collection` and `Network` components.
- Citations for all networks are now stored in this package until they are reliably provided together with the data
- It is now possible to pass multiple allowed variables when filtering sensors.

Version 1.0.1
=============
- Remove setuptools_scm from dependencies
- Remove option to pip install ismn[plot] for now.

Version 1.0
===========
- Rewrite package, objects for Networks, Stations, Sensors etc.
- Update ISMN_Interface to use new components, similar behaviour to old interface.
- Add MetaVar and MetaData modules for ismn metadata handling.
- Data readers are now collected in filehandler class.
- Metadata for soil layers is now filter by the sensor depth.
- Support reading from zip archive and extracted folders.
- New python_metadata structure (no absolute paths, no pickle format).
- Drop support for old ceop format (new ceop_sep format is supported!).
- Move lookup tables to const.py
- Update docs and tests.
- Use Github Actions instead of Travis CI.

Version 0.4
===========
- Update package pyscaffold package structure
- Drop python 2 support
- Add read function to ISMN_Interface
- Add travis pypi deployment
- Fix get_min_max_timestamp_header_values

Version 0.3.2
=============
- Add function to initialise different network

Version 0.3.1
=============
- Set allow_pickle to True when loading metadata

Version 0.3
===========
- Update readme
- Added information about landcover and climate to metadata.

Version 0.2
===========
- Add additional authors.

Version 0.1
===========
- Moved code from pytesmo into this package.
