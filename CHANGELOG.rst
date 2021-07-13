=========
Changelog
=========

Unreleased
==========

-

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
