====
ismn
====

.. image:: https://travis-ci.org/TUW-GEO/ismn.svg?branch=master
    :target: https://travis-ci.org/TUW-GEO/ismn

.. image:: https://coveralls.io/repos/TUW-GEO/ismn/badge.png?branch=master
  :target: https://coveralls.io/r/TUW-GEO/ismn?branch=master

.. image:: https://badge.fury.io/py/ismn.svg
    :target: http://badge.fury.io/py/ismn

.. image:: https://zenodo.org/badge/101878880.svg
   :target: https://zenodo.org/badge/latestdoi/101878880

.. image:: https://readthedocs.org/projects/ismn/badge/?version=latest
   :target: http://ismn.readthedocs.org/

Readers for the data from the International Soil Moisture Database (ISMN).

Citation
========

If you use the software please cite it using the following DOI:

.. image:: https://zenodo.org/badge/101878880.svg
   :target: https://zenodo.org/badge/latestdoi/101878880

Installation
============

This package should be installable through pip:

.. code::

    pip install ismn

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
   pip install -r test-requirements.txt

This script adds ``$HOME/miniconda/bin`` temporarily to the ``PATH`` to do this
permanently add ``export PATH="$HOME/miniconda/bin:$PATH"`` to your ``.bashrc``
or ``.zshrc``

The second to last line in the example activates the ``ismn`` environment.

After that you should be able to run:

.. code::

    python setup.py test

to run the test suite.

Description
===========

ISMN data can be downloaded for free after registration from the `ISMN Website
<http://ismn.geo.tuwien.ac.at/>`_

In case of the ISMN, two different formats are provided:

* Variables stored in separate files (CEOP formatted)

	this format is supported 100% and should work with all examples

* Variables stored in separate files (Header+values)

	this format is supported 100% and should work with all examples

If you downloaded ISMN data in one of the supported formats in the past it can
be that station names are not recognized correctly because they contained the
'_' character which is supposed to be the separator. If you experience problems
because of this please download new data from the ISMN since this issue should
be fixed.

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
``source activate ismn``.

Guidelines
----------

If you want to contribute please follow these steps:

- Fork the ismn repository to your account
- Clone the repository
- make a new feature branch from the ismn master branch
- Add your feature
- Please include tests for your contributions in one of the test directories.
  We use py.test so a simple function called test_my_feature is enough
- submit a pull request to our master branch


Note
====

This project has been set up using PyScaffold 2.5.7. For details and usage
information on PyScaffold see http://pyscaffold.readthedocs.org/.
