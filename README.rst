====
ismn
====

.. image:: https://travis-ci.org/TUW-GEO/ismn.svg?branch=master
    :target: https://travis-ci.org/TUW-GEO/ismn

.. image:: https://coveralls.io/repos/TUW-GEO/ismn/badge.png?branch=master
  :target: https://coveralls.io/r/TUW-GEO/ismn?branch=master

.. image:: https://badge.fury.io/py/ismn.svg
    :target: http://badge.fury.io/py/ismn

.. image:: https://readthedocs.org/projects/ismn/badge/?version=latest
   :target: http://ismn.readthedocs.org/

Readers for the data from the International Soil Moisture Database (ISMN).

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


Installation
============

This package should be installable through pip::

    pip install ismn

Setup Development environment
-----------------------------

1. Install Miniconda_. This will give you the ``conda`` command in your shell.
2. Run ``conda env create -f environment.yml`` this will install all the
   dependencies listed in the ``environment.yml`` file in this repository.
   By default this will create a new conda enviroment with the name ``ismn``.
   This can be changed by editing the ``environment.yml`` file.

.. _Miniconda: http://conda.pydata.org/miniconda.html

Example installation script
---------------------------

The following script will install miniconda and setup the environment on a UNIX
like system. Miniconda will be installed into ``$HOME/miniconda``.

::

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

After that you should be able to run::

    python setup.py test

to run the test suite.


Citation
========

If you use the software please cite it using the following DOI:

.. image:: https://zenodo.org/badge/101878880.svg
   :target: https://zenodo.org/badge/latestdoi/101878880

Note
====

This project has been set up using PyScaffold 2.5.7. For details and usage
information on PyScaffold see http://pyscaffold.readthedocs.org/.
