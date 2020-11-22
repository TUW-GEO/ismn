# -*- coding: utf-8 -*-
from pkg_resources import get_distribution, DistributionNotFound

try:
    # Change here if project is renamed and does not equal the package name
    dist_name = __name__
    __version__ = get_distribution(dist_name).version
except DistributionNotFound:
    __version__ = 'unknown'
finally:
    del get_distribution, DistributionNotFound

import os

# path into src
src_path = os.path.join(os.path.dirname(__file__), '..')

from ismn.network_collection import NetworkCollection

# path to tests dir, if it exists
tests_path = os.path.join(src_path, '..', 'tests')
if not os.path.exists(tests_path):
    tests_path = 'unknown'

# path to testsdata dir, if it exists
testdata_path = os.path.join(src_path, '..', 'tests', 'test_data')
if not os.path.exists(testdata_path):
    testdata_path = 'unknown'

