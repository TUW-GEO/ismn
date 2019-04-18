# The MIT License (MIT)
#
# Copyright (c) 2019 TU Wien
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""
This module tests the ISMN interface.
"""

import os
import unittest

from ismn.interface_v2 import IsmnFile, IsmnFileCollection


class IsmnFileTest(unittest.TestCase):

    def setUp(self):
        """
        Setup tests.
        """
        self.filename = os.path.join(
            os.path.dirname(
                __file__), 'test_data', 'format_ceop_sep', 'SMOSMANIA',
            'SMOSMANIA_SMOSMANIA_Narbonne_sm_0.050000_0.050000_ThetaProbe-ML2X_20070101_20070131.stm')

    def test_metadata(self):
        """
        Test
        """
        IsmnFile(self.filename)


class IsmnFileCollectionTest(unittest.TestCase):

    def setUp(self):
        """
        Setup tests.
        """
        self.path = os.path.join(os.path.dirname(__file__),
                                 'test_data', 'format_ceop_sep')

    def test_metadata(self):
        """
        Test
        """
        fc = IsmnFileCollection(self.path)
        fc.summary()


if __name__ == '__main__':
    unittest.main()
