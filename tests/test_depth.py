# -*- coding: utf-8 -*-

import unittest
from ismn.const import *
from ismn.meta import Depth


class DepthTest(unittest.TestCase):
    def setUp(self):
        """
        Setup test data.
        """
        self.d = Depth(0, 0.05)
        assert str(self.d) == "0.0 to 0.05 [m]"
        assert tuple(self.d) == (0.0, 0.05)

    def test_neg_depth(self):
        other = Depth(0.0, -0.05)
        assert str(other) == "0.0 to -0.05 [m]"

        assert self.d != other
        assert other.encloses(self.d) == False
        assert other.enclosed(self.d) == False
        assert self.d.across0 == other.across0 == False
        assert self.d.is_profile == other.is_profile == True

        assert other.perc_overlap(self.d) == 0.0

    def test_attributes(self):
        """
        Test depth attributes.
        """
        assert self.d.start == 0
        assert self.d.end == 0.05

    def test_is_profile(self):
        """
        Test if depth represents a profile.
        """
        assert self.d.is_profile

    def test_equal(self):
        """
        Test depth equality.
        """
        other = Depth(0, 0.05)
        assert self.d == other

    def test_invalid(self):
        try:
            Depth(0.5, 0.1)
            raise AssertionError
        except DepthError:
            pass
        try:
            Depth(-0.5, -0.1)
            raise AssertionError
        except DepthError:
            pass

    def test_enclose(self):
        """
        Test if other depth encloses depth.
        """
        other = Depth(0, 0.05)
        assert self.d.encloses(other)
        assert other.enclosed(self.d)

        other = Depth(0, 0.1)
        assert not self.d.encloses(other)
        assert not other.enclosed(self.d)
        assert other.across0 == False

        other = Depth(-0.1, -0.2)
        assert not self.d.encloses(other)
        assert not self.d.enclosed(other)
        assert other.across0 == False

    def test_perc_overlap(self):
        other = Depth(0.05, 0.1)
        assert other.across0 == False
        assert self.d.overlap(other) == True
        assert self.d.perc_overlap(other) == 0.0

        other = Depth(0.03, 0.05)
        assert other.across0 == False
        assert self.d.overlap(other) == True
        assert self.d.perc_overlap(other) == round(0.02 / 0.05, 7)

        other = Depth(0, 0.05)
        assert other.across0 == False
        assert self.d.overlap(other) == True
        assert self.d.perc_overlap(other) == 1.0

        other = Depth(-0.01, -0.05)
        assert other.across0 == False
        assert self.d.overlap(other) == False
        assert self.d.perc_overlap(other) == -1

        other = Depth(-0.01, 0.01)
        assert other.across0 == True
        assert self.d.overlap(other) == True
        assert self.d.perc_overlap(other) == round(0.01 / 0.06, 7)


if __name__ == "__main__":
    unittest.main()
