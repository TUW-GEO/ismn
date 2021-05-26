# -*- coding: utf-8 -*-

from ismn.meta import MetaVar, MetaData
from ismn.components import Depth
import pytest
import unittest
import numpy as np

# todo: test negative depth
class Test_MetaVar(unittest.TestCase):
    def test_MetaVar(self):
        var = MetaVar("myvar", 1.1, Depth(0, 1))
        assert str(var) == "myvar (0.0 to 1.0 [m]): 1.1"
        assert tuple(var) == ("myvar", 1.1, 0, 1)
        assert var == var

        nvar = MetaVar("negmyvar", 1.1, Depth(0, -1))
        assert str(nvar) == "negmyvar (0.0 to -1.0 [m]): 1.1"
        assert tuple(nvar) == ("negmyvar", 1.1, -0, -1)
        assert nvar != var

        other = MetaVar("other", 99)
        assert str(other) == "other (no depth): 99"
        assert tuple(other) == ("other", 99, None, None)

        assert other != var


class Test_MetaData(unittest.TestCase):
    def setUp(self) -> None:
        vars = [
            MetaVar("first", 1, Depth(0, 1)),
            MetaVar("second", 0),
            MetaVar("neg", -99, Depth(-0.25, -1)),
            MetaVar("dup", "3rd", Depth(1, 3)),
        ]

        self.dat = MetaData(vars)

        self.other = MetaData([MetaVar("dup", "3rd", Depth(2, 4)), MetaVar("4", 4)])

    def test_format(self):
        df = self.dat.to_pd()
        assert df["first", "val"] == 1
        assert df["neg", "depth_from"] == -0.25
        assert df["dup", "depth_to"] == 3
        ddict = self.dat.to_dict()
        assert ddict["dup"] == [("3rd", 1, 3)]
        assert ddict["second"] == [(0, None, None)]

    def test_MetaData(self):
        assert len(self.dat) == 4
        assert "second" in self.dat
        assert self.dat[1] in self.dat
        assert tuple(self.dat[0]) == ("first", 1, 0.0, 1.0)
        assert tuple(self.dat[1]) == ("second", 0, None, None)

        assert self.dat["dup"] == self.dat[3]

        assert self.dat.keys() == ["first", "second", "neg", "dup"]

        assert (
            MetaData([MetaVar.from_tuple(("first", 1, 0.0, 1.0))])
            == self.dat[["first"]]
        )

    def test_best_meta(self):
        self.dat.merge(self.other, inplace=True)

        assert len(self.dat) == 6

        # no depths overlap
        best_meta_9_10 = self.dat.best_meta_for_depth(Depth(9, 10))
        assert sorted(best_meta_9_10.keys()) == sorted(["second", "4"])

        # all depths overlap
        best_meta_inf = self.dat.best_meta_for_depth(Depth(-np.inf, np.inf))
        assert len(best_meta_inf) == len(self.dat) - 1  # one duplicate removed
        assert sorted(best_meta_inf.keys()) == sorted(
            ["second", "4", "dup", "first", "neg"]
        )
        # both valzes for dup where equally good, so the first was kept
        assert best_meta_inf["dup"].depth.start == 1
        assert best_meta_inf["dup"].depth.end == 3

        # all but one dup and neg depth overlaps
        best_meta_015 = self.dat.best_meta_for_depth(Depth(0, 1.5))
        assert len(best_meta_015) == len(self.dat) - 2
        assert best_meta_015["dup"].depth.start == 1
        assert best_meta_015["dup"].depth.end == 3

        # both duplicate depths overlap, but one more --> keep second, drop neg
        best_meta_231 = self.dat.best_meta_for_depth(Depth(2, 3.1))
        assert (
            len(best_meta_231) == len(self.dat) - 3
        )  # one duplicate and first and neg
        assert best_meta_231["dup"].depth.start == 2
        assert best_meta_231["dup"].depth.end == 4

        # both duplicate depths overlap, equally good -> keep first
        best_meta_23 = self.dat.best_meta_for_depth(Depth(2, 3))
        assert len(best_meta_23) == len(self.dat) - 3  # one dup and first and neg
        assert best_meta_23["dup"].depth.start == 1.0
        assert best_meta_23["dup"].depth.end == 3.0

        # one matches perfectly
        best_meta_13 = self.dat.best_meta_for_depth(Depth(1, 3))
        assert len(best_meta_13) == len(self.dat) - 2  # one dup only, no neg
        assert best_meta_13["dup"].depth.start == 1.0
        assert best_meta_13["dup"].depth.end == 3.0

        # check with negative
        best_meta_neg = self.dat.best_meta_for_depth(Depth(-0.5, 2.0))
        # one dup was outside depth and is dropped, rest remains
        assert sorted(best_meta_neg.keys()) == sorted(
            ["first", "second", "dup", "4", "neg"]
        )
        assert best_meta_neg["dup"].depth.start == 1.0
        assert best_meta_neg["dup"].depth.end == 3.0

        # check with negative
        best_meta_only_neg = self.dat.best_meta_for_depth(Depth(-0.5, -1.0))
        # only keep meta without depths and for neg depth
        assert sorted(best_meta_only_neg.keys()) == sorted(["second", "neg", "4"])
        assert best_meta_only_neg["neg"].depth.start == -0.25
        assert best_meta_only_neg["neg"].depth.end == -1


if __name__ == "__main__":
    unittest.main()
