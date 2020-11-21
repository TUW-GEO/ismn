# -*- coding: utf-8 -*-

from ismn.meta import MetaVar, MetaData
from ismn.components import Depth
import pytest
import unittest
import numpy as np

def test_MetaVar():
    var = MetaVar('myvar', 1.1, Depth(0, 1))
    assert str(var) == "myvar (0_1[m]): 1.1"
    assert tuple(var) == ('myvar', 1.1, 0, 1)
    assert var == var

    other = MetaVar('other', 99)
    assert str(other) == "other (no depth): 99"
    assert tuple(other) == ('other', 99, None, None)

    assert other != var

class Test_MetaData(unittest.TestCase):
    def setUp(self) -> None:
        vars = [MetaVar('first', 1, Depth(0, 1)),
                MetaVar('second', 0),
                MetaVar('dup', '3rd', Depth(1, 3))]

        self.dat = MetaData(vars)

        self.other = MetaData([MetaVar('dup', '3rd', Depth(2,4)),
                               MetaVar('4', 4)])
        
    def test_MetaData(self):
        assert len(self.dat) == 3
        assert 'second' in self.dat
        assert self.dat[1] in self.dat
        assert tuple(self.dat[0]) == ('first', 1, 0., 1.)
        assert tuple(self.dat[1]) == ('second', 0, None, None)

        assert str(self.dat) == 'first (0.0_1.0[m]), second (no depth), dup (1.0_3.0[m])'

        assert self.dat['dup'] == self.dat[2]

        assert self.dat.keys() == ['first', 'second', 'dup']

        assert MetaData.from_dict({'first': (1,0,1)}) == self.dat[['first']]


    def test_best_meta(self):
        self.dat.merge(self.other, inplace=True)

        assert len(self.dat) == 5

        with pytest.raises(ValueError):
            self.dat.to_pd() # conversion to dict with two vars of same name not supported

        # no depths overlap
        best_meta_9_10 = self.dat.best_meta_for_depth(Depth(9,10))
        assert sorted(best_meta_9_10.keys()) == sorted(['second', '4'])

        # all depths overlap
        best_meta_inf = self.dat.best_meta_for_depth(Depth(-np.inf,np.inf))
        assert len(best_meta_inf) == len(self.dat)-1 # one duplicate removed
        assert sorted(best_meta_inf.keys()) == sorted(['second', '4', 'dup', 'first'])
        # both valzes for dup where equally good, so the first was kept
        assert best_meta_inf['dup'].depth.start == 1
        assert best_meta_inf['dup'].depth.end == 3

        # all but one depth overlaps
        best_meta_015 = self.dat.best_meta_for_depth(Depth(0,1.5))
        assert len(best_meta_015) == len(self.dat)-1
        assert best_meta_015['dup'].depth.start == 1
        assert best_meta_015['dup'].depth.end == 3

        # both duplicate depths overlap, but one more --> keep second
        best_meta_231 = self.dat.best_meta_for_depth(Depth(2,3.1))
        assert len(best_meta_231) == len(self.dat)-2 # one duplicate and first
        assert best_meta_231['dup'].depth.start == 2
        assert best_meta_231['dup'].depth.end == 4

        # both duplicate depths overlap, equally good -> keep first
        best_meta_23 = self.dat.best_meta_for_depth(Depth(2,3))
        assert len(best_meta_23) == len(self.dat)-2 # one dup and first
        assert best_meta_23['dup'].depth.start == 1.
        assert best_meta_23['dup'].depth.end == 3.
        
        # one matches perfectlyt
        best_meta_13 = self.dat.best_meta_for_depth(Depth(1,3))
        assert len(best_meta_13) == len(self.dat)-1 # one dup only
        assert best_meta_13['dup'].depth.start == 1.
        assert best_meta_13['dup'].depth.end == 3.
        


if __name__ == '__main__':
    test_MetaVar()

