# -*- coding: utf-8 -*-

from ismn.collections import IsmnFileCollection
import os
import pickle

def usecase_coll_zip():
    path = r"C:\Temp\delete_me\ismn\testdata_ceop.zip"
    coll = IsmnFileCollection(path, load_data=False)
    coll.store_metadata()
    nets = coll.get_networks()
    stats = coll.iter_stations(None)
    sens = coll.iter_sensors(nets[0], stats[0])


def usecase_coll_missing_file():
    path = r"C:\Temp\delete_me\risma"
    coll = IsmnFileCollection(path, load_data=False)
    nets = coll.get_networks()
    stats = coll.iter_stations(None)
    sens = coll.iter_sensors(nets[0], stats[0])

def usecase_coll_large():
    path = r"D:\data-read\ISMN\global_20191024"
    coll = IsmnFileCollection(path, load_data=False)
    coll.store_metadata()
    nets = coll.get_networks()
    stats = coll.iter_stations(None)
    sens = coll.iter_sensors()

def usecase_coll_small():
    path = r"C:\Temp\delete_me\ismn\testdata_ceop"
    coll = IsmnFileCollection(path, load_data=False)

if __name__ == '__main__':
    usecase_coll_small()
