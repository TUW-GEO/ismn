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

'''
Created on Aug 1, 2013

@author: Christoph Paulik

Updated on Dec 14, 2018

@author: Philip Buttinger philip.buttinger@geo.tuwien.ac.at
'''


import os
import glob
import ismn.readers as readers
import numpy as np
import logging


def collect_from_folder(rootdir):
    """
    function walks the rootdir directory and looks for network
    folders and ISMN datafiles. It collects metadata for every
    file found and returns a numpy.ndarray of metadata

    Parameters
    ----------
    rootdir : string
        root directory on filesystem where the ISMN data was unzipped to

    Returns
    -------
    metadata : numpy.ndarray
        structured numpy array which contains the metadata for one file per row
    """

    logging.basicConfig(filename=os.path.join(rootdir, 'python_metadata', 'metadata.log'),
                        level=logging.DEBUG)

    metadata_catalog = []
    for root, subFolders, files in os.walk(rootdir):
        subFolders.sort()
        files.sort()

        # read additional metadata from csv file
        filename_csv = glob.glob('{}/*.csv'.format(root))
        # default values, if there is no csv file available or it crashes for e.g saturation
        lc_2000, lc_2005, lc_2010, lc_insitu, climate_KG, climate_insitu, saturation, clay_fraction, sand_fraction, \
        silt_fraction, organic_carbon = [np.nan, np.nan, np.nan, '', '', '', np.nan, np.nan, np.nan, np.nan, np.nan]
        if len(filename_csv) > 0:
            path_csv = os.path.join(root, filename_csv[0])
            try:
                lc_2000, lc_2005, lc_2010, lc_insitu, climate_KG, climate_insitu, saturation, clay_fraction, sand_fraction, \
                silt_fraction, organic_carbon = readers.get_metadata_from_csv(path_csv)
            except:
                logging.info('Error occured while reading metadata from csv file ({})'.format(root))
        else:
            if any(filename.endswith('.stm') for filename in files):
                logging.info('No csv file available ({})'.format(root))
            else:
                continue

        # print root,subFolders,files
        for filename in files:
            if filename.endswith('.stm'):
                fullfilename = os.path.join(root, filename)
                try:
                    metadata = readers.get_metadata(fullfilename)
                except (readers.ReaderException, IOError) as e:
                    continue

                for i, variable in enumerate(metadata['variable']):

                    metadata_catalog.append((metadata['network'], metadata['station'],
                                             variable, metadata['depth_from'][
                                                 i], metadata['depth_to'][i],
                                             metadata['sensor'], metadata[
                                                 'longitude'], metadata['latitude'],
                                             metadata['elevation'], fullfilename,
                                             lc_2000, lc_2005, lc_2010, lc_insitu, climate_KG, climate_insitu,
                                             saturation, clay_fraction, sand_fraction, silt_fraction, organic_carbon))

    return np.array(metadata_catalog, dtype=np.dtype([('network', object), ('station', object), ('variable', object),
                                                      ('depth_from', np.float), ('depth_to', np.float),
                                                      ('sensor', object), ('longitude', np.float),
                                                      ('latitude', np.float),
                                                      ('elevation', np.float), ('filename', object),
                                                      ('landcover_2000', object), ('landcover_2005', object),
                                                      ('landcover_2010', object), ('landcover_insitu', object),
                                                      ('climate', object), ('climate_insitu', object),
                                                      ('saturation', object),
                                                      ('clay_fraction', object), ('sand_fraction', object),
                                                      ('silt_fraction', object), ('organic_carbon', object)]))
