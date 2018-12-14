# Copyright (c) 2018, TU Wien, Department of Geodesy and Geoinformation
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#    * Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above copyright
#      notice, this list of conditions and the following disclaimer in the
#      documentation and/or other materials provided with the distribution.
#    * Neither the name of TU Wien, Department of Geodesy and Geoinformation nor the
#      names of its contributors may be used to endorse or promote products
#      derived from this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
