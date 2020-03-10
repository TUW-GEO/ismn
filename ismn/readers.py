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
Created on Jul 31, 2013

@author: Christoph Paulik

Updated on Dec 14, 2018

@author: Philip Buttinger philip.buttinger@geo.tuwien.ac.at
'''

import os
import pandas as pd
from datetime import datetime
import numpy as np
import logging
import io


variable_lookup = {'sm': 'soil moisture',
                   'ts': 'soil temperature',
                   'su': 'soil suction',
                   'p': 'precipitation',
                   'ta': 'air temperature',
                   'fc': 'field capacity',
                   'wp': 'permanent wilting point',
                   'paw': 'plant available water',
                   'ppaw': 'potential plant available water',
                   'sat': 'saturation',
                   'si_h': 'silt fraction',
                   'sd': 'snow depth',
                   'sa_h': 'sand fraction',
                   'cl_h': 'clay fraction',
                   'oc_h': 'organic carbon',
                   'sweq': 'snow water equivalent',
                   'tsf': 'surface temperature',
                   'tsfq': 'surface temperature quality flag original'
                   }


class ReaderException(Exception):
    pass


class ISMNTSError(Exception):
    pass


class ISMNTimeSeries(object):
    """
    class that contains a time series of ISMN data read from one text file

    Attributes
    ----------
    network : string
        network the time series belongs to
    station : string
        station name the time series belongs to
    latitude : float
        latitude of station
    longitude : float
        longitude of station
    elevation : float
        elevation of station
    variable : list
        variable measured
    depth_from : list
        shallower depth of layer the variable was measured at
    depth_to : list
        deeper depth of layer the variable was measured at
    sensor : string
        sensor name
    data : pandas.DataFrame
        data of the time series
    """

    def __init__(self, data):

        for key in data:
            setattr(self, key, data[key])

    def __repr__(self):

        return '%s %s %.2f m - %.2f m %s measured with %s ' % (
            self.network,
            self.station,
            self.depth_from[0],
            self.depth_to[0],
            self.variable[0],
            self.sensor)

    def plot(self, *args, **kwargs):
        """
        wrapper for pandas.DataFrame.plot which adds title to plot
        and drops NaN values for plotting
        Returns
        -------
        ax : axes
            matplotlib axes of the plot

        Raises
        ------
        ISMNTSError
            if data attribute is not a pandas.DataFrame
        """
        if type(self.data) is pd.DataFrame:
            tempdata = self.data.dropna()
            tempdata = tempdata[tempdata.columns[0]]
            ax = tempdata.plot(*args, figsize=(15, 5), **kwargs)
            ax.set_title(self.__repr__())
            return ax
        else:
            raise ISMNTSError("data attribute is not a pandas.DataFrame")


def get_info_from_file(filename):
    """
    reads first line of file and splits filename
    this can be used to construct necessary metadata information
    for all ISMN formats

    Parameters
    ----------
    filename : string
        filename including path

    Returns
    -------
    header_elements : list
        first line of file split into list
    filename_elements : list
        filename without path split by _
    """
    with io.open(filename, mode='r', newline=None) as f:
        header = f.readline()
    header_elements = header.split()

    path, filen = os.path.split(filename)
    filename_elements = filen.split('_')

    return header_elements, filename_elements


def get_metadata_header_values(filename):
    """
    get metadata from ISMN textfiles in the format called
    Variables stored in separate files (CEOP formatted)

    Parameters
    ----------
    filename : string
        path and name of file

    Returns
    -------
    metadata : dict
        dictionary of metadata information
    """

    header_elements, filename_elements = get_info_from_file(filename)

    if len(filename_elements) > 9:
        sensor = '_'.join(filename_elements[6:len(filename_elements) - 2])
    else:
        sensor = filename_elements[6]

    if filename_elements[3] in variable_lookup:
        variable = [variable_lookup[filename_elements[3]]]
    else:
        variable = [filename_elements[3]]

    metadata = {'network': header_elements[1],
                'station': header_elements[2],
                'latitude': float(header_elements[3]),
                'longitude': float(header_elements[4]),
                'elevation': float(header_elements[5]),
                'depth_from': [float(header_elements[6])],
                'depth_to': [float(header_elements[7])],
                'variable': variable,
                'sensor': sensor}

    return metadata


def read_format_header_values(filename):
    """
    Reads ISMN textfiles in the format called
    Variables stored in separate files (Header + values)

    Parameters
    ----------
    filename : string
        path and name of file

    Returns
    -------
    time_series : ISMNTimeSeries
        ISMNTimeSeries object initialized with metadata and data from file
    """

    metadata = get_metadata_header_values(filename)

    data = pd.read_csv(filename, skiprows=1, delim_whitespace=True,
                       names=['date', 'time', metadata['variable'][0],
                              metadata['variable'][0] + '_flag',
                              metadata['variable'][0] + '_orig_flag'],
                       parse_dates=[[0, 1]])

    data.set_index('date_time', inplace=True)
    data = data.tz_localize('UTC')

    metadata['data'] = data

    return ISMNTimeSeries(metadata)


def get_metadata_ceop_sep(filename):
    """
    get metadata from ISMN textfiles in the format called
    Variables stored in separate files (CEOP formatted)

    Parameters
    ----------
    filename : string
        path and name of file

    Returns
    -------
    metadata : dict
        dictionary of metadata information
    """

    header_elements, filename_elements = get_info_from_file(filename)

    if len(filename_elements) > 9:
        sensor = '_'.join(filename_elements[6:len(filename_elements) - 2])
    else:
        sensor = filename_elements[6]

    if filename_elements[3] in variable_lookup:
        variable = [variable_lookup[filename_elements[3]]]
    else:
        variable = [filename_elements[3]]

    metadata = {'network': filename_elements[1],
                'station': filename_elements[2],
                'variable': variable,
                'depth_from': [float(filename_elements[4])],
                'depth_to': [float(filename_elements[5])],
                'sensor': sensor,
                'latitude': float(header_elements[7]),
                'longitude': float(header_elements[8]),
                'elevation': float(header_elements[9])
                }

    return metadata


def read_format_ceop_sep(filename):
    """
    Reads ISMN textfiles in the format called
    Variables stored in separate files (CEOP formatted)

    Parameters
    ----------
    filename : string
        path and name of file

    Returns
    -------
    time_series : ISMNTimeSeries
        ISMNTimeSeries object initialized with metadata and data from file
    """

    metadata = get_metadata_ceop_sep(filename)

    data = pd.read_csv(filename, delim_whitespace=True, usecols=[0, 1, 12, 13, 14],
                       names=['date', 'time',
                              metadata['variable'][0],
                              metadata['variable'][0] + '_flag',
                              metadata['variable'][0] + '_orig_flag'],
                       parse_dates=[[0, 1]])

    data.set_index('date_time', inplace=True)
    data = data.tz_localize('UTC')

    metadata['data'] = data

    return ISMNTimeSeries(metadata)


def get_metadata_ceop(filename):
    """
    get metadata from ISMN textfiles in the format called
    CEOP Reference Data Format

    Parameters
    ----------
    filename : string
        path and name of file

    Returns
    -------
    metadata : dict
        dictionary of metadata information
    """

    header_elements, filename_elements = get_info_from_file(filename)

    metadata = {'network': filename_elements[1],
                'station': header_elements[6],
                'variable': ['ts', 'sm'],
                'sensor': 'n.s',
                'depth_from': ['multiple'],
                'depth_to': ['multiple'],
                'latitude': float(header_elements[7]),
                'longitude': float(header_elements[8]),
                'elevation': float(header_elements[9])
                }

    return metadata


def read_format_ceop(filename):
    """
    Reads ISMN textfiles in the format called
    CEOP Reference Data Format

    Parameters
    ----------
    filename : string
        path and name of file

    Returns
    -------
    time_series : ISMNTimeSeries
        ISMNTimeSeries object initialized with metadata and data from file
    """
    metadata = get_metadata_ceop(filename)
    data = pd.read_csv(filename, delim_whitespace=True, usecols=[0, 1, 11, 12, 13, 14, 15],
                       names=['date', 'time', 'depth_from',
                              metadata['variable'][0],
                              metadata['variable'][0] + '_flag',
                              metadata['variable'][1],
                              metadata['variable'][1] + '_flag'],
                       na_values=['-999.99'],
                       parse_dates=[[0, 1]])

    data = data.set_index('date_time')
    data = data.tz_localize('UTC')
    date_index = data.index
    depth_index = data['depth_from']

    del data['depth_from']

    data.index = pd.MultiIndex.from_arrays([depth_index,
                                            depth_index.rename('depth_to'),
                                            date_index])
    data.index.names = ['depth_from', 'depth_to', 'date']

    data = data.sort_index(level=0)

    metadata['depth_from'] = np.unique(
        data.index.get_level_values(0).values).tolist()
    metadata['depth_to'] = np.unique(
        data.index.get_level_values(1).values).tolist()
    metadata['data'] = data

    return ISMNTimeSeries(metadata)


def tail(f, lines=1, _buffer=4098):
    """Tail a file and get X lines from the end

    Parameters
    ----------
    f: file like object
    lines: int
       lines from the end of the file to read
    _buffer: int
       buffer to use to step backwards in the file.

    References
    ----------
    Found at http://stackoverflow.com/a/13790289/1314882
    """
    # place holder for the lines found
    lines_found = []

    # block counter will be multiplied by buffer
    # to get the block size from the end
    block_counter = -1

    # loop until we find X lines
    while len(lines_found) < lines:
        try:
            f.seek(block_counter * _buffer, os.SEEK_END)
        except IOError:  # either file is too small, or too many lines requested
            f.seek(0)
            lines_found = f.readlines()
            break

        lines_found = f.readlines()

        # we found enough lines, get out
        if len(lines_found) > lines:
            break

        # decrement the block counter to get the
        # next X bytes
        block_counter -= 1

    return lines_found[-lines:]


def get_min_max_timestamp_header_values(filename):
    """
    Get minimum and maximum observation timestamp from header values format.
    """
    with io.open(filename, mode='r', newline=None) as fid:
        _ = fid.readline()
        first = fid.readline()
        last = tail(fid)[0]

    min_date = datetime.strptime(first[:16], '%Y/%m/%d %H:%M')
    max_date = datetime.strptime(last[:16], '%Y/%m/%d %H:%M')
    return min_date, max_date


def get_min_max_timestamp_ceop_sep(filename):
    """
    Get minimum and maximum observation timestamp from ceop_sep format.
    """
    with io.open(filename, mode='r', newline=None) as fid:
        first = fid.readline()
        last = tail(fid)[0]

    min_date = datetime.strptime(first[:16], '%Y/%m/%d %H:%M')
    max_date = datetime.strptime(last[:16], '%Y/%m/%d %H:%M')
    return min_date, max_date


def get_min_max_timestamp_ceop(filename):
    """
    Get minimum and maximum observation timestamp from ceop format.
    """
    with io.open(filename, mode='r', newline=None) as fid:
        first = fid.readline()
        last = tail(fid)[0]

    min_date = datetime.strptime(first[:16], '%Y/%m/%d %H:%M')
    max_date = datetime.strptime(last[:16], '%Y/%m/%d %H:%M')
    return min_date, max_date


def get_min_max_timestamp(filename):
    """
    Determine the file type and get the minimum and maximum observation
    timestamp

    """
    dicton = globals()
    func = dicton['get_min_max_timestamp_' + get_format(filename)]
    return func(filename)


def get_format(filename):
    """
    get's the file format from the length of
    the header and filename information

    Parameters
    ----------
    filename : string

    Returns
    -------
    methodname : string
        name of method used to read the detected format

    Raises
    ------
    ReaderException
        if filename or header parts do not fit one of the formats
    """
    header_elements, filename_elements = get_info_from_file(filename)
    if len(filename_elements) == 5 and len(header_elements) == 16:
        return 'ceop'
    if len(header_elements) == 15 and len(filename_elements) >= 9:
        return 'ceop_sep'
    if len(header_elements) < 14 and len(filename_elements) >= 9:
        return 'header_values'
    raise ReaderException(
        "This does not seem to be a valid ISMN filetype %s" % filename)


def read_data(filename):
    """
    reads ISMN data in any format

    Parameters
    ----------
    filename: string

    Returns
    -------
    timeseries: IMSNTimeSeries
    """
    dicton = globals()
    func = dicton['read_format_' + get_format(filename)]
    return func(filename)


def get_metadata_from_csv(filename):
    """
    reads ISMN metadata from csv file

    Parameters
    ----------
    filename: str, path to csv file

    Returns
    -------
    landcover_2000: int, cci landcover classification for station (year 2000)
    landcover_2005: int, cci landcover classification for station (year 2005)
    landcover_2010: int, cci landcover classification for station (year 2010)
    landcover_insitu: str, in situ landcover classification
    climate: str, Koeppen Geiger climate classification for station
    climate_insitu: str, in situ climate classification for station
    saturation: nd.array, saturation for all available depths
    clay_fraction: nd.array, clay fraction for all available depths (in % weight)
    sand_fraction: nd.array, sand fraction for all available depths (in % weight)
    silt_fraction: nd.array, silt fraction for all available depths (in % weight)
    organic_carbon: nd.array, organic carbon for all available depths (in % weight)
    """
    def read_field(fieldname):
        if fieldname in data.index:
            dt = list()
            for i, j in zip(np.atleast_1d(data.loc[fieldname]['depth_from[m]']),
                            np.atleast_1d(data.loc[fieldname]['depth_to[m]'])):
                dt.append(('{}m_{}m'.format(i, j), np.float))
            return np.array([tuple(np.atleast_1d(data.loc[fieldname]['value']))], dtype=np.dtype(dt))
        else:
            return np.nan

    # some stations don't come with correct format in csv file (missing header)
    try:
        data = pd.read_csv(filename, delimiter=";")
        data.set_index('quantity_name', inplace=True)
    except:
        # set columns manually
        logging.info('no header: {}'.format(filename))
        data = pd.read_csv(filename, delimiter=";", header=None)
        cols = list(data.columns.values)
        cols[0:7] = ['quantity_name', 'unit', 'depth_from[m]', 'depth_to[m]',
                     'value', 'description', 'quantity_source_name']
        data.columns = cols
        data.set_index('quantity_name', inplace=True)

    # read landcover classifications
    lc = data.loc[['land cover classification']][['value', 'quantity_source_name']]
    lc_dict = {'CCI_landcover_2000': np.nan, 'CCI_landcover_2005': np.nan,
               'CCI_landcover_2010': np.nan, 'insitu': ''}
    for key in lc_dict.keys():
        if key in lc['quantity_source_name'].values:
            if key != 'insitu':
                lc_dict[key] = np.int(lc.loc[lc['quantity_source_name'] == key]['value'].values[0])
            else:
                lc_dict[key] = lc.loc[lc['quantity_source_name'] == key]['value'].values[0]
                logging.info('insitu land cover classification available: {}'.format(filename))

    # read climate classifications
    cl = data.loc[['climate classification']][['value', 'quantity_source_name']]
    cl_dict = {'koeppen_geiger_2007': '', 'insitu': ''}
    for key in cl_dict.keys():
        if key in cl['quantity_source_name'].values:
            cl_dict[key] = cl.loc[cl['quantity_source_name'] == key]['value'].values[0]
            if key == 'insitu':
                logging.info('insitu climate classification available: {}'.format(filename))

    saturation = read_field('saturation')
    clay_fraction = read_field('clay fraction')
    sand_fraction = read_field('sand fraction')
    silt_fraction = read_field('silt fraction')
    organic_carbon = read_field('organic carbon')

    return lc_dict['CCI_landcover_2000'], lc_dict['CCI_landcover_2005'], lc_dict['CCI_landcover_2010'], \
           lc_dict['insitu'], cl_dict['koeppen_geiger_2007'], cl_dict['insitu'], saturation, clay_fraction, \
           sand_fraction, silt_fraction, organic_carbon


def get_metadata(filename):
    """
    reads ISMN metadata from any format

    Parameters
    ----------
    filename: string

    Returns
    -------
    metadata: dict
    """
    dicton = globals()
    func = dicton['get_metadata_' + get_format(filename)]
    return func(filename)
