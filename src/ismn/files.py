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

import os
import pandas as pd
from copy import deepcopy

from collections import OrderedDict
from ismn.base import IsmnRoot
from ismn.components import *
from ismn import tables
from ismn.meta import MetaVar, MetaData

from tempfile import gettempdir, TemporaryDirectory
from pathlib import Path, PurePosixPath
import warnings

class IsmnFileError(IOError):
    pass

class IsmnFile(object):
    """
    General base class for data and static metadata files in ismn archive.

    Parameters
    ----------
    root: IsmnRoot or str
        Base archive that contains the file to read
    file_path : Path or str
        Path to the file in the archive.
    temp_root : Path or str, optional (default : gettempdir())
        Root directory where a separate subdir for temporary files
        will be created (and deleted).
    """

    def __init__(self, root, file_path, temp_root=gettempdir()):

        if not isinstance(root, IsmnRoot):
            root = IsmnRoot(root)

        self.root = root
        self.file_path = self.root._clean_subpath(file_path)

        if self.file_path not in self.root:
            raise IOError(f'Archive does not contain file: {self.file_path}')

        if not os.path.exists(temp_root):
            os.makedirs(temp_root, exist_ok=True)

        self.temp_root = temp_root

    def close(self):
        self.root.close()

    def open(self):
        self.root.open()


class StaticMetaFile(IsmnFile):
    """
    Represents a csv file containing site specific static variables.
    These attributes shall be assigned to all sensors at that site.

    Parameters
    ----------
    root: IsmnRoot or str
        Archive that contains the file to read
    file_path : Path or str
        Subpath to the file in the root. No leading slash!
    temp_root : Path or str, optional (default : gettempdir())
        Root directory where a separate subdir for temporary files
        will be created (and deleted).
    """

    def __init__(self, root, file_path, temp_root=gettempdir()):

        super(StaticMetaFile, self).__init__(root, file_path, temp_root)

        if self.file_path.suffix.lower() != '.csv':
            raise IsmnFileError(f'CSV file expected for StaticMetaFile object')

    def _read_field(self, fieldname:str, new_name=None) -> np.array:
        """
        Extract a field from the loaded csv metadata
        """
        field_vars = []

        if fieldname in self.data.index:

            froms = np.atleast_1d(self.data.loc[fieldname]['depth_from[m]'])
            tos = np.atleast_1d(self.data.loc[fieldname]['depth_to[m]'])
            vals = np.atleast_1d(self.data.loc[fieldname]['value'])

            for d_from, d_to, val in zip(froms, tos, vals):
                d = Depth(d_from, d_to)
                name = new_name if new_name is not None else fieldname
                try:
                    val = float(val)
                except ValueError:
                    pass # value is actually a string, that's ok
                field_vars.append(MetaVar(name, val, d))

        return field_vars

    def _read_csv(self, csvfile:Path) -> pd.DataFrame:
        """ Load static metadata data frame from csv """
        try:
            data = pd.read_csv(csvfile, delimiter=";")
            data.set_index('quantity_name', inplace=True)
        except:
            # set columns manually
            logging.info('no header: {}'.format(csvfile))
            data = pd.read_csv(csvfile, delimiter=";", header=None)
            cols = list(data.columns.values)
            cols[:len(tables.CSV_COLS)] = tables.CSV_COLS # todo: not safe
            data.columns = cols
            data.set_index('quantity_name', inplace=True)

        return data

    def read_metadata(self):
        """
        Read csv file containing static variables into data frame.

        Returns
        -------
        data : MetaData
            Data read from csv file.
        """
        if self.root.zip:
            if not self.root.isopen: self.root.open()
            with TemporaryDirectory(prefix='ismn', dir=self.temp_root) as tempdir:
                extracted = self.root.extract_file(self.file_path, tempdir)
                self.data = self._read_csv(extracted)
        else:
            self.data = self._read_csv(self.root.path / self.file_path)

        # read landcover classifications
        lc = self.data.loc[['land cover classification']][['value', 'quantity_source_name']]

        lc_dict = {'CCI_landcover_2000': tables.CSV_META_TEMPLATE['lc_2000'],
                   'CCI_landcover_2005': tables.CSV_META_TEMPLATE['lc_2005'],
                   'CCI_landcover_2010': tables.CSV_META_TEMPLATE['lc_2010'],
                   'insitu': tables.CSV_META_TEMPLATE['lc_insitu']}

        cl_dict = {'koeppen_geiger_2007': tables.CSV_META_TEMPLATE['climate_KG'],
                   'insitu': tables.CSV_META_TEMPLATE['climate_insitu']}

        for key in lc_dict.keys():
            if key in lc['quantity_source_name'].values:
                if key != 'insitu':
                    lc_dict[key] = np.int(lc.loc[lc['quantity_source_name']
                                                 == key]['value'].values[0])
                else:
                    lc_dict[key] = lc.loc[lc['quantity_source_name']
                                          == key]['value'].values[0]
                    logging.info(f'insitu land cover classification available: {self.file_path}')

        # read climate classifications
        cl = self.data.loc[['climate classification']][['value', 'quantity_source_name']]
        for key in cl_dict.keys():
            if key in cl['quantity_source_name'].values:
                cl_dict[key] = cl.loc[cl['quantity_source_name'] == key]['value'].values[0]
                if key == 'insitu':
                    logging.info(f'insitu climate classification available: {self.file_path}')

        metavars = []

        metavars.append(MetaVar('lc_2000', lc_dict['CCI_landcover_2000']))
        metavars.append(MetaVar('lc_2005', lc_dict['CCI_landcover_2005']))
        metavars.append(MetaVar('lc_2010', lc_dict['CCI_landcover_2010']))
        metavars.append(MetaVar('lc_insitu', lc_dict['insitu']))

        metavars.append(MetaVar('climate_KG', cl_dict['koeppen_geiger_2007']))
        metavars.append(MetaVar('climate_insitu', cl_dict['insitu']))

        static_meta = {
            'saturation': self._read_field('saturation'),
            'clay_fraction': self._read_field('clay fraction', new_name='clay_fraction'),
            'sand_fraction': self._read_field('sand fraction', new_name='sand_fraction'),
            'silt_fraction': self._read_field('silt fraction', new_name='silt_fraction'),
            'organic_carbon': self._read_field('organic carbon', new_name='organic_carbon')}
        for name, vars in static_meta.items():
            if len(vars) > 0:
                metavars += vars
            else:
                metavars.append(MetaVar(name, tables.CSV_META_TEMPLATE[name]))

        return MetaData(metavars)


class DataFile(IsmnFile):

    """
    IsmnFile class represents a single ISMN data file.
    This represents only .stm data files not metadata csv files.

    Parameters
    ----------
    root : IsmnRoot or str
        Archive to the downloaded data.
    file_path : str
        Path in the archive to the ismn file. No leading slash!
    load_data : bool, optional
        If True data will be loaded during metadata reading.
    load_metadata : bool, optional (default: True)
        Load metadata during initialisation.
    static_meta : OrderedDict, optional (default: None)
        If the static meta for the file has been read before, the OrderedDict
        returned by StaticMetaFile.read_metadata() can be passed here directly.
        This can be used to avoid reading the same static meta file e.g for
        multiple sensors at a station. By the default, the static_meta is loaded
        from the according csv file for the passed data file.
    temp_root : Path or str, optional (default : gettempdir())
        Root directory where a separate subdir for temporary files
        will be created (and deleted).

    Attributes
    ----------
    filename : str
        Filename.
    file_type : str
        File type information (e.g. ceop).
    metadata : dict
        Metadata information.
    data : numpy.ndarray
        Data stored in file.
    static_meta : OrderedDict
        Static meta data loaded from station csv file.
    # todo: update attrs and methods

    Methods
    -------
    check_metadata(self, variable, min_depth=0, max_depth=0.1, filter_static_vars=None)
        Evaluate whether the file complies with the passed metadata requirements
    load_data()
        Load data from file.
    read_data()
        Read data in file.
    read_metadata()
        Read metadata from file name and first line of file.
    _get_metadata_ceop_sep()
        Get metadata in the file format called CEOP in separate files.
    _get_metadata_header_values()
        Get metadata file in the format called Header Values.
    _get_metadata_from_file(delim='_')
        Read first line of file and split filename.
        Information is used to collect metadata information for all
        ISMN formats.
    _read_format_ceop_sep()
        Read data in the file format called CEOP in separate files.
    _read_format_header_values()
        Read data file in the format called Header Values.
    _read_csv(names=None, usecols=None, skiprows=0)
        Read data.
    """

    def __init__(self, root, file_path, load_data=False,
                 load_metadata=True, static_meta=None,
                 temp_root=gettempdir()):

        super(DataFile, self).__init__(root, file_path, temp_root)

        self.file_type = 'undefined'

        self.metadata = {}
        if load_metadata:
            self.metadata = self.read_metadata(static_meta=static_meta,
                                               best_meta_for_sensor=True)

        self.data = None
        if load_data:
            self.load_data()

    def __getitem__(self, item):
        return self.metadata[item]

    def load_data(self):
        """
        Load data from file.
        """
        if self.data is None:

            if not self.root.isopen: self.open()

            if self.file_type == 'ceop':
                # todo: what is this format?
                # self._read_format_ceop()
                raise NotImplementedError
            elif self.file_type == 'ceop_sep':
                self._read_format_ceop_sep()
            elif self.file_type == 'header_values':
                self._read_format_header_values()
            else:
                raise IOError(f"Unknown file format found for: {self.file_path}")
                # logger.warning(f"Unknown file type: {self.file_path}")

    def read_data(self):
        """
        Read data in file. Load file if necessary.

        Returns
        -------
        data : pd.DataFrame
            File content.
        """
        self.load_data()
        return self.data

    def check_metadata(self, variable, min_depth=0, max_depth=0.1,
                        filter_static_vars=None) -> bool:
        """
        Evaluate whether the file complies with the passed metadata requirements

        Parameters
        ----------
        variable : str
            Name of the required variable measured, e.g. soil_moisture
        min_depth : float, optional (default: 0)
            Minimum depth that the measurement should have.
        max_depth : float, optional (default: 0.1)
            Maximum depth that the measurement should have.
        filter_static_vars: dict
            Additional metadata keys and values for which the file list is filtered
            e.g. {'lc_2010': 10} to filter for a landcover class.

        Returns
        -------
        valid : bool
            Whether the metadata complies with the passed conditions or not.
        """

        if min_depth is None:
            min_depth = -np.inf
        if max_depth is None:
            max_depth = np.info

        lc_cl = list(tables.CSV_META_TEMPLATE.keys())

        if not (self.metadata['variable'].val == variable):
            return False

        sensor_depth = self.metadata['sensor'].depth

        if not Depth(min_depth, max_depth).encloses(sensor_depth):
            return False

        if filter_static_vars:
            fil_lc_cl = [True]
            for k in filter_static_vars.keys():
                if k not in lc_cl:
                    raise ValueError(f"{k} is not a valid metadata variable, "
                                     f"select one of {lc_cl}")
                fil_lc_cl.append(self.metadata[k].val == filter_static_vars[k])

            if not all(fil_lc_cl):
                return False

        return True

    def get_formatted_metadata(self, format='list'):
        """
        Return metadata for file in different formats such as list, dataframe,
        dict and structured array.

        Parameters
        ----------
        format : str, optional (default: 'list')
            Name of the format, one of 'struct', 'dict', 'pandas', 'list'

        Returns
        -------
        formatted_meta : Any
            Metadata in the selected format.
        """

        if self.metadata is None:
            warnings.warn("No metadata yet loaded.")
            return None

        meta = deepcopy(self.metadata)

        depth = meta.pop('depth')

        meta['depth_from'], meta['depth_to'] = depth.start, depth.end

        meta['archive'] = self.root.path
        meta['filepath'] = self.file_path

        if format.lower() == 'struct':
            return np.array([tuple(meta.values())],
                            [(k, object) for k in meta.keys()])
        elif format.lower() == 'dict':
            return meta
        elif format.lower() == 'pandas':
            return pd.Series(meta)
        elif format.lower() == 'list':
            return tuple(list(meta.values()))
        else:
            raise NotImplementedError(f"Format {format} is not (yet) implemented")

    def read_metadata(self, static_meta=None, best_meta_for_sensor=True):
        """
        Read metadata from file name and first line of file.

        Parameters
        ----------
        static_meta : MetaData, optional (default: None)
            Static meta data for the file, can be passed as a paramter e.g. if
            it was already loaded before to reduce number of file accesses.
        best_meta_for_sensor : bool, optional (default: True)
            Compare the sensor depth to metadata that is available in multiple
            depth layers (e.g. static metadata variables). Find the variable
            for which the depth matches best with the sensor depth.
        """
        header_elements, filename_elements = self._get_metadata_from_file()

        if len(filename_elements) == 5 and len(header_elements) == 16:
            self.file_type = 'ceop'
            raise RuntimeError('CEOP format not supported')
        elif len(header_elements) == 15 and len(filename_elements) >= 9:
            metadata, depth = self._get_metadata_ceop_sep()
            self.file_type = 'ceop_sep'
        elif len(header_elements) < 14 and len(filename_elements) >= 9:
            metadata, depth = self._get_metadata_header_values()
            self.file_type = 'header_values'
        else:
            raise IOError(f"Unknown file format found for: {self.file_path}")
            #logger.warning(f"Unknown file type: {self.file_path} in {self.archive}")

        # metadata.add('depth', None, depth)

        if static_meta is None:
            static_meta = self._get_static_metadata_from_csv()

        metadata = metadata.merge(static_meta)

        if best_meta_for_sensor:
            depth = metadata['sensor'].depth
            metadata = metadata.best_meta_for_depth(depth)

        self.metadata = metadata

        return self.metadata

    def _get_static_metadata_from_csv(self):
        """
        Read static metadata from csv file in the same directory as the ismn
        data file.

        Returns
        -------
        static_meta : MetaData
            Dictionary of static metadata
        """
        csv = self.root.find_files(self.file_path.parent, '*.csv')

        try:
            if len(csv) == 0:
                raise IsmnFileError("Expected 1 csv file for station, found 0. "
                                     "Use empty static metadata.")
            else:
                if len(csv) > 1:
                    warnings.warn(f"Expected 1 csv file for station, found {len(csv)}. "
                                  f"Use first file in dir.")
                static_meta_file = StaticMetaFile(self.root, csv[0])
                static_meta = static_meta_file.read_metadata()
        except IsmnFileError:
            static_meta = MetaData.from_dict(tables.CSV_META_TEMPLATE)

        return static_meta

    def _get_metadata_ceop_sep(self):
        """
        Get metadata in the file format called CEOP in separate files.

        Returns
        -------
        metadata : dict
            Metadata information.
        """
        header_elements, filename_elements = self._get_metadata_from_file()

        if len(filename_elements) > 9:
            sensor = '_'.join(filename_elements[6:len(filename_elements) - 2])
        else:
            sensor = filename_elements[6]

        if filename_elements[3] in tables.VARIABLE_LUT:
            variable = tables.VARIABLE_LUT[filename_elements[3]]
        else:
            variable = filename_elements[3]

        depth = Depth(float(filename_elements[4]),
                      float(filename_elements[5]))

        metadata = MetaData([MetaVar('network', filename_elements[1]),
                             MetaVar('station', filename_elements[2]),
                             MetaVar('variable', variable, depth),
                             MetaVar('sensor', sensor, depth),
                             MetaVar('latitude', float(header_elements[7])),
                             MetaVar('longitude', float(header_elements[8])),
                             MetaVar('elevation', float(header_elements[9])),
                             ])

        return metadata, depth


    def _get_metadata_header_values(self):
        """
        Get metadata file in the format called Header Values.

        Returns
        -------
        metadata : dict
            Metadata information.
        """
        header_elements, filename_elements = self._get_metadata_from_file()

        if len(filename_elements) > 9:
            sensor = '_'.join(filename_elements[6:len(filename_elements) - 2])
        else:
            sensor = filename_elements[6]

        if filename_elements[3] in tables.VARIABLE_LUT:
            variable = tables.VARIABLE_LUT[filename_elements[3]]
        else:
            variable = filename_elements[3]

        depth = Depth(float(filename_elements[6]),
                      float(filename_elements[7]))

        metadata = MetaData([MetaVar('network', header_elements[1]),
                             MetaVar('station', header_elements[2]),
                             MetaVar('variable', variable, depth),
                             MetaVar('sensor', sensor, depth),
                             MetaVar('latitude', float(header_elements[3])),
                             MetaVar('longitude', float(header_elements[4])),
                             MetaVar('elevation', float(header_elements[5])),
                             ])

        return metadata, depth

    def _get_metadata_from_file(self, delim='_'):
        """
        Read first line of file and split filename.
        Information is used to collect metadata information for all
        ISMN formats.

        Parameters
        ----------
        delim : str, optional
            File basename delimiter.

        Returns
        -------
        header_elements : list[str]
            First line of file split into list
        file_basename_elements : list[str]
            File basename without path split by 'delim'
        """
        if self.root.zip:
            if not self.root.isopen: self.root.open()
            with TemporaryDirectory(prefix='ismn', dir=self.temp_root) as tempdir:
                filename = self.root.extract_file(self.file_path, tempdir)

                with filename.open(mode='r', newline=None) as f:
                    header = f.readline()
        else:
            filename = self.root.path / self.file_path

            with filename.open(mode='r', newline=None) as f:
                header = f.readline()

        header_elements = header.split()
        path, basename = os.path.split(filename)
        file_basename_elements = basename.split(delim)

        return header_elements, file_basename_elements

    def _read_format_ceop_sep(self):
        """
        Read data in the file format called CEOP in separate files.
        """
        var = self.metadata['variable']
        varname = var.name
        names = ['date', 'time', varname, varname + '_flag', varname + '_orig_flag']
        usecols = [0, 1, 12, 13, 14]

        self.data = self._read_csv(names, usecols)

    def _read_format_header_values(self):
        """
        Read data file in the format called Header Values.
        """
        var = self.metadata['variable']
        varname = var.name
        names = ['date', 'time', varname, varname + '_flag', varname + '_orig_flag']

        self.data = self._read_csv(names, skiprows=1)

    def _read_csv(self, names=None, usecols=None, skiprows=0):
        """
        Read data.

        Parameters
        ----------
        names : list, optional
            List of column names to use.
        usecols : list, optional
            Return a subset of the columns.

        Returns
        -------
        data : pandas.DataFrame
            Time series.
        """
        readf = lambda f: pd.read_csv(f, skiprows=skiprows, usecols=usecols,
                                      names=names, delim_whitespace=True,
                                      parse_dates=[[0, 1]])
        if self.root.zip:
            with TemporaryDirectory(prefix='ismn', dir=self.temp_root) as tempdir:
                filename = self.root.extract_file(self.file_path, tempdir)
                data = readf(filename)
        else:
            data = readf(self.root.path / self.file_path)

        data.set_index('date_time', inplace=True)

        return data

if __name__ == '__main__':
    import pickle
    archive =  r"C:\Temp\delete_me\ismn\testdata_ceop.zip"
    filepath = "FMI/SAA111/FMI_FMI_SAA111_sm_0.050000_0.050000_5TE_20101001_20201005.stm"
    nodat = DataFile(archive, filepath, load_data=False)
    nodat.read_metadata()

    meta = nodat.metadata.get_meta_for_depth(0,0.05)

    flag = nodat.check_metadata('soil_moisture', 0, 0.1, {'lc_2010':110})
    nodat.close()
    with open(r"C:\Temp\delete_me\dumpnodat.pkl", 'wb') as handle:
        pickle.dump(nodat, handle)

    with open(r"C:\Temp\delete_me\dumpnodat.pkl", 'rb') as handle:
        nodat = pickle.load(handle)

    nodat.open()

    data = nodat.load_data()

    dat = DataFile(archive, filepath, load_data=True)
    dat.close()
    with open(r"C:\Temp\delete_me\dumpdat.pkl", 'wb') as handle:
        pickle.dump(dat, handle)

    with open(r"C:\Temp\delete_me\dumpdat.pkl", 'rb') as handle:
        dat = pickle.load(handle)

    dat.data
    dat.load_data()
    print(nodat.data)
    print(dat.data)