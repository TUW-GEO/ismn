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
from ismn.base import IsmnRoot
from ismn.components import *
from ismn import const
from ismn.const import IsmnFileError
from ismn.meta import MetaVar, MetaData

from tempfile import gettempdir, TemporaryDirectory
from pathlib import Path, PurePosixPath


class IsmnFile(object):
    """
    General base class for data and static metadata files (station csv file)
    in ismn archive.

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

    def __init__(self,
                 root,
                 file_path,
                 temp_root=gettempdir()):

        if not isinstance(root, IsmnRoot):
            root = IsmnRoot(root)

        self.root = root
        self.file_path = self.root._clean_subpath(file_path)

        if self.file_path not in self.root:
            raise IOError(f'Archive does not contain file: {self.file_path}')

        if not os.path.exists(temp_root):
            os.makedirs(temp_root, exist_ok=True)

        self.temp_root = temp_root

        self.metadata = MetaData()

    def __repr__(self):
        return f"{self.__class__.__name__}({self.root.path / self.file_path})"

    def __getitem__(self, item:int):
        return [self.root, self.file_path][item]


    def check_metadata(self, variable=None, allowed_depth=None,
                       filter_meta_dict=None) -> bool:
        """
        Evaluate whether the file complies with the passed metadata requirements

        Parameters
        ----------
        variable : str, optional (default: None)
            Name of the required variable measured, e.g. soil_moisture
        allowed_depth : Depth, optional (default: None)
            Depth range that is allowed, depth in metadata must be within this
            range.
        filter_meta_dict: dict, optional (default: None)
            Additional metadata keys and values for which the file list is filtered
            e.g. {'station': 'stationname'} to filter for a station name.

        Returns
        -------
        valid : bool
            Whether the metadata complies with the passed conditions or not.
        """

        lc_cl = list(const.CSV_META_TEMPLATE.keys())

        if variable:
            if not (self.metadata['variable'].val == variable):
                return False

        if allowed_depth is not None:
            try:
                sensor_depth = self.metadata['instrument'].depth
            except AttributeError:
                sensor_depth = self.metadata['variable'].depth

            if not allowed_depth.encloses(sensor_depth):
                return False

        if filter_meta_dict:
            fil_lc_cl = [True]
            for k in filter_meta_dict.keys():
                vs = self.metadata[k]
                if isinstance(vs, MetaVar):
                    vs = [vs]

                eval = False
                ref_list = np.atleast_1d(filter_meta_dict[k]).tolist()
                for v in vs:
                    eval = v.val in ref_list
                    if eval is True:
                        break
                fil_lc_cl.append(eval)

            if not all(fil_lc_cl):
                return False

        return True


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
    load_metadata : bool, optional (default: True)
        Load metadata during initialisation.
    temp_root : Path or str, optional (default : gettempdir())
        Root directory where a separate subdir for temporary files
        will be created (and deleted).
    """

    def __init__(self,
                 root,
                 file_path,
                 load_metadata=True,
                 temp_root=gettempdir()):

        super(StaticMetaFile, self).__init__(root, file_path, temp_root)

        if self.file_path.suffix.lower() != '.csv':
            raise IsmnFileError(f'CSV file expected for StaticMetaFile object')

        if load_metadata:
            self.metadata = self.read_metadata()

    @staticmethod
    def _read_field(data, fieldname: str, new_name=None) -> np.array:
        """
        Extract a field from the loaded csv metadata
        """
        field_vars = []

        if fieldname in data.index:

            froms = np.atleast_1d(data.loc[fieldname]['depth_from[m]'])
            tos = np.atleast_1d(data.loc[fieldname]['depth_to[m]'])
            vals = np.atleast_1d(data.loc[fieldname]['value'])

            for d_from, d_to, val in zip(froms, tos, vals):
                d = Depth(d_from, d_to)
                name = new_name if new_name is not None else fieldname
                try:
                    val = float(val)
                except ValueError:
                    pass  # value is actually a string, that's ok
                field_vars.append(MetaVar(name, val, d))

        return field_vars

    def _read_csv(self, csvfile: Path) -> pd.DataFrame:
        """ Load static metadata data frame from csv """
        try:
            data = pd.read_csv(csvfile, delimiter=";")
            data.set_index('quantity_name', inplace=True)
        except:
            # set columns manually
            logging.info('no header: {}'.format(csvfile))
            data = pd.read_csv(csvfile, delimiter=";", header=None)
            cols = list(data.columns)
            cols[:len(const.CSV_COLS)] = const.CSV_COLS  # todo: not safe
            data.columns = cols
            data.set_index('quantity_name', inplace=True)

        return data

    def read_metadata(self) -> MetaData:
        """
        Read csv file containing static variables into data frame.

        Returns
        -------
        metadata : MetaData
            Static metadata read from csv file.
        """
        if self.root.zip:
            if not self.root.isopen: self.root.open()
            with TemporaryDirectory(prefix='ismn', dir=self.temp_root) as tempdir:
                extracted = self.root.extract_file(self.file_path, tempdir)
                data = self._read_csv(extracted)
        else:
            data = self._read_csv(self.root.path / self.file_path)

        # read landcover classifications
        lc = data.loc[['land cover classification']][['value', 'quantity_source_name']]

        lc_dict = {'CCI_landcover_2000': const.CSV_META_TEMPLATE['lc_2000'],
                   'CCI_landcover_2005': const.CSV_META_TEMPLATE['lc_2005'],
                   'CCI_landcover_2010': const.CSV_META_TEMPLATE['lc_2010'],
                   'insitu': const.CSV_META_TEMPLATE['lc_insitu']}

        cl_dict = {'koeppen_geiger_2007': const.CSV_META_TEMPLATE['climate_KG'],
                   'insitu': const.CSV_META_TEMPLATE['climate_insitu']}

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
        cl = data.loc[['climate classification']][['value', 'quantity_source_name']]
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
            'saturation': self._read_field(data, 'saturation'),
            'clay_fraction': self._read_field(data, 'clay fraction', const.VARIABLE_LUT['cl_h']),
            'sand_fraction': self._read_field(data, 'sand fraction', const.VARIABLE_LUT['sa_h']),
            'silt_fraction': self._read_field(data, 'silt fraction', const.VARIABLE_LUT['si_h']),
            'organic_carbon': self._read_field(data, 'organic carbon', const.VARIABLE_LUT['oc_h']),
        }

        for name, vars in static_meta.items():
            if len(vars) > 0:
                metavars += vars
            else:
                metavars.append(MetaVar(name, const.CSV_META_TEMPLATE[name]))

        metadata = MetaData(metavars)

        return metadata


class DataFile(IsmnFile):
    """
    IsmnFile class represents a single ISMN data file.
    This represents only .stm data files not metadata csv files.

    Parameters
    ----------
    root : IsmnRoot or str
        Archive to the downloaded data.
    file_path : str or Path
        Path in the archive to the ismn file. No leading slash!
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

    Methods
    -------
    check_metadata(self, variable, min_depth=0, max_depth=0.1, filter_static_vars=None)
        Evaluate whether the file complies with the passed metadata requirements
    read_data()
        Read data in file.
    read_metadata()
        Read metadata from file name and first line of file.
    """

    def __init__(self,
                 root,
                 file_path,
                 load_metadata=True,
                 temp_root=gettempdir()):

        super(DataFile, self).__init__(root, file_path, temp_root)

        self.file_type = 'undefined'

        self.metadata = None

        if load_metadata:
            self.metadata = self.read_metadata(best_meta_for_sensor=True)

    def _read_lines(self, filename):
        """
        Read fist and last line from file as list, skips empty lines.
        """
        with filename.open(mode='rb', newline=None) as f:
            lines = f.read().splitlines()
            headr = lines[0].split()

            last, scnd = [], []
            i = 1
            while (not last) or (not scnd):
                if not last:
                    last = lines[-i].split()
                if not scnd:
                    scnd = lines[i].split()
                i += 1

        headr = [s.decode('ascii') for s in headr]
        scnd = [s.decode('ascii') for s in scnd]
        last = [s.decode('ascii') for s in last]
        
        return headr, scnd, last

    def _get_metadata_ceop_sep(self, elements=None):
        """
        Get metadata in the file format called CEOP in separate files.

        Parameters
        ----------
        elements : dict, optional (default: None)
            Previously loaded elements can be passed here to avoid reading the
            file again.
        Returns
        -------
        metadata : dict
            Metadata information.
        """
        if elements:
            headr = elements['headr']
            last = elements['last']
            fname = elements['fname']
        else:
            headr, _, last, fname = self._get_elements_from_file()

        if len(fname) > 9:
            instr = '_'.join(fname[6:len(fname) - 2])
        else:
            instr = fname[6]

        if fname[3] in const.VARIABLE_LUT:
            variable = const.VARIABLE_LUT[fname[3]]
        else:
            variable = fname[3]

        timerange_from = pd.to_datetime(' '.join(headr[:2]))
        timerange_to = pd.to_datetime(' '.join(last[:2]))

        depth = Depth(float(fname[4]),
                      float(fname[5]))

        metadata = MetaData([MetaVar('network', fname[1]),
                             MetaVar('station', fname[2]),
                             MetaVar('variable', variable, depth),
                             MetaVar('instrument', instr, depth),
                             MetaVar('timerange_from', timerange_from),
                             MetaVar('timerange_to', timerange_to),
                             MetaVar('latitude', float(headr[7])),
                             MetaVar('longitude', float(headr[8])),
                             MetaVar('elevation', float(headr[9])),
                             ])

        return metadata, depth

    def _get_metadata_header_values(self, elements=None):
        """
        Get metadata file in the format called Header Values.

        Parameters
        ----------
        elements : dict, optional (default: None)
            Previously loaded elements can be passed here to avoid reading the
            file again.

        Returns
        -------
        metadata : dict
            Metadata information.
        """
        if elements:
            headr = elements['headr']
            scnd = elements['scnd']
            last = elements['last']
            fname = elements['fname']
        else:
            headr, scnd, last, fname = self._get_elements_from_file()

        if len(fname) > 9:
            instrument = '_'.join(fname[6:len(fname) - 2])
        else:
            instrument = fname[6]

        if fname[3] in const.VARIABLE_LUT:
            variable = const.VARIABLE_LUT[fname[3]]
        else:
            variable = fname[3]

        timerange_from = pd.to_datetime(' '.join(scnd[:2]))
        timerange_to = pd.to_datetime(' '.join(last[:2]))

        depth = Depth(float(headr[6]),
                      float(headr[7]))

        metadata = MetaData([MetaVar('network', headr[1]),
                             MetaVar('station', headr[2]),
                             MetaVar('variable', variable, depth),
                             MetaVar('instrument', instrument, depth),
                             MetaVar('timerange_from', timerange_from),
                             MetaVar('timerange_to', timerange_to),
                             MetaVar('latitude', float(headr[3])),
                             MetaVar('longitude', float(headr[4])),
                             MetaVar('elevation', float(headr[5])),
                             ])

        return metadata, depth

    def _get_elements_from_file(self, delim='_', only_basename_elements=False):
        """
        Read first line of file and split filename.
        Information is used to collect metadata information for all
        ISMN formats.

        Parameters
        ----------
        delim : str, optional
            File basename delimiter.
        only_basename_elements : bool, optional (default: False)
            Parse only the filename and not the file contents.

        Returns
        -------
        headr : list[str] or None
            First line of file split into list, None if only_filename is True
        secnd : list[str] or None
            Second line of file split into list, None if only_filename is True
        last : list[str] or None
            Last non empty line elements,  None if only_filename is True
        file_basename_elements : list[str], None if only_filename is True
            File basename without path split by 'delim'
        """
        # todo: reading last line is slow...
        if only_basename_elements:
            headr = None
            secnd = None
            last = None
        else:
            if self.root.zip:
                if not self.root.isopen: self.root.open()

                with TemporaryDirectory(prefix='ismn',
                                        dir=self.temp_root) as tempdir:
                    filename = self.root.extract_file(self.file_path, tempdir)
                    headr, secnd, last = self._read_lines(filename)

            else:
                filename = self.root.path / self.file_path
                headr, secnd, last = self._read_lines(filename)

        path, basename = os.path.split(filename)
        file_basename_elements = basename.split(delim)

        return headr, secnd, last, file_basename_elements

    def _read_format_ceop_sep(self) -> pd.DataFrame:
        """
        Read data in the file format called CEOP in separate files.
        """
        var = self.metadata['variable']
        varname = var.val
        names = ['date', 'time', varname, varname + '_flag', varname + '_orig_flag']
        usecols = [0, 1, 12, 13, 14]

        return self._read_csv(names, usecols)

    def _read_format_header_values(self) -> pd.DataFrame:
        """
        Read data file in the format called Header Values.
        """
        var = self.metadata['variable']
        varname = var.val
        names = ['date', 'time', varname, varname + '_flag', varname + '_orig_flag']

        return self._read_csv(names, skiprows=1)

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

    def read_data(self) -> pd.DataFrame:
        """
        Read data in file. Load file if necessary.

        Returns
        -------
        data : pd.DataFrame
            File content.
        """

        if not self.root.isopen: self.open()

        if self.file_type == 'ceop':
            # todo: what is this format?
            # self._read_format_ceop()
            raise NotImplementedError
        elif self.file_type == 'ceop_sep':
            return self._read_format_ceop_sep()
        elif self.file_type == 'header_values':
            return self._read_format_header_values()
        else:
            raise IOError(f"Unknown file format found for: {self.file_path}")
            # logger.warning(f"Unknown file type: {self.file_path}")

    def read_metadata(self, best_meta_for_sensor=True) -> MetaData:
        """
        Read metadata from file name and first line of file.

        Parameters
        ----------
        best_meta_for_sensor : bool, optional (default: True)
            Compare the sensor depth to metadata that is available in multiple
            depth layers (e.g. static metadata variables). Find the variable
            for which the depth matches best with the sensor depth.
        """
        try:
            headr, scnd, last, fname = self._get_elements_from_file()
        except Exception as e:
            raise IOError(f"Unknown file format found for: {self.file_path}")

        elements = dict(headr=headr, scnd=scnd, last=last, fname=fname)

        if len(fname) == 5 and len(headr) == 16:
            self.file_type = 'ceop'
            raise RuntimeError('CEOP format not supported')
        elif len(headr) == 15 and len(fname) >= 9:
            metadata, depth = self._get_metadata_ceop_sep(elements)
            self.file_type = 'ceop_sep'
        elif len(headr) < 14 and len(fname) >= 9:
            metadata, depth = self._get_metadata_header_values(elements)
            self.file_type = 'header_values'
        else:
            raise IOError(f"Unknown file format found for: {self.file_path}")
            # logger.warning(f"Unknown file type: {self.file_path} in {self.archive}")

        if best_meta_for_sensor:
            depth = metadata['instrument'].depth
            metadata = metadata.best_meta_for_depth(depth)

        self.metadata = metadata

        return self.metadata


if __name__ == '__main__':
    from timeit import timeit

    fileroot = "/home/wolfgang/code/ismn/tests/test_data/Data_seperate_files_20170810_20180809"
    filepath = "COSMOS/ARM-1/COSMOS_COSMOS_ARM-1_sm_0.000000_0.190000_Cosmic-ray-Probe_20170810_20180809.stm"

    data_meta_file = DataFile(fileroot, filepath)

    fileroot =  "/home/wolfgang/code/ismn/tests/test_data/Data_seperate_files_20170810_20180809"
    filepath = "COSMOS/ARM-1/COSMOS_COSMOS_ARM-1_static_variables.csv"

    station_meta_file = StaticMetaFile(fileroot, filepath)
