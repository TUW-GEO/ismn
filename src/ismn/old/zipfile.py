# -*- coding: utf-8 -*-
from ismn.readers import IsmnFile
from tempfile import TemporaryDirectory, gettempdir
import zipfile as zf
import numpy as np
from ismn.utils import scan_archive


import os

def extract_from_archive(zip, subdir, out_path):
    """
    Extract all files in a zip file under the passed subdir

    Parameters
    ----------
    zip : str
        Path to ismn zipfile
    subdir : str
        No leading /!
        Subdir to extract, or a single file in the zip archive to extract
    out_path : str
        Path where the extracted file(s) is/are stored
    """
    with zf.ZipFile(zip) as zi:
        filelist = np.array(zi.namelist())
        if subdir in filelist: # single file was passed
            zi.extract(member=subdir, path=out_path)
        else: # subdir was passed
            filterlist = filter(lambda x: x.startswith(subdir), filelist)
            zi.extractall(members=list(filterlist), path=out_path)


class IsmnFile(object):

    """
    IsmnFile class represents a single ISMN file. This represents only DATA files
    not metadata csv files.

    Parameters
    ----------
    filename : str
        Filename.

    load_data : bool, optional
        If True data will be loaded during metadata reading.

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
    load_data()
        Load data from file.
    read_data()
        Read data in file.
    _read_metadata()
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

    def __init__(self, filename, load_data=False, static_meta=None):

        if not os.path.isfile(filename):
            raise IOError('File does not exist: {}'.format(filename))

        self.filename = filename
        self.file_type = 'undefined'
        self.metadata = {}
        self.data = None

        self.static_meta = static_meta

        self._read_metadata()

        if load_data:
            self.load_data()

    def __getitem__(self, item):
        return self.metadata[item]

    def load_data(self):
        """
        Load data from file.
        """
        if self.data is None:
            if self.file_type == 'ceop':
                # self._read_format_ceop()
                raise NotImplementedError
            elif self.file_type == 'ceop_sep':
                self._read_format_ceop_sep()
            elif self.file_type == 'header_values':
                self._read_format_header_values()
            else:
                logger.warning("Unknown file type: {}".format(self.filename))

    def read_data(self):
        """
        Read data in file.

        Returns
        -------
        data : pandas.DataFrame
            File content.
        """
        self.load_data()
        return self.data

    def _read_metadata(self):
        """
        Read metadata from file name and first line of file.
        """
        header_elements, filename_elements = self._get_metadata_from_file()

        if len(filename_elements) == 5 and len(header_elements) == 16:
            self.file_type = 'ceop'
            raise RuntimeError('CEOP format not supported')
        elif len(header_elements) == 15 and len(filename_elements) >= 9:
            self.metadata = self._get_metadata_ceop_sep()
            self.file_type = 'ceop_sep'
        elif len(header_elements) < 14 and len(filename_elements) >= 9:
            self.metadata = self._get_metadata_header_values()
            self.file_type = 'header_values'
        else:
            logger.warning("Unknown file type: {}".format(self.filename))

        if self.static_meta is None:
            self.static_meta = self._get_static_metadata_from_csv()

        self.metadata.update(self.static_meta)


    def _get_static_metadata_from_csv(self):
        """
        Read static metadata from csv file in the same directory as the ismn
        data file.

        Returns
        -------
        static_meta : OrderedDict
            Dictionary of static metadata
        """

        csv_path = glob.glob(os.path.join(os.path.dirname(self.filename), '*.csv'))
        assert len(csv_path) == 1, "More than 1 csv file found in path"
        return StaticMetaFile(csv_path[0]).read_csv()

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

        if filename_elements[3] in variable_lookup:
            variable = variable_lookup[filename_elements[3]]
        else:
            variable = filename_elements[3]

        metadata = {'network': filename_elements[1],
                    'station': filename_elements[2],
                    'variable': variable,
                    'depth': Depth(float(filename_elements[4]),
                                   float(filename_elements[5])),
                    'sensor': sensor,
                    'latitude': float(header_elements[7]),
                    'longitude': float(header_elements[8]),
                    'elevation': float(header_elements[9])}

        return metadata

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

        if filename_elements[3] in variable_lookup:
            variable = variable_lookup[filename_elements[3]]
        else:
            variable = filename_elements[3]

        metadata = {'network': header_elements[1],
                    'station': header_elements[2],
                    'latitude': float(header_elements[3]),
                    'longitude': float(header_elements[4]),
                    'elevation': float(header_elements[5]),
                    'depth': Depth(float(header_elements[6]),
                                   float(header_elements[7])),
                    'variable': variable,
                    'sensor': sensor}

        return metadata


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
        header_elements : list
            First line of file split into list
        file_basename_elements : list
            File basename without path split by 'delim'
        """
        with io.open(self.filename, mode='r', newline=None) as f:
            header = f.readline()

        header_elements = header.split()
        path, basename = os.path.split(self.filename)
        file_basename_elements = basename.split(delim)

        return header_elements, file_basename_elements

    def _read_format_ceop_sep(self):
        """
        Read data in the file format called CEOP in separate files.
        """
        names = ['date', 'time', self.metadata['variable'],
                 self.metadata['variable'] + '_flag',
                 self.metadata['variable'] + '_orig_flag']
        usecols = [0, 1, 12, 13, 14]

        self.data = self._read_csv(names, usecols)

    def _read_format_header_values(self):
        """
        Read data file in the format called Header Values.
        """
        names = ['date', 'time', self.metadata['variable'],
                 self.metadata['variable'] + '_flag',
                 self.metadata['variable'] + '_orig_flag']

        self.data = self._read_csv(names, skiprows=1)

    @zipped
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
        data = pd.read_csv(self.filename, skiprows=skiprows, usecols=usecols,
                           names=names, delim_whitespace=True,
                           parse_dates=[[0, 1]])

        data.set_index('date_time', inplace=True)

        return data


class ISMNZipFile(object):

    def __init__(self, data_path, load_data=False,
                 temp_root=gettempdir(), static_meta=None):
        """
        Represents an ISMN data file that is stored inside a downloaded zip
        archive.

        Parameters
        ----------
        data_path : str, optional (default: None)
            If the root path to the downloaded data is given, filename is relative
            to this path. Eg. when reading from zipfile, archive it the zip, while
            filename is the path inside the zip. If archive is not set, filename
            is an absolute path.
        filename
        load_data
        temp_root
        """
        self.data_path = data_path
        self.files = {}

        assert os.path.exists(data_path), f"{self.data_path} not found"
        assert zf.is_zipfile(data_path), f"{self.data_path} is not a zip archive."


        super(ISMNZipFile, self).__init__(path_unziped,  load_data, static_meta)

    def _build_filelist(self, load_data):
        """
        Scan ismn archive and build file object list.
        Reuse static metadata if possible
        """

        archive = scan_archive(self.data_path)

        for net_dir, stat_dirs in archive.items():
            for stat_dir in stat_dirs:
                root = os.path.join(self.data_path, stat_dir)
                static_meta = None

                for filename in glob.glob(os.path.join(root, '*.stm')):
                    f = IsmnFile(filename, load_data, static_meta)
                    static_meta = f.static_meta

                    if f['network'] not in self.files.keys():
                        self.files[f['network']] = {}
                    if f['station'] not in self.files[f['network']]:
                        self.files[f['network']][f['station']] = []

                    self.files[f['network']][f['station']].append(f)

    def _unzip(self):
        with TemporaryDirectory(prefix='ismn', dir=self.temp_root) as tempdir:
            extract_from_archive(self.data_path, station_dir, tempdir)
            station_path = os.path.normpath(os.path.join(tempdir, station_dir))
            station_meta = self.get_station_meta(station_path)
            metadata_catalog += station_meta
