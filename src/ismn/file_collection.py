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

from ismn.base import IsmnRoot
from ismn.components import *
from ismn.tables import *
from ismn.filehandlers import DataFile
from ismn.meta import MetaData, MetaVar

import ismn

pkg_version = ismn.__version__

import pandas as pd
from tempfile import gettempdir
from pathlib import Path, PurePosixPath

from typing import Union
from multiprocessing import Pool, cpu_count

class MetadataError(IOError):
    pass


class ISMNError(Exception):
    pass


def build_filelist_from_csv(data_path, csv_path) -> pd.DataFrame:
    """
    Use a pre-build metadata csv file to load the filelist

    Parameters
    ----------
    data_path : str or Path
        Root path where the data is store. Paths in csv file are RELATIVE to
        this path.
    csv_path : str or Path
        Path to the metadata csv file.
    """
    print(f"Found existing ismn metadata in {csv_path}.")

    metadata_df = pd.read_csv(csv_path, index_col=0, header=[0,1])

    file_paths = metadata_df['file_path'].iloc[:, 0].values
    metadata_df = metadata_df.drop(columns=['file_path'])

    filelist = {'network': [], 'station': [], 'instrument': [],
                'variable': [], 'sensor_depth_from': [], 'sensor_depth_to': [],
                'timerange_from': [], 'timerange_to': [],
                'file_path': [],  'filehandler': [], # file info, to reload file
               }

    for i, row in metadata_df.iterrows():
        metavars = []
        level = list(row.index.names).index('name')
        for metavar_name in row.index.levels[level]:
            if metavar_name == 'file_path': continue
            var = row.loc[metavar_name]
            depth_from, depth_to = var['depth_from'], var['depth_to']

            if np.all(np.isnan(np.array([depth_from, depth_to]))):
                depth = None
            else:
                depth = Depth(depth_from, depth_to)

            metavar = MetaVar(metavar_name, var['val'], depth)
            metavars.append(metavar)

        metadata = MetaData(metavars)
        f = DataFile(root=data_path,
                     file_path=str(PurePosixPath(file_paths[i])),
                     load_metadata=False)
        f.metadata = metadata

        filelist['network'].append(f.metadata['network'].val)
        filelist['station'].append(f.metadata['station'].val)
        filelist['instrument'].append(f.metadata['instrument'].val)
        filelist['variable'].append(f.metadata['variable'].val)

        filelist['sensor_depth_from'].append(f.metadata['instrument'].depth.start)
        filelist['sensor_depth_to'].append(f.metadata['instrument'].depth.end)

        filelist['timerange_from'].append(pd.Timestamp(f.metadata['timerange_from'].val))
        filelist['timerange_to'].append(pd.Timestamp(f.metadata['timerange_to'].val))
        filelist['file_path'].append(str(PurePosixPath(f.file_path)))
        filelist['filehandler'].append(f)

    files = pd.DataFrame.from_dict(filelist)
    files.index = files.index.astype('int')

    return files


def _read_station_dir(
        root: Union[IsmnRoot, Path, str],  # Path for reading from zip, avoid serialisation error
        stat_dir: Union[Path, str],
        temp_root: Path) -> dict:
    """
    Parallelizable function to read metadata for files in station dir
    """

    if not isinstance(root, IsmnRoot):
        proc_root = True
        root = IsmnRoot(root)
    else:
        proc_root = False

    filelist = {}
    for var in ['network', 'station', 'instrument', 'variable',
                'sensor_depth_from', 'sensor_depth_to', 'timerange_from',
                'timerange_to', 'file_path', 'filehandler']:
        filelist[var] = []

    data_files = root.find_files(stat_dir, '*.stm')
    static_meta = None

    for file_path in data_files:
        try:
            f = DataFile(root, file_path,
                         static_meta=static_meta,
                         temp_root=temp_root)
        except IOError as e:
            logger.error(f'Error loading ismn file: {e}')
            continue

        f.metadata = f.metadata.best_meta_for_depth(
            Depth(f.metadata['instrument'].depth.start,
                  f.metadata['instrument'].depth.end))

        filelist['network'].append(f.metadata['network'].val)
        filelist['station'].append(f.metadata['station'].val)
        filelist['instrument'].append(f.metadata['instrument'].val)
        filelist['variable'].append(f.metadata['variable'].val)
        filelist['sensor_depth_from'].append(f.metadata['instrument'].depth.start)
        filelist['sensor_depth_to'].append(f.metadata['instrument'].depth.end)
        filelist['timerange_from'].append(f.metadata['timerange_from'].val)
        filelist['timerange_to'].append(f.metadata['timerange_to'].val)
        filelist['file_path'].append(str(PurePosixPath(f.file_path)))
        filelist['filehandler'].append(f)

    if proc_root:
        root.close()

    return filelist


def build_filelist_from_data(root, parallel=True, temp_root=gettempdir()):
    """
    Build the file list
    Iterate over all networks and station folders in the root directory
    or archive. For each ismn station the according static metadata is
    loaded and stored in the file handler together with the specific meta
    data for each file.
    In the file list, for faster filtering, some metadata information is
    stored in separate columns, i.e. network/station names, variables and
    depths.

    Parameters
    ----------
    data_root : IsmnRoot
        Root where the data is stored.
    parallel : bool, optional (default: True)
        Run metadata extraction in parallel (use all CPU cores).
    temp_root : str or Path
        Root path where temporary files are stored.
    """
    print(f"Processing metadata for all ismn stations in {root.path}."
          f" This may take a few minutes ...")

    process_stat_dirs = []
    for net_dir, stat_dirs in root.cont.items():
        process_stat_dirs += list(stat_dirs)

    root = root.path if root.zip else root
    args = [(root, d, temp_root) for d in process_stat_dirs]

    pool = Pool(1 if not parallel else cpu_count())
    res = pool.starmap_async(_read_station_dir, args)

    pool.close()
    pool.join()

    result = res.get()

    df = [pd.DataFrame.from_dict(r) for r in result]
    df = pd.concat(df, axis=0, ignore_index=True, sort=False)
    df = df.sort_values(['network', 'station']).reset_index(drop=True)

    return df

class IsmnFileCollection(object):
    """
    The IsmnFileCollection class contains a pandas data frame with all ismn files
    in the given data directory. The file list can be loaded from a previously
    stored item, or built by iterating over all files in the data root.
    This class also contains function to filter the file list for faster access
    to files for a certain network/variable etc.
    """

    def __init__(self, root, filelist):
        """
        Use the @classmethods to create this object!

        Parameters
        ----------
        root : IsmnRoot
            Root object where data is stored.
        filelist : pd.DataFarme, optional (Default: None)
            A pre-built filelist to use, if None is give one from data_path is created.
        """
        self.root = root
        self.files = filelist

    @classmethod
    def from_scratch(cls, data_root, parallel=True, temp_root=gettempdir()):
        """
        Parameters
        ----------
        data_root : IsmnRoot or str or Path
            Root path of ISMN files or path to metadata pkl file.
            i.e. path to the downloaded zip file or the extracted zip directory (faster)
            or a file list that contains these infos already.
        parallel : bool, optional (default: True)
            Speed up metadata collecting with multiple processes.
        temp_root : str or Path
            Temporary folder where extraced data is copied during reading from
            zip archive.
        """
        if isinstance(data_root, IsmnRoot):
            root = data_root
        else:
            root = IsmnRoot(data_root)

        if not os.path.exists(temp_root):
            os.makedirs(temp_root, exist_ok=True)

        filelist = build_filelist_from_data(
            root, parallel=parallel, temp_root=temp_root)

        return cls(root, filelist=filelist)

    @classmethod
    def from_metadata_csv(cls, data_root, meta_csv_file):
        """
        Load a previously created and stored filelist from pkl.

        Parameters
        ----------
        data_root : IsmnRoot or str or Path
            Path where the ismn data is stored, can also be a zip file
        meta_csv_file : str or Path
            Csv file where the metadata is stored.
        """
        if isinstance(data_root, IsmnRoot):
            data_root = data_root.path
        else:
            data_root = Path(data_root)

        filelist = build_filelist_from_csv(data_root, csv_path=meta_csv_file)

        return cls(IsmnRoot(data_root), filelist=filelist)

    def to_metadata_csv(self, meta_csv_file):
        """
        Write filehandle metadata from filelist to metdata csv that contains
        ALL metadata / variables of the filehander. Can be read back in as
        filelist with filehandlers using from_metadata_csv().

        Parameters
        ----------
        meta_csv_file : Path or str, optional (default: None)
            Directory where the csv file with the correct name is crated
        """

        dfs = []
        filehandlers = self.files['filehandler'].values
        file_paths = self.files['file_path'].values

        for i, (filehandler, file_path) in enumerate(zip(filehandlers, file_paths)):
            df = filehandler.metadata.to_pd(True)
            df['file_path'] = file_path
            df.index = [i]
            dfs.append(df)

        metadata_df = pd.concat(dfs, axis=0)
        metadata_df = metadata_df.fillna(np.nan)

        os.makedirs(Path(os.path.dirname(meta_csv_file)), exist_ok=True)
        metadata_df.to_csv(meta_csv_file)

    def filter_col_val(self, col, vals, return_index=False):
        """
        Filter the file list for certain values in the columns, except the
        filehandler column and the depth range.

        Parameters
        ----------
        col : str
            Column based on which the filtering is performed.
        vals : Any
            Value(s) that are allowed in col.
        return_index : bool, optional (default: False)
            Return only the index, no the filtered data frame.

        Returns
        -------
        filtered_filelist or filtered_index : pd.DataFrame or np.array
            Filtered file list or indices of included elements
        """

        if col not in self.files.columns:
            raise ValueError(f"Column {col} is not in file list.")

        if col in ['filehandler']:
            raise ValueError(f"Cannot filter based on column {col}.")

        mask = np.isin(self.files[col], np.atleast_1d(vals))
        idx = self.files.loc[mask].index.values

        if return_index:
            return idx
        else:
            return self.files.loc[idx, :]

    def filter_depth(self, min_depth=-np.inf, max_depth=np.inf, return_index=False,
                     only_consider_depth_from=False):
        """
        Filter filelist by depth_from and depth_to columns.

        Parameters
        ----------
        min_depth : float, optional (default: -np.inf)
            Return files below this depth, i.e. with depth_from greater than this
            are kept.
        max_depth : float, optional (default: np.inf)
            Return files above this depth, i.e. with depth_to smaller than this
            are kept. If only_consider_depth_from is selected, files with
            depth_to smaller than this are kept.
        return_index : bool, optional (default: False)
            Return only the index, no the filtered data frame.
        only_consider_depth_from : bool, optional (default: False)
            Only check whether the depth_from of the file is within the passed
            depth range.

        Returns
        -------
        filtered_filelist or filtered_index : pd.DataFrame or np.array
            Filtered file list or indices of included elements
        """
        mask_from = np.greater_equal(self.files['sensor_depth_from'], min_depth)

        if only_consider_depth_from:
            mask_to = np.less_equal(self.files['sensor_depth_from'], max_depth)
        else:
            mask_to = np.less_equal(self.files['sensor_depth_to'], max_depth)

        idx = self.files.loc[(mask_from & mask_to)].index.values

        if return_index:
            return idx
        else:
            return self.files.loc[idx, :]

    def filter_metadata(self, filter_dict: dict, filelist=None, return_index=False):
        """
        Filter file list by comparing file metadata to passed metadata dict.

        Parameters
        ----------
        filter_dict: dict
            Additional metadata keys and values for which the file list is filtered
            e.g. {'lc_2010': 10} to filter for a landcover class.
        filelist : pd.DataFrame, optional (default: None)
            The filelist to filter, if None is passed, self.files is used as
            for the other filter functions.
        return_index : bool, optional (default: False)
            Return only the index, no the filtered data frame.

        Returns
        -------
        filtered_filelist or filtered_index : pd.DataFrame or np.array
            Filtered file list or indices of included elements
        """
        # FIler based on the metadata in the filehander

        if filelist is None:
            filelist = self.files

        filehandlers = filelist['filehandler'].values

        mask = []
        for filehandler in filehandlers:
            flags = []
            for meta_key, meta_vals in filter_dict.items():
                meta_vals = np.atleast_1d(meta_vals)

                if meta_key not in filehandler.metadata.keys():
                    raise ValueError(f"{meta_key} is not a valid metadata variable")

                # check if the metadata val is one of the passed allowed values
                flag = filehandler.metadata[meta_key].val in meta_vals
                flags.append(flag)
            mask.append(all(flags))

        idx = filelist.loc[mask].index.values

        if return_index:
            return idx
        else:
            return filelist.loc[idx, :]

    def close(self):
        # close root and all filehandlers
        self.root.close()
        for f in self.files['filehandler'].values:
            f.close()

if __name__ == '__main__':
    root_path = r"C:\Temp\delete_me\ismn\testdata_ceop.zip"
    fc = IsmnFileCollection.from_scratch(root_path)
    fc.to_metadata_csv(r"C:\Temp\delete_me\ismn")

    fc = IsmnFileCollection.from_metadata_csv(root_path,
                                              r"C:\Temp\delete_me\ismn")




