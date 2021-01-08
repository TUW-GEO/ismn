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
import logging
import pandas as pd
from tempfile import gettempdir
from pathlib import Path, PurePosixPath
from tqdm import tqdm
from typing import Union
from multiprocessing import Pool, cpu_count
from operator import itemgetter

from ismn.base import IsmnRoot
from ismn.const import *
from ismn.filehandlers import DataFile, StaticMetaFile
from ismn.meta import MetaData, MetaVar, Depth


def _read_station_dir(
        root: Union[IsmnRoot, Path, str],  # Path for reading from zip, avoid serialisation error
        stat_dir: Union[Path, str],
        temp_root: Path) -> (dict, list):
    """
    Parallelizable function to read metadata for files in station dir
    """
    infos = []

    if not isinstance(root, IsmnRoot):
        proc_root = True
        root = IsmnRoot(root)
    else:
        proc_root = False

    csv = root.find_files(stat_dir, '*.csv')

    try:
        if len(csv) == 0:
            raise IsmnFileError("Expected 1 csv file for station, found 0. "
                                "Use empty static metadata.")
        else:
            if len(csv) > 1:
                infos.append(f"Expected 1 csv file for station, found {len(csv)}. "
                             f"Use first file in dir.")
            static_meta_file = StaticMetaFile(root, csv[0], load_metadata=True,
                                              temp_root=temp_root)
            station_meta = static_meta_file.metadata
    except IsmnFileError:
        station_meta = MetaData.from_dict(CSV_META_TEMPLATE)

    data_files = root.find_files(stat_dir, '*.stm')

    filelist = []

    for file_path in data_files:
        try:
            f = DataFile(root, file_path, temp_root=temp_root)
        except IOError as e:
            infos.append(f'Error loading ismn file: {e}')
            continue

        f.metadata.merge(station_meta, inplace=True)

        f.metadata = f.metadata.best_meta_for_depth(
            Depth(f.metadata['instrument'].depth.start,
                  f.metadata['instrument'].depth.end))

        network = f.metadata['network'].val
        station = f.metadata['station'].val

        filelist.append((network, station, f))

        infos.append(f"Processed file {file_path}")

    if proc_root:
        root.close()

    return filelist, infos


class IsmnFileCollection(object):
    """
    The IsmnFileCollection class contains a list of file handlers to access data
    in the given data directory. The file list can be loaded from a previously
    stored csv file, or built by iterating over all files in the data root.
    This class also contains function to load filehandlers for certain networks
    only.

    Attributes
    ----------
    root : IsmnRoot
        Root object where data is stored.
    filelist : OrderedDict
        A collection of filehandlers and network names
    temp_root : Path
        Temporary root dir.
    """

    def __init__(self,
                 root,
                 filelist,
                 temp_root=gettempdir()):
        """
        Parameters
        ----------
        root : IsmnRoot
            Root object where data is stored.
        filelist : OrderedDict
            A collection of filehandler stored in lists with network name as key.
        temp_root : Path or str, optional (default : gettempdir())
            Root directory where a separate subdir for temporary files
            will be created (and deleted).
        """
        self.root = root
        self.filelist = filelist
        self.temp_root = Path(temp_root)

        os.makedirs(self.temp_root, exist_ok=True)

    def __repr__(self):
        return f"{self.__class__.__name__} for {len(self.filelist.keys())} Networks"

    @classmethod
    def build_from_scratch(cls, data_root, parallel=True, log_path=None,
                           temp_root=gettempdir()):
        """
        Parameters
        ----------
        data_root : IsmnRoot or str or Path
            Root path of ISMN files or path to metadata pkl file.
            i.e. path to the downloaded zip file or the extracted zip directory (faster)
            or a file list that contains these infos already.
        parallel : bool, optional (default: True)
            Speed up metadata collecting with multiple processes.
        log_path : str or Path, optional (default: None)
            Path where the log file is created. If None is set, no log file
            will be written.
        temp_root : str or Path, (default: gettempdir())
            Temporary folder where extracted data is copied during reading from
            zip archive.
        """
        if isinstance(data_root, IsmnRoot):
            root = data_root
        else:
            root = IsmnRoot(data_root)

        os.makedirs(temp_root, exist_ok=True)

        if log_path is not None:
            log_file = os.path.join(log_path, f"{root.name}.log")
        else:
            log_file = None

        if log_file:
            os.makedirs(os.path.dirname(log_file), exist_ok=True)
            logging.basicConfig(filename=log_file, level=logging.INFO,
                                format='%(levelname)s %(asctime)s %(message)s',
                                datefmt='%Y-%m-%d %H:%M:%S')

        n_proc = 1 if not parallel else cpu_count()

        logging.info(f"Collecting metadata with {n_proc} processes.")

        print(f"Processing metadata for all ismn stations into folder {root.path}."
              f" This may take a few minutes, but is only done once ...")

        process_stat_dirs = []
        for net_dir, stat_dirs in root.cont.items():
            process_stat_dirs += list(stat_dirs)

        args = [(root.path if root.zip else root, d, temp_root)
                for d in process_stat_dirs]

        pbar = tqdm(total=len(args), desc='Files Processed:')

        fl_elements = []

        def update(r):
            net_stat_fh, infos = r
            for i in infos:
                logging.info(i)
            for elements in net_stat_fh:
                fl_elements.append(elements)
            pbar.update()

        with Pool(n_proc) as pool:
            for arg in args:
                pool.apply_async(_read_station_dir, arg, callback=update,
                                 error_callback=logging.error)

            pool.close()
            pool.join()

        fl_elements.sort(key=itemgetter(0, 1))  # sort by net name, stat name

        # to ensure alphabetical order... not sure if necessary, slow?
        filelist = OrderedDict([])
        for net, stat, fh in fl_elements:
            if net not in filelist.keys():
                filelist[net] = []
            filelist[net].append(fh)

        logging.info(f"All processes finished.")

        return cls(root, filelist=filelist)

    @classmethod
    def from_metadata_csv(cls, data_root, meta_csv_file, temp_root=gettempdir()):
        """
        Load a previously created and stored filelist from pkl.

        Parameters
        ----------
        data_root : IsmnRoot or str or Path
            Path where the ismn data is stored, can also be a zip file
        meta_csv_file : str or Path
            Csv file where the metadata is stored.
        temp_root : str or Path, optional (default: gettempdir())
            Temporary folder where extracted data is copied during reading from
            zip archive.
        """
        if isinstance(data_root, IsmnRoot):
            root = data_root
        else:
            root = IsmnRoot(data_root)

        print(f"Found existing ismn metadata in {meta_csv_file}.")

        metadata_df = pd.read_csv(meta_csv_file, index_col=0, header=[0, 1],
                                  low_memory=False)

        # parse date cols as datetime
        for col in ['timerange_from', 'timerange_to']:
            metadata_df[col, 'val'] = pd.to_datetime(metadata_df[col, 'val'])

        lvars = []
        for c in metadata_df.columns:
            if c[0] not in lvars:
                lvars.append(c[0])

        # we assume triples for all vars except these, so they must be at the end
        assert lvars[-2:] == ['file_path', 'file_type'], \
            "file_type and file_path must be at the end."
        lvars = lvars[:-2]

        filelist = OrderedDict([])

        for row in metadata_df.values:  # todo: slow!?? parallelise?
            metavars = []

            for j, metavar_name in enumerate(lvars):
                depth_from, depth_to, val = row[j * 3], row[j * 3 + 1], row[j * 3 + 2]

                if np.all(np.isnan(np.array([depth_from, depth_to]))):
                    depth = None
                else:
                    depth = Depth(depth_from, depth_to)

                metavar = MetaVar(metavar_name, val, depth)
                metavars.append(metavar)

            metadata = MetaData(metavars)
            f = DataFile(root=root,
                         file_path=str(PurePosixPath(row[-2])),
                         load_metadata=False,
                         temp_root=temp_root)

            f.metadata = metadata
            f.file_type = row[-1]

            network = f.metadata['network'].val

            if network not in filelist.keys():
                filelist[network] = []

            filelist[network].append(f)

        return cls(root, filelist=filelist)

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

        for i, filehandler in enumerate(self.iter_filehandlers()):
            df = filehandler.metadata.to_pd(True, dropna=False)
            df['file_path'] = str(PurePosixPath(filehandler.file_path))
            df['file_type'] = filehandler.file_type

            df.index = [i]
            dfs.append(df)
            i += 1

        dfs = pd.concat(dfs, axis=0, sort=True)
        cols_end = ['file_path', 'file_type']

        dfs = dfs[[c for c in dfs.columns if c[0] not in cols_end] +
                  [c for c in dfs.columns if c[0] in cols_end]]
        dfs = dfs.fillna(np.nan)

        os.makedirs(Path(os.path.dirname(meta_csv_file)), exist_ok=True)
        dfs.to_csv(meta_csv_file)

    def get_filehandler(self, idx):
        """
        Get the nth filehandler in a list of all filehandlers for all networks.
        e.g. if there are 2 networks, with 3 filehandlers/sensors each, idx=4
        will return the first filehandler of the second network.

        Parameters
        ----------
        idx int
            Index of filehandler to read.

        Returns
        -------
        filehandler : DataFile
            nth filehandler of all filehandlers in the sorted list.
        """
        fs = 0
        for net, files in self.filelist.items():
            l = len(files)
            if fs + l > idx:
                return files[idx - fs]
            else:
                fs += l

    def iter_filehandlers(self, networks=None):
        """
        Iterator over files for networks

        Parameters
        ----------
        networks : list, optional (default: None)
            Name of networks to get files for, or None to use all networks.

        Yields
        -------
        file : DataFile
            Filehandler with metadata
        """
        for net, files in self.filelist.items():
            if (networks is None) or (net in networks):
                for f in files:
                    yield f

    def close(self):
        # close root and all filehandlers
        self.root.close()
        for f in self.iter_filehandlers():
            f.close()