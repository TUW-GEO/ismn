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

from ismn.base import IsmnRoot
from ismn.tables import *
from ismn.filehandlers import DataFile, StaticMetaFile
from ismn.meta import MetaData, MetaVar, Depth

import ismn

pkg_version = ismn.__version__

import os
import logging
import pandas as pd
from tempfile import gettempdir
from pathlib import Path, PurePosixPath
from tqdm import tqdm
from typing import Union
from multiprocessing import Pool, cpu_count
from collections import defaultdict

"""
The file collection / file list / sensor list contains all information about 
the sensors to quickly rebuild a network collection.
"""

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
            static_meta_file = StaticMetaFile(root, csv[0], load_metadata=True)
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

        filelist.append((network, f))

        infos.append(f"Processed file {file_path}")

    if proc_root:
        root.close()

    return filelist, infos

class IsmnFileCollection(object):
    """
    The IsmnFileCollection class contains a pandas data frame with all ismn files
    in the given data directory. The file list can be loaded from a previously
    stored item, or built by iterating over all files in the data root.
    This class also contains function to filter the file list for faster access
    to files for a certain network/variable etc.
    """

    def __init__(self,
                 root,
                 filelist):
        """
        Use the @classmethods to create this object!

        Parameters
        ----------
        root : IsmnRoot
            Root object where data is stored.
        filelist : OrderedDict
            A collection of filehandler, grouped by network.
        """
        self.root = root
        self.files = filelist

    def __repr__(self):
        return f"{self.__class__.__name__} for {len(self.files.keys())} Networks"

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

        root = root.path if root.zip else root
        args = [(root, d, temp_root) for d in process_stat_dirs]

        pbar = tqdm(total=len(args), desc='Files Processed:')

        fl_elements = []

        def update(r):
            net_fh, infos = r
            for i in infos:
                logging.info(i)
            for f in net_fh:
                fl_elements.append(f)
            pbar.update()

        with Pool(n_proc) as pool:
            for arg in args:
                pool.apply_async(_read_station_dir, arg, callback=update,
                                 error_callback=logging.error)

            pool.close()
            pool.join()

        filelist = defaultdict(list)
        for i, j in fl_elements:
            filelist[i].append(j)

        filelist = OrderedDict(sorted(filelist.items()))

        logging.info(f"All processes finished.")

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

        print(f"Found existing ismn metadata in {meta_csv_file}.")

        metadata_df = pd.read_csv(meta_csv_file, index_col=0, header=[0, 1], low_memory=False)
        file_paths = metadata_df.loc[:, 'file_path'].values[:, 0]
        file_types = metadata_df.loc[:, 'file_type'].values[:, 0]
        metadata_df = metadata_df.drop(columns=['file_path', 'file_type'], level='name')

        metadata_df = metadata_df.reindex(sorted(metadata_df.columns), axis=1)

        level = list(metadata_df.columns.names).index('name')
        vars = list(metadata_df.columns.levels[level].values)
        if 'file_path' in vars:
            vars.pop(vars.index('file_path'))
        if 'file_type' in vars:
            vars.pop(vars.index('file_type'))

        filelist = OrderedDict([])

        for i, row in enumerate(metadata_df.values):  # todo: slow!! parallelise?
            metavars = []

            for j, metavar_name in enumerate(vars):
                if metavar_name == 'file_path': continue
                var_tup = row[j * 3:(j * 3) + 3]
                depth_from, depth_to, val = var_tup[0], var_tup[1], var_tup[2]

                if np.all(np.isnan(np.array([depth_from, depth_to]))):
                    depth = None
                else:
                    depth = Depth(depth_from, depth_to)

                metavar = MetaVar(metavar_name, val, depth)
                metavars.append(metavar)

            metadata = MetaData(metavars)
            f = DataFile(root=data_root,
                         file_path=str(PurePosixPath(file_paths[i])),
                         load_metadata=False)
            f.metadata = metadata
            f.file_type = file_types[i]

            network = f.metadata['network'].val

            if network not in filelist.keys():
                filelist[network] = []
            
            filelist[network].append(f)

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

        i = 0
        for network, filehandlers in self.files.items():
            for filehandler in filehandlers:
                df = filehandler.metadata.to_pd(True, dropna=False)
                df['file_path'] = str(PurePosixPath(filehandler.file_path))
                df['file_type'] = filehandler.file_type
                df.index = [i]
                dfs.append(df)
                i += 1

        metadata_df = pd.concat(dfs, axis=0)
        metadata_df = metadata_df.fillna(np.nan)

        os.makedirs(Path(os.path.dirname(meta_csv_file)), exist_ok=True)
        metadata_df.to_csv(meta_csv_file)

    def close(self):
        # close root and all filehandlers
        self.root.close()
        for f in self.files['filehandler'].values:
            f.close()

if __name__ == '__main__':
    root_path = r"D:\data-read\ISMN\global_20191024"
    fc = IsmnFileCollection.build_from_scratch(root_path)
    fc.to_metadata_csv(r"C:\Temp\test\test.csv")
    #files = IsmnFileCollection.from_metadata_csv(root_path, r"C:\Temp\test\test.csv")

    # fc.to_metadata_csv(r"C:\Temp\delete_me\ismn\testdata_ceop.csv")

    # fc = IsmnFileCollection.from_metadata_csv(root_path,
    #                                           r"C:\Temp\delete_me\ismn\testdata_ceop.csv")




