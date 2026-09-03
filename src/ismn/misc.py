import os
import warnings
from datetime import datetime
from glob import glob

from repurpose.process import parallel_process


def collect_stm_cov(data_path: str, n_proc=1, progressbar=False):
    """
    Open all .stm time series files in a directory (slow) and detect the
    latest end date across all files.

    Files are filtered to those whose 4th underscore-separated field is "sm"
    (e.g. FMI_FMI_SOD023_sm_0.050000_..._19500101_20260511.stm). The end
    date is read from the first whitespace-separated token of the last line
    of each file (format "YYYY/MM/DD HH:MM ...").

    Parameters
    ----------
    data_path : str
        Path where the .stm files are stored.
    n_proc : int, optional (default: 1)
        Number of parallel workers (threads, since this is I/O-bound).
    progressbar : bool, optional (default: False)
        Show progress bar when looping through files.

    Returns
    -------
    last_end : datetime
        The most recent end date found across all .stm files.
    """
    fl = glob(
        os.path.join(data_path, '**', '*_*_*_sm_*.stm'),
        recursive=True)

    if len(fl) == 0:
        raise ValueError(f"No matching .stm files found in {data_path}")

    def _func(f: str) -> datetime:
        fname = os.path.basename(f)
        stem = fname[:-4] if fname.lower().endswith('.stm') else fname
        parts = stem.split('_')

        # Defensive: confirm "sm" position even though glob should ensure it
        if len(parts) < 4 or parts[3] != 'sm':
            raise ValueError(
                f"Filename does not match expected pattern: {fname}")

        # End date: first token of the last line. Seek from the end of the
        # file instead of reading the whole thing — these can be large.
        with open(f, 'rb') as fh:
            fh.seek(0, os.SEEK_END)
            filesize = fh.tell()
            bufsize = min(1024, filesize)
            data = b''
            while True:
                fh.seek(-bufsize, os.SEEK_END)
                data = fh.read(bufsize)
                if data.count(b'\n') >= 2 or bufsize >= filesize:
                    break
                bufsize = min(bufsize * 2, filesize)
        text = data.decode('utf-8', errors='replace').rstrip()
        last_nl = text.rfind('\n')
        last_line = text[last_nl + 1:] if last_nl != -1 else text

        if not last_line:
            raise ValueError(f"Empty last line in {f}")
        return datetime.strptime(last_line.split()[0], '%Y/%m/%d')

    ends = parallel_process(
        _func,
        ITER_KWARGS=dict(f=fl),
        show_progress_bars=progressbar,
        backend='threading',
        n_proc=n_proc)

    return max(ends)


def write_overview(data_path: str, period_to: datetime,
                   product: str = "ISMN",
                   version: str = "v202505") -> str:
    """
    Write (or overwrite) an overview.yml file inside `data_path`, describing
    the dataset's product, version, and latest end date.

    Parameters
    ----------
    data_path : str
        The data directory where overview.yml will be written.
    product : str, optional
        Product name. Default "ISMN".
    version : str, optional
        Product version. Default "v202505".

    Returns
    -------
    out_path : str
        Absolute path of the written overview.yml file.
    """
    out_path = os.path.join(os.path.abspath(data_path), "overview.yml")

    content = (
        f"product: {product}\n"
        f"version: {version}\n"
        f"period_to: {period_to.strftime('%Y-%m-%d')}\n"
    )

    with open(out_path, 'w') as f:
        f.write(content)

    return out_path


if __name__ == "__main__":
    print(collect_stm_cov("/tmp/ISMN", n_proc=4, progressbar=True))
