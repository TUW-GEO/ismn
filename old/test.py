
import os
import timeit


def read_lines(filename):
    """
    Read fist and last line from file as list, skips empty lines.
    """
    with open(filename, mode='rb', newline=None) as f:
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

    print(headr, scnd, last)

def faster(filename):
    with open(filename, "rb") as f:
        lines = f.read().splitlines()
        headr = lines[0].split()

        scnd = []
        i = 1
        while (not scnd):
            if not scnd:
                scnd = lines[i].split()
        f.seek(-2, os.SEEK_END)
        while f.read(1) != b'\n':
            f.seek(-2, os.SEEK_CUR)
        last = f.readline()
    print(headr, scnd, last)



if __name__ == '__main__':
    filename = "/home/wolfgang/data-read/ismn/Data_separate_files_20090804_20201212_5712_zm79_20201212/SCAN/RogersFarm#1/SCAN_SCAN_RogersFarm#1_ts_1.016000_1.016000_Hydraprobe-Analog-(2.5-Volt)_19500101_20201220.stm"
    #read_lines(filename)
    faster(filename)

