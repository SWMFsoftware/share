#!/usr/bin/env python
"""Prepares a geospace run by downloading inputs.

It makes the IMF input file (IMF.dat) and also replaces the F10.7 value in the
PARAM.in file. This uses an *already prepared* PARAM.in files to work with.
Please run in the same directory as your PARAM.in file (run dir).

You need swmfpy installed.
Installation:
    cd share/Python/swmfpy
    pip install -U --user .

Usage:
    python prepare_geospace.py
    # To avoid F10.7 change in PARAM.in
    python prepare_geospace.py --no_write_f107
    # To avoid writing an IMF.dat solar wind file
    python prepare_geospace.py --no_write_imf
"""
__author__ = 'Qusai Al Shidi'
__email__ = 'qusai@umich.edu'

import argparse
import datetime as dt
import numpy as np
from swmfpy import write_imf_from_omni, paramin
from swmfpy.web import get_omni_data


def arguments():
    """Make arguments"""
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-s', '--no_write_imf',
        help="Don't write IMF.dat (s)olar wind data",
        action="store_true",
        )

    parser.add_argument(
        '-l', '--no_write_f107',
        help="Don't write F10.7 f(l)ux value in the PARAM.in files",
        action="store_true",
        )

    parser.add_argument(
        '-f', '--paramin',
        help="PARAM.in file to read and write F10.7 value.",
        type=str,
        default="PARAM.in"
        )

    return parser.parse_args()


def get_times(paramin_file='PARAM.in'):
    """Start and end time for data retrieval.

       Returns:
            Time data starts from param.in file.
    """

    time_command = paramin.read_command('#STARTTIME',
                                        paramin_file)
    start_time = dt.datetime(year=int(time_command[1]),
                             month=int(time_command[2]),
                             day=int(time_command[3]),
                             hour=int(time_command[4]),
                             minute=int(time_command[5]))

    time_command = paramin.read_command('#ENDTIME',
                                        paramin_file)
    end_time = dt.datetime(year=int(time_command[1]),
                           month=int(time_command[2]),
                           day=int(time_command[3]),
                           hour=int(time_command[4]),
                           minute=int(time_command[5]))

    return (start_time, end_time)


def _main():
    """The hidden main block"""
    args = arguments()
    times = get_times(args.paramin)

    if not args.no_write_imf:
        write_imf_from_omni(*times)

    iono_descriptions = [
        'TypeConductanceModel',
        'UseFullCurrent',
        'UseFakeRegion2',
        'F107Flux',
        'StarLightPedConductance',
        'PolarCapPedConductance',
        ]
    if not args.no_write_f107:
        omni_lores = get_omni_data(*times, resolution='low')
        f10_7 = np.average(omni_lores['f10_7'])
        iono_command = paramin.read_command('#IONOSPHERE', args.paramin)
        change = {}
        iono_command[3] = str(f10_7)
        change[iono_command[0]] = map(
            list,
            zip(iono_command[1:], iono_descriptions)
            )
        paramin.replace_command(change, args.paramin)


if __name__ == '__main__':
    _main()
