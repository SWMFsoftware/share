#!/usr/bin/env python
"""Prepares a geospace run by downloading inputs.

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
DESC = """It makes the IMF input file (IMF.dat) and also replaces the F10.7
value in the PARAM.in file.
This uses an *already prepared* PARAM.in files to work with. If start and end
times aren't specified then it reads it from the PARAM.in file.
Please run in the same directory as your PARAM.in file (run dir).
"""
__doc__ += DESC


def arguments():
    """Make arguments"""
    parser = argparse.ArgumentParser(description=DESC)

    metavar_time = ('y', 'm', 'd', 'h', 'm', 's')

    parser.add_argument(
        '-t0', '--start_time',
        help="""The start time of the solar wind.
This well replace PARAM.in values.
If times are not specified this reads the PARAM.in file.""",
        nargs=6,
        type=int,
        metavar=metavar_time,
        default=None
        )

    parser.add_argument(
        '-t1', '--end_time',
        help="""The end time of the solar wind.
This well replace PARAM.in values.
If times are not specified this reads the PARAM.in file.""",
        nargs=6,
        type=int,
        metavar=metavar_time,
        default=None
        )

    parser.add_argument(
        '-v', '--verbose',
        help="Output verbose messages.",
        action="store_true",
        )

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
        metavar='PARAM.in',
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
    start_time = dt.datetime(*list(map(int, time_command[1:6])))

    time_command = paramin.read_command('#ENDTIME',
                                        paramin_file)
    end_time = dt.datetime(*list(map(int, time_command[1:6])))

    return (start_time, end_time)


def vprint(arg, *msg):
    "Verbose print of message"
    if arg.verbose:
        print(*msg)


def write_times(times, paramin_file='PARAM.in'):
    """Changes the times in the paramin file
    """
    change = {"#STARTTIME": [], "#ENDTIME": []}
    time_desc = [
        'year',
        'month',
        'day',
        'hour',
        'minute',
        'second',
        ]

    def time_as_list(thetime):
        'return time as list of str'
        raw = [thetime.year,
               thetime.month,
               thetime.day,
               thetime.hour,
               thetime.minute,
               thetime.second,
               ]
        return list(map(str, raw))

    change['#STARTTIME'] = zip(time_as_list(times[0]),
                               time_desc)
    change['#ENDTIME'] = zip(time_as_list(times[1]),
                             time_desc)
    paramin.replace_command(change, paramin_file)


def _main():
    """The hidden main block"""
    args = arguments()

    # Time handling
    if args.start_time and args.end_time:
        start_time = dt.datetime(*args.start_time)
        end_time = dt.datetime(*args.end_time)
        times = (start_time, end_time)
        write_times(times, args.paramin)
    elif args.start_time or args.end_time:
        raise ValueError('Please specify both start and end times')
    else:
        times = get_times(args.paramin)
    vprint(args, 'Start time:', times[0], '| End Time:', times[1])

    if not args.no_write_imf:
        write_imf_from_omni(*times)
        vprint(args, 'IMF.dat written')

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
        vprint(args, 'F10.7:', f10_7)
        iono_command = paramin.read_command('#IONOSPHERE', args.paramin)
        change = {}
        iono_command[3] = str(f10_7)
        change[iono_command[0]] = map(
            list,
            zip(iono_command[1:], iono_descriptions)
            )
        paramin.replace_command(change, args.paramin)
        vprint(args, 'Wrote F10.7 to PARAM.in file')


if __name__ == '__main__':
    _main()
