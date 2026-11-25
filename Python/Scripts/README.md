Python Scripts
==============

Collection of Python scripts to make using [SWMF](https://github.com/SWMFsoftware/SWMF) easier. Note most of these require a recent version of Python __3__. Supercomputers typically have this already.

Make sure [swmfpy](https://github.com/SWMFsoftware/swmfpy) is installed.

```bash
$ pip install swmfpy || pip install --user swmfpy
```

Table of Contents:
------------------

- [prepare_geospace.py](#prepare_geospacepy)
- [plot_indices.py](#plot_indices)
- [download_sme_kp.py](#download_sme_kppy)
- [swmf_datafile.py](#swmf_datafilepy)

prepare_geospace.py
-------------------

This script is to help with the inputs of the geospace model.

### Pre-steps

If you already have a run directory and `PARAM.in` file then skip ahead.

In the base directory of SWMF run:

```bash
SWMF$ ./Config.pl -install
# Compile then make the run directory if worked
SWMF$ make -j test_swpc_compile && make rundir
```

### Set up

It would be useful to set up an environment variable to this script location.
```bash
SWMF$ echo 'export PYSCRIPTS='"$(realpath share/Python/Scripts)" >> ~/.profile
SWMF$ . ~/.profile
```

Go to your run directory and copy a sensible `PARAM.in`:

```bash
SWMF$ cd run
SWMF/run$ cp Param/SWPC/PARAM.in_SWPC_v2_init PARAM.in
```

### Running

Then the final steps copy the script and the run directory into a work directory that your supercomputer allows before jobs:

```bash
SWMF$ cp run /some/work/dir/
SWMF$ cd /some/work/dir/
/some/work/dir$ python3 "$PYSCRIPTS"/prepare_geospace.py --start_time 2014 2 3 4 5 6 --end_time 2014 3 4 5 6 7
# Your PARAM.in will be overwritten with those values
# Then submit job
```

### Options

To find out the options of `prepare_geospace.py` then run:

```bash
/some/work/dir$ ./prepare_geospace.py --help
# Help message output
```

Options can be useful to script this file with shell scripting, for example, when submitting jobs.

plot_indices.py
---------------

A script to plot the global indices of your geospace run (SYM-H, cross polar cap
potential, AL.. etc.). This uses `swmfpy` to get the observed indices to
compare.

To run simply type on your command line:

```bash
python plot_indices.py -h
```

A simply way of using this script is to make sure you know where your
`log_****.log` and corresponding `geo_****.log` are and running:

```bash
# replace with your file names
python plot_indices.py log_YYYYMMMDD.log geo_YYYMMMDD.log
```

download_sme_kp.py
------------------

A script to download the Kp index and the SME/SML/SMU (SuperMAG Electrojet)
indices for a given UTC time period, saving each to its own CSV file. It only
uses the Python standard library, so no extra packages need to be installed.

Kp is fetched from GFZ Potsdam's web service, which needs no registration.
SME/SML/SMU are SuperMAG-exclusive products, so fetching them requires a free
SuperMAG "logon" username (just register with an email address at
https://supermag.jhuapl.edu, no API key needed).

Run once per time period:

```bash
python download_sme_kp.py --start 2024-05-10T00:00:00 --end 2024-05-13T00:00:00 \
    --userid YOUR_SUPERMAG_LOGON --out period1
```

This writes `period1_kp.csv` (3-hour Kp) and `period1_sme.csv` (1-minute
SME/SML/SMU). The SuperMAG username can also be set once via the
`SUPERMAG_USERID` environment variable instead of `--userid`; if neither is
given, only the Kp file is downloaded. Run `python download_sme_kp.py --help`
for all options.

swmf_datafile.py
----------------

A module to read and write SWMF data files (ASCII, `real4` and `real8`
binary formats) efficiently with `numpy`. Import it into your own Python
code to access the snapshots:

```python
import swmf_datafile as swmf

data = swmf.read_file("file.out")   # returns a dictionary, see below
swmf.show_data(data)                # print a summary of the dictionary
swmf.write_file(data, "file2.out")  # write it back out
```

After reading, the dictionary contains entries such as `head`, `step`,
`time`, `ndim`, `cart`, `npar`, `nvar`, `dims`, `pars`, `name`, `coord` and
`state` (see the top of `swmf_datafile.py` for the full description). For
writing, the dictionary must at least contain `name`, `coord` and `state`.

`swmf_datafile_test.py` is a self-contained example/test script that builds
sample data dictionaries, writes them out in all supported formats, and reads
them back, printing the content at each step. Run it directly to see example
usage, or compare its input/output to understand the data format:

```bash
python swmf_datafile_test.py
```

`swmf_datafile.py` can also be run directly from the command line to
inspect, select snapshots from, and convert the format of existing SWMF data
files:

```bash
# Check the content of file.outs
swmf_datafile.py file.outs

# Show the header info for the first snapshot
swmf_datafile.py -v -s 0 file.outs

# Extract the first snapshot and then every 10th from 30th to 200th
swmf_datafile.py -s 0,30:200:10 file.outs selected.outs

# Convert the file into an ascii file with 3 digit floats
swmf_datafile.py -f 10.3f file.outs ascii.out

# Convert the file into a single precision binary file
swmf_datafile.py -f real4 file.outs singleprec.outs

# Extract the last two snapshots (note the = sign!) and convert to real8 binary
swmf_datafile.py -f real8 -s=-2: file.outs doubleprec.out
```

Run `swmf_datafile.py --help` for all options.
