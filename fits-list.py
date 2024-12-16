#!/usr/bin/env python

# Copyright 2024 Martin Junius
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# ChangeLog
# Version 0.1 / 2024-12-16
#       First version

import sys
import argparse
import os

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# See docs at https://docs.astropy.org/en/stable/io/fits/
from astropy.io import fits

# Local modules
from verbose import verbose, warning, error
from csvoutput import csv_output

VERSION = "0.1 / 2024-12-16"
AUTHOR  = "Martin Junius"
NAME    = "fits-list"



# Command line options
class Options:
    """
    Global options
    """
    hdr_list = ["OBJECT", "DATE-OBS", "RA", "DEC", "CENTALT", "CENTAZ", "IMAGETYP", "EXPOSURE"]
                    # -H --header
    csv = False     # -C --csv
    output = None   # -o --output



def process_fits(file: str):
    """
    Process single FITS file

    :param file: file name
    :type file: str
    """
    with fits.open(file) as hdul:
        if ic.enabled:
            hdul.info()
        elif verbose.enabled:
            verbose(f"FITS file {file}")
        hdr = hdul[0].header
        ic(list(hdr.keys()))

        value_verbose = [ f"{h}={hdr.get(h)}" for h in Options.hdr_list ]
        value_list = [ hdr.get(h) for h in Options.hdr_list ]

        verbose(f"  {", ".join(value_verbose)}")
        if Options.csv:
            csv_output.add_row(value_list)



def process_file_or_dir(name: str):
    """
    Process single FITS file or traverse directory

    :param name: file or directory name
    :type name: str
    """
    if os.path.isfile(name):
        process_fits(name)

    elif os.path.isdir(name):
        for dir, subdir_list, file_list in os.walk(name):
            verbose(f"found directory {dir}")
            for file in file_list:
                file = os.path.join(dir, file)
                if file.lower().endswith(".fits") or file.lower().endswith(".fit"):
                    process_fits(file)

    else:
        error(f"no such file or directory {name}")



def init_csv_output():
    """
    Initialize CSV output
    """
    if Options.csv:
        csv_output.set_default_locale()
        csv_output.add_fields(Options.hdr_list)

def write_csv_output():
    """
    Write CSV output
    """
    if Options.csv:
        csv_output.write(Options.output)



def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = "List FITS file image headers",
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("-H", "--header", help=f"show image headers in HEADER list, (default {",".join(Options.hdr_list)})")
    arg.add_argument("-C", "--csv", action="store_true", help="output CSV list")
    arg.add_argument("-o", "--output", help="write CSV to file OUTPUT (default: stdout)")
    arg.add_argument("fits", nargs="+", help="FITS filename or directory")

    args = arg.parse_args()

    if args.debug:
        ic.enable()
        ic(sys.version_info)
        ic(args)
    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()
    # ... more options ...
    if args.header:
        Options.hdr_list = args.header.split(",")     
    Options.csv = args.csv
    Options.output = args.output

    # ... the action starts here ...
    init_csv_output()
    for fits in args.fits:
        # quick hack: Windows PowerShell adds a stray " to the end of dirname 
        # if it ends with a backslash \ AND contains a space!!!
        # see here https://bugs.python.org/issue39845
        process_file_or_dir(fits.rstrip("\""))
    write_csv_output()



if __name__ == "__main__":
    main()