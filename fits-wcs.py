#!/usr/bin/env python

# Copyright 2024-2025 Martin Junius
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
from astropy.time import Time
import astropy.units as u
import numpy as np

# Local modules
from verbose import verbose, warning, error, message
from csvoutput import csv_output

VERSION = "0.1 / 2025-09-17"
AUTHOR  = "Martin Junius"
NAME    = "fits-wcs"



# Command line options
class Options:
    """
    Global options
    """
    hdr_list = ["NAXIS1", "NAXIS2", "DATE-OBS", "CTYPE1", "CTYPE2", "EQUINOX", "LONPOLE", "LATPOLE",
                "CUNIT1", "CUNIT2", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2",
                "CD1_1", "CD1_2", "CD2_1", "CD2_2" ]
                        # -H --header
    csv = False         # -C --csv
    output = None       # -o --output
    list = False        # -L --list
    set_locale = False  # -l --locale


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

        verbose(f"{", ".join(value_verbose)}")
        if Options.csv:
            csv_output.add_row(value_list)

        # Skip files without solution
        if not hdr.get("CRVAL1"):
            warning(f"no WCS solution in FITS headers")
            return None, None

        # See https://danmoser.github.io/notes/gai_fits-imgs.html for FITS WCS header
        obstime    = Time(hdr["DATE-OBS"])
        center_ra  = hdr["CRVAL1"] * u.degree
        center_dec = hdr["CRVAL2"] * u.degree
        center_x   = hdr["CRPIX1"]
        center_y   = hdr["CRPIX2"]
        image_w    = hdr["NAXIS1"]
        image_h    = hdr["NAXIS2"]
        cd1_1      = hdr["CD1_1"]
        cd1_2      = hdr["CD1_2"]
        cd2_1      = hdr["CD2_1"]
        cd2_2      = hdr["CD2_2"]

        # transformation matrix
        # CD1_1 =  CDELT1 * cos (CROTA2)
        # CD1_2 = -CDELT2 * sin (CROTA2)
        # CD2_1 =  CDELT1 * sin (CROTA2)
        # CD2_2 =  CDELT2 * cos (CROTA2)
        cdelt1 = np.sqrt(np.square(cd1_1) + np.square(cd2_1)) * u.degree
        cdelt2 = np.sqrt(np.square(cd1_2) + np.square(cd2_2)) * u.degree
        scale = (cdelt1 + cdelt2) / 2   # yields same value as reported by astrometry.net
        crota2_a = np.atan2(cd2_1, cd1_1) * u.radian
        crota2_b = -np.atan2(cd1_2, cd2_2) * u.radian
        ic(cdelt1, cdelt2, scale, crota2_a, crota2_b)

        ra = center_ra.to(u.hourangle)
        dec = center_dec.to(u.degree)
        crota2 = crota2_a.to(u.degree)
        scale = scale.to(u.arcsec)

        verbose(f"obstime={obstime}")
        verbose(f"center RA={ra:.4f}, DEC={dec:.4f}, pos={center_x}({image_w})/{center_y}({image_h}), {scale:.2f}/pxl")
        verbose(f"transformation matrix: {cd1_1:+.8f} {cd1_2:+.8f}")
        verbose(f"                       {cd2_1:+.8f} {cd2_2:+.8f}")
        verbose(f"rotation angle={crota2:.2f}")

        if Options.list:
            print(f"{obstime}  {ra:.4f}  {dec:+.4f}  {crota2:.4f}")

        return obstime, crota2



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
            obstime_1st  = None
            rot_1st      = None
            obstime_last = None
            rot_last     = None
            for file in file_list:
                file = os.path.join(dir, file)
                if file.lower().endswith(".fits") or file.lower().endswith(".fit"):
                    obstime, rot = process_fits(file)
                    if obstime != None:
                        if obstime_1st == None:
                            obstime_1st = obstime
                            rot_1st     = rot
                        obstime_last = obstime
                        rot_last     = rot

        if obstime_1st != None and obstime_last != None:
            delta_time = obstime_last - obstime_1st
            delta_hour  = delta_time.to_value("sec") / 3600 * u.hour
            delta_rot  = rot_last     - rot_1st
            rel_rot    = delta_rot / delta_hour
            if Options.list:
                print(f"Rotation: {rel_rot:.3f}")

    else:
        error(f"no such file or directory {name}")



def init_csv_output():
    """
    Initialize CSV output
    """
    if Options.csv:
        csv_output.add_fields(Options.hdr_list)
        csv_output.set_float_format("%.6f")

def write_csv_output():
    """
    Write CSV output
    """
    if Options.csv:
        csv_output.write(Options.output, Options.set_locale)



def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = "List FITS file WCS headers",
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("-H", "--header", help=f"show image headers in HEADER list, \"+\" adds, (default {",".join(Options.hdr_list)})")
    arg.add_argument("-L", "--list", action="store_true", help="list header keywords")
    arg.add_argument("-C", "--csv", action="store_true", help="output CSV list")
    arg.add_argument("-l", "--locale", action="store_true", help="set locale for CSV output")
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
        h = args.header
        if h.startswith("+"):
            h = h.lstrip("+")
            h = h.lstrip(",")
            Options.hdr_list.extend(h.split(","))
        else:
            Options.hdr_list = h.split(",")
    Options.csv = args.csv
    Options.output = args.output
    Options.list = args.list
    Options.set_locale = args.locale

    init_csv_output()
    for fits in args.fits:
        # quick hack: Windows PowerShell adds a stray " to the end of dirname 
        # if it ends with a backslash \ AND contains a space!!!
        # see here https://bugs.python.org/issue39845
        process_file_or_dir(fits.rstrip("\""))
    write_csv_output()



if __name__ == "__main__":
    main()