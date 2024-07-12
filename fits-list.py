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
# Version 0.0 / 2024-07-12
#       First version

import sys
import argparse

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
# from astropy.coordinates import SkyCoord  # High-level coordinates
# from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
# from astropy.coordinates import Angle, Latitude, Longitude  # Angles
# import astropy.units as u
# from astropy.time        import Time, TimeDelta

# See docs at https://docs.astropy.org/en/stable/io/fits/
from astropy.io import fits

# Local modules
from verbose import verbose, warning, error

VERSION = "0.0 / 2024-07-12"
AUTHOR  = "Martin Junius"
NAME    = "fits-list"



# Command line options
class Options:
    hdr_list = ["DATE-LOC", "IMAGETYP", "EXPOSURE"]     # -H --header



def process_fits(file):
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
    


def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = "List FITS file image headers",
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("-H", "--header", help=f"show image headers in HEADER list, (default {",".join(Options.hdr_list)})")
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

    # ... the action starts here ...
    for fits in args.fits:
        process_fits(fits)



if __name__ == "__main__":
    main()