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
# Version 0.0 / 2024-xx-xx
#       TEXT

import sys
import argparse

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
import astropy.units as u
from astropy.time        import Time, TimeDelta

# Local modules
from verbose import verbose, warning, error

VERSION = "0.0 / 2024-xx-xx"
AUTHOR  = "Martin Junius"
NAME    = "astropy-template"



# Command line options
class Options:
    name = "abc"        # -n --name
    int  = 99           # -i --int



def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = "Generic python script template",
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("-n", "--name", help="example option name")
    arg.add_argument("-i", "--int", type=int, help="example option int")
    # nargs="+" for min 1 name argument
    arg.add_argument("name", nargs="*", help="filename")

    args = arg.parse_args()

    if args.debug:
        ic.enable()
        ic(sys.version_info)
        ic(args)
    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()
    # ... more options ...
    if args.name:
        Options.name = args.name
    if args.int:
        Options.int  = args.int
        
    # ... the action starts here ...



if __name__ == "__main__":
    main()