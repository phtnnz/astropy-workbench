#!/usr/bin/env python

# Copyright 2026 Martin Junius
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
# Version 0.0 / 2026-11-09
#       Read Autoslew POX files

VERSION     = "0.0 / 2025-11-09"
AUTHOR      = "Martin Junius"
NAME        = "pox"
DESCRIPTION = "Read Autoslew POX files"

import sys
import argparse
import re
from dataclasses import dataclass

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.coordinates import Angle
from astropy.coordinates import errors
from astropy.time        import Time, TimeDelta
import astropy.units as u
import astropy.constants as const
import numpy as np

# Local modules
from verbose import verbose, warning, error, message
from astroutils import location_to_string, get_location     # ...

DEFAULT_LOCATION = "M49"



# Command line options
class Options:
    loc: EarthLocation = None       # -l --location



class POX:
    # From ASA Tool NINA plugin:
    # public int Number { get; set; }
    # public string DateObs { get; set; }
    # public string TimeObs { get; set; }
    # public double ExpTime { get; set; }
    # public double TelescopeRA { get; set; }
    # public double SolvedRA { get; set; }
    # public double TelescopeDec { get; set; }
    # public double SolvedDec { get; set; }
    # public int PierSide { get; set; }
    pass

class POXList:
    def __init__(self):
        self.obstime = list()
        self.ra_telescope = list()
        self.ra_solved = list()
        self.dec_telescope = list()
        self.dec_solved = list()


    @classmethod
    def from_pox_file(cls, filename: str):
        self = cls()

        verbose(f"processing {filename}")
        with open(filename, "r") as file:
            total_number = None
            number = 0
            file_iter = iter(file)

            total_number = int( next(file_iter).strip("\n\"\'") )
            ic(total_number)

            while ( line := next(file_iter, "").strip("\n\"\'") ):
                # Number line
                ic(line)
                m = re.match(r'Number (\d+)', line)
                if not m:
                    warning(f"format error, reading stopped: {line}")
                    return self
                number = int(m.group(1))
                ic(number)

                # Obstime
                obstime = Time( next(file_iter, "").strip("\n\"\'") )
                # Time, ignored
                next(file_iter, "")
                # Exposure time, ignored
                next(file_iter, "")
                # RA telescope
                ra_telescope = float( next(file_iter, "").strip("\n\"\'") )
                # RA solved
                ra_solved = float( next(file_iter, "").strip("\n\"\'") )
                # DEC telescope
                dec_telescope = float( next(file_iter, "").strip("\n\"\'") )
                # DEC solved
                dec_solved = float( next(file_iter, "").strip("\n\"\'") )
                # Side of pier, ignored
                next(file_iter, "")

                # Check separator
                line = next(file_iter, "").strip("\n\"\'")
                m = re.match(r'\*+', line)
                if not m:
                    warning(f"format error, reading stopped: {line}")
                    return self

                ic(obstime, ra_telescope, ra_solved, dec_telescope, dec_solved)
                self.obstime.append(obstime)
                self.ra_telescope.append(ra_telescope)
                self.ra_solved.append(ra_solved)
                self.dec_telescope.append(dec_telescope)
                self.dec_solved.append(dec_solved)
        # /with
        return self


    def to_skycoord(self) -> tuple[SkyCoord, SkyCoord]:
        return ( SkyCoord(self.ra_telescope, self.dec_telescope, unit=(u.hour, u.deg)) ,
                 SkyCoord(self.ra_solved, self.dec_solved, unit=(u.hour, u.deg))        )



def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = DESCRIPTION,
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("-l", "--location", help=f"coordinates, named location or MPC station code, default {DEFAULT_LOCATION}")
    arg.add_argument("file", nargs="+", help="POX data file(s)")
    args = arg.parse_args()

    if args.debug:
        ic.enable()
        ic(sys.version_info, sys.path, args)
    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()

    # Observer location
    loc = get_location(args.location if args.location else DEFAULT_LOCATION)
    ic(loc, loc.to_geodetic())
    Options.loc = loc
    verbose(f"location {location_to_string(loc)}")

    # POX files
    for filename in args.file:
        pox = POXList.from_pox_file(filename)

        coord, coord_solved = pox.to_skycoord()
        ic(coord, coord_solved)



if __name__ == "__main__":
    main()
