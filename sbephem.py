#!/usr/bin/env python

# Copyright 2025 Martin Junius
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
# Version 0.0 / 2025-11-22
#       Get ephemeris for solar system objects

VERSION     = "0.0 / 2025-11-22"
AUTHOR      = "Martin Junius"
NAME        = "sbephem"
DESCRIPTION = "Ephemeris for solar system objects"

import sys
import argparse
from typing import Tuple, Any

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
from sbpy.data import Ephem

# Local modules
from verbose import verbose, warning, error, message
from astroutils import location_to_string, get_location     # ...

DEFAULT_LOCATION = "M49"



# Command line options
class Options:
    loc: EarthLocation = None       # -l --location



def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = DESCRIPTION,
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("-l", "--location", help=f"coordinates, named location or MPC station code, default {DEFAULT_LOCATION}")
    arg.add_argument("-f", "--file", help="read list of objects from file")
    arg.add_argument("-t", "--time", help="time for ephemeris")
    arg.add_argument("object", nargs="*", help="object name")

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

    # Observation time
    time = Time(args.time) if args.time else Time.now()
    ic(time)
    verbose(f"time {time} ({time.scale.upper()})")

    # Objects
    objects = []
    ##ADD: read from file
    if args.object:
        objects.append(args.object)
    ic(objects)
    if not objects:
        error("no objects from file or command line")

    for obj in objects:
        verbose(f"ephemeris for {obj}")
        epochs = {"start":  time,
                  "step":   5 * u.min,
                  "stop":   time + 1 * u.hour
                  }
        eph = Ephem.from_horizons(obj, location=loc, epochs=epochs)
        ic(eph._table.info)
        print(eph["targetname", "epoch", "solar_presence", "lunar_presence", "RA", "DEC", 
                  "RA*cos(Dec)_rate", "DEC_rate", "AZ", "EL", "Tmag", "Nmag", "velocityPA"])
        eph = Ephem.from_mpc(obj, location=loc, epochs=epochs)
        ic(eph._table.info)
        print(eph["Targetname", "Date", "RA", "Dec", "V", "Proper motion", "Direction", "Azimuth", "Altitude"])



if __name__ == "__main__":
    main()
