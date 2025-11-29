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
from astropy.units import Quantity, Magnitude
import astropy.units as u
import numpy as np
from sbpy.data import Ephem
from sbpy.data import Obs
from astroquery.mpc import MPC
from astroquery.exceptions import EmptyResponseError, InvalidQueryError

# Local modules
from verbose import verbose, warning, error, message
from astroutils import location_to_string, get_location     # ...
from neoutils import Exposure, exposure_from_ephemeris, id_type_from_name

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
    arg.add_argument("-J", "--jpl", action="store_true", help="use JPL Horizons ephemeris, default MPC")
    arg.add_argument("--obs", action="store_true", help="output MPC obs")
    arg.add_argument("--clear", action="store_true", help="clear MPC cache")
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
    if args.file:
        with open(args.file, "r") as file:
            for line in file:
                objects.append(line.strip())
    if args.object:
        objects.extend(args.object)
    ic(objects)
    if not objects:
        error("no objects from file or command line")

    # Clear astroquery cache
    if args.clear:
        MPC.clear_cache()    

    for obj in objects:
        verbose(f"object {obj}")

        # # Object ephemeris
        # try:
        #     eph = MPC.get_ephemeris(obj, location=loc, start=time, step=5*u.min, number=12)
        #     print(eph["Date", "RA", "Dec", "V", "Proper motion", "Direction", "Azimuth", "Altitude"])
        # except InvalidQueryError as err:
        #     warning(f"query MPC ephemeris for failed {obj}!")
        #     eph = None

        # # Observations
        # try:
        #     ##FIXME: must specify id_type
        #     obs = MPC.get_observations(obj, id_type="comet number")
        #     print(obs)
        # except (EmptyResponseError, ValueError) as err:
        #     warning(f"query MPC observations failed for {obj}!")
        #     obs = None

        # Get ephemerides via sbpy
        epochs = {"start":  time,
                  "step":   5 * u.min,
                  "stop":   time + 1 * u.hour
                  }
        if args.jpl:
            eph = Ephem.from_horizons(obj, location=loc, epochs=epochs)
            ic(eph.field_names)
            mag_col = "Tmag" if "Tmag" in eph.field_names else "V"
            mag = eph[mag_col][0]
            print(eph["targetname", "epoch", "solar_presence", "lunar_presence", "RA", "DEC", 
                    "RA*cos(Dec)_rate", "DEC_rate", "AZ", "EL", mag_col, "velocityPA"])
            exp = exposure_from_ephemeris(eph, "RA*cos(Dec)_rate,DEC_rate", mag)
            print(exp)
        else:
            # eph = Ephem.from_mpc(obj, location=loc, epochs=epochs)
            eph = Ephem.from_mpc(obj, location=loc, epochs=epochs, 
                                ra_format={'sep': ':', 'unit': 'hourangle', 'precision': 1}, 
                                dec_format={'sep': ':', 'precision': 1})
            ic(eph.field_names)
            mag = eph["V"][0]
            print(eph["Targetname", "Date", "RA", "Dec", "V", "Proper motion", "Direction", "Azimuth", "Altitude"])
            exp = exposure_from_ephemeris(eph, "Proper motion", mag)
            print(exp)

            if args.obs:
                obs = Obs.from_mpc(obj, id_type=id_type_from_name(obj))
                print(obs)
                # # Handle masked entries
                # for i in range(-1, -10, -1):
                #     mag = obs["mag"][i].unmasked
                #     if mag > Magnitude(0):
                #         break



if __name__ == "__main__":
    main()
