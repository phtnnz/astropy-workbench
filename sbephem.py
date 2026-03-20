#!/usr/bin/env python

# Copyright 2025-2026 Martin Junius
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
# Version 0.1 / 2025-12-23
#       Somewhat usable now ;-)
# Version 0.2 / 2026-01-03
#       Common column names for ephemeris tables
# Version 0.3 / 2026-03-15
#       Use functions from neoephem, with slightly different user interface

VERSION     = "0.3 / 2026-03-15"
AUTHOR      = "Martin Junius"
NAME        = "sbephem"
DESCRIPTION = "Ephemeris for solar system objects"

import sys
import argparse
import re

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.time import Time
import astropy.units as u
from astroquery.mpc import MPC
from sbpy.data import Obs
from sbpy.data.core import QueryError

# Local modules
from verbose    import verbose, warning, error, message
from neoconfig  import config
from neoephem   import get_ephem_jpl, get_ephem_mpc, get_local_circumstances, get_dec_limits

##FIXME: use config
DEFAULT_LOCATION = config.code



def id_type_from_name(name: str) -> str:
    id_type_regex = {   "asteroid numer":         "^[1-9][0-9]*$",
                        "asteroid designation":   "^20[0-9]{2}[ _][A-Z]{2}[0-9]{0,3}$",
                        "comet number":           "^[0-9]{1,3}[PIA]$",
                        "comet designation":      "^[PDCXAI]/20[0-9]{2}[ _][A-Z]{2}[0-9]{0,3}$"
                    }

    for id in id_type_regex.keys():
        m = re.match(id_type_regex[id], name)
        if m:
            ic(name, id)
            return id
    ## Default None or "asteroid designation"?
    return None



def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = DESCRIPTION,
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("-l", "--location", help=f"coordinates, named location or MPC station code, default {DEFAULT_LOCATION}")
    arg.add_argument("-f", "--file", help="read list of objects from file")
    arg.add_argument("-t", "--time", help="start time for ephemeris (1h, 5min steps)")
    arg.add_argument("-J", "--jpl", action="store_true", help="use JPL Horizons ephemeris, default MPC")
    arg.add_argument("-a", "--allnight", action="store_true", help="ephemeris for midnight +/- 8h (30min steps)")
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

    # Set min altitude for ephemeris to 0 deg
    config.min_alt = 0
    config.elev_min = 0

    # Observer location and local circumstances
    local = get_local_circumstances(args.location if args.location else DEFAULT_LOCATION)

    # Override config DEC limits
    min_dec, max_dec = get_dec_limits(local, config.min_alt*u.deg)
    config.min_dec = int(min_dec.degree)
    config.max_dec = int(max_dec.degree)

    # Observation time
    time = Time(args.time) if args.time else Time.now()
    ic(time)
    verbose(f"time {time.iso} ({time.scale.upper()})")

    if not args.allnight:
        local.epochs = {"start":  time,
                        "step":   5 * u.min,
                        "stop":   time + 1 * u.hour
                        }
        ##FIXME: better solution?
        local.naut_dusk = time
        local.naut_dawn = time + 1 * u.hour
    ic(local.epochs)
    verbose.print_lines(local)
    verbose("listing ephemeris for altitude>0 and after nautical dusk/before nautical dawn only!")

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

        # Get ephemerides
        if args.jpl:
            obj_data = get_ephem_jpl(objects, local, type)
        else:
            obj_data = get_ephem_mpc(objects, local, type)

        for obj, edata in obj_data.items():
            verbose.print_lines(edata.ephem["Targetname", "Obstime", "RA", "DEC", "Mag", 
                                        "Motion", "PA", "Az", "Alt", "Moon_dist", "Moon_alt"])

            ##CHECK: doesn't work properly
            if args.obs:
                try:
                    obs = Obs.from_mpc(obj, id_type=id_type_from_name(obj))
                    print(obs)
                    # # Handle masked entries
                    # for i in range(-1, -10, -1):
                    #     mag = obs["mag"][i].unmasked
                    #     if mag > Magnitude(0):
                    #         break
                except QueryError as e:
                    warning(f"MPC observations for {obj} failed")
                except ConnectionError as e:
                    warning(f"MPC request failed: {e}")



if __name__ == "__main__":
    main()
