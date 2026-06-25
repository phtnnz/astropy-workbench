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
# Version 0.4 / 2026-06-23
#       Moved and adapted to new directory structure under neoop/

VERSION     = "0.4 / 2026-06-23"
AUTHOR      = "Martin Junius"
NAME        = "neo-sbephem"
DESCRIPTION = "Ephemeris for solar system objects"

import sys
import argparse
import re

from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.time import Time
import astropy.units as u
from astroquery.mpc import MPC

# Local modules
from utils.verbose import verbose, warning, error, message
from neo.config    import config
from neo.ephem     import edata_add_ephem_jpl, edata_add_ephem_mpc, get_local_circumstances, get_dec_limits
from neo.exposure  import edata_add_exposure
from neo.classes   import EphemData
from mpc.observations import get_obs_from_mpc, get_lastobs_from_mpc, get_notseen_from_mpc

##FIXME: use config
DEFAULT_LOCATION = config.code



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
    arg.add_argument("--lastobs", action="store_true", help="output MPC obs last row")
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

        edata = EphemData("-", obj)

        # Get ephemerides
        if args.jpl:
            edata_add_ephem_jpl(edata, local)
        else:
            edata_add_ephem_mpc(edata, local)
        edata_add_exposure(edata, local)
        ic(edata)

        if edata.ephem:
            verbose.print_lines2(edata.ephem["Targetname", "Obstime", "RA", "DEC", "Mag", 
                                            "Motion", "PA", "Az", "Alt", "Moon_dist", "Moon_alt"])
            verbose("NEO exposure", edata.exposure)

        if args.obs:
            obs = get_obs_from_mpc(obj)
            verbose.print_lines2(obs)
        
        if args.lastobs:
            lastobs = get_lastobs_from_mpc(obj)
            verbose.print_lines(lastobs)
            notseen = get_notseen_from_mpc(obj)
            verbose(f"notseen = {notseen:.2f}")



if __name__ == "__main__":
    main()
