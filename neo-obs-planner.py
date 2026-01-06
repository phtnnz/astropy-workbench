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
# Version 0.1 / 2026-01-05
#       Copy of neoephem 0.1, importing neoephem functions, moved/adapted
#       functions from neocp 0.7 to neoutils, new options -s --start / -e --end

VERSION     = "0.1 / 2026-01-05"
AUTHOR      = "Martin Junius"
NAME        = "neo-obs-planner"
DESCRIPTION = "Plan NEO observations"

import sys
import argparse

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
# from astropy.coordinates import Angle
# from astropy.coordinates import errors
from astropy.time        import Time, TimeDelta
from astropy.table       import Table
# from astropy.units       import Quantity, Magnitude
import astropy.units as u
import numpy as np

from sbpy.data import Ephem
from sbpy.data import Obs
from sbpy.data.core import QueryError

from astroquery.mpc import MPC
from astroquery.exceptions import EmptyResponseError, InvalidQueryError

from astroplan import Observer

# Local modules
from verbose import verbose, warning, error, message
from astroutils import location_to_string, get_location
from neoclasses import Exposure, EphemTimes, EphemData, LocalCircumstances
from neoutils import process_obj_ephm_data, sort_obj_ephm_data
from neoconfig import config
from neoephem import get_ephem_jpl, get_ephem_mpc, get_local_circumstances

DEFAULT_LOCATION = config.code



def obs_planner_1(obj_data: dict[str, EphemData], local: LocalCircumstances) -> dict:
    # Start planner at naut. dusk / start time from options
    next_start_time = local.naut_dusk

    message("----------------------------------------------------------------------------------------")
    message("Object      Mag      Time start ephemeris   / end ephemeris                   Max motion")
    message("                     Time before            / after meridian               Moon distance")
    message("                     Time start exposure    / end exposure")
    message("                     # x Exp = total exposure time")
    message("                     RA, DEC, Alt, Az")

    for obj, edata in obj_data.items():
        exp_start = None
        exp_end = None

        etimes = edata.times
        start = etimes.start
        end = etimes.end
        total_time = edata.exposure.total_time
        alt_start = etimes.alt_start
        alt_end = etimes.alt_end
        before = etimes.before
        after = etimes.after

        message("----------------------------------------------------------------------------------------")
        message(f"{obj:9s}  {edata.mag}  {start}/{end}  {edata.motion:5.1f}")

        # check overlap with previous object
        if next_start_time > start:
            start = next_start_time
        if alt_start != None and next_start_time > alt_start:
            alt_start = next_start_time
        if alt_end != None and next_start_time > alt_end:
            alt_start = None
            alt_end = None
        if before != None and next_start_time > before:
            before = None
        if after != None and next_start_time > after:
            after = next_start_time
        if next_start_time > end:
            end = None

        ic(next_start_time.iso)
        ic(obj, start, end, total_time, alt_start, alt_end, before, after)

        # try to fit in alt_start ... alt_end interval
        if alt_start != None:
            # before meridian
            if before != None and alt_start + total_time <= before:
                exp_start = alt_start
                exp_end = alt_start + total_time
            # after meridian
            elif after != None and end != None and after + total_time <= end:
                exp_start = after
                exp_end = after + total_time
            # no meridian passing
            elif end != None and alt_start + total_time <= end:
                exp_start = alt_start
                exp_end = alt_start + total_time
        
        # no slot found, try to fit in start ... end interval
        if exp_start == None:
            # before meridian
            if before != None and start + total_time <= before:
                exp_start = start
                exp_end = start + total_time
            # after meridian
            elif after != None and end != None and after + total_time <= end:
                exp_start = after
                exp_end = after + total_time
            # no meridian passing
            elif end != None and start + total_time <= end:
                exp_start = start
                exp_end = start + total_time

        ic(exp_start, exp_end)

        if exp_start != None and exp_end != None:
            next_start_time = exp_end
            etimes.plan_start = exp_start
            etimes.plan_end = exp_end

        # Skip, if failed to allocate total_time
        if exp_start == None:
            message(f"SKIPPED: can't allocate exposure time {total_time:.2f} ({start} -- {end})")
            continue


    # end for
    message("----------------------------------------------------------------------------------------")



def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = DESCRIPTION,
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("-l", "--location", help=f"coordinates, named location or MPC station code, default {DEFAULT_LOCATION}")
    arg.add_argument("-f", "--file", help="read list of objects from file")
    arg.add_argument("-s", "--start", help="start time (UTC) (default naut. dusk)")
    arg.add_argument("-e", "--end", help="end time (UTC) (default naut. dawn)")

    arg.add_argument("-J", "--jpl", action="store_true", help="use JPL Horizons ephemeris, default MPC")
    arg.add_argument("--clear", action="store_true", help="clear MPC cache")
    arg.add_argument("object", nargs="*", help="object name")

    args = arg.parse_args()

    if args.debug:
        ic.enable()
        ic(sys.version_info, sys.path, args)
    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()

    # Observer location and local circumstances
    loc = get_location(args.location if args.location else DEFAULT_LOCATION)
    ic(loc, loc.to_geodetic())
    local = get_local_circumstances(loc)

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

    # Override start/end time
    if args.start:
        local.naut_dusk = Time(args.start)
    if args.end:
        local.naut_dawn = Time(args.end)
    verbose.print_lines(local)

    # Get ephemerides
    if args.jpl:
        obj_data = get_ephem_jpl(objects, local)
    else:
        obj_data = get_ephem_mpc(objects, local)
    # Process objects
    obj_data = process_obj_ephm_data(obj_data)
    verbose(f"original object sequence: {", ".join(obj_data.keys())}")

    for obj, edata in obj_data.items():
        verbose.print_lines(edata.ephem["Targetname", "Obstime", "RA", "DEC", "Mag", 
                                       "Motion", "PA", "Az", "Alt", "Moon_dist", "Moon_alt"])
        verbose(edata.exposure)
        verbose(edata.times)

    obj_data = sort_obj_ephm_data(obj_data)
    verbose(f"sorted object sequence: {", ".join(obj_data.keys())}")

    # Run obs planner
    obs_planner_1(obj_data, local)



if __name__ == "__main__":
    main()
