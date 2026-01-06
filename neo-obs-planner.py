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
import csv

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
from neoutils import process_obj_ephm_data, sort_obj_ephm_data, get_row_for_time
from neoconfig import config
from neoephem import get_ephem_jpl, get_ephem_mpc, get_local_circumstances

DEFAULT_LOCATION = config.code
TIME_FORMAT = "%Y-%m-%d %H:%M:%S"


# Command line options
class Options:
    """
    Command line options
    """
    csv: bool = False           # -C --csv
    output: str = None          # -o --output



def f_time(time: Time|None) -> str:
    return time.strftime(TIME_FORMAT) if time != None else "-" * 19



def obs_planner_1(obj_data: dict[str, EphemData], local: LocalCircumstances) -> dict:
    # Start planner at naut. dusk / start time from options
    next_start_time = local.naut_dusk

    message("----------------------------------------------------------------------------------")
    message("Object      Mag      Time start ephemeris/ end ephemeris                Max motion")
    message("                     Time before         / after meridian            Moon distance")
    message("                     Time start exposure / end exposure")
    message("                     # x Exp = total exposure time")
    message("                     RA, DEC, Alt, Az")

    ##FIXME: separate CSV output
    csv_rows = list()

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

        message("----------------------------------------------------------------------------------")
        message(f"{obj:9s}  {edata.mag}  {f_time(start)} / {f_time(end)}  {edata.motion:5.1f}")

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

        # Skip, if failed to allocate total_time
        if exp_start == None:
            message(f"SKIPPED: can't allocate exposure time {total_time:.2f} ({start} -- {end})")
            continue

        # Ephemeris row best matching start time
        row = get_row_for_time(edata.ephem, exp_start)
        ic(row)

        # Skip, if percentage of total exposure time is less than threshold
        if edata.exposure.percentage < config.min_perc_required:
            message(f"SKIPPED: only {edata.exposure.percentage:.0f}% of required total exposure time (< {config.min_perc_required}%)")
            continue

        # Skip, if moon distance is too small
        moon_dist = row["Moon_dist"][0]
        min_moon_dist = config.min_moon_dist * u.degree
        if moon_dist < min_moon_dist:
            message(f"SKIPPED: moon distance {moon_dist:.0f} < {min_moon_dist:.0f}")
            continue

        # Good to go!
        # Remember end of exposure
        next_start_time = exp_end
        etimes.plan_start = exp_start
        etimes.plan_end = exp_end

        ra, dec = row["RA"][0], row["DEC"][0]
        alt, az = row["Alt"][0], row["Az"][0]

        total = f"{edata.exposure.number} x {edata.exposure.single:2.0f} = {edata.exposure.total:3.1f} ({edata.exposure.percentage:.0f}%) / total {edata.exposure.total_time:3.1f}"
        message(f"                     {f_time(before)} / {f_time(after)}             {moon_dist:3.0f}")
        message(f"                     {f_time(exp_start)} / {f_time(exp_end)}")
        message(f"                     {total}")
        message(f"                     RA {ra:.4f}, DEC {dec:.4f}, Alt {alt:.0f}, Az {az:.0f}")

        ##FIXME: store output data in EphemData, separate function###########################################
        # CSV output:
        #   start time, end time, 
        #   target=id, observation date (YYYY-MM-DD), time ut (HH:MM), ra, dec, exposure, number, filter (L),
        #   type, mag, nobs, arc, notseen, total
        csv_row = { "target": obj,
                    # "+0000" to enforce UTC
                    "obstime": exp_start.strftime("%Y-%m-%d %H:%M:%S+0000"),
                    "ra": float(ra.value),
                    "dec": float(dec.value),
                    "exposure": float(edata.exposure.single.value),
                    "number": edata.exposure.number,
                    "filter": "L",
                    "start time": str(exp_start),
                    "end time": str(exp_end),
                    "type": "",
                    "mag": float(edata.mag.value),
                    "nobs": "",
                    "arc": "",
                    "notseen": "",
                    "total": total
                    }
        ic(csv_row)
        csv_rows.append(csv_row)
        #####################################################################################################


    # end for
    message("----------------------------------------------------------------------------------")

    # Output to CSV file for nina-create-sequence2
    if Options.csv:
        verbose(f"planned objects for nina-create-sequence2: {Options.output}")
        # csv_row is the last object, if any were found
        if csv_row:
            ##FIXME: improve csvoutput module to cover this usage
            fieldnames = csv_row.keys()
            ic(fieldnames)
            if Options.output:
                with open(Options.output, "w", newline="") as csvfile:
                    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                    writer.writeheader()
                    writer.writerows(csv_rows)
            else:
                writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(csv_rows)
        else:
            warning("no objects, no CSV output")



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
    arg.add_argument("-o", "--output", help="write CSV to OUTPUT file")
    arg.add_argument("-C", "--csv", action="store_true", help="use CSV output format")

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

    Options.csv    = args.csv
    Options.output = args.output

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
        # verbose(edata.exposure)
        # verbose(edata.times)

    obj_data = sort_obj_ephm_data(obj_data)
    verbose(f"sorted object sequence: {", ".join(obj_data.keys())}")

    # Run obs planner
    obs_planner_1(obj_data, local)



if __name__ == "__main__":
    main()
