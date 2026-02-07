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
# Version 0.1 / 2026-01-06
#       Copy of neoephem 0.1, importing neoephem functions, moved/adapted
#       functions from neocp 0.7 to neoutils, new options -s --start / -e --end /,
#       -C --csv / -o --output for CSV output, first working NEO obs planning
# Version 0.2 / 2026-01-06
#       Somewhat usable now, moved CSV output to separate function

VERSION     = "0.2 / 2026-01-06"
AUTHOR      = "Martin Junius"
NAME        = "neo-obs-planner"
DESCRIPTION = "NEO observation planner"

import sys
import argparse
import csv

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.time import Time
import astropy.units as u
from astroquery.mpc import MPC

# Local modules
from verbose    import verbose, warning, error, message
from astroutils import get_location
from neoconfig  import config
from neoclasses import EphemData, LocalCircumstances
from neoutils   import obj_data_add_times, sort_obj_data, get_row_for_time, motion_limit
from neoephem   import get_ephem_jpl, get_ephem_mpc, get_local_circumstances
from neoplot    import plot_objects



DEFAULT_LOCATION = config.code
TIME_FORMAT = "%Y-%m-%d %H:%M:%S"
TIME_FORMAT_TZ = "%Y-%m-%d %H:%M:%S+0000"



def f_time(time: Time|None, add_tz: bool=False) -> str:
    return time.strftime(TIME_FORMAT_TZ if add_tz else TIME_FORMAT) if time != None else "-" * 19



def obs_planner_1(obj_data: dict[str, EphemData], local: LocalCircumstances) -> dict:
    # Start planner at naut. dusk / start time from options
    next_start_time = local.naut_dusk
    objects  = []

    message("----------------------------------------------------------------------------------")
    message("Object     Mag       Time start ephemeris/ end ephemeris                Max motion")
    message("                     Time before         / after meridian            Moon distance")
    message("                     Time start exposure / end exposure")
    message("                     # x Exp = total exposure time")
    message("                     RA, DEC, Alt, Az")

    for obj, edata in obj_data.items():
        exp_start = None
        exp_end = None

        etimes = edata.times
        start = etimes.start
        end = etimes.end
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

        # check valid exposure data, skip if None = object is too fast
        if edata.exposure:
            total_time = edata.exposure.total_time
        else:
            message(f"SKIPPED: object too fast (>{motion_limit():.1f})")
            continue
       
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
        moon_dist = row["Moon_dist"]
        min_moon_dist = config.min_moon_dist * u.degree
        if moon_dist < min_moon_dist:
            message(f"SKIPPED: moon distance {moon_dist:.0f} < {min_moon_dist:.0f}")
            continue

        # Good to go!
        # Remember end of exposure
        next_start_time = exp_end
        etimes.plan_start = exp_start
        etimes.plan_end = exp_end

        ra, dec = row["RA"].to(u.hourangle), row["DEC"]
        alt, az = row["Alt"], row["Az"]
        edata.ra, edata.dec = ra, dec
        objects.append(obj)

        message(f"                     {f_time(before)} / {f_time(after)}             {moon_dist:3.0f}")
        message(f"                     {f_time(exp_start)} / {f_time(exp_end)}")
        message(f"                     {edata.exposure}")
        message(f"                     RA {ra:.4f}, DEC {dec:.4f}, Alt {alt:.0f}, Az {az:.0f}")

    # end for
    message("----------------------------------------------------------------------------------")
    message(f"{len(objects)} object(s) planned: {", ".join(objects)}")



def obs_csv_output(obj_data: dict[str, EphemData], output: str) -> None:
    csv_rows = list()

    # Traverse objects, only those with valid plan_start time
    for obj, edata in obj_data.items():
        if edata.times.plan_start != None:
            # CSV output:
            #   start time, end time, 
            #   target=id, observation date (YYYY-MM-DD), time ut (HH:MM), ra, dec, exposure, number, filter (L),
            #   type, mag, nobs, arc, notseen, total
            # RA output  = hourangle !!!
            # DEC output = degree
            csv_row = { "target": obj,
                        "obstime": f_time(edata.times.plan_start, add_tz=True),
                        "ra": float(edata.ra.value),
                        "dec": float(edata.dec.value),
                        "exposure": float(edata.exposure.single.value),
                        "number": edata.exposure.number,
                        "filter": "L",
                        "start time": f_time(edata.times.plan_start),
                        "end time": f_time(edata.times.plan_end),
                        "type": "",
                        "mag": float(edata.mag.value),
                        "nobs": "",
                        "arc": "",
                        "notseen": "",
                        "total": str(edata.exposure)
                        }
            ic(csv_row)
            csv_rows.append(csv_row)

    # Output to CSV file for nina-create-sequence2
    verbose(f"planned objects for nina-create-sequence2: {output}")
    # csv_row is the last object, if any were found
    if csv_row:
        ##FIXME: improve csvoutput module to cover this usage
        fieldnames = csv_row.keys()
        ic(fieldnames)
        if output:
            with open(output, "w", newline="") as csvfile:
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
    arg.add_argument("-P", "--plot", action="store_true", help="create altitude and sky plot with objects")

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
    obj_data = obj_data_add_times(obj_data)
    verbose(f"original object sequence: {", ".join(obj_data.keys())}")

    for obj, edata in obj_data.items():
        verbose.print_lines(edata.ephem["Targetname", "Obstime", "RA", "DEC", "Mag", 
                                       "Motion", "PA", "Az", "Alt", "Moon_dist", "Moon_alt"])
        # verbose(edata.exposure)
        # verbose(edata.times)

    obj_data = sort_obj_data(obj_data)
    verbose(f"sorted object sequence: {", ".join(obj_data.keys())}")

    # Run obs planner
    obs_planner_1(obj_data, local)
    if args.csv:
        obs_csv_output(obj_data, args.output)

    # Plot objects and Moon
    if args.plot:
        verbose("altitude and sky plot for objects")
        plot_objects(obj_data,
                      ##FIXME##
                      "tmp/neos-plot.png",
                      loc)



if __name__ == "__main__":
    main()
