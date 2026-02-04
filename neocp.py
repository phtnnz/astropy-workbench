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
# Version 0.0 / 2025-08-25
#       First attempt at parsing MPC NEOCP ephemerides
# Version 0.1 / 2025-09-04
#       Retrieval from MPC, altitude and sky plots
# Version 0.2 / 2025-09-12
#       Parse NEOCP and PCCP lists for additional data, separate altitude/sky plots
# Version 0.3 / 2025-09-17
#       Added CSV output for nina-create-sequence2
# Version 0.4 / 2025-09-18
#       Major improvement to planning, avoid overlap of observations,
#       altitude/sky plot only for planned objects
# Version 0.5 / 2025-09-25
#       Using config files neocp.json for all parameters
# Version 0.6 / 2025-10-03
#       Added -M --mag-limit option to override config settings,
#       added -p --prefix option for cached data, check against
#       motion limit derived from minimum exposure time
# Version 0.7 / 2025-10-15
#       Joined altitude/sky plot, new option -P --plot
# Version 0.8 / 2026-02-03
#       Started code refactoring and adaption to new data structures

VERSION = "0.8 / 2026-02-03"
AUTHOR  = "Martin Junius"
NAME    = "neocp"

import sys
import os
import argparse
import csv

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.coordinates import Angle, EarthLocation
import astropy.units as u
from astropy.units import Quantity, Magnitude
from astropy.time import Time
from astropy.table import QTable, Row

# Astroplan
from astroplan import Observer

# Local modules
from verbose import verbose, warning, error, message
from astroutils import mpc_station_location, location_to_string
from neoconfig import config
from neoplot import plot_objects
from mpcneocp import mpc_query_neocp_ephemerides, mpc_query_neocp_list, parse_html_ephemerides, parse_neocp_list, convert_text_ephemerides, print_ephemerides
from neoclasses import Exposure



# Requests timeout
TIMEOUT = config.requests_timeout
# Exposure times / s
EXP_TIMES = config.exposure_times



# Command line options
class Options:
    """
    Command line options
    """
    csv: bool = False           # -C --csv
    output: str = None          # -o --output
    code: str = config.code     # -l --location
    loc: EarthLocation = mpc_station_location(code)



##FIXME: replace wrappers with new code ##
from neoutils import Ephem
from neoutils import single_exp, motion_limit, exposure_calc

from neoutils import max_motion as _max_motion
from neoutils import flip_times as _flip_times
from neoutils import max_alt_time as _max_alt_time
from neoutils import get_row_for_time as _get_row_for_time
from neoutils import opt_alt_times as _opt_alt_times

## Wrapper for new functions ##
def max_motion(qt: QTable) -> Quantity:
    ephem = Ephem.from_table(qt)
    return _max_motion(ephem, "motion")

def flip_times(qt: QTable) -> tuple[Time, Time]:
    ephem = Ephem.from_table(qt)
    return _flip_times(ephem, "obstime", "az")

def max_alt_times(qt: QTable) -> tuple[Time, Time]:
    ephem = Ephem.from_table(qt)
    time_max = _max_alt_time(ephem, "obstime", "alt")
    # Return tuple for compatibility with flip_times()
    return time_max, time_max

def get_row_for_time(qt: QTable, t: Time) -> Row:
    ephem = Ephem.from_table(qt)
    return _get_row_for_time(ephem, t, "obstime")

def opt_alt_times(qt: QTable, alt: Angle) -> tuple[Time, Time]:
    ephem = Ephem.from_table(qt)
    return _opt_alt_times(ephem, alt, "obstime", "alt")



def get_times_from_eph(ephemerides: dict) -> dict:
    """
    Compile complete list of times (before, after, start, end, alt_first, alt_last)
    for entire ephemerides dictionary

    Parameters
    ----------
    ephemerides : dict
        Ephemerides dictionary {id: qtable}

    Returns
    -------
    dict
        Times dictionary {id: {name: value}}
    """
    times_list = {}

    for id, qt in ephemerides.items():
        time_before, time_after = flip_times(qt)
        if not time_before:     # No meridian passing
                ##CHECK: better solution than this hack?
                time_before, time_after = max_alt_times(qt)
        time_0 = qt["obstime"][0]
        time_1 = qt["obstime"][-1]
        time_alt0, time_alt1 = opt_alt_times(qt, config.opt_alt * u.degree)
        if not time_alt0:
            time_alt0, time_alt1 = time_0, time_1

        times_list[id] = {}
        times_list[id]["before"]     = time_before
        times_list[id]["after"]      = time_after
        times_list[id]["first"]      = time_0
        times_list[id]["last"]       = time_1
        times_list[id]["alt_first"]  = time_alt0
        times_list[id]["alt_last"]   = time_alt1
        ic(id, qt, times_list[id])

    return times_list



def process_objects(ephemerides: dict, neocp_list: dict, pccp_list: dict, times_list: dict) -> list:
    """
    Main planning for NEOCP observations

    Parameters
    ----------
    ephemerides : dict
        Ephemerides dictionary for all bjects
    neocp_list : dict
        NEOCP list data dictionary for all objects
    pccp_list : dict
        PCCP list data dictionary for all objects
    times_list : dict
        Times from ephemerides dictionary for all objects

    Returns
    -------
    list
        List of planned objects
    """
    ic(ephemerides.keys(), neocp_list.keys(), pccp_list.keys())

    message("-----------------------------------------------------------------------------------------------------------------------")
    message("             Score      MagV #Obs      Arc NotSeen  Time start ephemeris   / end ephemeris                  Max motion")
    message("                                                    Time before            / after meridian              Moon distance")
    message("                                                    Time start exposure    / end exposure")
    message("                                                    # x Exp = total exposure time")
    message("                                                    RA, DEC, Alt, Az")

    csv_row = None
    csv_rows = []
    objects  = []
    prev_time_end_exp = None

    for id, qt in ephemerides.items():
        item    = neocp_list[id]
        time    = times_list[id]
        type    = "PCCP" if id in pccp_list else "NEOCP"
        ic(id, qt, item, time, type)
        score   = item["score"]
        mag     = item["mag"]
        nobs    = item["nobs"]
        arc     = item["arc"]
        notseen = item["notseen"]
        time_before, time_after       = time["before"], time["after"]
        time_first, time_last         = time["first"], time["last"]
        time_alt_first, time_alt_last = time["alt_first"], time["alt_last"]
        max_m = max_motion(qt)

        message("-----------------------------------------------------------------------------------------------------------------------")
        message(f"{id}  {type:5s} {score:3d}  {mag}  {nobs:3d}  {arc:5.2f}  {notseen:4.1f}  {time_first}/{time_last}  {max_m:5.1f}")

        # Get exposure data
        exposure = exposure_calc(max_m, mag)
        if not exposure:                             # Object too fast
            message(f"SKIPPED: object too fast (>{motion_limit():.1f})")
            continue

        n_exp = exposure.number
        exp = exposure.single
        total_exp = exposure.total
        total_time = exposure.total_time
        perc_of_required = exposure.percentage        
        ic(n_exp, exp, total_exp, total_time, perc_of_required)

        ##### Plan exposure start and end time ... #####
        # Make sure to start after previous exposure
        if prev_time_end_exp != None and time_first < prev_time_end_exp:
            time_first = prev_time_end_exp
        ic(prev_time_end_exp, time_first)

        # 1st try: start at time_alt_first or time_first
        time_start_exp = time_alt_first if time_alt_first >= time_first else time_first
        time_end_exp   = time_start_exp + total_time
        ic("1st try - alt first", time_start_exp, time_end_exp)
        # ok if end <= before or start >= after and start >= previous

        # 2nd try: start at time_before - total = before passing meridian
        if time_end_exp > time_before and time_start_exp <= time_after:
            time_start_exp = time_before - total_time
            time_end_exp   = time_before
            ic("2nd try - before meridian", time_start_exp, time_end_exp)
        # ok if start >= first

        # 3rd try: start at time_after = after passing meridian
        if time_start_exp < time_first:
            time_start_exp = time_after
            time_end_exp   = time_after + total_time
            ic("3rd try - after meridian", time_start_exp, time_end_exp)
        # ok if end <= last
        ic(time_start_exp, time_end_exp)

        ##### Skip object for various reasons ... #####
        # Skip, if below threshold for # obs
        if nobs < config.min_n_obs:
            message(f"SKIPPED: only {nobs} obs (< {config.min_n_obs})")
            continue

        # Skip, if not seen for more than threshold days
        max_notseen = config.max_notseen * u.day
        if notseen > max_notseen:
            message(f"SKIPPED: not seen for {notseen:.1f} (> {max_notseen:.1f})")
            continue

        # Skip, if percentage of total exposure time is less than threshold
        if perc_of_required < config.min_perc_required:
            message(f"SKIPPED: only {perc_of_required:.0f}% of required total exposure time (< {config.min_perc_required}%)")
            continue

        # Skip, if arc is less than threshold
        min_arc = config.min_arc * u.day
        if arc < min_arc:
            message(f"SKIPPED: arc {arc:.2f} too small (< {min_arc})")
            continue

        # Skip, if failed to allocate total_time
        if time_end_exp > time_last:
            message(f"SKIPPED: can't allocate exposure time {total_time:.2f} ({time_first} -- {time_last})")
            continue

        # Table row best matching time_start_exp
        row = get_row_for_time(qt, time_start_exp)
        ic(row)
        moon_dist = row["moon_dist"][0] ##!!!
        ic(moon_dist)
        # Skip, if moon distance is too small
        min_moon_dist = config.min_moon_dist * u.degree
        if moon_dist < min_moon_dist:
            message(f"SKIPPED: moon distance {moon_dist:.0f} < {min_moon_dist:.0f}")
            continue

        ##### Good to go! #####
        # Remember end of exposure
        prev_time_end_exp = time_end_exp
        # Append to list of planned objects
        objects.append(id)
        ra, dec = row["ra"][0], row["dec"][0] ##!!!
        alt, az = row["alt"][0], row["az"][0] ##!!!

        message(f"                                                    {time_before}/{time_after}             {moon_dist:3.0f}")
        message(f"                                                    {time_start_exp}/{time_end_exp}")
        total = f"{n_exp} x {exp:2.0f} = {total_exp:3.1f} ({perc_of_required:.0f}%) / total {total_time:3.1f}"
        message(f"                                                    {total}")
        message(f"                                                    RA {ra:.4f}, DEC {dec:.4f}, Alt {alt:.0f}, Az {az:.0f}")

        # CSV output:
        #   start time, end time, 
        #   target=id, observation date (YYYY-MM-DD), time ut (HH:MM), ra, dec, exposure, number, filter (L),
        #   type, mag, nobs, arc, notseen, total
        csv_row = { "target": id,
                    # "observation date": time_start_exp.strftime("%Y-%m-%d"),
                    # "time ut": time_start_exp.strftime("%H:%M"),
                    # "+0000" to enforce UTC
                    "obstime": time_start_exp.strftime("%Y-%m-%d %H:%M:%S+0000"),
                    "ra": float(ra.value),
                    "dec": float(dec.value),
                    "exposure": float(exp.value),
                    "number": n_exp,
                    "filter": "L",
                    "start time": str(time_start_exp),
                    "end time": str(time_end_exp),
                    "type": type,
                    "mag": float(mag.value),
                    "nobs": nobs,
                    "arc": float(arc.value),
                    "notseen": float(notseen.value),
                    "total": total
                    }
        ic(csv_row)
        csv_rows.append(csv_row)

        ##MJ: only 1st object for debugging
        # return

    # Return list of planned objects
    message("-----------------------------------------------------------------------------------------------------------------------")
    message(f"{len(objects)} object(s) planned: {" ".join(objects)}")

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

    # Return planned objects
    return objects



def sort_by_flip_time(ephemerides: dict) -> dict:
    """
    Sort ephemerides dictionary by time before meridian/max altitude time

    Parameters
    ----------
    ephemerides : dict
        Ephemerides dictionary

    Returns
    -------
    dict
        Sorted ephemerides dictionary
    """
    time_dict = {}

    for id, qt in ephemerides.items():
        t, _ = flip_times(qt)
        if not t:
            t, _ = max_alt_times(qt)
        time_dict[id] = t
    
    # Sort dict by time (item[0] = id, item[1] = time)
    time_sorted = { id: time for id, time in sorted(time_dict.items(), key=lambda item: item[1]) }

    # Return table_dict sorted by time
    return { id: ephemerides[id] for id in time_sorted.keys() }



def main():
    prefix = Time.now().strftime("%Y%m%d")

    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = "Parse NEOCP ephemerides",
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("-l", "--location", help="MPC station code")
    arg.add_argument("-U", "--update-neocp", action="store_true", help="update NEOCP data from MPC")
    arg.add_argument("-P", "--plot", action="store_true", help="create altitude and sky plot with objects")
    arg.add_argument("-o", "--output", help="write CSV to OUTPUT file")
    arg.add_argument("-C", "--csv", action="store_true", help="use CSV output format")
    arg.add_argument("-M", "--mag-limit", help="override mag_limit from config")
    arg.add_argument("-p", "--prefix", help=f"prefix for cached MPD data, default {prefix}")

    args = arg.parse_args()

    if args.debug:
        ic.enable()
        ic(sys.version_info, sys.path)
        ic(args)
    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()

    Options.csv    = args.csv
    Options.output = args.output or os.path.join(config.downloads, f"{prefix}-neocp-plan.csv")
    
    if args.mag_limit:
        config.mag_limit = float(args.mag_limit)
    if args.location:
        loc = mpc_station_location(args.location)
        Options.code = args.location
        Options.loc  = loc
        ic(loc, loc.to_geodetic())
    verbose(f"location: {Options.code} {location_to_string(Options.loc)}")

    # Get next midnight
    time = Time(Time.now(), location=Options.loc)
    observer = Observer(location=Options.loc, description=Options.code)
    midnight = observer.midnight(Time.now(), which="next")
    min_time = midnight - 12 * u.hour
    max_time = midnight + 12 * u.hour
    ic(midnight, min_time, max_time)

    if args.prefix:
        prefix = args.prefix
        if args.update_neocp:
            error("don't use --prefix with --update-neocp")

    local_eph   = os.path.join(config.downloads, f"{prefix}-{config.local_eph}")
    local_neocp = os.path.join(config.downloads, f"{prefix}-{config.local_neocp}")
    local_pccp  = os.path.join(config.downloads, f"{prefix}-{config.local_pccp}")
    ic(prefix, local_eph, local_neocp, local_pccp)

    if args.update_neocp:
        verbose(f"download ephemerides from {config.url_neocp_query}")
        mpc_query_neocp_ephemerides(config.url_neocp_query, local_eph, Options.loc, Options.code)
        verbose(f"download NEOCP list from {config.url_neocp_list}")
        mpc_query_neocp_list(config.url_neocp_list, local_neocp)
        verbose(f"download PCCP list from {config.url_pccp_list}")
        mpc_query_neocp_list(config.url_pccp_list, local_pccp)

    try:
    # Parse ephemerides
        verbose(f"processing {local_eph}")
        with open(local_eph, "r") as file:
            content = file.readlines()
            ephemerides_txt = parse_html_ephemerides(content)
            ephemerides = convert_text_ephemerides(ephemerides_txt, min_time, max_time)
            times = get_times_from_eph(ephemerides)

        # Parse lists
        verbose(f"processing {local_neocp}")
        with open(local_neocp, "r") as file:
            content = file.readlines()
            neocp_list = parse_neocp_list(content)

        verbose(f"processing {local_pccp}")
        with open(local_pccp, "r") as file:
            content = file.readlines()
            pccp_list = parse_neocp_list(content)

    except FileNotFoundError as e:
        error(e)

    ephemerides = sort_by_flip_time(ephemerides)

    verbose("processing objects:", " ".join(ephemerides.keys()))
    print_ephemerides(ephemerides)
    objects = process_objects(ephemerides, neocp_list, pccp_list, times)

    # Plot objects and Moon
    if args.plot:
        verbose("altitude and sky plot for objects")
        plot_objects(ephemerides, objects, 
                     os.path.join(config.downloads, f"{prefix}-neocp-plot.png"),
                     Options.loc)



if __name__ == "__main__":
    main()
