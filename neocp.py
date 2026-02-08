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
# Version 0.9 / 2026-02-07
#       Refactoring ready, removed old code
# Version 0.10 / 2026-02-08
#       Unified handling of ephemeris column names

VERSION = "0.10 / 2026-02-08"
AUTHOR  = "Martin Junius"
NAME    = "neocp"

import sys
import os
import argparse

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
from mpcneocp import mpc_query_neocp_ephemerides, mpc_query_neocp_list, parse_html_ephemerides
from mpcneocp import parse_neocp_list, obj_data_from_text_ephemerides, obj_data_add_neocp_list
from neoclasses import EphemData
from neoutils import obj_data_add_times, sort_obj_data, verbose_obj_data
from neoutils import motion_limit, max_alt_time, get_row_for_time, obj_data_csv_output



# Command line options
class Options:
    """
    Command line options
    """
    csv: bool = False           # -C --csv
    output: str = None          # -o --output
    code: str = config.code     # -l --location
    loc: EarthLocation = mpc_station_location(code)



def obs_planner_neocp(obj_data: dict[str, EphemData]) -> None:
    """
    New version of process_objects()
    """
    message("-----------------------------------------------------------------------------------------------------------------------")
    message("             Score      MagV #Obs      Arc NotSeen  Time start ephemeris   / end ephemeris                  Max motion")
    message("                                                    Time before            / after meridian              Moon distance")
    message("                                                    Time start exposure    / end exposure")
    message("                                                    # x Exp = total exposure time")
    message("                                                    RA, DEC, Alt, Az")

    objects  = []
    prev_time_end_exp = None

    for obj, edata in obj_data.items():
        max_m = edata.motion
        eph = edata.ephem
        ic(obj, eph, max_m)
        # NEOCP/PCCP list data
        type    = edata.neocp.type
        score   = edata.neocp.score
        mag     = edata.neocp.mag
        nobs    = edata.neocp.nobs
        arc     = edata.neocp.arc
        notseen = edata.neocp.notseen
        ic(obj, type, score, mag, nobs, arc, notseen)
        # Times
        time_before, time_after       = edata.times.before, edata.times.after
        time_first, time_last         = edata.times.start, edata.times.end
        time_alt_first, time_alt_last = edata.times.alt_start, edata.times.alt_end
        # Previously handled by get_times_from_eph()
        if not time_alt_first:
            time_alt_first, time_alt_last = edata.times.start, edata.times.end
        if not time_before:
            time_before = time_after = edata.times.max_alt
        ic(time_before, time_after, time_first, time_last, time_alt_first, time_alt_last)

        message("-----------------------------------------------------------------------------------------------------------------------")
        message(f"{obj}  {type:5s} {score:3d}  {mag}  {nobs:3d}  {arc:5.2f}  {notseen:4.1f}  {time_first}/{time_last}  {max_m:5.1f}")

        # Get exposure data
        exposure = edata.exposure
        if not exposure:                             # Object too fast
            message(f"SKIPPED: object too fast (>{motion_limit():.1f})")
            continue

        # Skip, if ephemeris is only 1 line
        if len(eph) < 2:
            message(f"SKIPPED: only {len(eph)} line(s) of ephemeris data")
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
        row = get_row_for_time(eph, time_start_exp)
        ic(row)
        moon_dist = row["Moon_dist"]
        ic(moon_dist)
        # Skip, if moon distance is too small
        min_moon_dist = config.min_moon_dist * u.degree
        if moon_dist < min_moon_dist:
            message(f"SKIPPED: Moon distance {moon_dist:.0f} < {min_moon_dist:.0f}")
            continue

        ##### Good to go! #####
        # Remember end of exposure
        prev_time_end_exp = time_end_exp
        edata.times.plan_start = time_start_exp
        edata.times.plan_end = time_end_exp
        # Append to list of planned objects
        objects.append(obj)
        ra, dec = row["RA"], row["DEC"]
        alt, az = row["Alt"], row["Az"]
        edata.ra, edata.dec = ra, dec

        message(f"                                                    {time_before}/{time_after}             {moon_dist:3.0f}")
        message(f"                                                    {time_start_exp}/{time_end_exp}")
        total = f"{n_exp} x {exp:2.0f} = {total_exp:3.1f} ({perc_of_required:.0f}%) / total {total_time:3.1f}"
        message(f"                                                    {total}")
        message(f"                                                    RA {ra:.4f}, DEC {dec:.4f}, Alt {alt:.0f}, Az {az:.0f}")

        ##MJ: only 1st object for debugging
        # return

    # Return list of planned objects
    message("-----------------------------------------------------------------------------------------------------------------------")
    message(f"{len(objects)} object(s) planned: {" ".join(objects)}")



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
            obj_data = obj_data_from_text_ephemerides(ephemerides_txt, min_time, max_time)
            # obj_data = obj_data_add_times(obj_data)
            ##FIXME: using old sort order for testing
            obj_data = obj_data_add_times(obj_data, use_old_sort=True)

        # Parse lists
        verbose(f"processing {local_neocp}")
        with open(local_neocp, "r") as file:
            content = file.readlines()
            neocp_list = parse_neocp_list(content)
            obj_data = obj_data_add_neocp_list(obj_data, neocp_list)

        verbose(f"processing {local_pccp}")
        with open(local_pccp, "r") as file:
            content = file.readlines()
            pccp_list = parse_neocp_list(content)
            obj_data = obj_data_add_neocp_list(obj_data, pccp_list, is_pccp=True)

    except FileNotFoundError as e:
        error(e)

    obj_data = sort_obj_data(obj_data)
    verbose("planning objects:", " ".join(obj_data.keys()))
    verbose_obj_data(obj_data)    
    obs_planner_neocp(obj_data)

    # Output CSV plan
    if Options.csv:
        obj_data_csv_output(obj_data, Options.output)

    # Plot objects and Moon
    if args.plot:
        verbose("altitude and sky plot for objects")
        plot_objects(obj_data,
                     os.path.join(config.downloads, f"{prefix}-neocp-plot.png"),
                     Options.loc)



if __name__ == "__main__":
    main()
