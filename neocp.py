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
# Version 0.11 / 2026-02-24
#       Added log file output
# Version 0.12 / 2026-03-05
#       Added -f --force option to include objects in observation plan even if they fail checks
# Version 0.13 / 2026-05-18
#       Refactored for EphemDataList

VERSION = "0.13 / 2026-05-18"
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
import astropy.units as u

# Local modules
from verbose    import verbose, warning, error, message
from neoconfig  import config
from neoplot    import edata_list_plot
from mpcneocp   import neocp_get_edata_list
from neoclasses import EphemData, EphemDataList, LocalCircumstances
from neoutils   import motion_limit, get_row_for_time, fmt_time
from neoephem   import get_local_circumstances, edata_add_exposure
from neoutils   import edata_list_add_times, get_row_for_time, motion_limit, fmt_time, edata_list_csv_output
import neofiles

DEFAULT_LOCATION = config.code



def obs_planner_neocp(edata_list: EphemDataList) -> None:
    """
    New version of process_objects()
    """
    message("-----------------------------------------------------------------------------------------------------------------")
    message("             Score      MagV #Obs      Arc NotSeen  Time start ephemeris/ end ephemeris                Max motion")
    message("                                                    Time before         / after meridian            Moon distance")
    message("                                                    Time start exposure / end exposure")
    message("                                                    # x Exp = total exposure time")
    message("                                                    RA, DEC, Alt, Az")

    objects  = []
    prev_time_end_exp = None

    edata: EphemData
    for edata in edata_list:
        obj = edata.obj

        max_m = edata.motion
        eph = edata.ephem
        ic(obj, eph, max_m)
        # NEOCP/PCCP list data
        type    = edata.type
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

        message("-----------------------------------------------------------------------------------------------------------------")
        message(f"{obj:7s}  {type:5s} {score:3d}  {mag}  {nobs:3d}  {arc:5.2f}  {notseen:4.1f}  {fmt_time(time_first)} / {fmt_time(time_last)}  {max_m:5.1f}")

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

        # Table row best matching time_start_exp
        row = get_row_for_time(eph, time_start_exp)
        ic(row)
        moon_dist = row["Moon_dist"]
        ic(moon_dist)

        if not edata.force:
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

            # Skip, if moon distance is too small
            min_moon_dist = config.min_moon_dist * u.degree
            if moon_dist < min_moon_dist:
                message(f"SKIPPED: Moon distance {moon_dist:.0f} < {min_moon_dist:.0f}")
                continue
        # /if

        # Skip, if failed to allocate total_time
        if time_end_exp > time_last:
            message(f"SKIPPED: can't allocate exposure time {total_time:.2f} ({time_first} -- {time_last})")
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

        message(f"                                                    {fmt_time(time_before)} / {fmt_time(time_after)}             {moon_dist:3.0f}")
        message(f"                                                    {fmt_time(time_start_exp)} / {fmt_time(time_end_exp)}")
        total = f"{n_exp} x {exp:2.0f} = {total_exp:3.1f} ({perc_of_required:.0f}%) / total {total_time:3.1f}"
        message(f"                                                    {total}")
        message(f"                                                    RA {ra:.4f}, DEC {dec:.4f}, Alt {alt:.0f}, Az {az:.0f}")

        ##MJ: only 1st object for debugging
        # return

    # Return list of planned objects
    message("-----------------------------------------------------------------------------------------------------------------")
    message(f"{len(objects)} object(s) planned: {" ".join(objects)}")



def main():
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
    arg.add_argument("-p", "--prefix", help=f"prefix for cached MPD data, default {neofiles.prefix}")
    arg.add_argument("-f", "--force", help=f"skip checks for FORCE objects, include in observation plan")

    args = arg.parse_args()

    if args.debug:
        ic.enable()
        ic(sys.version_info, sys.path)
        ic(args)
    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()

    if args.mag_limit:
        config.mag_limit = float(args.mag_limit)

    # location and local circumstances
    local = get_local_circumstances(args.location if args.location else DEFAULT_LOCATION)
    verbose.print_lines(local)

    if args.prefix:
        if args.update_neocp:
            error("don't use --prefix with --update-neocp")
        neofiles.set_prefix(args.prefix)

    forced_objs = []
    if args.force:
        forced_objs = args.force.split(",")

    local_eph   = neofiles.path(config.local_eph)
    local_neocp = neofiles.path(config.local_neocp)
    local_pccp  = neofiles.path(config.local_pccp)
    ic(neofiles.prefix, local_eph, local_neocp, local_pccp)



    edata_list = neocp_get_edata_list(args.update_neocp, local)

    # Get exposure data from mag and motion
    edata_list.process(edata_add_exposure, local)

    ##FIXME: using old sort order for testing
    edata_list_add_times(edata_list, use_old_sort=True)
    edata_list.sort_by_time()
    ic(edata_list)

    verbose("planning objects:", " ".join(edata_list.objects()))
    verbose("forced objects:", " ".join(forced_objs))
    ##FIXME: use process() to implement
    # for obj in forced_objs:
    #     if obj in obj_data:
    #         obj_data[obj].force = True
    edata_list.verbose_ephem()

    log_file = neofiles.path("obs-planner-neocp.log")
    with verbose.logfile(log_file):
        # NEOCP planner
        verbose(f"obs-planner-neocp {fmt_time(neofiles.now)} {neofiles.now.scale.upper()}")
        obs_planner_neocp(edata_list)

        # Output CSV plan
        if args.csv:
            ##FIXME: output file name depending on mode
            edata_list_csv_output(edata_list, args.output or neofiles.path("neo-obs-plan.csv"))

        # Plot objects and Moon
        if args.plot:
            plot_file = neofiles.path("neocp-plot.png")
            verbose(f"altitude and sky plot for objects: {plot_file}")
            edata_list_plot(edata_list, plot_file, local.loc)



if __name__ == "__main__":
    main()
