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
# Version 0.3 / 2026-03-02
#       Output to neo_obs_data_dir, added log file
# Version 0.4 / 2026-03-15
#       Added -m --min-alt option
# Version 0.5 / 2026-05-18
#       Planner output aligned with neocp output
# Version 1.0 / 2026-05-25
#       Joint planner for NEOCP, SBWOBS and other objects, added all code and
#       options from neocp.py
# Version 1.1 / 2026-06-13
#       Removed -U --update-neocp option, no local files cached anymore
# Version 2.0 / 2026-06-16
#       Refactored version in new subdirectory neoop/
# Version 2.1 / 2026-06-16
#       New option --verbose-ephem for ephemerides output, fixed mag limits

VERSION     = "2.1 / 2026-06-20"
AUTHOR      = "Martin Junius"
NAME        = "neo-obs-planner"
DESCRIPTION = "NEOCP/NEO observation planner"

import sys
import argparse

from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.time import Time
import astropy.units as u
from astroquery.mpc import MPC

# Local modules
from utils.verbose import verbose, warning, error, message
from astro.utils   import fmt_time
from neo.config    import config
from neo.classes   import EphemData, EphemDataList, LocalCircumstances
from neo.utils     import edata_list_add_times, edata_list_csv_output
from neo.exposure  import motion_limit, edata_add_exposure
from mpc.ephemdata import edata_add_ephem_mpc
from neo.plot      import edata_list_plot
from jpl.sbwobs    import sbwobs_get_edata_list
from mpc.neocp     import neocp_get_edata_list
import neo.files

DEFAULT_LOCATION = config.code



def obs_planner_1(edata_list: EphemDataList, local: LocalCircumstances) -> None:
    # Start planner at naut. dusk / start time from options
    next_start_time = local.naut_dusk
    objects = list()
    skipped = list()

    message("-------------------------------------------------------------------------------------------------------------------")
    message("              Score       Mag #Obs      Arc NotSeen  Time start ephemeris/ end ephemeris                 Max motion")
    message("       /Uncertainty                                  Time before         / after meridian             Moon distance")
    message("                                                     Time start exposure / end exposure")
    message("                                                     # x Exp = total exposure time")
    message("                                                     RA, DEC, Alt, Az")

    edata: EphemData
    for edata in edata_list:
        obj = edata.obj
        eph = edata.ephem

        exp_start = None
        exp_end = None

        etimes = edata.times
        start = etimes.start
        end = etimes.end
        alt_start = etimes.alt_start
        alt_end = etimes.alt_end
        before = etimes.before
        after = etimes.after

        message("-------------------------------------------------------------------------------------------------------------------")
        type    = edata.type
        if edata.neocp:
            score   = edata.neocp.score
            mag     = edata.neocp.mag
            nobs    = edata.neocp.nobs
            arc     = edata.neocp.arc
            notseen = edata.neocp.notseen
        elif edata.dlx:
            score   = edata.dlx.uncertainty
            mag     = edata.mag
            nobs    = None
            arc     = None
            last    = edata.dlx.last_obs
            notseen = (Time.now() - last).to(u.day)
        else:
            score   = None
            mag     = edata.mag
            nobs    = None
            arc     = None
            notseen = None

        s_score = f"{score:3d}" if score != None else "   "
        s_nobs = f"{nobs:3d}" if nobs != None else "   "
        s_arc = f"{arc:5.2f}" if arc != None else "       "
        s_notseen = f"{notseen:4.1f}" if notseen != None else "      "
        message(f"{obj:9s} {type:5s} {s_score}  {mag}  {s_nobs}  {s_arc}  {s_notseen}  {fmt_time(start)} / {fmt_time(end)}   {edata.motion:5.1f}")

        # check overlap with previous object
        if next_start_time > start:
            start = next_start_time
        if alt_start and next_start_time > alt_start:
            alt_start = next_start_time
        if alt_end and next_start_time > alt_end:
            alt_start = None
            alt_end = None
        if before and next_start_time > before:
            before = None
        if after and next_start_time > after:
            after = next_start_time
        if next_start_time > end:
            end = None

        # check valid exposure data, skip if None = object is too fast
        ##FIXME: meanwhile this should be called earlier
        if edata.exposure:
            total_time = edata.exposure.total_time
        else:
            message(f"SKIPPED: object too fast (>{motion_limit():.1f})")
            skipped.append(obj)
            continue

        # Skip, if ephemeris is only 1 line
        if len(eph) < 2:
            message(f"SKIPPED: only {len(eph)} line(s) of ephemeris data")
            skipped.append(obj)
            continue

        ic(next_start_time.iso)
        ic(obj, start, end, total_time, alt_start, alt_end, before, after)

        # try to fit in alt_start ... alt_end interval
        if alt_start:
            # before meridian
            if before and alt_start + total_time <= before:
                exp_start = alt_start
                exp_end = alt_start + total_time
            # after meridian
            elif after and end and after + total_time <= end:
                exp_start = after
                exp_end = after + total_time
            # no meridian passing
            elif end and alt_start + total_time <= end:
                exp_start = alt_start
                exp_end = alt_start + total_time
        
        # no slot found, try to fit in start ... end interval
        if not exp_start:
            # before meridian
            if before and start + total_time <= before:
                exp_start = start
                exp_end = start + total_time
            # after meridian
            elif after and end and after + total_time <= end:
                exp_start = after
                exp_end = after + total_time
            # no meridian passing
            elif end and start + total_time <= end:
                exp_start = start
                exp_end = start + total_time

        ic(exp_start, exp_end)

        # Skip, if failed to allocate total_time
        if not exp_start:
            message(f"SKIPPED: can't allocate exposure time {total_time:.2f} ({start} -- {end})")
            skipped.append(obj)
            continue

        # Ephemeris row best matching start time
        row = eph.get_row_for_time(exp_start)
        ic(row)

        if not edata.force:
            ##### Skip object for various reasons ... #####
            # Skip, if percentage of total exposure time is less than threshold
            if edata.exposure.percentage < config.min_perc_required:
                message(f"SKIPPED: only {edata.exposure.percentage:.0f}% of required total exposure time (< {config.min_perc_required}%)")
                skipped.append(obj)
                continue

            # Skip, if moon distance is too small
            moon_dist = row["Moon_dist"]
            min_moon_dist = config.min_moon_dist * u.degree
            if moon_dist < min_moon_dist:
                message(f"SKIPPED: moon distance {moon_dist:.0f} < {min_moon_dist:.0f}")
                skipped.append(obj)
                continue

            # Skip, if below threshold for # obs
            if not nobs is None and nobs < config.min_n_obs:
                message(f"SKIPPED: only {nobs} obs (< {config.min_n_obs})")
                skipped.append(obj)
                continue

            # Skip, if not seen for more than threshold days
            max_notseen = config.max_notseen * u.day
            if edata.neocp and not notseen is None and notseen > max_notseen:
                message(f"SKIPPED: not seen for {notseen:.1f} (> {max_notseen:.1f})")
                skipped.append(obj)
                continue

            # Skip, if arc is less than threshold
            min_arc = config.min_arc * u.day
            if not arc is None and arc < min_arc:
                message(f"SKIPPED: arc {arc:.2f} too small (< {min_arc})")
                skipped.append(obj)
                continue
        # /if

        # Good to go!
        # Remember end of exposure
        next_start_time = exp_end
        etimes.plan_start = exp_start
        etimes.plan_end = exp_end

        ra, dec = row["RA"].to(u.hourangle), row["DEC"]
        alt, az = row["Alt"], row["Az"]
        edata.ra, edata.dec = ra, dec
        objects.append(obj)

        message(f"{'':53s}{fmt_time(before)} / {fmt_time(after)}              {moon_dist:3.0f}")
        message(f"{'':53s}{fmt_time(exp_start)} / {fmt_time(exp_end)}")
        message(f"{'':53s}{edata.exposure}")
        message(f"{'':53s}RA {ra:.4f}, DEC {dec:.4f}, Alt {alt:.0f}, Az {az:.0f}")

    # end for
    message("-------------------------------------------------------------------------------------------------------------------")
    message(f"{len(objects)} object(s) planned: {", ".join(objects)}")
    message(f"{len(skipped)} object(s) skipped: {", ".join(skipped)}")



def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = DESCRIPTION,
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("--verbose-ephem", action="store_true", help="verbose ephemerides")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("-l", "--location", help=f"coordinates, named location or MPC station code, default {DEFAULT_LOCATION}")
    arg.add_argument("-f", "--file", help="read list of objects from file")
    arg.add_argument("-s", "--start", help="start time (UTC) (default naut. dusk)")
    arg.add_argument("-e", "--end", help="end time (UTC) (default naut. dawn)")
    arg.add_argument("-o", "--output", help="write CSV to OUTPUT file")
    arg.add_argument("-C", "--csv", action="store_true", help="use CSV output format")
    arg.add_argument("-P", "--plot", action="store_true", help="create altitude and sky plot with objects")

    arg.add_argument("--clear", action="store_true", help="clear MPC cache")

    arg.add_argument("-M", "--mag-limit", type=float, help="override *mag_limits from config (see below)")
    arg.add_argument("--neocp-mag-limit", type=float, help=f"override neocp_mag_limit from config ({config.neocp_mag_limit})")
    arg.add_argument("--sbwobs-mag-limit", type=float, help=f"override sbwobs_mag_limit from config ({config.sbwobs_mag_limit})")

    arg.add_argument("-m", "--min-alt", type=float, help="override min_alt/elev_min from config")
    arg.add_argument("--neocp", action="store_true", help="observable NEOCP objects")
    arg.add_argument("--sbwobs", action="store_true", help="observable objects from JPL WOBS service")
    arg.add_argument("--asteroids", action="store_true", help=f"sbwobs: get asteroids default={config.sb_kind}")
    arg.add_argument("--neo", action="store_true", help=f"sbwobs: get NEOs default={config.sb_group}")
    arg.add_argument("--pha", action="store_true", help=f"sbwobs: get PHAs")
    arg.add_argument("--comets", action="store_true", help=f"sbwobs: get comets (overrides asteroids options)")
 
    arg.add_argument("-p", "--prefix", help=f"prefix for planner data, default {neo.files.prefix}")
    arg.add_argument("--force", help=f"skip checks for FORCE objects, include in observation plan")

    arg.add_argument("object", nargs="*", help="object name")

    args = arg.parse_args()

    if args.debug:
        ic.enable()
        ic(sys.version_info, sys.path, args)
    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()

    # override defaults from config
    if args.mag_limit:
        config.neocp_mag_limit = args.mag_limit
        config.sbwobs_mag_limit = args.mag_limit
    if args.neocp_mag_limit:
        config.neocp_mag_limit = args.neocp_mag_limit
    if args.sbwobs_mag_limit:
        config.sbwobs_mag_limit = args.sbwobs_mag_limit
    if not args.min_alt is None:
        config.min_alt = args.min_alt
    if args.asteroids:
        config.sb_kind = "a"
    if args.neo:
        config.sb_group = "neo"
    if args.pha:
        config.sb_group = "pha"
    if args.comets:
        config.sb_kind = "c"
        config.sb_group = None

    # Observer location and local circumstances
    local = LocalCircumstances.from_location(args.location if args.location else DEFAULT_LOCATION)

    # Update DEC limits in config
    min_dec, max_dec = local.get_dec_limits(config.min_alt*u.deg)
    config.min_dec = int(min_dec.degree)
    config.max_dec = int(max_dec.degree)

    # Prefix for cached NEOCP data
    if args.prefix:
        neo.files.set_prefix(args.prefix)

    # Forced objects
    forced_objs = []
    if args.force:
        forced_objs = args.force.split(",")


    edata_list = EphemDataList()

    # Objects from NEOCP
    if args.neocp:
        edata_list.extend( neocp_get_edata_list(local) )

    # Objects from SBWOBS
    if args.sbwobs:
        edata_list.extend( sbwobs_get_edata_list(local) )

    # Objects from file / command line
    type = "-"
    objects = list()
    if args.file:
        with open(args.file, "r") as file:
            for line in file:
                objects.append(line.strip())
    if args.object:
        objects.extend(args.object)
    if objects:
        edata_list = edata_list.append_objects(objects)

    if not edata_list:
        error("no objects from file or command line")
    
    ic(edata_list)

    # Clear astroquery cache
    if args.clear:
        MPC.clear_cache()    

    # Override start/end time
    if args.start:
        local.naut_dusk = Time(args.start)
    if args.end:
        local.naut_dawn = Time(args.end)
    verbose.print_lines(local)

    # Get ephemerides from MPC (JPL doesn't provide moon phase/distance)
    edata_list.process(edata_add_ephem_mpc, local)

    # Get exposure data from mag and motion
    edata_list.process(edata_add_exposure, local)
    ic(edata_list)

    # Process only objects with ephemeris and exposure data
    edata_list = EphemDataList([ edata for edata in edata_list if edata.ephem and edata.exposure ])

    # Process objects
    edata_list_add_times(edata_list)
    verbose(f"original object sequence: {edata_list.objects_str()}")
    # for edata in edata_list:
    #     verbose.print_lines(edata.ephem["Targetname", "Obstime", "RA", "DEC", "Mag", 
    #                                    "Motion", "PA", "Az", "Alt", "Moon_dist", "Moon_alt"])
    #     verbose(edata.exposure)
    #     verbose(edata.times)
    verbose.print_lines(edata_list)

    edata_list.sort_by_time()
    verbose(f"sorted object sequence: {edata_list.objects_str()}")

    verbose("forced objects:", ", ".join(forced_objs))
    if args.verbose_ephem:
        ##FIXME: use process() to implement
        edata: EphemData
        for edata in edata_list:
            if edata.obj in forced_objs:
                edata.force = True
            verbose("")
            verbose.print_lines2(edata.ephem["Targetname", "Obstime", "RA", "DEC", "Mag", 
                                            "Motion", "PA", "Az", "Alt", "Moon_dist", "Moon_alt"])
        # edata_list.verbose_ephem()

    verbose("")
    log_file = neo.files.path("obs-planner-1.log")
    with verbose.logfile(log_file):
        # NEOCP planner
        verbose(f"obs-planner-1 {fmt_time(neo.files.now)} {neo.files.now.scale.upper()}")

        # Run obs planner
        obs_planner_1(edata_list, local)
        if args.csv:
            ##FIXME: output file name depending on mode
            edata_list_csv_output(edata_list, args.output or neo.files.path("neo-obs-plan.csv"))

        # Plot objects and Moon
        if args.plot:
            plot_file = neo.files.path("neo-obs-plot.png")
            verbose(f"altitude and sky plot for objects: {plot_file}")
            edata_list_plot(edata_list, plot_file, local.loc)



if __name__ == "__main__":
    main()
