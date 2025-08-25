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
# Version 0.0 / 2025-08-25
#       First attempt at parsing MPC NEOCP ephemerids

VERSION = "0.0 / 2025-08-25"
AUTHOR  = "Martin Junius"
NAME    = "neocp"

import sys
import argparse
import csv
import re
from datetime import datetime, timezone
from zoneinfo import ZoneInfo
# Required on Windows
import tzdata

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.coordinates import SkyCoord, AltAz
import astropy.units as u
from astropy.time        import Time
import numpy as np
from astropy.table import QTable

# Astroplan
# import matplotlib.pyplot as plt
# from astroplan import FixedTarget, Observer
# from astroplan.plots import plot_airmass, plot_altitude, plot_sky, plot_finder_image, plot_parallactic

# Local modules
from verbose import verbose, warning, error
from astroutils import get_location, get_coord, coord_to_string, location_to_string



# Command line options
class Options:
    pass


def convert_all_to_qtable(eph_dict: dict) -> dict:
    qtable_dict = {}
    for id, eph in eph_dict.items():
        qtable_dict[id] = eph_to_qtable(id, eph)
    return qtable_dict


def eph_to_qtable(id: str, eph: list) -> QTable:
    ic(id)
    qt = QTable()
    qt.meta["comments"] = [ f"NEOCP temporary designation: {id}" ]
    # Initialize columns
    qt["datetime"]  = Time("2000-01-01 00:00")
    qt["coord"]     = SkyCoord(0, 0, unit=(u.hour, u.degree))
    qt["mag"]       = 0 * u.mag
    qt["motion"]    = 0 * u.arcsec / u.min
    qt["altaz"]     = AltAz(0 * u.degree, 0 * u.degree)
    qt["moon_dist"] = 0 * u.degree
    qt["moon_alt"]  = 0 * u.degree

    for line in eph:
        # Date       UT      R.A. (J2000) Decl.  Elong.  V        Motion     Object     Sun         Moon        Uncertainty
        #             h m                                      "/min   P.A.  Azi. Alt.  Alt.  Phase Dist. Alt.        
        # 2025 08 26 0400   02 48 19.2 -26 44 25 114.9  19.5    1.79  239.8  064  +81   -17    0.09  133  -38
        # ^0         ^11    ^18        ^29              ^47   ^53            ^67  ^72                ^93  ^98
        time      = Time( line[0:10].replace(" ", "-") + " " + line[11:13]+":"+line[13:15] )
        coord     = SkyCoord( line[18:28], line[29:38], unit=(u.hour, u.deg) )
        mag       = float(line[47:51]) * u.mag
        motion    = float(line[53:59].strip()) * u.arcsec / u.min
        altaz     = AltAz(float(line[67:70]) * u.degree, float(line[72:75]) * u.degree)
        moon_dist = float(line[93:96]) * u.degree
        moon_alt  = float(line[98:101]) * u.degree
        ic(time, coord, mag, motion, altaz, moon_dist, moon_alt)
        qt.add_row([ time, coord, mag, motion, altaz, moon_dist, moon_alt ])

    ic(qt)
    # qt.write(sys.stdout, format="ascii")




def parse_neocp_eph(content: list) -> dict:
    neocp_id = None
    neocp_eph = {}

    content_iter = iter(content)
    while (line := next(content_iter, None)) != None:
        line = line.strip()
        # ic(line)

        # Observatory code
        m = re.match(r"observatory code ([0-9A-Z]{3}).", line)
        if m:
            ic(line)
            code = m.group(1)
            ic(code)
    
        # <hr> marks start of NEOCP ephemerids
        m = re.search(r"<p></p><hr><p>", line)
        if m:
            ic(line)

            line = next(content_iter, None).strip()
            m = re.match(r"</p><p><b>(.+)</b>", line)
            if m:
                ic(line)
                neocp_id = m.group(1)
                ic(neocp_id)

                # Read lines until </pre>, ephemerids data starts with date
                neocp_eph[neocp_id] = []
                while (line := next(content_iter, None).strip()) != None:
                    m = re.match(r"</pre>", line)
                    if m:
                        ic(line)
                        break
                    m = re.match(r"</p><pre>Date", line)
                    if m:
                        line = line[9:]
                        ic(line)
                        line = next(content_iter, None).rstrip()
                        ic(line)
                    m = re.match(r"\d\d\d\d \d\d \d\d \d\d", line)
                    if m:
                        line = line[0:100]
                        ic(line)
                        neocp_eph[neocp_id].append(line)

                # Early return, just the 1st entry in the list for debugging
                # return neocp_eph

                # m = re.match(r"", line)
                # if m:
                #     ic(line)

    return neocp_eph



def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = "Parse NEOCP ephemerids",
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("resultpage", help="NEOCP query results page (.htm)")

    args = arg.parse_args()

    if args.debug:
        ic.enable()
        ic(sys.version_info)
        ic(args)
    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()

    with open(args.resultpage, "r") as file:
        content = file.readlines()
        eph_dict = parse_neocp_eph(content)
        table_dict = convert_all_to_qtable(eph_dict)



if __name__ == "__main__":
    main()
