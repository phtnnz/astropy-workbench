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
#       First attempt at parsing MPC NEOCP ephemerides
# Version 0.1 / 2025-09-04
#       Retrieval from MPC, altitude and sky plots

VERSION = "0.1 / 2025-09-04"
AUTHOR  = "Martin Junius"
NAME    = "neocp"

import sys
import argparse
import csv
import re
import requests
from datetime import datetime, timezone
from zoneinfo import ZoneInfo
# Required on Windows
import tzdata

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.coordinates import SkyCoord, AltAz, Angle, EarthLocation, FK5
import astropy.units as u
from astropy.units import Quantity, Magnitude
from astropy.time        import Time
import numpy as np
from astropy.table import QTable

# Astroplan
import matplotlib.pyplot as plt
from astroplan import FixedTarget, Observer
from astroplan.plots import plot_altitude, plot_sky

# Local modules
from verbose import verbose, warning, error
from astroutils import get_location, get_coord, coord_to_string, location_to_string


# Requests timeout
TIMEOUT = 60

# Exposure times / s
EXP_TIMES = [ 5, 10, 15, 20, 30, 45, 60 ]


# MPC pages:
# NEOCP form
#       https://minorplanetcenter.net/iau/NEO/toconfirm_tabular.html
# NEOCP list
#       https://minorplanetcenter.net/iau/NEO/neocp.txt
# NEOCP query ephemerides
#       https://cgi.minorplanetcenter.net/cgi-bin/confirmeph2.cgi
#
# Example request for M49
# W=a&mb=-30&mf=20.5&dl=-90&du=%2B40&nl=75&nu=100&sort=d&Parallax=1&obscode=M49&long=&lat=&alt=&int=1&start=0&raty=a&mot=m&dmot=p&out=f&sun=n&oalt=26

URL_NEOCP_PAGE  = "https://minorplanetcenter.net/iau/NEO/toconfirm_tabular.html"
URL_NEOCP_QUERY = "https://cgi.minorplanetcenter.net/cgi-bin/confirmeph2.cgi"
URL_NEOCP_LIST  = "https://minorplanetcenter.net/iau/NEO/neocp.txt"
LOCAL_QUERY = "downloads/NEOCP-ephemerides.html"
LOCAL_LIST  = "downloads/NEOCP-list.txt"



# Command line options
class Options:
    mag_limit = 20.5        # -L --mag-limit
    arcsec_tolerance = 3    # Max trail tolerance in arcsec
    resolution = 1.33       # arcsec / pixel resolution of camera/telescope
    code = "M49"            # -l --location
    loc = get_location(code)



def neocp_query_ephemeris(filename: str) -> None:
    ##FIXME: get from config files
    url = URL_NEOCP_QUERY
    data = { 
        "W":	"a",                    # all objects with ...
        "mb":	"-30",                  # max brightness (V)
        "mf":	str(Options.mag_limit), # min brightness (V)
        "dl":	"-90",      ##   # min DEC
        "du":	"+40",      ##   # max DEC
        "nl":	"75",       ##   # min NEO score
        "nu":	"100",           # max NEO score
        "sort":	"d",
        "Parallax":	"1",                # 0=geocentric, 1=code, 2=Lon/Lat/Alt
        "obscode":	Options.code,       # observatory code
        "long":	"",
        "lat":	"",
        "alt":	"",
        "int":	"1",             # interval 0=1h, 1=30m, 2=10m, 3=1m
        "start":"1",        ##   # start now + X hours
        "raty":	"a",             # format: h=trunc. sexagesimal, a=full, d=decimal
        "mot":	"m",             # motion: s="/sec, m="/min, h="/hr, d=degree/day
        "dmot":	"p",             # motion: p=total, r=separate RA/DEC coord. mo tion, s=separate RA/DEC sky motion
        "out":	"f",             # f=full output, b=brief output
        "sun":	"n",             # suppress output: x=never, s=sunset/rise, c=civil, n=nautical, a=astronomical
        "oalt":	"26"        ##   # suppress below min altitude
    }

    ic(url, data)
    response = requests.get(url, params=data, timeout=TIMEOUT)
    ic(response.status_code)
    if response.status_code != 200:
        error(f"query to {url} failed")

    with open(filename, mode="w", encoding=response.encoding) as file:
        file.write(response.text)



def neocp_query_list(filename: str) -> None:
    ##FIXME: get from config files
    url = URL_NEOCP_LIST

    ic(url)
    response = requests.get(url, timeout=TIMEOUT)
    ic(response.status_code)
    if response.status_code != 200:
        error(f"query to {url} failed")

    with open(filename, mode="w", encoding=response.encoding) as file:
        file.write(response.text)



def qtable_to_altaz(id: str, qt: QTable) -> AltAz:
    altaz = AltAz(alt=qt["alt"], az=qt["az"], obstime=qt["obstime"], location=Options.loc)
    # Quick hack to get a proper label from plot_altitude
    altaz.name = id
    ic(altaz)
    return altaz



def plot_alt_objects(table_dict: dict) -> None:
    # Get next midnight
    time = Time(Time.now(), location=Options.loc)
    observer = Observer(location=Options.loc, description=Options.code)
    midnight = observer.midnight(Time.now(), which="next")
    ic(midnight)

    # Intervals around midnight
    time_interval = midnight + np.linspace(-8, 8, 160)*u.hour
    moon          = observer.moon_altaz(time_interval)
    # Quick hack to get a proper label from plot_altitude
    moon.name     = "Moon"

    # Plot all NEOCP objects in dict
    for id, qt in table_dict.items():
        altaz = qtable_to_altaz(id, qt)
        plot_altitude(altaz, observer, altaz.obstime, style_kwargs=dict(fmt="o"))

    plot_altitude(moon, observer, moon.obstime, brightness_shading=True, style_kwargs=dict(fmt="y--"))
    plt.legend(bbox_to_anchor=(1.0, 1.02))

    plt.savefig("tmp/plot.png", bbox_inches="tight")
    plt.close()



def plot_sky_objects(table_dict: dict) -> None:
    # Get next midnight
    time = Time(Time.now(), location=Options.loc)
    observer = Observer(location=Options.loc, description=Options.code)
    midnight = observer.midnight(Time.now(), which="next")
    ic(midnight)

    # Intervals around midnight
    time_interval = midnight + np.linspace(-5, 5, 11)*u.hour
    moon          = observer.moon_altaz(time_interval)
    # Quick hack to get a proper label from plot_altitude
    moon.name     = "Moon"

    # Plot all NEOCP objects in dict
    for id, qt in table_dict.items():
        altaz = qtable_to_altaz(id, qt)
        plot_sky(altaz, observer, altaz.obstime)

    plot_sky(moon, observer, moon.obstime, style_kwargs=dict(color="y", marker="x"))
    plt.legend(bbox_to_anchor=(1.48, 1.11))

    plt.savefig("tmp/plot.png", bbox_inches="tight")
    plt.close()



def max_motion(qt: QTable) -> Quantity:
    max_motion = -1 * u.arcsec / u.min
    for row in qt:
        if row["motion"] > max_motion:
            max_motion = row["motion"]
    return max_motion

def exp_time_from_motion(motion: Quantity) -> Quantity:
    exp = Options.arcsec_tolerance * u.arcsec / motion
    exp_max = EXP_TIMES[0]
    for exp1 in EXP_TIMES:
        if exp1 > exp.to(u.s).value:
            break
        exp_max = exp1
    return exp_max * u.s



def process_objects(table_dict: dict) -> None:
    for id, qt in table_dict.items():
        time0 = qt["obstime"][0]
        alt0  = qt["alt"][0]
        max_m = max_motion(qt)
        exp   = exp_time_from_motion(max_m)
        ic(id, time0, alt0, max_m, exp)



def sort_by_alt_max_time(table_dict: dict) -> dict:
    time_dict = {}

    for id, qt in table_dict.items():
        max_alt = -1 * u.degree
        max_time = None
        for row in qt:
            if row["alt"] > max_alt:
                max_alt = row["alt"]
                max_time = row["obstime"]
        ic(max_alt, max_time)
        time_dict[id] = max_time
    
    # Sort dict by time (item[0] = id, item[1] = time)
    time_sorted = { id: time for id, time in sorted(time_dict.items(), key=lambda item: item[1]) }

    # Return table_dict sorted by time
    return { id: table_dict[id] for id in time_sorted.keys() }



def print_table_dict(table_dict: dict) -> None:
    for id, qt in table_dict.items():
        print("===================================================================================================================")
        print(f"NEOCP {id} ephemerides")
        print(qt)



def convert_all_to_qtable(eph_dict: dict) -> dict:
    qtable_dict = {}
    for id, eph in eph_dict.items():
        qt = eph_to_qtable(id, eph)
        if len(qt["mag"]) == 0:
            verbose(f"skipping NEOCP {id=} (empty)")
            continue
        mag = qt["mag"][0]
        if mag > Options.mag_limit * u.mag:
            verbose(f"skipping NEOCP {id=} {mag=}")
            continue
        qtable_dict[id] = qt
    return qtable_dict


def eph_to_qtable(id: str, eph: list) -> QTable:
    ic(id)
    qt = QTable()
    qt.meta["comments"] = [ f"NEOCP temporary designation: {id}" ]
    # Initialize empty columns with proper Quantity type, same as .add_column()
    # A bit ugly, but I found no other to handle this
    qt["obstime"]   = Time("2000-01-01 00:00")
    qt["ra"]        = 0 * u.hourangle
    qt["dec"]       = 0 * u.degree
    qt["mag"]       = 0 * u.mag
    qt["motion"]    = 0 * u.arcsec / u.min
    qt["pa"]        = 0 * u.degree
    qt["alt"]       = 0 * u.degree
    qt["az"]        = 0 * u.degree
    qt["moon_dist"] = 0 * u.degree
    qt["moon_alt"]  = 0 * u.degree

    for line in eph:
        # Date       UT      R.A. (J2000) Decl.  Elong.  V        Motion     Object     Sun         Moon        Uncertainty
        #             h m                                      "/min   P.A.  Azi. Alt.  Alt.  Phase Dist. Alt.        
        # 2025 08 26 0400   02 48 19.2 -26 44 25 114.9  19.5    1.79  239.8  064  +81   -17    0.09  133  -38
        # ^0         ^11    ^18        ^29              ^46   ^52     ^60    ^67  ^72                ^91  ^96
        time      = Time(line[0:10].replace(" ", "-") + " " + line[11:13]+":"+line[13:15])
        ra        = Angle(line[18:28], unit=u.hourangle)
        dec       = Angle(line[29:38], unit=u.degree)
        mag       = Magnitude(line[46:50], unit=u.mag)
        motion    = Quantity(line[52:58].strip(), unit=u.arcsec / u.min)
        pa        = Angle(line[60:65], unit=u.degree)
        # NEOCP ephemerides counts az from S=0, Astropy N=0
        az        = Angle(line[67:70], unit=u.degree) - 180*u.degree
        az.wrap_at(360 * u.degree, inplace=True)
        alt       = Angle(line[72:75], unit=u.degree)
        moon_dist = Angle(line[91:94], unit=u.degree)
        moon_alt  = Angle(line[96:99], unit=u.degree)
        ic(time, ra, dec, mag, motion, alt, az, moon_dist, moon_alt)
        qt.add_row([ time, ra, dec, mag, motion, pa, alt, az, moon_dist, moon_alt ])

        # # Test: compare alt/az in ephemerides to values computed from ra/dec
        # coord=SkyCoord(ra, dec, frame=FK5, equinox="J2000", obstime=time)
        # altaz=coord.transform_to(AltAz(obstime=time, location=Options.loc))
        # ic(alt, az, altaz.alt, altaz.az)                                 

    ic(qt)
    # qt.write(sys.stdout, format="ascii")
    return qt




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
    
        # <hr> marks start of NEOCP ephemerides
        m = re.search(r"<p>(?:</p>)?<hr><p>", line)
        if m:
            ic(line)

            line = next(content_iter, None).strip()
            m = re.match(r"(?:</p>)?<p><b>(.+)</b>", line)
            if m:
                ic(line)
                neocp_id = m.group(1)
                ic(neocp_id)

                # Read lines until </pre>, ephemerides data starts with date
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

                ##DEBUG##
                # Early return, just the 1st entry in the list for debugging
                # return neocp_eph
    return neocp_eph



def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = "Parse NEOCP ephemerides",
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("-L", "--mag-limit", help=f"set mag limit for NEOCP, default {Options.mag_limit}")
    arg.add_argument("-l", "--location", help="MPC station code")
    arg.add_argument("-U", "--update-neocp", action="store_true", help="update NEOCP data from MPC")
    arg.add_argument("-P", "--plot_objects", action="store_true", help="create plot with objects")
    arg.add_argument("-S", "--sky-plot", action="store_true", help="create sky plot (default altitude ploat)")

    args = arg.parse_args()

    if args.debug:
        ic.enable()
        ic(sys.version_info)
        ic(args)
    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()

    if args.location:
        loc = get_location(args.location)
        Options.code = args.location
        Options.loc  = loc
        ic(loc, loc.to_geodetic())
    verbose(f"location: {Options.code} {location_to_string(Options.loc)}")

    if args.update_neocp:
        neocp_query_ephemeris(LOCAL_QUERY)
        neocp_query_list(LOCAL_LIST)

    with open(LOCAL_QUERY, "r") as file:
        content = file.readlines()
        eph_dict = parse_neocp_eph(content)
        table_dict = convert_all_to_qtable(eph_dict)
        print_table_dict(table_dict)
        # table_dict_sorted = sort_by_alt_max_time(table_dict)
        # process_objects(table_dict_sorted)

        # Plot objects and Moon
        if args.plot_objects or args.sky_plot:
            if args.sky_plot:
                plot_sky_objects(table_dict)
            else:
                plot_alt_objects(table_dict)



if __name__ == "__main__":
    main()
