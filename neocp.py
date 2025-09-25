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
# Version 0.2 / 2025-09-12
#       Parse NEOCP and PCCP lists for additional data, separate altitude/sky plots
# Version 0.3 / 2025-09-17
#       Added CSV output for nina-create-sequence2
# Version 0.4 / 2025-09-18
#       Major improvement to planning, avoid overlap of observations,
#       altitude/sky plot only for planned objects

VERSION = "0.4 / 2025-09-18"
AUTHOR  = "Martin Junius"
NAME    = "neocp"

import sys
import argparse
import csv
import re
import requests
from datetime import datetime, timezone
from zoneinfo import ZoneInfo
from typing import Tuple
from itertools import pairwise
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
from astropy.table import QTable, Row

# Astroplan
import matplotlib.pyplot as plt
from astroplan import Observer
from astroplan.plots import plot_altitude, plot_sky

# Local modules
from verbose import verbose, warning, error
from astroutils import get_location, location_to_string


# Requests timeout
TIMEOUT = 60

# Exposure times / s
EXP_TIMES = [ 5, 10, 15, 20, 30, 45, 60 ]


# MPC pages:
# NEOCP form
#       https://minorplanetcenter.net/iau/NEO/toconfirm_tabular.html
# PCCP form
#       https://minorplanetcenter.net/iau/NEO/pccp_tabular.html
# NEOCP list (contains also PCCP objects)
#       https://minorplanetcenter.net/iau/NEO/neocp.txt
# PCCP list
#       https://minorplanetcenter.net/iau/NEO/pccp.txt
#
# NEOCP/PCCP query ephemerides
#       https://cgi.minorplanetcenter.net/cgi-bin/confirmeph2.cgi
#
# Example request for M49
# W=a&mb=-30&mf=20.5&dl=-90&du=%2B40&nl=75&nu=100&sort=d&Parallax=1&obscode=M49&long=&lat=&alt=&int=1&start=0&raty=a&mot=m&dmot=p&out=f&sun=n&oalt=26

URL_NEOCP_QUERY = "https://cgi.minorplanetcenter.net/cgi-bin/confirmeph2.cgi"
LOCAL_EPH = "NEOCP-ephemerides.html"

URL_NEOCP_LIST = "https://minorplanetcenter.net/iau/NEO/neocp.txt"
LOCAL_NEOCP  = "NEOCP-list.txt"

URL_PCCP_LIST = "https://minorplanetcenter.net/iau/NEO/pccp.txt"
LOCAL_PCCP   = "PCCP-list.txt"

DOWNLOADS = "./downloads/"


# Command line options
class Options:
    """
    Command line options
    """
    csv: bool = False       # -C --csv
    output: str = None      # -o --output

    ##FIXME: use config file
    mag_limit = 20.5
    pixel_tolerance = 2     # Max trail tolerance in pixels
    resolution = 1.33       # arcsec / pixel resolution of camera/telescope
    code = "M49"            # -l --location
    loc = get_location(code)
    min_alt = 26            # min altitude
    dead_time_slew_center = 90 * u.s
    dead_time_af          = 100 * u.s
    dead_time_image       = 1.5 * u.s
    dead_time_guiding     = 30 * u.s
    safety_margin         = 5 * u.min
    min_n_obs = 4
    max_notseen = 3 * u.day
    opt_alt = 50 * u.degree # min altitude for optimal results
    min_n_exp = 15          # min number of exposures
    max_n_exp = 60          # max number of exposures
    min_perc_required = 50  # min percentage of required total exposure time
    min_moon_dist = 50 * u.degree     # min distance from moon



def mpc_query_ephemerides(url: str, filename: str) -> None:
    """
    Retrieve NEOCP ephemerides from MPC

    Parameters
    ----------
    url : str
        MPC web service URL
    filename : str
        Local file name in cache
    """
    ##FIXME: get query params from config file
    ic(url, filename)

    # Compute min/max DEC from min altitude and latitude
    lat = Options.loc.lat.degree
    min_dec = -90 + Options.min_alt + lat
    max_dec = +90 - Options.min_alt + lat
    min_dec = -90 if min_dec < -90 else int(min_dec)
    max_dec = +90 if max_dec > +90 else int(max_dec)
    ic(lat, min_dec, max_dec)

    data = { 
        "W":	"a",                # all objects with ...
        "mb":	"-30",              # max brightness (V)
        "mf":	Options.mag_limit,  # min brightness (V)
        "dl":	min_dec,            # min DEC
        "du":	max_dec,            # max DEC
        "nl":	"0",                # min NEO score
        "nu":	"100",              # max NEO score
        "sort":	"d",
        "Parallax":	"1",            # 0=geocentric, 1=code, 2=Lon/Lat/Alt
        "obscode":	Options.code,   # observatory code
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
    # For some reason GET works, but POST doesn't?
    response = requests.get(url, params=data, timeout=TIMEOUT)
    ic(response.status_code)
    if response.status_code != 200:
        error(f"query to {url} failed")

    with open(filename, mode="w", encoding=response.encoding) as file:
        file.write(response.text)



def mpc_query_list(url: str, filename: str) -> None:
    """
    Retrieve NEOCP/PCCP txt list from MPC

    Parameters
    ----------
    url : str
        MPC web service URL
    filename : str
        Local file name in cache
    """
    ic(url, filename)
    response = requests.get(url, timeout=TIMEOUT)
    ic(response.status_code)
    if response.status_code != 200:
        error(f"query to {url} failed")

    with open(filename, mode="w", encoding=response.encoding) as file:
        file.write(response.text)



def qtable_to_altaz(id: str, qt: QTable) -> AltAz:
    """
    Convert table rows "alt"/"az" to AltAz object

    Parameters
    ----------
    id : str
        NEOCP id = temporary designation
    qt : QTable
        Table with ephemeris, including alt/az

    Returns
    -------
    AltAz
        AltAz coordinates object for altitude/sky plot
    """
    altaz = AltAz(alt=qt["alt"], az=qt["az"], obstime=qt["obstime"], location=Options.loc)
    # Quick hack to get a proper label for plot_altitude()
    altaz.name = id
    ic(altaz)
    return altaz


def plot_alt_objects(ephemerides: dict, objects: list, filename: str) -> None:
    """
    Generate altitude plot

    Parameters
    ----------
    ephemerides : dict
        Dictionary with ephemerides for all objects
    objects : list
        Subset list of objects for plot
    filename : str
        File name for generated PNG
    """
    # Get next midnight
    time = Time(Time.now(), location=Options.loc)
    observer = Observer(location=Options.loc, description=Options.code)
    midnight = observer.midnight(Time.now(), which="next")
    ic(midnight)

    # Intervals around midnight
    time_interval = midnight + np.linspace(-8, 8, 160)*u.hour
    moon          = observer.moon_altaz(time_interval)
    # Quick hack to get a proper label for plot_altitude()
    moon.name     = "Moon"

    # Plot all NEOCP objects
    for id in objects:
        qt = ephemerides[id]
        altaz = qtable_to_altaz(id, qt)
        plot_altitude(altaz, observer, altaz.obstime, style_kwargs=dict(fmt="o"))

    plot_altitude(moon, observer, moon.obstime, brightness_shading=True, style_kwargs=dict(fmt="y--"))
    plt.legend(bbox_to_anchor=(1.0, 1.02))

    plt.savefig(f"tmp/{filename}", bbox_inches="tight")
    plt.close()


def plot_sky_objects(ephemerides: dict, objects: list, filename: str) -> None:
    """
    Generate sky plot

    Parameters
    ----------
    ephemerides : dict
        Dictionary with ephemerides for all objects
    objects : list
        Subset list of objects for plot
    filename : str
        File name for generated PNG
    """
    # Get next midnight
    time = Time(Time.now(), location=Options.loc)
    observer = Observer(location=Options.loc, description=Options.code)
    midnight = observer.midnight(Time.now(), which="next")
    ic(midnight)

    # Intervals around midnight
    time_interval = midnight + np.linspace(-5, 5, 11)*u.hour
    moon          = observer.moon_altaz(time_interval)
    # Quick hack to get a proper label for plot_sky()
    moon.name     = "Moon"

    # Plot all NEOCP objects
    for id in objects:
        qt = ephemerides[id]
        altaz = qtable_to_altaz(id, qt)
        plot_sky(altaz, observer, altaz.obstime)

    plot_sky(moon, observer, moon.obstime, style_kwargs=dict(color="y", marker="x"))
    plt.legend(bbox_to_anchor=(1.48, 1.11))

    plt.savefig(f"tmp/{filename}", bbox_inches="tight")
    plt.close()



def max_motion(qt: QTable) -> Quantity:
    """
    Get max value for "motion" from table

    Parameters
    ----------
    qt : QTable
        Ephemeris table

    Returns
    -------
    Quantity
        Max motion as arcsec/min
    """
    max_motion = -1 * u.arcsec / u.min
    for row in qt:
        if row["motion"] > max_motion:
            max_motion = row["motion"]
    return max_motion

def exp_time_from_motion(motion: Quantity) -> Quantity:
    """
    Compute exposure time depending on motion value,

    Parameters
    ----------
    motion : Quantity
        Motion value as arcsec/min

    Returns
    -------
    Quantity
        Exposure time as secs
    """
    exp = Options.pixel_tolerance * Options.resolution * u.arcsec / motion
    exp_max = EXP_TIMES[0]
    for exp1 in EXP_TIMES:
        if exp1 > exp.to(u.s).value:
            break
        exp_max = exp1
    return exp_max * u.s



def is_east(az: Angle) -> bool:
    """
    Test for east azimut position

    Parameters
    ----------
    az : Angle
        Azimut angle

    Returns
    -------
    bool
        True if east, False if west
    """
    az180     = 180 * u.degree
    return True if az >= 0 and az < az180 else False

def is_west(az: Angle) -> bool:
    """
    Test for west azimut position

    Parameters
    ----------
    az : Angle
        Azimut angle

    Returns
    -------
    bool
        True if west, False if east
    """
    az180     = 180 * u.degree
    az360     = 360 * u.degree
    return True if az >= az180 and az < az360 else False

def flip_times(qt: QTable) -> Tuple[Time, Time]:
    """
    Get time before and after meridian passing from ephemeris table

    Parameters
    ----------
    qt : QTable
        Ephemeris table

    Returns
    -------
    Tuple[Time, Time]
        Time before meridian, time after meridian passing
        None, None if object doesn't pass meridian
    """
    prev_time = None
    prev_az   = None

    for row in qt:
        time = row["obstime"]
        az   = row["az"]
        # ic(time, az)
        if not prev_az == None:
            if is_east(prev_az) and is_west(az):     # South flip
                return (prev_time, time)
            if is_west(prev_az) and is_east(az):     # North flip
                return (prev_time, time)
        prev_time = time
        prev_az   = az

    # No meridian passing found
    return None, None



def max_alt_times(qt: QTable) -> Tuple[Time, Time]:
    """
    Get time of maximum altitude from ephemeris table

    Parameters
    ----------
    qt : QTable
        Ephemeris table

    Returns
    -------
    Tuple[Time, Time]
        Time of max altitude, Time of max altitude
        (twice for return tuple to be compatible with flip_times())
    """
    max_alt = -90 * u.degree
    time_max = None
    for row in qt:
        if row["alt"] > max_alt:
            max_alt = row["alt"]
            time_max = row["obstime"]
    # Return tuple for compatibility with flip_times()
    return time_max, time_max



def get_row_for_time(qt: QTable, t: Time) -> Row:
    """
    Get row from ephemeris table closest to specified time

    Parameters
    ----------
    qt : QTable
        Ephemeris table
    t : Time
        Time for retrieving row

    Returns
    -------
    Row
        Corresponding row from ephemeris table
        None, if not found
    """
    for r1, r2 in pairwise(qt):
        if r1["obstime"] <= t and t <= r2["obstime"]:
            return r1
    # Not matching interval found
    return None



def opt_alt_times(qt: QTable, alt: Angle) -> Tuple[Time, Time]:
    """
    Get times for object above specified altitude

    Parameters
    ----------
    qt : QTable
        Ephemeris table
    alt : Angle
        Altitude angle

    Returns
    -------
    Tuple[Time, Time]
        Start time above altitude, end time above altitude
    """
    time_alt0 = None
    time_alt1 = None

    for row in qt:
        if time_alt0 == None and row["alt"] >= alt:
            time_alt0 = row["obstime"]
        if time_alt0 != None and row["alt"] >= alt:
            time_alt1 = row["obstime"]
        if time_alt1 != None and row["alt"] < alt:
            break
    
    return time_alt0, time_alt1



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
                time_before, time_after = max_alt_times(qt)
        time_0 = qt["obstime"][0]
        time_1 = qt["obstime"][-1]
        time_alt0, time_alt1 = opt_alt_times(qt, Options.opt_alt)
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

    verbose("             Score      MagV #Obs      Arc NotSeen  Time start ephemeris   / end ephemeris                  Max motion")
    verbose("                                                    Time before            / after meridian")
    verbose("                                                    Time start exposure    / end exposure")
    verbose("                                                    # x Exp   = total exposure time")
    verbose("                                                    RA, DEC")

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

        # Calculate single exposure, number of exposures, total exposure, total time
        max_m = max_motion(qt)
        exp   = exp_time_from_motion(max_m)         # Single exposure

        rel_brightness = 10 ** (0.4 * (mag.value - 18))
        total_exp = 4 * u.min * rel_brightness      # Total exposure
        n_exp = int(total_exp / exp) + 1            # Number of exposures
        perc_of_required = 100.                     # Percentage actual / total exposure
        if n_exp < Options.min_n_exp:
            perc_of_required = Options.min_n_exp / n_exp * 100
            n_exp = Options.min_n_exp
        if n_exp > Options.max_n_exp:
            perc_of_required = Options.max_n_exp / n_exp * 100
            n_exp = Options.max_n_exp
        total_exp = (n_exp * exp).to(u.min)
        total_time = (  total_exp 
                      + Options.dead_time_slew_center + Options.dead_time_slew_center 
                      + Options.dead_time_guiding + Options.safety_margin
                      + n_exp * Options.dead_time_image )
        ic(n_exp, exp, total_exp, total_time, perc_of_required)

        verbose("----------------------------------------------------------------------------------------------------------------------")
        verbose(f"{id}  {type:5s} {score:3d}  {mag}  {nobs:3d}  {arc:5.2f}  {notseen:4.1f}  {time_first}/{time_last}  {max_m:4.1f}")

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
        if time_end_exp > time_before and time_start_exp < time_after:
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
        if nobs < Options.min_n_obs:
            warning(f"{id}: SKIPPED: only {nobs} obs (< {Options.min_n_obs})")
            continue

        # Skip, if not seen for more than threshold days
        if notseen > Options.max_notseen:
            warning(f"{id}: SKIPPED: not seen for {notseen:.1f} (> {Options.max_notseen:.1f})")
            continue

        # Skip, if percentage of total exposure time is less than threshold
        if perc_of_required < Options.min_perc_required:
            warning(f"{id}: SKIPPED: only {perc_of_required:.0f}% of required total exposure time (< {Options.min_perc_required}%)")
            continue

        # Skip, if arc is less than threshold
        if arc < 0.05 * u.day:
            warning(f"{id}: SKIPPED: arc {arc:.2f} too small (< 0.05 d)")
            continue

        # Skip, if failed to allocate total_time
        if time_end_exp > time_last:
            warning(f"{id}: SKIPPED: can't allocate exposure time {total_time} ({time_first} -- {time_last})")
            continue

        # Table row best matching time_start_exp
        row = get_row_for_time(qt, time_start_exp)
        ic(row)
        moon_dist = row["moon_dist"]
        # Skip, if moon distance is too small
        if moon_dist < Options.min_moon_dist:
            warning(f"{id}: SKIPPED: moon distance {moon_dist:.0f} < {Options.min_moon_dist}")
            continue

        ##### Good to go! #####
        # Remember end of exposure
        prev_time_end_exp = time_end_exp
        # Append to list of planned objects
        objects.append(id)
        ra, dec = row["ra"], row["dec"]
        alt, az = row["alt"], row["az"]

        verbose(f"                                                    {time_before}/{time_after}")
        verbose(f"                                                    {time_start_exp}/{time_end_exp}")
        total = f"{n_exp} x {exp:2.0f} = {total_exp:3.1f} ({perc_of_required:.0f}%) / total {total_time:3.1f}"
        verbose(f"                                                    {total}")
        verbose(f"                                                    RA {ra:.4f}, DEC {dec:.4f}, Alt {alt:.0f}, Az {az:.0f}")

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

    # Output to CSV file for nina-create-sequence2
    fieldnames = [  "start time", "end time", 
                    "target", "obstime", "ra", "dec", "exposure", "number", "filter",
                    "type", "mag", "nobs", "arc", "notseen", "total" ]
    if Options.csv:
        ##FIXME: use improved csvoutput module
        if Options.output:
            with open(Options.output, "w", newline="") as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(csv_rows)
        else:
            writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(csv_rows)

    # Return list of planned objects
    return objects



def sort_by_flip_time(table_dict: dict) -> dict:
    time_dict = {}

    for id, qt in table_dict.items():
        t, _ = flip_times(qt)
        if not t:
            t, _ = max_alt_times(qt)
        time_dict[id] = t
    
    # Sort dict by time (item[0] = id, item[1] = time)
    time_sorted = { id: time for id, time in sorted(time_dict.items(), key=lambda item: item[1]) }

    # Return table_dict sorted by time
    return { id: table_dict[id] for id in time_sorted.keys() }



def print_table_dict(table_dict: dict) -> None:
    for id, qt in table_dict.items():
        print("===================================================================================================================")
        print(f"NEOCP {id} ephemerides")
        print(qt)
    print("===================================================================================================================")



def convert_eph_list_to_qtable(eph_dict: dict) -> dict:
    qtable_dict = {}
    for id, eph in eph_dict.items():
        qt = eph_to_qtable(id, eph)
        if len(qt["mag"]) == 0:
            verbose(f"skipping NEOCP {id=} (empty)")
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
        # NEOCP ephemerides count azimuth from S=0, Astropy from N=0
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




def parse_html_eph(content: list) -> dict:
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



def parse_neocp_list(content: list) -> dict:
    neocp_list = {}

    for line in content:
        # Format example:
        # ZTF105X  76 2025 09 10.2  22.6070  +0.8481 18.2 Updated Sept. 11.49 UT          26   1.25 24.5  0.038
        # P12e56E 100 2025 08 30.4  23.2926 -12.0924 18.9 Updated Sept. 11.20 UT           6   0.09 23.0 12.007
        # H452509  63 2016 06 04.3  11.5743  +5.6542 34.1 Updated Sept. 10.20 UT           5   0.09 25.6 3386.166
        # ^0      ^8  ^12           ^26     ^34      ^43  ^48                     ???    ^79  ^84   ^90  ^95
        # Temp.   Score             RA      DEC      MagV                         Note   Nobs       H    Not Seen/dys
        # Desig.   	  Discovery                           Updated                             Arc
        line = line.rstrip()
        ic(line)
        vals = {}
        vals["score"] = int(line[8:11])
        vals["mag"]   = float(line[43:47]) * u.mag
        vals["nobs"]  = int(line[79:82])
        # Time difference between 1st and last observation in days
        vals["arc"]   = float(line[84:89]) * u.day
        vals["notseen"] = float(line[95:]) * u.day

        id = line[0:7]
        ic(id, vals)
        neocp_list[id] = vals

    return neocp_list



def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = "Parse NEOCP ephemerides",
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("-l", "--location", help="MPC station code")
    arg.add_argument("-U", "--update-neocp", action="store_true", help="update NEOCP data from MPC")
    arg.add_argument("-A", "--alt-plot", action="store_true", help="create altitude plot with objects")
    arg.add_argument("-S", "--sky-plot", action="store_true", help="create sky plot with objects")
    arg.add_argument("-o", "--output", help="write CSV to OUTPUT file")
    arg.add_argument("-C", "--csv", action="store_true", help="use CSV output format")

    args = arg.parse_args()

    if args.debug:
        ic.enable()
        ic(sys.version_info, sys.path)
        ic(args)
    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()

    Options.csv    = args.csv
    Options.output = args.output

    if args.location:
        loc = get_location(args.location)
        Options.code = args.location
        Options.loc  = loc
        ic(loc, loc.to_geodetic())
    verbose(f"location: {Options.code} {location_to_string(Options.loc)}")

    prefix = Time.now().strftime("%Y%m%d-")
    local_eph   = DOWNLOADS + prefix + LOCAL_EPH
    local_neocp = DOWNLOADS + prefix + LOCAL_NEOCP
    local_pccp  = DOWNLOADS + prefix + LOCAL_PCCP
    ic(prefix, local_eph, local_neocp, local_pccp)

    if args.update_neocp:
        verbose(f"download ephemerides from {URL_NEOCP_QUERY}")
        mpc_query_ephemerides(URL_NEOCP_QUERY, local_eph)
        verbose(f"download NEOCP list from {URL_NEOCP_LIST}")
        mpc_query_list(URL_NEOCP_LIST, local_neocp)
        verbose(f"download PCCP list from {URL_PCCP_LIST}")
        mpc_query_list(URL_PCCP_LIST, local_pccp)

    # Parse ephemerides
    verbose(f"processing {local_eph}")
    with open(local_eph, "r") as file:
        content = file.readlines()
        ephemerides_txt = parse_html_eph(content)
        ephemerides = convert_eph_list_to_qtable(ephemerides_txt)
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

    ephemerides = sort_by_flip_time(ephemerides)

    verbose("processing objects")
    print_table_dict(ephemerides)
    objects = process_objects(ephemerides, neocp_list, pccp_list, times)

    # Plot objects and Moon
    if args.sky_plot:
        verbose("sky plot for objects")
        plot_sky_objects(ephemerides, objects, "neocp-plot-sky.png")
    if args.alt_plot:
        verbose("altitude plot for objects")
        plot_alt_objects(ephemerides, objects, "neocp-plot-alt.png")



if __name__ == "__main__":
    main()
