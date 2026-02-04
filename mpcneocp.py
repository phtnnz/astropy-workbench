#!/usr/bin/env python

# Copyright 2026 Martin Junius
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
# Version 0.0 / 2026-02-03
#       NEOCP ephemeris and list queries moved from neocp.py

VERSION = "0.0 / 2026-02-03"
AUTHOR  = "Martin Junius"
NAME    = "mpcneocp"

import re
import requests

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.coordinates import Angle, EarthLocation
import astropy.units as u
from astropy.units import Quantity, Magnitude
from astropy.time import Time
from astropy.table import QTable

# Local modules
from verbose import verbose, warning, error, message
from neoconfig import config



# Requests timeout
TIMEOUT = config.requests_timeout
# Exposure times / s
EXP_TIMES = config.exposure_times

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

def mpc_query_neocp_ephemerides(url: str, filename: str, loc: EarthLocation, code: str) -> None:
    """
    Retrieve NEOCP ephemerides from MPC

    Parameters
    ----------
    url : str
        MPC web service URL
    filename : str
        Local file name in cache
    loc : EarthLocation
        Observer location
    code : str
        MPC station code
    """
    ic(url, filename, code)

    # Compute min/max DEC from min altitude and latitude
    lat = loc.lat.degree
    mag_limit = config.mag_limit
    min_alt = config.min_alt

    min_dec = -90 + min_alt + lat
    max_dec = +90 - min_alt + lat
    min_dec = -90 if min_dec < -90 else int(min_dec)
    max_dec = +90 if max_dec > +90 else int(max_dec)
    ic(lat, min_dec, max_dec)

    data = { 
        "W":	"a",            # all objects with ...
        "mb":	"-30",          # max brightness (V)
        "mf":	mag_limit,      # min brightness (V)
        "dl":	min_dec,        # min DEC
        "du":	max_dec,        # max DEC
        "nl":	"0",            # min NEO score
        "nu":	"100",          # max NEO score
        "sort":	"d",
        "Parallax":	"1",        # 0=geocentric, 1=code, 2=Lon/Lat/Alt
        "obscode": code,        # observatory code
        "long":	"",
        "lat":	"",
        "alt":	"",
        "int":	"1",            # interval 0=1h, 1=30m, 2=10m, 3=1m
        "start":"1",            # start now + X hours
        "raty":	"a",            # format: h=trunc. sexagesimal, a=full, d=decimal
        "mot":	"m",            # motion: s="/sec, m="/min, h="/hr, d=degree/day
        "dmot":	"p",            # motion: p=total, r=separate RA/DEC coord. mo tion, s=separate RA/DEC sky motion
        "out":	"f",            # f=full output, b=brief output
        "sun":	"n",            # suppress output: x=never, s=sunset/rise, c=civil, n=nautical, a=astronomical
        "oalt":	min_alt         # suppress below min altitude
    }

    ic(url, data)
    # For some reason GET works, but POST doesn't?
    response = requests.get(url, params=data, timeout=TIMEOUT)
    ic(response.status_code)
    if response.status_code != 200:
        error(f"query to {url} failed")

    with open(filename, mode="w", encoding=response.encoding) as file:
        file.write(response.text)



def mpc_query_neocp_list(url: str, filename: str) -> None:
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



def print_ephemerides(ephemerides: dict) -> None:
    """
    Print ephemerides for all objects

    Parameters
    ----------
    ephemerides : dict
        Ephemerides dictionary
    """
    for id, qt in ephemerides.items():
        verbose("===================================================================================================================")
        verbose(f"NEOCP {id} ephemerides")
        verbose.print_lines(qt)
    verbose("===================================================================================================================")



def convert_text_ephemerides(eph_dict: dict, min_time: Time, max_time: Time) -> dict:
    """
    Convert ephemerides from plain text format to table for all objects

    Parameters
    ----------
    eph_dict : dict
        Ephemerides dictionary in plain text format
    min_time : Time
        Lower limit for ephemeris time
    max_time : Time
        Upper limit for ephemeris time

    Returns
    -------
    dict
        Ephemerides dictionary in table form for further processing
    """
    qtable_dict = {}
    for id, eph in eph_dict.items():
        qt = convert_text_ephemeris1(id, eph, min_time, max_time)
        if len(qt["mag"]) == 0:
            verbose(f"skipping NEOCP {id=} (empty)")
            continue
        qtable_dict[id] = qt
    return qtable_dict


def convert_text_ephemeris1(id: str, eph: list, min_time: Time, max_time: Time) -> QTable:
    """
    Convert ephemeris in plain text format to table

    Parameters
    ----------
    id : str
        Object id
    eph : list
        Ephemeris as list of plain text lines from MPC query results
    min_time : Time
        Lower limit for ephemeris time
    max_time : Time
        Upper limit for ephemeris time

    Returns
    -------
    QTable
        Ephemeris table
    """
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
        ic(time, min_time, max_time)
        if time < min_time or time > max_time:
            ic("skipping")
            continue
        ic(ra, dec, mag, motion, alt, az, moon_dist, moon_alt)
        qt.add_row([ time, ra, dec, mag, motion, pa, alt, az, moon_dist, moon_alt ])

        # # Test: compare alt/az in ephemerides to values computed from ra/dec
        # coord=SkyCoord(ra, dec, frame=FK5, equinox="J2000", obstime=time)
        # altaz=coord.transform_to(AltAz(obstime=time, location=Options.loc))
        # ic(alt, az, altaz.alt, altaz.az)                                 

    ic(qt)
    return qt



def parse_html_ephemerides(content: list) -> dict:
    """
    Parse HTML page with the result of the MPC NEOCP ephemerides query

    Parameters
    ----------
    content : list
        Web page content

    Returns
    -------
    dict
        Ephemerides dictionary in plain text format
    """
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
    """
    Parse plain text HTML page of NEOCP/PCCP list

    Parameters
    ----------
    content : list
        Web page content

    Returns
    -------
    dict
        Dictionary of relevant data from NEOCP list
    """
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
