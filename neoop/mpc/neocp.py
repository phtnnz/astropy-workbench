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
# Version 0.1 / 2026-02-08
#       Refactoring complete
# Version 0.2 / 2026-05-18
#       Refactoring for EphemDataList
# Version 1.0 / 2026-06-13
#       Removed locally cached files

VERSION = "1.0 / 2026-06-13"
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
from verbose import verbose, warning, error
from neoconfig import config
from neoclasses import EphemData, Ephem, NEOCPListData, EphemDataList, LocalCircumstances
from neoutils import get_mag0, max_motion
from neoephem import get_dec_limits
import neofiles


# Requests timeout
TIMEOUT = config.requests_timeout

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

def mpc_query_neocp_ephemerides(url: str, local: LocalCircumstances) -> None:
    loc = local.loc
    code = local.code
    ic(url, code)

    # Compute min/max DEC from min altitude and latitude
    mag_limit = config.neocp_mag_limit
    min_alt   = config.min_alt
    min_dec, max_dec = get_dec_limits(local, min_alt * u.deg)
    ic(min_dec, max_dec)
    min_dec, max_dec = int(min_dec.value), int(max_dec.value)

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

    return response.text



def mpc_query_neocp_list(url: str) -> None:
    """
    Retrieve NEOCP/PCCP txt list from MPC

    Parameters
    ----------
    url : str
        MPC web service URL
    filename : str
        Local file name in cache
    """
    ic(url)
    response = requests.get(url, timeout=TIMEOUT)
    ic(response.status_code)
    if response.status_code != 200:
        error(f"query to {url} failed")

    return response.text



def edata_list_from_text_ephemerides(eph_text: dict[str, list[str]], local: LocalCircumstances) -> EphemDataList:
    edata_list = EphemDataList()

    for obj, lines in eph_text.items():
        qt = convert_text_ephemeris1(obj, lines, local)
        if len(qt) == 0:
            verbose(f"skipping NEOCP {obj=} (empty)")
            continue
        eph1 = Ephem.from_table(qt)
        mag = get_mag0(eph1)
        motion = max_motion(eph1)

        edata = EphemData("-", obj, None, eph1, None, None, mag, motion)
        edata_list.append(edata)

    return edata_list



def edata_list_add_neocp_list(edata_list: EphemDataList, neocp_list: dict[str, NEOCPListData], is_pccp: bool=False) -> EphemDataList:
    for edata in edata_list:
        obj = edata.obj
        neocp = neocp_list.get(obj)
        if neocp:
            neocp.type = "PCCP" if is_pccp else "NEOCP"
            edata.neocp = neocp
            edata.type = neocp.type
        elif not is_pccp:
            warning(f"no NEOCP / PCCP list data for {obj}")
    return edata_list



def convert_text_ephemeris1(id: str, eph: list[str], local: LocalCircumstances) -> QTable:
    """
    Convert ephemeris in plain text format to table

    Parameters
    ----------
    id : str
        Object id
    eph : list
        Ephemeris as list of plain text lines from MPC query results
    local : LocalCircumstances
        Location data &c.
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
    qt["Targetname"] = id
    qt["Obstime"]    = Time("2000-01-01 00:00")
    qt["RA"]         = 0 * u.hourangle
    qt["DEC"]        = 0 * u.degree
    qt["Mag"]        = 0 * u.mag
    qt["Motion"]     = 0 * u.arcsec / u.min
    qt["PA"]         = 0 * u.degree
    qt["Alt"]        = 0 * u.degree
    qt["Az"]         = 0 * u.degree
    qt["Moon_dist"]  = 0 * u.degree
    qt["Moon_alt"]   = 0 * u.degree

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
        ic(time, local.naut_dusk, local.naut_dawn)
        if time < local.naut_dusk or time > local.naut_dawn:
            ic("skipping")
            continue
        ic(ra, dec, mag, motion, alt, az, moon_dist, moon_alt)
        qt.add_row([ id, time, ra, dec, mag, motion, pa, alt, az, moon_dist, moon_alt ])

        # # Test: compare alt/az in ephemerides to values computed from ra/dec
        # coord=SkyCoord(ra, dec, frame=FK5, equinox="J2000", obstime=time)
        # altaz=coord.transform_to(AltAz(obstime=time, location=Options.loc))
        # ic(alt, az, altaz.alt, altaz.az)                                 

    ic(qt)
    return qt



def parse_html_ephemerides(content: str) -> dict[str, list[str]]:
    """
    Parse HTML page with the result of the MPC NEOCP ephemerides query

    Parameters
    ----------
    content : str
        Web page content

    Returns
    -------
    dict
        Dictionary with ephemerides in list of plain text format
    """
    neocp_id = None
    neocp_text_eph = {}

    content_iter = iter(content.splitlines())
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
                neocp_id = m.group(1).strip()
                ic(neocp_id)

                # Read lines until </pre>, ephemerides data starts with date
                neocp_text_eph[neocp_id] = []
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
                        neocp_text_eph[neocp_id].append(line)

                ##DEBUG##
                # Early return, just the 1st entry in the list for debugging
                # return neocp_eph
    return neocp_text_eph



def parse_neocp_list(content: str) -> dict[str, NEOCPListData]:
    """
    Parse plain text HTML page of NEOCP/PCCP list

    Parameters
    ----------
    content : str
        Web page content

    Returns
    -------
    dict[str, NEOCPListData]
        Dict of relevant data from NEOCP list
    """
    neocp_list = dict()

    for line in content.splitlines():
        # Format example:
        # ZTF105X  76 2025 09 10.2  22.6070  +0.8481 18.2 Updated Sept. 11.49 UT          26   1.25 24.5  0.038
        # P12e56E 100 2025 08 30.4  23.2926 -12.0924 18.9 Updated Sept. 11.20 UT           6   0.09 23.0 12.007
        # H452509  63 2016 06 04.3  11.5743  +5.6542 34.1 Updated Sept. 10.20 UT           5   0.09 25.6 3386.166
        # ^0      ^8  ^12           ^26     ^34      ^43  ^48                     ???    ^79  ^84   ^90  ^95
        # Temp.   Score             RA      DEC      MagV                         Note   Nobs       H    Not Seen/dys
        # Desig.   	  Discovery                           Updated                             Arc

        # Arc = time difference between 1st and last observation in days

        line = line.rstrip()
        ic(line)
        obj = line[0:7].strip()
        data = NEOCPListData(type="NEOCP", 
                             score=int(line[8:11]),              
                             mag=float(line[43:47]) * u.mag,
                             nobs=int(line[79:82]),
                             arc=float(line[84:89]) * u.day,
                             notseen=float(line[95:]) * u.day)
        ic(obj, data)
        neocp_list[obj] = data

    return neocp_list



def neocp_get_edata_list(local: LocalCircumstances) -> EphemDataList:
    verbose(f"download ephemerides from {config.url_neocp_query}")
    content_ephem = mpc_query_neocp_ephemerides(config.url_neocp_query, local)
    ephemerides_txt = parse_html_ephemerides(content_ephem)
    edata_list = edata_list_from_text_ephemerides(ephemerides_txt, local)

    verbose(f"download NEOCP list from {config.url_neocp_list}")
    content_neocp = mpc_query_neocp_list(config.url_neocp_list)
    neocp_list = parse_neocp_list(content_neocp)
    edata_list_add_neocp_list(edata_list, neocp_list)

    verbose(f"download PCCP list from {config.url_pccp_list}")
    content_pccp  = mpc_query_neocp_list(config.url_pccp_list)
    pccp_list = parse_neocp_list(content_pccp)
    edata_list_add_neocp_list(edata_list, pccp_list, is_pccp=True)

    verbose(f"NEOCP objects ({edata_list.len()}): {edata_list.objects_str()}")
    return edata_list
