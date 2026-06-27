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
# Version 0.1 / 2026-06-22
#       Moved MPC functions from jpl.sbwobs

VERSION     = "0.1 / 2026-06-22"
AUTHOR      = "Martin Junius"
NAME        = "mpc.lastobs"
DESCRIPTION = "Retrieve lists of last observations from MPC"

import requests
import re

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.time import Time
import astropy.units as u

# Local modules
from utils.verbose import verbose, warning, error, message
from neo.config import config
from neo.classes import MPCDLxData, EphemData



def mpc_query_customize(url: str, list_type: str) -> str:
    ic(url)

    timeout = config.requests_timeout

    data = { 
        "ra1":          "00 00",                    # RA limits
        "ra2":          "24 00",
        "dec1":         f"{config.min_dec:+02d} 00",# DEC limits
        "dec2":         f"{config.max_dec:+02d} 00",
        "elong1":       60,                         # Elongation limits
        "elong2":       180,
        "mag1":         0,                          # V mag limits
        "mag2":         config.sbwobs_mag_limit,
        "to":           1,                          # Single-opposition unnumbered objects
        "wh":           list_type                   # "DLN" | "DLU"
        # Valid list types:
        # *DLN   Dates Of Last Observation Of NEOs
        #  DLNR  Dates Of Last Observation Of NEOs (R.A. order)
        #  BNR   Bright NEO Recovery Opportunities
        #  FNR   Faint NEO Recovery Opportunities
        # *DLU   Dates Of Last Observation Of Unusual Minor Planets
        #  DLD   Dates Of Last Observation Of Distant Objects
        #  DLDR  Dates Of Last Observation Of Distant Objects (R.A. order)
        # * = used with this script
    }

    ic(url, data)
    verbose(f"query {url}")
    response = requests.get(url, params=data, timeout=timeout)
    ic(response.status_code)
    if response.status_code != 200:
        error(f"query to {url} failed")

    return response.text



def parse_date(date: str) -> str:
    months = {"jan": 1, "feb": 2, "mar": 3, "apr": 4, "may": 5, "jun": 6,
              "jul": 7, "aug": 8, "sep": 9, "oct": 10, "nov": 11, "dec": 12}
    m = re.match(r"(\d{4}) ([a-z]{3}). ([ 0-9]{2})$", date.lower())
    if m:
        year  = int(m.group(1))
        month = months[m.group(2)]
        day   = int(m.group(3))
        return f"{year:04d}-{month:02d}-{day:02d}"
    else:
        raise ValueError(date)



def parse_txt_line(txt_line: str, list_type: str) -> MPCDLxData:
    ic(txt_line)
##DLU format examples
#         2025 WA             Amo   1.3/+16/18.3/147/07.42  *2025 Nov. 16  C23  17.8 c  7    1  2.50 0.60   5
#         2025 VD6            Amo   3.1/+30/20.3/168/04.63  *2025 Nov. 16  958  19.7 G  8    3  1.69 0.40  18
#         ^9                  ^29  ^34                  57^ ^^60           ^74  ^79     ^87^90  ^95
    if list_type == "DLU":
        # Fix one char shift if "Motn" is 3 digits
        if txt_line[57:59] != "  ":
            if txt_line[58:60] == "  ":
                txt_line = txt_line[:58] + txt_line[59:]
            else:
                warning("bad format in data txt_line")
                warning(txt_line)
                return None
        designation = txt_line[ 9:27].strip()
        type        = txt_line[29:32]
        currently   = txt_line[34:58].strip()
        marker      = txt_line[59:60].strip()
        last_obs    = parse_date(txt_line[60:72])
        code        = txt_line[74:77]
        mag         = txt_line[79:84].strip()
        filter      = txt_line[84:85].strip()
        uncertainty = txt_line[87:88]
        arc         = txt_line[90:93].strip()
        ic(designation, type, currently, last_obs, code, mag, uncertainty, arc)

##DLN format examples
#  Designation                         Currently              Last obs.  Code   Mag.   U  Arc      Currently + 7d           Currently + 30d          Currently + 60 d
#                                   R.A./Dec/  V /El./Motn                                       R.A./Dec/  V /El./Motn   R.A./Dec/  V /El./Motn   R.A./Dec/  V /El./Motn
# Sample:
#         2025 WA             Amo   1.3/+16/18.3 /147/07.42  2025 Nov. 16  C23  17.8 c  7    1   6.4/+01/18.3/141/08.40   9.1/-08/22.2/121/00.19   8.9/-04/23.3/151/00.30
#         ^9                  ^29  ^34                   58^ ^60           ^74  ^79     ^87^90  ^95
    elif list_type == "DLN":
        # Fix one char shift if "Motn" is 3 digits
        if txt_line[58:60] != "  ":
            if txt_line[59:61] == "  ":
                txt_line = txt_line[:59] + txt_line[60:]
            else:
                warning("bad format in data line")
                warning(txt_line)
                return None
        designation = txt_line[ 9:27].strip()
        type        = txt_line[29:32]
        currently   = txt_line[34:59].strip()
        marker      = ""
        last_obs    = parse_date(txt_line[60:72])
        code        = txt_line[74:77]
        mag         = txt_line[79:84].strip()
        filter      = txt_line[84:85].strip()
        uncertainty = txt_line[87:88]
        arc         = txt_line[90:93].strip()
        ic(designation, type, currently, last_obs, code, mag, uncertainty, arc)

    else:
        error(f"unknown list type {list_type}")

##Common
    ra, dec, mag1, elongation, motion = currently.split("/")
    # ic(ra, dec, mag1, elongation, motion)
    dlx  = MPCDLxData(
        designation,
        type,
        float(ra) * u.hourangle,
        float(dec) * u.deg,
        float(mag1) if mag1.strip() else 0,
        float(elongation) * u.deg,
        float(motion) * u.arcsec / u.min,
        marker,
        Time(last_obs),
        code,
        float(mag) if mag else 0,
        filter,
        # Uncertainty may be empty ' '
        int(uncertainty) if uncertainty.isdigit() else 0,
        float(arc) * u.day
    )
    ic(dlx)
    return dlx



def mpc_parse_customize(content: str, list_type: str) -> dict[str, EphemData]:
    objects = dict()

    in_unnumbered1 = False
    for line in content.splitlines():
        m = re.match(r'<h2>(.*)</h2>', line)
        if m:
            ic(m)
            h2 = m.group(1)
            in_unnumbered1 = True if h2 == "One-Opposition Unnumbered Objects" else False
            continue
        if not in_unnumbered1:
            continue
        m = re.match(r'<input type="checkbox" name="Obj" VALUE=".+">(.*)$', line)
        if not m:
            continue
        ic(m)
        txt_line = m.group(1)
        # print(txt_line)
        dlx: MPCDLxData = parse_txt_line(txt_line, list_type)
        if dlx:
            edata = EphemData("", dlx.designation, None, None, None, None, dlx.vmag, None, dlx=dlx)
            ic(edata)
            objects[dlx.designation] = edata

    return objects



def mpc_query_lastobs(url: str) -> str:
    ic(url)

    timeout = config.requests_timeout

    ic(url)
    verbose(f"query {url}")
    response = requests.get(url, timeout=timeout)
    ic(response.status_code)
    if response.status_code != 200:
        error(f"query to {url} failed")

    return response.text



def mpc_parse_lastobs(content: str) -> dict[str, EphemData]:
    objects = dict()

    for txt_line in content.splitlines():
        dlx = parse_txt_line(txt_line, "DLN")       # same format als customize/DLN
        if dlx:
            # Check mag limit
            mag = dlx.last_mag or 99.99
            if mag > config.sbwobs_mag_limit:
                continue
            # Check DEC limit
            dec = dlx.dec.value or 0.0
            if dec < config.min_dec or dec > config.max_dec:
                continue
            ##FIXME: check last obs data here?
            edata = EphemData("", dlx.designation, None, None, None, None, dlx.vmag, None, dlx=dlx)
            objects[dlx.designation] = edata

    return objects



