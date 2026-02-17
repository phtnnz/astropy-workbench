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
# Version 0.0 / 2025-11-14
#       Retrieve "What's Observable?" from JPL
# Version 0.1 / 2025-11-18
#       Added MPC "Dates Of Last Observation Of NEOs", "Dates Of Last Observation 
#       Of Unusual Minor Planets", output observable objects also in "Unusual" list
# Version 0.2 / 2025-12-28
#       New options --asteroids / --comets / --neo / --pha for sbwobs
#       object selection
# Version 0.3 / 2026-01-03
#       New option -M --mag-limit to override values in config, filter out last obs
#       older than 14 days or uncertainty < 3

VERSION     = "0.3 / 2026-01-03"
AUTHOR      = "Martin Junius"
NAME        = "sbwobs"
DESCRIPTION = "Retrieve observable NEOs/comets from JPL/MPC"

import sys
import argparse
import requests
import json
import re
from typing import Tuple, Any

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
# from astropy.coordinates import AltAz, EarthLocation, SkyCoord
# from astropy.coordinates import Angle
# from astropy.coordinates import errors
from astropy.time import Time
import astropy.units as u
# import astropy.constants as const
# import numpy as np

# Local modules
from verbose import verbose, warning, error, message
from astroutils import location_to_string, get_location
from neoconfig import config
from neoutils import fmt_time

DEFAULT_LOCATION = "M49"



def jpl_query_sbwobs(url: str, start: Time=None, end: Time=None) -> str:
    """Retrieve What's Observable from JPL via API, see
    https://ssd.jpl.nasa.gov/tools/sbwobs.html#/
    https://ssd-api.jpl.nasa.gov/doc/sbwobs.html

    Parameters
    ----------
    url : str
        JPL API URL
    start : Time, optional
        Observation time start (obs-time), by default None = now
    end : Time, optional
        Observation end time (end-time), by default None

    Returns
    -------
    str
        Query result, to be parsed as JSON
    """
    ic(url, start, end)

    timeout = config.requests_timeout

    data = { 
        "mpc-code":     config.mpc_code,        # MPC observatory code
        "obs-time":     fmt_time(Time.now()),   # Date/time of the observation
        "elev-min":     config.elev_min,        # Minimum altitude
        "vmag-max":     config.vmag_max,        # Max V mag = minimum brightness
        "output-sort":  config.output_sort,     # Sort records

        "sb-ns":        config.sb_ns,           # Numbered (n) ./. unnumbered (u)
        "sb-kind":      config.sb_kind,         # Asteroids (a) ./. comets (c)
        # "sb-group":     sb_group              # NEOs (neo) ./. PHA (pha)
    }
    if config.sb_group:
        data["sb-group"] = config.sb_group
    if start:
        data["obs-time"] = fmt_time(start)
    if end:
        data["obs-end"] = fmt_time(end)

    ic(url, data)
    verbose(f"query {url}")
    response = requests.get(url, params=data, timeout=timeout)
    ic(response.status_code)
    if response.status_code != 200:
        error(f"query to {url} failed")

    return response.text



def jpl_parse_sbwobs(text: str) -> dict:
    """
    Parse JSON retrieve from JPL What's Observable

    Parameters
    ----------
    text : str
        String to load JSON from

    Returns
    -------
    dict
        Parsed dictionary, key is object designation
    """
    obj = json.loads(text)
    
    ic(obj.keys())
    # Top-level keys: 'signature', 'obs_constraints', 'sb_constraints', 'location', 
    # 'fields', 'sort_by', 'data', 'total_objects', 'shown_objects'
    fields = obj.get("fields") or error(f"can't get fields list")
    data   = obj.get("data")   or error(f"can't get data list")
    ic(fields)
    # 'Designation'             object name
    # 'Full name'               full name
    # 'Rise time'               HH:MM rise time > min-elev
    # 'Transit time'            HH:MM transit time
    # 'Set time'                HH:MM set time < min elev
    # 'Max. time observable'    HH:MM observable time    
    # 'R.A.'                    RA (deviates from ephemeris!)
    # 'Dec.'                    DEC (")
    # 'Vmag'                    V magnitude
    # 'Helio. range (au)'
    # 'Topo.range (au)'
    # 'Object-Observer-Sun (deg)'
    # 'Object-Observer-Moon (deg)'
    # 'Galactic latitude (deg)'
    objects = {}
    for d in data:
        # Convert to dictionary
        object1 = { field: d[idx] for idx, field in enumerate(fields) }
        designation = object1.get("Designation")
        ic(designation, object1)
        objects[designation] = object1
    return objects



def to_string(obj1: dict, obj2: dict) -> str:
    return f"{obj1.get("Designation"):11s} {obj1.get("Rise time"):6s} {obj1.get("Transit time"):6s} {obj1.get("Set time"):6s}  {float(obj1.get("Vmag")):4.1f}  {obj2.get("Uncertainty")}  {obj2.get("Last OBS")}"

def to_string1(obj1: dict) -> str:
    return f"{obj1.get("Designation"):11s} {obj1.get("Rise time"):6s} {obj1.get("Transit time"):6s} {obj1.get("Set time"):6s}  {obj1.get("Vmag")}"



def mpc_query_customize(url: str, filename: str, list_type: str) -> None:
    ic(url, filename)

    timeout = config.requests_timeout

    data = { 
        "ra1":          "00 00",                    # RA limits
        "ra2":          "24 00",
        "dec1":         f"{config.min_dec:+02d} 00",# DEC limits
        "dec2":         f"{config.max_dec:+02d} 00",
        "elong1":       60,                         # Elongation limits
        "elong2":       180,
        "mag1":         0,                          # V mag limits
        "mag2":         config.vmag_max,
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

    if filename:
        with open(filename, mode="w", encoding=response.encoding) as file:
            file.write(response.text)
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



def parse_txt_line(txt_line: str, list_type: str) -> dict:
    ic(txt_line)
##DLU
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
        # ic(designation, type, currently, last_obs, code, mag, uncertainty, arc)

##DLN
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
        # ic(designation, type, currently, last_obs, code, mag, uncertainty, arc)

    else:
        error(f"unknown list type {list_type}")

##Common
    ra, dec, mag1, elongation, motion = currently.split("/")
    # ic(ra, dec, mag1, elongation, motion)
    object  = { "Designation":   designation,
                "Type":          type,
                "RA":            float(ra),
                "DEC":           float(dec),
                "VMag":          float(mag1) if mag1.strip() else "",
                "Elongation":    float(elongation),
                "Motion":        float(motion),
                "Marker":        marker,
                "Last OBS":      last_obs,
                "MPC Code":      code,
                "Last mag":      float(mag) if mag else "",
                "Filter":        filter,
                "Uncertainty":   uncertainty,
                "Arc":           float(arc)
                }
    ic(object)
    return object



def mpc_parse_customize(content: str, filename: str, list_type: str) -> dict:
    if not content:
        with open(filename) as file:
            content = file.read()

    objects = {}

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

        object1 = parse_txt_line(txt_line, list_type)
        if object1:
            objects[ object1.get("Designation") ] = object1

    return objects



def mpc_query_lastobs(url: str, filename: str = None) -> None:
    ic(url, filename)

    timeout = config.requests_timeout

    ic(url)
    verbose(f"query {url}")
    response = requests.get(url, timeout=timeout)
    ic(response.status_code)
    if response.status_code != 200:
        error(f"query to {url} failed")

    if filename:
        with open(filename, mode="w", encoding=response.encoding) as file:
            file.write(response.text)
    return response.text



def mpc_parse_lastobs(content: str, filename: str = None) -> dict:
    if not content:
        with open(filename) as file:
            content = file.read()

    objects = {}
    for txt_line in content.splitlines():
        object1 = parse_txt_line(txt_line, "DLN")       # same format als customize/DLN
        if object1:
            # Check mag limit
            mag = object1.get("Last mag") or 99.99
            if mag > config.vmag_max:
                continue
            # Check DEC limit
            dec = object1.get("DEC") or 0.0
            if dec < config.min_dec or dec > config.max_dec:
                continue
            objects[ object1.get("Designation") ] = object1

    return objects



def object_filter(key: str, obj1: dict, obj2: dict) -> bool:
    ic(key, obj1, obj2)

    if obj1:
        # not yet checking anything here
        pass
    if obj2:
        last_obs = obj2.get("Last OBS")
        if last_obs:
            t_obs = Time(last_obs)
            t_now = Time.now()
            ic(t_obs, t_now)
            ##FIXME: add config option
            if t_obs < t_now - 14 * u.day:
                verbose(f"{key}: last obs {last_obs} too old, not included")
                return False

        uncertainty = obj2.get("Uncertainty")
        if uncertainty:
            ##FIXME: add config option
            if int(uncertainty) < 3:
                verbose(f"{key}: uncertainty < 3, not included")
                return False

    # Default is not filtered
    return True



def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = DESCRIPTION,
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("--asteroids", action="store_true", help=f"get asteroids default={config.sb_kind}")
    arg.add_argument("--neo", action="store_true", help=f"get NEOs default={config.sb_group}")
    arg.add_argument("--pha", action="store_true", help=f"get PHAs")
    arg.add_argument("--comets", action="store_true", help=f"get comets (overrides asteroid options)")
    arg.add_argument("-o", "--output", help="write object list to OUTPUT")
    arg.add_argument("-M", "--mag-limit", help="override mag_limit from config")

    args = arg.parse_args()

    if args.debug:
        ic.enable()
        ic(sys.version_info, sys.path, args)
    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()

    # override defaults from config
    if args.mag_limit:
        config.mag_limit = float(args.mag_limit)
        config.vmag_max  = float(args.mag_limit)
    if args.asteroids:
        config.sb_kind = "a"
    if args.neo:
        config.sb_group = "neo"
    if args.pha:
        config.sb_group = "pha"
    if args.comets:
        config.sb_kind = "c"
        config.sb_group = None

    # get sbwobs objects from JPL
    query = jpl_query_sbwobs(config.sbwobs_url)
    objs1 = jpl_parse_sbwobs(query)
    keys1 = objs1.keys()
    verbose(f"WOBS objects ({len(keys1)}): {", ".join(sorted(keys1))}")

    if args.comets:
        # Comets
        keys_selected = sorted(keys1)
        verbose("---------------------------------------")
        verbose("Designation Rise   Trans  Set     Vmag")
        verbose("---------------------------------------")
        for key in keys_selected:
            verbose(to_string1(objs1.get(key)))
        verbose("---------------------------------------")

    else:
        # Asteroids
        query = mpc_query_customize(config.customize_url, None, "DLU")   # "DLN" | "DLU"
        objs2 = mpc_parse_customize(query, None, "DLU")
        keys2 = objs2.keys()
        verbose(f"DLU objects ({len(keys2)}): {", ".join(sorted(keys2))}")

        # Filter DLU for "Last OBS"
        keys2_filtered = [ k for k in keys2 if object_filter(k, objs1.get(k), objs2.get(k)) ]

        # Intersection 1 & 2: observable objects also in DLU list
        keys_selected = sorted(keys1 & keys2_filtered)
        verbose(f"WOBS & DLU objects ({len(keys_selected)}): {", ".join(keys_selected)}")

        verbose("-----------------------------------------------------")
        verbose("Designation Rise   Trans  Set     Vmag  U  Last Obs")
        verbose("-----------------------------------------------------")
        for key in keys_selected:
            verbose(to_string(objs1.get(key), objs2.get(key)))
        verbose("-----------------------------------------------------")

    # DLN is a subset of DLU, tested on 2025-11-18, not used
    # mpc_query_customize(config.customize_url, "tmp/dln.html", "DLN")   # "DLN" | "DLU"
    # objs3 = mpc_parse_customize(None, "tmp/dln.html", "DLN")
    # keys3 = sorted(objs3.keys())
    # verbose(f"DLN objects ({len(keys3)}): {", ".join(keys3)}")

    # mpc_query_lastobs(config.lastobs_url, "tmp/lastobs.txt")
    # objs4 = mpc_parse_lastobs(None, "tmp/lastobs.txt")
    # keys4 = sorted(objs4.keys())
    # verbose(f"Last obs objects ({len(keys4)}): {", ".join(keys4)}")

    if args.output:
        with open(args.output, "w") as file:
            file.writelines([line + "\n" for line in keys_selected])



if __name__ == "__main__":
    main()
