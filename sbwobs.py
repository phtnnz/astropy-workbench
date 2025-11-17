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
# Version 0.0 / 2025-11-14
#       Retrieve "What's Observable?" from JPL

VERSION     = "0.0 / 2025-11-14"
AUTHOR      = "Martin Junius"
NAME        = "sbwobs"
DESCRIPTION = "Retrieve JPL What's Observable"

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
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.coordinates import Angle
from astropy.coordinates import errors
from astropy.time        import Time, TimeDelta
import astropy.units as u
import astropy.constants as const
import numpy as np

# Local modules
from verbose import verbose, warning, error, message
from astroutils import location_to_string, get_location     # ...
from jsonconfig import JSONConfig, config

DEFAULT_LOCATION = "M49"

CONFIGFILE = "sbwobs.json"
config = JSONConfig(CONFIGFILE)
config.set_error_on_missing()



# Command line options
class Options:
    pass



def jpl_query_sbwobs(url: str, filename: str) -> None:
    """
    Retrieve What's Observable from JPL

    Parameters
    ----------
    url : str
        JPL API URL
    filename : str
        Local file name in cache
    """
    ic(url, filename)

    timeout = config.requests_timeout
    ##FIXME: get from command line, twilight times from astroplan start/end -> obs-time/obs-end
    time_obs = Time.now().strftime("%Y-%m-%d")

    data = { 
        "mpc-code":     config.mpc_code,    # MPC observatory code
        "obs-time":     time_obs,           # Date/time of the observation
        "elev-min":     config.elev_min,    # Minimum altitude
        "vmag-max":     config.vmag_max,    # Max V mag = minimum brightness
        "output-sort":  config.output_sort, # Sort records

        "sb-ns":        config.sb_ns,       # Numbered (n) ./. unnumbered (u)
        "sb-kind":      config.sb_kind,     # Asteroids (a) ./. comets (c)
        "sb-group":     config.sb_group     # NEOs (neo) ./. PHA (pha)
    }

    ic(url, data)
    response = requests.get(url, params=data, timeout=timeout)
    ic(response.status_code)
    if response.status_code != 200:
        error(f"query to {url} failed")

    with open(filename, mode="w", encoding=response.encoding) as file:
        file.write(response.text)



def jpl_parse_sbwobs(text: str, filename: str) -> dict:
    """
    Parse JSON retrieve from JPL What's Observable

    Parameters
    ----------
    text : str
        String to load JSON from
    filename : str
        File name to load JSON from

    Returns
    -------
    dict
        Parsed dictionary, key is object designation
    """
    if text:
        obj = json.loads(text)
    else:
        with open(filename) as file:
            obj = json.load(file)
    
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
        "wh":           list_type
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
    response = requests.get(url, params=data, timeout=timeout)
    ic(response.status_code)
    if response.status_code != 200:
        error(f"query to {url} failed")

    with open(filename, mode="w", encoding=response.encoding) as file:
        file.write(response.text)



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

##FIXME: move the following part of the parsing to a separate function
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
            ic(txt_line)
            designation = txt_line[ 9:27].strip()
            type        = txt_line[29:32]
            currently   = txt_line[34:58].strip()
            marker      = txt_line[59:60].strip()
            last_obs    = txt_line[60:72]
            code        = txt_line[74:77]
            mag         = txt_line[79:84].strip()
            filter      = txt_line[84:85].strip()
            uncertainty = txt_line[87:88]
            arc         = txt_line[90:93].strip()
            # ic(designation, type, currently, last_obs, code, mag, uncertainty, arc)

##DLN
#         2025 WA             Amo   1.3/+16/18.3 /147/07.42  2025 Nov. 16  C23  17.8 c  7    1   6.4/+01/18.3/141/08.40   9.1/-08/22.2/121/00.19   8.9/-04/23.3/151/00.30
#         ^9                  ^29  ^34                   58^ ^60           ^74  ^79     ^87^90  ^95
        if list_type == "DLN":
            # Fix one char shift if "Motn" is 3 digits
            if txt_line[58:60] != "  ":
                if txt_line[59:61] == "  ":
                    txt_line = txt_line[:59] + txt_line[60:]
                else:
                    warning("bad format in data line")
                    warning(txt_line)
            ic(line)
            designation = txt_line[ 9:27].strip()
            type        = txt_line[29:32]
            currently   = txt_line[34:59].strip()
            marker      = ""
            last_obs    = txt_line[60:72]
            code        = txt_line[74:77]
            mag         = txt_line[79:84].strip()
            filter      = txt_line[84:85].strip()
            uncertainty = txt_line[87:88]
            arc         = txt_line[90:93].strip()
            # ic(designation, type, currently, last_obs, code, mag, uncertainty, arc)

##Common
        ra, dec, mag1, elongation, motion = currently.split("/")
        # ic(ra, dec, mag1, elongation, motion)
        object1 = {"Designation":   designation,
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
        ic(object1)
        objects[designation] = object1



def mpc_query_lastobs(url: str, filename: str) -> None:
    ic(url, filename)

    timeout = config.requests_timeout

    ic(url)
    response = requests.get(url, timeout=timeout)
    ic(response.status_code)
    if response.status_code != 200:
        error(f"query to {url} failed")

    with open(filename, mode="w", encoding=response.encoding) as file:
        file.write(response.text)



def mpc_parse_lastobs(content: str, filename: str) -> dict:
    if not content:
        with open(filename) as file:
            content = file.read()

#  Designation                         Currently              Last obs.  Code   Mag.   U  Arc      Currently + 7d           Currently + 30d          Currently + 60 d
#                                   R.A./Dec/  V /El./Motn                                       R.A./Dec/  V /El./Motn   R.A./Dec/  V /El./Motn   R.A./Dec/  V /El./Motn
# Sample:
#         2025 WC             Apo   0.4/-47/16.6 /100/421.97  2025 Nov. 17  H21  18.0 G  5    1  16.2/+03/31.8/024/00.12  16.2/+01/32.3/032/00.26  16.9/-08/31.2/046/00.63
#         2025 WA             Amo   1.2/+16/18.3 /147/07.21  2025 Nov. 16  C23  17.8 c  7    1   6.3/+02/18.2/142/08.64   9.1/-08/22.2/121/00.19   8.9/-04/23.3/151/00.30
#         ^9                  ^29  ^34                   58^ ^60           ^74  ^79     ^87^90  ^95

    objects = {}
    for line in content.splitlines():
        # Fix one char shift if "Motn" is 3 digits
        if line[58:60] != "  ":
            if line[59:61] == "  ":
                line = line[:59] + line[60:]
            else:
                warning("bad format in data line")
                warning(line)
        ic(line)
        designation = line[ 9:27].strip()
        type        = line[29:32]
        currently   = line[34:59].strip()
        last_obs    = line[60:72]
        code        = line[74:77]
        mag         = line[79:84].strip()
        filter      = line[84:85].strip()
        uncertainty = line[87:88]
        arc         = line[90:93].strip()
        # ic(designation, type, currently, last_obs, code, mag, uncertainty, arc)
        ra, dec, mag1, elongation, motion = currently.split("/")
        # ic(ra, dec, mag1, elongation, motion)
        object1 = {"Designation":   designation,
                   "Type":          type,
                   "RA":            float(ra),
                   "DEC":           float(dec),
                   "VMag":          float(mag1),
                   "Elongation":    float(elongation),
                   "Motion":        float(motion),
                   "Last OBS":      last_obs,
                   "MPC Code":      code,
                   "Last mag":      float(mag) if mag else "",
                   "Filter":        filter,
                   "Uncertainty":   uncertainty,
                   "Arc":           float(arc)
                   }
        ic(object1)
        objects[designation] = object1



def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = DESCRIPTION,
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")

    args = arg.parse_args()

    if args.debug:
        ic.enable()
        ic(sys.version_info, sys.path, args)
    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()

    ##TEST##
    # jpl_query_sbwobs(config.sbwobs_url, "tmp/sbwobs.json")
    # jpl_parse_sbwobs(None, "tmp/sbwobs.json")

    # mpc_query_customize(config.customize_url, "tmp/sbwobs.html", "DLU")   # "DLN" | "DLU"
    # mpc_parse_customize(None, "tmp/sbwobs.html", "DLU")

    mpc_query_customize(config.customize_url, "tmp/sbwobs.html", "DLN")   # "DLN" | "DLU"
    mpc_parse_customize(None, "tmp/sbwobs.html", "DLN")

    # mpc_query_lastobs(config.lastobs_url, "tmp/sbwobs.txt")
    # mpc_parse_lastobs(None, "tmp/sbwobs.txt")



if __name__ == "__main__":
    main()
