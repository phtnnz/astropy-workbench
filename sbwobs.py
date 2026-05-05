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
# Version 0.4 / 2026-02-17
#       New option -l --location, some refactoring

VERSION     = "0.4 / 2026-02-17"
AUTHOR      = "Martin Junius"
NAME        = "sbwobs"
DESCRIPTION = "Retrieve observable NEOs/comets from JPL/MPC"

import sys
import argparse
import requests
import json
import re

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.coordinates import Angle
from astropy.time import Time
import astropy.units as u
from astropy.units import Quantity, Magnitude

# Local modules
from verbose import verbose, warning, error, message
from neoconfig import config
from neoutils import fmt_time
from neoephem import get_local_circumstances
from neoclasses import LocalCircumstances, JPLWObsData, MPCDLxData, EphemData



def jpl_query_sbwobs(url: str, local: LocalCircumstances) -> str:
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
    ic(url, local)

    timeout = config.requests_timeout
    start = local.naut_dusk
    end   = local.naut_dawn
    ##FIXME: use lon/lat/alt if MPC code is None

    data = { 
        "mpc-code":     local.code,             # MPC observatory code
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



def jpl_parse_sbwobs(text: str) -> dict[str, EphemData]:
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
        obj1 = { field: d[idx] for idx, field in enumerate(fields) }
        designation = obj1.get("Designation")
        ic(designation, obj1)
        vmag = obj1.get("Vmag")
        if vmag[-1].isupper():
            vmag = vmag[:-1]
        wobs = JPLWObsData(
            obj1.get('Designation'),
            obj1.get('Full name'),
            obj1.get('Rise time'),
            obj1.get('Transit time'),
            obj1.get('Set time'),
            obj1.get('Max. time observable'),
            Angle(obj1.get('R.A.'), unit=u.hourangle),
            Angle(obj1.get('Dec.').replace("\'", " ").replace("\"", " "), unit=u.deg),
            Magnitude(vmag),
            float(obj1.get('Helio. range (au)')) * u.au,
            float(obj1.get('Topo.range (au)')) * u.au,
            float(obj1.get('Object-Observer-Sun (deg)')) * u.deg,
            float(obj1.get('Object-Observer-Moon (deg)')) * u.deg,
            float(obj1.get('Galactic latitude (deg)')) * u.deg,
            # Type from config
            config.sb_group if config.sb_group else "comet" if config.sb_kind == "c" else "-"
        )
        edata = EphemData(wobs.type, wobs.designation, None, None, None, None, wobs.vmag, None, wobs=wobs)
        ic(edata)
        objects[designation] = edata
    return objects



def to_string(edata: EphemData) -> str:
    wobs = edata.wobs
    dlx = edata.dlx
    if dlx:
        return f"{wobs.type.upper():5s}  {wobs.designation:11s} {wobs.rise_time:6s} {wobs.transit_time:6s} {wobs.set_time:6s}  {float(wobs.vmag.value):4.1f}  {dlx.uncertainty}  {str(dlx.last_obs):10.10s}"
    else:
        return str(wobs)



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
        dlx: MPCDLxData = parse_txt_line(txt_line, list_type)
        if dlx:
            edata = EphemData("", dlx.designation, None, None, None, None, dlx.vmag, None, dlx=dlx)
            ic(edata)
            objects[dlx.designation] = edata

    return objects



## Currently not used ##
def mpc_query_lastobs(url: str) -> str:
    ic(url)

    timeout = config.requests_timeout

    ic(url)
    verbose(f"query {url}")
    response = requests.get(url, timeout=timeout)
    ic(response.status_code)
    if response.status_code != 200:
        error(f"query to {url} failed")




## Currently not used ##
def mpc_parse_lastobs(content: str) -> dict[str, MPCDLxData]:
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



def object_filter(key: str, edata: EphemData) -> bool:
    wobs = edata.wobs
    dlx = edata.dlx
    ic(key, wobs, dlx)

    if wobs:
        # not yet checking anything here
        pass
    if dlx:
        last_obs = dlx.last_obs
        if last_obs:
            t_obs = last_obs
            t_now = Time.now()
            ic(t_obs, t_now)
            ##FIXME: add config option
            if t_obs < t_now - 14 * u.day:
                verbose(f"{key}: last obs {last_obs} too old, not included")
                return False

        uncertainty = dlx.uncertainty
        if uncertainty:
            ##FIXME: add config option
            if int(uncertainty) < 3:
                verbose(f"{key}: uncertainty < 3, not included")
                return False

    # Default is not filtered
    return True



def sbwobs_get_obj_edata(local: LocalCircumstances) -> dict[str, EphemData]:
    # get sbwobs objects from JPL
    query = jpl_query_sbwobs(config.sbwobs_url, local)
    obj_edata1 = jpl_parse_sbwobs(query)
    keys1 = obj_edata1.keys()
    verbose(f"WOBS objects ({len(keys1)}): {", ".join(sorted(keys1))}")

    if config.sb_kind == "c":
        # Comets
        keys_selected = sorted(keys1)
        obj_edata = obj_edata1

        verbose("---------------------------------------------")
        verbose("Type   Designation Rise   Trans  Set     Vmag")
        verbose("---------------------------------------------")
        for key in keys_selected:
            verbose(to_string(obj_edata1.get(key)))
        verbose("---------------------------------------------")
    else:
        # Asteroids
        query = mpc_query_customize(config.customize_url, "DLU")   # "DLN" | "DLU"
        obj_edata2 = mpc_parse_customize(query, "DLU")
        keys2 = obj_edata2.keys()
        verbose(f"DLU objects ({len(keys2)}): {", ".join(sorted(keys2))}")

        ## DLN is a subset of DLU, tested on 2025-11-18, not used ##
        # mpc_query_customize(config.customize_url, "DLN")   # "DLN" | "DLU"
        # mpc_parse_customize(content, "DLN")

        ## LastObs list currently not used ##
        # mpc_query_lastobs(config.lastobs_url)
        # mpc_parse_lastobs(content)

        # Filter DLU for "Last OBS"
        # None: data from WOBS not yet used in object_filter()
        keys2_filtered = [ k for k in keys2 if object_filter(k, obj_edata2.get(k)) ]

        # Intersection 1 & 2: observable objects also in DLU list
        keys_selected = sorted(keys1 & keys2_filtered)
        verbose(f"WOBS & DLU objects ({len(keys_selected)}): {", ".join(keys_selected)}")
        # Build new dict
        obj_edata = dict()
        for k in keys_selected:
            edata = obj_edata1.get(k)
            # Copy MPCDLxData from DLU data
            edata.dlx = obj_edata2.get(k).dlx
            obj_edata[k] = edata

        verbose("------------------------------------------------------------")
        verbose("Type   Designation Rise   Trans  Set     Vmag  U  Last Obs")
        verbose("------------------------------------------------------------")
        for key in keys_selected:
            verbose(to_string(obj_edata.get(key)))
        verbose("------------------------------------------------------------")

    return obj_edata



def sbwobs_get_objects(local: LocalCircumstances) -> list[str]:
    # wrapper for sbwobs_get_obj_edata()
        return sbwobs_get_obj_edata(local).keys()



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
    arg.add_argument("-l", "--location", help=f"coordinates, named location or MPC station code, default {config.mpc_code}")

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

    # Observer location and local circumstances
    local = get_local_circumstances(args.location if args.location else config.mpc_code)

    keys_selected = sbwobs_get_objects(local)

    if args.output:
        with open(args.output, "w") as file:
            file.writelines([line + "\n" for line in keys_selected])



if __name__ == "__main__":
    main()
