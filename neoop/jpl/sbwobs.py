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
# Version 1.0 / 2026-06-16
#       Moved and adapted to new directory structure under neoop/
# Version 1.1 / 2026-06-22
#       Moved MPC functions to new module mpc.lastobs

VERSION     = "1.1 / 2026-06-22"
AUTHOR      = "Martin Junius"
NAME        = "jpl.sbwobs"
DESCRIPTION = "Retrieve observable NEOs/comets from JPL"

import requests
import json

from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.coordinates import Angle
from astropy.time import Time
import astropy.units as u
from astropy.units import Quantity, Magnitude

# Local modules
from utils.verbose import verbose, warning, error, message
from neo.config import config
from neo.utils import fmt_time
from neo.classes import LocalCircumstances, JPLWObsData, MPCDLxData, EphemData, EphemDataList
from mpc.lastobs import mpc_query_customize, mpc_parse_customize



def jpl_query_verbose_error(text: str) -> None:
    obj = json.loads(text)
    verbose(f"ERROR: {obj['code']} {obj['message']}")
    verbose(f"ERROR: info {obj['moreInfo']}")



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
        "elev-min":     config.min_alt,         # Minimum altitude
        "vmag-max":     config.sbwobs_mag_limit,# Max V mag = minimum brightness
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
        jpl_query_verbose_error(response.text)
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
            if t_obs < t_now - config.max_last_obs * u.day:
                verbose(f"{key}: last obs {last_obs} too old, not included")
                return False

        uncertainty = dlx.uncertainty
        if uncertainty:
            min_uncertainty = config.min_uncertainty
            if int(uncertainty) < min_uncertainty:
                verbose(f"{key}: uncertainty < 3, not included")
                return False

    # Default is not filtered
    return True



def sbwobs_get_edata_list(local: LocalCircumstances) -> EphemDataList:
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
            verbose(obj_edata1.get(key))
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
            verbose(obj_edata.get(key))
        verbose("------------------------------------------------------------")

    return EphemDataList.from_dict(obj_edata)



def sbwobs_get_objects(local: LocalCircumstances) -> list[str]:
    # wrapper for sbwobs_get_obj_edata()
    return sbwobs_get_edata_list(local).objects()



