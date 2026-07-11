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
# Version 1.2 / 2026-06-26
#       Get last obs for comets, using mpc.observations, filter accordingly
# Version 1.3 / 2026-07-11
#       Provides from_sbwobs class method for EphemDataList

# Usage
#       import jpl.sbwobs

VERSION     = "1.3 / 2026-07-11"
AUTHOR      = "Martin Junius"
NAME        = "jpl.sbwobs"
DESCRIPTION = "Retrieve observable NEOs/comets from JPL"

import requests
import json
from time import sleep

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
from astro.utils import fmt_time
from neo.classes import LocalCircumstances, WObsData, EphemData, EphemDataList, EphemDataDict, Obs
from mpc.lastobs import mpc_query_customize, mpc_parse_customize, mpc_query_lastobs, mpc_parse_lastobs



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



def jpl_parse_sbwobs(text: str) -> EphemDataDict:
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
    objects = EphemDataDict()
    for d in data:
        # Convert to dictionary
        obj1 = { field: d[idx] for idx, field in enumerate(fields) }
        designation = obj1.get("Designation")
        ic(designation, obj1)
        vmag = obj1.get("Vmag")
        if vmag[-1].isupper():
            vmag = vmag[:-1]
        wobs = WObsData(
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
        last_obs = wobs.last_obs
        if last_obs:
            t_obs = last_obs
            t_now = Time.now()
            ic(t_obs, t_now)
            if t_obs < t_now - config.max_last_obs * u.day:
                verbose(f"{wobs.designation}: last obs {last_obs.iso} too old, not included")
                return False

    if dlx:
        last_obs = dlx.last_obs
        if last_obs:
            t_obs = last_obs
            t_now = Time.now()
            ic(t_obs, t_now)
            if t_obs < t_now - config.max_last_obs * u.day:
                verbose(f"{key}: last obs {last_obs.iso} too old, not included")
                return False

        uncertainty = dlx.uncertainty
        if uncertainty:
            min_uncertainty = config.min_uncertainty
            if int(uncertainty) < min_uncertainty:
                verbose(f"{key}: uncertainty < 3, not included")
                return False

    # Default is not filtered
    return True



def edata_comet_add_prefix(edata: EphemData) -> EphemData:
    if not edata.wobs or edata.type != "comet":
        return None
    if edata.wobs.full_name[1] == "/":
        # Comet full name contains X/ prefix
        edata.obj = edata.wobs.designation = f"{edata.wobs.full_name[0]}/{edata.obj}"
    return edata



def edata_add_last_obs(edata: EphemData, local: LocalCircumstances) -> EphemData:
    if not edata.wobs:
        return None
    obs = Obs.from_object(edata.obj)
    edata.wobs.last_obs = obs.get_last_obs()
    verbose(f"{edata.obj}: last obs {edata.wobs.last_obs.iso}")
    ic(edata.obj, edata.wobs.last_obs.iso)
    return edata



@classmethod
def edata_list_from_sbwobs(cls, local: LocalCircumstances, list_type: str="DLU") -> EphemDataList:
    # get sbwobs objects from JPL
    query = jpl_query_sbwobs(config.sbwobs_url, local)
    obj_edata1: EphemDataDict = jpl_parse_sbwobs(query)
    keys1 = obj_edata1.keys()
    verbose(f"WOBS {obj_edata1}")

    if config.sb_kind == "c":
        # Comets
        for key, edata in obj_edata1.items():
            edata = obj_edata1.get(key)
            edata_comet_add_prefix(edata)
            edata_add_last_obs(edata, local)
            ##FIXME: config
            sleep(0.25) # avoid rapid fire to MPC servers

        # Filter list
        keys1_filtered = [ k for k in keys1 if object_filter(k, obj_edata1.get(k)) ]
        keys_selected = sorted(keys1_filtered)

        # Build new dict
        obj_edata = dict()
        for k in keys_selected:
            obj_edata[k] = obj_edata1.get(k)

        verbose("----------------------------------------------------------")
        verbose("Type   Designation  Rise   Trans  Set     Vmag  Last Obs")
        verbose("----------------------------------------------------------")
        for key, edata in obj_edata.items():
            verbose(edata)
        verbose("----------------------------------------------------------")
    else:
        # Asteroids
        if list_type == "LASTOBS":
            # LastObs
            query = mpc_query_lastobs(config.lastobs_url)
            obj_edata2 = mpc_parse_lastobs(query)
        else:
            # DLU/DLN
            query = mpc_query_customize(config.customize_url, "DLU")   # "DLN" | "DLU"
            obj_edata2 = mpc_parse_customize(query, "DLU")

        keys2 = obj_edata2.keys()
        verbose(f"{list_type} {obj_edata2}")

        # Filter MPC list
        keys2_filtered = [ k for k in keys2 if object_filter(k, obj_edata2.get(k)) ]

        # Intersection 1 & 2: observable objects also in MPC list
        keys_selected = sorted(keys1 & keys2_filtered)
        # Build new dict
        obj_edata = EphemDataDict()
        for k in keys_selected:
            edata = obj_edata1.get(k)
            # Copy MPCDLxData from MPC data
            edata.dlx = obj_edata2.get(k).dlx
            obj_edata[k] = edata
        verbose(f"WOBS & {list_type} {obj_edata}")

        verbose("------------------------------------------------------------")
        verbose("Type   Designation Rise   Trans  Set     Vmag  U  Last Obs")
        verbose("------------------------------------------------------------")
        for key in keys_selected:
            verbose(obj_edata.get(key))
        verbose("------------------------------------------------------------")

    return EphemDataList.from_dict(obj_edata)

EphemDataList.from_sbwobs = edata_list_from_sbwobs


def sbwobs_get_objects(local: LocalCircumstances, list_type: str="DLU") -> list[str]:
    # wrapper for EphemDataList.from_sbwobs()
    edata_list = EphemDataList.from_sbwobs(local, list_type)
    return edata_list.objects()
