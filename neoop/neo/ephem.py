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
# Version 0.1 / 2026-01-04
#       Copy of sbephem 0.2

VERSION     = "0.1 / 2026-01-04"
AUTHOR      = "Martin Junius"
NAME        = "neoephem"
DESCRIPTION = "Ephemeris for solar system objects"

import sys
import argparse
import re

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.coordinates import Angle
from astropy.time import Time
import astropy.units as u
import numpy as np

from sbpy.data import Ephem
from sbpy.data.core import QueryError, FieldError

from astroquery.mpc import MPC

from astroplan import Observer

# Local modules
from verbose import verbose, warning, error, message
from astroutils import get_location
from neoclasses import Exposure, EphemTimes, EphemData, EphemDataList, LocalCircumstances
from neoutils import exposure_calc, max_motion, get_mag0, motion_limit
from neoconfig import config

DEFAULT_LOCATION = config.code



def _rename_columns_mpc(eph: Ephem) -> None:
    """Rename MPC ephemeris table to common column names

    Args:
        eph (Ephem): Ephemeris table
    """
    eph.table.rename_columns(("Date",    "Dec",      "V",             "Proper motion", "Direction", 
                               "Azimuth", "Altitude", "Moon distance", "Moon altitude" ),
                              # -->
                              ("Obstime", "DEC",      "Mag",           "Motion",        "PA",        
                               "Az",      "Alt",      "Moon_dist",     "Moon_alt"      ))


def _rename_columns_jpl(eph: Ephem) -> None:
    """Rename JPL ephemeris table to common column names

    Args:
        eph (Ephem): Ephemeris table
    """

    mag_col = "Tmag" if "Tmag" in eph.field_names else "V"
    eph.table.rename_columns(("targetname", "epoch",   "AZ", "EL",  mag_col, "velocityPA"),
                             # -->
                             ("Targetname", "Obstime", "Az", "Alt", "Mag",   "PA" ))



def edata_add_ephem_mpc(edata: EphemData, local: LocalCircumstances) -> None:
    if edata.ephem:
        verbose(f"already got ephemeris for {edata.obj}")
        return

    min_alt = config.min_alt

    obj = edata.obj
    try:
        verbose(f"{obj} ephemeris from MPC")
        ic(local.epochs)
        eph = Ephem.from_mpc(obj, location=local.loc, epochs=local.epochs, 
                                ra_format={'sep': ':', 'unit': 'hourangle', 'precision': 1}, 
                                dec_format={'sep': ':', 'precision': 1} )
        # Rename columns to common names
        _rename_columns_mpc(eph)
        ic(eph.field_names)

        mask = (eph["Alt"] > min_alt * u.deg) & (eph["Obstime"] >= local.naut_dusk) & (eph["Obstime"] <= local.naut_dawn)
        eph1 = eph[mask]
        if len(eph1) == 0:
            warning(f"skipping empty ephemeris for {obj}")
            return
        mag = get_mag0(eph1)
        motion = max_motion(eph1)

        # Copy to EphemData
        if edata.wobs:
            edata.type = edata.wobs.type.upper()
        edata.obj = obj
        edata.ephem = eph1
        edata.mag = mag
        edata.motion = motion

    except (QueryError, FieldError) as e:
        warning(f"MPC ephemeris for {obj} failed")



def edata_add_ephem_jpl(edata: EphemData, local: LocalCircumstances) -> None:
    min_alt = config.min_alt

    obj = edata.obj
    try:
        verbose(f"{obj} ephemeris from JPL")
        eph = Ephem.from_horizons(obj, location=local.loc, epochs=local.epochs)
        # Compute total motion from RA/DEC rates
        eph["Motion"] = np.sqrt( np.square(eph["RA*cos(Dec)_rate"]) + np.square(eph["DEC_rate"]) )
        # Rename columns to common names
        _rename_columns_jpl(eph)
        ic(eph.field_names)

        mask = (eph["Alt"] > min_alt * u.deg) & (eph["Obstime"] >= local.naut_dusk) & (eph["Obstime"] <= local.naut_dawn)
        eph1 = eph[mask]
        if len(eph1) == 0:
            warning(f"skipping empty ephemeris for {obj}")
            return
        ##Quick hack: missing Moon_dist, Moon_alt in JPL ephemeris?
        eph1["Moon_dist"] = 180 * u.degree
        eph1["Moon_alt"]  = -90 * u.degree
        mag = get_mag0(eph1)
        motion = max_motion(eph1)

        # Copy to EphemData
        if edata.wobs:
            edata.type = edata.wobs.type.upper()
        edata.obj = obj
        edata.ephem = eph1
        edata.mag = mag
        edata.motion = motion

    except (QueryError, FieldError) as e:
        warning(f"JPL ephemeris for {obj} failed")



def edata_add_exposure(edata: EphemData, local: LocalCircumstances) -> None:
    obj = edata.obj
    if not edata.ephem:
        return
    if edata.motion != None and edata.mag != None:
        edata.exposure = exposure_calc(edata.motion, edata.mag)
        ic(edata.obj, edata.exposure)
    if not edata.exposure:
        warning(f"exposure calculation for {obj} failed, too fast?")
        if edata.motion != None:
            warning(f"motion={edata.motion:.2f}, limit={motion_limit():.2f}")



def get_local_circumstances(location: str) -> LocalCircumstances:
    """Get local circumentances: location, observer, dusk, dawn, epochs parameter

    Args:
        location (str): observer location

    Returns:
        LocalCircumstances: local circumstances data
    """
    loc = get_location(location)
    ic(loc, loc.to_geodetic())
    # MPC station code
    m = re.search(r'^([0-9A-Z]{3})$', location)
    if m:
        code = m.group(1)
    else:
        code = None
    ic(code)
 
    # Astroplan
    observer = Observer(location=loc, description=loc.info.name)
    ic(observer)

    # Observation times for upcoming night
    time = Time.now()
    ic(time)

    midnight = observer.midnight(time, which="next")
    twilight_evening = observer.twilight_evening_nautical(time, which="next")
    twilight_morning = observer.twilight_morning_nautical(time, which="next")
    if twilight_evening > twilight_morning:
        twilight_evening = observer.twilight_evening_nautical(time, which="previous")
    ic(midnight.iso, twilight_evening.iso, twilight_morning.iso)

    # Round midnight time to nearest 30 min
    rem, day = np.modf(midnight.jd)
    n_round = 24 * 2    # 24 h / 30 min
    rem = round(rem*n_round) / n_round
    jd1 = day + rem
    midnight1 = Time(jd1, format="jd")
    ic(day, rem, midnight1.iso)
    epochs = {"start":  midnight1 - 8 * u.hour,
              "step":   30 * u.min,
              "stop":   midnight1 + 9 * u.hour
             }
    ic(epochs)

    return LocalCircumstances(loc, observer, twilight_evening, twilight_morning, epochs, code)



def get_dec_limits(local: LocalCircumstances, min_alt: Angle) -> tuple[Angle, Angle]:
    """Compute min/max DEC from min altitude and latitude

    Parameters
    ----------
    local : LocalCircumstances
        Observer location related data
    min_alt : Angle
        Minimum altitude

    Returns
    -------
    tuple[Angle, Angle]
        Minimum DEC, maximum DEC
    """
    lat = local.loc.lat
    min_dec = -90*u.deg + min_alt + lat
    max_dec = +90*u.deg - min_alt + lat
    if min_dec < -90*u.deg:
        min_dec = -90*u.deg 
    if max_dec > +90*u.deg:
        max_dec = +90*u.deg 
    ic(lat, min_alt, min_dec, max_dec)
    return Angle(min_dec), Angle(max_dec)
