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
from astropy.coordinates import AltAz, EarthLocation
from astropy.time        import Time
import astropy.units as u
import numpy as np

from sbpy.data import Ephem
from sbpy.data.core import QueryError

from astroquery.mpc import MPC

from astroplan import Observer

# Local modules
from verbose import verbose, warning, error, message
from astroutils import get_location
from neoclasses import Exposure, EphemTimes, EphemData, LocalCircumstances
from neoutils import exposure_calc, max_motion, get_mag0
from neoconfig import config

DEFAULT_LOCATION = config.code



def _rename_columns_mpc(eph: Ephem) -> None:
    """Rename MPC ephemeris table to common column names

    Args:
        eph (Ephem): Ephemeris table
    """
    eph._table.rename_columns(("Date",    "Dec",      "V",             "Proper motion", "Direction", 
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



def get_ephem_mpc(objects: list, local: LocalCircumstances) -> dict[str, EphemData]:
    """Get ephemerides for object list from MPC

    Args:
        objects (list): object list
        local (LocalCircumstances): local circumstances for upcoming observation night

    Returns:
        dict[EphemData]: object dict with ephemerides data
    """
    min_alt = config.min_alt
    obj_data = {}

    for obj in objects:
        try:
            verbose(f"{obj} ephemeris from MPC")
            # eph = Ephem.from_mpc(obj, location=loc, epochs=epochs)
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
                continue
            mag = get_mag0(eph1)
            motion = max_motion(eph1)
            exp = exposure_calc(motion, mag)
            data = EphemData(obj, None, eph1, None, exp, mag, motion)
            obj_data[obj] = data
        except QueryError as e:
            warning(f"MPC ephemeris for {obj} failed")

    return obj_data



def get_ephem_jpl(objects: list, local: LocalCircumstances) -> dict[str, EphemData]:
    """Get ephemerides for object list from JPL

    Args:
        objects (list): object list
        local (LocalCircumstances): local circumstances for upcoming observation night

    Returns:
        dict[EphemData]: object dict with ephemerides data
    """
    min_alt = config.min_alt
    obj_data = {}

    for obj in objects:
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
                continue
            ##Quick hack: missing Moon_dist, Moon_alt in JPL ephemeris?
            eph1["Moon_dist"] = 180 * u.degree
            eph1["Moon_alt"]  = -90 * u.degree
            mag = get_mag0(eph1)
            motion = max_motion(eph1)
            exp = exposure_calc(motion, mag)
            data = EphemData(obj, None, eph1, None, exp, mag, motion)
            obj_data[obj] = data
        except QueryError as e:
            warning(f"MPC ephemeris for {obj} failed")

    return obj_data



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



def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = DESCRIPTION,
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("-l", "--location", help=f"coordinates, named location or MPC station code, default {DEFAULT_LOCATION}")
    arg.add_argument("-f", "--file", help="read list of objects from file")
    arg.add_argument("-J", "--jpl", action="store_true", help="use JPL Horizons ephemeris, default MPC")
    arg.add_argument("--clear", action="store_true", help="clear MPC cache")
    arg.add_argument("object", nargs="*", help="object name")

    args = arg.parse_args()

    if args.debug:
        ic.enable()
        ic(sys.version_info, sys.path, args)
    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()

    # Observer location and local circumstances
    local = get_local_circumstances(args.location if args.location else DEFAULT_LOCATION)

    # Objects
    objects = []
    if args.file:
        with open(args.file, "r") as file:
            for line in file:
                objects.append(line.strip())
    if args.object:
        objects.extend(args.object)
    ic(objects)
    if not objects:
        error("no objects from file or command line")

    # Clear astroquery cache
    if args.clear:
        MPC.clear_cache()    

    verbose.print_lines(local)

    if args.jpl:
        obj_data = get_ephem_jpl(objects, local)
    else:
        obj_data = get_ephem_mpc(objects, local)

    ic(obj_data)
    for obj, data in obj_data.items():
        verbose.print_lines(data.ephem["Targetname", "Obstime", "RA", "DEC", "Mag", 
                                       "Motion", "PA", "Az", "Alt", "Moon_dist", "Moon_alt"])
        verbose(data.exposure)



if __name__ == "__main__":
    main()
