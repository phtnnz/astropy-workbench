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
# Version 0.0 / 2025-11-22
#       Get ephemeris for solar system objects
# Version 0.1 / 2025-12-23
#       Somewhat usable now ;-)
# Version 0.2 / 2026-01-03
#       Common column names for ephemeris tables

VERSION     = "0.2 / 2026-01-03"
AUTHOR      = "Martin Junius"
NAME        = "sbephem"
DESCRIPTION = "Ephemeris for solar system objects"

import sys
import argparse
import re

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
# from astropy.coordinates import Angle
# from astropy.coordinates import errors
from astropy.time        import Time, TimeDelta
from astropy.table       import Table
# from astropy.units       import Quantity, Magnitude
import astropy.units as u
import numpy as np

from sbpy.data import Ephem
from sbpy.data import Obs
from sbpy.data.core import QueryError

from astroquery.mpc import MPC
from astroquery.exceptions import EmptyResponseError, InvalidQueryError

from astroplan import Observer

# Local modules
from verbose import verbose, warning, error, message
from astroutils import location_to_string, get_location     # ...
from neoutils import Exposure, exposure_from_ephemeris

##FIXME: use config
DEFAULT_LOCATION = "M49"



# Command line options
class Options:
    loc: EarthLocation = None       # -l --location



def rename_columns_mpc(eph: Ephem) -> None:
    """Rename MPC ephemeris table to common column names

    Args:
        table (Table): Ephemeris table
    """
    eph._table.rename_columns(("Date",    "Dec",      "V",             "Proper motion", "Direction", 
                               "Azimuth", "Altitude", "Moon distance", "Moon altitude" ),
                              ("Obstime", "DEC",      "Mag",           "Motion",        "PA",        
                               "Az",      "Alt",      "Moon_dist",     "Moon_alt"      ))


def rename_columns_jpl(eph: Ephem) -> None:
    """Rename JPL ephemeris table to common column names

    Args:
        table (Table): Ephemeris table
    """

    mag_col = "Tmag" if "Tmag" in eph.field_names else "V"
    eph.table.rename_columns(("targetname", "epoch",   "AZ", "EL",  mag_col, "velocityPA"),
                             ("Targetname", "Obstime", "Az", "Alt", "Mag",   "PA" ))



def id_type_from_name(name: str) -> str:
    id_type_regex = {   "asteroid numer":         "^[1-9][0-9]*$",
                        "asteroid designation":   "^20[0-9]{2}[ _][A-Z]{2}[0-9]{0,3}$",
                        "comet number":           "^[0-9]{1,3}[PIA]$",
                        "comet designation":      "^[PDCXAI]/20[0-9]{2}[ _][A-Z]{2}[0-9]{0,3}$"
                    }

    for id in id_type_regex.keys():
        m = re.match(id_type_regex[id], name)
        if m:
            ic(name, id)
            return id
    ## Default None or "asteroid designation"?
    return None





def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = DESCRIPTION,
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("-l", "--location", help=f"coordinates, named location or MPC station code, default {DEFAULT_LOCATION}")
    arg.add_argument("-f", "--file", help="read list of objects from file")
    arg.add_argument("-t", "--time", help="start time for ephemeris (1h, 5min steps)")
    arg.add_argument("-J", "--jpl", action="store_true", help="use JPL Horizons ephemeris, default MPC")
    arg.add_argument("-a", "--allnight", action="store_true", help="ephemeris for midnight +/- 8h (30min steps)")
    arg.add_argument("--obs", action="store_true", help="output MPC obs")
    arg.add_argument("--clear", action="store_true", help="clear MPC cache")
    arg.add_argument("object", nargs="*", help="object name")

    args = arg.parse_args()

    if args.debug:
        ic.enable()
        ic(sys.version_info, sys.path, args)
    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()

    # Observer location
    loc = get_location(args.location if args.location else DEFAULT_LOCATION)
    ic(loc, loc.to_geodetic())
    Options.loc = loc
    verbose(f"location {location_to_string(loc)}")

    observer = Observer(location=loc, description=loc.info.name)
    ic(observer)

    # Observation time
    time = Time(args.time) if args.time else Time.now()
    ic(time)
    verbose(f"time {time.iso} ({time.scale.upper()})")

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
    verbose(f"midnight {midnight.iso} / rounded {midnight1.iso} ({time.scale.upper()})")
    verbose(f"nautical twilight {twilight_evening.iso} / {twilight_morning.iso} ({time.scale.upper()})")

    if args.allnight:
        epochs = {"start":  midnight1 - 8 * u.hour,
                  "step":   30 * u.min,
                  "stop":   midnight1 + 9 * u.hour
                  }
    else:
        epochs = {"start":  time,
                  "step":   5 * u.min,
                  "stop":   time + 1 * u.hour
                  }
    ic(epochs)

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

    for obj in objects:
        verbose(f"object {obj}")

        # # Object ephemeris
        # try:
        #     eph = MPC.get_ephemeris(obj, location=loc, start=time, step=5*u.min, number=12)
        #     print(eph["Date", "RA", "Dec", "V", "Proper motion", "Direction", "Azimuth", "Altitude"])
        # except InvalidQueryError as err:
        #     warning(f"query MPC ephemeris for failed {obj}!")
        #     eph = None

        # # Observations
        # try:
        #     ##FIXME: must specify id_type
        #     obs = MPC.get_observations(obj, id_type="comet number")
        #     print(obs)
        # except (EmptyResponseError, ValueError) as err:
        #     warning(f"query MPC observations failed for {obj}!")
        #     obs = None

        # Get ephemerides via sbpy
        if args.jpl:
            eph = Ephem.from_horizons(obj, location=loc, epochs=epochs)
            # Compute total motion from RA/DEC rates
            eph["Motion"] = np.sqrt( np.square(eph["RA*cos(Dec)_rate"]) + np.square(eph["DEC_rate"]) )
            # Rename columns to common names
            rename_columns_jpl(eph)
            ic(eph.field_names)

            mag = eph["Mag"][0]
            ##FIXME: get min altitude from config
            mask = (eph["Alt"] > 25 * u.deg) & (eph["Obstime"] > twilight_evening) & (eph["Obstime"] < twilight_morning)
            eph1 = eph[mask]
            message.print_lines(eph1["Targetname", "Obstime", "RA", "DEC", "Mag", "Motion", "PA", "Az", "Alt"])

            exp = exposure_from_ephemeris(eph, "Motion", mag)
            message(exp)
        else:
            try:
                # eph = Ephem.from_mpc(obj, location=loc, epochs=epochs)
                eph = Ephem.from_mpc(obj, location=loc, epochs=epochs, 
                                        ra_format={'sep': ':', 'unit': 'hourangle', 'precision': 1}, 
                                        dec_format={'sep': ':', 'precision': 1} )
                # Rename columns to common names
                rename_columns_mpc(eph)
                ic(eph.field_names)

                mag = eph["Mag"][0]
                ##FIXME: get min altitude from config
                mask = (eph["Alt"] > 25 * u.deg) & (eph["Obstime"] > twilight_evening) & (eph["Obstime"] < twilight_morning)
                eph1 = eph[mask]
                message.print_lines(eph1["Targetname", "Obstime", "RA", "DEC", "Mag", 
                                         "Motion", "PA", "Az", "Alt", "Moon_dist", "Moon_alt"])

                exp = exposure_from_ephemeris(eph, "Motion", mag)
                message(exp)
            except QueryError as e:
                warning(f"MPC ephemeris for {obj} failed")

            try:
                if args.obs:
                    obs = Obs.from_mpc(obj, id_type=id_type_from_name(obj))
                    print(obs)
                    # # Handle masked entries
                    # for i in range(-1, -10, -1):
                    #     mag = obs["mag"][i].unmasked
                    #     if mag > Magnitude(0):
                    #         break
            except QueryError as e:
                warning(f"MPC observations for {obj} failed")
            except ConnectionError as e:
                warning(f"MPC request failed: {e}")



if __name__ == "__main__":
    main()
