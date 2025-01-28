#!/usr/bin/env python

# Copyright 2024-2025 Martin Junius
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
# Version 0.1 / 2025-01-27
#       First steps with astroplan

import sys
import argparse
from datetime import datetime, timezone
from zoneinfo import ZoneInfo
# Required on Windows
import tzdata

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, CartesianRepresentation
from astropy.coordinates import ICRS, GCRS, PrecessedGeocentric, FK4, FK5, HADec  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
from astropy.coordinates import errors
import astropy.units as u
from astropy.time        import Time, TimeDelta
import numpy as np

# Astroplan
from astroplan import Observer

# Local modules
from verbose import verbose, warning, error
from astroutils import ra_from_lst_ha, ra_dec_to_string, angle_to_string, hourangle_to_string, get_location

VERSION = "0.1 / 2025-01-27"
AUTHOR  = "Martin Junius"
NAME    = "rise-n-set"



# Command line options
class Options:
    pass



def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = "Rise and set times",
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("-t", "--time", help="time (UTC) for computation, default now")
    arg.add_argument("-l", "--location", help="coordinates, named location or MPC station code")

    args = arg.parse_args()

    if args.debug:
        ic.enable()
        ic(sys.version_info)
        ic(args)
    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()

    # Default location is Hakos, Namibia for personal reasons ;-)
    loc = EarthLocation(lon=16.36167*u.deg , lat=-23.23639*u.deg, height=1853*u.m)
    tz  = ZoneInfo("Africa/Windhoek")
    if args.location:
        loc = get_location(args.location)
    ic(loc, loc.to_geodetic())

    if args.time:
        time = Time(args.time, location=loc)
    else:
        time = Time(Time.now(), location=loc)
    ic(time)

    observer = Observer(location=loc, 
                        # pressure=1000*u.hPa, 
                        # temperature=20*u.deg_C,
                        # relative_humidity=0.4, 
                        description=loc.info.name,
                        # timezone="Africa/Windhoek"
                        )
    ic(observer)
    ic(observer.moon_altaz(time))

    # To compare with timeanddate.com / LunaSolCal &c.:
    # pressure=0
    # horizon=-0.8333*u.deg 
    # is strictly required!
    # See https://astroplan.readthedocs.io/en/stable/faq/precision.html#how-are-sunrise-and-sunset-defined 

    format = "%Y-%m-%d %H:%M:%S %Z"

    verbose(f"time:              {time.to_datetime(timezone=timezone.utc).strftime(format)}")

    dt_moon_rise = observer.moon_rise_time(time, horizon=-0.8333*u.deg).to_datetime(timezone=timezone.utc)
    dt_moon_set = observer.moon_set_time(time, horizon=-0.8333*u.deg).to_datetime(timezone=timezone.utc)
    dt_moon_rise_local = dt_moon_rise.astimezone(tz)
    dt_moon_set_local = dt_moon_set.astimezone(tz)
    ic(dt_moon_rise, dt_moon_set, dt_moon_rise_local, dt_moon_set_local)
    verbose(f"nearest moon rise: {dt_moon_rise.strftime(format)}  {dt_moon_rise_local.strftime(format)}")
    verbose(f"nearest moon set:  {dt_moon_set.strftime(format)}  {dt_moon_set_local.strftime(format)}")

    dt_sun_rise = observer.sun_rise_time(time, horizon=-0.8333*u.deg).to_datetime(timezone=timezone.utc)
    dt_sun_set = observer.sun_set_time(time, horizon=-0.8333*u.deg).to_datetime(timezone=timezone.utc)
    dt_sun_rise_local = dt_sun_rise.astimezone(tz)
    dt_sun_set_local = dt_sun_set.astimezone(tz)
    ic(dt_sun_rise, dt_sun_set, dt_sun_rise_local, dt_sun_set_local)
    verbose(f"nearest sun rise:  {dt_sun_rise.strftime(format)}  {dt_sun_rise_local.strftime(format)}")
    verbose(f"nearest sun set:   {dt_sun_set.strftime(format)}  {dt_sun_set_local.strftime(format)}")



if __name__ == "__main__":
    main()
