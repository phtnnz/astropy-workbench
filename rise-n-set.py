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
#       Sun/Moon rise/set with astroplan
# Version 0.2 / 2025-06-21
#       Added civil/nautical/astronomical dusk/dawn

VERSION = "0.2 / 2025-06-21"
AUTHOR  = "Martin Junius"
NAME    = "rise-n-set"

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
from astropy.coordinates import EarthLocation
import astropy.units as u
from astropy.time        import Time

# Astroplan
from astroplan import Observer

# Local modules
from verbose import verbose, warning, error
from astroutils import get_location



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

    dt_c_dusk = observer.twilight_evening_civil(time).to_datetime(timezone=timezone.utc)
    dt_n_dusk = observer.twilight_evening_nautical(time).to_datetime(timezone=timezone.utc)
    dt_a_dusk = observer.twilight_evening_astronomical(time).to_datetime(timezone=timezone.utc)
    dt_c_dusk_local = dt_c_dusk.astimezone(tz)
    dt_n_dusk_local = dt_n_dusk.astimezone(tz)
    dt_a_dusk_local = dt_a_dusk.astimezone(tz)
    ic(dt_a_dusk, dt_n_dusk, dt_c_dusk)
    ic(dt_a_dusk_local, dt_n_dusk_local, dt_c_dusk_local)
    verbose(f"nearest civil dusk {dt_c_dusk.strftime(format)}  {dt_c_dusk_local.strftime(format)}")
    verbose(f"nearest naut. dusk {dt_n_dusk.strftime(format)}  {dt_n_dusk_local.strftime(format)}")
    verbose(f"nearest astr. dusk {dt_a_dusk.strftime(format)}  {dt_a_dusk_local.strftime(format)}")

    dt_c_dawn = observer.twilight_morning_civil(time).to_datetime(timezone=timezone.utc)
    dt_n_dawn = observer.twilight_morning_nautical(time).to_datetime(timezone=timezone.utc)
    dt_a_dawn = observer.twilight_morning_astronomical(time).to_datetime(timezone=timezone.utc)
    dt_c_dawn_local = dt_c_dawn.astimezone(tz)
    dt_n_dawn_local = dt_n_dawn.astimezone(tz)
    dt_a_dawn_local = dt_a_dawn.astimezone(tz)
    ic(dt_a_dawn, dt_n_dawn, dt_c_dawn)
    ic(dt_a_dawn_local, dt_n_dawn_local, dt_c_dawn_local)
    verbose(f"nearest civil dawn {dt_c_dawn.strftime(format)}  {dt_c_dawn_local.strftime(format)}")
    verbose(f"nearest naut. dawn {dt_n_dawn.strftime(format)}  {dt_n_dawn_local.strftime(format)}")
    verbose(f"nearest astr. dawn {dt_a_dawn.strftime(format)}  {dt_a_dawn_local.strftime(format)}")


if __name__ == "__main__":
    main()
