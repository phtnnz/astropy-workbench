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
# Version 0.1 / 2025-01-30
#       First steps with astroplan
# Version 0.2 / 2025-08-06
#       Somewhat complete version plotting altitude and sky for
#       multiple objects
# Version 0.3 / 2025-08-12
#       Added support for CSV list with objects, improved plots

VERSION = "0.3 / 2025-08-12"
AUTHOR  = "Martin Junius"
NAME    = "plot-sky"

import sys
import argparse
import csv
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
import numpy as np

# Astroplan
import matplotlib.pyplot as plt
from astroplan import FixedTarget, Observer
from astroplan.plots import plot_airmass, plot_altitude, plot_sky, plot_finder_image, plot_parallactic

# Local modules
from verbose import verbose, warning, error
from astroutils import get_location, get_coord, coord_to_string, location_to_string



# Command line options
class Options:
    pass



def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = "Plot objects and sky",
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("-t", "--time", help="time (UTC) for computation, default now")
    arg.add_argument("-l", "--location", help="coordinates, named location or MPC station code")
    arg.add_argument("-q", "--query-simbad", action="store_true", help="query Simbad for OBJECT name")
    arg.add_argument("-A", "--altitude", action="store_true", help="plot altitude (default: sky)")
    arg.add_argument("-f", "--file", help="read list of objects from file")
    arg.add_argument("object", nargs="*", help="object name (-q required) or [name=]\"RA DEC\" coordinates")

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
    verbose(f"location: {location_to_string(loc)}")

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

    # .strftime() format
    format = "%Y-%m-%d %H:%M:%S %Z"

    verbose(f"time: {time.to_datetime(timezone=timezone.utc).strftime(format)}")
    midnight = observer.midnight(time, which="next")
    ic(midnight)
    verbose(f"next midnight: {midnight.to_datetime(timezone=timezone.utc).strftime(format)}")
    time = midnight

    # Intervals around midnight
    time_interval_full = time + np.linspace(-8, 8, 160)*u.hour
    moon_vals_full = observer.moon_altaz(time_interval_full)
    time_interval_pm5 = time + np.linspace(-5, 5, 11)*u.hour
    moon_vals_pm5 = observer.moon_altaz(time_interval_pm5)
    # Quick hack to set label in legend
    moon_vals_full.name = "Moon"
    moon_vals_pm5.name = "Moon"
    ic(time_interval_pm5, moon_vals_pm5)

    # plot objects from file
    if args.file:
        with open(args.file, newline="") as file:
            line = file.readline()
            sep = ";" if ";" in line else ","
            file.seek(0)
            reader = csv.DictReader(file, delimiter=sep)
            for row in reader:
                name = row.get("Name") or row.get("name") or row.get("target")
                ra   = row.get("RA") or row.get("ra") or row.get("RAJ2000")
                dec  = row.get("DE") or row.get("DEC") or row.get("de") or row.get("dec") or row.get("DEJ2000")
                obj  = f"{ra} {dec}"
                ic(name, ra, dec, obj)

                coord = get_coord(obj, args.query_simbad)
                if not coord:
                    warning(f"{obj} not a coordinate or not found (--query-simbad)")
                    continue

                verbose(f"object: {name}={coord_to_string(coord, short=True)}")
                target = FixedTarget(name=name, coord=coord)
                ic(target)

                if args.altitude:
                    # Must use fmt="", not marker="none" to avoid warnings from plot_date()!
                    plot_altitude(target, observer, time_interval_full, style_kwargs=dict(fmt=""))
                else:
                    plot_sky(target, observer, time_interval_pm5)

    # plot objects from command line
    for obj in args.object:
        if "=" in obj:
            (name, obj) = obj.split("=")
            verbose(f"object name: {name}")
        else:
            name = obj
        ic(name, obj)
        coord = get_coord(obj, args.query_simbad)
        if not coord:
            warning(f"{obj} not a coordinate or not found (--query-simbad)")
            continue

        target = FixedTarget(name=name, coord=coord)
        verbose(f"object coord: {coord_to_string(coord)}")
        ic(target)

        if args.altitude:
            # Must use fmt="", not marker="none" to avoid warnings from plot_date()!
            plot_altitude(target, observer, time_interval_full, style_kwargs=dict(fmt=""))
        else:
            plot_sky(target, observer, time_interval_pm5)

    # Add moon plot and legend
    if args.altitude:
        plot_altitude(moon_vals_full, observer, time_interval_full, brightness_shading=True, style_kwargs=dict(fmt="y--"))
        # Set legend for last curve
        plt.legend(bbox_to_anchor=(1.0, 1.02))
        # plt.legend(loc='upper right').get_texts()[-1].set_text("Moon")
        # plt.tight_layout()
    else:
        plot_sky(moon_vals_pm5, observer, time_interval_pm5, style_kwargs=dict(color="y", marker="x"))
        plt.legend(bbox_to_anchor=(1.48, 1.11))

    plt.savefig("tmp/plot.png", bbox_inches="tight")
    plt.close()


    # ic("plot_airmass")
    # plot_airmass(target, observer, time, brightness_shading=True, altitude_yaxis=True)
    # plt.savefig("tmp/plot-airmass.png", bbox_inches="tight")
    # plt.close()

    # Doesn't work anymore, see https://github.com/astropy/astroplan/pull/591
    # and https://github.com/astropy/astroplan/pull/622 seems to be fixed but not yet
    # in the PyPI package?
    # ax, hdu = plot_finder_image(target)
    # plt.savefig("tmp/plot-sky4.png", bbox_inches="tight")
    # plt.close()

    # not working inside VSCode terminal
    # plt.show()



if __name__ == "__main__":
    main()
