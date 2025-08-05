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

VERSION = "0.1 / 2025-01-30"
AUTHOR  = "Martin Junius"
NAME    = "plot-sky"

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
import matplotlib.pyplot as plt
from astroplan import FixedTarget, Observer
from astroplan.plots import plot_airmass, plot_altitude, plot_sky, plot_finder_image, plot_parallactic

# Local modules
from verbose import verbose, warning, error
from astroutils import get_location, get_coord, coord_to_string



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
    arg.add_argument("object", nargs="+", help="object name or \"RA DEC\" coordinates")

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

    # .strftime() format
    format = "%Y-%m-%d %H:%M:%S %Z"

    verbose(f"time: {time.to_datetime(timezone=timezone.utc).strftime(format)}")

    # test object
    # obj = "03:57:25.611 -46:11:07.57" # SN 2024abfo precise position
    # name = "SN 2024abfo"
    # coord = SkyCoord(obj, unit=(u.hour, u.deg))
    # ic(obj, coord)

    for obj in args.object:
        coord = get_coord(obj, args.query_simbad)
        target = FixedTarget(name=obj, coord=coord)
        verbose(f"object: {coord_to_string(coord)}")
        ic(target)

        midnight = observer.midnight(time, which="next")
        ic(midnight)
        verbose(f"next midnight: {midnight.to_datetime(timezone=timezone.utc).strftime(format)}")
        time = midnight

        time_interval50 = time + np.linspace(-4, 5, 50)*u.hour
        time_interval10 = time + np.linspace(-4, 5, 10)*u.hour
        moon_vals50 = observer.moon_altaz(time_interval50)
        moon_vals10 = observer.moon_altaz(time_interval10)
        ic(moon_vals10)

        ic("plot_altitude")
        plot_altitude(target, observer, time_interval50, brightness_shading=True)
        plot_altitude(moon_vals50, observer, time_interval50)
        plt.legend(loc='lower left').get_texts()[1].set_text("Moon")
        plt.tight_layout()
        verbose("altitude plot")
        plt.savefig("tmp/plot-altitude.png", bbox_inches="tight")
        plt.close()

        # ic("plot_airmass")
        # plot_airmass(target, observer, time, brightness_shading=True, altitude_yaxis=True)
        # plt.savefig("tmp/plot-airmass.png", bbox_inches="tight")
        # plt.close()

        ic("plot_sky")
        plot_sky(target, observer, time_interval10)
        plot_sky(moon_vals10, observer, time_interval10)
        plt.legend(loc='lower left').get_texts()[1].set_text("Moon")
        verbose("sky plot")
        plt.savefig("tmp/plot-sky.png", bbox_inches="tight")
        plt.close()

        # Doesn't work anymore, see https://github.com/astropy/astroplan/pull/591
        # ax, hdu = plot_finder_image(target)
        # plt.savefig("tmp/plot-sky4.png", bbox_inches="tight")
        # plt.close()

        # not working inside VSCode terminal
        # plt.show()



if __name__ == "__main__":
    main()
