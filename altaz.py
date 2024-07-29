#!/usr/bin/env python

# Copyright 2024 Martin Junius
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
# Version 0.1 / 2024-07-29
#       Calculate Alt/Az coordinates for objects

import sys
import argparse
from datetime import datetime, timezone

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
import astropy.units as u
from astropy.time        import Time, TimeDelta
import numpy as np

# Local modules
from verbose import verbose, warning, error


VERSION = "0.1 / 2024-07-29"
AUTHOR  = "Martin Junius"
NAME    = "altaz"



# Command line options
class Options:
    pass



def convert_to_alt_az(coord, utc_time):
    """ Convert RA/DEC to local Alt/Az """
    # Based on https://docs.astropy.org/en/stable/generated/examples/coordinates/plot_obs-planning.html

    # IAS Hakos, Namibia
    loc  = EarthLocation(lat=-23.233*u.deg, lon=16.350*u.deg , height=1834*u.m)
    time = Time(utc_time, location=loc)
    lst  = time.sidereal_time("mean")
    hour_angle = lst - coord.ra 
    ic(loc, time, lst, coord.ra, hour_angle)
    verbose(f"  location: lat={loc.lat:.3f}, lon={loc.lon:.3f}, height={loc.height:.0f}")

    altaz = coord.transform_to( AltAz(obstime=time, location=loc) )
    ic(altaz)
    verbose(f"  Az={altaz.az:.3f}, Alt={altaz.alt:.3f}")

    # Calculate parallactic angle https://en.wikipedia.org/wiki/Parallactic_angle 
    # Based on https://github.com/lsst-ts/ts_observatory_control/blob/develop/python/lsst/ts/observatory/control/utils/utils.py
    # Eqn (14.1) of Meeus' Astronomical Algorithms
    # q = 0     N axis aligned with Meridian
    # q < 0     N axis turned q degrees left
    # q > 0     N axis turned q degrees right
    H = hour_angle.radian
    q = Angle( np.arctan2( np.sin(H),
                           (np.tan(loc.lat.radian) * np.cos(coord.dec.radian) 
                            - np.sin(coord.dec.radian) * np.cos(H))             ), u.rad).to(u.deg)
    ic(q)
    verbose(f"  parallactic angle={q:.3f}")



def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = "Generic python script template",
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("-T", "--time", help="ISO format time incl. timezone, e.g. 2024-07-09T22:00:00+02:00, default now")
    arg.add_argument("object", nargs="+", help="object name(s)")

    args = arg.parse_args()

    if args.debug:
        ic.enable()
        ic(sys.version_info)
        ic(args)
    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()
    # ... more options ...
    if args.time:
        loc_time = datetime.fromisoformat(args.time)
    else:
        loc_time = datetime.now().astimezone()
    utc_time = loc_time.astimezone(timezone.utc)
    ic(loc_time, utc_time)

    # ... the action starts here ...
    for obj in args.object:
        coord = SkyCoord.from_name(obj)
        verbose(f"{obj}: RA/DEC = {coord.to_string("hmsdms")}")
        convert_to_alt_az(coord, utc_time)



if __name__ == "__main__":
    main()