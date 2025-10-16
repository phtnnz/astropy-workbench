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
# Version 0.1 / 2025-10-16
#       Started some solar eclipse related calculations, compute separation
#       and eclipse phase type for location/time, search for contact times
#       in +/- 2h interval around max eclipse

VERSION = "0.1 / 2025-10-16"
AUTHOR  = "Martin Junius"
NAME    = "eclipse"

import sys
import argparse
from typing import Tuple, Any

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
import astropy.constants as const
import numpy as np

from astropy.coordinates import solar_system_ephemeris
from astropy.coordinates import get_body

# Local modules
from verbose import verbose, warning, error
from astroutils import altaz_to_string, location_to_string, get_location

# Moon equatorial radius
R_moon = 1738.1 * u.km
# Sun equatorial radius, IAU
R_sun = const.R_sun
S_sun = np.asin(R_sun / (1 * u.au))
# See SEML and https://iopscience.iop.org/article/10.3847/1538-4365/ac1279
S_sun2 = 959.95 * u.arcsec
R_sun2 = np.sin(S_sun2) * 1 * u.au



# Command line options
class Options:
    pass



def test_tse2026() -> Tuple[EarthLocation, Time]:
    """
    Get location and time for test case: solar eclipse 12 Aug 2026, Burgos, Spain

    Returns
    -------
    Tuple[EarthLocation, Time]
        Location and time of max eclipse
    """

    # Test case TSE 2026 @ Burgos, Spain
    #
    # 42° 21' 01.68" N	  <—>  	42.35047°
    # 3° 41' 21.68" W	  <—>  	-3.68935°
    # 891.0m
    #
    # Maximum eclipse (MAX) : 2026/08/12 18:29:17.6
    loc = EarthLocation(lon=-3.68935*u.degree, lat=42.35047*u.degree, height=891*u.m)
    time = Time("2026-08-12 18:29:17.6", location=loc)
    ic(loc, time)
    return loc, time



def sun_and_moon_series(loc: EarthLocation, time: Time) -> Tuple[Time, Time, Time, Time, Time]:
    ic(loc)

    verbose("searching for contact times ... (takes some time)")
    # Tests
    # # +/- 2 h, 1 min intervals
    # time_interval = time + np.linspace(-2, 2, 4*60+1)*u.hour
    # # +/- 2 min, 1 sec intervals
    # time_interval = time + np.linspace(-2, 2, 4*60+1)*u.min
    # # +/- 1 min, 0.1 sec intervals
    # time_interval = time + np.linspace(-1, 1, 2*60*10+1)*u.min

    # +/- 2 h, 0.5 sec intervals
    time_interval = time + np.linspace(-2, 2, 4*60*60*2+1)*u.hour
    ic(time, time_interval)

    sun  = get_body("sun", time_interval, loc)
    moon = get_body("moon", time_interval, loc)
    ic(sun, moon)
    sun_dist  = sun.distance.to(u.km)
    moon_dist = moon.distance.to(u.km)
    ic(sun_dist, moon_dist)
    sun_altaz  = sun.transform_to( AltAz(obstime=time_interval, location=loc) )
    moon_altaz = moon.transform_to( AltAz(obstime=time_interval, location=loc) )
    ic(sun_altaz, moon_altaz)
    sep = sun.separation(moon).to(u.arcmin)
    ic(sep)
    sun_size = (2 * np.asin(R_sun2 / sun_dist)).to(u.arcmin)
    moon_size = (2 * np.asin(R_moon / moon_dist)).to(u.arcmin)
    ic(sun_size, moon_size)
    partial_sep = (sun_size + moon_size).to(u.arcmin) / 2
    total_sep   = (sun_size - moon_size).to(u.arcmin) / 2
    ic(partial_sep, total_sep)

    last_type = None
    last_time = None
    last_sep  = 180 * u.degree
    max_found = False
    C1        = None
    C2        = None
    MAX       = None
    C3        = None
    C4        = None

    # Search for changes in phase type
    ##FIXME: A=annular doesn't work!
    for t, sep1, p_sep, t_sep in zip(time_interval, sep, partial_sep, total_sep):
        type = type_from_sep(sep1, p_sep, t_sep)
        if not last_type:
            last_type = type
            last_time = t

        if sep1 < last_sep:
            last_sep = sep1
        elif not max_found:
            verbose(last_time, type, "MAX")
            MAX = last_time
            max_found = True

        if type=="P" and last_type=="-":
            verbose(t, type, "C1")
            C1 = t
        if type=="T" and last_type=="P":
            verbose(t, type, "C2")
            C2 = t
        if type=="P" and last_type=="T":
            verbose(last_time, last_type, "C3")
            C3 = last_time
        if type=="-" and last_type=="P":
            verbose(last_time, last_type, "C4")
            C4 = last_time

        last_type = type
        last_time = t 

    return C1, C2, MAX, C3, C4



def sun_and_moon(loc: EarthLocation, time: Time) -> Tuple[Angle, Angle, Angle, Angle]:
    ic(loc)

    sun  = get_body("sun", time, loc)
    moon = get_body("moon", time, loc)
    ic(sun, moon)

    sun_dist  = sun.distance.to(u.km)
    moon_dist = moon.distance.to(u.km)
    verbose(f"sun distance {sun_dist:_.3f} = {sun_dist.to(u.au):_.3f}")
    verbose(f"moon distance {moon_dist:_.3f} = {moon_dist.to(u.au):_.3f}")

    ic(R_moon, R_sun.to(u.km), R_sun2.to(u.km), (1*u.au).to(u.km))

    sun_size = 2 * np.asin(R_sun2 / sun_dist)
    verbose(f"sun angular size (Quaglia) {sun_size.to(u.degree):.3f} = {sun_size.to(u.arcmin):.3f} = {sun_size.to(u.arcsec):.3f}")
    moon_size = 2 * np.asin(R_moon / moon_dist)
    verbose(f"moon angular size {moon_size.to(u.degree):.3f} = {moon_size.to(u.arcmin):.3f} = {moon_size.to(u.arcsec):.3f}")

    sun_altaz  = sun.transform_to( AltAz(obstime=time, location=loc) )
    moon_altaz = moon.transform_to( AltAz(obstime=time, location=loc) )
    ic(sun_altaz, moon_altaz)
    verbose(f"sun position {altaz_to_string(sun_altaz.alt, sun_altaz.az)}")
    verbose(f"moon position {altaz_to_string(moon_altaz.alt, moon_altaz.az)}")

    # Separation between sun and moon
    sep = sun.separation(moon).to(u.arcmin)
    ic(sep)

    # / 2 for radii
    partial_sep = (sun_size + moon_size).to(u.arcmin) / 2
    total_sep   = (sun_size - moon_size).to(u.arcmin) / 2
    ratio       = moon_size / sun_size
    ic(partial_sep, total_sep, ratio)
    type = type_from_sep(sep, partial_sep, total_sep)
    ic(type)
    verbose(f"moon/sun ratio {ratio:.4f}")
    verbose(f"moon-sun separation {sep:.2f}, phase {type}")

    # return separation for minimum search and root finding
    #  min=MAX  root=C1/C4       root=C2/C3 T   root=C2/C3 A
    return sep, sep-partial_sep, sep+total_sep, sep-total_sep



def type_from_sep(sep: Angle, partial_sep: Angle, total_sep: Angle) -> str:
    type = "-"
    if sep <= partial_sep:
        type = "P"
    if total_sep <= 0 and sep <= -total_sep:
        type = "T"
    if total_sep > 0 and sep <= total_sep:
        type = "A"
    # ic(sep, partial_sep, total_sep, type)
    return type



def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = "Solar eclipse calculations",
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("-t", "--time", help="time (UTC), default now")
    arg.add_argument("-l", "--location", help="coordinates, named location or MPC station code")
    arg.add_argument("-e", "--ephemeris", help="use EPHEMERIS, default \"de440\"")
    arg.add_argument("--tse2026", action="store_true", help="test case TSE 12 Aug 2026")

    args = arg.parse_args()

    if args.debug:
        ic.enable()
        ic(sys.version_info)
        ic(args)
    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()

    # Location and time
    loc = None
    time = None
    if args.tse2026:
        loc, time = test_tse2026()
    if args.location:
        loc = get_location(args.location)
    if loc == None:
        error("no location specified")
    ic(loc, loc.to_geodetic())
    verbose(f"location {location_to_string(loc)}")
    if args.time:
        time = Time(args.time, location=loc)
    elif not time:
        time = Time(Time.now(), location=loc)
    ic(time)
    verbose(f"time UTC {time}")

    # Ephemeris
    ephemeris = args.ephemeris or "de440"
    verbose(f"using {ephemeris} ephemeris")
    solar_system_ephemeris.set(ephemeris)

    sun_and_moon(loc, time)
    sun_and_moon_series(loc, time)



if __name__ == "__main__":
    main()
