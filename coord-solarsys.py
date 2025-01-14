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
# Version 0.1 / 2025-01-14
#       Get coordinates for solar system bodies

import sys
import argparse
import re

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

from astropy.coordinates import solar_system_ephemeris
from astropy.coordinates import get_body_barycentric, get_body


# Local modules
from verbose import verbose, warning, error
from mpclocation import mpc_station_location

VERSION = "0.1 / 2025-01-14"
AUTHOR  = "Martin Junius"
NAME    = "coord-solarsys"



# Command line options
class Options:
    pass



def ra_from_lst_ha(lst: Angle, ha: Angle):
    ra = lst - ha
    ra.wrap_at(24*u.hourangle, inplace=True)
    ic(lst, ha, ra)
    return ra


def ra_dec_to_string(ra: Angle, dec: Angle):
    return f"RA={ra.to_string(unit=u.hour, precision=2)} DEC={dec.to_string(unit=u.degree, precision=2)}"

def angle_to_string(a: Angle):
    return f"{a.to_string(unit=u.degree, precision=2)}"

def hourangle_to_string(a: Angle):
    return f"{a.to_string(unit=u.hour, precision=2)}"


def coord_to_altaz(obj: str, loc: EarthLocation, time: Time):
    verbose(f"object {obj}, time {time}")

    # GCRS coord 
    coord = get_body(obj, time, loc)
    ic(coord)
    verbose(f"GCRS coord                       {ra_dec_to_string(coord.ra, coord.dec)}")

    coord_pgc = coord.transform_to(PrecessedGeocentric(equinox=Time(time, format="jyear"), obstime=time))
    ic(coord_pgc)
    verbose(f"PrecessedGeocentric JNOW coord   {ra_dec_to_string(coord_pgc.ra, coord_pgc.dec)}")

    # Transform to topocentric HADec
    hadec_jnow = coord.transform_to(HADec(obstime=time, location=loc))
    hadec_jnow_w_refraction = coord.transform_to(HADec(obstime=time, location=loc, 
                                                 pressure=1000*u.hPa, temperature=20*u.deg_C,
                                                 relative_humidity=0.4, obswl=0.54*u.micron))
    ic(hadec_jnow, hadec_jnow_w_refraction)
    
    lst_mean     = time.sidereal_time("mean")
    ra_mean      = ra_from_lst_ha(lst_mean, hadec_jnow.ha)
    ic(lst_mean, ra_mean)
    verbose(f"LST={hourangle_to_string(lst_mean)}, HA={hourangle_to_string(hadec_jnow.ha)}")
    verbose(f"Topocentric JNOW coord           {ra_dec_to_string(ra_mean, hadec_jnow.dec)}")

    # Transform to topocentric AltAz
    altaz = hadec_jnow.transform_to( AltAz(obstime=time, location=loc) )
    ic(altaz)
    verbose(f"Alt={angle_to_string(altaz.alt)} Az={angle_to_string(altaz.az)}")

    # Calculate parallactic angle https://en.wikipedia.org/wiki/Parallactic_angle 
    # Based on https://github.com/lsst-ts/ts_observatory_control/blob/develop/python/lsst/ts/observatory/control/utils/utils.py
    # Eqn (14.1) of Meeus' Astronomical Algorithms
    # q = 0     N axis aligned with Meridian
    # q < 0     N axis turned q degrees left
    # q > 0     N axis turned q degrees right
    H = hadec_jnow.ha.radian
    q = Angle( np.arctan2( np.sin(H),
                           (np.tan(loc.lat.radian) * np.cos(coord.dec.radian) 
                            - np.sin(hadec_jnow.dec.radian) * np.cos(H))             ), u.rad).to(u.deg)
    ic(q)
    verbose(f"parallactic angle={angle_to_string(q)}")



def get_location(name: str) -> EarthLocation:
    """
    Try to interpret location name as address, site name, MPC code

    :param name: location name
    :type name: str
    :return: location object
    :rtype: EarthLocation
    """
    loc = None
    m = re.match(r'^([0-9.]+) ([+-]?[0-9.]+) ([0-9.]+)$', name)
    if m:
        (lon, lat, height) = [ float(v) for v in m.groups() ]
        loc = EarthLocation(lon=lon*u.degree, lat=lat*u.degree, height=height*u.m)
        verbose(f"location {lon=} {lat=} {height=}")

    if loc == None:
        try:
            loc = EarthLocation.of_site(name)
        except errors.UnknownSiteException as e:
            verbose(f"location {name} not in astropy database")
            loc = None

    if loc == None:
        try:
            loc = mpc_station_location(name)
        except (LookupError, ValueError):
            verbose(f"location {name} not an MPC station code")
            loc = None

    if loc == None:
        ic(EarthLocation.get_site_names())
        error(f"named location {name} not found")

    return loc



def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = "Astropy solar system bodies coordinates",
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("-t", "--time", help="time (UTC) for JNOW coordinates, default now")
    arg.add_argument("-l", "--location", help="coordinates, named location or MPC station code")
    arg.add_argument("object", nargs="+", help="solar system body")

    args = arg.parse_args()

    if args.debug:
        ic.enable()
        ic(sys.version_info)
        ic(args)
    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()

    # Default location is Hakos, Namibia for personal reasons ;-)
    loc      = EarthLocation(lon=16.36167*u.deg , lat=-23.23639*u.deg, height=1825*u.m)
    if args.location:
        loc = get_location(args.location)
    ic(loc, loc.to_geodetic())

    if args.time:
        time = Time(args.time, location=loc)
    else:
        time = Time(Time.now(), location=loc)
    ic(time)

    # Ephemeris, TODO: add option for JPL
    solar_system_ephemeris.set('builtin')

    for obj in args.object:
        verbose(f"object {obj}")
        verbose(f"time (UTC) {time}")
        coord_to_altaz(obj, loc, time)            



if __name__ == "__main__":
    main()
