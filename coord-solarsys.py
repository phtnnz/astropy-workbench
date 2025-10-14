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
# Version 0.2 / 2025-01-27
#       Use astroutils module
# Version 0.3 / 2025-10-14
#       Output distance and sun/moon angular size, taking recent discussions
#       about the solar radius into account (IAS vs. Quaglia et al.)

VERSION = "0.3 / 2025-10-14"
AUTHOR  = "Martin Junius"
NAME    = "coord-solarsys"

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
import astropy.constants as const
import numpy as np

from astropy.coordinates import solar_system_ephemeris
from astropy.coordinates import get_body_barycentric, get_body

# Local modules
from verbose import verbose, warning, error
from astroutils import ra_from_lst_ha, ra_dec_to_string, angle_to_string, hourangle_to_string, get_location, location_to_string

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



def coord_to_altaz(obj: str, loc: EarthLocation, time: Time):
    # GCRS coord 
    loc0 = EarthLocation.from_geocentric(0, 0, 0, unit=u.m)
    ic(loc0, loc)
    try:
        coord = get_body(obj, time, loc)
        coord0 = get_body(obj, time, loc0)      # center of Earth
        ic(coord, coord0)
    except KeyError:
        error(f"get_object failed, supported bodies: {",".join(solar_system_ephemeris.bodies)}")
    distance = coord.distance.to(u.km)
    verbose(f"distance {distance:_.3f} = {distance.to(u.au):_.3f}")
    ic(R_moon, R_sun.to(u.km), R_sun2.to(u.km), (1*u.au).to(u.km))
    if obj.lower() == "moon":
        verbose(f"moon radius {R_moon.to(u.km)}")
        size = 2 * np.asin(R_moon / distance)
        verbose(f"moon angular size {size.to(u.degree):.3f} = {size.to(u.arcmin):.3f} = {size.to(u.arcsec):.3f}")
    if obj.lower() == "sun":
        verbose(f"sun radius IAU {R_sun.to(u.km)} = {S_sun.to(u.arcsec):.2f} @ 1 AU")
        verbose(f"sun radius Quaglia et al. 2021 {R_sun2.to(u.km):.0f} = {S_sun2.to(u.arcsec):.2f} @ 1 AU")
        size = 2 * np.asin(R_sun / distance)
        verbose(f"sun angular size (IAU) {size.to(u.degree):.3f} = {size.to(u.arcmin):.3f} = {size.to(u.arcsec):.3f}")
        size = 2 * np.asin(R_sun2 / distance)
        verbose(f"sun angular size (Quaglia) {size.to(u.degree):.3f} = {size.to(u.arcmin):.3f} = {size.to(u.arcsec):.3f}")
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



def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = "Astropy solar system bodies coordinates",
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("-t", "--time", help="time (UTC) for JNOW coordinates, default now")
    arg.add_argument("-l", "--location", help="coordinates, named location or MPC station code")
    arg.add_argument("-e", "--ephemeris", help="use EPHEMERIS (e.g. JPL \"de430\"), default \"builtin\"")
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
    loc = EarthLocation(lon=16.36167*u.deg , lat=-23.23639*u.deg, height=1853*u.m)
    if args.location:
        loc = get_location(args.location)
    ic(loc, loc.to_geodetic())
    verbose(f"location: {location_to_string(loc)}")


    if args.time:
        time = Time(args.time, location=loc)
    else:
        time = Time(Time.now(), location=loc)
    ic(time)

    # Ephemeris
    ephemeris = args.ephemeris or "builtin"
    verbose(f"using {ephemeris} ephemeris")
    solar_system_ephemeris.set(ephemeris)

    for obj in args.object:
        verbose(f"object {obj}")
        verbose(f"time (UTC) {time}")
        coord_to_altaz(obj, loc, time)            



if __name__ == "__main__":
    main()
