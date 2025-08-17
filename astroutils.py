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
#       Utility functions moved to this module

import re

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
from astropy.coordinates import errors
import astropy.units as u
from astropy.time        import Time, TimeDelta
import numpy as np

# Local modules
from verbose import verbose, warning, error
from mpclocation import mpc_station_location
from querysimbad import query_simbad


VERSION = "0.1 / 2025-01-27"
AUTHOR  = "Martin Junius"
NAME    = "astroutils"



def ra_from_lst_ha(lst: Angle, ha: Angle) -> Angle:
    """
    Compute RA from local sidereal time and hour angle

    :param lst: local sidereal time
    :type lst: Angle
    :param ha: hour angle
    :type ha: Angle
    :return: RA
    :rtype: Angle
    """
    ra = lst - ha
    ra.wrap_at(24*u.hourangle, inplace=True)
    ic(lst, ha, ra)
    return ra



def coord_to_string(coord: SkyCoord, short: bool=False) -> str:
    """
    Format SkyCoord as string

    :param coord: coordinates object
    :type coord: SkyCoord
    :return: _description_
    :rtype: str
    """
    return ra_dec_to_string(coord.ra, coord.dec, short)

def ra_dec_to_string(ra: Angle, dec: Angle, short: bool=False) -> str:
    """
    Format RA and DEC as string

    :param ra: RA
    :type ra: Angle
    :param dec: DEC
    :type dec: Angle
    :return: formatted string "RA=... DEC=..."
    :rtype: str
    """
    if short:
        return f"{ra.to_string(unit=u.hour, precision=0)} {dec.to_string(unit=u.degree, precision=0)}"
    else:
        return f"RA={ra.to_string(unit=u.hour, precision=2)} DEC={dec.to_string(unit=u.degree, precision=2)}"

def angle_to_string(a: Angle, decimal: bool=False) -> str:
    """
    Format angle as string

    :param a: angle
    :type a: Angle
    :param decimal: flag for decimal output, default False
    :type decimal: bool
    :return: formatted string angle in degrees, minutes, seconds
    :rtype: str
    """
    return f"{a.to_string(unit=u.degree, decimal=decimal, precision=2)}"

def hourangle_to_string(a: Angle) -> str:
    """
    Format hour angle as string

    :param a: hour angle
    :type a: Angle
    :return: formatted string hour angle in hours, minutes, seconds
    :rtype: str
    """
    return f"{a.to_string(unit=u.hour, precision=2)}"



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



def location_to_string(loc: EarthLocation) -> str:
    return f"lon={loc.to_geodetic().lon.to_string(unit=u.degree, precision=1)} lat={loc.to_geodetic().lat.to_string(unit=u.degree, precision=1)} height={loc.to_geodetic().height.to_string(precision=0)}"



def get_coord(name: str, simbad=False) -> SkyCoord:
    """
    Get coordinates for object, either by converting "RA DEC" string, or querying Simbad

    :param name: object name or "RA DEC" string
    :type name: str
    :param simbad: query simbad for object name, defaults to False
    :type simbad: bool, optional
    :return: coordinates object
    :rtype: SkyCoord
    """
    # Try ICRS coords first
    try:
        coord = SkyCoord(name, unit=(u.hour, u.deg))
    except ValueError:
        ic("not ra/dec coordinates")
        coord = None
        pass

    # Try to query Simbad for object name
    if not coord and simbad:
        coord = query_simbad(name, w_velocity=False)

    ic(name, coord)
    return coord;



if __name__ == "__main__":
    error("no main() function")
