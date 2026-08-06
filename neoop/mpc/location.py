#!/usr/bin/env python

# Copyright 2024-2026 Martin Junius
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
# Version 0.1 / 2024-01-09
#       EarthLocation from MPC longitude/parallax data
# Version 0.2 / 2025-01-24
#       Formatted verbose output for main, add station name to
#       EarthLocation.info.name
# Version 1.0 / 2026-06-16
#       Moved and adapted to new directory structure under neoop/

import re
import numpy as np

# AstroPy
from astropy.coordinates import EarthLocation  # High-level coordinates
import astropy.constants as C
import astropy.units as u
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
from astropy.coordinates import errors, name_resolve

# Astroquery
from astroquery.mpc import MPC

from icecream import ic
# Disable debugging
ic.disable()

# Local modules
from utils.verbose import verbose, warning, error


VERSION = "1.0 / 2026-06-16"
AUTHOR  = "Martin Junius"
NAME    = "mpc.location"



# R_earth = C.R_earth.to_value(u.m) # astropy.contants is not precise enough!
                                    # 6378100.0 m
R_earth = 6378137.0                 # GRS 80/WGS 84 value (Wikipedia)
                                    # https://en.wikipedia.org/wiki/World_Geodetic_System
# R_earth = 6378136.6               # UNITED STATES NAVAL OBSERVATORY, CIRCULAR NO. 179, 2005
                                    # The IAU Resolutions on Astronomical Reference Systems,
                                    # Time Scales, and Earth Rotation Models
                                    # https://aa.usno.navy.mil/downloads/Circular_179.pdf



def mpc_station_location(station: str) -> EarthLocation:
    """
    Get geocentric location for MPC station code

    :param station: MPC station code
    :type station: str
    :return: location object
    :rtype: EarthLocation
    """
    (long, rho_cos_phi, rho_sin_phi, name) = MPC.get_observatory_location(station)
    ic(long, rho_cos_phi, rho_sin_phi, name)

    loc = mpc_parallax_to_location(long, rho_cos_phi, rho_sin_phi)
    loc.info.name = name
    return loc



def mpc_parallax_to_location(longitude: Longitude, rho_cos_phi: float, rho_sin_phi: float) -> EarthLocation:
    """
    Convert MPC station location data to geocentric location

    :param longitude: longitude from station data
    :type longitude: Longitude
    :param rho_cos_phi: rho*cos(phi) parallax in Earth radii
    :type rho_cos_phi: float
    :param rho_sin_phi: rho*sin(phi) parallax in Earth radii
    :type rho_sin_phi: float
    :return: location object
    :rtype: EarthLocation
    """
    lon = longitude.to_value(u.radian)
    ic(R_earth, lon)

    # compute cartesian components of geocentric position
    x = R_earth * rho_cos_phi * np.cos(lon) * u.m
    y = R_earth * rho_cos_phi * np.sin(lon) * u.m
    z = R_earth * rho_sin_phi               * u.m
    ic(x, y, z)
    return EarthLocation.from_geocentric(x, y, z)



def get_location(name: str) -> EarthLocation:
    """
    Try to interpret location name as lon/lat coordinates, MPC station code,
    site name, or openstreetmap address

    Parameters
    ----------
    name : str
        Location text

    Returns
    -------
    EarthLocation
        Astropy location object
    """
    loc = None
    m = re.match(r'^(-?[0-9.]+) ([+-]?[0-9.]+) ([0-9.]+)$', name)
    if m:
        (lon, lat, height) = [ float(v) for v in m.groups() ]
        loc = EarthLocation(lon=lon*u.degree, lat=lat*u.degree, height=height*u.m)
        verbose(f"location {lon=} {lat=} {height=}")

    if loc == None:
        try:
            loc = mpc_station_location(name)
        except (LookupError, ValueError):
            verbose(f"location {name} not an MPC station code")
            loc = None

    if loc == None:
        try:
            loc = EarthLocation.of_site(name)
        except errors.UnknownSiteException as e:
            verbose(f"location {name} not in astropy database")
            loc = None

    if loc == None:
        try:
            loc = EarthLocation.of_address(name)
        except name_resolve.NameResolveError as e:
            verbose(f"location {name} not a recognized address")
            loc = None

    if loc == None:
        error(f"named location {name} not found")

    return loc
