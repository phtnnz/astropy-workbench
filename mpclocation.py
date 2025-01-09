#!/usr/bin/env python

# Copyright 2025 Martin Junius
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

import sys
import argparse
import numpy as np

# AstroPy
from astropy.coordinates import EarthLocation  # High-level coordinates
import astropy.constants as C
import astropy.units as u
from astropy.coordinates import Angle, Latitude, Longitude  # Angles

# Astroquery
from astroquery.mpc import MPC

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()
# Local modules
from verbose import verbose, warning, error


VERSION = "0.1 / 2024-01-09"
AUTHOR  = "Martin Junius"
NAME    = "mpclocation"



def test_m49():
    """
    Run test case M49
    """
    loc1 = EarthLocation(lat=-23.23639*u.deg, lon=16.36167*u.deg , height=1825*u.m) # M49, Hakos, Namibia
    ic(loc1, loc1.to_geodetic())

    loc2 = mpc_station_location("M49")
    ic(loc2, loc2.to_geodetic())

    # https://projectpluto.com/parallax.htm
    #
    # M49 !+016.361720000  -23.236548535  1853.495   IAS Remote Observatory, Hakos
    # M49  16.361720.919630-0.392206IAS Remote Observatory, Hakos
    # Longitude   16.361720000 = +16 21 42.19199
    # Latitude  -23.236548535 = -23 14 11.57472
    # Altitude 1853.49468 meters above the WGS84 ellipsoid (_not_ above sea level/geoid)
    # Parallax constants 0.91963000000 -0.39220600000
    # In meters: 5865526.12931 -2501543.60022
    # xyz in Earth radii -0.7300066 -0.5592940 -0.3922060
    # xyz in meters      -4656081.99331 -3567253.45962 -2501543.60022
    # This point is somewhere in Namibia
    lon_p = 16.361720000
    lat_p = -23.236548535
    alt_p = 1853.49468
    ic("projectpluto results in comparison")
    ic(lon_p, lat_p, alt_p)



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

    return mpc_parallax_to_location(long, rho_cos_phi, rho_sin_phi)



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
    # R_earth = C.R_earth.to_value(u.m) # not precise enough!
    R_earth = 6378137.0                 # GRS 80/WG S84 value (Wikipedia)
                                        # https://en.wikipedia.org/wiki/World_Geodetic_System
    lon = longitude.to_value(u.radian)
    ic(R_earth, lon)

    # compute cartesian components of geocentric position
    x = R_earth * rho_cos_phi * np.cos(lon) * u.m
    y = R_earth * rho_cos_phi * np.sin(lon) * u.m
    z = R_earth * rho_sin_phi               * u.m
    ic(x, y, z)
    return EarthLocation.from_geocentric(x, y, z)



### Test run as a command line script ###
def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = "Convert MPC station location to EarthLocation",
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("--test_m49", action="store_true", help="MPC station M49 test case")
    arg.add_argument("station", nargs="*", help="MPC station code")

    args = arg.parse_args()

    if args.debug:
        ic.enable()
        ic(args, sys.version_info, sys.path)
    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()

    if args.test_m49:
        ic.enable()
        test_m49()
    else:
        for station in args.station:
            pass



if __name__ == "__main__":
    main()
