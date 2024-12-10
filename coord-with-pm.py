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
# Version 0.1 / 2024-12-10
#       SkyCoord for object at particular date/time with proper motion

import sys
import argparse

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
from astropy.time        import Time, TimeDelta
import astropy.units as u

# Local modules
from verbose import verbose, warning, error
from querysimbad import query_simbad


VERSION = "0.1 / 2024-12-10"
AUTHOR  = "Martin Junius"
NAME    = "coord-with-pm"



# Command line options
class Options:
    pass



# Proper motion data from Simbad
def test_gaia_bh3_1(t: Time):
    name = "Gaia DR3 4318465066420528000"
    verbose(f"----- {name} test case 1 - Simbad -----")
    object_with_proper_motion(name, t)


# Proper motion data from Wikipedia article
def test_gaia_bh3_2(t: Time):
    name = "Gaia DR3 4318465066420528000"
    verbose(f"----- {name} test case 2 - Wikipedia -----")
    coord = SkyCoord.from_name(name)
    coord_pm = SkyCoord(ra=coord.ra, dec=coord.dec, 
                      # proper motion from Wikipedia https://en.wikipedia.org/wiki/Gaia_BH3
                      # distance must set, else apply_space_motion() will throw an "erfa" warning
                      distance=591*u.pc, radial_velocity=-333.2*u.km/u.s,
                      pm_ra_cosdec=-28.317*u.mas/u.yr, pm_dec=-155.221*u.mas/u.yr,
                      obstime="J2000") 
    coord_with_proper_motion(name, coord_pm, t)



def object_with_proper_motion(name: str, t: Time):
    coord = query_simbad(name)
    coord_with_proper_motion(name, coord, t)


def coord_with_proper_motion(name: str, coord: SkyCoord, t: Time):
    ic(coord)
    verbose(f"{name}: catalog position (ICRS) = {coord.to_string("hmsdms")}")
    verbose(f"{name}: proper motion dRA={coord.pm_ra_cosdec}, dDEC={coord.pm_dec}")

    # coord_j2000 = coord.fk5
    # ic(coord_j2000)
    # verbose(f"{name} catalog position (FK5 J2000) = {coord_j2000.to_string("hmsdms")}")

    ic(t)
    # coord_now = coord_j2000.apply_space_motion(new_obstime=t)
    coord_now = coord.apply_space_motion(new_obstime=t)
    ic(coord_now)
    verbose(f"{name}: position now (ICRS)     = {coord_now.to_string("hmsdms")}")



def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = "Calculate current position of object with proper motion",
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("--date", help="calculate position for DATE/time (default now)")
    arg.add_argument("--test-gaia-bh3", action="store_true", help="test with object Gaia BH3")
    arg.add_argument("object", nargs="*", help="object name (queried from simbad)")

    args = arg.parse_args()

    if args.debug:
        ic.enable()
        ic(sys.version_info)
        ic(args)
    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()

    if args.date:
        t = Time(args.date)
    else:
        t = Time.now()
    if args.test_gaia_bh3:
        verbose.enable()
        test_gaia_bh3_1(t)
        test_gaia_bh3_2(t)
    else:
        for obj in args.object:
            object_with_proper_motion(obj, t)



if __name__ == "__main__":
    main()
