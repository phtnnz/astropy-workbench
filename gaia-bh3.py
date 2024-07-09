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
# Version 0.0 / 2024-07-09
#       SkyCoord handling of Gaia BH3 with proper motion

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
import astropy.units as u
from astropy.time        import Time, TimeDelta

# Local modules
from verbose import verbose, warning, error

VERSION = "0.0 / 2024-07-09"
AUTHOR  = "Martin Junius"
NAME    = "gaia-bh3"



# Command line options
class Options:
    pass



def gaia_bh3():
    name = "Gaia DR3 4318465066420528000"
    bh3 = SkyCoord.from_name(name)
    ic(bh3)
    bh3_pm = SkyCoord(ra=bh3.ra, dec=bh3.dec, 
                      # proper motion from Wikipedia https://en.wikipedia.org/wiki/Gaia_BH3
                      distance=591*u.pc, radial_velocity=-333.2*u.km/u.s,
                      pm_ra_cosdec=-28.317*u.mas/u.yr, pm_dec=-155.221*u.mas/u.yr,
                      obstime="J2000") 
    ic(bh3_pm)
    verbose(f"Gaia BH3 ({name}) catalog position (ICRS) = {bh3_pm.to_string("hmsdms")}")
 
    now = Time.now()
    ic(now)
    bh3_now = bh3_pm.apply_space_motion(new_obstime=now)
    ic(bh3_now)
    verbose(f"Gaia BH3 ({name}) position now (ICRS) = {bh3_now.to_string("hmsdms")}")



def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = "Calculate current position of Gaia BH3",
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")

    args = arg.parse_args()

    if args.debug:
        ic.enable()
        ic(sys.version_info)
        ic(args)
    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()
    # ... more options ...
        
    # ... the action starts here ...
    gaia_bh3()



if __name__ == "__main__":
    main()