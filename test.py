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
# Version 0.0 / 2024-11-17
#       Test with SkyCoord

import sys
import argparse
from datetime import datetime, UTC

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.coordinates import ICRS, ITRS, Galactic, FK4, FK5, HADec, CIRS  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
import astropy.units as u
from astropy.time        import Time, TimeDelta

# Local modules
from verbose import verbose, warning, error

VERSION = "0.0 / 2024-xx-xx"
AUTHOR  = "Martin Junius"
NAME    = "astropy-template"



# Command line options
class Options:
    name = "abc"        # -n --name
    int  = 99           # -i --int



def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = "Astropy SkyCoord tests",
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("object", nargs="*", help="sky coord \"RA DEC\"")

    args = arg.parse_args()

    if args.debug:
        ic.enable()
        ic(sys.version_info)
        ic(args)
    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()

    if args.object:
        print("nothing happens ... ;-)")
        return
        for obj in args.object:
            pass

    # NINA:
    #
    # 2024-11-22T22:06:20.8901|INFO|SequenceItem.cs|
    #   Run|208|Starting Category: Telescope, Item: Center, Coordinates RA: 03:57:26; Dec: -46° 11' 08"; Epoch: J2000
    # 2024-11-22T22:06:20.9889|INFO|TelescopeVM.cs|
    #   SlewToCoordinatesAsync|926|Slewing from RA: 13:10:34; Dec: -82° 10' 33"; Epoch: JNOW to RA: 03:58:15; Dec: -46° 06' 45"; Epoch: JNOW - Alt: 51° 00' 18"; Az: 135° 44' 08"
    # 2024-11-22T22:06:58.2396|INFO|ImageSolver.cs|
    #   Solve|54|Platesolve successful: Coordinates: RA: 03:56:49; Dec: -46° 10' 35"; Epoch: J2000
    # 2024-11-22T22:06:58.2445|INFO|CenteringSolver.cs|
    #   Center|99|Centering Solver - Scope Position: RA: 03:58:15; Dec: -46° 06' 45"; Epoch: JNOW; Offset: RA: 00:00:00; Dec: 00° 00' 00"; Distance: 00° 00' 00"; Bearing: 00° 00' 00"; Centering Coordinates: RA: 03:58:15; Dec: -46° 06' 45"; Epoch: JNOW; Solved: RA: 03:57:38; Dec: -46° 06' 12"; Epoch: JNOW; Separation RA: 00:00:36; Dec: -00° 00' 33"; Distance: 00° 06' 17"; Bearing: 84° 58' 06"; Threshold: 0.5
    #
    # 2024-11-22T22:06:59.2706|INFO|TelescopeVM.cs|
    #   Sync|800|Syncing scope from RA: 03:58:15; Dec: -46° 06' 45"; Epoch: JNOW to RA: 03:57:38; Dec: -46° 06' 12"; Epoch: JNOW
    # 2024-11-22T22:07:09.2757|INFO|CenteringSolver.cs|
    #   Center|132|Slewing to target after sync. Current Position: RA: 03:57:38; Dec: -46° 06' 12"; Epoch: JNOW; Target coordinates: RA: 03:58:15; Dec: -46° 06' 45"; Epoch: JNOW; Offset RA: 00:00:00; Dec: 00° 00' 00"; Distance: 00° 00' 00"; Bearing: 00° 00' 00"
    # 2024-11-22T22:07:09.2790|INFO|TelescopeVM.cs|
    #   SlewToCoordinatesAsync|926|Slewing from RA: 03:57:38; Dec: -46° 06' 12"; Epoch: JNOW to RA: 03:58:15; Dec: -46° 06' 45"; Epoch: JNOW - Alt: 51° 08' 03"; Az: 135° 49' 13"
    # 2024-11-22T22:07:30.8207|INFO|ImageSolver.cs|
    #   Solve|54|Platesolve successful: Coordinates: RA: 03:57:26; Dec: -46° 11' 07"; Epoch: J2000
    # 2024-11-22T22:07:30.8318|INFO|CenteringSolver.cs|
    #   Center|99|Centering Solver - Scope Position: RA: 03:58:15; Dec: -46° 06' 45"; Epoch: JNOW; Offset: RA: 00:00:00; Dec: 00° 00' 00"; Distance: 00° 00' 00"; Bearing: 00° 00' 00"; Centering Coordinates: RA: 03:58:15; Dec: -46° 06' 45"; Epoch: JNOW; Solved: RA: 03:58:15; Dec: -46° 06' 44"; Epoch: JNOW; Separation RA: -00:00:00; Dec: -00° 00' 01"; Distance: 00° 00' 02"; Bearing: -68° 54' 32"; Threshold: 0.5
    #
    # Autoslew:
    # 22:06:21.084: Interface SlewToTargetAsynch
    # 22:06:21.084: Starting Slew to RA 03.97h  DE -46.11d  EPOCH 0000
    # 22:06:58.259: Interface setting TargetRightAscension to 3.96 False
    # 22:06:58.260: Interface setting TargetDeclination to -46.10 True
    # 22:06:58.260: Interface SynchToTarget
    # 22:07:09.282: Interface setting TargetRightAscension to 3.97 False
    # 22:07:09.283: Interface setting TargetDeclination to -46.11 True
    # 22:07:09.283: Interface SlewToTargetAsynch
    # 22:07:09.283: Starting Slew to RA 03.97h  DE -46.11d  EPOCH 0000


    obj = "03:57:25.611 -46:11:07.57" # SN 2024abfo

    coord = SkyCoord(obj, unit=(u.hour, u.deg))
    ic(coord)
    print(f"coord ICRS      {coord.to_string("hmsdms")}")
    print(f"Autoslew M49    03h58m14.52s       -46d06m45s")
    print()

    # From https://www.cloudynights.com/topic/861776-asiair-small-program-to-convert-from-j2000-to-jnow/
    # _fk5 = FK5(equinox=Time(Time(datetime.now(UTC), scale='utc').jd, format="jd", scale="utc"))
    # coord_now = SkyCoord(object, frame=FK4, unit=(u.hourangle, u.deg)).transform_to(_fk5)
    # ic(_fk5, coord_now)
    # print(f"coord now={coord_now.to_string("hmsdms")}")

    # ic(coord_now.ra - coord.ra, coord_now.dec - coord.dec)

    coord1 = coord.transform_to(FK5(equinox="J2000"))
    ic(coord1)
    print(f"coord FK5 J2000 {coord1.to_string("hmsdms")}")

    coord2 = coord.transform_to(FK5(equinox=f"J{Time.now().jyear}"))
    ic(coord2)
    print(f"coord FK5 now   {coord2.to_string("hmsdms")}")

    loc  = EarthLocation(lat=-23.23639*u.deg, lon=16.36167*u.deg , height=1825*u.m)
    time = Time(Time.now(), location=loc)

    hadec1 = coord.transform_to(HADec(obstime=time, location=loc))
    hadec2 = coord.transform_to(HADec(obstime=time, location=loc, 
                                     pressure=1000*u.hPa, temperature=20*u.deg_C,
                                     relative_humidity=0.4, obswl=0.54*u.micron))
    lst  = time.sidereal_time("mean")
    ra1 = lst - hadec1.ha
    if ra1 >= 24*u.hourangle:
        ra1 -= 24*u.hourangle
    ra2 = lst - hadec2.ha
    if ra2 >= 24*u.hourangle:
        ra2 -= 24*u.hourangle
    ic(loc, time, lst, ra1, ra2)
    ic(hadec1, ra1.to_string(unit=u.hour), hadec1.dec.to_string(unit=u.degree))
    ic(hadec2, ra2.to_string(unit=u.hour), hadec2.dec.to_string(unit=u.degree))
    print(f"lst {lst.to_string(unit=u.hour)}")
    print(f"ha {hadec1.ha.to_string(unit=u.degree)}")
    print(f"hadec1 now       {ra1.to_string(unit=u.hour)} {hadec1.dec.to_string(unit=u.degree)}")
    print(f"hadec2 now refr  {ra2.to_string(unit=u.hour)} {hadec2.dec.to_string(unit=u.degree)}")


if __name__ == "__main__":
    main()
