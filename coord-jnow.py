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
# Version 0.1 / 2024-12-12
#       SkyCoord transformations
# Version 0.2 / 2025-01-08
#       Output ICRS, FK5/J2000, FK4/B1950, FK5/JNOW, GCRS, topocentrc LST, 
#       hourangle, parallactic angle, JNOW, AltAz, and with refraction
#       New options -t --time, -j --j2000, -q --query-simbad
#       Query Simbad for object

import sys
import argparse

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, CartesianRepresentation
from astropy.coordinates import ICRS, GCRS, FK4, FK5, HADec  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
import astropy.units as u
from astropy.time        import Time, TimeDelta
import numpy as np

# Local modules
from verbose import verbose, warning, error
from querysimbad import query_simbad

VERSION = "0.2 / 2025-01-08"
AUTHOR  = "Martin Junius"
NAME    = "coord-jnow"



# Command line options
class Options:
    j2000 = False           # -j --j2000
    query_simbad = False    # -q --query-simbad



def test_sn2024abfo():
    ### Test case with position of SN 2024abfo in comparison to NINA / Autoslew ###

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

    obj1     = "03:57:25.611 -46:11:07.57" # SN 2024abfo precise position
    obj      = "03:57:25.6 -46:11:07.6" # SN 2024abfo NINA sequencer position
    jnow     = "03:58:15 -46:06:45"    # NINA Epoch: JNOW - Alt: 51° 00' 18"; Az: 135° 44' 08"
    altaz    = "51:00:18 135:44:08"
    autoslew = "03h58m14.52s -46d06m45s"   # Autoslew position with Ep: "real"
    loc      = EarthLocation(lat=-23.23639*u.deg, lon=16.36167*u.deg , height=1825*u.m) # Hakos, Namibia
    time     = Time("2024-11-22 20:06:20", location=loc)
    ic(obj, jnow, autoslew, loc, time)
    verbose(f"NINA JNOW: {jnow}")
    verbose(f"NINA AltAz: {altaz}")
    coord_to_jnow_altaz(obj, loc, time)



def ra_from_lst_ha(lst: Angle, ha: Angle):
    ra = lst - ha
    if ra >= 24*u.hourangle:
        ra -= 24*u.hourangle
    ic(lst, ha, ra)
    return ra


def ra_dec_to_string(ra: Angle, dec: Angle):
    return f"RA={ra.to_string(unit=u.hour, precision=2)} DEC={dec.to_string(unit=u.degree, precision=2)}"

def angle_to_string(a: Angle):
    return f"{a.to_string(unit=u.degree, precision=2)}"

def hourangle_to_string(a: Angle):
    return f"{a.to_string(unit=u.hour, precision=2)}"


def coord_to_jnow_altaz(obj: str, loc: EarthLocation, time: Time):
    if Options.query_simbad:
        # Query object name
        coord = query_simbad(obj)
        ic(coord)
        verbose(f"ICRS coord {ra_dec_to_string(coord.ra, coord.dec)}")
        coord_j2000 = coord.transform_to(FK5(equinox="J2000"))
    elif Options.j2000:
        # FK5/J2000 coord
        coord = SkyCoord(obj, unit=(u.hour, u.deg), frame=FK5, equinox="J2000")
        coord_j2000 = coord
    else:
        # ICRS coord
        coord = SkyCoord(obj, unit=(u.hour, u.deg))
        ic(coord)
        verbose(f"ICRS coord {ra_dec_to_string(coord.ra, coord.dec)}")
        coord_j2000 = coord.transform_to(FK5(equinox="J2000"))

    ic(coord_j2000)
    verbose(f"FK5 J2000 coord {ra_dec_to_string(coord_j2000.ra, coord_j2000.dec)}")

    coord_b1950 = coord.transform_to(FK4(equinox="B1950"))
    ic(coord_b1950)
    verbose(f"FK4 B1950 coord {ra_dec_to_string(coord_b1950.ra, coord_b1950.dec)}")

    coord_jnow = coord.transform_to(FK5(equinox=Time(time, format="jyear")))
    ic(coord_jnow)
    verbose(f"FK5 JNOW coord {ra_dec_to_string(coord_jnow.ra, coord_jnow.dec)}")

    # Transform to geocentric GCRS
    coord_gcrs = coord.transform_to(GCRS(obstime=time))
    ic(coord_gcrs)
    verbose(f"GCRS coord {ra_dec_to_string(coord_gcrs.ra, coord_gcrs.dec)}")

    # Transform to topocentric HADec
    hadec_jnow = coord.transform_to(HADec(obstime=time, location=loc))
    hadec_jnow_w_refraction = coord.transform_to(HADec(obstime=time, location=loc, 
                                                 pressure=1000*u.hPa, temperature=20*u.deg_C,
                                                 relative_humidity=0.4, obswl=0.54*u.micron))
    ic(hadec_jnow, hadec_jnow_w_refraction)
    
    lst_mean     = time.sidereal_time("mean")
    lst_apparent = time.sidereal_time("apparent")
    ra_mean      = ra_from_lst_ha(lst_mean, hadec_jnow.ha)
    ra_apparent  = ra_from_lst_ha(lst_apparent, hadec_jnow.ha)
    ic(lst_mean, lst_apparent, ra_mean, ra_apparent)
    verbose(f"Topocentric mean LST={hourangle_to_string(lst_mean)}, hourangle={hourangle_to_string(hadec_jnow.ha)}")
    verbose(f"Topocentric JNOW coord (lst mean) {ra_dec_to_string(ra_mean, hadec_jnow.dec)}")
    # verbose(f"Topocentric JNOW coord (lst apparent) {ra_dec_to_string(ra_apparent, hadec_jnow.dec)}")

    ra_mean_w_refraction = ra_from_lst_ha(lst_mean, hadec_jnow_w_refraction.ha)
    ra_apparent_w_refraction = ra_from_lst_ha(lst_apparent, hadec_jnow_w_refraction.ha)
    ic(ra_mean_w_refraction, ra_apparent_w_refraction)
    verbose(f"Topocentric JNOW coord (lst mean) w/REFRACTION {ra_dec_to_string(
        ra_mean_w_refraction, hadec_jnow_w_refraction.dec)}")
    # verbose(f"Topocentric JNOW coord (lst apparent) w/REFRACTION {ra_dec_to_string(
    #     ra_apparent_w_refraction, hadec_jnow_w_refraction.dec)}")

    # Transform to topocentric AltAz
    altaz = hadec_jnow.transform_to( AltAz(obstime=time, location=loc) )
    altaz_w_refraction = hadec_jnow_w_refraction.transform_to( AltAz(obstime=time, location=loc) )
    ic(altaz, altaz_w_refraction)
    verbose(f"Topocentric Alt={angle_to_string(altaz.alt)} Az={angle_to_string(altaz.az)}")
    verbose(f"Topocentric w/REFRACTION Alt={angle_to_string(altaz_w_refraction.alt)} Az={angle_to_string(altaz_w_refraction.az)}")

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
    verbose(f"Parallactic angle={angle_to_string(q)}")



def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = "Astropy SkyCoord transformations",
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("--test-sn2024abfo", action="store_true", help="test case SN 2024abfo")
    arg.add_argument("-j", "--j2000", action="store_true", help="object coordinates FK5/J2000, default ICRS")
    arg.add_argument("-t", "--time", help="time (UTC) for JNOW coordinates, default now")
    arg.add_argument("-q", "--query-simbad", action="store_true", help="query Simbad for OBJECT name")
    arg.add_argument("object", nargs="*", help="sky coord \"RA DEC\" or object name (-q)")

    args = arg.parse_args()

    if args.debug:
        ic.enable()
        ic(sys.version_info)
        ic(args)
    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()

    ##FIXME: loc / time from command line
    loc      = EarthLocation(lat=-23.23639*u.deg, lon=16.36167*u.deg , height=1825*u.m) # Hakos, Namibia

    Options.j2000 = args.j2000
    Options.query_simbad = args.query_simbad

    if args.time:
        time = Time(args.time, location=loc)
    else:
        time = Time(Time.now(), location=loc)
    ic(loc, time)

    if args.test_sn2024abfo:
        verbose.enable()
        test_sn2024abfo()
    elif args.object:
        for obj in args.object:
            verbose(f"object {obj}")
            verbose(f"time (UTC) {time}")
            coord_to_jnow_altaz(obj, loc, time)            



if __name__ == "__main__":
    main()
