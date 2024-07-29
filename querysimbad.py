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
# Version 0.1 / 2024-07-14
#       Query Simbad for object coordinates, distance, proper motion, velocity
# Version 0.2 / 2024-07-29
#       Renamed, can be used as a module

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

# Astroquery
from astroquery.simbad import Simbad

# Local modules
from verbose import verbose, warning, error

VERSION = "0.2 / 2024-07-29"
AUTHOR  = "Martin Junius"
NAME    = "querysimbad"



def query_simbad(obj):
    verbose(f"query object {obj}")

    simbad = Simbad()
    simbad.add_votable_fields("distance", "propermotions", "velocity")
    ic(simbad.get_votable_fields())

    result = simbad.query_object(obj)
    ic(type(result))
    ic(result)
    ic(result.info)

    id       = result["MAIN_ID"][0]
    id_u     = result["MAIN_ID"].unit
    ra       = result["RA"][0]
    dec      = result["DEC"][0]
    ra_u     = result["RA"].unit
    dec_u    = result["DEC"].unit
    if str(ra_u) != '"h:m:s"' or str(dec_u) != '"d:m:s"':
        error(f"query_simbad: RA/DEC units ({str(ra_u)}/{str(dec_u)}) must be hour/degree (h:m:s/d:m:s)")
    dist     = result["Distance_distance"][0]
    dist_u   = result["Distance_unit"][0]
    # No proper unit in Distance_distance column, work-around
    if not dist_u or dist_u == "pc":
        dist_u = u.pc
    else:
        error(f"query_simbad: unknown distance unit {dist_u}")
    pmra     = result["PMRA"][0]
    pmdec    = result["PMDEC"][0]
    radvel   = result["RVZ_RADVEL"][0]
    pmra_u   = result["PMRA"].unit
    pmdec_u  = result["PMDEC"].unit
    radvel_u = result["RVZ_RADVEL"].unit
    ic(id, id_u, ra, ra_u, dec, dec_u, dist, dist_u, pmra, pmra_u, pmdec, pmdec_u, radvel, radvel_u)

    coord = SkyCoord(ra=ra, dec=dec, unit=(u.hour, u.deg),
                    distance=dist*dist_u, radial_velocity=-radvel*radvel_u,
                    pm_ra_cosdec=pmra*pmra_u, pm_dec=pmdec*pmdec_u,
                    obstime="J2000") 
    ic(coord, coord.to_string("hmsdms"))

    return coord



def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = "Query Simbad for OBJECT coordinates, distance, proper motion, velocity",
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("object", nargs="+", help="object name")

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
    for obj in args.object:
        query_simbad(obj)



if __name__ == "__main__":
    main()