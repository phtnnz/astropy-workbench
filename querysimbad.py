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
# Version 0.1 / 2024-07-14
#       Query Simbad for object coordinates, distance, proper motion, velocity
# Version 0.2 / 2024-07-29
#       Renamed, can be used as a module
# Version 0.3 / 2024-12-16
#       Added docstrings
# Version 0.4 / 2025-01-30
#       Major rework for new table format

import sys
import argparse

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.coordinates import SkyCoord  # High-level coordinates
import astropy.units as u
import numpy as np

# Astroquery
from astroquery.simbad import Simbad

# Local modules
from verbose import verbose, warning, error

VERSION = "0.4 / 2025-01-30"
AUTHOR  = "Martin Junius"
NAME    = "querysimbad"



def query_simbad(obj: str, w_velocity: bool=True) -> SkyCoord:
    """
    Query Simbad and return SkyCoord with position, distance, proper motion, radial velocity

    :param obj: object name
    :type obj: str
    :param w_velocity: flag to include distance/velocity
    :type obj: bool, default: True
    :return: coordinates of object
    :rtype: SkyCoord
    """
    verbose(f"query object {obj}")

    simbad = Simbad()
    simbad.add_votable_fields("mesdistance", "propermotions", "velocity")
    ic(simbad.get_votable_fields())

    result = simbad.query_object(obj)
    ic(type(result))
    ic(result)
    ic(result.info)

    # Object not found
    if not result:
        return None
    
    id       = result["main_id"][0]
    id_u     = result["main_id"].unit
    ra       = result["ra"][0]
    dec      = result["dec"][0]
    ra_u     = result["ra"].unit
    dec_u    = result["dec"].unit
    ic(id, ra, dec)
    ic(id_u, ra_u, dec_u)

    dist     = result["mesdistance.dist"][0]
    dist_u   = result["mesdistance.unit"][0]
    # No proper unit in the.unit column, work-around
    if not dist_u or dist_u.strip() == "pc":
        dist_u = u.pc
    elif dist_u.strip() == "Mpc":
        dist_u = u.Mpc
    elif dist_u.strip() == "kpc":
        dist_u = u.kpc
    else:
        error(f"query_simbad: unknown distance unit {dist_u}")
    ic(dist, dist_u)

    pmra     = result["pmra"][0]
    pmdec    = result["pmdec"][0]
    radvel   = result["rvz_radvel"][0]
    pmra_u   = result["pmra"].unit
    pmdec_u  = result["pmdec"].unit
    radvel_u = result["rvz_radvel"].unit
    ic(pmra, pmra_u, pmdec, pmdec_u, radvel, radvel_u)

    if np.ma.is_masked(pmra):
        ic("no propermotions/velocity data")
        w_velocity = False

    if w_velocity:
        coord = SkyCoord(ra=ra*ra_u, dec=dec*dec_u,
                        distance=dist*dist_u, radial_velocity=radvel*radvel_u,
                        pm_ra_cosdec=pmra*pmra_u, pm_dec=pmdec*pmdec_u,
                        obstime="J2000") 
    else:
        coord = SkyCoord(ra=ra*ra_u, dec=dec*dec_u, distance=dist*dist_u) 

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

    for obj in args.object:
        coord = query_simbad(obj)
        verbose(f"{obj}: position (ICRS) = {coord.to_string("hmsdms")}")




if __name__ == "__main__":
    main()