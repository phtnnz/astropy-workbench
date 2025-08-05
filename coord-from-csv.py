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
# Version 0.1 / 2025-06-30
#       Coord conversion using positions from CSV file, based on coord-jnow.py 0.3

VERSION = "0.1 / 2025-06-30"
AUTHOR  = "Martin Junius"
NAME    = "coord-from-csv"

import sys
import argparse
import csv

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, CartesianRepresentation
from astropy.coordinates import ICRS, GCRS, PrecessedGeocentric, FK4, FK5, HADec  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
from astropy.coordinates import errors, get_constellation
import astropy.units as u
from astropy.time        import Time, TimeDelta
import numpy as np

# Local modules
from verbose import verbose, warning, error, message
from csvoutput import csv_output

from astroutils import ra_from_lst_ha, ra_dec_to_string, angle_to_string, hourangle_to_string, get_location



# Command line options
class Options:
    """
    Global options
    """
    csv = False         # -C --csv
    output = None       # -o --output
    set_locale = False  # -l --locale



def coord_to_jnow_altaz(ra: float, dec: float, loc: EarthLocation, date_obs: str) -> None:
    time = Time(date_obs, location=loc)
    ic(time)

    # FK5/J2000 coord
    coord = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame=FK5, equinox="J2000")

    # constellation = str(get_constellation(coord))
    # constellation_s = str(get_constellation(coord, short_name=True))
    # verbose(f"constellation {constellation} ({constellation_s})")

    verbose(f"FK5 J2000 coord         {ra_dec_to_string(coord.ra, coord.dec)}")

    # Transform to topocentric HADec
    hadec_jnow = coord.transform_to(HADec(obstime=time, location=loc))
    ic(hadec_jnow)
    
    lst_mean     = time.sidereal_time("mean")
    ra_mean      = ra_from_lst_ha(lst_mean, hadec_jnow.ha)
    ic(lst_mean, ra_mean)
    # verbose(f"LST={hourangle_to_string(lst_mean)}, HA={hourangle_to_string(hadec_jnow.ha)}")
    verbose(f"Topocentric JNOW coord  {ra_dec_to_string(ra_mean, hadec_jnow.dec)}")

    # Transform to topocentric AltAz
    altaz = coord.transform_to( AltAz(obstime=time, location=loc) )
    ic(altaz)
    verbose(f"Alt={angle_to_string(altaz.alt, decimal=True)} Az={angle_to_string(altaz.az, decimal=True)}")

    ra_j2000  = float(coord.ra.degree)
    dec_j2000 = float(coord.dec.degree)
    ra_jnow   = float(ra_mean.degree)
    dec_jnow  = float(hadec_jnow.dec.degree)
    alt       = float(altaz.alt.degree)
    az        = float(altaz.az.degree)

    # message(f"{date_obs=} {ra_j2000=}, {dec_j2000=}, {ra_jnow=}, {dec_jnow=}, {alt=}, {az=}")
    message(f"{date_obs=}, {ra_j2000=}, {dec_j2000=}, {alt=}, {az=}")

    if Options.csv:
        csv_output(fields=[ "date_obs", "ra_j2000", "dec_j2000", "alt", "az" ])
        csv_output(row=[ date_obs, ra_j2000, dec_j2000, alt, az ])



def process_csv(loc: EarthLocation, csvfile: str) -> None:
    verbose(f"processing CSV file {csvfile}")
    ic(loc)

    with open(csvfile) as file:
        reader = csv.DictReader(file)
        for row in reader:
            ic(row)
            # If WCS headers are present, use them
            ra  = row.get("CRVAL1") or row.get("RA")
            dec = row.get("CRVAL2") or row.get("DEC")
            if not ra or not dec:
                warning("headers CRVAL1/CRVAL2/RA/DEC not found")
                continue
            # UTC time stamp
            date_obs = row.get("DATE-OBS")
            if not date_obs:
                warning("header DATE-OBS not found")
                continue

            ic(date_obs, ra, dec)
            coord_to_jnow_altaz(ra, dec, loc, date_obs)



def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = "Astropy SkyCoord transformations from CSV file",
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("-L", "--location", help="named location or MPC station code")
    arg.add_argument("-C", "--csv", action="store_true", help="output CSV list")
    arg.add_argument("-l", "--locale", action="store_true", help="set locale for CSV output")
    arg.add_argument("-o", "--output", help="write CSV to file OUTPUT (default: stdout)")

    arg.add_argument("csvfile", nargs="+", help="input CSV file")

    args = arg.parse_args()

    if args.debug:
        ic.enable()
        ic(sys.version_info)
        ic(args)
    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()

    Options.csv = args.csv
    Options.output = args.output
    Options.set_locale = args.locale

    # Default location is Hakos, Namibia for personal reasons ;-)
    loc      = EarthLocation(lon=16.36167*u.deg , lat=-23.23639*u.deg, height=1853*u.m)
    if args.location:
        loc = get_location(args.location)
    ic(loc, loc.to_geodetic())

    for csvfile in args.csvfile:
        process_csv(loc, csvfile)

    if Options.csv:
        csv_output.set_float_format("%.6f")
        csv_output.write(Options.output, Options.set_locale)



if __name__ == "__main__":
    main()
