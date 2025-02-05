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
# Version 0.1 / 2025-01-05
#       Query VizieR catalog
# Version 0.2 / 2025-01-06
#       New option, -f --ra-dec-float to output RA/DEC as float degrees,
#       --replace-comma to replace "," in fields with the given string
# Version 0.3 / 2025-01-27
#       Added --object option, not that useful, because .query_object()
#       doesn't seem to support wildcards * or %
#       Added -m --match (regex) option
#       Added -l --locale option for CSV output
# Version 0.4 / 2025-02-05
#       Added --constellation option, adds constellation full and short name
#       to CSV output

VERSION = "0.4 / 2025-02-05"
AUTHOR  = "Martin Junius"
NAME    = "queryvizier"

import sys
import argparse
import re

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.coordinates import SkyCoord, get_constellation
import astropy.units as u
from astropy.table import Table, Row

# Astroquery
from astroquery.vizier import Vizier

# Local modules
from verbose import verbose, warning, error
from csvoutput import csv_output



class Options:
    row_limit = -1              # -n --row-limit
    ra_dec_float = False        # -f --ra-dec-float
    replace_comma = None        # --replace-comma
    object = ""                 # -O --object, "" = query all catalog entries
    match = None                # -m --match
    constellation = False       # --constellation



def _replace_ra_dec(val, col: str, ra: float, dec: float):
    if Options.ra_dec_float:
        if col == "RAJ2000": return ra
        if col == "DEJ2000": return dec
    if Options.replace_comma and isinstance(val, str):
        ic(val, type(val))
        return val.replace(",", Options.replace_comma)
    return val


def query_vizier(cat: str, cols: list=None):
    # Process TABLE=OUT in cols list
    if cols:
        table_cols = [ col.split("=")[0] if "=" in col else col for col in cols ]
        out_cols   = [ col.split("=")[1] if "=" in col else col for col in cols ]
        ic(table_cols, out_cols)

    # Query VizieR
    vizier = Vizier(catalog=cat,
                    columns=["**"],                     # "*" = default columns, "**" = all columns
                    row_limit = Options.row_limit
                   )
    meta = vizier.get_catalog_metadata()
    ic(vizier, meta)
    verbose(f"catalog title: {meta["title"][0]}")
    verbose(f"      authors: {meta["authors"][0]}")

    result = vizier.query_object(Options.object)
    ic(result)
    for name in result.keys():
        table = result[name]
        if not cols:
            table_cols = out_cols = table.keys()
        if Options.constellation:
            out_cols.extend(["constellation", "constellation_short"])
        ic(table, table_cols, out_cols)
        verbose(f"table columns = {table_cols}")
        verbose(f"output columns = {out_cols}")
        csv_output(fields=out_cols)
        for row in table:
            ic(row)
            # Match string
            if Options.match:
                if not any(isinstance(field, str) and re.search(Options.match, field) for field in row):
                    continue
            # Special handling for coordinates
            coord = SkyCoord(row["RAJ2000"], row["DEJ2000"], unit=(u.hour, u.degree))
            ra = float(coord.ra.degree)
            dec = float(coord.dec.degree)
            ic(coord, ra, dec)
            values = [ _replace_ra_dec(row[col], col, ra, dec)
                       for col in table_cols ]          # Selected columns
            if Options.constellation:
                values.append(str(get_constellation(coord)))
                values.append(str(get_constellation(coord, short_name=True)))
            ic(values)
            csv_output(row=values)



def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = "Query VizieR catalog",
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("--columns", help="columns to retrieve, comma-separated, default all")
    arg.add_argument("-n", "--row-limit", help="number of rows to retrieve, default unlimited")
    arg.add_argument("-f", "--ra-dec-float", action="store_true", help="output RA/DEC as float degrees")
    arg.add_argument("-C", "--csv", action="store_true", help="CSV output")
    arg.add_argument("-l", "--locale", action="store_true", help="set locale for CSV output")
    arg.add_argument("-o", "--output", help="output file")
    arg.add_argument("--replace-comma", help="replace \",\" in field with REPLACE_COMMA")
    arg.add_argument("--object", help="query specific object, default \"\"=all, no wildcards")
    arg.add_argument("-m", "--match", help="output rows containing regex MATCH only")
    arg.add_argument("--constellation", action="store_true", help="get and output constellation as an extra column")

    arg.add_argument("catalog", help="catalog name")

    args = arg.parse_args()

    if args.debug:
        ic.enable()
        ic(args, sys.version_info, sys.path)
    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()

    Options.ra_dec_float = args.ra_dec_float
    
    if args.row_limit:
        Options.row_limit = int(args.row_limit)
    Options.replace_comma = args.replace_comma
    columns = None
    if args.columns:
        columns = args.columns.split(",")
    if args.object:
        Options.object = args.object
    Options.match = args.match
    Options.constellation = args.constellation

    verbose(f"query catalog {args.catalog}")
    query_vizier(args.catalog, columns)

    if args.csv:
        csv_output.write(args.output, set_locale=args.locale)


if __name__ == "__main__":
    main()
