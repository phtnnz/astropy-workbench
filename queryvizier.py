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
# Version 0.1 / 2025-01-05
#       Query VizieR catalog

import sys
import argparse

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.coordinates import SkyCoord  # High-level coordinates
import astropy.units as u
from astropy.table import Table, Row

# Astroquery
from astroquery.vizier import Vizier

# Local modules
from verbose import verbose, warning, error
from csvoutput import csv_output

VERSION = "0.1 / 2025-01-05"
AUTHOR  = "Martin Junius"
NAME    = "queryvizier"



def query_vizier(cat: str, cols: list=None, row_limit: int=-1):
    # Process TABLE=OUT in cols list
    if cols:
        table_cols = [ col.split("=")[0] if "=" in col else col for col in cols ]
        out_cols   = [ col.split("=")[1] if "=" in col else col for col in cols ]
        ic(table_cols, out_cols)

    # Query VizieR
    vizier = Vizier(catalog=cat,
                    columns=["**"],                     # "*" = default columns, "**" = all columns
                    row_limit = row_limit
                   )
    meta = vizier.get_catalog_metadata()
    ic(vizier, meta)
    result = vizier.query_object("")                    # "" = query all catalog entries
    ic(result)
    for name in result.keys():
        table = result[name]
        if not cols:
            table_cols = out_cols = table.keys()
        ic(table, table_cols, out_cols)
        verbose(f"table columns = {table_cols}")
        verbose(f"output columns = {out_cols}")
        csv_output(fields=out_cols)
        for row in table:
            ic(row)
            # values = list(row)                        # Convert row to list
            values = [ row[col] for col in table_cols ] # Selected columns
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
    arg.add_argument("--row-limit", help="number of rows to retrieve, default unlimited")
    arg.add_argument("-C", "--csv", action="store_true", help="CSV output")
    arg.add_argument("-o", "--output", help="output file")

    arg.add_argument("catalog", help="catalog name")

    args = arg.parse_args()

    if args.debug:
        ic.enable()
        ic(args, sys.version_info, sys.path)
    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()

    cat = args.catalog
    columns = None
    row_limit = int(args.row_limit or -1)
    if args.columns:
        columns = args.columns.split(",")

    verbose(f"query catalog {cat}, {row_limit=}")
    query_vizier(cat, columns, row_limit)

    if args.csv:
        csv_output.write(args.output, set_locale=False)


if __name__ == "__main__":
    main()