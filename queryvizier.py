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

# Astroquery
from astroquery.vizier import Vizier

# Local modules
from verbose import verbose, warning, error

VERSION = "0.1 / 2025-01-05"
AUTHOR  = "Martin Junius"
NAME    = "queryvizier"



def query_vizier(cat: str):
    # vizier = Vizier(catalog=cat)
    vizier = Vizier(catalog=cat,
                    columns=["*", 'Names'] ) ##FIXME: get from command line
    meta = vizier.get_catalog_metadata()
    ic(vizier, meta, meta.keys())
    result = vizier.query_object("")        # query all catalog entries
    ic(result)
    for name in result.keys():
        table = result[name]
        ic(table)
    pass



def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = "Query VizieR catalog",
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("catalog", help="catalog name")

    args = arg.parse_args()

    if args.debug:
        ic.enable()
        ic(args, sys.version_info, sys.path)
    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()
    cat = args.catalog
    verbose(f"query catalog {cat}")
    query_vizier(cat)



if __name__ == "__main__":
    main()