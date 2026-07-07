#!/usr/bin/env python

# Copyright 2026 Martin Junius
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
# Version 0.1 / 2026-06-21
#       Query previous NEOCP objects
# Version 0.2 / 2026-06-23
#       Added -f --file option to read from plan CSV file

VERSION     = "0.2 / 2026-06-23"
AUTHOR      = "Martin Junius"
NAME        = "neo-list-prev"
DESCRIPTION = "Query previous NEOCP list from MPC"

import sys
import argparse
import csv

from icecream import ic
# Disable debugging
ic.disable()

# Local modules
from utils.verbose import verbose, warning, error, message
from neo.config import config
from mpc.prevneocp import mpc_query_prev_neocp, mpc_parse_prev_neocp



def objects_from_csv(csvfile: str) -> list[str]:
    verbose(f"processing CSV file {csvfile}")
    objects = list()

    with open(csvfile) as file:
        reader = csv.DictReader(file)
        for row in reader:
            # Make field names lower case
            row = { key.lower():val for key, val in row.items() if key}
            ic(row)

            # Get object
            target = row.get("object") or row.get("name") or row.get("target")
            if not target:
                error("can't find target name in CSV data")
            obj_type = row.get("type")

            if not obj_type or obj_type=="NEOCP" or obj_type=="PCCP":
                objects.append(target)

    return objects



def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = DESCRIPTION,
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("-f", "--file", help="read objects from CSV FILE")
    arg.add_argument("object", nargs="*", help="object name")

    args = arg.parse_args()

    if args.debug:
        ic.enable()
        ic(sys.version_info, sys.path, args)
    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()

    content = mpc_query_prev_neocp(config.prev_neocp_url)
    prev_objects = mpc_parse_prev_neocp(content)

    objects = list()
    if args.file:
        objects.extend(objects_from_csv(args.file))
    if args.object:
        objects.extend(args.object)

    if objects:    
        for obj in objects:
            if obj in prev_objects:
                ic(prev_objects[obj])
                message(prev_objects[obj])
            else:
                message(f"{obj:7s}   not in previous NEOCP list")
    else:
        for obj in sorted(prev_objects.keys()):
            message(prev_objects[obj])



if __name__ == "__main__":
    main()
