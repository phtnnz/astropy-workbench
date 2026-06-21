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

VERSION     = "0.1 / 2026-06-21"
AUTHOR      = "Martin Junius"
NAME        = "neo-list-prev"
DESCRIPTION = "Query previous NEOCP list from MPC"

import sys
import argparse
from icecream import ic
# Disable debugging
ic.disable()

# Local modules
from utils.verbose import verbose, warning, error, message
from neo.config import config
from mpc.prevneocp import mpc_query_prev_neocp, mpc_parse_prev_neocp



def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = DESCRIPTION,
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("object", nargs="*", help="object name")

    args = arg.parse_args()

    if args.debug:
        ic.enable()
        ic(sys.version_info, sys.path, args)
    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()

    content = mpc_query_prev_neocp(config.prev_neocp_url)
    objects = mpc_parse_prev_neocp(content)

    if args.object:
        for obj in args.object:
            if obj in objects:
                ic(objects[obj])
                verbose(objects[obj])
            else:
                verbose(f"{obj:7s}   not in previous NEOCP list")
    else:
        for obj in sorted(objects.keys()):
            verbose(objects[obj])



if __name__ == "__main__":
    main()
