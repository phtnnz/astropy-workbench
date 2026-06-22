#!/usr/bin/env python

# Copyright 2025-2026 Martin Junius
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
# Version 0.1 / 2026-06-22
#       Query JPL SBWOBS and MPC dates of last observation

VERSION     = "0.1 / 2026-06-22"
AUTHOR      = "Martin Junius"
NAME        = "neo-sbwobs"
DESCRIPTION = "Retrieve observable NEOs/comets from JPL/MPC"

import sys
import argparse

from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
import astropy.units as u

# Local modules
from utils.verbose import verbose, warning, error, message
from neo.config import config
from neo.ephem import get_local_circumstances, get_dec_limits
from jpl.sbwobs import sbwobs_get_objects


def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = DESCRIPTION,
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("--asteroids", action="store_true", help=f"get asteroids default={config.sb_kind}")
    arg.add_argument("--neo", action="store_true", help=f"get NEOs default={config.sb_group}")
    arg.add_argument("--pha", action="store_true", help=f"get PHAs")
    arg.add_argument("--comets", action="store_true", help=f"get comets (overrides asteroid options)")
    arg.add_argument("-o", "--output", help="write object list to OUTPUT")
    arg.add_argument("-M", "--mag-limit", help="override mag_limit from config")
    arg.add_argument("-l", "--location", help=f"coordinates, named location or MPC station code, default {config.code}")
    arg.add_argument("--dln", action="store_true", help=f"use DLN list (default: DLU)")
    arg.add_argument("--lastobs", action="store_true", help=f"use LastObs list (default: DLU)")


    args = arg.parse_args()

    if args.debug:
        ic.enable()
        ic(sys.version_info, sys.path, args)
    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()

    # override defaults from config
    if args.mag_limit:
        config.mag_limit = float(args.mag_limit)
        config.vmag_max  = float(args.mag_limit)
    if args.asteroids:
        config.sb_kind = "a"
    if args.neo:
        config.sb_group = "neo"
    if args.pha:
        config.sb_group = "pha"
    if args.comets:
        config.sb_kind = "c"
        config.sb_group = None

    # Observer location and local circumstances
    local = get_local_circumstances(args.location if args.location else config.code)

    # Update DEC limits in config
    min_dec, max_dec = get_dec_limits(local, config.min_alt*u.deg)
    config.min_dec = int(min_dec.degree)
    config.max_dec = int(max_dec.degree)

    # MPC list type
    list_type = "DLU"
    if args.dln:
        list_type = "DLN"
    if args.lastobs:
        list_type = "LASTOBS"
    keys_selected = sbwobs_get_objects(local, list_type)

    if args.output:
        with open(args.output, "w") as file:
            file.writelines([line + "\n" for line in keys_selected])



if __name__ == "__main__":
    main()
