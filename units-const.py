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
# Version 0.1 / 2025-02-10
#       Test with astropy.units und .constants

VERSION = "0.1 / 2025-02-10"
AUTHOR  = "Martin Junius"
NAME    = "units-const"

import sys
import argparse
import numpy as np

# AstroPy
import astropy.constants as const
import astropy.units as u

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()
# Local modules
from verbose import message, verbose, warning, error



def test():
    message("c\n", const.c, sep="")
    message("eps0\n", const.eps0, sep="")
    message("mu0\n", const.mu0, sep="")

    c2_1 = const.c ** 2
    message("c^2  = ", c2_1)

    c2_2 = 1 / (const.eps0 * const.mu0)
    message("1 / (eps0 * mu0)  = ", c2_2)
    c2_2a = c2_2.to(u.m ** 2 / u.s ** 2)
    message("1 / (eps0 * mu0)  = ", c2_2a)

    eps0 = const.eps0.to(u.A * u.s / (u.V * u.m))
    mu0  = const.mu0.to(u.V * u.s / (u.A * u.m))
    c2_2b = 1 / (eps0 * mu0)
    message("eps0  = ", eps0)
    message("mu0   = ", mu0)
    message("1 / (eps0 * mu0)  = ", c2_2b)

    

### Test run as a command line script ###
def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = "astropy.units / .constants playground",
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")

    args = arg.parse_args()

    if args.debug:
        ic.enable()
        ic(args, sys.version_info, sys.path)
    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()

    test()


if __name__ == "__main__":
    main()
