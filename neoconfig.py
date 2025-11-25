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
# Version 0.1 / 2025-11-25
#       New global config module for NEO modules

VERSION     = "0.1 / 2025-11-25"
AUTHOR      = "Martin Junius"
NAME        = "neoconfig"
DESCRIPTION = "Global NEO config module"

import sys
import argparse

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# Local modules
from verbose import verbose, warning, error, message
from jsonconfig import JSONConfig, config



CONFIGFILE = "neo-config.json"
config = JSONConfig(CONFIGFILE)
config.set_error_on_missing()



def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = DESCRIPTION,
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")

    args = arg.parse_args()

    if args.debug:
        ic.enable()
        ic(sys.version_info, sys.path, args)
    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()

    config.info()
    json_keys = config.get_keys()
    for k in json_keys:
        verbose(f"{k:22s} : {config.get(k)}")



if __name__ == "__main__":
    main()
