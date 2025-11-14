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
# Version 0.0 / 2025-11-14
#       Retrieve "What's Observable?" from JPL

VERSION     = "0.0 / 2025-11-14"
AUTHOR      = "Martin Junius"
NAME        = "jplsbwobs"
DESCRIPTION = "Retrieve JPL What's Observable"

import sys
import argparse
import re
import requests
from typing import Tuple, Any

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.coordinates import Angle
from astropy.coordinates import errors
from astropy.time        import Time, TimeDelta
import astropy.units as u
import astropy.constants as const
import numpy as np

# Local modules
from verbose import verbose, warning, error, message
from astroutils import location_to_string, get_location     # ...
from jsonconfig import JSONConfig, config

DEFAULT_LOCATION = "M49"

CONFIGFILE = "jplsbwobs.json"
config = JSONConfig(CONFIGFILE)
config.set_error_on_missing()



# Command line options
class Options:
    pass



def jpl_query_sbwobs(url: str, filename: str) -> None:
    """
    Retrieve What's Observable from JPL

    Parameters
    ----------
    url : str
        JPL API URL
    filename : str
        Local file name in cache
    """
    ic(url, filename)

    timeout = config.requests_timeout
    ##FIXME: get from command line, twilight times from astroplan start/end -> obs-time/obs-end
    time_obs = Time.now().strftime("%Y-%m-%d")

    data = { 
        "mpc-code":     config.mpc_code,    # MPC observatory code
        "obs-time":     time_obs,           # Date/time of the observation
        "elev-min":     config.elev_min,    # Minimum altitude
        "vmag-max":     config.vmag_max,    # Max V mag = minimum brightness
        "output-sort":  config.output_sort, # Sort records

        "sb-ns":        config.sb_ns,       # Numbered (n) ./. unnumbered (u)
        "sb-kind":      config.sb_kind,     # Asteroids (a) ./. comets (c)
        "sb-group":     config.sb_group     # NEOs (neo) ./. PHA (pha)
    }

    ic(url, data)
    response = requests.get(url, params=data, timeout=timeout)
    ic(response.status_code)
    if response.status_code != 200:
        error(f"query to {url} failed")

    with open(filename, mode="w", encoding=response.encoding) as file:
        file.write(response.text)



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

    ##TEST##
    jpl_query_sbwobs(config.sbwobs_url, "tmp/sbwobs.json")


if __name__ == "__main__":
    main()
