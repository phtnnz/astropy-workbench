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
# Version 0.0 / 2025-06-29
#       Test for plate solving with astrometry.net

VERSION = "0.0 / 2025-06-29"
AUTHOR  = "Martin Junius"
NAME    = "solve"

CONFIGFILE = "astrometry.json"

import sys
import argparse

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()
# Local modules
from verbose import verbose, warning, error
from jsonconfig import JSONConfig

# AstroPy
from astropy.coordinates import SkyCoord, get_constellation
import astropy.units as u
from astropy.table import Table, Row

# Astroquery
from astroquery.astrometry_net import AstrometryNet
from astroquery.exceptions import TimeoutError


config = JSONConfig(CONFIGFILE)



# Command line options
class Options:
    pass



def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = "Plate solver with astrometry.net",
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("filename", nargs="+", help="filename")

    args = arg.parse_args()

    if args.debug:
        ic.enable()
        ic(sys.version_info, sys.path, args)
    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()

    config.info()
    apikey = config.get("astrometry.net", "API key")
    ic(apikey)

    filename = args.filename[0]
    verbose(f"processing {filename}")

    # From https://astroquery.readthedocs.io/en/latest/astrometry_net/astrometry_net.html
    ast = AstrometryNet()
    ast.api_key = apikey

    # ic(ast.show_allowed_settings())
    
    submission_id = None
    while True:
        try:
            if not submission_id:
                wcs_header = ast.solve_from_image(filename, submission_id=submission_id, 
                                                  crpix_center=True, publicly_visible="n",
                                                  allow_commercial_use="n", allow_modifications="n")
            else:
                wcs_header = ast.monitor_submission(submission_id, solve_timeout=120)
        except TimeoutError as e:
            submission_id = e.args[1]
        else:
            # got a result, so terminate
            print("")                   # Added newline after "Solving ..." progress indicator
            break

    ic(wcs_header)
    if wcs_header:
        # Sucessful
        center_ra  = wcs_header["CRVAL1"]
        center_dec = wcs_header["CRVAL2"]
        center_x   = wcs_header["CRPIX1"]
        center_y   = wcs_header["CRPIX2"]
        image_w    = wcs_header["IMAGEW"]
        image_h    = wcs_header["IMAGEH"]

        scale = None
        cpu_time = None
        for comment in wcs_header["COMMENT"]:
            if comment.startswith("scale: "):
                scale = comment[len("scale: "):]
            if comment.startswith("cpu time: "):
                cpu_time = comment[len("cpu time: "):]

        verbose(f"center RA={center_ra} DEC={center_dec} pos={center_x}({image_w})/{center_y}({image_h}) {scale} cpu={cpu_time}")
    else:
        # Code to execute when solve fails
        warning("plate solve failed")



if __name__ == "__main__":
    main()
