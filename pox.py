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
# Version 0.0 / 2026-11-09
#       Read Autoslew POX files

VERSION     = "0.0 / 2025-11-09"
AUTHOR      = "Martin Junius"
NAME        = "pox"
DESCRIPTION = "Read Autoslew POX files"

import sys
import argparse
import re
from dataclasses import dataclass

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time        import Time
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt

# Local modules
from verbose import verbose, warning, error, message
from astroutils import location_to_string, get_location     # ...

DEFAULT_LOCATION = "M49"



# Command line options
class Options:
    loc: EarthLocation = None       # -l --location



class POX:
    # From ASA Tool NINA plugin:
    # public int Number { get; set; }
    # public string DateObs { get; set; }
    # public string TimeObs { get; set; }
    # public double ExpTime { get; set; }
    # public double TelescopeRA { get; set; }
    # public double SolvedRA { get; set; }
    # public double TelescopeDec { get; set; }
    # public double SolvedDec { get; set; }
    # public int PierSide { get; set; }
    pass

class POXList:
    def __init__(self):
        self.obstime = list()
        self.ra_telescope = list()
        self.ra_solved = list()
        self.dec_telescope = list()
        self.dec_solved = list()
        self.coord_telescope = None
        self.coord_solved = None
        self.altaz_telescope = None
        self.altaz_solved = None


    @classmethod
    def from_pox_file(cls, filename: str):
        self = cls()

        verbose(f"processing {filename}")
        with open(filename, "r") as file:
            total_number = None
            number = 0
            file_iter = iter(file)

            total_number = int( next(file_iter).strip("\n\"\'") )
            ic(total_number)

            while ( line := next(file_iter, "").strip("\n\"\'") ):
                # Number line
                ic(line)
                m = re.match(r'Number (\d+)', line)
                if not m:
                    warning(f"format error, reading stopped: {line}")
                    return self
                number = int(m.group(1))
                ic(number)

                # Obstime
                obstime = Time( next(file_iter, "").strip("\n\"\'") )
                # Time, ignored
                next(file_iter, "")
                # Exposure time, ignored
                next(file_iter, "")
                # RA telescope
                ra_telescope = float( next(file_iter, "").strip("\n\"\'") )
                # RA solved
                ra_solved = float( next(file_iter, "").strip("\n\"\'") )
                # DEC telescope
                dec_telescope = float( next(file_iter, "").strip("\n\"\'") )
                # DEC solved
                dec_solved = float( next(file_iter, "").strip("\n\"\'") )
                # Side of pier, ignored
                next(file_iter, "")

                # Check separator
                line = next(file_iter, "").strip("\n\"\'")
                m = re.match(r'\*+', line)
                if not m:
                    warning(f"format error, reading stopped: {line}")
                    return self

                ic(obstime, ra_telescope, ra_solved, dec_telescope, dec_solved)
                self.obstime.append(obstime)
                self.ra_telescope.append(ra_telescope)
                self.ra_solved.append(ra_solved)
                self.dec_telescope.append(dec_telescope)
                self.dec_solved.append(dec_solved)
        # /with
        return self


    def to_skycoord(self) -> tuple[SkyCoord, SkyCoord]:
        self.coord_telescope = SkyCoord(self.ra_telescope, self.dec_telescope, unit=(u.hour, u.deg))
        self.coord_solved = SkyCoord(self.ra_solved, self.dec_solved, unit=(u.hour, u.deg))
        return self.coord_telescope, self.coord_solved


    def to_altaz(self, loc: EarthLocation) -> tuple[AltAz, AltAz]:
        if not self.coord_telescope or not self.coord_solved:
            self.to_skycoord()
        self.altaz_telescope = self.coord_telescope.transform_to( AltAz(obstime=self.obstime, location=loc) )
        self.altaz_solved = self.coord_solved.transform_to( AltAz(obstime=self.obstime, location=loc) )
        return self.altaz_telescope, self.altaz_solved



def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = DESCRIPTION,
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("-l", "--location", help=f"coordinates, named location or MPC station code, default {DEFAULT_LOCATION}")
    arg.add_argument("file", nargs="+", help="POX data file(s)")
    args = arg.parse_args()

    if args.debug:
        ic.enable()
        ic(sys.version_info, sys.path, args)
    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()

    # Observer location
    loc = get_location(args.location if args.location else DEFAULT_LOCATION)
    ic(loc, loc.to_geodetic())
    Options.loc = loc
    verbose(f"location {location_to_string(loc)}")

    # POX files
    for idx, filename in enumerate(args.file):
        pox = POXList.from_pox_file(filename)

        coord, coord_solved = pox.to_skycoord()
        ic(coord, coord_solved)
        altaz, altaz_solved = pox.to_altaz(loc)
        ic(altaz, altaz_solved)

        sep = altaz.separation(altaz_solved)
        ic(sep.to(u.arcmin))
        pa = altaz.position_angle(altaz_solved)
        ic(pa.to(u.degree))

        # ##Test
        # altaz = AltAz(az=[0, 90, 180, 270]*u.degree, alt=[30, 40, 50, 60]*u.degree)
        # altaz_solved = AltAz(az=[0+15, 90+15, 180, 270]*u.degree, alt=[30-10, 40, 50+5, 60+5]*u.degree)
        # ic(altaz, altaz_solved)
        # sep = altaz.separation(altaz_solved)
        # ic(sep.to(u.arcmin))
        # pa = altaz.position_angle(altaz_solved)
        # ic(pa.to(u.degree))
        # ##/Test

        fig = plt.figure()
        ax = fig.add_subplot(111, polar=True)
        ax.set_ylim(90, 0)
        ax.set_theta_zero_location("N")
        
        r1, theta1  = altaz.alt.degree, altaz.az.radian
        ic(r1, theta1)
        ax.plot(theta1, r1, "o", zorder=0)

        r2, theta2  = altaz_solved.alt.degree, altaz_solved.az.radian
        ic(r2, theta2)
        ax.plot(theta2, r2, "o", zorder=1)

        dtheta = pa.radian
        dr = sep.arcmin
        ic(dtheta, dr)
        # dr, dtheta -> x, y for Az 0 = N !
        # arrows point from the solved coordinates in the direction of the needed correction,
        # length is the correction in arcmin, not displayed to scale
        ax.quiver(theta2, r2, -dr*np.sin(theta2 - dtheta), +dr*np.cos(theta2 - dtheta), zorder=2)

        plotfile = f"tmp/plot_{idx+1:02d}.png"
        verbose(f"plot saved to {plotfile}")
        plt.savefig(plotfile, bbox_inches="tight")
        plt.close()



if __name__ == "__main__":
    main()
