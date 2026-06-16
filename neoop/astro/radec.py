#!/usr/bin/env python

# Copyright 2023-2025 Martin Junius
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
# Version 0.1 / 2024-07-29
#       Module for parsing ra/dec coordinates in various format, not depending on astropy
# Version 0.2 / 2025-10-03
#       Fixed sign handling for DEC
# Version 0.3 / 2025-10-15
#       Rewritten using AstroPy's SkyCoord

VERSION = "0.3 / 2025-12-15"
AUTHOR  = "Martin Junius"
NAME    = "radec"

import argparse

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
import astropy.units as u

# Local modules
from verbose import verbose, error



# Coord class derived from SkyCoord
class Coord(SkyCoord):
    """Coord(inate) class based on AstroPy SkyCoord"""

    def __init__(self, ra: any, dec: any) -> None:
        """
        Initialize Coord object

        Parameters
        ----------
        ra : any
            RA coordinate, any input accepted by SkyCoord (hourangle)
        dec : any
            DEC coordinate, any input accepted by SkyCoord (degree)
        """
        ic(ra, dec)
        super().__init__(ra, dec, unit=(u.hourangle, u.degree))
        ic(self.ra, self.dec)

        self.ra_h,  self.ra_m,  self.ra_s  = self._split_hms(self.ra)
        self.dec_sign, self.dec_d, self.dec_m, self.dec_s = self._split_signed_dms(self.dec)
        _, self.dec_pm, self.dec_neg = self._dec_sign(self.dec_sign)

        ic(self.ra_h,  self.ra_m,  self.ra_s,
            self.dec_sign, self.dec_d, self.dec_m, self.dec_s,
            self.dec_pm, self.dec_neg)


    def _split_signed_dms(self, a: Angle) -> tuple[int, int, int, float]:
        """
        Internal: get signed_dms tuple for DEC

        Parameters
        ----------
        a : Angle
            DEC angle

        Returns
        -------
        tuple[int, int, int, float]
            sign, degrees (positive), minutes, seconds
        """
        sgn, d, m, s = a.signed_dms
        return int(sgn), int(d), int(m), float(s)


    def _split_hms(self, a: Angle) -> tuple[int, int, float]:
        """
        Internal: get hms tuple for RA

        Parameters
        ----------
        a : Angle
            RA hourangle

        Returns
        -------
        tuple[int, int, float]
            hours, minutes, seconds
        """
        h, m, s = a.hms
        return int(h), int(m), float(s)


    def _dec_sign(self, v: float) -> tuple[int, str, bool]:
        """
        Internal: get sign representations for DEC

        Parameters
        ----------
        v : float
            DEC decimal angle or -1/+1

        Returns
        -------
        tuple[int, str, bool]
            DEC sign tuple (-1, "-", True) if negative / (1, "+", False) if positive
        """
        return (-1, "-", True) if v < 0 else (+1, "+", False)


    def to_string(self, format: str="hmsdms") -> str:
        """
        Convert coordinate to string

        Parameters
        ----------
        format : str, optional
            Format string, by default "hmsdms"
            Supported: "hmsdms", "decimal", " ", "mpc", "mpc1"

        Returns
        -------
        str
            String representation of coordinates
        """
        if format == "decimal":
            return f"{self.ra.degree:.7f} {self.dec.degree:.7f}"
        if format == " ":
            return f"{self.ra_h:02d} {self.ra_m:02d} {self.ra_s:06.3f} {self.dec_pm}{self.dec_d:02d} {self.dec_m:02d} {self.dec_s:06.3f}"
        if format == "mpc":
            return f"{self.ra_h:02d} {self.ra_m:02d} {self.ra_s:06.3f}{self.dec_pm}{self.dec_d:02d} {self.dec_m:02d} {self.dec_s:05.2f}"
        if format == "mpc1":
            return f"{self.ra_h:02d} {self.ra_m:02d} {self.ra_s:05.2f} {self.dec_pm}{self.dec_d:02d} {self.dec_m:02d} {self.dec_s:04.1f} "
        return f"{self.ra_h:02d}h{self.ra_m:02d}m{self.ra_s:06.3f}s {self.dec_pm}{self.dec_d:02d}d{self.dec_m:02d}m{self.dec_s:06.3f}s"


    def __repr__(self) -> str:
        """
        String representation

        Returns
        -------
        str
            String representation, format "hmsdms"
        """
        return self.to_string()



### Test run as a command line script ###
def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = "RA/DEC parser",
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("ra", help="ra")
    arg.add_argument("dec", help="dec")

    args = arg.parse_args()

    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()
    if args.debug:
        ic.enable()

    print(f"{args.ra=} {args.dec=}")
    try:
        coord1 = Coord(args.ra, args.dec)
    except ValueError as e:
        error(f"illegal RA/DEC values")

    print(f"{coord1 = }\n{coord1.to_string(format='decimal') = }")
    print("Regression with coord1 decimal values ...")
    coord2 = Coord(coord1.ra, coord1.dec)
    print(f"{coord2 = }")
    print(f"{coord2.to_string(format='hmsdms')  = }")
    print(f"{coord2.to_string(format='decimal') = }")
    print(f"{coord2.to_string(format=' ')       = }")
    print(f"{coord2.to_string(format='mpc')     = }")
    print(f"{coord2.to_string(format='mpc1')    = }")


    

if __name__ == "__main__":
    main()
