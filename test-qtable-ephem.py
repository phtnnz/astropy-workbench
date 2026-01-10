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
# Version 0.0 / 2025-xx-xx
#       TEXT

VERSION     = "0.0 / 2025-xx-xx"
AUTHOR      = "Martin Junius"
NAME        = "template"
DESCRIPTION = "mj's Python template"

import sys
import argparse

# AstroPy
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.coordinates import Angle
from astropy.coordinates import errors
from astropy.time  import Time, TimeDelta
from astropy.table import QTable, Row
import astropy.units as u
import astropy.constants as const
import numpy as np
from sbpy.data import Ephem
from icecream import ic



def main():

    qt = QTable()
    qt["Obstime"]   = Time("2000-01-01 00:00")
    qt["RA"]        = 0 * u.hourangle
    qt["DEC"]       = 0 * u.degree
    qt.add_row([ Time("2026-01-01"), 10*u.hourangle,  12.34*u.deg ])
    qt.add_row([ Time("2026-01-02"), 12*u.hourangle, -56.78*u.deg ])

    ic(qt)
    ic(qt["RA"][0])
    for row in qt:
        ic(row)
        ic(row["RA"])

    eph = Ephem.from_table(qt)
    ic(eph)
    ic(eph["RA"][0])
    for row in eph:
        ic(row)
        ic(row["RA"])



if __name__ == "__main__":
    main()
