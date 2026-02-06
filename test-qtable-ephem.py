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
from astropy.time  import Time
from astropy.table import Table, QTable, Row
import astropy.units as u
from sbpy.data import Ephem
from icecream import ic



def get0(tlike: QTable|Row, col: str) -> any:
    if isinstance(tlike, Table) or isinstance(tlike, Ephem):
        return tlike[col][0]
    if isinstance(tlike, Row):
        return tlike[col]
    raise ValueError(f"{type(tlike)}: not a Table or Row")



def main():
    qt = QTable()
    qt["Obstime"]   = Time("2000-01-01 00:00")
    qt["RA"]        = 0 * u.hourangle
    qt["DEC"]       = 0 * u.degree
    qt.add_row([ Time("2026-01-01"), 10*u.hourangle,  12.34*u.deg ])
    qt.add_row([ Time("2026-01-02"), 12*u.hourangle, -56.78*u.deg ])

    ic("QTable")
    ic(qt, len(qt), qt["RA"][0])
    for row in qt:
        ic(row, row["RA"])

    eph = Ephem.from_table(qt)
    ic("Ephem")
    ic(eph, len(eph), eph["RA"][0])
    for row in eph:
        ic(row, row["RA"], get0(row, "RA"))
        ##HACK: bind get0() method to row object, NOT RECOMMENDED
        row.get0 = get0.__get__(row)
        ic(row.get0("RA"))
        ##HACK: get a proper Row object
        row0 = row._table[0]
        ic(row0)



if __name__ == "__main__":
    main()
