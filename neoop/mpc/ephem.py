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
# Version 0.1 / 2026-06-27
#       New ephemeris handling using astroquery directly,
#       Ephem class defined in neo.classes

VERSION     = "0.1 / 2026-06-27"
AUTHOR      = "Martin Junius"
NAME        = "mpc.ephem"
DESCRIPTION = "MPC ephemeris and observations"

from dataclasses import dataclass
from typing import Self

from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
import astropy.units as u

# Local modules
from neo.local import LocalCircumstances

# AstroPy
from astropy.units import Quantity, Magnitude
from astropy.table import Table, QTable
import astropy.units as u
import numpy as np
from astroquery.mpc import MPC



@dataclass
class Ephem:
    table: QTable = None

    def _rename_columns(self) -> None:
        self.table.rename_columns(("Date",    "Dec",      "V",             "Proper motion", "Direction", 
                                   "Azimuth", "Altitude", "Moon distance", "Moon altitude" ),
                                  # -->
                                  ("Obstime", "DEC",      "Mag",           "Motion",        "PA",        
                                   "Az",      "Alt",      "Moon_dist",     "Moon_alt"      ))


    def _convert_columns(self) -> None:
        self.table["RA"] = self.table["RA"].to(u.hourangle)


    def get_ephemeris(self, obj: str, local: LocalCircumstances) -> Self:
        # For compatibility with sbpy.data.Ephem
        start = local.epochs.get("start")
        step  = local.epochs.get("step")
        stop  = local.epochs.get("stop")
        if stop:
            number = int((stop - start) / step) + 1
        else:
            number = 10
        ic(start, step, stop, number)
        
        table = MPC.get_ephemeris(obj, location=local.loc,
                                  start=start, step=step, number=number)
        table["Targetname"] = obj
        # Adopted from sbpy.data: convert Table returned from MPC.get_ephemeris
        # to QTable with Quantity when indexing
        self.table = QTable(table, meta={**table.meta})
        self._rename_columns()
        self._convert_columns()
        return self


    @classmethod
    def from_object(cls, obj: str, local: LocalCircumstances) -> Self:
        ephem = cls()
        ephem.get_ephemeris(obj, local)
        return ephem
    

    @classmethod
    def from_table(cls, table: Table) -> Self:
        ephem = cls()
        ephem.table = QTable(table, meta={**table.meta})
        return ephem


    def get_mag0(self, column: str="Mag") -> Magnitude:
        return self.table[column][0]


    def get_max_motion(self, column: str="Motion") -> Quantity:
        max_m = -1 * u.arcsec / u.min
        if "," in column:
            # Separate RA*cos(DEC), DEC motion columns
            col1, col2 = column.split(",")
        else:
            # Single column with proper motion
            col1 = column
            col2 = None
        for row in self.table:
            if col2:
                motion = np.sqrt( np.square(row[col1]) + np.square(row[col2]) )
            else:
                motion = row[col1]
            if motion > max_m:
                max_m = motion

        return max_m.to(u.arcsec / u.min)


    def __getitem__(self, item) -> any:
        return self.table[item]
    

    def __setitem__(self, item, value) -> None:
        self.table[item] = value
    

    def __len__(self) -> int:
        return len(self.table)
