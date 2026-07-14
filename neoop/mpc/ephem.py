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
#       New ephemeris handling using astroquery directly

VERSION     = "0.1 / 2026-06-27"
AUTHOR      = "Martin Junius"
NAME        = "mpc.ephem"
DESCRIPTION = "MPC ephemeris and observations"

from dataclasses import dataclass
from typing import Self
from itertools import pairwise

from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.units import Quantity, Magnitude
from astropy.table import Table, QTable, Row
from astropy.time import Time
from astropy.coordinates import Angle, AltAz
import astropy.units as u
import numpy as np
from astroquery.mpc import MPC

# Local modules
from astro.utils import is_east, is_west
from neo.local import LocalCircumstances



@dataclass
class Ephem:
    table: QTable = None

    def __getitem__(self, item) -> any:
        return self.table[item]

    def __setitem__(self, item, value) -> None:
        self.table[item] = value

    def __len__(self) -> int:
        return len(self.table)


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


    def get_flip_times(self, col_obstime: str="Obstime", col_az: str="Az") -> tuple[Time, Time]:
        prev_time = None
        prev_az   = None

        for row in self.table:
            time = row[col_obstime]
            az   = row[col_az]
            # ic(time, az)
            if not prev_az == None:
                if is_east(prev_az) and is_west(az):     # South flip
                    return (prev_time, time)
                if is_west(prev_az) and is_east(az):     # North flip
                    return (prev_time, time)
            prev_time = time
            prev_az   = az

        # No meridian passing found
        return None, None


    def get_opt_alt_times(self, alt: Angle, col_obstime: str="Obstime", col_alt: str="Alt") -> tuple[Time, Time]:
        time_alt0 = None
        time_alt1 = None

        for row in self.table:
            if time_alt0 == None and row[col_alt] >= alt:
                time_alt0 = row[col_obstime]
            if time_alt0 != None and row[col_alt] >= alt:
                time_alt1 = row[col_obstime]
            if time_alt1 != None and row[col_alt] < alt:
                break
        
        return time_alt0, time_alt1


    def get_max_alt_time(self, col_obstime: str="Obstime", col_alt: str="Alt") -> Time:
        max_alt = -90 * u.degree
        time_max = None
        for row in self.table:
            if row[col_alt] > max_alt:
                max_alt = row[col_alt]
                time_max = row[col_obstime]
        return time_max


    def get_row_for_time(self, t: Time, col_obstime: str="Obstime") -> Row:
        for r1, r2 in pairwise(self.table):
            if r1[col_obstime] <= t and t <= r2[col_obstime]:
                return r1
        
        # No matching interval found
        return None


    def to_altaz(self, name: str, local: LocalCircumstances, col_obstime: str="Obstime", col_alt: str="Alt", col_az: str="Az") -> AltAz:
        """
        Convert ephemeris cols "Alt"/"Az" to AltAz object

        Parameters
        ----------
        id : str
            NEOCP id = temporary designation
        eph : Ephem
            Ephemeris, including alt/az

        Returns
        -------
        AltAz
            AltAz coordinates object for altitude/sky plot
        """
        altaz = AltAz(alt=self[col_alt], az=self[col_az], obstime=self[col_obstime], location=local.loc)
        altaz.name = name
        ic(altaz)
        return altaz
