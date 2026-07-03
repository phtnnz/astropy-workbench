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
# Version 0.1 / 2026-06-25
#       Get observations from MPC database

VERSION     = "0.1 / 2026-06-25"
AUTHOR      = "Martin Junius"
NAME        = "mpc.observations"
DESCRIPTION = "Retrieve MPC observations data"

import re
from dataclasses import dataclass
from typing import Self

from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.time import Time
from astropy.table import Row
from astropy.table import Table, QTable
import astropy.units as u

from astroquery.mpc import MPC



@dataclass
class Obs:
    table: QTable = None

    def _rename_columns_mpc(self) -> None:
        self.table.rename_columns(("Date",    "Dec",      "V",             "Proper motion", "Direction", 
                                   "Azimuth", "Altitude", "Moon distance", "Moon altitude" ),
                                  # -->
                                  ("Obstime", "DEC",      "Mag",           "Motion",        "PA",        
                                   "Az",      "Alt",      "Moon_dist",     "Moon_alt"      ))


    def _id_type_from_name(name: str) -> str:
        id_type_regex = {   "asteroid number":        r'^[1-9][0-9]*$',
                            "asteroid designation":   r'^\d{4}[ _][A-Z]{1,2}\d{0,3}$',
                            "comet number":           r'^[0-9]{1,3}[PIA]$',
                            "comet designation":      r'^[PDCXAI]\/\d{4}[ _][A-Z]{1,2}\d{0,3}$'
                        }

        for id, regex in id_type_regex.items():
            m = re.match(regex, name)
            if m:
                ic(name, id)
                return id
        ## Default None or "asteroid designation"?
        return None


    def get_observations(self, obj: str) -> Self:
        table = MPC.get_observations(obj)
        # table["Targetname"] = obj
        # # table is already a QTable
        # self.table = QTable(table, meta={**table.meta})
        # self._rename_columns_mpc()
        self.table = table
        return self


    @classmethod
    def from_object(cls, obj: str) -> Self:
        obs = cls()
        obs.get_observations(obj)
        return obs
    

    def __getitem__(self, item) -> any:
        return self.table[item]
    

    def __len__(self) -> int:
        return len(self.table)


    def get_last_row_from_mpc(self) -> Row:
        # # Handle masked entries
        # for i in range(-1, -10, -1):
        #     mag = obs["mag"][i].unmasked
        #     if mag > Magnitude(0):
        #         break
        return self.table[-1]


    def get_last_obs(self) -> Time:
            jd = self.table[-1].get("epoch")
            time = Time(jd, format="jd")
            time.format = "iso"
            return time
