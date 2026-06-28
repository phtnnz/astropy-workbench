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

VERSION     = "0. / 2026-06-25"
AUTHOR      = "Martin Junius"
NAME        = "mpc.observations"
DESCRIPTION = "Retrieve MPC observations data"

import re

from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.time import Time
from astropy.table import Row
from astropy.units import Quantity
import astropy.units as u

# Local modules
from utils.verbose import verbose, warning, error, message
from neo.classes import Obs
from astro.astroutils import time_jd_as_iso



def get_obs_from_mpc(obj: str) -> Obs:
    obs = Obs.from_object(obj)
    return obs



def get_last_row_from_mpc(obj: str) -> Row:
    obs = get_obs_from_mpc(obj)
    # # Handle masked entries
    # for i in range(-1, -10, -1):
    #     mag = obs["mag"][i].unmasked
    #     if mag > Magnitude(0):
    #         break
    return obs.table[-1]



def get_last_obs_from_mpc(obj: str) -> Quantity:
    lastrow = get_last_row_from_mpc(obj)
    return time_jd_as_iso(lastrow.get("epoch"))
