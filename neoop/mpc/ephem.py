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

from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
import astropy.units as u

# Local modules
from utils.verbose import verbose, warning
from neo.classes import Ephem, EphemData, LocalCircumstances
from neo.config import config


def edata_add_ephem_mpc(edata: EphemData, local: LocalCircumstances) -> None:
    if edata.ephem:
        verbose(f"already got ephemeris for {edata.obj}")
        return

    min_alt = config.min_alt

    obj = edata.obj
    verbose(f"{obj} ephemeris from MPC")
    eph = Ephem.from_object(obj, local)
    mask = (eph["Alt"] > min_alt * u.deg) & (eph["Obstime"] >= local.naut_dusk) & (eph["Obstime"] <= local.naut_dawn)
    eph1 = Ephem(eph[mask])
    if len(eph1) == 0:
        warning(f"skipping empty ephemeris for {obj}")
        return
    mag = eph1.get_mag0()
    motion = eph1.get_max_motion()

    # Copy to EphemData
    if edata.wobs:
        edata.type = edata.wobs.type.upper()
    edata.obj = obj
    edata.ephem = eph1
    edata.mag = mag
    edata.motion = motion
