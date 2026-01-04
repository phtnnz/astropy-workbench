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
# Version 0.1 / 2026-01-04
#       Dataclasses from neoephem, neoutils

VERSION     = "0.1 / 2026-01-04"
AUTHOR      = "Martin Junius"
NAME        = "neoclasses"
DESCRIPTION = "Dataclasses for ephemeris/planning"

from dataclasses import dataclass

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.coordinates import EarthLocation
from astropy.time        import Time
from sbpy.data import Ephem
from astroplan import Observer
from astropy.units import Quantity

# Local objects
from astroutils import location_to_string



# Dataclasses
@dataclass
class Exposure:
    number: int             # number of exposure
    single: Quantity        # single exposure time
    total: Quantity         # total net exposure time
    total_time: Quantity    # total gross exposure time incl. overhead
    percentage: float       # percentage of required total exposure time



@dataclass
class EphemTimes:
    start: Time             # ephemeris start time
    end: Time               # ephemeris end time
    before: Time            # time before meridian
    after: Time             # time after meridian
    alt_start: Time         # start time of optimal altitude
    alt_end: Time           # end time of optimal altitude
    plan_start: Time        # planned start time
    plan_end: Time          # planned end time



@dataclass
class EphemData:
    obj: str                # object name
    sort_time: Time         # time for sorting objects
    ephem: Ephem            # ephemeris of object
    times: EphemTimes       # ephemeris/planned times
    exposure: Exposure      # exposure data object



@dataclass
class LocalCircumstances:
    loc: EarthLocation      # location
    observer: Observer      # astroplan observer
    naut_dusk: Time         # nautical dusk
    naut_dawn: Time         # nautical dawn
    epochs: dict            # epochs parameter for Ephem.from_mpc()/from_jpb()

    def __str__(self):
        return f"location {location_to_string(self.loc)}\nnautical twilight {self.naut_dusk.iso} / {self.naut_dawn.iso} ({self.naut_dusk.scale.upper()})"
