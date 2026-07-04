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
#       LocalCircumstandes functions

VERSION     = "0.1 / 2026-06-27"
AUTHOR      = "Martin Junius"
NAME        = "neo.local"
DESCRIPTION = "Location and local circumstances"

import re
from dataclasses import dataclass
from typing import Self

from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.time import Time
from astropy.coordinates import EarthLocation, Angle

import astropy.units as u
import numpy as np
from astroplan import Observer

# Local modules
from mpc.location import get_location



@dataclass
class LocalCircumstances:
    """Observer location and time data"""
    loc: EarthLocation          # location
    observer: Observer          # astroplan observer
    naut_dusk: Time             # nautical dusk
    naut_dawn: Time             # nautical dawn
    epochs: dict                # epochs parameter for Ephem.from_mpc()/from_jpb()
    code: str = None            # MPC station code
    midnight: Time = None       # midnight

    def _loc_to_string(self) -> str:
        return f"lon={self.loc.to_geodetic().lon.to_string(unit=u.degree, precision=1)} lat={self.loc.to_geodetic().lat.to_string(unit=u.degree, precision=1)} height={self.loc.to_geodetic().height.to_string(precision=0)}"


    def __str__(self) -> str:
        return f"location {self._loc_to_string()} code {self.code if self.code else "---"}\nnautical twilight {self.naut_dusk.iso} / {self.naut_dawn.iso} ({self.naut_dusk.scale.upper()})"


    @classmethod
    def from_location(cls, location: str) -> Self:
        loc = get_location(location)
        ic(loc, loc.to_geodetic())
        # MPC station code
        m = re.search(r'^([0-9A-Z]{3})$', location)
        if m:
            code = m.group(1)
        else:
            code = None
        ic(code)
    
        # Astroplan
        observer = Observer(location=loc, description=loc.info.name)
        ic(observer)

        # Observation times for upcoming night
        time = Time.now()
        ic(time)

        midnight = observer.midnight(time, which="next")
        twilight_evening = observer.twilight_evening_nautical(time, which="next")
        twilight_morning = observer.twilight_morning_nautical(time, which="next")
        if twilight_evening > twilight_morning:
            twilight_evening = observer.twilight_evening_nautical(time, which="previous")
        ic(midnight.iso, twilight_evening.iso, twilight_morning.iso)

        # Round midnight time to nearest 30 min
        rem, day = np.modf(midnight.jd)
        n_round = 24 * 2    # 24 h / 30 min
        rem = round(rem*n_round) / n_round
        jd1 = day + rem
        midnight1 = Time(jd1, format="jd")
        ic(day, rem, midnight1.iso)
        epochs = {"start":  midnight1 - 8 * u.hour,
                "step":   30 * u.min,
                "stop":   midnight1 + 9 * u.hour
                }
        ic(epochs)

        return cls(loc, observer, twilight_evening, twilight_morning, epochs, code, midnight1)


    def get_dec_limits(self, min_alt: Angle) -> tuple[Angle, Angle]:
        lat = self.loc.lat
        min_dec = -90*u.deg + min_alt + lat
        max_dec = +90*u.deg - min_alt + lat
        if min_dec < -90*u.deg:
            min_dec = -90*u.deg 
        if max_dec > +90*u.deg:
            max_dec = +90*u.deg 
        ic(lat, min_alt, min_dec, max_dec)
        return Angle(min_dec), Angle(max_dec)
