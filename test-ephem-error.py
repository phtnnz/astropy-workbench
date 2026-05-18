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
import astropy.units as u
from sbpy.data import Ephem
from icecream import ic
from astropy.coordinates import EarthLocation
from astropy.time import Time



def main():
    loc = EarthLocation(lat=-23.23639*u.deg, lon=16.36167*u.deg , height=1825*u.m) # M49, Hakos, Namibia
    obj = "2026 JH2"

    eph = Ephem.from_mpc(obj, location=loc)
    ic(eph)

    epochs = {"start":  Time('2026-05-19 00:00:00.000'),
              "step":   30 * u.min,
              "stop":   Time('2026-05-19 06:00:00.000')
             }
    eph = Ephem.from_mpc(obj, location=loc, epochs=epochs)
    ic(eph)



if __name__ == "__main__":
    main()
