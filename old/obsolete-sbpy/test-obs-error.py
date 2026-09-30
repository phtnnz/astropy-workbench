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
#
# https://github.com/NASA-Planetary-Science/sbpy/issues/451

VERSION     = "0.0 / 2025-xx-xx"
AUTHOR      = "Martin Junius"
NAME        = "template"
DESCRIPTION = "mj's Python template"

import sys
import argparse

# AstroPy
import astropy.units as u
from sbpy.data import Ephem, Obs
from icecream import ic
from astropy.coordinates import EarthLocation
from astropy.time import Time



def main():
    loc = EarthLocation(lat=-23.236549*u.deg, lon=16.361720*u.deg , height=1853*u.m) # M49, Hakos, Namibia

    # Ok
    obs = Obs.from_mpc("2026 MP1", id_type='asteroid designation')
    ic(obs)

    # Ok
    obs = Obs.from_mpc("C/2026 L1", id_type='comet designation')
    ic(obs)

    # Fails 
    # "20001 results. Please restrict query to return less than 20001 results."
    obs = Obs.from_mpc("2026 MP1")
    ic(obs)



if __name__ == "__main__":
    main()
