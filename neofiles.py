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
# Version 0.0 / 2026-03-02
#       New module for NEO obs planner file handling
#
# Usage:
#       import neofiles
#       neofiles.now
#       neofiles.prefix
#       neofiles.set_prefix(prefix)
#       neofiles.path(filename)

VERSION     = "0.1 / 2026-03-02"
AUTHOR      = "Martin Junius"
NAME        = "neofiles"
DESCRIPTION = "NEO file handling"

import os

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# AstroPy & friends
import astropy.units as u
from astropy.coordinates import Angle
from astropy.units import Quantity, Magnitude
from astropy.time import Time
from astropy.table import QTable, Row
import numpy as np
from sbpy.data import Ephem

# Local modules
from verbose import verbose, warning, error, message
from neoconfig import config



now = Time.now()
prefix = now.strftime("%Y%m%d")
neo_data_dir = config.neo_obs_data_dir or config.downloads



def set_prefix(new_prefix: str) -> None:
    now = Time.now()
    prefix = new_prefix



def path(filename: str) -> str:
    return os.path.join(neo_data_dir, f"{prefix}-{filename}")
