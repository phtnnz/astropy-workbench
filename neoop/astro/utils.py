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
# Version 0.1 / 2026-07-04
#       Astropy utility functions

VERSION     = "0.1 / 2026-07-04"
AUTHOR      = "Martin Junius"
NAME        = "astro.utils"
DESCRIPTION = "Astropy utility functions"

from icecream import ic
# Disable debugging
ic.disable()

# AstroPy & friends
import astropy.units as u
from astropy.coordinates import Angle
from astropy.time import Time



def is_east(az: Angle) -> bool:
    """
    Test for east azimut position

    Parameters
    ----------
    az : Angle
        Azimut angle

    Returns
    -------
    bool
        True if east, False if west
    """
    az180     = 180 * u.degree
    return True if az >= 0 and az < az180 else False



def is_west(az: Angle) -> bool:
    """
    Test for west azimut position

    Parameters
    ----------
    az : Angle
        Azimut angle

    Returns
    -------
    bool
        True if west, False if east
    """
    az180     = 180 * u.degree
    az360     = 360 * u.degree
    return True if az >= az180 and az < az360 else False



TIME_FORMAT    = "%Y-%m-%d %H:%M:%S"
TIME_FORMAT_TZ = "%Y-%m-%d %H:%M:%S+0000"

def fmt_time(time: Time|None, add_tz: bool=False) -> str:
    """Format time

    Parameters
    ----------
    time : Time | None
        Time value
    add_tz : bool, optional
        Add timezone "+0000", by default False

    Returns
    -------
    str
        Formatted time
    """
    return time.strftime(TIME_FORMAT_TZ if add_tz else TIME_FORMAT) if time != None else "-" * 19
