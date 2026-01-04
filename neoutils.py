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
# Version 0.1 / 2025-11-25
#       New module with NEO utility functions
#
# Usage:
#       from neoutils import ...

VERSION     = "0.1 / 2025-11-25"
AUTHOR      = "Martin Junius"
NAME        = "neoutils"
DESCRIPTION = "NEO utility functions"

import re
from typing import Tuple

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
import astropy.units as u
from astropy.units import Quantity, Magnitude
from astropy.time import Time, TimeDelta
from astropy.table import QTable, Row
import numpy as np

# Local modules
from verbose import verbose, warning, error, message
from neoconfig import config
from neoclasses import Exposure



def single_exp(motion: Quantity) -> Quantity:
    """
    Compute exposure time depending on motion value,

    Parameters
    ----------
    motion : Quantity
        Motion value as arcsec/min

    Returns
    -------
    Quantity
        Exposure time as secs or None, if too fast
    """
    exp_times = config.exposure_times
    exp = config.pixel_tolerance * config.resolution * u.arcsec / motion
    #            ^ pixels                 ^ arcsec/pixel          ^ arcsec/min
    exp = exp.to(u.s).value
    exp_min = exp_times[0]
    if exp < 0.9 * exp_min: # allow a bit of tolerance
        return None
    for exp1 in exp_times:
        if exp1 > exp:
            break
        exp_min = exp1
    return exp_min * u.s



def motion_limit() -> Quantity:
    """
    Get motion limit

    Returns
    -------
    Quantity
        Motion limit derived from minimum exposure time
    """
    exp_times = config.exposure_times
    return config.pixel_tolerance * config.resolution * u.arcsec / (exp_times[0] * u.s).to(u.min)



def total_exp(max_motion: Quantity, mag: Magnitude) -> Exposure:
    """
    Calculate number of exposures, single exposure time, total exposure time, total time, percentage
    from object motion and magnitude

    Parameters
    ----------
    max_motion : Quantity
        Object motion
    mag : Magnitude
        Object magnitude

    Returns
    -------
    Exposure
        Exposure object
    """
    exp   = single_exp(max_motion)     # Single exposure / s
    if exp == None:                                     # Object too fast
        return None

    min_n_exp = config.min_n_exp
    max_n_exp = config.max_n_exp
    base_mag  = config.base_mag
    base_exp  = config.base_exp

    rel_brightness = 10 ** (0.4 * (mag.value - base_mag))
    total_exp = base_exp * u.s * rel_brightness  # Total exposure
    n_exp = int(total_exp / exp) + 1                    # Number of exposures
    ic(base_mag, base_exp, mag.value, rel_brightness, total_exp, n_exp)
    perc_of_required = 100.                             # Percentage actual / total exposure
    if n_exp < min_n_exp:
        perc_of_required = min_n_exp / n_exp * 100
        n_exp = min_n_exp
    if n_exp > max_n_exp:
        perc_of_required = max_n_exp / n_exp * 100
        n_exp = max_n_exp
    total_exp = (n_exp * exp).to(u.min)
    total_time = (  total_exp 
                    + config.dead_time_slew_center * u.s 
                    + config.dead_time_af * u.s
                    + config.dead_time_guiding * u.s  
                    + config.safety_margin * u.s
                    + n_exp * config.dead_time_image * u.s )
    ic(n_exp, exp, total_exp, total_time, perc_of_required)

    return Exposure(n_exp, exp, total_exp, total_time, perc_of_required)



def max_motion(ephemeris: QTable, column: str="motion") -> Quantity:
    """
    Get max value for motion column(s) from ephemeris table

    Parameters
    ----------
    ephemeris : QTable
        Ephemeris table
    column : str, optional
        Motion column name(s), by default "motion", comma separated for RA/DEC motion

    Returns
    -------
    Quantity
        Max motion of object
    """
    max_m = -1 * u.arcsec / u.min
    if "," in column:
        # Separate RA*cos(DEC), DEC motion columns
        col1, col2 = column.split(",")
    else:
        # Single column with proper motion
        col1 = column
        col2 = None
    for row in ephemeris:
        if col2:
            motion = np.sqrt( np.square(row[col1]) + np.square(row[col2]) )
        else:
            motion = row[col1]
        if motion > max_m:
            max_m = motion

    return max_m



def exposure_from_ephemeris(ephemeris: QTable, column: str, mag: Magnitude) -> Exposure:
    max_m = max_motion(ephemeris, column)
    if max_m == None:
        return None
    
    ic(max_m, mag)
    return total_exp(max_m, mag)



def id_type_from_name(name: str) -> str:
    id_type_regex = {   "asteroid numer":         "^[1-9][0-9]*$",
                        "asteroid designation":   "^20[0-9]{2}[ _][A-Z]{2}[0-9]{0,3}$",
                        "comet number":           "^[0-9]{1,3}[PIA]$",
                        "comet designation":      "^[PDCXAI]/20[0-9]{2}[ _][A-Z]{2}[0-9]{0,3}$"
                    }

    for id in id_type_regex.keys():
        m = re.match(id_type_regex[id], name)
        if m:
            ic(name, id)
            return id
    ## Default None or "asteroid designation"?
    return None


