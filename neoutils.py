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
# Version 0.1 / 2025-11-25
#       New module with NEO utility functions
#
# Usage:
#       from neoutils import ...

VERSION     = "0.1 / 2025-11-25"
AUTHOR      = "Martin Junius"
NAME        = "neoutils"
DESCRIPTION = "NEO utility functions"

from typing import Tuple
from dataclasses import dataclass

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
import astropy.units as u
from astropy.units import Quantity, Magnitude
from astropy.time import Time, TimeDelta
from astropy.table import QTable, Row

# Local modules
from verbose import verbose, warning, error, message
from neoconfig import config



@dataclass
class Exposure:
    number: int
    single: Quantity
    total: Quantity
    total_time: Quantity
    percentage: float



def single_exp_time_from_motion(motion: Quantity) -> Quantity:
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
    exp   = single_exp_time_from_motion(max_motion)     # Single exposure / s
    if exp == None:                                     # Object too fast
        return None

    min_n_exp = config.min_n_exp
    max_n_exp = config.max_n_exp

    rel_brightness = 10 ** (0.4 * (mag.value - config.base_mag))
    total_exp = config.base_exp * u.s * rel_brightness  # Total exposure
    n_exp = int(total_exp / exp) + 1                    # Number of exposures
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
    Get max value for motion column from ephemeris table

    Parameters
    ----------
    ephemeris : QTable
        Ephemeris table
    column : str, optional
        Motion column name, by default "motion"

    Returns
    -------
    Quantity
        Max motion of object
    """
    max_motion = -1 * u.arcsec / u.min
    for row in ephemeris:
        if row[column] > max_motion:
            max_motion = row[column]
    return max_motion

