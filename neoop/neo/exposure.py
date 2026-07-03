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
# Version 1.0 / 2026-06-20
#       Exposure functions from utils.py
#
# Usage:
#       from neo.exposure import ...

VERSION = "1.0 / 2026-06-20"
AUTHOR      = "Martin Junius"
NAME        = "neo.exposure"
DESCRIPTION = "NEO exposure functions"

from icecream import ic
# Disable debugging
ic.disable()

# AstroPy & friends
import astropy.units as u
from astropy.units import Quantity, Magnitude

# Local modules
from utils.verbose import verbose, warning
from neo.config import config
from neo.classes import Exposure, EphemData, LocalCircumstances



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



def exposure_calc(max_motion: Quantity, mag: Magnitude) -> Exposure:
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
    exp   = single_exp(max_motion)                      # Single exposure / s
    if exp == None:                                     # Object too fast
        return None

    min_n_exp = config.min_n_exp
    max_n_exp = config.max_n_exp
    base_mag  = config.base_mag
    base_exp  = config.base_exp

    rel_brightness = 10 ** (0.4 * (mag.value - base_mag))
    total_exp = base_exp * u.s * rel_brightness         # Total exposure
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
    total_time = ( total_exp 
                   + config.dead_time_slew_center * u.s 
                   + config.dead_time_af * u.s
                   + config.dead_time_guiding * u.s  
                   + config.safety_margin * u.s
                   + n_exp * config.dead_time_image * u.s )
    ic(n_exp, exp, total_exp, total_time, perc_of_required)

    return Exposure(n_exp, exp, total_exp, total_time, perc_of_required)



def edata_add_exposure(edata: EphemData, local: LocalCircumstances) -> EphemData:
    obj = edata.obj
    if not edata.ephem:
        return
    if edata.motion != None and edata.mag != None:
        edata.exposure = exposure_calc(edata.motion, edata.mag)
        ic(edata.obj, edata.exposure)
    if not edata.exposure:
        warning(f"exposure calculation for {obj} failed, too fast?")
        if edata.motion != None:
            warning(f"motion={edata.motion:.2f}, limit={motion_limit():.2f}")

    return edata
