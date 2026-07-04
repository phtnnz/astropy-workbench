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

from dataclasses import dataclass
from typing import Self

from icecream import ic
# Disable debugging
ic.disable()

# AstroPy & friends
import astropy.units as u
from astropy.units import Quantity, Magnitude

# Local modules
from neo.config import config



@dataclass
class Exposure:
    """Exposure data"""
    number: int                 # number of exposure
    single: Quantity            # single exposure time
    total: Quantity             # total net exposure time
    total_time: Quantity        # total gross exposure time incl. overhead
    percentage: float           # percentage of required total exposure time

    def __str__(self):
        return f"{self.number} x {self.single:2.0f} = {self.total:3.1f} ({self.percentage:.0f}%) / total {self.total_time:3.1f}"


    @classmethod
    def single_exp(cls, motion: Quantity) -> Quantity:
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


    @classmethod
    def motion_limit(cls) -> Quantity:
        """
        Get motion limit

        Returns
        -------
        Quantity
            Motion limit derived from minimum exposure time
        """
        exp_times = config.exposure_times
        return config.pixel_tolerance * config.resolution * u.arcsec / (exp_times[0] * u.s).to(u.min)


    @classmethod
    def from_motion_mag(cls, max_motion: Quantity, mag: Magnitude) -> Self:
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
        exp1 = cls.single_exp(max_motion)                      # Single exposure / s
        if exp1 == None:                                     # Object too fast
            return None

        min_n_exp = config.min_n_exp
        max_n_exp = config.max_n_exp
        base_mag  = config.base_mag
        base_exp  = config.base_exp

        rel_brightness = 10 ** (0.4 * (mag.value - base_mag))
        total_exp = base_exp * u.s * rel_brightness         # Total exposure
        n_exp = int(total_exp / exp1) + 1                    # Number of exposures
        ic(base_mag, base_exp, mag.value, rel_brightness, total_exp, n_exp)
        perc_of_required = 100.                             # Percentage actual / total exposure
        if n_exp < min_n_exp:
            perc_of_required = min_n_exp / n_exp * 100
            n_exp = min_n_exp
        if n_exp > max_n_exp:
            perc_of_required = max_n_exp / n_exp * 100
            n_exp = max_n_exp
        total_exp = (n_exp * exp1).to(u.min)
        total_time = ( total_exp 
                    + config.dead_time_slew_center * u.s 
                    + config.dead_time_af * u.s
                    + config.dead_time_guiding * u.s  
                    + config.safety_margin * u.s
                    + n_exp * config.dead_time_image * u.s )
        ic(n_exp, exp1, total_exp, total_time, perc_of_required)

        return cls(n_exp, exp1, total_exp, total_time, perc_of_required)
