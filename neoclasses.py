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
# Version 0.2 / 2026-02-04
#       Added NEOCPListData

VERSION     = "0.2 / 2026-02-04"
AUTHOR      = "Martin Junius"
NAME        = "neoclasses"
DESCRIPTION = "Dataclasses for ephemeris/planning"

from dataclasses import dataclass

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.coordinates import EarthLocation, Angle
from astropy.time        import Time
from sbpy.data import Ephem
from astroplan import Observer
from astropy.units import Quantity, Magnitude

# Local modules
from astroutils import location_to_string



# Dataclasses
@dataclass
class NEOCPListData:
    """Extra data from NEOCP / PCCP list"""
    type: str                   # "PCCP" or "NEOCP"
    score: int                  # NEOCP score
    mag: Magnitude              # Magnitude
    nobs: int                   # Number of observations
    arc: Quantity               # Orbit arc (days last - first observation)
    notseen: Quantity           # Not seen for x days

    def __str__(self):
        return f"{self.type} {self.score} {self.mag} #{self.nobs} {self.arc} {self.notseen}"



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



@dataclass
class EphemTimes:
    """Times for planner"""
    start: Time                 # ephemeris start time
    end: Time                   # ephemeris end time
    before: Time                # time before meridian
    after: Time                 # time after meridian
    alt_start: Time             # start time of optimal altitude
    alt_end: Time               # end time of optimal altitude
    plan_start: Time            # planned start time
    plan_end: Time              # planned end time
    max_alt: Time = None        # time of max altitude



@dataclass
class EphemData:
    """All object data incl. ephemeris"""
    obj: str                    # object name
    sort_time: Time             # time for sorting objects
    ephem: Ephem                # ephemeris of object
    times: EphemTimes           # ephemeris/planned times
    exposure: Exposure          # exposure data object
    mag: Magnitude              # magnitude of object
    motion: Quantity            # max. motion of object
    ra: Angle = None            # planned RA
    dec: Angle = None           # planned DEC
    neocp: NEOCPListData = None # data from NEOCP list


@dataclass
class LocalCircumstances:
    """Observer location and time data"""
    loc: EarthLocation          # location
    observer: Observer          # astroplan observer
    naut_dusk: Time             # nautical dusk
    naut_dawn: Time             # nautical dawn
    epochs: dict                # epochs parameter for Ephem.from_mpc()/from_jpb()
    code: str = None            # MPC station code

    def __str__(self) -> str:
        return f"location {location_to_string(self.loc)} code {self.code if self.code else "---"}\nnautical twilight {self.naut_dusk.iso} / {self.naut_dawn.iso} ({self.naut_dusk.scale.upper()})"


@dataclass
class JPLWObs:
    """Data from JPL WObs servive"""
    designation: str            # 'Designation'             object name
    full_name: str              # 'Full name'               full name
    rise_time: str              # 'Rise time'               HH:MM rise time > min-elev
    transit_time: str           # 'Transit time'            HH:MM transit time
    set_time: str               # 'Set time'                HH:MM set time < min elev
    max_time_obs: str           # 'Max. time observable'    HH:MM observable time    
    ra: Angle                   # 'R.A.'                    RA (deviates from ephemeris!)
    dec: Angle                  # 'Dec.'                    DEC (")
    vmag: Magnitude             # 'Vmag'                    V magnitude
    helio_range: Quantity       # 'Helio. range (au)'
    topo_range: Quantity        # 'Topo.range (au)'
    obj_obs_sun: Angle          # 'Object-Observer-Sun (deg)'
    obj_obs_moon: Angle         # 'Object-Observer-Moon (deg)'
    galatic_lat: Angle          # 'Galactic latitude (deg)'

    def __str__(self) -> str:
        return f"{self.designation:11s} {self.rise:time:6s} {self.transit_time:6s} {self.set_time:6s}  {self.vmag}"


@dataclass
class MPCDLx:
    """Data from MPC DLU / DLN lists"""
    designation: str
    type: str
    ra: Angle
    dec: Angle
    vmag: Magnitude
    elongation: Angle
    motion: Quantity
    marker: str
    last_obs: Time
    code: str
    last_mag: Magnitude
    filter: str
    uncertainty: int
    arc: Quantity
