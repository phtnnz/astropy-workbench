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
# Version 0.3 / 2026-02-24
#       Added JPLWObsData, MPCDLxData
# Version 0.4 / 2026-05-16
#       Added EphemDataList
# Version 1.0 / 2026-06-16
#       Moved and adapted to new directory structure under neoop/

VERSION = "1.0 / 2026-06-16"
AUTHOR      = "Martin Junius"
NAME        = "neoclasses"
DESCRIPTION = "Dataclasses for ephemeris/planning"

from dataclasses import dataclass
from typing import Self
from collections.abc import Callable

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
from astro.astroutils import location_to_string
from utils.verbose import verbose



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
class JPLWObsData:
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
    # Extra
    type: str                   #  Object type: neo, pha, comet

    def __str__(self) -> str:
        return f"{self.type.upper()}  {self.designation:11s} {self.rise_time:6s} {self.transit_time:6s} {self.set_time:6s}  {self.vmag.value:4.1f}"



@dataclass
class MPCDLxData:
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



@dataclass
class EphemData:
    """All object data incl. ephemeris"""
    type: str                   # object type: NEO, PHA, COMET, NEOCP, PCCP
    obj: str                    # object name
    sort_time: Time = None      # time for sorting objects
    ephem: Ephem = None         # ephemeris of object
    times: EphemTimes = None    # ephemeris/planned times
    exposure: Exposure = None   # exposure data object
    mag: Magnitude = None       # magnitude of object
    motion: Quantity = None     # max. motion of object
    ra: Angle = None            # planned RA
    dec: Angle = None           # planned DEC
    neocp: NEOCPListData = None # data from NEOCP list
    wobs: JPLWObsData = None    # data from JPL SBWOBS service
    dlx: MPCDLxData = None      # data from MPC DLU/DLN lists
    force: bool = False         # force observation of this object (e.g. from --force option)

    def __str__(self) -> str:
        wobs = self.wobs
        dlx = self.dlx
        if dlx:
            return f"{wobs.type.upper():5s}  {wobs.designation:11s} {wobs.rise_time:6s} {wobs.transit_time:6s} {wobs.set_time:6s}  {float(wobs.vmag.value):4.1f}  {dlx.uncertainty}  {str(dlx.last_obs):10.10s}"
        else:
            return str(wobs)



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



class EphemDataList(list):
    @classmethod
    def from_dict(cls, obj_edata: dict[str, EphemData]) -> Self:
        return cls(obj_edata.values())
    
    @classmethod
    def from_objects(cls, objects: list[str]) -> Self:
        return cls([ EphemData("-", obj) for obj in objects ])

    def objects(self) -> list[str]:
        return [ edata.obj for edata in self ]

    def objects_str(self) -> str:
        return ", ".join(self.objects())
    
    def len(self) -> int:
        return len(self)
    
    def append_objects(self, objects: list[str]) -> Self:
        self.extend([ EphemData("-", obj) for obj in objects ])
        return self

    def sort_by_time(self) -> None:
        return self.sort(key=lambda item: item.sort_time)

    def __str__(self) -> str:
        return "\n".join(f"{edata.type} {edata.obj} {edata.sort_time.iso if edata.sort_time else ''}" for edata in self)

    def verbose_ephem(self) -> None:
        for edata in self:
            if edata.ephem:
                verbose("===================================================================================================================")
                verbose(f"{edata.obj} ephemeris")
                verbose.print_lines2(edata.ephem)
        verbose("===================================================================================================================")

    def process(self, func: Callable, local: LocalCircumstances) -> None:
        for edata in self:
            func(edata, local)



@dataclass
class PrevNEOCPData:
    trk_sub: str
    designation: str = None
    comment: str = None
    date: Time = None
    mpec_no: str = None
    mpec_url: str = None

    def __str__(self) -> str:
        if self.designation:
            if self.mpec_no:
                return f"{self.trk_sub:7s} = {self.designation:10s}  MPEC {self.mpec_no:9s}  {self.mpec_url}"
            else:
                return f"{self.trk_sub:7s} = {self.designation:10s}"
        else:
            return f"{self.trk_sub:7s}   {self.comment}"
