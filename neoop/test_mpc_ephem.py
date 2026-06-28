#!/usr/bin/env python

# Copyright 2026 Martin Junius
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
# Version 0.1 / 2026-06-27
#       Test mpc.ephem

from icecream import ic

from mpc.ephem import Ephem, edata_add_ephem_mpc
from neo.local import get_local_circumstances
from neo.classes import EphemData

ic.enable()
obj = "C/2026 L1"

# Observer location and local circumstances
local = get_local_circumstances("M49")
ic(local)

eph = Ephem()
eph.get_ephemeris(obj, local)
ic(eph)

eph = Ephem.from_object(obj, local)
ic(eph)

motion = eph["Motion"]
ic(motion, motion[0])

edata = EphemData("-", obj, sort_time=None, ephem=None, times=None, exposure=None, mag=None, motion=None)
edata_add_ephem_mpc(edata, local)
ic(edata)
print(edata.ephem["Targetname", "Obstime", "RA", "DEC", "Mag", 
                  "Motion", "PA", "Az", "Alt", "Moon_dist", "Moon_alt"])
