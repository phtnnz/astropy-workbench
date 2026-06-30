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
# Version 1.1 / 2026-06-28
#       Added Ephem class (replacing sbpy)
# Version 1.2 / 2026-06-28
#       Added Obs class (replacing sbpy)
# Version 2.0 / 2026-06-30
#       Now just a wrapper for all the various classes

VERSION     = "2.0 / 2026-06-30"
AUTHOR      = "Martin Junius"
NAME        = "neo.classes"
DESCRIPTION = "Dataclasses for ephemeris/planning"



# One-stop shopping ;-) for all classes
from neo.local           import LocalCircumstances
from mpc.observations    import Obs
from mpc.ephem           import Ephem
from mpc.ephemdata       import EphemData, EphemDataList, EphemTimes, Exposure, NEOCPListData, JPLWObsData, MPCDLxData
from mpc.prevneocp       import PrevNEOCPData
