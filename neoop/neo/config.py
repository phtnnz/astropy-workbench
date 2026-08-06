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
#       New global config module for NEO modules
# Version 1.0 / 2026-06-16
#       Moved and adapted to new directory structure under neoop/
#
# Usage:
#       from neoconfig import config

VERSION     = "1.0 / 2026-06-16"
AUTHOR      = "Martin Junius"
NAME        = "neoconfig"
DESCRIPTION = "Global NEO config module"

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# Local modules
from utils.jsonconfig import JSONConfig, config



CONFIGFILE = "neo-config.json"
config = JSONConfig(CONFIGFILE)
config.set_error_on_missing()
