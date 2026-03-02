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
# Version 0.0 / 2026-03-02
#       New module for NEO obs planner file handling
#
# Usage:
#       import neofiles
#       neofiles.now
#       neofiles.prefix
#       neofiles.set_prefix(prefix)
#       neofiles.path(filename)

VERSION     = "0.1 / 2026-03-02"
AUTHOR      = "Martin Junius"
NAME        = "neofiles"
DESCRIPTION = "NEO file handling"

import os

# The following libs must be installed with pip

# AstroPy & friends
from astropy.time import Time

# Local modules
from neoconfig import config
from verbose import warning



now = Time.now()
prefix = now.strftime("%Y%m%d")
neo_obs_data_dir = config.neo_obs_data_dir

# If subdir doesn't exist, create it, NINA will create it only with the 1st frame
if not os.path.isdir(neo_obs_data_dir):
    warning(f"directory {neo_obs_data_dir} doesn't exist, creating it")
    os.makedirs(neo_obs_data_dir)



def set_prefix(new_prefix: str) -> None:
    now = Time.now()
    prefix = new_prefix



def path(filename: str) -> str:
    return os.path.join(neo_obs_data_dir, f"{prefix}-{filename}")
