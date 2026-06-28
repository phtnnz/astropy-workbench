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
# Version 0.1 / 2026-06-21
#       Query previous NEOCP objects

VERSION     = "0.1 / 2026-06-21"
AUTHOR      = "Martin Junius"
NAME        = "mpc.prevneocp"
DESCRIPTION = "Retrieve previous NEOCP from MPC"

import re
import requests
from icecream import ic
# Disable debugging
ic.disable()

# Local modules
from utils.verbose import verbose, warning, error, message
from neo.config import config
from neo.classes import PrevNEOCPData



def mpc_query_prev_neocp(url: str) -> str:
    ic(url)

    timeout = config.requests_timeout

    ic(url)
    verbose(f"query {url}")
    response = requests.get(url, timeout=timeout)
    ic(response.status_code)
    if response.status_code != 200:
        error(f"query to {url} failed")

    return response.text



def mpc_parse_prev_neocp(content: str) -> dict[str, PrevNEOCPData]:
    objects = dict()

    in_prev = False
    for line in content.splitlines():
        m = re.match(r'<h2>(.*)</h2>', line)
        if m:
            ic(m)
            h2 = m.group(1)
            in_prev = True if h2 == "Previous NEO Confirmation Page Objects" else False
            continue
        if line.startswith("</ul>"):
            in_prev = False
        if not in_prev or not line.startswith("<li>"):
            continue

        # <li> 2019 NG<sub>2</sub> = C1FFEF5 (June 21.43 UT)
        # <li> 2026 MH<sub>1</sub> = P12nxW2 (June 21.34 UT)   [see <a href="/mpec/K26/K26M57.html"><i>MPEC</i> 2026-M57</a>]
        m = re.match(r'<li>\s*(.+?)\s+=\s+(\w+?)\s+\((.+?)\)\s*(\[.+\])?$', line)
        if m:
            prov_id, trk_sub, date, mpec = m.groups()
            prov_id = prov_id.replace("<sub>", "").replace("</sub>", "")
            mpec_url = None
            mpec_no  = None
            if mpec:
                m1 = re.search(r'href="(.+?)".*</i>\s*(.+)</a>', mpec)
                if m1:
                    mpec_url = config.mpc_url + m1.group(1)
                    mpec_no  = m1.group(2)
            ic(trk_sub, prov_id, date, mpec_no, mpec_url)
            objects[trk_sub] = PrevNEOCPData(trk_sub, prov_id, None, date, mpec_no, mpec_url)
            continue

        # <li> CE74892 was not confirmed (Feb. 21.55 UT)
        # <li> A11yweM was not a minor planet (Feb. 20.95 UT)
        m = re.match(r'<li>\s*(\w+?)\s+([\w ]+?)\s+\((.+?)\)$', line)
        if m:
            trk_sub, comment, date = m.groups()
            ic(trk_sub, comment, date)
            objects[trk_sub] = PrevNEOCPData(trk_sub, None, comment, date)
            continue
            
        if not m:
            warning("can't match line:", line)


    return objects
