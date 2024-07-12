#!/usr/bin/env python

# Copyright 2023 Martin Junius
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
# Version 0.0 / 2024-07-12
#       locale aware global CSV output class

import csv
import locale

VERSION = "0.0 / 2024-07-12"
AUTHOR  = "Martin Junius"
NAME    = "csvoutput"



class CSVOutput:
    _cache  = []
    _fields = None
    _locale = None

    def set_default_locale():
        locale.setlocale(locale.LC_ALL, "")

    def add_csv_row(data):
        CSVOutput._cache.append(data)

    def add_csv_fields(fields):
        CSVOutput._fields = fields

    def write_csv(file):
        with open(file, 'w', newline='', encoding="utf-8") as f:
            ##FIXME: check  locale.RADIXCHAR
            if locale.localeconv()['decimal_point'] == ",":
                # Use ; as the separator and quote all fields for easy import in "German" Excel
                writer = csv.writer(f, dialect="excel", delimiter=";", quoting=csv.QUOTE_ALL)
            else:
                writer = csv.writer(f, dialect="excel")
            if CSVOutput._fields:
                writer.writerow(CSVOutput._fields)
            writer.writerows(CSVOutput._cache)



def val(val):
    return locale.format_string("%.3f", val)

