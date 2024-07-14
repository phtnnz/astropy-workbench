#!/usr/bin/env python

# Copyright 2024 Martin Junius
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

# Usage
#   from csvoutput import CSVOutput
#   CSVOutput.set_default_locale(locale="")
#   CSVOutput.set_float_format(fmt="%.3f")
#   CSVOutput.add_row(list)
#   CSVOutput.add_fields(list)
#   CSVOutput.write(file="", set_locale=True)       file="" uses stdout

# ChangeLog
# Version 0.0 / 2024-07-12
#       locale aware global CSV output class


import csv
import locale
import sys

VERSION = "0.0 / 2024-07-12"
AUTHOR  = "Martin Junius"
NAME    = "csvoutput"



DEFAULT_FLOAT_FORMAT = "%.3f"

class CSVOutput:
    _cache      = []
    _fields     = None
    _float_fmt  = None    # format for floats if set


    def set_default_locale(loc=""):
        locale.setlocale(locale.LC_ALL, loc)
        if locale.localeconv()['decimal_point'] == ",":
            # Set float format to automatically format all float values as strings using locale
            if not CSVOutput._float_fmt:
                CSVOutput.set_float_format()


    def set_float_format(fmt=DEFAULT_FLOAT_FORMAT):
        CSVOutput._float_fmt = fmt


    def _fmt(v):
        return locale.format_string(CSVOutput._float_fmt, v)


    def add_row(data):
        CSVOutput._cache.append(data)


    def add_fields(fields):
        CSVOutput._fields = fields


    def _write(f):            
        if locale.localeconv()['decimal_point'] == ",":
            # Use ; as the separator and quote all fields for easy import in "German" Excel
            writer = csv.writer(f, dialect="excel", delimiter=";", quoting=csv.QUOTE_ALL)
        else:
            writer = csv.writer(f, dialect="excel")
        if CSVOutput._fields:
            writer.writerow(CSVOutput._fields)
        
        # Write rows one by one the enable conversion of float values
        for row in CSVOutput._cache:
            if CSVOutput._float_fmt:
                row = [ CSVOutput._fmt(v) if type(v) == float else v   for v in row ]
            writer.writerow(row)


    def write(file=None, set_locale=True):
        if set_locale:
            CSVOutput.set_default_locale()

        if file:
            with open(file, 'w', newline='', encoding="utf-8") as f:
                CSVOutput._write(f)
        else:
                CSVOutput._write(sys.stdout)
