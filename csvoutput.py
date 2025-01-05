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
#   from csvoutput import csv_output
#   csv_output.set_default_locale(locale="")
#   csv_output.set_float_format(fmt="%.3f")
#   csv_output.add_row([a, b, c, ...])
#   csv_output.add_fields([name1, name2, ...])
#   csv_output.write(file="", set_locale=True)       file="" uses stdout
#   csv_output(a, b, c, ...)
#   csv_output(row=[a, b, c, ...])
#   csv_output(fields=[a, b, c, ...])
#
# Compatibility with version 1.0
#   from csvoutput import csv_output as CSVOutput

# ChangeLog
# Version 1.0 / 2024-07-12
#       locale aware global CSV output class
# Version 2.0 / 2024-11-21
#       Reworked as a proper csv_output object, added new interface


import csv
import locale
import sys
import typing

VERSION = "2.0 / 2024-11-21"
AUTHOR  = "Martin Junius"
NAME    = "csvoutput"


DEFAULT_FLOAT_FORMAT = "%.3f"

class CSVOutput:
    """CSV output class"""

    def __init__(self):
        self._cache      = []
        self._fields     = None
        self._float_fmt  = None    # format for floats if set


    def __call__(self, *args, **kwargs):
        fields = kwargs.get("fields")
        if fields:
            self.add_fields(fields)
        row = list(args) or kwargs.get("row")
        if row:
            self.add_row(row)


    def set_default_locale(self, loc: str=""):
        locale.setlocale(locale.LC_ALL, loc)
        if locale.localeconv()['decimal_point'] == ",":
            # Set float format to automatically format all float values as strings using locale
            if not self._float_fmt:
                self.set_float_format()


    def set_float_format(self, fmt: str=DEFAULT_FLOAT_FORMAT):
        self._float_fmt = fmt


    def _fmt(self, v: float):
        return locale.format_string(self._float_fmt, v)


    def add_row(self, data: list):
        self._cache.append(data)


    def add_fields(self, fields: list):
        self._fields = fields


    def _write(self, f: typing.TextIO):            
        if locale.localeconv()['decimal_point'] == ",":
            # Use ; as the separator and quote all fields for easy import in "German" Excel
            writer = csv.writer(f, dialect="excel", delimiter=";", quoting=csv.QUOTE_ALL)
        else:
            writer = csv.writer(f, dialect="excel")
        if self._fields:
            writer.writerow(self._fields)
        
        # Write rows one by one the enable conversion of float values
        for row in self._cache:
            if self._float_fmt:
                row = [ self._fmt(v) if type(v) == float else v   for v in row ]
            writer.writerow(row)


    def write(self, file: str=None, set_locale: bool=True):
        if set_locale:
            self.set_default_locale()

        if file:
            with open(file, 'w', newline='', encoding="utf-8") as f:
                self._write(f)
        else:
                self._write(sys.stdout)



# Global object
csv_output = CSVOutput()
