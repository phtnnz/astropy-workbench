#!/usr/bin/env python

# Copyright 2024-2025 Martin Junius
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
#   csv_output.set_force_float(force=True)
#   csv_output.add_row([a, b, c, ...])
#   csv_output.add_fields([name1, name2, ...])
#   csv_output.write(file="", set_locale=True)       file="" uses stdout
#   csv_output(a, b, c, ...)
#   csv_output([a, b, c, ...])
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
# Version 2.1 / 2024-12-16
#       Added docstrings
# Version 2.2 / 2025-11-21
#       Changed csv_output() to accepts both tuples and lists, 
#       added .set_force_float() to convert any string that looks
#       like a float to a float and formats it according to locale
# Version 2.3 / 2025-11-21
#       Converted to numpy docstring format

import csv
import locale
import sys
from typing import TextIO, Any

VERSION = "2.3 / 2025-11-21"
AUTHOR  = "Martin Junius"
NAME    = "csvoutput"


DEFAULT_FLOAT_FORMAT = "%f"

class CSVOutput:
    """
    CSV output class
    """

    def __init__(self) -> None:
        """
        Create CSV output object
        """
        self._cache       = []
        self._fields      = None
        self._float_fmt   = None    # format for floats if set
        self._force_float = False


    def __call__(self, *args, **kwargs) -> None:
        """
        Make CSV output object callable

        Default behavior like .add_row()

        Optional keyword args
        row=[ data1, data2, date3, ... ]
        fields=[ name1, name2, name3, ... ]
        """
        if args:
            if isinstance(args[0], list):
                self.add_row(args[0])
            else:
                self.add_row(list(args))
        fields = kwargs.get("fields")
        if fields:
            self.add_fields(fields)
        row = kwargs.get("row")
        if row:
            self.add_row(row)
        


    def set_default_locale(self, loc: str="") -> None:
        """
        Set default locale for CSV output

        Parameters
        ----------
        loc : str, optional
            locale name, by default "" = system locale
        """
        locale.setlocale(locale.LC_ALL, loc)
        if locale.localeconv()['decimal_point'] == ",":
            # Set float format to automatically format all float values as strings using locale
            if not self._float_fmt:
                self.set_float_format()


    def set_force_float(self, force: bool=True) -> None:
        """
        Set force to float conversion flag

        Parameters
        ----------
        force : bool
            Force flag
        """
        self._force_float = force


    def set_float_format(self, fmt: str=DEFAULT_FLOAT_FORMAT) -> None:
        """
        Set format for float numbers

        Parameters
        ----------
        fmt : str, optional
            %-style format string, by default DEFAULT_FLOAT_FORMAT
        """
        self._float_fmt = fmt


    def _fmt(self, v: Any) -> str:
        """
        Format fields in row, optionally forcing interpretation as float

        Parameters
        ----------
        v : Any
            Value

        Returns
        -------
        str
            Formatted value as string
        """
        if self._force_float:
            try:
                return locale.format_string(self._float_fmt, float(v))
            except ValueError:
                return str(v)
        elif isinstance(v, float):
            return locale.format_string(self._float_fmt, float(v))
        else:
            return str(v)
    

    def add_row(self, data: list) -> None:
        """
        Add data row to CSV output

        Parameters
        ----------
        data : list
            List of values [a, b, c, ...]
        """
        self._cache.append(data)


    def add_fields(self, fields: list) -> None:
        """
        Add field names header to CSV output

        Parameters
        ----------
        fields : list
            List of field names [txt1, txt2, ...]
        """
        self._fields = fields


    def _write(self, f: TextIO) -> None:
        """
        Internal, use csv.writer to output CSV data

        Parameters
        ----------
        f : TextIO
            File descriptor
        """
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
                row = [ self._fmt(v) for v in row ]
            writer.writerow(row)


    def write(self, file: str=None, set_locale: bool=True) -> None:
        """
        Write CSV output to file or stdout

        Parameters
        ----------
        file : str, optional
            Filename, by default None = write to stdout
        set_locale : bool, optional
            Automatically set default system locale, by default True
        """
        if set_locale:
            self.set_default_locale()

        if file:
            with open(file, 'w', newline='', encoding="utf-8") as f:
                self._write(f)
        else:
                self._write(sys.stdout)



"""
Global CSVOutput object
"""
csv_output = CSVOutput()
