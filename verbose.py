#!/usr/bin/env python

# Copyright 2023-2026 Martin Junius
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
# Version 0.1 / 2023-11-04
#       First version of verbose module
# Version 0.2 / 2023-12-18
#       Added warning(), error() with abort
# Version 1.0 / 2024-01-06
#       Version bumped to 1.0
# Version 1.1 / 2024-08-08
#       Output "exiting" message only if verbose is enabled
#       errno is now global for verbose, warning, error
# Version 1.2 / 2024-12-15
#       Added docstrings
# Version 1.3 / 2024-12-16
#       Added message(), just like print(), but can be disabled
# Version 1.4 / 2025-12-17
#       Added .print_lines(), splits multi-line string representations
# Version 2.0 / 2026-02-20
#       Added context manager and logfile output
#
#       Usage:  from verbose import message, verbose, warning, error
#               message(print-like-args)
#               verbose(print-like-args)
#               warning(print-like-args)
#               error(print-like-args)
#               .print_lines(line(s), ...)
#               .enable(flag=True)
#               .disable()
#               .enabled
#               .set_prog(name)         global for all objects
#               .set_errno(errno)       relevant only for error()
#
#               with verbose.logfile(LOGFILE):
#                   ...

VERSION = "2.0 / 2026-02-20"
AUTHOR  = "Martin Junius"
NAME    = "verbose"

import argparse
import sys
from typing import TextIO
from contextlib import contextmanager



class Verbose:
    """Class for verbose-stype log output objects"""
    _progname: str = ""          # global program name
    _log_file: TextIO = None     # log file
    _errno: int = 1              # exit code, 1 for generic errors

    def __init__(self, flag: bool=False, prefix: str=None, abort: bool=False):
        """Constructor

        Parameters
        ----------
        flag : bool, optional
            Output enabled, by default False
        prefix : str, optional
            Output prefix, by default None
        abort : bool, optional
            Error abort after output, by default False
        """
        self.enabled = flag
        self.prefix = prefix
        self.abort = abort


    def __call__(self, *args, **kwargs):
        """Make Verbose objects callable like print()
        """
        if not self.enabled:
            return
        preargs = ()
        if Verbose._progname:
            preargs = (f"{Verbose._progname}:", )
        if self.prefix:
            preargs = preargs + (f"{self.prefix}:", )
        if preargs:
            args = preargs + args
        print(*args, **kwargs)
        if Verbose._log_file:
            print(*args, file=Verbose._log_file, **kwargs)
        if self.abort:
            self._exit()


    def print_lines(self, *args, **kwargs) -> None:
        """Multi-line output using default string representation
        """
        for arg in args:
            for line in str(arg).splitlines():
                self.__call__(line, **kwargs)


    def enable(self, flag: bool=True):
        """Enable output

        Parameters
        ----------
        flag : bool, optional
            Enable flag, by default True
        """
        self.enabled = flag


    def disable(self):
        """Disable output
        """
        self.enabled = False


    def set_prog(self, name: str=""):
        """Set program name for output

        Parameters
        ----------
        name : str, optional
            Program name, by default ""
        """
        Verbose._progname = name


    def set_errno(self, errno: int):
        """Set errno for error exit

        Parameters
        ----------
        errno : int
            System error number
        """
        Verbose._errno = errno


    def _exit(self):
        """Internal: exit script
        """
        if verbose.enabled:
            if Verbose._progname:
                print(Verbose._progname + ": ", end="")
            print(f"exiting ({Verbose._errno})")
        sys.exit(Verbose._errno)


    @contextmanager
    def logfile(self, filename: str):
        """Context manager for addtional log file output

        Parameters
        ----------
        filename : str
            Log file name
        """
        try:
            Verbose._log_file = open(filename, mode="w")
            yield None
        finally:
            Verbose._log_file.close()
            Verbose._log_file = None


"""Global callable objects
"""
message = Verbose(True)
verbose = Verbose()
warning = Verbose(True, "WARNING")
error   = Verbose(True, "ERROR", True)



# Test
def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = "Test script for verbose module",
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")

    args = arg.parse_args()

    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()
    if args.debug:
        message("Nothing to debug ;-)")

    message("Just", "a", "message at the beginning")
    warning("Just a first warning ;-)")
    verbose("Test", "1", "for", "verbose()")
    verbose("Test", "2", "for more", "verbose()", "with some formatting {:04d}".format(11+12))
    verbose("Changing progname")
    verbose.set_prog(NAME+"2")
    warning("A", "warning", "message", " --- but no abort here!")
    warning.set_prog(NAME+"3")
    warning("Another change to progname occurred")
    verbose.print_lines("abc", "def", "line1\nline2\nline3")

    test_log = "_verbose_test.log"
    verbose.set_prog(NAME)
    with verbose.logfile(test_log):
        verbose(f"writing to log file {test_log}")
        verbose("yet another line in log file")
        warning("Warning in log file")
        verbose("end of log file")
    error.set_errno(99)
    error("Error test", "for Verbose module")

    

if __name__ == "__main__":
    main()
