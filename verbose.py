#!/usr/bin/env python

# Copyright 2023-2024 Martin Junius
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
#
#       Usage:  from verbose import message, verbose, warning, error
#               message(print-like-args)
#               verbose(print-like-args)
#               warning(print-like-args)
#               error(print-like-args)
#               .enable(flag=True)
#               .disable()
#               .enabled
#               .set_prog(name)         global for all objects
#               .set_errno(errno)       relevant only for error()

import argparse
import sys

VERSION = "1.3 / 2024-12-16"
AUTHOR  = "Martin Junius"
NAME    = "verbose"



class Verbose:
    """
    Class for verbose-style objects
    """
    progname = None             # global program name
    errno    = 1                # exit code, 1 for generic errors

    def __init__(self, flag: bool=False, prefix: str=None, abort: bool=False):
        """
        Create verbose-style object

        :param flag: enable output flag, defaults to False
        :type flag: bool, optional
        :param prefix: output prefix (in addition to program name), defaults to None
        :type prefix: str, optional
        :param abort: abort after output flag, defaults to False
        :type abort: bool, optional
        """
        self.enabled = flag
        self.prefix = prefix
        self.abort = abort

    def __call__(self, *args, **kwargs):
        """
        Make verbose-style object callable, all parameters passed to print()
        """
        if not self.enabled:
            return
        if Verbose.progname:
            print(Verbose.progname + ": ", end="")
        if self.prefix:
            print(self.prefix + ": ", end="")
        print(*args, **kwargs)
        if self.abort:
            self._exit()

    def enable(self, flag: bool=True):
        """
        Enable (default) or disable (flag=False) output

        :param flag: enable output flag, defaults to True
        :type flag: bool, optional
        """
        self.enabled = flag

    def disable(self):
        """
        Disable output
        """
        self.enabled = False

    def set_prog(self, name: str):
        """
        Set program name prefix

        :param name: program name
        :type name: str
        """
        Verbose.progname = name

    def set_errno(self, errno: int):
        """
        Set global errno for abort exit()

        :param errno: error code
        :type errno: int
        """
        Verbose.errno = errno

    def _exit(self):
        """
        Internal, exit program
        """
        if verbose.enabled:
            if Verbose.progname:
                print(Verbose.progname + ": ", end="")
            print(f"exiting ({Verbose.errno})")
        sys.exit(Verbose.errno)


message = Verbose(True)
verbose = Verbose()
warning = Verbose(True, "WARNING")
error   = Verbose(True, "ERROR", True)




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
    error.set_errno(99)
    error("Error test", "for Verbose module")

    

if __name__ == "__main__":
    main()