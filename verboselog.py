#!/usr/bin/env python

# Copyright 2023-2025 Martin Junius
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
# Version 0.1 / 2026-01-02
#       Based on verbose 1.4, using logging for output
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

import argparse
import sys
import logging

VERSION = "0.1 / 2026-01-02"
AUTHOR  = "Martin Junius"
NAME    = "verboselog"



class Verbose:
    """
    Class for verbose-style objects using logging
    """
    progname: str = ""              # global program name
    errno: int    = 1               # exit code, 1 for generic errors
    logger: logging.Logger = None   # global logger object

    def __init__(self, flag: bool, level: int=logging.INFO, abort: bool=False):
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
        self.level = level
        self.abort = abort


    def __call__(self, *args, **kwargs):
        """
        Make verbose-style object callable, all parameters passed to print()
        """
        if not self.enabled:
            return

        # Make sure that logger exists
        if not Verbose.logger:
            self.set_prog()

        # print(*args, **kwargs)
        msg = " ".join(args)
        Verbose.logger.log(self.level, msg, **kwargs)
        if self.abort:
            self._exit()


    def print_lines(self, *args, **kwargs) -> None:
        """
        Print multi-line string representation of object(s)
        """
        for arg in args:
            for line in str(arg).splitlines():
                self.__call__(line, **kwargs)


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

    def set_prog(self, name: str=""):
        """
        Set program name prefix

        :param name: program name
        :type name: str
        """
        if name:
            format = "%(asctime)s %(name)s:%(levelname)s: %(message)s"
        else:
            format = "%(asctime)s %(levelname)s: %(message)s"
        logging.basicConfig(encoding='utf-8', level=logging.DEBUG,
                            format=format, datefmt='%Y-%m-%d %H:%M:%S')
        Verbose.logger = logging.getLogger(name)
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
        Verbose.logger.error(f"exiting ({Verbose.errno})")
        sys.exit(Verbose.errno)


message = Verbose(True, logging.INFO)
verbose = Verbose(False, logging.INFO)
warning = Verbose(True, logging.WARNING)
error   = Verbose(True, logging.ERROR, True)



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
    verbose("Test", "2", "for more", "verbose()", f"with some formatting {11+12=:04d}")
    verbose("Changing progname")
    verbose.set_prog(NAME+"2")
    warning("A", "warning", "message", "--- but no abort here!")
    warning.set_prog(NAME+"3")
    warning("Another change to progname occurred")
    verbose.print_lines("abc", "def", "line1\nline2\nline3")
    error.set_errno(99)
    error("Error test", "for Verbose module")

    

if __name__ == "__main__":
    main()