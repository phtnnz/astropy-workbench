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
# Version 0.1 / 2023-12-18
#       First version of JSONConfig module
# Version 0.2 / 2024-01-23
#       Rewritten for multiple config files, global config object
#
#       Usage:  from jsonconfig import config
#               from jsonconfig import JSONConfig
#               class MyConfig(JSONConfig)
# Version 0.3 / 2024-06-19
#       Search also in .config
# Version 0.4 / 2024-07-23
#       Added code for getting the path of Windows' Documents directory
# Version 0.5 / 2024-08-28
#       New method .info(), verbose output of config filename and top-level keys
# Version 0.6 / 2025-06-29
#       Method .get() now supports nested keys, e.g. config.get("main", "sub", "setting")

import os
import sys
import argparse
import json
# Windows specific
import ctypes.wintypes

# The following libs must be installed with pip
from icecream import ic
# Disable debugging

# Local modules
from verbose import verbose, warning, error



VERSION = "0.6 / 2025-06-29"
AUTHOR  = "Martin Junius"
NAME    = "JSONConfig"


CONFIG     = ".config"
CONFIGDIR  = "astro-python"
CONFIGFILE = "astro-python-config.json"

#ic.enable()
ic(CONFIGDIR, CONFIGFILE)



class JSONConfig:
    """ JSONConfig base class """

    def __init__(self, file, warn=True, err=True):
        ic("config init", file)
        self.config = {}
        self.read_config(file, warn, err)


    def read_config(self, file, warn=True, err=True):
        ic(file)
        file1 = self.search_config(file)
        if(file1):
            self.configfile = file1
            json  = self.read_json(file1)
            # Merge with existing config
            self.config = self.config | json
        elif err:
            error(f"config {file} not found")
        elif warn:
            warning(f"config {file} not found")


    def info(self):
        verbose(f"config file {self.configfile}")
        verbose("config keys:", " ".join( [k for k in self.config.keys() if not k.startswith("#")] ))

    def search_config(self, file):
        # If full path use as is
        if os.path.isfile(file):
            return file
        
        # Search config file in current directory, LOCALAPPDATA, APPDATA
        searchpath = []

        path = os.path.join(os.path.curdir, CONFIG)
        if os.path.isdir(path):
            searchpath.append(path)

        path = os.path.join(os.path.curdir, CONFIG, CONFIGDIR)
        if os.path.isdir(path):
            searchpath.append(path)

        appdata = os.environ.get('LOCALAPPDATA')
        if appdata:
            path = os.path.join(appdata, CONFIGDIR)
            if os.path.isdir(path):
                searchpath.append(path)
        else:
            ic("environment LOCALAPPDATA not set!")

        appdata = os.environ.get('APPDATA')
        if appdata:
            path = os.path.join(appdata, CONFIGDIR)
            if os.path.isdir(path):
                searchpath.append(path)
        else:
            ic("environment APPDATA not set!")

        # Add Python search path to list
        searchpath.extend([ os.path.join(d, CONFIG) for d in sys.path ])

        # Search for config file in searchpath list
        for path in searchpath:
            ic(path)
            file1 = os.path.join(path, file)
            if os.path.isfile(file1):
                ic(file1)
                return file1

        return None


    def read_json(self, file):
        with open(file, 'r') as f:
            return json.load(f)


    def write_json(self, file):
        with open(file, 'w') as f:
            json.dump(self.config, f, indent = 2)


    def get(self, *keys):
        cf = self.config
        for k in keys:
            cf = cf.get(k)
            if not cf:
                return None
        return cf


    def get_keys(self):
        return self.config.keys()


    def get_json(self):
        # For backwards compatibility
        return self.config


    def get_documents_path(self):
        # Windows hack to get path of Documents folder, which might reside on other drives than C:
        CSIDL_PERSONAL = 5       # My Documents
        SHGFP_TYPE_CURRENT = 0   # Get current, not default value
        buf = ctypes.create_unicode_buffer(ctypes.wintypes.MAX_PATH)
        ctypes.windll.shell32.SHGetFolderPathW(None, CSIDL_PERSONAL, None, SHGFP_TYPE_CURRENT, buf)
        return buf.value



# Global config object
config = JSONConfig(CONFIGFILE, False, False)



def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = "Test for module",
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="debug messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("-c", "--config", help="read CONFIG file")

    args = arg.parse_args()

    verbose.set_prog(NAME)
    if args.verbose:
        verbose.enable()
    if args.debug:
        ic.enable()
        ic(sys.version_info, sys.path)
    if args.config:
        config.read_config(args.config)

    print("JSON config keys =", ", ".join(config.get_keys()))
    print("Documents path =", config.get_documents_path())



if __name__ == "__main__":
    main()
