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
# Version 0.1 / 2023-06-24
#       Added -n option, read filter name from CSV
# Version 0.2 / 2023-06-26
#       Use new target template "Target NEO v2.json"
#       Use new base template "NEW v4 nautical.json"
#       Works with nautical dusk/dawn time, too, see caveat below
# Version 0.3 / 2023-07-01
#       Added -N option, append number of frames to target name
# Version 0.4 / 2023-07-03
#       New target template "./NINA-Templates-IAS/Base NEO nautical.json"
#       New base template "./NINA-Templates-IAS/Base NEO nautical.json"
#       Removed -l option
# Version 1.0 / 2023-07-04
#       Cleaner handling of SelectedProvider references in completed sequence
#       Process External Script with TARGET
#       General clean-up
# Version 1.1 / 2024-06-28
#       Changes for Remote3, added -3 --remote3 option
# --- nina-create-sequence2 ---
# Version 1.2 / 2024-07-25
#       Use JSONConfig, verbose modules
#       Removed a lot of options (now handled in config)
#       New option -A --debug-print-attr
#       New config setting "container", add target items to this container in target area
# Version 1.3 / 2024-08-06
#       Use new radec.Coord
#       More alternatives and defaults for CSV data field names
# Version 1.4 / 2024-09-02
#       Added --date option, default today, use this date for output and template format
#       New "subdir" config setting, see below
# Version 1.5 / 2025-08-20
#       Added -l --list-targets option, just output the list of targets
# Version 1.6 / 2025-09-27
#       Added "autofocus_first_target_only" setting to config


# See here https://www.newtonsoft.com/json/help/html/SerializingJSON.htm for the JSON serializing used in N.I.N.A

# Entries in nina-create-sequence.json:
# "<NAME>": {
#     "timezone":  "<IANA TIMEZONE NAME>"
#     "template":  "<BASE SEQUENCE TEMPLATE>",
#     "target":    "<SINGLE TARGET TEMPLATE>",
#     "container": "<CONTAINER NAME OR EMPTY>",
#     "format":    "<TARGETNAME {x}>",
#     "output":    "<OUTPUT FILENAME {x}>.json",
#     "subdir":    "_asteroids_{1}",            (optional)
#     "autofocus_first_target_only": "yes"      (optional)
# }
#
# Format placeholders:
# 0=target, 1=date, 2=seq, 3=number
# "subdir" is optional, or can be blank ""

VERSION = "1.6 / 2025-09-27"
AUTHOR  = "Martin Junius"
NAME    = "nina-create-sequence2"

import sys
import os
import argparse
import json
import csv
from datetime import datetime, timezone, timedelta, date
from zoneinfo import ZoneInfo
import copy

# The following libs must be installed with pip
import tzdata
from icecream import ic
# Disable debugging
ic.disable()

# Local modules
from verbose          import verbose, warning, error, message
from jsonconfig       import JSONConfig, config
from radec            import Coord



DEFAULT_FILTER_NAMES = [ "L", "R", "G", "B", "Ha", "OIII", "SII"]


CONFIG = "nina-create-sequence.json"

class SequenceConfig(JSONConfig):
    """ JSON Config for creating NINA sequences """

    def __init__(self, file=None):
        super().__init__(file)

  

config = SequenceConfig(CONFIG)



class Options:
    """ Command line options """
    debug_print_attr = False            # -A --debug-print-attr
    date = date.today().isoformat()     # --date
    list_targets = False                # -l --list-targets



class TargetData:
    """ Hold data to update N.I.N.A template """

    def __init__(self, name: str, target: str, coord: Coord, time: datetime, 
                 number: int, exposure: float, filter: str="L", binning: str="2x2", subdir: str=None,
                 no_autofocus: bool=False):
        self.name = name
        self.targetname = target
        self.coord = coord
        self.time = time
        self.number = number
        self.exposure = exposure
        self.filter = filter
        self.binning = 0
        if str(binning).startswith("1"):
            self.binning = 1
        if str(binning).startswith("2"):
            self.binning = 2
        self.subdir = subdir
        self.no_autofocus = no_autofocus



class NINABase:
    """ Base class for NINASequence and NINATarget """

    def __init__(self):
        self.obj = None


    def read_json(self, file):
        with open(file, 'r') as f:
            data = json.load(f)
        self.obj = data


    def write_json(self, file):
        with open(file, 'w') as f:
            json.dump(self.obj, f, indent = 2)


    def traverse(self, func=None, param=None):
        self.traverse_obj(self.obj, ">", 1, func, param)


    def traverse_obj(self, obj, indent, level, func=None, param=None):
        if Options.debug_print_attr:
            print(indent, "KEYS =", ", ".join(obj.keys()))

        if func:
            func(self, obj, indent + ">", param)

        for k, val in obj.items():
            if k=="$id" or k=="$type" or k=="$ref" or k=="Name":
                self.print_attr(obj, k, indent+" >")

            if type(val) is dict:
                """dict"""
                self.traverse_obj(val, indent + " >", level + 1, func, param)
            elif type(val) is str:
                """str"""
            elif type(obj[k]) is list:
                """list"""
                for val1 in val:
                    self.traverse_obj(val1, indent + " >", level + 1, func, param)
            else:
                """rest"""


    def set_prefix(self, prefix):
        self.id_prefix = prefix


    def add_prefix_to_id(self, obj, indent, param=None):
        # change "$id" and "$ref"
        self.__add_prefix(obj, "$id")
        self.__add_prefix(obj, "$ref")


    def __add_prefix(self, obj, key):
        if key in obj:
            # obj[key] = str(self.id_prefix*1000 + int(obj[key]))
            # JSON serializing in NINA uses pure numeric ids, but strings
            # work as well and are collision free
            obj[key] = f"{self.id_prefix:04d}_{int(obj[key]):04d}"


    def process_provider(self, obj, indent, dict):
        # search for Provider {...} and change additional occurences to reference
        self.print_attr(obj, "SelectedProvider", "")
        if "SelectedProvider" in obj.keys():
            prov = obj["SelectedProvider"]
            if "$type" in prov.keys():
                type = prov["$type"]
                id   = prov["$id"]
                if type in dict.keys():
                    # already exists
                    ref = dict[type]
                    obj["SelectedProvider"] = { "$ref": ref }

                else:
                    # 1st occurence, don't touch
                    dict[type] = id
        self.print_attr(obj, "SelectedProvider", "")


    def print_attr(self, obj, name, indent):
        if Options.debug_print_attr:
            if name in obj:
                print(indent, name, "=", obj[name])




class NINATarget(NINABase):
    """ Holds data read from N.I.N.A JSON template for target """

    def __init__(self):
        self.name            = None # template name
        self.targetname      = None # astro target name
        # instance variables created by process_data(), referencing the data objects, used by update_target_data()
        self.target          = None
        self.coord           = None
        self.waitfortime     = None
        self.conditions0     = None
        self.filter          = None
        self.exposure        = None
        self.binning         = None
        self.timecondition   = None
        self.script_w_target = None
        self.autofocus       = None
        # Instance variables created by process_data(), referencing the inner containers
        # Items in container
        # [0] Pre-imaging (slew, AF, WaitForTime)
        # [1] Imaging (exposure loop)
        # [2] Post-imaging (create flag &c.)
        self.container_items = None
        self.container_conditions = None
        self.container_triggers = None

        NINABase.__init__(self)


    def update_target_data(self, data: TargetData):
        # NEW NAME --> obj["Name"]
        self.obj["Name"] = data.name
        self.name = data.name

        # NEW TARGET --> target["TargetName"]
        self.target["TargetName"] = data.targetname

        # NEW COORD --> coord["..."]
        # Coordinates in WaitForAltitude don't need to be updated, handled automatically by N.I.N.A when loading
        self.coord["RAHours"]   = int(data.coord.ra_h)
        self.coord["RAMinutes"] = int(data.coord.ra_m)
        self.coord["RASeconds"] = float(data.coord.ra_s)
        # NINA: for negative DEC values, both NegativeDec is set to true *and* DecDegrees is negative!
        #       DEC = -01:29:39 -> "NegativeDec": true, "DecDegrees": -1
        #       DEC = -00:29:39 -> "NegativeDec": true, "DecDegrees": 0
        self.coord["NegativeDec"] = data.coord.dec_neg
        self.coord["DecDegrees"] = int(data.coord.dec_d) * data.coord.dec_sign
        self.coord["DecMinutes"] = int(data.coord.dec_m)
        self.coord["DecSeconds"] = float(data.coord.dec_s)

        # NEW TIME --> waitfortime["..."]
        if self.waitfortime:
            if not data.time:
                error("WaitForTime in sequence, but no start time in CSV data")
            self.waitfortime["Hours"]   = int(data.time.hour)
            self.waitfortime["Minutes"] = int(data.time.minute)
            self.waitfortime["Seconds"] = int(data.time.second)

        # NEW NUMBER OF EXPOSURES --> conditions0["Iterations"]
        self.conditions0["Iterations"] = int(data.number)

        # NEW FILTER --> filter["_name"]
        self.filter["_name"] = data.filter

        # NEW EXPOSURETIME --> exposure["ExposureTime"]
        self.exposure["ExposureTime"] = float(data.exposure)
        # NEW BINNING --> exposure["Binning"]["X|Y"]
        if data.binning > 0:
            self.exposure["Binning"]["X"] = data.binning
            self.exposure["Binning"]["Y"] = data.binning

        # Update target for External Script
        if self.script_w_target:
            ## FIXME: quick'n'dirty to work with new _asteroids_YYYY-MM-DD subdirectories
            ##        replace nina-flag-ready.bat         
            script = self.script_w_target["Script"]
            newtarget = os.path.join(data.subdir, data.targetname)
            self.script_w_target["Script"] = script.replace("\"TARGET\"", f"\"{newtarget}\"")

        # Disable autofocus
        if self.autofocus:
            if data.no_autofocus:
                verbose("skipping autofocus")
                # hack, convert to Annotation, del won't work because this needs the container
                # and we only have the item here
                self.autofocus["$type"] = "NINA.Sequencer.SequenceItem.Utility.Annotation, NINA.Sequencer"
                self.autofocus["Text"] = f"{NAME}: removed Run Autofocus"
            ic(self.autofocus)
            pass



    def process_data(self):
        self.name = self.obj["Name"]
        verbose("NINATarget(process_data):", "name =", self.name)

        self.target = self.obj["Target"]
        self.targetname = self.target["TargetName"]
        self.coord = self.target["InputCoordinates"]

        items = self.obj["Items"]

        self.container_items      = []
        self.container_conditions = []
        self.container_triggers   = []

        for item in items["$values"]:
            type = item["$type"]
            # print("type =", type)
            if "Container.SequentialContainer" in type:
                self.container_items.append(item["Items"]["$values"])
                self.container_conditions.append(item["Conditions"]["$values"])
                self.container_triggers.append(item["Triggers"]["$values"])

        self.__process_container_list()


    def __process_container_list(self):
        for container in self.container_items + self.container_conditions + self.container_triggers:
            for item in container:
                if "WaitForTime" in item["$type"]:
                    self.waitfortime = item

                if "SmartExposure" in item["$type"]:
                    # NEW NUMBER OF EXPOSURES --> conditions0["Iterations"]
                    self.conditions0 = item["Conditions"]["$values"][0]

                    items = item["Items"]
                    for item2 in items["$values"]:
                        if "SwitchFilter" in item2["$type"]:
                            self.filter = item2["Filter"]
                        if "TakeExposure"  in item2["$type"]:
                            self.exposure = item2

                if "TimeCondition" in item["$type"]:
                    self.timecondition = item

                if "ExternalScript" in item["$type"]:
                    self.script_w_target = item

                if "RunAutofocus" in item["$type"]:
                    self.autofocus = item


    def add_parent(self, id):
        # verbose("NINATarget(add_parent):", "id =", id)
        self.obj["Parent"] = { "$ref": id }


    def set_expanded(self, flag):
        # verbose("NINATarget(set_expanded):", "flag =", flag)
        self.obj["IsExpanded"] = flag




class NINASequence(NINABase):
    """ Holds data read from N.I.N.A JSON sequence """

    def __init__(self):
        self.name         = None # template name
        # instance variables created by process_data(), referencing the data objects, used by update_target_data()
        self.start_list   = None
        self.targets_list = None
        self.end_list     = None
        self.start_id     = None
        self.targets_id   = None
        self.end_id       = None

        # This dictionary holds the ids of the various time SelectedProvider{}, to be replaced with reference
        # on further occurences
        self.provider_dict = {}

        NINABase.__init__(self)


    def process_data(self):
        self.name = self.obj["Name"]
        verbose("NINASequence(process_data):", "name =", self.name)

        items = self.obj["Items"]["$values"]
        self.area_list = items
        for item in items:
            for k in item.keys():
                if k=="$id" or k=="$type" or k=="$ref" or k=="Name":
                    self.print_attr(item, k, ">")

        self.start_list   = self.area_list[0]["Items"]["$values"]
        self.targets_list = self.area_list[1]["Items"]["$values"]
        self.end_list     = self.area_list[2]["Items"]["$values"]
        self.start_id     = self.area_list[0]["$id"]
        self.targets_id   = self.area_list[1]["$id"]
        self.end_id       = self.area_list[2]["$id"]
        # ic(self.start_id, self.start_list)
        # ic(self.targets_id, self.targets_list)
        # ic(self.end_id, self.end_list)


    def search_container(self, container):
        """ Search for named container in target area (targets_list) """
        if container:
            verbose("NINASequence(search_container):", "name =", container)
            for item in self.targets_list:
                type = item["$type"]
                name = item["Name"]
                if "Container.SequentialContainer" in type and name==container:
                    self.targets_list = item["Items"]["$values"]
                    self.targets_id   = item["$id"]
                    if ic.enabled:
                        ic("NEW in container")
                        ic(self.targets_id, self.targets_list)                    
                    break
                


    def append_target(self, target):
        # verbose("NINASequence(append_target):", "name =", target.name)

        self.targets_list.append(target.obj)


    def process_csv(self, target_tmpl: NINATarget, 
                    file: str, target_format: str, tzname: str, 
                    subdir: str=None, autofocus1: bool=False):
        tz_local = ZoneInfo(tzname)

        with open(file, newline='') as f:
            reader = csv.DictReader(f)
            first_target = True
            for count, row1 in enumerate(reader):
                # Make field names lower case
                row = { key.lower():val for key, val in row1.items() if key}
                # Number in sequence, default 0
                seq = int(row.get("#") or count + 1)

                # Target name, must not contain [/:"]
                target = row.get("object") or row.get("name") or row.get("target")
                if not target:
                    error("can't find target name in CSV data")
                target = target.replace("/", "").replace(":", "").replace("\"", "")

                # Date / time UTC and local
                ## FIXME: should be way more flexible
                time_utc = time_local = None
                # Try ISO format first
                obstime = row.get("obstime") or row.get("observation time")
                if obstime:
                    time_utc = datetime.fromisoformat(obstime)
                    time_local = time_utc.astimezone(tz_local)
                # NEO Planner CSV output
                else:
                    date1 = row.get("observation date")
                    time1 = row.get("time ut")
                    if date1 and time1:
                        time_utc = datetime.fromisoformat(date1.replace(" ", "-") + "T" + time1 + ":00+00:00")
                        # Python 3.9 doesn't like the "Z" timezone declaration, thus +00:00
                        # convert to local time zone (now configurable)
                        time_local = time_utc.astimezone(tz_local)

                # Use various field names for RA/DEC coordinates in CSV data
                ra  = row.get("ram")  or row.get("ra")
                dec = row.get("decm") or row.get("dec.") or row.get("dec") or row.get("de")
                # End marker
                if target=="Azelfafage" or target=="0" or not target or not ra or not dec:
                    break
                try:
                    coord = Coord(ra, dec)
                except ValueError:
                    error(f"invalid coordinates {ra=} {dec=}")

                # Default 30s exposure time
                exp = float(row.get("exposure time") or row.get("exposure") or row.get("exp") or 30)
                # Default 5 frames
                number = int(row.get("no images") or row.get("number") or 5)

                # Default filter L
                filter = row.get("filter") or "L"
                if filter:
                    for fn in DEFAULT_FILTER_NAMES:
                        if filter.startswith(fn):
                            filter = fn
                            break

                # Format target name for templates:
                # 0=target, 1=date, 2=seq, 3=number
                # (from config)
                # date1 = time_local.date() if time_local else ""
                # Use today or --date
                date1 = Options.date
                formatted_target = target_format.format(target, date1, seq, number)
                ic(target, formatted_target)

                # Replace target with formatted target, as NINA currently only supports $$TARGETNAME$$
                # in filename templates under Options > Imaging
                target = formatted_target

                # Just output the target list
                if Options.list_targets:
                    print(target)
                    continue

                message("------------------------------------------------------------------")
                message(f"#{seq:03d} {target:32s} {coord}")

                if time_utc:
                    message(f"UT={time_utc} / local {time_local}")
                message(f"{number:d}x{exp:.1f}s filter={filter}")

                # AF only for first target if option "autofocus_first_target_only" is set
                no_af = autofocus1 and not first_target
                # default for binning
                data = TargetData(formatted_target, target, coord, time_local, number, exp, filter, 
                                  subdir=subdir, no_autofocus=no_af)

                # create deep copy of target object, update with data read from CSV
                target_new = copy.deepcopy(target_tmpl)
                target_new.update_target_data(data)

                ### append to main sequence targets
                # the following changes are necessary to append target_new to the list of the target area
                target_new.set_prefix(seq)
                target_new.traverse(NINABase.add_prefix_to_id)
                # add parent ref to target object
                target_new.add_parent(self.targets_id)
                # collapse view
                target_new.set_expanded(False)
                self.append_target(target_new)

                first_target = False

            if not Options.list_targets:
                message("------------------------------------------------------------------")

        # update all SelectedProvider {...} with references
        self.traverse(NINABase.process_provider, self.provider_dict)




def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = "Create/populate multiple N.I.N.A target templates/complete sequence with data from NEO Planner CSV",
        epilog      = "Version: " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="debug messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")
    arg.add_argument("-A", "--debug-print-attr", action="store_true", help="extra debug output")
    arg.add_argument("-D", "--destination-dir", help="output dir for created sequence")
    arg.add_argument("-o", "--output", help="output .json file")
    arg.add_argument("-n", "--no-output", action="store_true", help="dry run, don't create output files")
    arg.add_argument("-l", "--list-targets", action="store_true", help="list targets only")
    arg.add_argument("-S", "--setting", help="use template/target SETTING from config")
    arg.add_argument("--date", help=f"use DATE for generating sequence (default {Options.date})")
    arg.add_argument("filename", nargs="+", help="CSV target data list")
   
    args = arg.parse_args()

    if args.debug:
        ic.enable()
        ic(sys.version_info, sys.path)
    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()

    Options.debug_print_attr = args.debug_print_attr
    Options.list_targets = args.list_targets

    if args.date:
        Options.date = args.date

    if args.setting:
        if not args.setting in config.get_keys():
            error(f"setting {args.setting} not in config")
        setting = config.get(args.setting)
    else:
        error(f"must supply a setting with --setting, valid:",
              ", ".join([k for k in config.get_keys() if not k.startswith("#")]))
    ic(args.setting, setting)

    target_template = setting["target"]
    verbose("processing target template", target_template)
    sequence_template = setting["template"]
    verbose("processing sequence template", sequence_template)
    target_format = setting["format"]
    verbose("target format (0=target, 1=date, 2=seq, 3=number)", target_format)
    output_format = setting["output"]
    verbose("output format (1=date)", output_format)
    container = setting["container"]
    verbose(f"add target items to container '{container}', empty=target area")
    tzname = setting["timezone"]
    verbose("timezone", tzname)
    subdir = setting.get("subdir")
    verbose("subdir (1=date)", subdir)
    if subdir:
        subdir = subdir.format("", Options.date)
    autofocus_first_target_only = setting.get("autofocus_first_target_only") or False
    verbose(f"autofocus first target only = {autofocus_first_target_only}")

    if args.destination_dir:
        destination_dir = args.destination_dir
    else:
        destination_dir = os.path.join(config.get_documents_path(), "N.I.N.A")
    verbose("destination directory", destination_dir)
    if args.output:
        output = args.output
    else:
        output = output_format.format("", Options.date)
    verbose("output file", output)

    target = NINATarget()
    target.read_json(target_template)
    target.process_data()

    sequence = NINASequence()
    sequence.read_json(sequence_template)
    sequence.process_data()
    sequence.search_container(container)

    for f in args.filename:
        verbose("processing CSV file", f)
        sequence.process_csv(target, f, target_format, tzname, subdir, autofocus_first_target_only)

    output_path = os.path.join(destination_dir, output)
    if not args.no_output and not args.list_targets:
        verbose("writing JSON sequence", output_path)
        sequence.write_json(output_path)


   
if __name__ == "__main__":
   main()
