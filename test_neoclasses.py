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
# Version 0.0 / 2026-05-16
#       Test EphemDataList

# Tests for neoclasses


# EphemDataList
from neoclasses import EphemData, EphemDataList
test_obj_edata = {
    "1999 A": EphemData("TEST", "1999 A"),
    "2000 B": EphemData("TEST", "2000 B"),
    "2001 C": EphemData("TEST", "2001 C"),
    "(1234)": EphemData("TEST", "(1234)")
}

test_objects = [ "1999 A", "2000 B", "2001 C", "(1234)"]

test_list1 = EphemDataList.from_dict(test_obj_edata)
print(test_list1)
test_list2 = EphemDataList.from_objects(test_objects)
print(test_list2)
