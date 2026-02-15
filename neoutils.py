#!/usr/bin/env python

# Copyright 2025-2026 Martin Junius
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
# Version 0.1 / 2025-11-25
#       New module with NEO utility functions
# Version 0.2 / 2026-02-04
#       Ephemeris column names can be passed as optional parameters
#
# Usage:
#       from neoutils import ...

VERSION     = "0.2 / 2026-02-04"
AUTHOR      = "Martin Junius"
NAME        = "neoutils"
DESCRIPTION = "NEO utility functions"

import sys
import csv
from itertools import pairwise

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# AstroPy & friends
import astropy.units as u
from astropy.coordinates import Angle
from astropy.units import Quantity, Magnitude
from astropy.time import Time
from astropy.table import QTable, Row
import numpy as np
from sbpy.data import Ephem

# Local modules
from verbose import verbose, warning, error, message
from neoconfig import config
from neoclasses import Exposure, EphemData, EphemTimes



def single_exp(motion: Quantity) -> Quantity:
    """
    Compute exposure time depending on motion value,

    Parameters
    ----------
    motion : Quantity
        Motion value as arcsec/min

    Returns
    -------
    Quantity
        Exposure time as secs or None, if too fast
    """
    exp_times = config.exposure_times
    exp = config.pixel_tolerance * config.resolution * u.arcsec / motion
    #            ^ pixels                 ^ arcsec/pixel          ^ arcsec/min
    exp = exp.to(u.s).value
    exp_min = exp_times[0]
    if exp < 0.9 * exp_min: # allow a bit of tolerance
        return None
    for exp1 in exp_times:
        if exp1 > exp:
            break
        exp_min = exp1
    return exp_min * u.s



def motion_limit() -> Quantity:
    """
    Get motion limit

    Returns
    -------
    Quantity
        Motion limit derived from minimum exposure time
    """
    exp_times = config.exposure_times
    return config.pixel_tolerance * config.resolution * u.arcsec / (exp_times[0] * u.s).to(u.min)



def exposure_calc(max_motion: Quantity, mag: Magnitude) -> Exposure:
    """
    Calculate number of exposures, single exposure time, total exposure time, total time, percentage
    from object motion and magnitude

    Parameters
    ----------
    max_motion : Quantity
        Object motion
    mag : Magnitude
        Object magnitude

    Returns
    -------
    Exposure
        Exposure object
    """
    exp   = single_exp(max_motion)                      # Single exposure / s
    if exp == None:                                     # Object too fast
        return None

    min_n_exp = config.min_n_exp
    max_n_exp = config.max_n_exp
    base_mag  = config.base_mag
    base_exp  = config.base_exp

    rel_brightness = 10 ** (0.4 * (mag.value - base_mag))
    total_exp = base_exp * u.s * rel_brightness         # Total exposure
    n_exp = int(total_exp / exp) + 1                    # Number of exposures
    ic(base_mag, base_exp, mag.value, rel_brightness, total_exp, n_exp)
    perc_of_required = 100.                             # Percentage actual / total exposure
    if n_exp < min_n_exp:
        perc_of_required = min_n_exp / n_exp * 100
        n_exp = min_n_exp
    if n_exp > max_n_exp:
        perc_of_required = max_n_exp / n_exp * 100
        n_exp = max_n_exp
    total_exp = (n_exp * exp).to(u.min)
    total_time = ( total_exp 
                   + config.dead_time_slew_center * u.s 
                   + config.dead_time_af * u.s
                   + config.dead_time_guiding * u.s  
                   + config.safety_margin * u.s
                   + n_exp * config.dead_time_image * u.s )
    ic(n_exp, exp, total_exp, total_time, perc_of_required)

    return Exposure(n_exp, exp, total_exp, total_time, perc_of_required)



def get_mag0(eph: Ephem, column: str="Mag") -> Magnitude:
    """Get magnitude from ephemeris (first row)

    Parameters
    ----------
    eph : Ephem
        ephemeris
    column : str, optional
        column name, by default "Mag"

    Returns
    -------
    Magnitude
        Magnitude quantity
    """
    return eph[column][0]



def get_table(eph: Ephem) -> QTable:
    """Hack for proper iteration: get actual table from Ephem object

    Parameters
    ----------
    eph : Ephem
        Ephemeris object

    Returns
    -------
    QTable
        Table object
    """
    # Hack for iterating over Ephem object:
    # for row in eph yields single line QTable objects
    # for row in get_table(eph) yields Row object as expected
    # Reported as a bug to sbpy, may change with future implementations
    return eph._table



def max_motion(eph: Ephem, column: str="Motion") -> Quantity:
    """
    Get max value for motion column(s) from ephemeris table

    Parameters
    ----------
    eph : QTable
        Ephemeris table
    column : str, optional
        Motion column name(s), by default "motion", comma separated for RA/DEC motion

    Returns
    -------
    Quantity
        Max motion of object
    """
    max_m = -1 * u.arcsec / u.min
    if "," in column:
        # Separate RA*cos(DEC), DEC motion columns
        col1, col2 = column.split(",")
    else:
        # Single column with proper motion
        col1 = column
        col2 = None
    for row in get_table(eph):
        if col2:
            motion = np.sqrt( np.square(row[col1]) + np.square(row[col2]) )
        else:
            motion = row[col1]
        if motion > max_m:
            max_m = motion

    return max_m.to(u.arcsec / u.min)



def is_east(az: Angle) -> bool:
    """
    Test for east azimut position

    Parameters
    ----------
    az : Angle
        Azimut angle

    Returns
    -------
    bool
        True if east, False if west
    """
    az180     = 180 * u.degree
    return True if az >= 0 and az < az180 else False

def is_west(az: Angle) -> bool:
    """
    Test for west azimut position

    Parameters
    ----------
    az : Angle
        Azimut angle

    Returns
    -------
    bool
        True if west, False if east
    """
    az180     = 180 * u.degree
    az360     = 360 * u.degree
    return True if az >= az180 and az < az360 else False

def flip_times(eph: Ephem, col_obstime: str="Obstime", col_az: str="Az") -> tuple[Time, Time]:
    """
    Get time before and after meridian passing from ephemeris table

    Parameters
    ----------
    eph : Ephem
        Ephemeris table
    col_obstime : str
        Name of obstime column
    col_az : str
        Name of azimuth column

    Returns
    -------
    Tuple[Time, Time]
        Time before meridian, time after meridian passing
        None, None if object doesn't pass meridian
    """
    prev_time = None
    prev_az   = None

    for row in get_table(eph):
        time = row[col_obstime]
        az   = row[col_az]
        # ic(time, az)
        if not prev_az == None:
            if is_east(prev_az) and is_west(az):     # South flip
                return (prev_time, time)
            if is_west(prev_az) and is_east(az):     # North flip
                return (prev_time, time)
        prev_time = time
        prev_az   = az

    # No meridian passing found
    return None, None



def opt_alt_times(eph: Ephem, alt: Angle, col_obstime: str="Obstime", col_alt: str="Alt") -> tuple[Time, Time]:
    """
    Get times for object above specified altitude

    Parameters
    ----------
    eph : Ephem
        Ephemeris table
    alt : Angle
        Altitude angle
    col_obstime : str
        Name of obstime column
    col_alt : str
        Name of altitude column

    Returns
    -------
    tuple[Time, Time]
        Start time above altitude, end time above altitude
    """
    time_alt0 = None
    time_alt1 = None

    for row in get_table(eph):
        if time_alt0 == None and row[col_alt] >= alt:
            time_alt0 = row[col_obstime]
        if time_alt0 != None and row[col_alt] >= alt:
            time_alt1 = row[col_obstime]
        if time_alt1 != None and row[col_alt] < alt:
            break
    
    return time_alt0, time_alt1



def max_alt_time(eph: Ephem, col_obstime: str="Obstime", col_alt: str="Alt") -> Time:
    """
    Get time of maximum altitude from ephemeris table

    Parameters
    ----------
    eph : Ephem
        Ephemeris table
    col_obstime : str
        Name of obstime column
    col_alt : str
        Name of altitude column

    Returns
    -------
    Time
        Time of max altitude
    """
    max_alt = -90 * u.degree
    time_max = None
    for row in get_table(eph):
        if row[col_alt] > max_alt:
            max_alt = row[col_alt]
            time_max = row[col_obstime]
    return time_max



def get_row_for_time(eph: Ephem, t: Time, col_obstime: str="Obstime") -> Row:
    """
    Get row from ephemeris table closest to specified time

    Parameters
    ----------
    eph : Ephem
        Ephemeris table
    t : Time
        Time for retrieving row
    col_obstime : str
        Name of obstime column

    Returns
    -------
    Row
        Corresponding row from ephemeris table
        None, if not found
    """
    for r1, r2 in pairwise(get_table(eph)):
        if r1[col_obstime] <= t and t <= r2[col_obstime]:
            return r1
    
    # No matching interval found
    return None



def ephem_data_add_times(edata: EphemData, col_obstime: str="Obstime", col_alt: str="Alt", col_az: str="Az", use_old_sort: bool=False) -> None:
    """Add EphemTimes data to EphemData object

    Parameters
    ----------
    edata : EphemData
        EphemData object
    col_obstime : str, optional
        observation time column name, by default "Obstime"
    col_alt : str, optional
        Altitude column name, by default "Alt"
    col_az : str, optional
        Azimuth column name, by default "Az"
    use_old_sort : bool, optional
        Use old sort order (neocp.py), by default False

    Raises
    ------
    ValueError
        Empty ephemeris
    """
    eph = edata.ephem
    # Ephemeris could be empty
    if len(eph) == 0:
        raise ValueError(f"ephemeris for object {edata.obj} is empty")

    edata.times = EphemTimes(eph[col_obstime][0], eph[col_obstime][-1],
                             None, None, None, None, None, None)
    etimes = edata.times
    etimes.before, etimes.after = flip_times(eph, col_obstime, col_az)
    etimes.alt_start, etimes.alt_end = opt_alt_times(eph, config.opt_alt * u.deg, col_obstime, col_alt)

    etimes.max_alt = max_alt_time(eph, col_obstime, col_alt)
    edata.sort_time = etimes.max_alt
    if use_old_sort:
        # Old sort order of neocp.py
        if etimes.before:
            edata.sort_time = etimes.before
    else:
        # New sort order of neo-obs-planner
        if etimes.alt_start != None:
            edata.sort_time = etimes.alt_start



def obj_data_add_times(obj_data: dict[str, EphemData], col_obstime: str="Obstime", col_alt: str="Alt", col_az: str="Az", use_old_sort: bool=False) -> dict[str, EphemData]:
    """Fill EphemTimes for all objects

    Args:
        obj_data (dict[str, EphemData]): objects dict

    Returns:
        dict[str, EphemData]: processed objects dict
    """
    for obj in obj_data.keys():
        ephem_data_add_times(obj_data[obj], col_obstime, col_alt, col_az, use_old_sort)
    return obj_data



def sort_obj_data(obj_data: dict[str, EphemData]) -> dict[str, EphemData]:
    """Sort objects dict by sort_time attribute

    Args:
        obj_data (dict[str, EphemData]): objects dict

    Returns:
        dict[str, EphemData]: sorted objects dict
    """
    return { obj: edata for obj, edata in sorted(obj_data.items(), key=lambda item: item[1].sort_time) }



def verbose_obj_data(obj_data: dict[str, EphemData]) -> None:
    """Verbose output of object data dict

    Parameters
    ----------
    obj_data : dict[str, EphemData]
        Object data dict
    """
    for obj in obj_data:
        verbose("===================================================================================================================")
        verbose(f"{obj} ephemeris")
        verbose.print_lines(obj_data[obj].ephem)
    verbose("===================================================================================================================")



TIME_FORMAT = "%Y-%m-%d %H:%M:%S"
TIME_FORMAT_TZ = "%Y-%m-%d %H:%M:%S+0000"

def fmt_time(time: Time|None, add_tz: bool=False) -> str:
    """Format time

    Parameters
    ----------
    time : Time | None
        Time value
    add_tz : bool, optional
        Add timezone "+0000", by default False

    Returns
    -------
    str
        Formatted time
    """
    return time.strftime(TIME_FORMAT_TZ if add_tz else TIME_FORMAT) if time != None else "-" * 19



def obj_data_csv_output(obj_data: dict[str, EphemData], output: str) -> None:
    csv_rows = list()

    # Traverse objects, only those with valid plan_start time
    for obj, edata in obj_data.items():
        if edata.times.plan_start:
            # CSV output:
            #   start time, end time, 
            #   target=id, observation date (YYYY-MM-DD), time ut (HH:MM), ra, dec, exposure, number, filter (L),
            #   type, mag, nobs, arc, notseen, total
            # RA output  = hourangle !!!
            # DEC output = degree
            csv_row = { "target": obj,
                        "obstime": fmt_time(edata.times.plan_start, add_tz=True),
                        "ra": float(edata.ra.value),
                        "dec": float(edata.dec.value),
                        "exposure": float(edata.exposure.single.value),
                        "number": edata.exposure.number,
                        "filter": "L",
                        "start time": fmt_time(edata.times.plan_start),
                        "end time": fmt_time(edata.times.plan_end),
                        "type": edata.neocp.type if edata.neocp else "",
                        "mag": float(edata.mag.value),
                        "nobs": int(edata.neocp.nobs) if edata.neocp else "",
                        "arc": float(edata.neocp.arc.value) if edata.neocp else "",
                        "notseen": float(edata.neocp.notseen.value) if edata.neocp else "",
                        "total": str(edata.exposure)
                        }
            ic(csv_row)
            csv_rows.append(csv_row)

    # Output to CSV file for nina-create-sequence2
    verbose(f"planned objects for nina-create-sequence2: {output}")
    # csv_row is the last object, if any were found
    if csv_row:
        ##FIXME: improve csvoutput module to cover this usage
        fieldnames = csv_row.keys()
        ic(fieldnames)
        if output:
            with open(output, "w", newline="") as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(csv_rows)
        else:
            writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(csv_rows)
    else:
        warning("no objects, no CSV output")
