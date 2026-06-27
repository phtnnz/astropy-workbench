#!/usr/bin/env python

# Copyright 2024-2026 Martin Junius
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
# Version 0.1 / 2025-01-27
#       Utility functions moved to this module
# Version 0.2 / 2025-10-18
#       Added time_jd_as_iso(), added address search via of_address() to get_location()
# Version 1.0 / 2026-06-16
#       Moved and adapted to new directory structure under neoop/

from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.coordinates import Angle
import astropy.units as u
from astropy.time        import Time
import numpy as np

# Local modules
from utils.verbose import error



VERSION = "1.0 / 2026-06-16"
AUTHOR  = "Martin Junius"
NAME    = "astroutils"



def ra_from_lst_ha(lst: Angle, ha: Angle) -> Angle:
    """
    Compute RA from local sidereal time and hour angle

    :param lst: local sidereal time
    :type lst: Angle
    :param ha: hour angle
    :type ha: Angle
    :return: RA
    :rtype: Angle
    """
    ra = lst - ha
    ra.wrap_at(24*u.hourangle, inplace=True)
    ic(lst, ha, ra)
    return ra



def coord_to_string(coord: SkyCoord, short: bool=False) -> str:
    """
    Format SkyCoord as string

    :param coord: coordinates object
    :type coord: SkyCoord
    :return: _description_
    :rtype: str
    """
    return ra_dec_to_string(coord.ra, coord.dec, short)

def ra_dec_to_string(ra: Angle, dec: Angle, short: bool=False) -> str:
    """
    Format RA and DEC as string

    :param ra: RA
    :type ra: Angle
    :param dec: DEC
    :type dec: Angle
    :return: formatted string "RA=... DEC=..."
    :rtype: str
    """
    if short:
        return f"{ra.to_string(unit=u.hour, precision=0)} {dec.to_string(unit=u.degree, precision=0)}"
    else:
        return f"RA={ra.to_string(unit=u.hour, precision=2)} DEC={dec.to_string(unit=u.degree, precision=2)}"

def angle_to_string(a: Angle, decimal: bool=False) -> str:
    """
    Format angle as string

    :param a: angle
    :type a: Angle
    :param decimal: flag for decimal output, default False
    :type decimal: bool
    :return: formatted string angle in degrees, minutes, seconds
    :rtype: str
    """
    return f"{a.to_string(unit=u.degree, decimal=decimal, precision=2)}"


def altaz_to_string(alt: Angle, az: Angle, decimal: bool=True) -> str:
    """
    Format alt/az as string

    Parameters
    ----------
    alt : Angle
        altitude
    az : Angle
        azimuth
    decimal : bool, optional
        decimal output, by default True

    Returns
    -------
    str
        string representation of alt/az
    """
    return f"alt={angle_to_string(alt, decimal)} az={angle_to_string(az, decimal)}"


def hourangle_to_string(a: Angle) -> str:
    """
    Format hour angle as string

    :param a: hour angle
    :type a: Angle
    :return: formatted string hour angle in hours, minutes, seconds
    :rtype: str
    """
    return f"{a.to_string(unit=u.hour, precision=2)}"



def location_to_string(loc: EarthLocation) -> str:
    return f"lon={loc.to_geodetic().lon.to_string(unit=u.degree, precision=1)} lat={loc.to_geodetic().lat.to_string(unit=u.degree, precision=1)} height={loc.to_geodetic().height.to_string(precision=0)}"



def time_jd_as_iso(jd: np.float64) -> Time:
    """
    Create Time object from JD value, set format to "iso"

    Parameters
    ----------
    jd : np.float64
        JD value

    Returns
    -------
    Time
        Astropy Time object
    """
    time = Time(jd, format="jd")
    time.format = "iso"
    return time
