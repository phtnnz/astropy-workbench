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
# Version 0.0 / 2026-02-03
#       plot functions moved to this module

VERSION = "0.0 / 2026-02-03"
AUTHOR  = "Martin Junius"
NAME    = "neoplot"

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.coordinates import AltAz, EarthLocation
import astropy.units as u
from astropy.time import Time
import numpy as np
from astropy.table import QTable

# Astroplan
import matplotlib.pyplot as plt
from astroplan import Observer
from astroplan.plots import plot_altitude, plot_sky

# Local modules



def qtable_to_altaz(id: str, qt: QTable, loc: EarthLocation) -> AltAz:
    """
    Convert table rows "alt"/"az" to AltAz object

    Parameters
    ----------
    id : str
        NEOCP id = temporary designation
    qt : QTable
        Table with ephemeris, including alt/az

    Returns
    -------
    AltAz
        AltAz coordinates object for altitude/sky plot
    """
    altaz = AltAz(alt=qt["alt"], az=qt["az"], obstime=qt["obstime"], location=loc)
    # Quick hack to get a proper label for plot_altitude()
    altaz.name = id
    ic(altaz)
    return altaz



def plot_objects(ephemerides: dict, objects: list, filename: str, loc: EarthLocation) -> None:
    """
    Generate altitude and sky plot

    Parameters
    ----------
    ephemerides : dict
        Dictionary with ephemerides for all objects
    objects : list
        Subset list of objects for plot
    filename : str
        File name for generated PNG
    """
    # Get next midnight
    time = Time(Time.now(), location=loc)
    observer = Observer(location=loc, description=loc.info.name)
    midnight = observer.midnight(Time.now(), which="next")
    ic(midnight)

    # Intervals around midnight
    time_interval = midnight + np.linspace(-8, 8, 160)*u.hour
    moon          = observer.moon_altaz(time_interval)
    # Quick hack to get a proper label for plot_altitude()
    moon.name     = "Moon"

    # Subplots
    fig = plt.figure(figsize=(15, 6))
    ax1 = plt.subplot(1, 2, 1)

    # Plot altitude for all NEOCP objects
    for id in objects:
        qt = ephemerides[id]
        altaz = qtable_to_altaz(id, qt)
        plot_altitude(altaz, observer, altaz.obstime, ax1, style_kwargs=dict(fmt="o"))
    plot_altitude(moon, observer, moon.obstime, ax1, brightness_shading=True, style_kwargs=dict(fmt="y--"))
    plt.legend(bbox_to_anchor=(1.0, 1.015))

    # Hourly intervals around midnight
    time_interval = midnight + np.linspace(-5, 5, 11)*u.hour
    moon          = observer.moon_altaz(time_interval)
    # Quick hack to get a proper label for plot_sky()
    moon.name     = "Moon"

    # Subplot for altitude
    ax2 = plt.subplot(1, 2, 2, projection='polar')

    # Plot sky for all NEOCP objects
    for id in objects:
        qt = ephemerides[id]
        altaz = qtable_to_altaz(id, qt)
        plot_sky(altaz, observer, altaz.obstime, ax2)
    plot_sky(moon, observer, moon.obstime, ax2, style_kwargs=dict(color="y", marker="x"))
    # plt.legend(bbox_to_anchor=(1.32, 1.15))

    plt.subplots_adjust(wspace=0.3)
    plt.savefig(filename, bbox_inches="tight")
    plt.close()
