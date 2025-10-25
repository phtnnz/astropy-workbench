#!/usr/bin/env python

# Copyright 2025 Martin Junius
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
# Version 0.1 / 2025-10-25
#       Plot field rotation over alt/az

VERSION     = "0.1 / 2025-10-25"
AUTHOR      = "Martin Junius"
NAME        = "field-rot"
DESCRIPTION = "Plot field rotation"

import sys
import argparse
from typing import Tuple, Any

# The following libs must be installed with pip
from icecream import ic
# Disable debugging
ic.disable()

# AstroPy
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.coordinates import Angle
from astropy.coordinates import errors
from astropy.time        import Time, TimeDelta
import astropy.units as u
import astropy.constants as const
from astropy.units import Quantity
import numpy as np

from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

# Local modules
from verbose import verbose, warning, error, message
from astroutils import location_to_string, get_location     # ...

DEFAULT_LOCATION = "M49"



# Command line options
class Options:
    loc: EarthLocation = None       # -l --location



def field_rot(loc: EarthLocation, alt: Angle, az: Angle) -> Quantity:
    """
    Compute field rotation
    See also https://de.wikipedia.org/wiki/Bildfelddrehung

    Parameters
    ----------
    loc : EarthLocation
        Observer location
    alt : Angle
        Altitude angle
    az : Angle
        Azimuth angle

    Returns
    -------
    Quantity
        Field rotation
    """
    ic(loc, alt, az)
    sday  = 1 * u.sday
    omega = 360 * u.deg / sday.to(u.hour)
    rot   = omega * np.cos(loc.lat) * np.cos(az) / np.cos(alt)
    ic(omega, loc.lat, rot)
    return rot



def plot(loc: EarthLocation) -> None:
    """
    3D plot of field rotation over azimuth and altitude for the given location

    Parameters
    ----------
    loc : EarthLocation
        Observer location
    """
    fig = plt.figure(figsize=(15, 10))
    ax  = fig.add_subplot(projection='3d')

    az  = np.linspace(0, 360, 72+1) * u.deg
    alt = np.linspace(0, 85,  17+1) * u.deg
    ic(az, alt)

    az, alt = np.meshgrid(az, alt)
    rot = field_rot(loc, alt, az)
    ic(rot)

    ax.set_xticks([0, 60, 120, 180, 240, 300, 360])
    surf = ax.plot_surface(az, alt, rot, rstride=1, cstride=1, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)
    fig.colorbar(surf, shrink=0.4, aspect=10, pad=0.05)
    ax.set_title(f"Field Rotation Rate for {loc.info.name}", fontsize=20)
    ax.set_xlabel("Azimuth deg", fontsize=16)
    ax.set_ylabel("Altitude deg", fontsize=16)
    ax.set_zlabel("Field rotation deg / h", fontsize=16)

    plt.savefig("tmp/plot.png", bbox_inches="tight")
    plt.close()



def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = DESCRIPTION,
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")

    arg.add_argument("-l", "--location", help=f"coordinates, named location or MPC station code, default {DEFAULT_LOCATION}")

    args = arg.parse_args()

    if args.debug:
        ic.enable()
        ic(sys.version_info, sys.path, args)
    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()

    # Observer location
    loc = get_location(args.location if args.location else DEFAULT_LOCATION)
    ic(loc, loc.to_geodetic())
    Options.loc = loc
    verbose(f"location {location_to_string(loc)}")

    # Tests
    field_rot(loc, 45 * u.deg,   0 * u.deg)
    field_rot(loc, 45 * u.deg,  90 * u.deg)
    field_rot(loc, 60 * u.deg, 180 * u.deg)
    field_rot(loc, 60 * u.deg, 270 * u.deg)

    # 3D plot
    plot(loc)



if __name__ == "__main__":
    main()
