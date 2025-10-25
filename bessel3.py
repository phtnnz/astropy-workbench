#!/usr/bin/env python

# Copyright 2025 Martin Junius, Uwe Pilz (VdS)
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
# Version 0.0 / 2025-10-18
#       Included Uwe Pilz' code mostly verbatim
# Version 0.1 / 2025-10-19
#       Major overhaul to fit into the astropy-workbench environment,
#       using astropy and numpy, a lot of renaming for better understanding
# Version 0.2 / 2025-10-21
#       Further rework, output clean-up, new option -0 --positive-mag-only
# Version 0.3 / 2025-10-23
#       Further rework and clean-up, variables renamed in accordance with [ESAA]
# Version 0.4 / 2025-10-23
#       Added sun altitude using astropy get_sun(),
#       new option -T --totality for listing
#
# See [ESAA] Explanatory Supplement to the Astronomical Almanac, 3rd Edtion
# Chapter 11 - Eclipses of the Sun and Moon


VERSION     = "0.4 / 2024-10-23"
AUTHOR      = "Martin Junius"
NAME        = "bessel3"
DESCRIPTION = "Solar eclipse local circumstances"

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
import astropy.units as unit
import astropy.constants as const
import numpy as np
from numpy.polynomial import Polynomial

from astropy.coordinates import solar_system_ephemeris
from astropy.coordinates import get_sun

# Local modules
from verbose import verbose, warning, error, message
from astroutils import location_to_string, get_location     # ...

DEFAULT_LOCATION = "Burgos, Spain"



# Earth equatorial radius
R_earth = 6378137.0 * unit.m        # GRS 80/WGS 84 value (Wikipedia)
                                    # https://en.wikipedia.org/wiki/World_Geodetic_System
R_earth_pol = 6356752.0 * unit.m    # https://en.wikipedia.org/wiki/Earth
f_earth = (R_earth - R_earth_pol) / R_earth
# Moon equatorial radius
R_moon = 1738100.0 * unit.m
R_moon = 0.2725076 * R_earth        # IAU 1982
R_moon = 0.272281  * R_earth        # smaller value from https://eclipse.gsfc.nasa.gov/SEpubs/20080801/TP214149b.pdf



# Command line options
class Options:
    loc: EarthLocation = None       # -l --location
    list: bool = False              # -L --list
    pos_mag: bool = False           # -0 --positive-mag-only
    totality: bool = False          # -T --totality



def test_tse2026() -> Tuple[EarthLocation, Time]:
    """
    Get location and time for test case: solar eclipse 12 Aug 2026, Burgos, Spain

    Returns
    -------
    Tuple[EarthLocation, Time]
        Location and time of max eclipse
    """
    # Test case TSE 2026 @ Burgos, Spain
    #
    # 42° 21' 01.68" N	  <—>  	42.35047°
    # 3° 41' 21.68" W	  <—>  	-3.68935°
    # 891.0m
    #
    # Maximum eclipse (MAX) : 2026/08/12 18:29:17.6
    loc = EarthLocation(lon=-3.68935*unit.degree, lat=42.35047*unit.degree, height=891*unit.m)
    time = Time("2026-08-12 18:29:17.6", scale="ut1", location=loc)
    ic(loc, time, time.jd)
    return loc, time



# A major part of the following code ist originally by Uwe Pilz
# VdS-Journal 95: Lokaler Verlauf einer Finsternis
# https://fg-astrophysik.vdsastro.de/prg95
# (no license provided)
# Heavily modified and using numpy, astropy
#
# sofiBessel3 : lokaler Verlauf einer Sonnenfinsternis
# ohne Berücksichtigung der lokalen Sonnenhöhe

# Besselian elements for TSE 12 Aug 2026 from
# https://eclipse.gsfc.nasa.gov/SEsearch/SEdata.php?Ecl=20260812
#
# x, y   - Cartesian coordinates of the lunar shadow axis in the Fundamental Plane 
#          (in units of Earth's equatorial radius)
# L1, L2 - Radii of the Moon's penumbral and umbral/antumbral shadows in the Fundamental Plane 
#          (in units of Earth's equatorial radius)
# d      - Declination of the Moon's shadow axis on the celestial sphere
# µ      - Hour angle of the Moon's shadow axis on the celestial sphere
# f1, f2 - Angles of the penumbral and umbral/antumbral shadow cones with respect to the axis of the lunar shadow

bessel_T0    = 18
time_T0      = Time("2026-08-12 18:00:00", scale="tt")

bessel_x     = Polynomial([  0.4755140,  0.5189249, -0.0000773, -0.0000080 ])
bessel_y     = Polynomial([  0.7711830, -0.2301680, -0.0001246,  0.0000038 ])
bessel_d     = Polynomial([ 14.7966700, -0.0120650, -0.0000030,  0.0       ])
bessel_l1    = Polynomial([  0.5379550,  0.0000939, -0.0000121,  0.0       ])
bessel_l2    = Polynomial([ -0.0081420,  0.0000935, -0.0000121,  0.0       ])
bessel_mu    = Polynomial([ 88.7477900, 15.0030900,  0.0000000,  0.0       ])

bessel_tanf1 = 0.0046141
bessel_tanf2 = 0.0045911

# all variables with _p ("prime") suffix are derivatives
bessel_x_p   = bessel_x.deriv()
bessel_y_p   = bessel_y.deriv()
bessel_d_p   = bessel_d.deriv()
bessel_mu_p  = bessel_mu.deriv()



# Trigonometry using degree, as in Meeus
def sin(w):         return np.sin( np.deg2rad(w) )
def cos(w):         return np.cos( np.deg2rad(w) )
def tan(w):         return np.tan( np.deg2rad(w) )
def atan2(y, x):    return np.rad2deg( np.atan2(y, x) )
def asin(x):        return np.rad2deg( np.asin(x) )
def acos(x):        return np.rad2deg( np.acos(x) )
def atan(x):        return np.rad2deg( np.atan(x) )
def sqrt(x):        return np.sqrt(x)
def sq(x):          return x*x



def fundamental_plane(t: float, rho_sin_phi: float, rho_cos_phi: float, longitude: float, delta_t: float) \
    -> Tuple[float, float, float, float, float, float, float, float, float, float, float, float]:
    """
    Calculate coordinates and derivatives for observer location at specified time (TT)

    Parameters
    ----------
    t : float
        TT - Terrestial Time
    rho_sin_phi_p : float
        Observer location, rho*sin(phi') part in earth radii
    rho_cos_phi_p : float
        Observer location, rho*cos(phi') part in earth radii
    longitude : float
        Longitude part in degrees
    delta_t : float
        Delta T in seconds

    Returns
    -------
    Tuple[float, float, float, float, float, float, float, float, float]
        u, v : observer position relative to shadow axis x, y
        u_p, v_p : derivatives of U, V
        L1 : penumbra size at zeta
        L2 : umbra size at zeta
        xi, eta, zeta : position of observer over fundamental plane
    """
    x    = bessel_x(t)
    y    = bessel_y(t)
    d    = bessel_d(t)
    mu   = bessel_mu(t)
    l1   = bessel_l1(t)
    l2   = bessel_l2(t)
    # derivatives
    x_p  = bessel_x_p(t)
    y_p  = bessel_y_p(t)
    d_p  = bessel_d_p(t)
    mu_p = bessel_mu_p(t)

    # local coordinates in fundamental plane
    # https://de.wikipedia.org/wiki/Besselsche_Elemente
    # \theta = ( \mu - 1.002738 * 360° / 86400 s * \Delta T ) + \lambda
    # Sidereal day on Earth is approximately 86164.0905 s
    # 1.002738 = 86400 s / 86164.09 s
    #
    # local coordinates in fundamental plane
    # longitude = \lambda
    # theta = local hourangle corrected for TT
    theta = mu - 360 / 86164.0905 * delta_t  + longitude
    # [ESAA] (11.41)
    xi    = rho_cos_phi * sin(theta)
    eta   = rho_sin_phi * cos(d) - rho_cos_phi * cos(theta) * sin(d)
    zeta  = rho_sin_phi * sin(d) + rho_cos_phi * cos(theta) * cos(d)

    # new straight forward, proper derivatives, easier to understand ;-)
    # np.deg2rad because we're doing the calculations with degrees, not radians
    xi_p  = np.deg2rad( rho_cos_phi * cos(theta) * mu_p )
    eta_p = np.deg2rad(-rho_sin_phi * sin(d) * d_p                    # == -zeta * d_p
                       -rho_cos_phi * cos(theta) * cos(d) * d_p
                       +rho_cos_phi * sin(theta) * mu_p * sin(d) )    # == xi * mu_p * sin(d)

    # [ESAA] (11.119)
    # m = [u, v, 0] = [x-xi, y-eta, 0] = r - rho
    # n = \dot m = [u_p, v_p, 0]
    # max eclipse at m*n = 0
    # C1/C4 at u^2 + v^2 - L1^2 = 0
    # C2/C3 at u^2 + v^2 - L2^2 = 0
    u   = x   - xi                   # coordinates relative to shadow axis
    v   = y   - eta 
    u_p = x_p - xi_p                 # derivates of U, V
    v_p = y_p - eta_p
    L1  = l1  - zeta * bessel_tanf1  # penumbra size at zeta
    L2  = l2  - zeta * bessel_tanf2  # umbra size at zeta

    return u, v, u_p, v_p, L1, L2, xi, eta, zeta



def mag_and_pos_angle(u: float, v: float, L1: float, L2: float, xi: float, eta: float) -> Tuple[float, float, float, float]:
    """
    Compute magnitude, moon/sun size ration, position angles

    Parameters
    ----------
    u : float
        Position relative to shadow axis x
    v : float
        Position relative to shadow axis y
    L1 : float
        Size of penumbra
    L2 : float
        Size of umbra
    xi : float
        Observer position on fundamental plane x
    eta : float
        Observer position on fundamental plane y

    Returns
    -------
    Tuple[float, float, float, float]
        M_1 : magnitude
        M_2 : moon/sun size ratio = magnitude on center line
        Q : position angle from celestial north in degrees
        V : position angle from vertex ("up") in degrees
    """
    # Magnitude
    m   = sqrt(sq(u) + sq(v))
    M_1 = (L1 - m) / (L1 + L2)

    # Moon/sun size ration
    M_2 = (L1 - L2) / (L1 + L2)

    # [ESAA] (11.120)
    # Position angle on fundamental plane, in relation to celestial north
    # m = [L*sin(Q), L*cos(Q), 0]
    Q = atan2(u, v)
    if Q < 0: Q += 360

    # [ESAA] (11.121)
    # Position angle in relation to vertex, "up"
    # Parallatic angle tan(C) = xi/eta
    C = atan2(xi, eta)
    V = Q - C
    if V < 0: V += 360

    ic(M_1, M_2, C, V, Q)
    return M_1, M_2, V, Q



def geocentric(loc: EarthLocation) -> Tuple[float, float]:
    """
    Get geocentric parameters from location object

    Parameters
    ----------
    loc : EarthLocation
        Location of observer

    Returns
    -------
    Tuple[float, float]
        rho_sin_phi : rho*sin(phi') parameter in earth radii
        rho_cos_phi : rho*cos(phi') parameter in earth radii
    """
    # Var 1: geocentric from EarthLocation X, Y, Z
    phi_p = atan( loc.z.value / sqrt( sq(loc.x.value) + sq(loc.y.value) ) )
    rho = sqrt( sq(loc.x.value) + sq(loc.y.value) + sq(loc.z.value) ) / R_earth.value
    rho_sin_phi = rho * sin(phi_p)
    rho_cos_phi = rho * cos(phi_p)
    ic(phi_p, rho, rho_sin_phi, rho_cos_phi)

    # # Var 2: geocentric coordinates
    # # https://de.wikipedia.org/wiki/Besselsche_Elemente
    # e2 = 1 - sq(R_earth_pol) / sq(R_earth)
    # C = 1 / sqrt(1 - e2 * sq(sin(latitude)))
    # S = (1 - e2) * C
    # phi_p = atan((R_earth.value*S + height) / (R_earth.value*C + height) * tan(latitude) )
    # rho = (R_earth.value*C + height) * cos(latitude) / (R_earth.value * cos(phi_p))
    # rho_sin_phi = rho * sin(phi_p)
    # rho_cos_phi = rho * cos(phi_p)
    # ic(e2, C, S, phi_p, rho, rho_sin_phi, rho_cos_phi)

    return rho_sin_phi, rho_cos_phi



def sun_alt(t: Time, loc: EarthLocation) -> Angle:
    """
    Get altitude of the Sun for time and locations

    Parameters
    ----------
    t : Time
        Time of observation
    loc : EarthLocation
        Observer location

    Returns
    -------
    Angle
        Sun's altitude
    """
    sun = get_sun(t)
    altaz = sun.transform_to( AltAz(obstime=t, location=loc) )
    alt = altaz.alt
    ic(t, sun, alt)

    return alt



def bessel3(delta_t: float, loc: EarthLocation) -> None:
    """
    Compute local circumstances at location

    Parameters
    ----------
    delta_t : float
        Delta T value for the calculation
    loc : EarthLocation
        Observer position
    """
    longitude = loc.lon.value
    latitude  = loc.lat.value
    height    = loc.height.value
    ic(longitude, latitude, height)

    # geocentric coordinates
    rho_sin_phi, rho_cos_phi = geocentric(loc)

    # [ESAA] (11.119)
    # iterate towards m*n = [u, v, 0]*[u_p, v_p, 0] = 0
    # start at T0
    t = 0
    for _ in range(5):  # 5 iterations are enough
        u, v, u_p, v_p, L1, L2, xi, eta, _ = fundamental_plane(t, rho_sin_phi, rho_cos_phi, longitude, delta_t)
        ic(t, u, v, u_p, v_p, L1, L2)
        t_corr = -(u * u_p + v * v_p) / (sq(u_p) + sq(v_p)) # -m*n / n^2
        ic(t_corr)
        t = t + t_corr

    t_max        = t
    time_max_tt  = time_T0 + t_max * unit.hour
    time_max_ut1 = time_max_tt.ut1
    ic(t_max, time_max_tt, time_max_ut1)

    M_1, M_2, Q, V = mag_and_pos_angle(u, v, L1, L2, xi, eta)
    message(f"Magnitude: {M_1*100:.1f}%")
    message(f"Moon/sun size ratio: {M_2:.3f}")
    message(f"Time of MAX eclipse (UT1): {time_max_ut1}")
    message(f"Position angle at MAX (north): {V:.1f} deg")
    message(f"Position angle at MAX (zenith): {Q:.1f} deg")
    alt = sun_alt(time_max_ut1, loc)
    message(f"Sun altitude at MAX {alt:.2f}")


    if Options.list:
        # Results for 2.5h centered around tmax
        message("========================================================")
        message("Time                     magnitude  pos angle     sun")
        message("(UT1)                               north  up     alt")
        message("-----------------------  ---------  -----  -----  ------")

        if Options.totality:
            # Higher time resolution +/- 4.5 min around MAX
            t_range = np.linspace(-0.075, 0.075, 270+1)
        else:
            # +/- 1.25 h around MAX
            t_range = np.linspace(-1.25, 1.25, 150+1)

        for t in t_max + t_range:
            time_tt  = time_T0 + t * unit.hour
            time_ut1 = time_tt.ut1
            u, v, _, _, L1, L2, xi, eta, _ = fundamental_plane(t, rho_sin_phi, rho_cos_phi, longitude, delta_t)
            M_1, M_2, Q, V = mag_and_pos_angle(u, v, L1, L2, xi, eta)
            if M_1 >= 0 or not Options.pos_mag:
                alt = sun_alt(time_ut1, loc).degree
                message(f"{time_ut1}  {M_1*100:5.1f}%    {V:4.0f}°   {Q:4.0f}°  {alt:5.1f}°")



def main():
    """
    Main function of script
    """
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = DESCRIPTION,
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")

    arg.add_argument("-L", "--list", action="store_true", help="list time and magnitude centered around MAX")
    arg.add_argument("-0", "--positive-mag-only", action="store_true", help="list positive magnitude only")
    arg.add_argument("-T", "--totality", action="store_true", help="list with higher time resolution around MAX")
    arg.add_argument("-t", "--time", help=f"time (UT1), default T0={time_T0.ut1}")
    arg.add_argument("-l", "--location", help=f"coordinates, named location or MPC station code, default {DEFAULT_LOCATION}")
    arg.add_argument("--tse2026", action="store_true", help="test case TSE 12 Aug 2026")

    args = arg.parse_args()

    if args.debug:
        ic.enable()
        ic(sys.version_info, sys.path, args)
    if args.verbose:
        verbose.set_prog(NAME)
        verbose.enable()

    Options.list = args.list
    Options.pos_mag = args.positive_mag_only
    Options.totality = args.totality

    # Location and time
    loc = None
    time = None
    if args.tse2026:
        loc, time = test_tse2026()
    if args.location:
        loc = get_location(args.location)
    if loc == None:
        error("no location specified")
    ic(loc, loc.to_geodetic())
    Options.loc = loc
    verbose(f"location {location_to_string(loc)}")
    if args.time:
        time = Time(args.time, scale="ut1", location=loc)
    elif not time:
        time = Time(time_T0.ut1, location=loc)

    # Delta T = TT - UT1
    # https://de.wikipedia.org/wiki/Delta_T
    # https://en.wikipedia.org/wiki/%CE%94T_(timekeeping)
    # https://en.wikipedia.org/wiki/Leap_second
    # https://www.iers.org/IERS/EN/DataProducts/tools/timescales/timescales.html
    # (slight difference in UT1, UTC - UT1 ~ 0.1 s)
    delta_t = ((time.tt.jd - time.ut1.jd) * unit.day).to(unit.s)
    # alternate calculation (37 = leap seconds in 2025)
    delta_t2 = ((time.utc.jd - time.ut1.jd) * unit.day).to(unit.s) + (37 + 32.184) * unit.s
    ic(time, delta_t, delta_t2)

    verbose(f"time {time} (UT1), {time.tt} (TT), Delta T={delta_t:.2f}")
    verbose(f"time T0 {time_T0} (TT), {time_T0.ut1} (UT1)")

    # Run calculation for local circumstances
    bessel3(delta_t.value, loc)


if __name__ == "__main__":
    main()
