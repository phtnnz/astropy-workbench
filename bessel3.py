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
# Version 0.5 / 2025-10-27
#       Started to include calculation for central line
# Version 0.6 / 2025-10-30
#       Added duration and width for central line, new options -C --central-line,
#       -P --plot, changed -0 --positive-mag-only to -p
# Version 0.7 / 2025-11-05
#       Added output of greatest eclipse (for -C --central-line)
#
# See [ESAA] Explanatory Supplement to the Astronomical Almanac, 3rd Edtion
# Chapter 11 - Eclipses of the Sun and Moon


VERSION     = "0.7 / 2024-11-05"
AUTHOR      = "Martin Junius"
NAME        = "bessel3"
DESCRIPTION = "Solar eclipse local circumstances and central line"

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
from astropy.units       import Quantity
import astropy.units as unit
import astropy.constants as const
from astropy.coordinates import solar_system_ephemeris
from astropy.coordinates import get_sun

import numpy as np
from numpy.polynomial    import Polynomial

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# SciPy
from scipy import optimize

# Local modules
from verbose import verbose, warning, error, message
from astroutils import location_to_string, get_location     # ...

DEFAULT_LOCATION = "Burgos, Spain"



# Earth equatorial and polar radius
# IERS (2003) values (Wikipedia)
# https://en.wikipedia.org/wiki/Earth
# https://en.wikipedia.org/wiki/Earth_ellipsoid
R_earth = 6378136.6 * unit.m                    # a
R_earth_pol = 6356751.9 * unit.m                # b
f_earth = (R_earth - R_earth_pol) / R_earth     # flattening
e2_earth = 1 -  R_earth_pol**2 / R_earth**2     # ellipticity squared
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
# https://eclipsewise.com/solar/SEprime/2001-2100/SE2026Aug12Tprime.html
#
# x, y   - Cartesian coordinates of the lunar shadow axis in the Fundamental Plane 
#          (in units of Earth's equatorial radius)
# L1, L2 - Radii of the Moon's penumbral and umbral/antumbral shadows in the Fundamental Plane 
#          (in units of Earth's equatorial radius)
# d      - Declination of the Moon's shadow axis on the celestial sphere
# µ      - Hour angle of the Moon's shadow axis on the celestial sphere
# f1, f2 - Angles of the penumbral and umbral/antumbral shadow cones with respect to the axis of the lunar shadow

eclipse_name = "TSE 12 Aug 2026"
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
    # Sidereal day (unit.sday) on Earth is approximately 86164.0905 s
    # 1.002738 = 86400 s / 86164.09 s
    #
    # local coordinates in fundamental plane
    # longitude = \lambda
    # theta = local hourangle corrected for TT
    sday  = unit.sday.to(unit.s)
    theta = mu - 360 / sday * delta_t  + longitude
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



def fundamental_plane_xy(t: float) -> float:
    """
    Compute distance of (x, y) from earth center

    Parameters
    ----------
    t : float
        Time in hours relative to T0

    Returns
    -------
    float
        Distance of (x, y)
    """
    x    = bessel_x(t)
    y    = bessel_y(t)

    return sqrt(sq(x) + sq(y))



def find_greatest_eclipse() -> float:
    """
    Compute time of greatest eclipse

    Returns
    -------
    float
        Time of greatest eclipse relative to T0 in hours
    """
    # Find mininum of (x, y) distance, which is greatest eclipse
    func    = lambda x: fundamental_plane_xy(x[0])
    # +/- 1 hour around T0
    sol_max = optimize.minimize(func, x0=(0), bounds=[(-1, 1)])
    t_max  = sol_max.x[0]
    ic(sol_max, t_max)

    return t_max



def solve_quadrant(sin_theta: float, cos_theta: float) -> float:
    """
    Find angle in one of the four quadrants

    Parameters
    ----------
    sin_theta : float
        sin value
    cos_theta : float
        cos value

    Returns
    -------
    float
        Angle in degrees
    """
    if (sin_theta>=0 and cos_theta>=0): return  asin(sin_theta)
    if (sin_theta <0 and cos_theta>=0): return  asin(sin_theta)
    if (sin_theta <0 and cos_theta <0): return -acos(cos_theta)
    if (sin_theta>=0 and cos_theta <0): return  acos(cos_theta)
    # not reached



def fundamental_plane2(t: float, delta_t: float) -> Tuple[Quantity, Quantity, Quantity, Quantity, Quantity]:
    """
    Calculation for central line on fundamental plane

    Parameters
    ----------
    t : float
        Time relative to T0 in TT
    delta_t : float
        Delta T in s

    Returns
    -------
    Tuple[float, float, float, float, float]
        longitude : longitude coord of central eclipse at T in degrees
        latitude : latitude coord of central eclipse at T in degrees
        M_2 : moon/sun size ration at T
        duration : duration of total/annular eclipse at T in s
        width : width of umbra at T in km
    """
    x    = bessel_x(t)
    y    = bessel_y(t)
    d    = bessel_d(t)
    mu   = bessel_mu(t)
    l1   = bessel_l1(t)
    l2   = bessel_l2(t)
    ic(x, y, d, mu, l1, l2)

    # derivatives
    x_p  = bessel_x_p(t)
    y_p  = bessel_y_p(t)
    d_p  = bessel_d_p(t)
    mu_p = bessel_mu_p(t)

    # Auxiliary Besselians [ESAA] (11.61)
    rho_1 = sqrt( 1 - e2_earth * sq(cos(d)) )
    rho_2 = sqrt( 1 - e2_earth * sq(sin(d)) )
    sin_d_1 = sin(d) / rho_1
    cos_d_1 = sqrt(1 - e2_earth) * cos(d) / rho_1
    sin_d_1_d_2 = e2_earth * sin(d) * cos(d) / rho_1
    cos_d_1_d_2 = sqrt(1 - e2_earth) / rho_1 / rho_2
    ic(rho_1, rho_2, sin_d_1, cos_d_1, sin_d_1_d_2, cos_d_1_d_2)

    # Convert point(xi, eta, 0) on fundamental plain [ESAA] 11.3.3.3
    ## For central line ##
    xi   = x
    eta  = y
    ic(xi, eta)
    #####################
    eta_1 = eta / rho_1                                         # (11.55)
    if 1 - sq(xi) - sq(eta_1) < 0:                              # no result
        return None, None, None, None, None
    zeta_1 = sqrt(1 - sq(xi) - sq(eta_1))                       # (11.56)
    ic(eta_1, zeta_1)

    phi_1 = asin(eta_1 * cos_d_1 + zeta_1*sin_d_1)              # (11.59) 2nd row
    sin_theta = xi / cos(phi_1)                                 #         1st row
    cos_theta = (-eta_1*sin_d_1 + zeta_1*cos_d_1) / cos(phi_1)  #         3rd row
    ic(phi_1, sin_theta, cos_theta)

    # Geodesic coordinates
    theta = solve_quadrant(sin_theta, cos_theta)
    # longitude = lambda
    # theta = local hourangle corrected for TT
    sday  = unit.sday.to(unit.s)
    ic(sday, theta, mu, delta_t)
    longitude = theta - mu * unit.deg + 360 / sday * delta_t * unit.deg
    if longitude > 180 * unit.deg:
        longitude -= 360 * unit.deg
    if longitude < -180 * unit.deg:
        longitude += 360 * unit.deg
    # latitude = phi
    latitude = atan( 1 / sqrt(1 - e2_earth) * tan(phi_1) )      # (11.52) divide eqs
    zeta = rho_2 * (zeta_1*cos_d_1_d_2 - eta_1*sin_d_1_d_2)     # (11.60)
    ic(longitude, latitude, zeta)

    L1  = l1  - zeta * bessel_tanf1  # penumbra size at zeta
    L2  = l2  - zeta * bessel_tanf2  # umbra size at zeta
    M_2 = (L1 - L2) / (L1 + L2)      # magnitude at central line = moon/sun size ratio
    ic(L1, L2, M_2)

    # [ESAA] 11.3.5.5
    # Central line, duration of central eclipse, and width of path
    # Duration
    x_dot    = x_p
    y_dot    = y_p
    d_dot    = np.deg2rad(d_p)                                  # d, mu are in degree
    mu_dot   = np.deg2rad(mu_p)
    xi_dot   = mu_dot * ( -y*sin(d) + zeta*cos(d) )             # (11.99)
    eta_dot  = mu_dot * x * sin(d) - d_dot * zeta
    n2       = sq(x_dot - xi_dot) + sq(y_dot - eta_dot)
    n        = sqrt(n2)
    L2       = abs(L2)                                          # < 0 for TSE, > 0 for ASE
    duration = 2 * L2 / n * unit.hour                           # (11.100)
    ic(n2, n, duration.to(unit.s))

    # Width
    width = 2 * L2 / sqrt(sq(zeta) +                            # (11.101)
                          sq(xi  / n * (x_dot - xi_dot ) +
                             eta / n * (y_dot - eta_dot)  ) ) * R_earth.to(unit.km)
    ic(width)

    return longitude, latitude, M_2, duration.to(unit.s), width



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
    # # See also [ESAA] (11.47-49)
    # C = 1 / sqrt(1 - e2_earth * sq(sin(latitude)))
    # S = (1 - e2) * C
    # phi_p = atan( (R_earth.value*S + height) / (R_earth.value*C + height) * tan(latitude) )
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



def central_line(delta_t: float, plot: bool=False) -> None:
    """
    List and plot coordinates of central line

    Parameters
    ----------
    delta_t : float
        Delta T in s
    """
    message("===================================================================================")
    message("Time (UT1)               longitude      latitude      magnitude  duration  width   ")
    message("-----------------------  -------------  ------------  ---------  --------  --------")

    lons = []
    lats = []
    times = []

    ##TEST##
    # t_range = [0]

    # T0 +/- 2h, 1 min intervals
    t_range = np.linspace(-2, 2, 4*60+1)
    for t in t_range:
        time_tt  = time_T0 + t * unit.hour
        time_ut1 = time_tt.ut1
        lon, lat, M_2, duration, width = fundamental_plane2(t, delta_t)
        if lon != None:
            message(f"{time_ut1}  {lon:9.4f}  {lat:8.4f}  {M_2:.5f}    {duration:5.1f}   {width:5.1f}")
            times.append(time_ut1)
            lons.append(lon.value)
            lats.append(lat.value)

    # Greatest eclipse
    t_max = find_greatest_eclipse()
    lon, lat, M_2, duration, width = fundamental_plane2(t_max, delta_t)
    time_tt  = time_T0 + t_max * unit.hour
    time_ut1 = time_tt.ut1
    message("----- greatest eclipse ------------------------------------------------------------")
    message(f"{time_ut1}  {lon:9.4f}  {lat:8.4f}  {M_2:.5f}    {duration:5.1f}   {width:5.1f}")
    
    if plot:
        ic(lons, lats)
        plot_path(times, lons, lats)



def plot_path(times: list[Time], lons: list[Quantity], lats: list[Quantity]) -> None:
    """
    Plot eclipse path on the Earth

    Parameters
    ----------
    times : list[Time]
        Time labels
    lons : list[Quantity]
        Longitude coordinates
    lats : list[Quantity]
        Latitude coordinates
    """
    verbose("plotting eclipse central line")
    fig = plt.figure(figsize=(10, 10), dpi=300)
    ax  = fig.add_subplot()

    c_lon = lons[len(lons) // 2]
    c_lat = lats[len(lats) // 2]

    # lon_0, lat_0 are the center point of the projection.
    map = Basemap(projection='ortho', lon_0=c_lon, lat_0=c_lat, ax=ax)
    # map.drawcoastlines()
    map.fillcontinents(color="wheat", lake_color="lightblue")

    # Draw parallels and meridians.
    map.drawparallels(np.linspace(-60, 60, 5))
    map.drawmeridians(np.linspace(-180, 150, 12))
    map.drawmapboundary(fill_color="lightblue")
    map.drawcountries()

    ax.set_title(eclipse_name, fontsize=20)
    # plt.title(eclipse_name)

    # Plot eclipse path
    x, y = map(lons, lats)
    map.plot(x, y, color="red", linewidth=5)

    # Add time annotations
    format = "%H:%M UT1"
    t = times[0]
    lon = lons[0]
    lat = lats[0]
    plt.annotate(t.strftime(format), xy=map(lon, lat), ha="center", va="center", color="black", fontsize=10)
    t = times[-1]
    lon = lons[-1]
    lat = lats[-1]
    plt.annotate(t.strftime(format), xy=map(lon, lat), ha="center", va="center", color="black", fontsize=10)

    plt.savefig("tmp/plot.png", bbox_inches="tight")
    plt.close()



def local_circumstances(delta_t: float, loc: EarthLocation) -> None:
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
                message(f"{time_ut1}  {M_1*100:5.1f}%    {V:4.0f}°  {Q:4.0f}°  {alt:5.1f}°")



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

    arg.add_argument("-C", "--central-line", action="store_true", help="list coordinates of central line")
    arg.add_argument("-P", "--plot", action="store_true", help="plot central line")
    arg.add_argument("-L", "--list", action="store_true", help="list time and magnitude centered around MAX")
    arg.add_argument("-p", "--positive-mag-only", action="store_true", help="list positive magnitude only")
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

    # Debug some constants
    ic(R_earth, R_earth_pol, f_earth, e2_earth, R_moon)
    ic(1 - f_earth)
    ic(sqrt(1-e2_earth))

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

    if args.central_line:
        # List coordinates for central line
        central_line(delta_t.value, args.plot)
    else:
        # Run calculation for local circumstances
        local_circumstances(delta_t.value, loc)



if __name__ == "__main__":
    main()
