#!/usr/bin/env python

# Copyright 2025 Martin Junius, Uwe Pilz (VdS) and contributors
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
#       Major overhault to fit into the astropy-workbench environment,
#       using astropy and numpy, a lot of renaming for better understanding

VERSION     = "0.1 / 2024-10-19"
AUTHOR      = "Martin Junius"
NAME        = "bessel3"
DESCRIPTION = "Eclipse local circumstances"

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
import numpy as np
from numpy.polynomial import Polynomial

# Local modules
from verbose import verbose, warning, error, message
from astroutils import location_to_string, get_location     # ...

DEFAULT_LOCATION = "Burgos, Spain"



# Earth equatorial radius
R_earth = 6378137.0 * u.m           # GRS 80/WGS 84 value (Wikipedia)
                                    # https://en.wikipedia.org/wiki/World_Geodetic_System
R_earth_pol = 6356752.0 * u.m       # https://en.wikipedia.org/wiki/Earth
f_earth = (R_earth - R_earth_pol) / R_earth
# Moon equatorial radius
R_moon = 1738100.0 * u.m
R_moon = 0.2725076 * R_earth        # IAU 1982
R_moon = 0.272281  * R_earth        # smaller value from https://eclipse.gsfc.nasa.gov/SEpubs/20080801/TP214149b.pdf



# Command line options
class Options:
    loc: EarthLocation = None       # -l --location
    list: bool = False              # -L --list



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
    loc = EarthLocation(lon=-3.68935*u.degree, lat=42.35047*u.degree, height=891*u.m)
    time = Time("2026-08-12 18:29:17.6", location=loc)
    ic(loc, time, time.jd)
    return loc, time



# A major part of the following code ist originally by Uwe Pilz
# VdS-Journal 95: Lokaler Verlauf einer Finsternis
# https://fg-astrophysik.vdsastro.de/prg95
# (no license provided)
# Modified to use numpy and astropy
#
# sofiBessel3 : lokaler Verlauf einer Sonnenfinsternis
# ohne Berücksichtigung der lokalen Sonnenhöhe

# Besselian element for TSE 12 Aug 2026 from
# https://eclipse.gsfc.nasa.gov/SEsearch/SEdata.php?Ecl=20260812
#
# x, y   - Cartesian coordinates of the lunar shadow axis in the Fundamental Plane 
#          (in units of Earth's equatorial radius)
# L1, L2 - Radii of the Moon's penumbral and umbral/antumbral shadows in the Fundamental Plane 
#          (in units of Earth's equatorial radius)
# d      - Declination of the Moon's shadow axis on the celestial sphere
# µ      - Hour angle of the Moon's shadow axis on the celestial sphere
# f1, f2 - Angles of the penumbral and umbral/antumbral shadow cones with respect to the axis of the lunar shadow

bessel_T0 = 18
bessel_x        = [  0.4755140,  0.5189249, -0.0000773, -0.0000080 ]
bessel_y        = [  0.7711830, -0.2301680, -0.0001246,  0.0000038 ]
bessel_d        = [ 14.7966700, -0.0120650, -0.0000030,  0.0       ]
bessel_l1       = [  0.5379550,  0.0000939, -0.0000121,  0.0       ]
bessel_l2       = [ -0.0081420,  0.0000935, -0.0000121,  0.0       ]
bessel_mu       = [ 88.7477900, 15.0030900,  0.0000000,  0.0       ]
bessel_tanf1    = 0.0046141
bessel_tanf2    = 0.0045911

polynomial_x    = Polynomial(bessel_x)
polynomial_y    = Polynomial(bessel_y)
polynomial_d    = Polynomial(bessel_d)
polynomial_l1   = Polynomial(bessel_l1)
polynomial_l2   = Polynomial(bessel_l2)
polynomial_mu   = Polynomial(bessel_mu) 

polynomial_x_p  = polynomial_x.deriv()
polynomial_y_p  = polynomial_y.deriv()
polynomial_d_p  = polynomial_d.deriv()
polynomial_mu_p = polynomial_mu.deriv()



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



# Ergebnisse der Fundamentalebene für einen Zeitpunkt
def calc_on_fundamental_plane(t: float, rho_sin_phi_p: float, rho_cos_phi_p: float, longitude: float, delta_t: float) -> Tuple[float, float, float, float, float, float, float, float, float]:
    """
    Caculate coordinates and derivatives of observer location at specified time (TT)

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
        U, V : observer position relative to shadow axis x, y
        U_p, V_p : derivatives of U, V
        l1_zeta : penumbra size at zeta
        l2_zeta : umbra size at zeta
        d : declination of shadow axis in degrees
        theta : hourangle of shadow axis at observer longitude in degrees
    """
    x    = polynomial_x(t)
    y    = polynomial_y(t)
    d    = polynomial_d(t)
    mu   = polynomial_mu(t)
    l1   = polynomial_l1(t)
    l2   = polynomial_l2(t)
    # derivatives
    x_p  = polynomial_x_p(t)
    y_p  = polynomial_y_p(t)
    d_p  = polynomial_d_p(t)
    mu_p = polynomial_mu_p(t)

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
    xi    = rho_cos_phi_p * sin(theta)
    eta   = rho_sin_phi_p * cos(d) - rho_cos_phi_p * cos(theta) * sin(d)
    zeta  = rho_sin_phi_p * sin(d) + rho_cos_phi_p * cos(theta) * cos(d)
    # original derivatives approximation
    # xi_p  = np.deg2rad(bessel_mu[1] * rho_cos_phi_p * cos(theta))
    # eta_p = np.deg2rad(bessel_mu[1] * xi * sin(d) - zeta * bessel_d[1])
    # straight forward, proper derivatives, easier to understand ;-)
    xi_p  = np.deg2rad( rho_cos_phi_p * cos(theta) * mu_p )
    eta_p = np.deg2rad(-rho_sin_phi_p * sin(d) * d_p                    # == -zeta * d_p
                       -rho_cos_phi_p * cos(theta) * cos(d) * d_p
                       +rho_cos_phi_p * sin(theta) * mu_p * sin(d) )    # == xi * mu_p * sin(d)

    U       = x   - xi                  # coordinates relative to shadow axis
    V       = y   - eta
    U_p     = x_p - xi_p                # derivates of U, V
    V_p     = y_p - eta_p
    l1_zeta = l1 - zeta * bessel_tanf1  # penumbra size at zeta
    l2_zeta = l2 - zeta * bessel_tanf2  # umbra size at zeta

    return U, V, U_p, V_p, l1_zeta, l2_zeta, d, theta



# Ergebnisse für die Ausgabe
def ergebnisberechnung(t: float, U: float, V: float, L1S: float, L2S: float, D: float, H: float, latitude: float, delta_t: float) -> Tuple[float, float, float, float, float, float, float]:
    # Zeit des lokalen Maximums
    TD = bessel_T0 + t
    # Magnitude G
    m = sqrt(U * U + V * V)
    G = (L1S - m) / (L1S + L2S)
    # Durchmesserverhältnis
    A = (L1S - L2S) / (L1S + L2S)
    # Positionswinklel zur maximalen Verfinsterung, vom Nordpol der Sonne
    Pm = atan2(U , V)
    if Pm < 0:
        Pm = Pm + 360.0
    # Positionswinkel Zm bezogen auf die Zenit-Richtung
    sinH = sin(D) * sin(latitude) + cos(D) * cos(latitude) * cos(H)
    h = asin(sinH)
    sinq = (cos(latitude) * sin(H)) / cos(h)
    q = asin(sinq)
    Zm = Pm - q
    # Positionswinkel
    P = atan2(U, V)
    if P < 0:
        P = P + 360
    UT = TD - delta_t / 3600
    UTh = int(UT)
    UTm0 = 60 * (UT - UTh)
    UTm = int(UTm0)
#    UTs = int(60 * (UTm0 - UTm) + 0.5)
    UTs = 60 * (UTm0 - UTm)
    return UTh, UTm, UTs, G, A, Zm, P



def bessel3(delta_t: float) -> None:
    longitude = Options.loc.lon.value
    latitude  = Options.loc.lat.value
    height    = Options.loc.height.value
    ic(longitude, latitude, height)

    # geocentric from EarthLocation X, Y, Z
    phi_p = atan( Options.loc.z.value / sqrt( sq(Options.loc.x.value) + sq(Options.loc.y.value) ) )
    rho = sqrt( sq(Options.loc.x.value) + sq(Options.loc.y.value) + sq(Options.loc.z.value) ) / R_earth.value
    rho_sin_phi_p = rho * sin(phi_p)
    rho_cos_phi_p = rho * cos(phi_p)
    ic(phi_p, rho_sin_phi_p, rho_cos_phi_p)

    # geocentric coordinates
    # https://de.wikipedia.org/wiki/Besselsche_Elemente
    # e2 = 1 - sq(R_earth_pol) / sq(R_earth)
    # C = 1 / sqrt(1 - e2 * sq(sin(latitude)))
    # S = (1 - e2) * C
    # phi_p = atan((R_earth.value*S + height) / (R_earth.value*C + height) * tan(latitude) )
    # rho = (R_earth.value*C + height) * cos(latitude) / (R_earth.value * cos(phi_p))
    # rho_sin_phi_p = rho * sin(phi_p)
    # rho_cos_phi_p = rho * cos(phi_p)
    # ic(e2, C, S, phi_p, rho, rho_sin_phi_p, rho_cos_phi_p)

    # geocentric coordinates from prog95_3.py, slight difference 0.1° for phi_p (U)
    # rechtwinklige geozentrische Koordinaten
    # ratio_earth = R_earth_pol.value / R_earth.value
    # ratio_hoehe = height / R_earth.value
    # ic(ratio_earth, ratio_hoehe)
    # phi_p         = atan(ratio_earth * tan(latitude))
    # rho_sin_phi_p = ratio_earth * sin(phi_p) + ratio_hoehe * sin(latitude)
    # rho_cos_phi_p =               cos(phi_p) + ratio_hoehe * cos(latitude)
    # ic(phi_p, rho_sin_phi_p, rho_cos_phi_p)

    # Start at t = 0, relative to T0 from besselians
    # 1) Ergebnisse zum Maximums-Zeitpunkt
    t = 0
    for i in range(5):  # 5 Iterationen genügen
        U, V, a, b, L1S, L2S, D, H = calc_on_fundamental_plane(t, rho_sin_phi_p, rho_cos_phi_p, longitude, delta_t)
        ic(t, a, b, U, V, L1S, L2S, D, H)
        # Zeitpunkt der maximalen Finsternis per Iteration
        tm = -(U * a + V * b) / (sq(a) + sq(b))
        t = t + tm
    UTh, UTm, UTs, G, A, Zm, P = ergebnisberechnung(t, U, V, L1S, L2S, D, H, latitude, delta_t)
    print("Uhrzeit des Maximums:", UTh, "h", UTm, "m", UTs, "s UT")  # dynamische Zeit
    print("Magnitude: ", int(1000 * G) / 10, "%")
    print("Verhältnis Durchmesser Mond/Sonne: ", int(1000 * A) / 1000)
    print("Positionswinkel, bezogen auf Nord: ", int(P+0.5), "Grad")
    print("Positionswinkel, bezogen auf Zenit: ", int(Zm+0.5), "Grad")
    tmax = t

    if Options.list:
        # Results for 2.5h centered around tmax
        print("UT             Magnitude Nord-  Zenit- Posw.")
        for t in tmax + np.linspace(-1.25, 1.25, 150+1):
            U, V, a, b, L1S, L2S, D, H = calc_on_fundamental_plane(t, rho_sin_phi_p, rho_cos_phi_p, longitude, delta_t)
            UTh, UTm, UTs, G, A, Zm, P = ergebnisberechnung(t, U, V, L1S, L2S, D, H, latitude, delta_t)
            # if G >= 0:
            #     print("%2.0f h %02.0f m %02.0f s  %5.1f%%     %3.0f° %3.0f°" % (UTh, UTm, UTs, 100*G, P+0.5, Zm+0.5))
            print(f"{UTh:02.0f}h {UTm:02.0f}m {UTs:04.1f}s  {G*100:.1f}%    {P:.0f}°   {Zm:.0f}°")

#/ Uwe Pilz, Februar 2025



def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = DESCRIPTION,
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")

    arg.add_argument("-L", "--list", action="store_true", help="list time and magnitude centered around MAX")
    arg.add_argument("-t", "--time", help="time (UTC), default now")
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
        time = Time(args.time, location=loc)
    elif not time:
        time = Time(Time.now(), location=loc)

    # Delta T = TT - UT1
    # https://de.wikipedia.org/wiki/Delta_T
    # https://en.wikipedia.org/wiki/%CE%94T_(timekeeping)
    # https://en.wikipedia.org/wiki/Leap_second
    # https://www.iers.org/IERS/EN/DataProducts/tools/timescales/timescales.html
    # (slight difference in UT1, UTC - UTC1 ~ 0.1 s)
    # delta_t = ((time.tt.jd - time.ut1.jd) * u.day).to(u.s)
    ##FIXME: using Delta T = TT - UTC here, as it's also used to convert TT -> UTC
    delta_t = ((time.tt.jd - time.utc.jd) * u.day).to(u.s)
    # alternate calculation (37 = leap seconds in 2025)
    delta_t2 = ((time.utc.jd - time.ut1.jd) * u.day).to(u.s) + (37 + 32.184) * u.s
    ic(time, delta_t, delta_t2)

    verbose(f"time UTC {time}, Delta T={delta_t:.2f}")

    # Run calculation for local circumstances
    bessel3(delta_t.value)


if __name__ == "__main__":
    main()
