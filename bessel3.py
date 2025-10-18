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
# Version 0.0 / 2025-10-18
#       Included Uwe Pilz' code mostly verbatim

VERSION     = "0.0 / 2024-xx-xx"
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



# Command line options
class Options:
    loc: EarthLocation = None       # -l --location



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



##### Following code taken from Uwe Pilz, Februar 2025 #####
##### https://fg-astrophysik.vdsastro.de/prg95         #####
##### No license provided                              #####
##### Modified to use AstroPy & Friends                #####

# VdS-Journal 95: Lokaler Verlauf einer Finsternis
# Die Berechnung erfolgt minütlich. Vor Benutzung müssen die sog. Bessel-Elemente der 
# gewünschten Finsternis eingesetzt werden. Die Datumsangabe muss aktiviert 
# ("einkommentiert") werden.

# Geographische Koordinaten werden astronomisch richtig ausgegebe, also Ost = negativ.
# Das ist umgekehrt zur Google-Maps-Anzeige.

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

T0 = 18
bessel_x       = [  0.4755140,  0.5189249, -0.0000773, -0.0000080 ]
bessel_y       = [  0.7711830, -0.2301680, -0.0001246,  0.0000038 ]
bessel_d       = [ 14.7966700, -0.0120650, -0.0000030,  0.0       ]
bessel_l1      = [  0.5379550,  0.0000939, -0.0000121,  0.0       ]
bessel_l2      = [ -0.0081420,  0.0000935, -0.0000121,  0.0       ]
bessel_u       = [ 88.7477900, 15.0030900,  0.0000000,  0.0       ]
bessel_tanf1   = 0.0046141
bessel_tanf2   = 0.0045911

polynomial_x   = Polynomial(bessel_x)
polynomial_y   = Polynomial(bessel_y)
polynomial_d   = Polynomial(bessel_d)
polynomial_l1  = Polynomial(bessel_l1)
polynomial_l2  = Polynomial(bessel_l2)
polynomial_u   = Polynomial(bessel_u)

polynomial_x_p = polynomial_x.deriv()
polynomial_y_p = polynomial_y.deriv()


# Trigonometry using degree, as in Meeus
def sin(w):         return np.sin( np.deg2rad(w) )
def cos(w):         return np.cos( np.deg2rad(w) )
def tan(w):         return np.tan( np.deg2rad(w) )
def atan2(y, x):    return np.rad2deg( np.atan2(y, x) )
def asin(x):        return np.rad2deg( np.asin(x) )
def acos(x):        return np.rad2deg( np.acos(x) )
def atan(x):        return np.rad2deg( np.atan(x) )


# Ergebnisse der Fundamentalebene für einen Zeitpunkt
def magnitudePoswinkel(t: float, pSinPS: float, pCosPS: float, longitude: float, delta_t: float) -> Tuple[float, float, float, float, float, float, float, float, float]:
    X  = polynomial_x(t)
    Y  = polynomial_y(t)
    D  = polynomial_d(t)
    M  = polynomial_u(t)
    L1 = polynomial_l1(t)
    L2 = polynomial_l2(t)

    # derivatives
    X_p = polynomial_x_p(t)
    Y_p = polynomial_y_p(t)

    # Koordinaten bezüglich der Fundamentalebene

    # https://de.wikipedia.org/wiki/Besselsche_Elemente
    # theta = ( \mu - 1.002738 * 360° / 86400 s * Delta_T ) + \Lambda
    # Sidereal day on Earth is approximately 86164.0905 s
    # 1.002738 = 86400 s / 86164.09 s
    # theta = local hourangle corrected for TT
    #
    # H: theta, M: \mu, longitude = \Lambda
    H     = M - 360 / 86164.0905 * delta_t  + longitude
    xi    = pCosPS * sin(H)
    eta   = pSinPS * cos(D) - pCosPS * cos(H) * sin(D)
    zeta  = pSinPS * sin(D) + pCosPS * cos(H) * cos(D)
    # stündliche Änderungen
    xi_p  = np.deg2rad(bessel_u[1] * pCosPS * cos(H))
    eta_p = np.deg2rad(bessel_u[1] * xi * sin(D) - zeta * bessel_d[1])

    U     = X - xi
    V     = Y - eta
    a     = X_p - xi_p
    b     = Y_p - eta_p
    L1S   = L1 - zeta * bessel_tanf1
    L2S   = L2 - zeta * bessel_tanf2
    n2    = a * a + b * b

    return a, b, U, V, n2, L1S, L2S, D, H



# Ergebnisse für die Ausgabe
def ergebnisberechnung(t: float, U: float, V: float, L1S: float, L2S: float, D: float, H: float, latitude: float, delta_t: float) -> Tuple[float, float, float, float, float, float, float]:
    # Zeit des lokalen Maximums
    TD = T0 + t
    # Magnitude G
    m = np.sqrt(U * U + V * V)
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



def bessel3(longitude: float, latitude: float, hoehe: float, delta_t: float) -> None:
    ic(longitude, latitude, hoehe)

    # Hauptprogramm
    t = 0

    # rechtwinklige geozentrische Koordinaten
    U = atan(0.99664719 * tan(latitude))
    pSinPS = 0.99664719 * sin(U) + hoehe / 6378140 * sin(latitude)
    pCosPS = cos(U) + hoehe / 6378140 * cos(latitude)

    # 1) Ergebnisse zum Maximums-Zeitpunkt
    for i in range(5):  # 5 Iterationen genügen
        a, b, U, V, n2, L1S, L2S, D, H = magnitudePoswinkel(t, pSinPS, pCosPS, longitude, delta_t)
        # Zeitpunkt der maximalen Finsternis per Iteration
        tm = -(U * a + V * b) / n2
        t = t + tm
        UTh, UTm, UTs, G, A, Zm, P = ergebnisberechnung(t, U, V, L1S, L2S, D, H, latitude, delta_t)
    print("Uhrzeit des Maximums:", UTh, "h", UTm, "m", UTs, "s UT")  # dynamische Zeit
    print("Magnitude: ", int(1000 * G) / 10, "%")
    print("Verhältnis Durchmesser Mond/Sonne: ", int(1000 * A) / 1000)
    print("Positionswinkel, bezogen auf Nord: ", int(P+0.5), "Grad")
    print("Positionswinkel, bezogen auf Zenit: ", int(Zm+0.5), "Grad")

    # 2) Ergebnisse aller 10 Minuten
    tmax = t
    t = tmax - 3 
    print("UT             Magnitude Nord-  Zenit- Posw.")
    for i in range(480):  # 3 h aller Minuten
        a, b, U, V, n2, L1S, L2S, D, H = magnitudePoswinkel(t, pSinPS, pCosPS, longitude, delta_t)
        UTh, UTm, UTs, G, A, Zm, P = ergebnisberechnung(t, U, V, L1S, L2S, D, H, latitude, delta_t)
        if G >= 0:
            print("%2.0f h %02.0f m %02.0f s  %5.1f%%     %3.0f° %3.0f°" % (UTh, UTm, UTs, 100*G, P+0.5, Zm+0.5))
            # print(UTh,"h", UTm,"m", UTs,"s  ", int(1000 * G) / 10, "%  ", int(P+0.5),"°  ", int(Zm+0.5), "°")
        t = t + 1 / 60

#/ Uwe Pilz, Februar 2025



def main():
    arg = argparse.ArgumentParser(
        prog        = NAME,
        description = DESCRIPTION,
        epilog      = "Version " + VERSION + " / " + AUTHOR)
    arg.add_argument("-v", "--verbose", action="store_true", help="verbose messages")
    arg.add_argument("-d", "--debug", action="store_true", help="more debug messages")

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
    # Delta T from astropy Time
    # checked with https://www.iers.org/IERS/EN/DataProducts/tools/timescales/timescales.html
    # slight difference in UT1
    delta_t = ((time.tt.jd - time.ut1.jd) * u.day).to(u.s)
    # approximation, see https://de.wikipedia.org/wiki/Delta_T
    y = time.to_datetime().year + (time.to_datetime().month - 0.5) / 12
    delta_t2 = 67.62 + 0.3645 * (y - 2015) + 0.0039755 * (y - 2015)**2
    ic(time, delta_t, y, delta_t2)
    verbose(f"time UTC {time}, Delta T={delta_t:.2f}")

    # Run calculation for local circumstances
    bessel3(loc.lon.value, loc.lat.value, loc.height.value, delta_t.value)


if __name__ == "__main__":
    main()
