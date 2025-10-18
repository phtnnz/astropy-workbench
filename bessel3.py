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

from math import *

# Besselsche Elemente
# yr=2026; mnt=8; day=12
T0 = 18
bessel_x     = [ 0.4755140, 0.5189249, -0.0000773, -0.0000080]
bessel_y     = [ 0.7711830, -0.2301680, -0.0001246, 0.0000038]
bessel_d     = [ 14.7966700, -0.0120650, -0.0000030, 0.0]
bessel_l1    = [ 0.5379550, 0.0000939, -0.0000121, 0.0]
bessel_l2    = [ -0.0081420, 0.0000935, -0.0000121, 0.0]
bessel_u     = [ 88.7477900, 15.0030900, 0.0000000, 0.0]
bessel_tanf1 = 0.0046141
bessel_tanf2 = 0.0045911


G# Winkelfunktionen in Gradmaß, wie im Meeus
def Sin(w):    return sin(pi * w / 180)
def Cos(w):    return cos(pi * w / 180)
def Tan(w):    return tan(pi * w / 180)
def Atan2(y, x):    return 180 * atan2(y, x) / pi
def Asin(x):    return 180 * asin(x) / pi
def Acos(x):    return 180 * acos(x) / pi
def Atan(x):    return 180 * atan(x) / pi
def sq(x):    return x * x


def solveQuadrant(sinX, cosX):
    if sinX >= 0 and cosX >= 0:
        return Asin(sinX)
    if sinX < 0 and cosX >= 0:
        return Asin(sinX)
    if sinX < 0 and cosX < 0:
        return -Acos(cosX)
    if sinX >= 0 and cosX < 0:
        return Acos(cosX)


# Ergebnisse der Fundamentalebene für einen Zeitpunkt
def magnitudePoswinkel(t: float, pSinPS: float, pCosPS: float, Lambda: float, delta_t: float) -> Tuple[float, float, float, float, float, float, float, float, float]:
    X = 0
    Y = 0
    D = 0
    L1 = 0
    L2 = 0
    M = 0
    for j in range(3, -1, -1):
        X = X * t + bessel_x[j]  # die Koordinaten
        Y = Y * t + bessel_y[j]
        D = D * t + bessel_d[j]  # Deklination ...
        M = M * t + bessel_u[j]  # ... und Stundenwinkel im orstfesten Äquatorsystem
        L1 = L1 * t + bessel_l1[j]  # die Radien der Schattenkegel, Halbschatten ...
        L2 = L2 * t + bessel_l2[j]  # ... und Kernschatten

    # die Ableitungen / Veränderungen des Ortes i.f. Fundamental-E.
    XS = bessel_x[1] + 2 * bessel_x[2] * t + 3 * bessel_x[3] * t * t
    YS = bessel_y[1] + 2 * bessel_y[2] * t + 3 * bessel_y[3] * t * t

    # Koordinaten bezüglich der Fundamentalebene
    H = M - Lambda - 0.00417807 * delta_t
    xi = pCosPS * Sin(H)
    eta = pSinPS * Cos(D) - pCosPS * Cos(H) * Sin(D)
    zeta = pSinPS * Sin(D) + pCosPS * Cos(H) * Cos(D)
    # stündliche Änderungen
    xiS = 0.01745329 * bessel_u[1] * pCosPS * Cos(H)
    etaS = 0.01745329 * (bessel_u[1] * xi * Sin(D) - zeta * bessel_d[1])

    U = X - xi
    V = Y - eta
    a = XS - xiS
    b = YS - etaS
    L1S = L1 - zeta * bessel_tanf1
    L2S = L2 - zeta * bessel_tanf2
    n2 = a * a + b * b
    return a, b, U, V, n2, L1S, L2S, D, H


# Ergebnisse für die Ausgabe
def ergebnisberechnung(t: float, U: float, V: float, L1S: float, L2S: float, D: float, H: float, phi: float, delta_t: float) -> Tuple[float, float, float, float, float, float, float]:
    # Zeit des lokalen Maximums
    TD = T0 + t
    # Magnitude G
    m = sqrt(U * U + V * V)
    G = (L1S - m) / (L1S + L2S)
    # Durchmesserverhältnis
    A = (L1S - L2S) / (L1S + L2S)
    # Positionswinklel zur maximalen Verfinsterung, vom Nordpol der Sonne
    Pm = Atan2(U , V)
    if Pm < 0:
        Pm = Pm + 360.0
    # Positionswinkel Zm bezogen auf die Zenit-Richtung
    sinH = Sin(D) * Sin(phi) + Cos(D) * Cos(phi) * Cos(H)
    h = Asin(sinH)
    sinq = (Cos(phi) * Sin(H)) / Cos(h)
    q = Asin(sinq)
    Zm = Pm - q
    # Positionswinkel
    P = Atan2(U, V)
    if P < 0:
        P = P + 360
    UT = TD - delta_t / 3600
    UTh = int(UT)
    UTm0 = 60 * (UT - UTh)
    UTm = int(UTm0)
#    UTs = int(60 * (UTm0 - UTm) + 0.5)
    UTs = 60 * (UTm0 - UTm)
    return UTh, UTm, UTs, G, A, Zm, P



def bessel3(Lambda: float, phi: float, hoehe: float, delta_t: float) -> None:
    ic(Lambda, phi, hoehe)

    # Hauptprogramm
    t = 0

    # rechtwinklige geozentrische Koordinaten
    U = Atan(0.99664719 * Tan(phi))
    pSinPS = 0.99664719 * Sin(U) + hoehe / 6378140 * Sin(phi)
    pCosPS = Cos(U) + hoehe / 6378140 * Cos(phi)

    # 1) Ergebnisse zum Maximums-Zeitpunkt
    for i in range(5):  # 5 Iterationen genügen
        a, b, U, V, n2, L1S, L2S, D, H = magnitudePoswinkel(t, pSinPS, pCosPS, Lambda, delta_t)
        # Zeitpunkt der maximalen Finsternis per Iteration
        tm = -(U * a + V * b) / n2
        t = t + tm
        UTh, UTm, UTs, G, A, Zm, P = ergebnisberechnung(t, U, V, L1S, L2S, D, H, phi, delta_t)
    print("Uhrzeit des Maximums:", UTh, "h", UTm, "m", UTs, "s UT")  # dynamische Zeit
    print("Magnitude: ", int(1000 * G) / 10, "%")
    print("Verhältnis Durchmesser Mond/Sonne: ", int(1000 * A) / 1000)
    print("Positionswinkel, bezogen auf Nord: ", int(P+0.5), "Grad")
    print("Positionswinkel, bezogen auf Zenit: ", int(Zm+0.5), "Grad")

    # 2) Ergebnisse aller 10 Minuten
    # tmax = t
    # t = tmax - 3 
    # print("UT             Magnitude Nord-  Zenit- Posw.")
    # for i in range(480):  # 3 h aller Minuten
    #     a, b, U, V, n2, L1S, L2S, D, H = magnitudePoswinkel(t)
    #     UTh, UTm, UTs, G, A, Zm, P = ergebnisberechnung()
    #     if G >= 0:
    #         print("%2.0f h %02.0f m %02.0f s  %5.1f%%     %3.0f° %3.0f°" % (UTh, UTm, UTs, 100*G, P+0.5, Zm+0.5))
    #         # print(UTh,"h", UTm,"m", UTs,"s  ", int(1000 * G) / 10, "%  ", int(P+0.5),"°  ", int(Zm+0.5), "°")
    #     t = t + 1 / 60

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
    bessel3(-loc.lon.value, loc.lat.value, loc.height.value, delta_t.value)


if __name__ == "__main__":
    main()
