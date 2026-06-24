# NEO / NEOCP Observation Planning

## NEO Obs Planner

```
usage: neo-obs-planner [-h] [-v] [--verbose-ephem] [-d] [-l LOCATION] [-f FILE] [-s START] [-e END] [-o OUTPUT] [-C] [-P] [--clear] [-M MAG_LIMIT]
                       [--neocp-mag-limit NEOCP_MAG_LIMIT] [--sbwobs-mag-limit SBWOBS_MAG_LIMIT] [-m MIN_ALT] [--neocp] [--sbwobs] [--asteroids] [--neo]
                       [--pha] [--comets] [-p PREFIX] [--force FORCE]
                       [object ...]

NEOCP/NEO observation planner

positional arguments:
  object                object name

options:
  -h, --help            show this help message and exit
  -v, --verbose         verbose messages
  --verbose-ephem       verbose ephemerides
  -d, --debug           more debug messages
  -l LOCATION, --location LOCATION
                        coordinates, named location or MPC station code, default M49
  -f FILE, --file FILE  read list of objects from file
  -s START, --start START
                        start time (UTC) (default naut. dusk)
  -e END, --end END     end time (UTC) (default naut. dawn)
  -o OUTPUT, --output OUTPUT
                        write CSV to OUTPUT file
  -C, --csv             use CSV output format
  -P, --plot            create altitude and sky plot with objects
  --clear               clear MPC cache
  -M MAG_LIMIT, --mag-limit MAG_LIMIT
                        override *mag_limits from config (see below)
  --neocp-mag-limit NEOCP_MAG_LIMIT
                        override neocp_mag_limit from config (20.5)
  --sbwobs-mag-limit SBWOBS_MAG_LIMIT
                        override sbwobs_mag_limit from config (19.5)
  -m MIN_ALT, --min-alt MIN_ALT
                        override min_alt/elev_min from config
  --neocp               observable NEOCP objects
  --sbwobs              observable objects from JPL WOBS service
  --asteroids           sbwobs: get asteroids default=a
  --neo                 sbwobs: get NEOs default=neo
  --pha                 sbwobs: get PHAs
  --comets              sbwobs: get comets (overrides asteroids options)
  -p PREFIX, --prefix PREFIX
                        prefix for planner data, default 20260624
  --force FORCE         skip checks for FORCE objects, include in observation plan

Version 2.1 / 2026-06-20 / Martin Junius
```

Retrieve lists and ephemerides for upcoming night
```
neo-obs-planner.py -v --neocp --sbwobs 
```
```--neocp``` = get NEOCP objects, ```--sbwobs``` get observable "unusual" NEO objects

Use options ```-CP``` to create CSV plan output and graphic plot
```
neo-obs-planner.py -v --neocp --sbwobs -CP
```

Use option ```-M``` to limit magnitude for *both* NEOCP and NEO around full moon phase
```
neo-obs-planner.py -v --neocp --sbwobs -M 19.5
```

Output saved to ./neo-obs-data/


## N.I.N.A.

Create sequence:
```
nina-create-sequence2.py -v --setting remote3-neo .\neo-obs-data\YYYYMMDD-neo-obs-plan.csv
```

## Example

```
.\neo-obs-planner.py -v --neocp --sbwobs -CP -M 19.5                                        
neo-obs-planner: download ephemerides from https://cgi.minorplanetcenter.net/cgi-bin/confirmeph2.cgi
neo-obs-planner: download NEOCP list from https://minorplanetcenter.net/iau/NEO/neocp.txt
neo-obs-planner: download PCCP list from https://minorplanetcenter.net/iau/NEO/pccp.txt
neo-obs-planner: NEOCP objects (1): A11E1GK
neo-obs-planner: query https://ssd-api.jpl.nasa.gov/sbwobs.api
neo-obs-planner: WOBS objects (39): 2001 FF90, 2005 NE21, 2006 KD40, 2007 QA2, 2007 TO74, 2011 LA19, 2011 LJ1, 2011 SZ21, 2012 LE11, 2016 BN1, 2017 KN4, 2018 JB3, 2018 LH5, 2018 NC15, 2019 FH, 2019 GU2, 2020 TV8, 2022 LN3, 2022 PM1, 2022 UC6, 2025 WC4, 2026 BS12, 2026 HL1, 2026 KC3, 2026 KX, 2026 LS, 2026 MC2, 2026 MD, 2026 MD2, 2026 MJ1, 2026 MM1, 2026 MN, 2026 MN2, 2026 MP1, 2026 MW, 2026 MW1, 2026 MX2, 2026 MY, 2026 MZ1
neo-obs-planner: query https://cgi.minorplanetcenter.net/cgi-bin/customize.cgi
neo-obs-planner: DLU objects (17): 2010 DD20, 2017 BH30, 2019 UN13, 2019 YB, 2022 DX29, 2022 ET, 2022 WT27, 2024 MW, 2024 UC21, 2025 CZ1, 2025 DW, 2026 BS12, 2026 MN, 2026 MP1, 2026 MW1, 2026 MX2, 2026 MZ1
neo-obs-planner: 2026 BS12: last obs 2026-04-12 00:00:00.000 too old, not included
neo-obs-planner: 2025 DW: last obs 2025-02-20 00:00:00.000 too old, not included
neo-obs-planner: 2025 CZ1: last obs 2025-02-07 00:00:00.000 too old, not included
neo-obs-planner: 2024 UC21: last obs 2024-10-28 00:00:00.000 too old, not included
neo-obs-planner: 2024 MW: last obs 2024-06-28 00:00:00.000 too old, not included
neo-obs-planner: 2022 WT27: last obs 2022-11-23 00:00:00.000 too old, not included
neo-obs-planner: 2022 ET: last obs 2022-03-07 00:00:00.000 too old, not included
neo-obs-planner: 2022 DX29: last obs 2022-02-26 00:00:00.000 too old, not included
neo-obs-planner: 2019 YB: last obs 2019-12-17 00:00:00.000 too old, not included
neo-obs-planner: 2019 UN13: last obs 2019-10-31 00:00:00.000 too old, not included
neo-obs-planner: 2017 BH30: last obs 2017-01-30 00:00:00.000 too old, not included
neo-obs-planner: 2010 DD20: last obs 2010-02-18 00:00:00.000 too old, not included
neo-obs-planner: WOBS & DLU objects (5): 2026 MN, 2026 MP1, 2026 MW1, 2026 MX2, 2026 MZ1
neo-obs-planner: ------------------------------------------------------------
neo-obs-planner: Type   Designation Rise   Trans  Set     Vmag  U  Last Obs
neo-obs-planner: ------------------------------------------------------------
neo-obs-planner: NEO    2026 MN     19:00  23:09  03:18   19.3  6  2026-06-23
neo-obs-planner: NEO    2026 MP1    22:20  02:42  07:04*  18.3  8  2026-06-23
neo-obs-planner: NEO    2026 MW1    22:10  02:41  07:13*  17.1  8  2026-06-23
neo-obs-planner: NEO    2026 MX2    14:42* 19:34  00:26   19.5  8  2026-06-24
neo-obs-planner: NEO    2026 MZ1    17:11* 21:37  02:03   18.1  7  2026-06-24
neo-obs-planner: ------------------------------------------------------------
neo-obs-planner: location lon=16d21m42.2s lat=-23d14m11.6s height=1853. m code M49
neo-obs-planner: nautical twilight 2026-06-24 17:10:58.809 / 2026-06-25 04:43:13.825 (UTC)
neo-obs-planner: already got ephemeris for A11E1GK
neo-obs-planner: 2026 MN ephemeris from MPC
neo-obs-planner: 2026 MP1 ephemeris from MPC
neo-obs-planner: 2026 MW1 ephemeris from MPC
neo-obs-planner: 2026 MX2 ephemeris from MPC
neo-obs-planner: 2026 MZ1 ephemeris from MPC
neo-obs-planner: original object sequence: A11E1GK, 2026 MN, 2026 MP1, 2026 MW1, 2026 MX2, 2026 MZ1
neo-obs-planner: NEOCP A11E1GK 2026-06-24 22:30:00.000
neo-obs-planner: NEO 2026 MN 2026-06-24 21:00:00.000
neo-obs-planner: NEO 2026 MP1 2026-06-25 00:30:00.000
neo-obs-planner: NEO 2026 MW1 2026-06-25 00:00:00.000
neo-obs-planner: NEO 2026 MX2 2026-06-24 17:30:00.000
neo-obs-planner: NEO 2026 MZ1 2026-06-24 19:00:00.000
neo-obs-planner: sorted object sequence: 2026 MX2, 2026 MZ1, 2026 MN, A11E1GK, 2026 MW1, 2026 MP1
neo-obs-planner: forced objects: 
neo-obs-planner: 
neo-obs-planner: obs-planner-1 2026-06-24 18:14:05 UTC
neo-obs-planner: -------------------------------------------------------------------------------------------------------------------
neo-obs-planner:               Score       Mag #Obs      Arc NotSeen  Time start ephemeris/ end ephemeris                 Max motion
neo-obs-planner:        /Uncertainty                                  Time before         / after meridian             Moon distance
neo-obs-planner:                                                      Time start exposure / end exposure
neo-obs-planner:                                                      # x Exp = total exposure time
neo-obs-planner:                                                      RA, DEC, Alt, Az
neo-obs-planner: -------------------------------------------------------------------------------------------------------------------
neo-obs-planner: 2026 MX2  NEO     8  19.3 mag                 0.8 d  2026-06-24 17:30:00 / 2026-06-25 00:00:00     6.7 arcsec / min
neo-obs-planner: SKIPPED: moon distance 19 deg < 50 deg
neo-obs-planner: -------------------------------------------------------------------------------------------------------------------
neo-obs-planner: 2026 MZ1  NEO     7  18.0 mag                 0.8 d  2026-06-24 17:30:00 / 2026-06-25 02:00:00    19.3 arcsec / min
neo-obs-planner: SKIPPED: moon distance 41 deg < 50 deg
neo-obs-planner: -------------------------------------------------------------------------------------------------------------------
neo-obs-planner: 2026 MN   NEO     6  19.2 mag                 1.8 d  2026-06-24 19:30:00 / 2026-06-25 03:00:00     5.2 arcsec / min
neo-obs-planner:                                                      2026-06-24 23:00:00 / 2026-06-24 23:30:00               64 deg
neo-obs-planner:                                                      2026-06-24 21:00:00 / 2026-06-24 21:21:55
neo-obs-planner:                                                      30 x 30 s = 15.0 min (120%) / total 21.9 min
neo-obs-planner:                                                      RA 18.4283 hourangle, DEC -2.1150 deg, Alt 46 deg, Az 67 deg
neo-obs-planner: -------------------------------------------------------------------------------------------------------------------
neo-obs-planner: A11E1GK   NEOCP 100  18.8 mag    6   0.35 d   0.2 d  2026-06-24 18:30:00 / 2026-06-25 04:00:00     2.3 arcsec / min
neo-obs-planner:                                                      2026-06-24 23:00:00 / 2026-06-24 23:30:00               62 deg
neo-obs-planner:                                                      2026-06-24 23:30:00 / 2026-06-25 00:06:55
neo-obs-planner:                                                      30 x 60 s = 30.0 min (333%) / total 36.9 min
neo-obs-planner:                                                      RA 18.5715 hourangle, DEC -61.7631 deg, Alt 51 deg, Az 177 deg
neo-obs-planner: -------------------------------------------------------------------------------------------------------------------
neo-obs-planner: 2026 MW1  NEO     8  17.1 mag                 1.8 d  2026-06-24 22:30:00 / 2026-06-25 04:30:00    35.8 arcsec / min
neo-obs-planner:                                                      2026-06-25 02:30:00 / 2026-06-25 03:00:00              107 deg
neo-obs-planner:                                                      2026-06-25 00:06:55 / 2026-06-25 00:16:10
neo-obs-planner:                                                      53 x  2 s = 1.8 min (100%) / total 9.3 min
neo-obs-planner:                                                      RA 21.9478 hourangle, DEC -16.2922 deg, Alt 52 deg, Az 87 deg
neo-obs-planner: -------------------------------------------------------------------------------------------------------------------
neo-obs-planner: 2026 MP1  NEO     8  18.3 mag                 1.8 d  2026-06-24 22:30:00 / 2026-06-25 04:30:00    27.0 arcsec / min
neo-obs-planner:                                                      2026-06-25 02:30:00 / 2026-06-25 03:00:00              111 deg
neo-obs-planner:                                                      2026-06-25 00:30:00 / 2026-06-25 00:43:06
neo-obs-planner:                                                      64 x  5 s = 5.3 min (100%) / total 13.1 min
neo-obs-planner:                                                      RA 22.0174 hourangle, DEC -9.4758 deg, Alt 49 deg, Az 77 deg
neo-obs-planner: -------------------------------------------------------------------------------------------------------------------
neo-obs-planner: 4 object(s) planned: 2026 MN, A11E1GK, 2026 MW1, 2026 MP1
neo-obs-planner: 2 object(s) skipped: 2026 MX2, 2026 MZ1
neo-obs-planner: planned objects for nina-create-sequence2: neo-obs-data\20260624-neo-obs-plan.csv
neo-obs-planner: altitude and sky plot for objects: neo-obs-data\20260624-neo-obs-plot.png
```

```
.\nina-create-sequence2.py -v --setting remote3-neo .\neo-obs-data\20260624-neo-obs-plan.csv
nina-create-sequence2: processing target template D:/Users/mj/Documents/N.I.N.A/Templates/NINA-Templates-IAS-Common/Target NEO.template.json
nina-create-sequence2: processing sequence template D:/Users/mj/Documents/N.I.N.A/Templates/NINA-Templates-IAS-Common/Base Remote3 NAUTICAL.json
nina-create-sequence2: target format (0=target, 1=date, 2=seq, 3=number) {1} {2:03d} {0} (n{3:03d})
nina-create-sequence2: output format (1=date) NEO-{1}.json
nina-create-sequence2: add target items to container '', empty=target area
nina-create-sequence2: timezone Africa/Windhoek
nina-create-sequence2: subdir (1=date) _asteroids_{1}
nina-create-sequence2: autofocus first target only = False
nina-create-sequence2: destination directory D:\Users\mj\Documents\N.I.N.A
nina-create-sequence2: output file NEO-2026-06-24.json
nina-create-sequence2: NINATarget(process_data): name = Target NEO
nina-create-sequence2: NINASequence(process_data): name = Base Remote3 NAUTICAL
nina-create-sequence2: processing CSV file .\neo-obs-data\20260624-neo-obs-plan.csv
nina-create-sequence2: ------------------------------------------------------------------
nina-create-sequence2: #001 2026-06-24 001 2026 MN (n030)    18h25m41.800s -02d06m54.000s
nina-create-sequence2: UT=2026-06-24 21:00:00+00:00 / local 2026-06-24 23:00:00+02:00
nina-create-sequence2: 30x30.0s filter=L
nina-create-sequence2: ------------------------------------------------------------------
nina-create-sequence2: #002 2026-06-24 002 A11E1GK (n030)    18h34m17.500s -61d45m47.000s
nina-create-sequence2: UT=2026-06-24 23:30:00+00:00 / local 2026-06-25 01:30:00+02:00
nina-create-sequence2: 30x60.0s filter=L
nina-create-sequence2: ------------------------------------------------------------------
nina-create-sequence2: #003 2026-06-24 003 2026 MW1 (n053)   21h56m52.000s -16d17m32.000s
nina-create-sequence2: UT=2026-06-25 00:06:55+00:00 / local 2026-06-25 02:06:55+02:00
nina-create-sequence2: 53x2.0s filter=L
nina-create-sequence2: ------------------------------------------------------------------
nina-create-sequence2: #004 2026-06-24 004 2026 MP1 (n064)   22h01m02.800s -09d28m33.000s
nina-create-sequence2: UT=2026-06-25 00:30:00+00:00 / local 2026-06-25 02:30:00+02:00
nina-create-sequence2: 64x5.0s filter=L
nina-create-sequence2: ------------------------------------------------------------------
nina-create-sequence2: writing JSON sequence D:\Users\mj\Documents\N.I.N.A\NEO-2026-06-24.json
```

```
.\nina-create-sequence2.py -l --setting remote3-neo .\neo-obs-data\20260624-neo-obs-plan.csv
2026-06-24 001 2026 MN (n030)   NEO   18 25 41.80 -02 06 54.0  19.2
2026-06-24 002 A11E1GK (n030)   NEOCP 18 34 17.50 -61 45 47.0  18.8
2026-06-24 003 2026 MW1 (n053)  NEO   21 56 52.00 -16 17 32.0  17.1
2026-06-24 004 2026 MP1 (n064)  NEO   22 01 02.80 -09 28 33.0  18.3
```
