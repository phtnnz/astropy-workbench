# NEO / NEOCP Observation Planning

## NEOCP

Update lists and ephemeris for upcoming night:
```
neocp -v -U
```

Plot and output CSV plan:
```
neocp -v -P -C
```

Output in ./neocp-data/


## NEOs

Find observable "unusual" NEOs:
```
sbwobs.py -v -o .\tmp\20260203-neos.list
```

Plan observations:
```
neo-obs-planner.py -v -f .\tmp\20260203-neos.list
```

... with altitude/sky plot and output to CSV plan:
```
neo-obs-planner.py -v -f .\tmp\20260203-neos.list -P -C -o .\tmp\20260203-neos.csv
```

## N.I.N.A.

Create sequence:
```
nina-create-sequence2.py -v .\neocp-data\20260208-neocp-plan.csv --setting remote3-neo
```

## Example

```
> .\neocp.py -v -C -P -U
neocp: location: M49 lon=16d21m42.2s lat=-23d14m11.6s height=1853. m
neocp: download ephemerides from https://cgi.minorplanetcenter.net/cgi-bin/confirmeph2.cgi
neocp: download NEOCP list from https://minorplanetcenter.net/iau/NEO/neocp.txt
neocp: download PCCP list from https://minorplanetcenter.net/iau/NEO/pccp.txt
neocp: processing neocp-data\20260208-NEOCP-ephemerides.html
neocp: skipping NEOCP obj='C19RNP5' (empty)
neocp: processing neocp-data\20260208-NEOCP-list.txt
neocp: processing neocp-data\20260208-PCCP-list.txt
neocp: planning objects: CE5W292 P22liAW P22lpqR ZTF10BK ZTF10BL ZTF10BN CE5ZET2 ST26B21 CE60272 TF26B20 C45TCR1 ST26B19 ST26B20 P22lbY0 ZTF10BM C45TCY1 CE519G2
[...]
neocp: -----------------------------------------------------------------------------------------------------------------------
neocp:              Score      MagV #Obs      Arc NotSeen  Time start ephemeris   / end ephemeris                  Max motion
neocp:                                                     Time before            / after meridian              Moon distance
neocp:                                                     Time start exposure    / end exposure
neocp:                                                     # x Exp = total exposure time
neocp:                                                     RA, DEC, Alt, Az
neocp: -----------------------------------------------------------------------------------------------------------------------
neocp: CE5W292  NEOCP  79  19.3 mag   38  18.39 d   0.1 d  2026-02-08 18:30:00.000/2026-02-08 22:00:00.000    5.6 arcsec / min
neocp:                                                     2026-02-08 19:00:00.000/2026-02-08 19:30:00.000             135 deg
neocp:                                                     2026-02-08 18:30:00.000/2026-02-08 18:49:25.500
neocp:                                                     37 x 20 s = 12.3 min (100%) / total 19.4 min
neocp:                                                     RA 5.7147 hourangle, DEC 24.6986 deg, Alt 41 deg, Az 16 deg
neocp: -----------------------------------------------------------------------------------------------------------------------
neocp: P22liAW  NEOCP  10  18.9 mag   38  21.40 d   0.1 d  2026-02-08 18:30:00.000/2026-02-09 00:00:00.000    4.5 arcsec / min
neocp:                                                     2026-02-08 20:00:00.000/2026-02-08 20:30:00.000             115 deg
neocp:                                                     2026-02-08 18:49:25.500/2026-02-08 19:05:34.000
neocp:                                                     19 x 30 s = 9.5 min (100%) / total 16.1 min
neocp:                                                     RA 6.4575 hourangle, DEC -8.0153 deg, Alt 62 deg, Az 61 deg
neocp: -----------------------------------------------------------------------------------------------------------------------
neocp: P22lpqR  NEOCP  81  20.3 mag   23  24.83 d   0.6 d  2026-02-08 18:30:00.000/2026-02-09 00:30:00.000    8.0 arcsec / min
neocp:                                                     2026-02-08 20:30:00.000/2026-02-08 21:00:00.000             109 deg
neocp:                                                     2026-02-08 19:05:34.000/2026-02-08 19:33:14.000
neocp:                                                     60 x 20 s = 20.0 min (60%) / total 27.7 min
neocp:                                                     RA 7.2844 hourangle, DEC 3.9089 deg, Alt 51 deg, Az 50 deg
neocp: -----------------------------------------------------------------------------------------------------------------------
neocp: ZTF10BK  NEOCP 100  17.6 mag    6   0.05 d   0.4 d  2026-02-08 18:30:00.000/2026-02-09 01:30:00.000  211.1 arcsec / min
neocp: SKIPPED: object too fast (>79.8 arcsec / min)
neocp: -----------------------------------------------------------------------------------------------------------------------
neocp: ZTF10BL  NEOCP 100  17.1 mag    4   0.01 d   0.5 d  2026-02-08 19:30:00.000/2026-02-09 00:30:00.000   37.4 arcsec / min
neocp: SKIPPED: arc 0.01 d too small (< 0.05 d)
neocp: -----------------------------------------------------------------------------------------------------------------------
neocp: ZTF10BN  NEOCP 100  19.5 mag    4   0.04 d   0.4 d  2026-02-08 19:00:00.000/2026-02-09 02:00:00.000    6.6 arcsec / min
neocp: SKIPPED: arc 0.04 d too small (< 0.05 d)
neocp: -----------------------------------------------------------------------------------------------------------------------
neocp: CE5ZET2  NEOCP 100  20.4 mag   11   0.14 d   0.3 d  2026-02-08 19:30:00.000/2026-02-09 01:30:00.000   14.0 arcsec / min
neocp: SKIPPED: only 27% of required total exposure time (< 30%)
neocp: -----------------------------------------------------------------------------------------------------------------------
neocp: ST26B21  NEOCP 100  19.1 mag    4   0.03 d   0.3 d  2026-02-08 19:00:00.000/2026-02-09 03:00:00.000  174.9 arcsec / min
neocp: SKIPPED: object too fast (>79.8 arcsec / min)
neocp: -----------------------------------------------------------------------------------------------------------------------
neocp: CE60272  NEOCP  98  20.3 mag   15   0.27 d   0.2 d  2026-02-08 20:00:00.000/2026-02-09 02:30:00.000   13.3 arcsec / min
neocp:                                                     2026-02-08 23:00:00.000/2026-02-08 23:30:00.000              80 deg
neocp:                                                     2026-02-08 22:30:00.000/2026-02-08 22:47:40.000
neocp:                                                     60 x 10 s = 10.0 min (30%) / total 17.7 min
neocp:                                                     RA 9.7505 hourangle, DEC 14.7217 deg, Alt 47 deg, Az 30 deg
neocp: -----------------------------------------------------------------------------------------------------------------------
neocp: TF26B20  NEOCP 100  19.3 mag    5   0.43 d   0.9 d  2026-02-08 20:00:00.000/2026-02-09 03:00:00.000    9.5 arcsec / min
neocp:                                                     2026-02-08 23:30:00.000/2026-02-09 00:00:00.000              73 deg
neocp:                                                     2026-02-08 22:47:40.000/2026-02-08 23:08:24.500
neocp:                                                     53 x 15 s = 13.2 min (100%) / total 20.7 min
neocp:                                                     RA 9.9854 hourangle, DEC 6.7336 deg, Alt 56 deg, Az 31 deg
neocp: -----------------------------------------------------------------------------------------------------------------------
neocp: C45TCR1  NEOCP  98  18.2 mag   16   0.31 d   0.0 d  2026-02-08 20:30:00.000/2026-02-09 03:30:00.000  498.2 arcsec / min
neocp: SKIPPED: object too fast (>79.8 arcsec / min)
neocp: -----------------------------------------------------------------------------------------------------------------------
neocp: ST26B19  NEOCP  98  20.1 mag   10   0.17 d   0.3 d  2026-02-08 22:00:00.000/2026-02-09 03:00:00.000   19.4 arcsec / min
neocp: SKIPPED: only 15% of required total exposure time (< 30%)
neocp: -----------------------------------------------------------------------------------------------------------------------
neocp: ST26B20  NEOCP 100  19.5 mag    4   0.04 d   0.4 d  2026-02-08 21:30:00.000/2026-02-09 03:30:00.000   91.4 arcsec / min
neocp: SKIPPED: object too fast (>79.8 arcsec / min)
neocp: -----------------------------------------------------------------------------------------------------------------------
neocp: P22lbY0  PCCP   35  19.8 mag   62  10.80 d   0.5 d  2026-02-08 22:30:00.000/2026-02-09 03:30:00.000    0.4 arcsec / min
neocp:                                                     2026-02-09 01:30:00.000/2026-02-09 02:00:00.000              54 deg
neocp:                                                     2026-02-08 23:08:24.500/2026-02-08 23:36:06.000
neocp:                                                     21 x 60 s = 21.0 min (100%) / total 27.7 min
neocp:                                                     RA 12.1337 hourangle, DEC 18.3575 deg, Alt 32 deg, Az 48 deg
neocp: -----------------------------------------------------------------------------------------------------------------------
neocp: ZTF10BM  NEOCP 100  18.0 mag    4   0.00 d   0.4 d  2026-02-08 21:00:00.000/2026-02-09 03:30:00.000  107.4 arcsec / min
neocp: SKIPPED: object too fast (>79.8 arcsec / min)
neocp: -----------------------------------------------------------------------------------------------------------------------
neocp: C45TCY1  NEOCP  96  17.6 mag    9   0.31 d   0.1 d  2026-02-09 00:30:00.000/2026-02-09 03:30:00.000   64.2 arcsec / min
neocp:                                                     2026-02-09 02:00:00.000/2026-02-09 02:30:00.000              63 deg
neocp:                                                     2026-02-09 00:30:00.000/2026-02-09 00:39:40.000
neocp:                                                     60 x  2 s = 2.0 min (41%) / total 9.7 min
neocp:                                                     RA 12.7024 hourangle, DEC 35.3364 deg, Alt 26 deg, Az 25 deg
neocp: -----------------------------------------------------------------------------------------------------------------------
neocp: CE519G2  NEOCP 100  19.5 mag    8   0.03 d   3.2 d  2026-02-09 01:00:00.000/2026-02-09 03:30:00.000   30.9 arcsec / min
neocp: SKIPPED: arc 0.03 d too small (< 0.05 d)
neocp: -----------------------------------------------------------------------------------------------------------------------
neocp: 7 object(s) planned: CE5W292 P22liAW P22lpqR CE60272 TF26B20 P22lbY0 C45TCY1
neocp: planned objects for nina-create-sequence2: neocp-data\20260208-neocp-plan.csv
neocp: altitude and sky plot for objects
```

```
> .\nina-create-sequence2.py -v .\neocp-data\20260208-neocp-plan.csv --setting remote3-neo
nina-create-sequence2: processing target template D:/Users/mj/Documents/N.I.N.A/Templates/NINA-Templates-IAS-Common/Target NEO (Discord).template.json
nina-create-sequence2: processing sequence template D:/Users/mj/Documents/N.I.N.A/Templates/NINA-Templates-IAS-Common/Base Remote3 NAUTICAL.json
nina-create-sequence2: target format (0=target, 1=date, 2=seq, 3=number) {1} {2:03d} {0} (n{3:03d})
nina-create-sequence2: output format (1=date) NEO-{1}.json
nina-create-sequence2: add target items to container '', empty=target area
nina-create-sequence2: timezone Africa/Windhoek
nina-create-sequence2: subdir (1=date) _asteroids_{1}
nina-create-sequence2: autofocus first target only = False
nina-create-sequence2: destination directory D:\Users\mj\Documents\N.I.N.A
nina-create-sequence2: output file NEO-2026-02-08.json
nina-create-sequence2: NINATarget(process_data): name = Target NEO (Discord)
nina-create-sequence2: NINASequence(process_data): name = Base Remote3 NAUTICAL
nina-create-sequence2: processing CSV file .\neocp-data\20260208-neocp-plan.csv
nina-create-sequence2: ------------------------------------------------------------------
nina-create-sequence2: #001 2026-02-08 001 CE5W292 (n037)    05h42m53.000s +24d41m55.000s
nina-create-sequence2: UT=2026-02-08 18:30:00+00:00 / local 2026-02-08 20:30:00+02:00
nina-create-sequence2: 37x20.0s filter=L
nina-create-sequence2: ------------------------------------------------------------------
nina-create-sequence2: #002 2026-02-08 002 P22liAW (n019)    06h27m26.900s -08d00m55.000s
nina-create-sequence2: UT=2026-02-08 18:49:25+00:00 / local 2026-02-08 20:49:25+02:00
nina-create-sequence2: 19x30.0s filter=L
nina-create-sequence2: ------------------------------------------------------------------
nina-create-sequence2: #003 2026-02-08 003 P22lpqR (n060)    07h17m03.800s +03d54m32.000s
nina-create-sequence2: UT=2026-02-08 19:05:34+00:00 / local 2026-02-08 21:05:34+02:00
nina-create-sequence2: 60x20.0s filter=L
nina-create-sequence2: ------------------------------------------------------------------
nina-create-sequence2: #004 2026-02-08 004 CE60272 (n060)    09h45m01.700s +14d43m18.000s
nina-create-sequence2: UT=2026-02-08 22:30:00+00:00 / local 2026-02-09 00:30:00+02:00
nina-create-sequence2: 60x10.0s filter=L
nina-create-sequence2: ------------------------------------------------------------------
nina-create-sequence2: #005 2026-02-08 005 TF26B20 (n053)    09h59m07.400s +06d44m01.000s
nina-create-sequence2: UT=2026-02-08 22:47:40+00:00 / local 2026-02-09 00:47:40+02:00
nina-create-sequence2: 53x15.0s filter=L
nina-create-sequence2: ------------------------------------------------------------------
nina-create-sequence2: #006 2026-02-08 006 P22lbY0 (n021)    12h08m01.400s +18d21m27.000s
nina-create-sequence2: UT=2026-02-08 23:08:24+00:00 / local 2026-02-09 01:08:24+02:00
nina-create-sequence2: 21x60.0s filter=L
nina-create-sequence2: ------------------------------------------------------------------
nina-create-sequence2: #007 2026-02-08 007 C45TCY1 (n060)    12h42m08.800s +35d20m11.000s
nina-create-sequence2: UT=2026-02-09 00:30:00+00:00 / local 2026-02-09 02:30:00+02:00
nina-create-sequence2: 60x2.0s filter=L
nina-create-sequence2: ------------------------------------------------------------------
nina-create-sequence2: writing JSON sequence D:\Users\mj\Documents\N.I.N.A\NEO-2026-02-08.json
'''
