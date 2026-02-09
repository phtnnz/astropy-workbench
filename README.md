# Astropy-Workbench

Python scripts using Astropy and friends

Copyright 2024-2026 Martin Junius

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.


## About

This is just my personal playground for working with astropy and friends,
not necessarily usable for anyone else.


## Installation

Clone the repository, create an virtual environment, and activate it (Windows command line):
```
> git clone https://github.com/phtnnz/astropy-workbench.git
> cd astropy-workbench
> python -m venv venv
> .\venv\Scripts\activate.bat
```

Install the required packages (at lot!):
```
> pip install -r requirements.txt
```


## VizieR

### astroquery.vizier Package
VizieR Query Tool

This package is for querying the VizieR service, primarily hosted at: https://vizier.cds.unistra.fr

This work makes use of the VizieR catalogue access tool, CDS,
Strasbourg, France.  The original description of the VizieR service was
published in A&AS 143, 23

### queryvizier.py

Query VizieR catalog

```
usage: queryvizier [-h] [-v] [-d] [--columns COLUMNS] [-n ROW_LIMIT] [-f] [-C] [-l] [-o OUTPUT] [--replace-comma REPLACE_COMMA] [--object OBJECT] [-m MATCH]
                   [--constellation]
                   catalog

Query VizieR catalog

positional arguments:
  catalog               catalog name

options:
  -h, --help            show this help message and exit
  -v, --verbose         verbose messages
  -d, --debug           more debug messages
  --columns COLUMNS     columns to retrieve, comma-separated, default all
  -n ROW_LIMIT, --row-limit ROW_LIMIT
                        number of rows to retrieve, default unlimited
  -f, --ra-dec-float    output RA/DEC as float degrees
  -C, --csv             CSV output
  -l, --locale          set locale for CSV output
  -o OUTPUT, --output OUTPUT
                        output file
  --replace-comma REPLACE_COMMA
                        replace "," in field with REPLACE_COMMA
  --object OBJECT       query specific object, default ""=all, no wildcards
  -m MATCH, --match MATCH
                        output rows containing regex MATCH only
  --constellation       get and output constellation as an extra column

Version 0.4 / 2025-02-05 / Martin Junius
```


## Catalogs

### PixInsight CSV catalog headers
Mandatory: id, alpha, delta  
(id=Name, alpha=RA float degrees, delta=DEC float degrees)

Common: magnitude, diameter, axisRatio, posAngle

Other: Common name, PGC, PGC2, Messier, NGC/IC, Spectral type, HD, HIP

### cat/dust-clouds.csv

https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=J/A%2bA/383/631&-out.max=50&-out.form=HTML%20Table&-out.add=_r&-out.add=_RAJ,_DEJ&-sort=_r&-oc.form=sexa

Prepared for PixInsight's AnnotateImage

```
> .\queryvizier.py -C -v --columns RAJ2000=alpha,DEJ2000=delta,Names=id,MajAxis=diameter --replace-comma=" / " -f -o cat\dust-clouds.csv 'J/A+A/383/631'
queryvizier: query catalog J/A+A/383/631
queryvizier: catalog title: Catalogue of dust clouds in the Galaxy
queryvizier:       authors: Dutra C.M.; Bica E.
queryvizier: table columns = ['RAJ2000', 'DEJ2000', 'Names', 'MajAxis']
queryvizier: output columns = ['alpha', 'delta', 'id', 'diameter']
```

### cat/CGs.csv

```
> .\queryvizier.py -C -v --columns Names,RAJ2000,DEJ2000 --replace-comma=" / " -m "CG\d+[A-Z]* DN" -l -o .\cat\CGs.csv 'J/A+A/383/631'
```

Catalog excerpt with cometary globules (CGs)

```
> .\queryvizier.py --constellation -C -o .\cat\CGs-constellation.csv -v --columns Names,RAJ2000,DEJ2000 --replace-comma=" / 
" -m "CG\d+[A-Z]* DN" -l 'J/A+A/383/631'
```

Catalog excerpt with cometary globules (CGs) and added constellation names


## References

https://arxiv.org/abs/astro-ph/0602086  
https://aa.usno.navy.mil/downloads/Circular_179.pdf

The IAU Resolutions on Astronomical Reference Systems, Time Scales, and Earth Rotation Models  
Explanation and Implementation  
by George H. Kaplan

http://www.iausofa.org/cookbooks.html  
http://www.iausofa.org/sofa_pn_c.pdf  
http://www.iausofa.org/sofa_ts_c.pdf

International Astronomical Union, Standards of Fundamental Astronomy

SOFA Tools for Earth Attitude  
Software version 18  
Document revision 1.7

SOFA Time Scale and Calendar Tools  
Software version 18  
Document revision 1.63
