# Astropy-Workbench

Python scripts using astropy and astroquery

Copyright 2024 Martin Junius

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


## VizieR

### astroquery.vizier Package
VizieR Query Tool

This package is for querying the VizieR service, primarily hosted at: https://vizier.cds.unistra.fr

Note: If the access to catalogues with VizieR was helpful for your research work, the following acknowledgment would be appreciated:

This research has made use of the VizieR catalogue access tool, CDS,
Strasbourg, France.  The original description of the VizieR service was
published in A&AS 143, 23

### queryvizier.py

Query VizieR catalog

```
usage: queryvizier [-h] [-v] [-d] [--columns COLUMNS] [-n ROW_LIMIT] [-f] [-C] [-o OUTPUT] [--replace-comma REPLACE_COMMA] catalog

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
  -o OUTPUT, --output OUTPUT
                        output file
  --replace-comma REPLACE_COMMA
                        replace "," in field with REPLACE_COMMA

Version 0.2 / 2025-01-06 / Martin Junius
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
