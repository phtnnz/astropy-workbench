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
