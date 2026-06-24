# NEO / NEOCP Observation Planning

## NEO Obs Planner

NEO-OBS-PLANNER USAGE ...

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
nina-create-sequence2.py -v --setting remote3-neo .\neo-obs-data\20260208-neo-obs-plan.csv
```

## Example

```
'''
