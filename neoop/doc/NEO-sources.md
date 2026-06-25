# Web sources for NEO and comet data

## MPC

### MPC Database Search
https://minorplanetcenter.net/db_search

Example:
https://minorplanetcenter.net/db_search/show_object?utf8=%E2%9C%93&object_id=C%2F2026+L1

Then, observations in 80 columns format:
https://minorplanetcenter.net/tmp2/C_2026_L1.txt


### MPC Web Service:
https://minorplanetcenter.net/search_db

Example:
```
url = 'https://minorplanetcenter.net/search_db'

params = {
  "table": "observations",      # orbits or observations
  "designation": "C/2026 L1"    # provId / designation of object, comets with CPI prefix
}

r = requests.post(url, params, auth = ('mpc_ws', 'mpc!!ws'))
```


### Dates Of Last Observation Of Unusual Minor Planets
https://minorplanetcenter.net/iau/lists/LastUnusual.html  
https://www.minorplanetcenter.net/iau/lists/Customize.html

API:
https://cgi.minorplanetcenter.net/cgi-bin/customize.cgi 


### NEO Confirmation Page
https://www.minorplanetcenter.net/iau/NEO/toconfirm_tabular.html

API:
https://cgi.minorplanetcenter.net/cgi-bin/confirmeph2.cgi


### Previous NEO Confirmation Page Objects
https://www.minorplanetcenter.net/iau/NEO/ToConfirm_PrevDes.html


## JPL

### What's Observable?
https://ssd.jpl.nasa.gov/tools/sbwobs.html#/

API:
https://ssd-api.jpl.nasa.gov/doc/sbwobs.html
