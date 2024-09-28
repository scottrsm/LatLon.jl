# LatLon.jl
Functions are provided to manipulate Lat/Lon geo-coordinates including functions
to work with NASA KML files.

A non-Haversine formula is used to compute the distance between points on the sphere
which is faster than the standard Haversine formula. The new formula 
also works for points on the whole sphere, not just points in the same region.

A Jupyter notebook is provided which examines recent solar eclipses in the US.


