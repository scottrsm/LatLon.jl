# LatLon.jl

[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://scottrsm.github.io/LatLon.jl/dev/)

Functions are provided to manipulate Lat/Lon geo-coordinates including functions
to work with NASA KML files.

A non-Haversine formula is used to compute the distance between points on the sphere:
the angle between the two points is recovered from both the cosine (dot product) and the
sine (cross product) of the angle with a two-argument arc-tangent, which is accurate
for the whole sphere, including very close and near-antipodal points.

A Jupyter notebook is provided which examines recent solar eclipses in the US.

## Documentation
- HTML (latest, built from `main`): https://scottrsm.github.io/LatLon.jl/dev/
- Markdown source: [docs/src/index.md](docs/src/index.md)
