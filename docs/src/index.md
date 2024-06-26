# LatLon.jl Documentation

# Overview
This module computes geo-distances between points and between sets.
Also provided are retrieve functions for geo-data from NASA via KML files.

The distance functions can be applied to the data from the geo-data
from NASA files. However, it should be noted that the distance functions 
assume a perfect sphere and do not, for instance, correct for the *polar
flattening* of most planetary objects.

**NOTE:** This module uses a constant *internal* structure, `MC`.
The values from this structure are used in several of the functions below.
In addition, some functions use values from this structure as default values
for some of the arguments.

# Exports
- Geo-Distance Functions:
    - geo\_dist
        - Computes the distance between two points on a sphere.
          The default sphere is the earth with the radius given in Kilometers.
    - geo\_midpoint
        - Computes the mid lat/lon between to lat/lon coordinates. 
    - latlon\_set\_dist
        - Computes the distance between two geo-sets on a sphere.
- Geo-Extraction Functions:
    - center\_latlon\_from\_NASA\_xml\_file
        - Retrieve lat/lon data from a NASA XML "center" file.
    - upath\_latlon\_from\_NASA\_xml\_file
        - Retrieve lat/lon data from a NASA XML "upath" file.

```@meta
CurrentModule = LatLon
```

## Geo-Distance Functions

```@docs
geo_dist
```

```@docs
geo_midpoint
```

```@docs
latlon_set_dist
```

## Geo-Extraction Functions

```@docs
center_latlon_from_NASA_xml_file
```

```@docs
upath_latlon_from_NASA_xml_file
```

## Index

```@index
```

