# LatLon.jl Documentation

# Overview
This module computes geo-distances between points and sets.
In addition, can also retrieve geo-data from NASA XML files.
The distance functions can be applied to the data from the geo-data
from NASA files. However, it should be noted that the distance functions 
assume a perfect sphere and do not, for instance, correct for the *polar
flattening* of most planetary objects.

# Exports
- Distance Functions:
    - geo\_dist
        - Computes the distance between two points on a sphere.
          The default sphere is the earth with distance in Kilometers.
    - geo\_midpoint
        - Computes the mid lat/lon between to lat/lon coordinates. 
    - latlon\_set\_dist
        - Computes the distance between two geo-sets on a sphere.
- Geo-extraction Functions:
    - center\_latlon\_from\_NASA\_xml\_file
        - Retrieve lat/lon data from a NASA XML "center" file.
    - upath\_latlon\_from\_NASA\_xml\_file
        - Retrieve lat/lon data from a NASA XML "upath" file.

```@meta
CurrentModule = LatLon
```

## Geo-distance Functions

```@docs
geo_dist
```

```@docs
geo_midpoint
```

```@docs
latlon_set_dist
```

## Geo-extraction Functions

```@docs
center_latlon_from_NASA_xml_file
```

```@docs
upath_latlon_from_NASA_xml_file
```

## Index

```@index
```

