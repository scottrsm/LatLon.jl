var documenterSearchIndex = {"docs":
[{"location":"#LatLon.jl-Documentation","page":"LatLon.jl Documentation","title":"LatLon.jl Documentation","text":"","category":"section"},{"location":"#Overview","page":"LatLon.jl Documentation","title":"Overview","text":"","category":"section"},{"location":"","page":"LatLon.jl Documentation","title":"LatLon.jl Documentation","text":"This module computes geo-distances between points and between sets. Also provided are retrieve functions for geo-data from NASA via KML files.","category":"page"},{"location":"","page":"LatLon.jl Documentation","title":"LatLon.jl Documentation","text":"The distance functions can be applied to the data from the geo-data from NASA files. However, it should be noted that the distance functions  assume a perfect sphere and do not, for instance, correct for the polar flattening of most planetary objects.","category":"page"},{"location":"","page":"LatLon.jl Documentation","title":"LatLon.jl Documentation","text":"NOTE: This module uses a constant internal structure, MC. The values from this structure are used in several of the functions below. In addition, some functions use values from this structure as default values for some of the arguments.","category":"page"},{"location":"#Exports","page":"LatLon.jl Documentation","title":"Exports","text":"","category":"section"},{"location":"","page":"LatLon.jl Documentation","title":"LatLon.jl Documentation","text":"Geo-Distance Functions:\ngeo_dist\nComputes the distance between two points on a sphere. The default sphere is the earth with the radius given in Kilometers.\ngeo_midpoint\nComputes the mid lat/lon between to lat/lon coordinates. \nlatlon_set_dist\nComputes the distance between two geo-sets on a sphere.\nGeo-Extraction Functions:\ncenter_latlon_from_NASA_xml_file\nRetrieve lat/lon data from a NASA XML \"center\" file.\nupath_latlon_from_NASA_xml_file\nRetrieve lat/lon data from a NASA XML \"upath\" file.","category":"page"},{"location":"","page":"LatLon.jl Documentation","title":"LatLon.jl Documentation","text":"CurrentModule = LatLon","category":"page"},{"location":"#Geo-Distance-Functions","page":"LatLon.jl Documentation","title":"Geo-Distance Functions","text":"","category":"section"},{"location":"","page":"LatLon.jl Documentation","title":"LatLon.jl Documentation","text":"geo_dist","category":"page"},{"location":"#LatLon.geo_dist","page":"LatLon.jl Documentation","title":"LatLon.geo_dist","text":"geo_dist(coord1::Vector{Float64}  , \n         coord2::Vector{Float64}  ,\n         radius=MC.radius::Float64 )\n\nComputes the distance between two points on a sphere represented as  two lon/lat vectors in signed degrees. That is, north latitude is  positive, south latitude is negative, while east longitude is positive  and west longitude is negative. The radius, R, defaults to the Earth's radius in Kilometers.\n\nDetails\n\nDot product, dp, of the geo-positions on the unit sphere gives the cosine of    the angle between the two points (on the \"great circle\").   This is the information contained in coord1 and coord2.   Once the cosine is computed we can easily retrieve the angle (in Radians) and    then it is easy to find the distance between the two points    – the length of the arc connecting them on a \"great circle\".\n\nThe Cartesian points of the lon/lat on the unit sphere are:\nbf v_1 = ( cos(theta_1)cos(phi_1) sin(theta_1)cos(phi_1) sin(phi_1) )\nbf v_2 = ( cos(theta_2)cos(phi_2) sin(theta_2)cos(phi_2) sin(phi_2) )\nThe dot product becomes:\ndp = bf v_1 cdot bf v_2 = cos(phi_1)cos(phi_2) left( cos(theta_1)cos(theta_2) + sin(theta_1)sin(theta_2) right) + sin(phi_1)sin(phi_2)\nSimplifying\ndp = cos(phi_1)cos(phi_2) cos(theta_1 - theta_2) + sin(phi_1)sin(phi_2)\ndp = left( cos(phi_1)cos(phi_2) + sin(phi_1)sin(phi_2) right) cos(theta_1 - theta_2) + (1 - cos(theta_1 - theta_2)) sin(phi_1)sin(phi_2)\ndp = cos(theta_1 - theta_2) cos(phi_1 - phi_2) + (1 - cos(theta_1 - theta_2)) sin(phi_1)sin(phi_2)\nSet A = cos(theta_1 - theta_2) B = cos(phi_1 - phi_2), then\ndp = A  B + (1 - A) sin(phi_1) sin(phi_2)\nThe angle between vectors in Radians is then\npsi = cos^-1(dp).\nDistance between the two points on the \"great circle\", Delta, is\nDelta = R psi\n\nArguments\n\ncoord1::Vector{Float64} – A 2-element vector: [lon, lat] in signed degrees.\ncoord2::Vector{Float64} – A 2-element vector: [lon, lat] in signed degrees.\nR::Float64              – The radius of the sphere (Default is Earth's radius in Kilometers.)\n\nReturn\n\nThe distance (in the units of the radius, R) via a \"great circle\" path.\n\n\n\n\n\n","category":"function"},{"location":"","page":"LatLon.jl Documentation","title":"LatLon.jl Documentation","text":"geo_midpoint","category":"page"},{"location":"#LatLon.geo_midpoint","page":"LatLon.jl Documentation","title":"LatLon.geo_midpoint","text":"geo_midpoint(coord1::Vector{Float{64}, \n             coord2::Vector{Float64}  )\n\nComputes and returns the mid-point of two points on the sphere as  a vector (lon/lat) in signed decimal degrees.\n\nArguments\n\ncoord1::Vector{Float64} - A 2-element vector: [lon, lat] in signed degrees.\ncoord2::Vector{Float64} - A 2-element vector: [lon, lat] in signed degrees.\n\nReturn\n\nThe 2-Vector representing the longitude and latitude.\n\n\n\n\n\n","category":"function"},{"location":"","page":"LatLon.jl Documentation","title":"LatLon.jl Documentation","text":"latlon_set_dist","category":"page"},{"location":"#LatLon.latlon_set_dist","page":"LatLon.jl Documentation","title":"LatLon.latlon_set_dist","text":"latlon_set_dist(coord1s::Vector{Float64}, \n                coord2s::Vector{Float64},\n                R=MC.radius::Float64     )\n\nComputes the set distance between two sets as represented by their lon/lat coordinates, in signed decimal degrees. It does this the hard way by computing the distance of all pairs of points between the two sets. The radius, R, defaults to the Earth's radius in Kilometers.\n\nArguments\n\ncoord1s::Matrix{Float64} – The coordinates of the first set.\ncoord2s::Matrix{Float64} – The coordinates of the second set.\nR::Float64               – Radius of sphere (Default is Earth's radius in Kilometers.)\n\nReturn\n\nThe minimum distance (in the same units as R) between the sets along with  the index of the points for each set representing the closest points from  each set.\n\nA Tuple: (distinkm, set1index, set2index)\n\n\n\n\n\n","category":"function"},{"location":"#Geo-Extraction-Functions","page":"LatLon.jl Documentation","title":"Geo-Extraction Functions","text":"","category":"section"},{"location":"","page":"LatLon.jl Documentation","title":"LatLon.jl Documentation","text":"center_latlon_from_NASA_xml_file","category":"page"},{"location":"#LatLon.center_latlon_from_NASA_xml_file","page":"LatLon.jl Documentation","title":"LatLon.center_latlon_from_NASA_xml_file","text":"center_latlon_from_NASA_xml_file(file::String)\n\nRetrieves the latitude and longitude from a \"center\" NASA KML file  as a matrix of Float64: 2xN. Here N is the number of points. Each column of this matrix (a 2-vector) is a lon/lat pair, in  signed decimal degrees.\n\nArguments\n\nfile::String – String representing a lon/lat file in NASA KML format.\n\nReturn\n\n::Matrix{Float64} – 2xN matrix of Lon/Lat pairs.\n\n\n\n\n\n","category":"function"},{"location":"","page":"LatLon.jl Documentation","title":"LatLon.jl Documentation","text":"upath_latlon_from_NASA_xml_file","category":"page"},{"location":"#LatLon.upath_latlon_from_NASA_xml_file","page":"LatLon.jl Documentation","title":"LatLon.upath_latlon_from_NASA_xml_file","text":"upath_latlon_from_NASA_xml_file(file::String)\n\nRetrieves the latitude and longitude of the annular ring representing the boundary of the \"total\" part of an eclipse from a \"upath\" NASA KML file as a matrix of Float64: 2xN. Here, N is the number of points. Each column of this matrix (a 2-vector) is a lon/lat pair, in  signed decimal degrees. The resulting set of points graphically sweeps out an annulus.\n\nArguments\n\nfile::String – String representing a  Lon/Lat file in NASA \"upath\" KML format.\n\nReturn\n\n::Matrix{Float64} – 2xN matrix of Lon/Lat pairs.\n\n\n\n\n\n","category":"function"},{"location":"#Index","page":"LatLon.jl Documentation","title":"Index","text":"","category":"section"},{"location":"","page":"LatLon.jl Documentation","title":"LatLon.jl Documentation","text":"","category":"page"}]
}
