module LatLon

import LightXML as LX

export geo_midpoint,  geo_dist,  latlon_set_dist
export center_latlon_from_NASA_xml_file, upath_latlon_from_NASA_xml_file


# Internal structure that contains constants used by the functions below.
struct ModuleConsts
	deg_2_rad::Float64 # Conversion factor from Degrees to Radians.
	rad_2_deg::Float64 # Conversion factor from Radians to Degrees.
	radius::Float64    # The radius of the sphere of interest.
end

# The single instance of the ModuleConsts structure.
# Note: The radius is set to the radius of the Earth in Kilometers.
const MC = ModuleConsts(π / 180.0, 180.0 / π, 6371.0)



"""
    center_latlon_from_NASA_xml_file(file::String)

Retrieves the latitude and longitude from a "center" NASA kml file 
as a matrix of `Float64: 2xN`. Here `N` is the number of points.
Each column of this matrix (a 2-vector) is a lon/lat pair, in 
signed decimal degrees.
"""
function center_latlon_from_NASA_xml_file(file::String)
    kml_doc = LX.parse_file(file)
    droot = LX.root(kml_doc)
    data = LX.content(droot["Document"][1]["Folder"][1]["Placemark"][1]["LineString"][1]["coordinates"][1])
    mat = stack(map(x -> split(x, ","), split(data, " ")))
	return map(x -> parse(Float64, x), mat)
end



"""
    upath_latlon_from_NASA_xml_file(file::String)

Retrieves the latitude and longitude of the annular ring representing
the boundary of the "total" part of an eclipse from a "upath" NASA kml file
as a matrix of `Float64: 2xN`. Here, `N` is the number of points.
Each column of this matrix (a 2-vector) is a lon/lat pair, in 
signed decimal degrees.
"""
function upath_latlon_from_NASA_xml_file(file::String)
    kml_doc = LX.parse_file(file)
    rt = LX.root(kml_doc)
    data = LX.content(rt["Document"][1]["Folder"][1]["Placemark"][1]["Polygon"][1]["outerBoundaryIs"][1]["LinearRing"][1]["coordinates"][1])
    mat = stack(map(x -> split(x, ","), split(data, " ")))
    return map(x -> parse(Float64, x), mat)
end



"""
	geo_midpoint(coord1::Vector{Float{64}, 
                 coord2::Vector{Float64}  )

Computes and returns the mid-point of two points on the sphere as 
a vector (lon/lat) in signed decimal degrees.

# Arguments
- coord1::Vector{Float64} - A 2-element vector: [lon, lat] in signed degrees.
- coord2::Vector{Float64} - A 2-element vector: [lon, lat] in signed degrees.

# Return
The 2-Vector representing the longitude and latitude.
"""
function geo_midpoint(coord1::Vector{Float64}, 
					  coord2::Vector{Float64} )

    # Convert angles to radians.
    θ1   = coord1[1] * MC.deg_2_rad
    θ2   = coord2[1] * MC.deg_2_rad
    ϕ1   = coord1[2] * MC.deg_2_rad
    ϕ2   = coord2[2] * MC.deg_2_rad

    # Get the cosine of the latitudes.
    cph1 = cos(ϕ1)
    cph2 = cos(ϕ2)

    cc1  = cos(θ1) * cph1
    cs1  = sin(θ1) * cph1
    cc2  = cos(θ2) * cph2
    cs2  = sin(θ2) * cph2

    # The Cartesian coordinates of the two points projected onto the unit sphere.
    p1   = [cc1, cs1, sin(ϕ1)]
    p2   = [cc2, cs2, sin(ϕ2)]

    # Create mid point (no need to divide by 2 as we will normalize)
    mp = p1 .+ p2

    # Project to the point on the surface of the sphere.
    mp /= sqrt(sum(mp .* mp))

    # Get the latitude, ϕ, and longitude, θ.
    ϕ  = asin(mp[3])
    θ  = atan(mp[2], mp[1])

    # Return the lon/lat vector.
    return MC.rad_2_deg .* [θ, ϕ]
end



"""
	geo_dist(coord1::Vector{Float64}  , 
             coord2::Vector{Float64}  ,
             radius=MC.radius::Float64 )

Computes the distance between two coordinates on a sphere represented as 
two lon/lat vectors in signed degrees. That is, north latitude is 
positive, south latitude is negative, while east longitude is positive 
and west longitude is negative. The radius, R, defaults to the Earth's radius
in Kilometers.

# Details
  Dot product, `dp`, of the geo-positions on the unit sphere gives the cosine of 
  the angle between the two points (on the "great circle").
  This is the information contained in `coord1` and `coord2`.
  Once the cosine is computed we can easily retrieve the angle (in Radians) and 
  then it is easy to find the distance between the two points 
  -- the length of the arc connecting them on a "great circle".
- The Cartesian points of the lon/lat on the unit sphere are:\n
  v1 = ( cos(θ1)cos(ϕ1), sin(θ1)cos(ϕ1), sin(ϕ1) )\n
  v2 = ( cos(θ2)cos(ϕ2), sin(θ2)cos(ϕ2), sin(ϕ2) )
- The dot product becomes:\n
  dp = v1 ⋅ v2 = cos(ϕ1)cos(ϕ2) ( cos(θ1)cos(θ2) + sin(θ1)sin(θ2) ) + sin(ϕ1)sin(ϕ2)
- Simplifying\n
  dp = cos(ϕ1)cos(ϕ2) cos(θ1 - θ2) + sin(ϕ1)sin(ϕ2)\n
  dp = ( cos(ϕ1)cos(ϕ2) + sin(ϕ1)sin(ϕ2) ) cos(θ1 - θ2) + (1 - cos(θ1 - θ2)) sin(ϕ1)sin(ϕ2)\n
  dp = cos(θ1 - θ2) cos(ϕ1 - ϕ2) + (1 - cos(θ1 - θ2)) sin(ϕ1)sin(ϕ2)

  Set A = cos(θ1 - θ2), B = cos(ϕ1 - ϕ2), then\n
  dp = A B + (1 - A) sin(ϕ1) sin(ϕ2)
- The angle between vectors in Radians is then
  ψ = acos(dp)
- Distance between the two points on the "great circle", `Δ`, is\n
  Δ = R ψ

# Arguments
- coord1::Vector{Float64} -- A 2-element vector: [lon, lat] in signed degrees.
- coord2::Vector{Float64} -- A 2-element vector: [lon, lat] in signed degrees.
- R::Float64              -- The radius of the sphere (Default is Earth's radius in Kilometers.)

# Return
The distance in (in the units of the radius, `R`) via a "great circle" path.
"""
function geo_dist(coord1::Vector{Float64}, 
                  coord2::Vector{Float64},
                  R=MC.radius::Float64    )
    # Longitude differences.
    θ1   = coord1[1] 
    θ2   = coord2[1] 
    dθ   = (θ1 - θ2) * MC.deg_2_rad

    # Latitude differences.
    ϕ1   = coord1[2] * MC.deg_2_rad
    ϕ2   = coord2[2] * MC.deg_2_rad
    dϕ   = ϕ1 - ϕ2 

    # Cosines of the differences between lat/lon angles.
    A = cos(dθ)
    B = cos(dϕ)

    # The dot product of the Cartesian coordinates of the two points 
    # that the lat/lon coordinates represent (on the unit sphere) is
    # the cosine of the angle between the two points. 
    # See the documentation section "Details" above.
    dp = A * B + (1 - A) * sin(ϕ1) * sin(ϕ2)
    
    # Retrieve the angle from the cosine of the angle between the two points in Radians.
    # Then easily compute the distance between the two points on the sphere  based on its radius.
    return R * acos(dp)
end


"""
	latlon_set_dist(coord1s::Vector{Float64}, 
                    coord2s::Vector{Float64},
                    R=MC.radius::Float64     )

Computes the set distance between two sets as represented by their lat/lon coordinates,
in singed decimal degrees.
It does this the hard way by computing the distance of all pairs of points
between the two sets. The radius, R, defaults to the Earth's radius in Kilometers.

# Arguments
- coord1s::Matrix{Float64} -- The coordinates of the first set.
- coord2s::Matrix{Float64} -- The coordinates of the second set.
- R::Float64               -- Radius of sphere (Default is Earth's radius in Kilometers.)

# Return
The minimum distance (in the same units as `R`) between the sets along with 
the index of the points for each set representing the closest points from 
each set.

A Tuple: (dist_in_km, set1_index, set2_index)
"""
function latlon_set_dist(coord1s::Matrix{Float64}, 
                         coord2s::Matrix{Float64},
                         R=MC.radius::Float64     )
    dmin = Inf
    min1 = 0
    min2 = 0
    _, N1 = size(coord1s)
    _, N2 = size(coord2s)
    @inbounds for i1 in 1:N1
        p1 = coord1s[:, i1]
        for i2 in 1:N2
            d = geo_dist(p1, coord2s[:, i2], R)
            if d < dmin
                dmin = d
                min1 = i1
                min2 = i2
            end
        end
    end
    return (dmin, min1, min2)
end


end # LatLon Module

