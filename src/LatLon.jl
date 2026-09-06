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



#= Parse the text of a KML `<coordinates>` element into a 2xN matrix of (lon, lat) columns.
   Tuples are separated by any whitespace; each tuple is "lon,lat" or "lon,lat,alt"
   (the altitude, if present, is dropped).
=#
function _parse_kml_coordinates(data::AbstractString)
    tuples = split(strip(data))
    mat = Matrix{Float64}(undef, 2, length(tuples))
    for (j, t) in enumerate(tuples)
        fields = split(t, ",")
        length(fields) >= 2 || throw(ArgumentError("KML coordinate tuple must have at least a longitude and a latitude: $(repr(t))"))
        mat[1, j] = parse(Float64, fields[1])
        mat[2, j] = parse(Float64, fields[2])
    end
    return mat
end


"""
    center_latlon_from_NASA_xml_file(file::String)

Retrieves the latitude and longitude from a "center" NASA KML file 
as a matrix of `Float64: 2xN`. Here `N` is the number of points.
Each column of this matrix (a 2-vector) is a lon/lat pair, in 
signed decimal degrees.

# Arguments
- file::String -- String representing a lon/lat file in NASA KML format.

# Return
::Matrix{Float64} -- 2xN matrix of Lon/Lat pairs.
"""
function center_latlon_from_NASA_xml_file(file::String)
    kml_doc = LX.parse_file(file)
    try
        droot = LX.root(kml_doc)
        data = LX.content(droot["Document"][1]["Folder"][1]["Placemark"][1]["LineString"][1]["coordinates"][1])
        return _parse_kml_coordinates(data)
    finally
        LX.free(kml_doc)
    end
end



"""
    upath_latlon_from_NASA_xml_file(file::String)

Retrieves the latitude and longitude of the annular ring representing
the boundary of the "total" part of an eclipse from a "upath" NASA KML file
as a matrix of `Float64: 2xN`. Here, `N` is the number of points.
Each column of this matrix (a 2-vector) is a lon/lat pair, in 
signed decimal degrees.
The resulting set of points graphically sweeps out an annulus.

# Arguments
- file::String -- String representing a  Lon/Lat file in NASA "upath" KML format.

# Return
::Matrix{Float64} -- 2xN matrix of Lon/Lat pairs.
"""
function upath_latlon_from_NASA_xml_file(file::String)
    kml_doc = LX.parse_file(file)
    try
        rt = LX.root(kml_doc)
        data = LX.content(rt["Document"][1]["Folder"][1]["Placemark"][1]["Polygon"][1]["outerBoundaryIs"][1]["LinearRing"][1]["coordinates"][1])
        return _parse_kml_coordinates(data)
    finally
        LX.free(kml_doc)
    end
end



"""
	geo_midpoint(coord1::AbstractVector{<:Real}, 
                 coord2::AbstractVector{<:Real})

Computes and returns the mid-point of two points on the sphere as 
a vector (lon/lat) in signed decimal degrees.

# Arguments
- coord1::AbstractVector{<:Real} - A 2-element vector: [lon, lat] in signed degrees.
- coord2::AbstractVector{<:Real} - A 2-element vector: [lon, lat] in signed degrees.

# Return
The 2-Vector representing the longitude and latitude.

Throws a `DomainError` if the points are antipodal (the mid-point is not unique).
"""
function geo_midpoint(coord1::AbstractVector{<:Real}, 
					  coord2::AbstractVector{<:Real})
    _check_coord(coord1)
    _check_coord(coord2)

    # Convert angles to radians.
    θ1   = float(coord1[1]) * MC.deg_2_rad
    θ2   = float(coord2[1]) * MC.deg_2_rad
    ϕ1   = float(coord1[2]) * MC.deg_2_rad
    ϕ2   = float(coord2[2]) * MC.deg_2_rad

    # Get the cosine of the latitudes.
    cph1 = cos(ϕ1)
    cph2 = cos(ϕ2)

    # The Cartesian coordinates of the two points projected onto the unit sphere,
    # summed to get the (un-normalized) mid point.
    x = cos(θ1) * cph1 + cos(θ2) * cph2
    y = sin(θ1) * cph1 + sin(θ2) * cph2
    z = sin(ϕ1) + sin(ϕ2)

    # Project to the point on the surface of the sphere.
    nrm = sqrt(x * x + y * y + z * z)
    if nrm < 1.0e-15
        throw(DomainError((coord1, coord2), "geo_midpoint is undefined for antipodal points"))
    end

    # Get the latitude, ϕ, and longitude, θ.
    ϕ  = asin(clamp(z / nrm, -1.0, 1.0))
    θ  = atan(y, x)

    # Return the lon/lat vector.
    return MC.rad_2_deg .* [θ, ϕ]
end

# A coordinate is a 2-element vector: [lon, lat].
function _check_coord(c::AbstractVector)
    length(c) == 2 || throw(DimensionMismatch("a lon/lat coordinate must be a 2-element vector; got length $(length(c))"))
    return nothing
end



"""
	geo_dist(coord1::AbstractVector{<:Real}, 
             coord2::AbstractVector{<:Real},
             R::Real=MC.radius              )

Computes the distance between two points on a sphere represented as 
two lon/lat vectors in signed degrees. That is, north latitude is 
positive, south latitude is negative, while east longitude is positive 
and west longitude is negative. The radius, R, defaults to the Earth's radius
in Kilometers.

**NOTE:** The distance function described below is **NOT** the usual *Haversine* formula.
It works for the entire sphere and, unlike a formula based on the arc-cosine of the
dot product, it is accurate for both very small and near-antipodal separations
(the angle is recovered with a two-argument arc-tangent of the sine and cosine of the
angle, which is well conditioned everywhere).

# Details
  The dot product, `dp`, of the geo-positions on the unit sphere gives the cosine of 
  the angle between the two points on the "great circle" connecting them, and the
  norm of their cross product gives the sine of that angle.
  Inputs `coord1` and `coord2` contain the positions of the two points on the 
  unit sphere with respect to latitude and longitude.
  Once the sine and cosine are computed we can retrieve the angle (in Radians) and 
  then it is easy to find the distance between the two points 
  -- the length of the "great circle" arc connecting them.
- The Cartesian points of the lon/lat on the unit sphere are:\n
  ``{\\bf v}_1 = ( \\cos(\\theta_1)\\cos(\\phi_1), \\sin(\\theta_1)\\cos(\\phi_1), \\sin(\\phi_1) )``\n
  ``{\\bf v}_2 = ( \\cos(\\theta_2)\\cos(\\phi_2), \\sin(\\theta_2)\\cos(\\phi_2), \\sin(\\phi_2) )``
- The dot product becomes:\n
  ``dp = {\\bf v}_1 {\\cdot} {\\bf v}_2 = \\cos(\\phi_1)\\cos(\\phi_2) \\left( \\cos(\\theta_1)\\cos(\\theta_2) + \\sin(\\theta_1)\\sin(\\theta_2) \\right) + \\sin(\\phi_1)\\sin(\\phi_2)``
- Simplifying\n
  ``dp = \\cos(\\phi_1)\\cos(\\phi_2) \\cos(\\theta_1 - \\theta_2) + \\sin(\\phi_1)\\sin(\\phi_2)``\n
  ``dp = \\left( \\cos(\\phi_1)\\cos(\\phi_2) + \\sin(\\phi_1)\\sin(\\phi_2) \\right) \\cos(\\theta_1 - \\theta_2) + (1 - \\cos(\\theta_1 - \\theta_2)) \\sin(\\phi_1)\\sin(\\phi_2)``\n
  ``dp = \\cos(\\theta_1 - \\theta_2) \\cos(\\phi_1 - \\phi_2) + (1 - \\cos(\\theta_1 - \\theta_2)) \\sin(\\phi_1)\\sin(\\phi_2)``
- The norm of the cross product, ``|{\\bf v}_1 \\times {\\bf v}_2|``, is the sine of the
  same angle. With ``\\Delta\\theta = \\theta_1 - \\theta_2`` it is\n
  ``\\sqrt{ \\left(\\cos(\\phi_2)\\sin(\\Delta\\theta)\\right)^2 + \\left(\\cos(\\phi_1)\\sin(\\phi_2) - \\sin(\\phi_1)\\cos(\\phi_2)\\cos(\\Delta\\theta)\\right)^2 }``
- Procedure to Compute Distance:
    - Set Intermediate Variables:\n
         ``A = \\cos(\\theta_1 - \\theta_2), \\, B = \\cos(\\phi_1 - \\phi_2)``\n
         ``dp = A \\, B + (1 - A) \\sin(\\phi_1) \\sin(\\phi_2)`` (the cosine of the angle)\n
         ``cp = |{\\bf v}_1 \\times {\\bf v}_2|`` (the sine of the angle)
    - The angle between vectors in Radians is\n
         ``\\psi = {\\rm atan}(cp, dp)``.
    - Distance between the two points on the "great circle" is\n
         ``\\Delta = R \\, \\psi``

# Arguments
- coord1::AbstractVector{<:Real} -- A 2-element vector: [lon, lat] in signed degrees.
- coord2::AbstractVector{<:Real} -- A 2-element vector: [lon, lat] in signed degrees.
- R::Real                        -- The radius of the sphere (Default is Earth's radius in Kilometers.)

# Return
The distance (in the units of the radius, `R`) via a "great circle" path.
"""
function geo_dist(coord1::AbstractVector{<:Real}, 
                  coord2::AbstractVector{<:Real},
                  R::Real=MC.radius              )
    _check_coord(coord1)
    _check_coord(coord2)

    # Longitude differences.
    θ1   = float(coord1[1])
    θ2   = float(coord2[1])
    dθ   = (θ1 - θ2) * MC.deg_2_rad

    # Latitudes and their difference.
    ϕ1   = float(coord1[2]) * MC.deg_2_rad
    ϕ2   = float(coord2[2]) * MC.deg_2_rad
    dϕ   = ϕ1 - ϕ2 

    # Cosines of the differences between lat/lon angles.
    A = cos(dθ)
    B = cos(dϕ)

    sϕ1 = sin(ϕ1)
    sϕ2 = sin(ϕ2)
    cϕ1 = cos(ϕ1)
    cϕ2 = cos(ϕ2)

    #= The dot product of the Cartesian coordinates of the two points 
       that the lat/lon coordinates represent (on the unit sphere) is
       the cosine of the angle between the two points. 
       See the documentation section "Details" above.
	=#
    dp = A * B + (1 - A) * sϕ1 * sϕ2

    #= The norm of the cross product of the same two points is the sine of the angle.
       The second component, cos(ϕ1)sin(ϕ2) - sin(ϕ1)cos(ϕ2)cos(dθ), is rewritten as
       sin(ϕ2 - ϕ1) + sin(ϕ1)cos(ϕ2)(1 - cos(dθ)) with 1 - cos(dθ) = 2 sin²(dθ/2),
       which avoids cancellation for nearby points.
	=#
    cx = cϕ2 * sin(dθ)
    cy = -sin(dϕ) + sϕ1 * cϕ2 * (2 * sin(dθ / 2)^2)
    cp = sqrt(cx * cx + cy * cy)
    
    # Retrieve the angle (in Radians) from its sine and cosine -- well conditioned for all separations.
    # Then compute the distance between the two points on the sphere based on its radius.
    return R * atan(cp, dp)
end


"""
	latlon_set_dist(coord1s::AbstractMatrix{<:Real}, 
                    coord2s::AbstractMatrix{<:Real},
                    R::Real=MC.radius               )

Computes the set distance between two sets as represented by their lon/lat coordinates,
in signed decimal degrees.
It does this the hard way by computing the distance of all pairs of points
between the two sets. The radius, R, defaults to the Earth's radius in Kilometers.

# Arguments
- coord1s::AbstractMatrix{<:Real} -- The coordinates of the first set: a `2xN1` matrix, each column is [lon, lat].
- coord2s::AbstractMatrix{<:Real} -- The coordinates of the second set: a `2xN2` matrix, each column is [lon, lat].
- R::Real                         -- Radius of sphere (Default is Earth's radius in Kilometers.)

# Input Contract
- `size(coord1s, 1) == size(coord2s, 1) == 2`
- Both sets are non-empty.

# Return
The minimum distance (in the same units as `R`) between the sets along with 
the index of the points for each set representing the closest points from 
each set.

A Tuple: (dist_in_km, set1_index, set2_index)
"""
function latlon_set_dist(coord1s::AbstractMatrix{<:Real}, 
                         coord2s::AbstractMatrix{<:Real},
                         R::Real=MC.radius               )
    size(coord1s, 1) == 2 || throw(DimensionMismatch("latlon_set_dist: `coord1s` must be a 2xN matrix (columns are [lon, lat]); got size $(size(coord1s))"))
    size(coord2s, 1) == 2 || throw(DimensionMismatch("latlon_set_dist: `coord2s` must be a 2xN matrix (columns are [lon, lat]); got size $(size(coord2s))"))
    (size(coord1s, 2) > 0 && size(coord2s, 2) > 0) || throw(ArgumentError("latlon_set_dist: both coordinate sets must be non-empty"))

    dmin = Inf
    min1 = 0
    min2 = 0
    _, N1 = size(coord1s)
    _, N2 = size(coord2s)
    for i1 in 1:N1
        p1 = @view coord1s[:, i1]
        for i2 in 1:N2
            d = geo_dist(p1, @view(coord2s[:, i2]), R)
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

