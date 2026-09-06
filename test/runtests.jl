using Test
using LatLon

const R_EARTH = 6371.0

# Independent reference: the Haversine formula (accurate at short range).
function haversine(c1, c2, R=R_EARTH)
	d2r = π / 180
	ϕ1 = c1[2] * d2r; ϕ2 = c2[2] * d2r
	dϕ = ϕ2 - ϕ1
	dθ = (c2[1] - c1[1]) * d2r
	a = sin(dϕ / 2)^2 + cos(ϕ1) * cos(ϕ2) * sin(dθ / 2)^2
	return 2R * asin(min(1.0, sqrt(a)))
end

@testset  "LatLon (Fidelity) " begin
    @test length(detect_ambiguities(LatLon)) == 0
end

@testset "KML parsing" begin
	proj_path = dirname(pathof(LatLon))
	c23 = center_latlon_from_NASA_xml_file(joinpath(proj_path, "../data/2023/center.kml"))
	u23 = upath_latlon_from_NASA_xml_file(joinpath(proj_path, "../data/2023/upath_hi.kml"))
	@test size(c23, 1) == 2 && size(c23, 2) > 100
	@test size(u23, 1) == 2 && size(u23, 2) > 100
	@test all(-180 .<= c23[1, :] .<= 180) && all(-90 .<= c23[2, :] .<= 90)

	# Whitespace and altitude fields in a coordinates element are handled.
	mktempdir() do dir
		f = joinpath(dir, "t.kml")
		write(f, """<?xml version="1.0" encoding="UTF-8"?>
		<kml><Document><Folder><Placemark><LineString><coordinates>
		  -1.5,2.25,0 
		  -3,4
		</coordinates></LineString></Placemark></Folder></Document></kml>""")
		@test center_latlon_from_NASA_xml_file(f) == [-1.5 -3.0; 2.25 4.0]
	end
end

@testset "Distance Between Sets" begin
	proj_path = dirname(pathof(LatLon))
	c23 = center_latlon_from_NASA_xml_file(joinpath(proj_path, "../data/2023/center.kml"))
	c24 = center_latlon_from_NASA_xml_file(joinpath(proj_path, "../data/2024/center.kml"))

  	dist = Int64(round(geo_dist(c23[:, 2350], c24[:, 2260]) * 5280 * 0.62137, digits=0))
	@test dist == 1027

	lonlat1 = [-98.1858158111572, 29.755859375]
	lonlat2 = [-100.949504375458, 31.771484375]
	@test  346.3308134412445 ≈ geo_dist(lonlat1, lonlat2) atol=1.0e-8

	# The closest pair between two small sets.
	A = [0.0 10.0 20.0; 0.0 0.0 0.0]
	B = [10.5 50.0; 0.0 0.0]
	d, i, j = latlon_set_dist(A, B)
	@test (i, j) == (2, 1)
	@test d ≈ geo_dist([10.0, 0.0], [10.5, 0.0])
	@test latlon_set_dist(A, B, 1.0)[1] ≈ 0.5 * π / 180

	# Any AbstractMatrix{<:Real} works.
	@test latlon_set_dist(view(A, :, 1:2), B)[2] == 2
	@test latlon_set_dist([0 10; 0 0], [10 50; 0 0])[1] ≈ 0.0 atol=1e-12

	# Contract violations.
	@test_throws DimensionMismatch latlon_set_dist(A', B)
	@test_throws DimensionMismatch latlon_set_dist(A, zeros(3, 2))
	@test_throws ArgumentError latlon_set_dist(zeros(2, 0), B)
end

@testset "geo_dist" begin
	# Identical points, symmetry, the equator, and the antipode.
	@test geo_dist([10.0, 20.0], [10.0, 20.0]) == 0.0
	@test geo_dist([0.0, 0.0], [90.0, 0.0]) ≈ R_EARTH * π / 2
	@test geo_dist([0.0, 0.0], [180.0, 0.0]) ≈ R_EARTH * π
	@test geo_dist([0.0, 90.0], [0.0, -90.0]) ≈ R_EARTH * π
	@test geo_dist([-98.0, 30.0], [10.0, -40.0]) == geo_dist([10.0, -40.0], [-98.0, 30.0])

	# Longitude wrap-around at ±180°.
	@test geo_dist([179.9, 0.0], [-179.9, 0.0]) ≈ haversine([179.9, 0.0], [-179.9, 0.0])

	# Agreement with Haversine from continental scale down to centimetres.
	c1 = [-98.0, 30.0]
	for dlat in (10.0, 1.0, 1e-2, 1e-4, 1e-5, 1e-6, 1e-7)
		c2 = [-98.0, 30.0 + dlat]
		@test geo_dist(c1, c2) ≈ haversine(c1, c2) rtol=1e-9
	end
	@test geo_dist([-98.0, 30.0], [-98.0 + 1e-7, 30.0]) > 0.0

	# Other element types and radii.
	@test geo_dist([0, 0], [90, 0]) ≈ R_EARTH * π / 2
	@test geo_dist(Float32[0, 0], Float32[90, 0]) ≈ R_EARTH * π / 2
	@test geo_dist([0.0, 0.0], [90.0, 0.0], 1) ≈ π / 2
	@test_throws DimensionMismatch geo_dist([0.0, 0.0, 0.0], [90.0, 0.0])
end

@testset "geo_midpoint" begin
	@test geo_midpoint([0.0, 0.0], [90.0, 0.0]) ≈ [45.0, 0.0]
	@test geo_midpoint([10.0, 20.0], [10.0, 20.0]) ≈ [10.0, 20.0]
	@test geo_midpoint([0.0, 89.0], [180.0, 89.0]) ≈ [90.0, 90.0]
	@test geo_midpoint([0, 0], [0, 90]) ≈ [0.0, 45.0]
	# The mid-point is equidistant from both ends.
	a = [-98.0, 30.0]; b = [10.0, -40.0]
	m = geo_midpoint(a, b)
	@test geo_dist(a, m) ≈ geo_dist(m, b)
	@test geo_dist(a, m) + geo_dist(m, b) ≈ geo_dist(a, b)
	@test_throws DomainError geo_midpoint([0.0, 0.0], [180.0, 0.0])
	@test_throws DimensionMismatch geo_midpoint([0.0], [180.0, 0.0])
end
