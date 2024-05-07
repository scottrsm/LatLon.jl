using Test
using LatLon


@testset  "LatLon (Fidelity) " begin
    @test length(detect_ambiguities(LatLon)) == 0
end

@testset "Distance Between Sets" begin
	proj_path = dirname(pathof(LatLon))
	c23 = center_latlon_from_NASA_xml_file(joinpath(proj_path, "../data/2023/center.kml"))
	c24 = center_latlon_from_NASA_xml_file(joinpath(proj_path, "../data/2024/center.kml"))

  	dist = Int64(round(geo_dist(c23[:, 2350], c24[:, 2260]) * 5280 * 0.62137, digits=0))
	@test dist == 1027

	lonlat1 = [-98.1858158111572, 29.755859375]
	lonlat2 = [-100.949504375458, 31.771484375]
	@test  346.3308134412445 ≈ geo_dist(lonlat1, lonlat2) atol=1.0e-10
end

