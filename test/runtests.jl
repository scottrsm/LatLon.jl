using Test
using LatLon


@testset  "LatLon (Fidelity) " begin
    @test length(detect_ambiguities(LatLon)) == 0
end

