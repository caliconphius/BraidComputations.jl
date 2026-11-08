@testset "Huxford–Schillewaert example computations" begin
    include(abspath(@__DIR__, "..", "examples", "huxford_schillewaert.jl"))
    @test length(hs_words_up_to(6)) == 1457
    @test all(w -> !hs_equation(w, 1), filter(!isempty, words))
end
