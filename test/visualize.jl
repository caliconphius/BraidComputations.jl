@testset "dependency-free SVG visualization" begin
    x = BraidGroup(4)([1, 2, 1, -3])
    partner_dir = mktempdir()
    for orientation in (:left_to_right, :top_to_bottom)
        image = braid_svg(x; orientation)
        @test startswith(image, "<?xml")
        @test occursin("<svg", image)
        @test occursin("σ1", image)
        @test occursin("σ3⁻", image)
        @test endswith(image, "</svg>\n")
        if orientation == :top_to_bottom
            @test occursin("translate", image)
        end
    end

    @test braid_svg(Int[], 3) |> x -> occursin("Braid on 3 strands", x)
    red = BraidGroup(4)(1)
    svg = crs_svg(red)
    outer = outer_crs_svg(red)
    @test occursin("A[1,1]", svg)
    @test occursin("<ellipse", outer)
    @test length(findall("<ellipse", outer)) >= 2
    @test occursin("empty CRS", crs_svg(BraidGroup(4)(Int[])))

    path = save_braid_svg(joinpath(partner_dir, "example.svg"), x)
    @test read(path, String) == braid_svg(x)
    path2 = save_crs_svg(joinpath(partner_dir, "crs.svg"), red)
    @test read(path2, String) == crs_svg(red)
    path3 = save_outer_crs_svg(joinpath(partner_dir, "outer.svg"), red)
    @test read(path3, String) == outer_crs_svg(red)
    @test_throws ArgumentError braid_svg(BraidGroup(3)(1); orientation=:diagonal)
end
