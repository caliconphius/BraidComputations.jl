using BraidComputations

mkpath(joinpath(@__DIR__, "..", "build"))
B = BraidGroup(4)
x = B([1, 2, 1, 1, 2, 1, 3])

path1 = save_braid_svg(joinpath(@__DIR__, "..", "build", "braid.svg"), x;
                       orientation=:left_to_right)
path2 = save_braid_svg(joinpath(@__DIR__, "..", "build", "braid-vertical.svg"), x;
                       orientation=:top_to_bottom)
path3 = save_crs_svg(joinpath(@__DIR__, "..", "build", "crs.svg"), x)
path4 = save_outer_crs_svg(joinpath(@__DIR__, "..", "build", "outer-crs.svg"), x)

for path in (path1, path2, path3, path4)
    println(path)
end
