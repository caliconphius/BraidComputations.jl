using Test
using Random
using BraidComputations
using Farey
using Farey.Monoids

function all_permutations(n)
    n == 1 && return [[1]]
    result = Vector{Vector{Int}}()
    for prefix in all_permutations(n - 1), position in 1:n
        values = copy(prefix)
        insert!(values, position, n)
        push!(result, values)
    end
    result
end

include("permutations.jl")
include("transducer.jl")
include("left_normal_forms.jl")
include("garside_operations.jl")
include("bkl.jl")
include("bkl_normal_forms.jl")
include("sliding_circuits.jl")
include("canonical_reduction.jl")
include("visualize.jl")
include("braids.jl")
include("words.jl")
include("huxford_schillewaert.jl")
