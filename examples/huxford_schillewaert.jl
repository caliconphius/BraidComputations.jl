using BraidComputations
using Farey.Monoids

# Explicit calculations in Huxford–Schillewaert, *Braid groups, elliptic
# curves, and resolving the quartic* (arXiv:2309.12999). These checks reproduce
# the algebraic computations; they are not a mechanized proof of classification.
# Related UI conventions appear in packages/Farey/working/kassel_morphs.jl
# (ϕ = G1i ∘ D1 and ψ = D1i ∘ G1), for later Farey/Christoffel comparison.

println("Section 2: presentations, quotients, and matrices")
B3 = BraidGroup(3)
B4 = BraidGroup(4)
α3 = B3([1, 2])
β3 = B3([1, 2, 1])
α4 = B4([1, 2, 3])
β4 = α4 * B4(1)

@assert α3^3 == β3^2
@assert is_central(α3^3)
@assert α4^4 == β4^3
@assert is_central(α4^4)

# R_*(σ1)=R_*(σ3)=σ1, R_*(σ2)=σ2.
R_images = Dict(1 => B3(1), 2 => B3(2), 3 => B3(1))
R(word) = prod([R_images[abs(i)]^sign(i) for i in word]; init=one(B3))
@assert R([1, 2, 1]) == R([2, 1, 2])
@assert R([1, 3]) == R([3, 1])

T = [1 1; 0 1]
L = [1 0; -1 1]
psl_equal(A, B) = A == B || A == -B
@assert T*L*T == L*T*L
minus_identity = Int[-1 0; 0 -1]
@assert psl_equal((T*L)^3, minus_identity)
@assert psl_equal((T*L*T)^2, minus_identity)

g = B3([2, 1, -2])
rg_α = g * α3 * inv(g)
rg_β = β3
@assert rg_α^3 == rg_β^2
@assert canonical_reduction_system(B3(1)).standard == ParabolicSystem([StandardParabolic(1, 1, 3)])

println("Section 3: periodic braids and outer canonical reduction systems")
@assert isempty(canonical_reduction_system(α3))
@assert nielsen_thurston_type(α3) == :periodic
@assert canonical_reduction_system(B3(1)).standard == ParabolicSystem([StandardParabolic(1, 1, 3)])

# Distinct explicit outer-reduction configurations in B4.
@assert outer_crs_partition(B4(1)*B4(3)) == [2, 2]
@assert outer_crs_partition(B4(1)) == [2, 1, 1]
@assert outer_crs_partition(B4([1, 2, 1, 1, 2, 1])) == [3, 1]
@show outer_crs_partition(B4(1)*B4(3))
@show outer_crs_partition(B4(1))
@show outer_crs_partition(B4([1, 2, 1, 1, 2, 1]))

println("Section 4: endomorphism forms")
@assert rg_α^3 == rg_β^2
@assert inv(g)*rg_α*g == α3

println("Section 5: irreducible homomorphisms and Gassner's splitting")
ψ1 = B4([1, 2])
ψ2 = B4([3, 2])
@assert ψ1*ψ2*ψ1 == ψ2*ψ1*ψ2
ψ_z3 = (ψ1*ψ2)^3
@show nielsen_thurston_type(ψ_z3)

x = B4([-1, 3])
y = inv(x) * inv(B4(2)) * x * B4(2)
σ1 = B4(1)
σ2 = B4(2)
@assert R(letters(x)) == one(B3)
@assert R(letters(y)) == one(B3)
@assert inv(σ1)*x*σ1 == x
@assert inv(σ2)*x*σ2 == x*y
@assert inv(σ1)*y*σ1 == y*inv(x)
@assert inv(σ2)*y*σ2 == y

# Conjugation in F2 by σ1σ2 induces f(x)=xy, f(y)=x^{-1}.
h = B4([1, 2])
@assert inv(h)*x*h == x*y
@assert inv(h)*y*h == inv(x)

println("Section 6: explicit F2 automorphism and bounded equation verification")
# Use signed integral letters ±1, ±2 for x±1, y±1. Products always reduce
# freely, independently of the braid/Garside normalizers.
function hs_reduce(word)
    result = Int[]
    for letter in word
        if !isempty(result) && result[end] == -letter
            pop!(result)
        else
            push!(result, letter)
        end
    end
    result
end
hs_prod(a, b) = hs_reduce(vcat(a, b))
hs_inv(w) = -reverse(w)
hs_pow(w, n::Integer) = n == 0 ? Int[] : n > 0 ? hs_prod(collect(Iterators.flatten(fill(w, n)))) : hs_pow(hs_inv(w), -n)

function hs_f(w, inverse=false)
    hs_reduce(mapreduce(letter -> begin
        if inverse
            sign(letter) > 0 ? (abs(letter) == 1 ? [-2] : [2, 1]) :
                               (abs(letter) == 1 ? [2] : [-1, -2])
        else
            abs(letter) == 1 ? (sign(letter) > 0 ? [1, 2] : [-2, -1]) : [-sign(letter)]
        end
    end, hs_prod, w; init=Int[]))
end

function hs_iterate_f(w, n::Integer)
    result = collect(w)
    step = n < 0
    for _ in 1:abs(n)
        result = hs_f(result, step)
    end
    result
end

c = [1, 2, -1, -2]
@assert hs_f(c) == c
for w in ([1], [2], [1, 2, -2, 1], c, [1, -1, 2, -2, 2])
    @assert hs_iterate_f(w, 6) == hs_reduce(vcat(c, w, hs_inv(c)))
end
@assert hs_iterate_f([1], -1) == [-2]
for w in ([1], [2], [-1], [-2], c, [1, 2, -1, 2, -2])
    @assert hs_f(hs_f(w, true), false) == hs_reduce(w)
    @assert hs_f(hs_f(w, false), true) == hs_reduce(w)
end
@assert hs_f([1], true) == [-2]

hs_equation(w, k) = hs_prod(hs_iterate_f(w, 6k+1), hs_iterate_f(w, -6k-1)) == hs_reduce(w)

function hs_words_up_to(maxlength)
    result = [Int[]]
    frontier = [Int[]]
    for _ in 1:maxlength
        next = Vector{Int}[]
        for w in frontier, letter in (1, -1, 2, -2)
            v = hs_reduce([w; letter])
            length(v) == length(w)+1 || continue
            push!(result, v)
            push!(next, v)
        end
        frontier = next
    end
    result
end

words = hs_words_up_to(6)
solutions = filter(w -> hs_equation(w, 0), words)
orbit_x = [hs_iterate_f([1], n) for n in -6:12]
solution_set = Set{Vector{Int}}(orbit_x)
@assert all(w -> isempty(w) || w in solution_set, solutions)
@assert !isempty(solutions)
@assert all(w -> isempty(w) || !hs_equation(w, 1), words)
@assert all(w -> isempty(w) || !hs_equation(w, -1), words)

println("  Reduced words checked (length ≤ 6): $(length(words))")
println("  Nonidentity equation solutions at k=0: $(sum(!isempty, solutions))")
println("  f-orbit of x for n=0:5: $([hs_iterate_f([1], n) for n in 0:5])")

println("All explicit Huxford–Schillewaert braid checks completed.")
