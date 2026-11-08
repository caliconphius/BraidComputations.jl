"""A signed Birman–Ko–Lee generator `a_{i,j}` on `n` strands."""
struct BKLAtom
    i::Int
    j::Int
    sign::Int
    function BKLAtom(i::Integer, j::Integer, sign::Integer=1)
        i == j && throw(ArgumentError("BKL generators require distinct indices"))
        sign in (-1, 1) || throw(ArgumentError("a BKL exponent must be +1 or -1"))
        lo, hi = minmax(Int(i), Int(j))
        new(lo, hi, Int(sign))
    end
end

bkl_atom(i::Integer, j::Integer, sign::Integer=1) = BKLAtom(i, j, sign)
Base.show(io::IO, a::BKLAtom) =
    print(io, "a($(a.i),$(a.j))", a.sign == 1 ? "" : "⁻")

"""
    artin_letters(a::BKLAtom)

Write the stored `i<j` generator (the original paper's `a_{j,i}`) as
`(σ_{j-1}…σ_{i+1})σ_i(σ_{i+1}^{-1}…σ_{j-1}^{-1})`. Adjacent BKL generators
therefore reduce to the ordinary Artin generators.
"""
function artin_letters(a::BKLAtom)
    result = collect(a.j-1:-1:a.i+1)
    push!(result, a.sign * a.i)
    append!(result, -(a.i+1:a.j-1))
    result
end

"""A permutation which is below the standard `n`-cycle in absolute order."""
struct BKLSimpleBraid
    permutation::Permutation
    function BKLSimpleBraid(permutation::Permutation)
        n = nstrands(permutation)
        c = ncycle_permutation(n)
        absolute_length(permutation) + absolute_length(inv(permutation) * c) == n - 1 ||
            throw(ArgumentError("permutation is not below the n-cycle in absolute order"))
        new(permutation)
    end
end

permutation(x::BKLSimpleBraid) = x.permutation
nstrands(x::BKLSimpleBraid) = nstrands(x.permutation)
absolute_length(x::BKLSimpleBraid) = absolute_length(x.permutation)
Base.isone(x::BKLSimpleBraid) = isone(x.permutation)
Base.copy(x::BKLSimpleBraid) = x

function Base.show(io::IO, x::BKLSimpleBraid)
    isone(x) && return print(io, "1")
    print(io, join(string.(bkl_letters(x)), "·"))
end

"""
    BKLStructure(B)

The classical braid group equipped with its Birman–Ko–Lee/dual Garside
structure. The constructor enumerates the simple permutations once. This
initial implementation is intended for small and moderate `n`: enumeration
checks all `n!` permutations, although only `Catalan(n)` are retained.
"""
struct BKLStructure
    group::BraidGroup
    delta::Permutation
    simples::Vector{Permutation}
    function BKLStructure(group::BraidGroup)
        delta = ncycle_permutation(group.n)
        candidates = _all_permutation_values(group.n)
        simples = [Permutation(values) for values in candidates if _below_ncycle(Permutation(values))]
        new(group, delta, simples)
    end
end

BKLStructure(n::Integer) = BKLStructure(BraidGroup(n))
nstrands(D::BKLStructure) = D.group.n
group(D::BKLStructure) = D.group
bkl_simple_permutations(D::BKLStructure) = copy(D.simples)

# BKL equation (5): δ = σₙ₋₁…σ₁, an n-cycle in descending order.
ncycle_permutation(n::Integer) = Permutation([n; collect(1:n-1)])

"""Cycles of a permutation, including fixed points by default."""
function cycles(p::Permutation; include_fixed=true)
    unseen = Set(1:nstrands(p))
    result = Vector{Vector{Int}}()
    while !isempty(unseen)
        start = minimum(unseen)
        cycle = Int[]
        current = start
        while current in unseen
            push!(cycle, current)
            delete!(unseen, current)
            current = p.images[current]
        end
        (include_fixed || length(cycle) > 1) && push!(result, cycle)
    end
    result
end

absolute_length(p::Permutation) = nstrands(p) - length(cycles(p))
cycles(x::BKLSimpleBraid; kwargs...) = cycles(x.permutation; kwargs...)

function _below_ncycle(p::Permutation)
    c = ncycle_permutation(nstrands(p))
    absolute_length(p) + absolute_length(inv(p) * c) == nstrands(p) - 1
end

is_bkl_simple(p::Permutation) = _below_ncycle(p)

function _all_permutation_values(n::Integer)
    n == 1 && return [[1]]
    result = Vector{Vector{Int}}()
    for prefix in _all_permutation_values(n - 1), position in 1:n
        values = copy(prefix)
        insert!(values, position, n)
        push!(result, values)
    end
    result
end

bkl_simple(D::BKLStructure, p::Permutation) = begin
    nstrands(p) == nstrands(D) || throw(ArgumentError("permutation has the wrong strand count"))
    p in D.simples || throw(ArgumentError("permutation is not BKL-simple"))
    BKLSimpleBraid(p)
end

function bkl_generators(D::BKLStructure)
    [bkl_atom(i, j) for i in 1:nstrands(D)-1 for j in i+1:nstrands(D)]
end

"""The signed word in target `δ = σ₁…σₙ₋₁`, returned as a dual normal form."""
struct BKLNormalForm
    group::BraidGroup
    delta_power::Int
    factors::Vector{Permutation}
    function BKLNormalForm(group::BraidGroup, delta_power::Integer,
                           factors::Vector{Permutation})
        all(is_bkl_simple, factors) ||
            throw(ArgumentError("a factor is not BKL-simple"))
        new(group, Int(delta_power), copy(factors))
    end
end

nstrands(x::BKLNormalForm) = x.group.n
group(x::BKLNormalForm) = x.group
delta_power(x::BKLNormalForm) = x.delta_power
infimum(x::BKLNormalForm) = x.delta_power
supremum(x::BKLNormalForm) = x.delta_power + length(x.factors)
canonical_length(x::BKLNormalForm) = length(x.factors)
bkl_factors(x::BKLNormalForm) = BKLSimpleBraid.(x.factors)
normal_form(x::BKLNormalForm) = (delta=delta_power(x), factors=bkl_factors(x))
Base.copy(x::BKLNormalForm) = BKLNormalForm(x.group, x.delta_power, copy(x.factors))
Base.isone(x::BKLNormalForm) = isempty(x.factors) && iszero(x.delta_power)

function Base.:(==)(a::BKLNormalForm, b::BKLNormalForm)
    nstrands(a) == nstrands(b) && a.delta_power == b.delta_power && a.factors == b.factors
end
Base.isequal(a::BKLNormalForm, b::BKLNormalForm) = a == b
Base.hash(x::BKLNormalForm, h::UInt) =
    hash((x.group.n, x.delta_power, Tuple(x.factors)), h)

bkl_delta(D::BKLStructure) = BKLNormalForm(D.group, 1, Permutation[])
δ(D::BKLStructure) = bkl_delta(D)
δ(B::BraidGroup) = δ(BKLStructure(B))
δ(n::Integer) = δ(BKLStructure(n))

function Base.show(io::IO, x::BKLNormalForm)
    isone(x) && return print(io, "1")
    pieces = String[]
    d = x.delta_power
    d != 0 && push!(pieces, d == 1 ? "δ" : "δ^($d)")
    append!(pieces, [join(string.(bkl_letters(factor)), "·") for factor in x.factors])
    print(io, join(pieces, " · "))
end

"""Test `a ≤ b` in absolute order, meaning `a` is a positive prefix of `b`."""
function absolute_prefix(a::Permutation, b::Permutation)
    nstrands(a) == nstrands(b) || throw(DimensionMismatch("permutations have different sizes"))
    absolute_length(a) + absolute_length(inv(a) * b) == absolute_length(b)
end

function bkl_left_gcd(D::BKLStructure, a::Permutation, b::Permutation)
    nstrands(a) == nstrands(D) == nstrands(b) ||
        throw(DimensionMismatch("permutation has the wrong strand count"))
    candidates = [p for p in D.simples if absolute_prefix(p, a) && absolute_prefix(p, b)]
    maximum_length = maximum(absolute_length, candidates)
    candidates = filter(p -> absolute_length(p) == maximum_length, candidates)
    length(candidates) == 1 ||
        throw(ErrorException("absolute-order meet was not unique"))
    first(candidates)
end

bkl_right_complement(D::BKLStructure, a::Permutation) = inv(a) * D.delta

function bkl_local_sliding(D::BKLStructure, a::Permutation, b::Permutation)
    s = bkl_left_gcd(D, bkl_right_complement(D, a), b)
    (left=a*s, right=inv(s)*b)
end

function is_bkl_left_weighted(D::BKLStructure, a::Permutation, b::Permutation)
    isone(bkl_left_gcd(D, bkl_right_complement(D, a), b))
end

"""Conjugation `τ^k(p) = δ^(-k) p δ^k`, the dual rotation."""
function bkl_flip(D::BKLStructure, p::Permutation, exponent::Integer=1)
    c = D.delta^exponent
    inv(c) * p * c
end

"""A deterministic reduced word in all BKL reflections `a_{i,j}`."""
function bkl_letters(p::Permutation)
    is_bkl_simple(p) || throw(ArgumentError("permutation is not BKL-simple"))
    result = BKLAtom[]
    current = p
    while !isone(current)
        n = nstrands(current)
        found = false
        for i in 1:n-1, j in i+1:n
            t = transposition(n, i, j)
            next = t * current
            if absolute_length(next) == absolute_length(current) - 1 && is_bkl_simple(next)
                push!(result, bkl_atom(i, j))
                current = next
                found = true
                break
            end
        end
        found || throw(ErrorException("failed to find a reduced BKL reflection"))
    end
    result
end

bkl_letters(x::BKLSimpleBraid) = bkl_letters(x.permutation)

function bkl_letters(x::BKLNormalForm)
    result = BKLAtom[]
    if x.delta_power > 0
        for _ in 1:x.delta_power, i in reverse(1:nstrands(x)-1)
            push!(result, bkl_atom(i, i+1))
        end
    else
        for _ in 1:-x.delta_power, i in 1:nstrands(x)-1
            push!(result, bkl_atom(i, i+1, -1))
        end
    end
    for factor in x.factors
        append!(result, bkl_letters(factor))
    end
    result
end

artin_letters(atoms::AbstractVector{BKLAtom}) = reduce(vcat, artin_letters.(atoms); init=Int[])
artin_letters(x::Union{BKLNormalForm,BKLSimpleBraid}) = artin_letters(bkl_letters(x))

Braid(x::Union{BKLNormalForm,BKLSimpleBraid}) = BraidGroup(nstrands(x))(artin_letters(x))
(B::BraidGroup)(x::Union{BKLNormalForm,BKLSimpleBraid}) = begin
    nstrands(x) == B.n || throw(ArgumentError("braid has the wrong strand count"))
    B(artin_letters(x))
end

function permutation(x::BKLNormalForm)
    result = ncycle_permutation(nstrands(x))^x.delta_power
    for factor in x.factors
        result *= factor
    end
    result
end

function is_bkl_left_normal(D::BKLStructure, x::BKLNormalForm)
    nstrands(x) == nstrands(D) || return false
    all(f -> !isone(f) && f != D.delta, x.factors) || return false
    all(i -> is_bkl_left_weighted(D, x.factors[i], x.factors[i+1]), 1:length(x.factors)-1)
end
