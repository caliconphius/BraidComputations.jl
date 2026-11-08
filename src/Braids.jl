"""The classical Artin braid group on `n ≥ 2` strands."""
struct BraidGroup
    n::Int
    function BraidGroup(n::Integer)
        n >= 2 || throw(ArgumentError("the current braid implementation needs n ≥ 2"))
        new(Int(n))
    end
end

nstrands(B::BraidGroup) = B.n

"""
    StandardBraid(permutation)

A standard (simple, non-repeating, permutation) braid. Its permutation is an
immutable representative in the bijection of book Lemma 9.1.10.
"""
struct StandardBraid
    permutation::Permutation
end

permutation(x::StandardBraid) = x.permutation
nstrands(x::StandardBraid) = nstrands(x.permutation)
inversion_number(x::StandardBraid) = inversion_number(x.permutation)
letters(x::StandardBraid) = letters(x.permutation)
Base.isone(x::StandardBraid) = isone(x.permutation)
Base.one(x::StandardBraid) = StandardBraid(identity_permutation(nstrands(x)))

function Base.show(io::IO, x::StandardBraid)
    isone(x) && return print(io, "1")
    print(io, join(("σ$(i)" for i in letters(x)), "·"))
end

standard_braid(p::Permutation) = StandardBraid(p)
standard_braid(B::BraidGroup, p::Permutation) = begin
    nstrands(p) == B.n || throw(ArgumentError("permutation has the wrong strand count"))
    StandardBraid(p)
end

standard_braid(B::BraidGroup, word::AbstractVector{<:Integer}) =
    standard_braid_from_word(B, word)

# Concrete constructors avoid hiding the permutation in a word constructor.
function _permutation_for_positive_word(n::Int, word)
    p = identity_permutation(n)
    for i in word
        p *= transposition(n, i)
    end
    p
end

function standard_braid_from_word(B::BraidGroup, word::AbstractVector{<:Integer})
    _check_positive_word(B, word)
    p = _permutation_for_positive_word(B.n, word)
    inversion_number(p) == length(word) ||
        throw(ArgumentError("the positive word is not reduced, so it is not a standard braid"))
    StandardBraid(p)
end

"""
    Braid

A braid in Δ-prefix right-greedy normal form. The represented element is
`Δ^delta_power * x₁ * … * xᵣ`, where every stored factor is a proper,
nontrivial standard braid. Use `normal_form` to inspect the factors.
"""
struct Braid
    group::BraidGroup
    delta_power::Int
    factors::Vector{Permutation}
    function Braid(group::BraidGroup, delta_power::Integer, factors::Vector{Permutation})
        all(p -> nstrands(p) == group.n, factors) ||
            throw(ArgumentError("factor has the wrong strand count"))
        new(group, Int(delta_power), copy(factors))
    end
end

nstrands(x::Braid) = x.group.n
group(x::Braid) = x.group
delta_power(x::Braid) = x.delta_power
infimum(x::Braid) = x.delta_power
supremum(x::Braid) = x.delta_power + length(x.factors)
canonical_length(x::Braid) = length(x.factors)
standard_factors(x::Braid) = StandardBraid.(x.factors)
normal_form(x::Braid) = (delta=delta_power(x), factors=standard_factors(x))
right_greedy_form(x::Braid) = x

function Base.:(==)(a::Braid, b::Braid)
    nstrands(a) == nstrands(b) && a.delta_power == b.delta_power && a.factors == b.factors
end
Base.isequal(a::Braid, b::Braid) = a == b
Base.hash(x::Braid, h::UInt) = hash((x.group.n, x.delta_power, Tuple(x.factors)), h)
Base.copy(x::Braid) = Braid(x.group, x.delta_power, copy(x.factors))
Base.isone(x::Braid) = isempty(x.factors) && iszero(x.delta_power)
Base.one(B::BraidGroup) = _delta_normal(B, 0, Permutation[])
Base.one(x::Braid) = one(x.group)

function _delta_normal(B::BraidGroup, delta::Integer, factors::Vector{Permutation})
    result = filter(!isone, copy(factors))
    Δperm = longest_permutation(B.n)
    d = Int(delta)
    while !isempty(result) && last(result) == Δperm
        pop!(result)
        d += 1
        # HΔ = Δ flip(H). Moving trailing Δ factors to the prefix preserves
        # right-greediness because flip is a positive automorphism.
        result = flip.(result)
    end
    any(==(Δperm), result) &&
        throw(ErrorException("internal error: non-trailing Δ in right-greedy factors"))
    Braid(B, d, result)
end

function (B::BraidGroup)(x::StandardBraid)
    nstrands(x) == B.n || throw(ArgumentError("braid has the wrong strand count"))
    x.permutation == longest_permutation(B.n) ? _delta_normal(B, 1, Permutation[]) :
        _delta_normal(B, 0, [x.permutation])
end

(B::BraidGroup)(p::Permutation) = B(standard_braid(B, p))

function (B::BraidGroup)(factors::Vector{Permutation}, delta::Integer=0)
    _delta_normal(B, delta, factors)
end

"""The positive half-twist `Δ`, returned in braid normal form."""
delta(B::BraidGroup) = B(StandardBraid(longest_permutation(B.n)))
Δ(B::BraidGroup) = delta(B)
Δ(n::Integer) = Δ(BraidGroup(n))

function standard_delta(B::BraidGroup)
    StandardBraid(longest_permutation(B.n))
end

generators(B::BraidGroup) = [B(i) for i in 1:B.n-1]

function _check_letter(B::BraidGroup, i::Integer)
    i != 0 && 1 <= abs(i) < B.n ||
        throw(ArgumentError("braid generators on $(B.n) strands are ±1,…,±$(B.n-1)"))
    nothing
end

function _check_positive_word(B::BraidGroup, word)
    for i in word
        _check_letter(B, i)
        i > 0 || throw(ArgumentError("expected a positive braid word"))
    end
    nothing
end

function _check_signed_word(B::BraidGroup, word)
    for i in word
        _check_letter(B, i)
    end
    nothing
end

"""Convert a braid to a signed Artin-generator word from its normal form."""
function letters(x::Braid)
    result = Int[]
    Δletters = letters(longest_permutation(nstrands(x)))
    if x.delta_power > 0
        for _ in 1:x.delta_power
            append!(result, Δletters)
        end
    else
        inverse_Δletters = -reverse(Δletters)
        for _ in 1:-x.delta_power
            append!(result, inverse_Δletters)
        end
    end
    for factor in x.factors
        append!(result, letters(factor))
    end
    result
end

function permutation(x::Braid)
    result = identity_permutation(nstrands(x))
    Δperm = longest_permutation(nstrands(x))
    isodd(x.delta_power) && (result = Δperm)
    for factor in x.factors
        result *= factor
    end
    result
end

function Base.show(io::IO, B::BraidGroup)
    print(io, "BraidGroup(", B.n, ")")
end

function Base.show(io::IO, x::Braid)
    isone(x) && return print(io, "1")
    pieces = String[]
    d = x.delta_power
    d != 0 && push!(pieces, d == 1 ? "Δ" : "Δ^($d)")
    append!(pieces, [join(("σ$(i)" for i in letters(factor)), "·") for factor in x.factors])
    print(io, join(pieces, " · "))
end

function Base.:*(a::Braid, b::StandardBraid)
    a * a.group(b)
end
function Base.:*(a::StandardBraid, b::Braid)
    b.group(a) * b
end
function Base.:*(a::StandardBraid, b::StandardBraid)
    nstrands(a) == nstrands(b) || throw(DimensionMismatch("braids have different strand counts"))
    B = BraidGroup(nstrands(a))
    B(a) * B(b)
end

Base.inv(x::StandardBraid) = BraidGroup(nstrands(x))(-reverse(letters(x)))
Base.inv(x::Braid) = x.group(-reverse(letters(x)))

function Base.:^(x::Braid, n::Integer)
    n < 0 && return inv(x)^(-n)
    Base.power_by_squaring(x, n)
end

"""Conjugation in the same orientation as Farey: `x^c = inv(c) * x * c`."""
Base.:^(x::Braid, c::Braid) = inv(c) * x * c

function _product_word(a::Braid, b::Braid)
    result = letters(a)
    append!(result, letters(b))
    result
end
