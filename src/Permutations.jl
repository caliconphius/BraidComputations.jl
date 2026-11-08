"""
    SymmetricGroup(n)

The symmetric group on `1:n`. Elements are permutations, and their product is
composition from right to left: `p * q` is `p ∘ q`.
"""
struct SymmetricGroup
    n::Int
    function SymmetricGroup(n::Integer)
        n >= 1 || throw(ArgumentError("the symmetric group needs n ≥ 1"))
        new(Int(n))
    end
end

"""
    Permutation(images)

A permutation in one-line notation. `images[i]` is the image of `i`. The input
is copied into immutable storage.
"""
struct Permutation
    images::Tuple{Vararg{Int}}
    function Permutation(images::AbstractVector{<:Integer})
        values = Int.(images)
        isempty(values) && throw(ArgumentError("a permutation needs at least one point"))
        sort(values) == collect(1:length(values)) ||
            throw(ArgumentError("images must be a permutation of 1:n"))
        new(Tuple(values))
    end
end

nstrands(p::Permutation) = length(p.images)
nstrands(S::SymmetricGroup) = S.n
images(p::Permutation) = collect(p.images)
Base.:(==)(a::Permutation, b::Permutation) = a.images == b.images
Base.isequal(a::Permutation, b::Permutation) = a == b
Base.hash(p::Permutation, h::UInt) = hash(p.images, h)
Base.copy(p::Permutation) = p

identity_permutation(n::Integer) = Permutation(collect(1:n))
Base.one(S::SymmetricGroup) = identity_permutation(S.n)
Base.one(p::Permutation) = identity_permutation(nstrands(p))
Base.isone(p::Permutation) = all(i -> p.images[i] == i, 1:nstrands(p))

function transposition(n::Integer, i::Integer)
    1 <= i < n || throw(ArgumentError("a simple transposition needs 1 ≤ i < n"))
    values = collect(1:n)
    values[i], values[i + 1] = values[i + 1], values[i]
    Permutation(values)
end

function transposition(n::Integer, i::Integer, j::Integer)
    1 <= i < j <= n || throw(ArgumentError("a transposition needs 1 ≤ i < j ≤ n"))
    values = collect(1:n)
    values[i], values[j] = values[j], values[i]
    Permutation(values)
end

(S::SymmetricGroup)(i::Integer) = transposition(S.n, abs(i))
(S::SymmetricGroup)(p::Permutation) = begin
    nstrands(p) == S.n || throw(ArgumentError("permutation has the wrong size"))
    p
end

generators(S::SymmetricGroup) = [S(i) for i in 1:S.n-1]

function _check_same_size(a::Permutation, b::Permutation)
    nstrands(a) == nstrands(b) || throw(DimensionMismatch("permutations have different sizes"))
    nothing
end

"""Compose permutations, returning `a ∘ b`."""
function Base.:*(a::Permutation, b::Permutation)
    _check_same_size(a, b)
    Permutation([a.images[b.images[i]] for i in 1:nstrands(a)])
end

Base.inv(p::Permutation) = Permutation(sortperm(images(p)))

function Base.:^(p::Permutation, n::Integer)
    n < 0 && return inv(p)^(-n)
    Base.power_by_squaring(p, n)
end

function Base.show(io::IO, S::SymmetricGroup)
    print(io, "SymmetricGroup(", S.n, ")")
end

function Base.show(io::IO, p::Permutation)
    print(io, "Permutation(", images(p), ")")
end

"""The inversion set `{(i,j) : i<j and p(i)>p(j)}`."""
function inversion_set(p::Permutation)
    result = Set{Tuple{Int,Int}}()
    for i in 1:nstrands(p)-1, j in i+1:nstrands(p)
        p.images[i] > p.images[j] && push!(result, (i, j))
    end
    result
end

inversion_number(p::Permutation) = length(inversion_set(p))

"""
    weak_order_meet(a, b)

The meet of `a` and `b` in the weak order used in *Word Processing in Groups*,
Chapter 9. Its inversion set is computed by the closure formula preceding the
proof of Proposition 9.1.8, not just by intersecting inversion sets.
"""
function weak_order_meet(a::Permutation, b::Permutation)
    _check_same_size(a, b)
    candidates = intersect(inversion_set(a), inversion_set(b))
    n = nstrands(a)
    result = Set{Tuple{Int,Int}}()

    # Conditions for (i,k) only refer to strictly shorter pairs, so distance
    # order implements the induction in the book's formula.
    for distance in 1:n-1, i in 1:n-distance
        k = i + distance
        (i, k) in candidates || continue
        if all(j -> (i, j) in result || (j, k) in result, i+1:k-1)
            push!(result, (i, k))
        end
    end

    # The empty inversion set is the identity. Reconstruct the permutation by
    # repeatedly adjoining a left descent to the reduced word.
    return _permutation_from_inversions(result, n)
end

"""Reconstruct a permutation from a valid inversion set."""
function permutation_from_inversions(inversions::Set{Tuple{Int,Int}}, n::Integer)
    _permutation_from_inversions(inversions, n)
end

# An inversion set orders the *positions* by their permutation values: for
# i<j, inversion (i,j) means that j precedes i in that order. Sort positions by
# this relation and assign their ranks as values.
function _permutation_from_inversions(inversions::Set{Tuple{Int,Int}}, n::Integer)
    positions = sort(collect(1:n); lt=(a, b) -> begin
        a == b && return false
        lo, hi = minmax(a, b)
        inverted = (lo, hi) in inversions
        a == lo ? !inverted : inverted
    end)
    values = zeros(Int, n)
    for (rank, position) in enumerate(positions)
        values[position] = rank
    end
    p = Permutation(values)
    inversion_set(p) == inversions ||
        throw(ArgumentError("the supplied set is not a permutation inversion set"))
    p
end

longest_permutation(n::Integer) = Permutation(collect(n:-1:1))
longest_permutation(S::SymmetricGroup) = longest_permutation(S.n)

"""Complementation by the longest element, as in book §9.1."""
complement(p::Permutation) = longest_permutation(nstrands(p)) * p

"""The involution induced by conjugation by Δ: `σᵢ ↦ σₙ₋ᵢ`."""
flip(p::Permutation) = longest_permutation(nstrands(p)) * p * longest_permutation(nstrands(p))

"""Whether the positive product of two simple elements is again simple."""
function is_simple_product(a::Permutation, b::Permutation)
    _check_same_size(a, b)
    inversion_number(a * b) == inversion_number(a) + inversion_number(b)
end

"""
    letters(p)

A deterministic reduced word in the adjacent transpositions. This implements
the descent argument in Lemma 9.1.5: every letter removes one inversion.
"""
function letters(p::Permutation)
    result = Int[]
    current = p
    while !isone(current)
        i = findfirst(i -> current.images[i] > current.images[i+1], 1:nstrands(current)-1)
        pushfirst!(result, i)
        current = current * transposition(nstrands(current), i)
    end
    result
end
