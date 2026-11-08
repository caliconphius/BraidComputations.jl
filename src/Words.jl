_generator_id(i::Integer) = UInt64(i)
_generator_id(g) = UInt64(Farey.Monoids.id(g))

function _generator_ids(generators)
    ids = _generator_id.(collect(generators))
    all(!iszero, ids) || throw(ArgumentError("generator IDs must be nonzero"))
    length(unique(ids)) == length(ids) ||
        throw(ArgumentError("Farey generators must have distinct IDs"))
    ids
end

"""
    farey_letters(word; generators=nothing)

Convert a `Farey.Monoids` free word to signed integer generator positions.
With no `generators`, return the signed underlying Farey IDs. Otherwise,
`generators` gives the correspondence `σᵢ ↔ generators[i]`. Inverse letters
are negative. This is deliberately an adapter: free-word equality is not
braid-group equality.
"""
function farey_letters(word; generators=nothing)
    ids = generators === nothing ? nothing : _generator_ids(generators)
    result = Int[]
    for generator in Farey.Monoids.eachgen(word)
        id = UInt64(Farey.Monoids.id(generator))
        id == 0 && continue
        position = ids === nothing ? Int(id) : findfirst(==(id), ids)
        position === nothing &&
            throw(ArgumentError("Farey generator ID $id is absent from the supplied generator list"))
        push!(result, generator.inv ? -position : position)
    end
    result
end

"""
    farey_word(word; generators=nothing)
    farey_word(x; generators=nothing)

Construct a `Farey.Monoids` free word from signed generator positions or from
one of this package's normal forms. By default, position `i` becomes
`FreeGen(i)`. With `generators`, position `i` instead uses `generators[i]`.
"""
function farey_word(word::AbstractVector{<:Integer}; generators=nothing)
    ids = generators === nothing ? nothing : _generator_ids(generators)
    result = Farey.Monoids.FreeGen(0)
    for i in word
        i == 0 && continue
        index = abs(Int(i))
        id = ids === nothing ? index : begin
            index <= length(ids) ||
                throw(ArgumentError("generator position $index exceeds the supplied generator list"))
            ids[index]
        end
        result *= Farey.Monoids.FreeGen(id; p=(i < 0 ? -1 : 1))
    end
    result
end

farey_word(x::Union{Braid,LeftNormalForm,StandardBraid,Permutation}; generators=nothing) =
    farey_word(letters(x); generators)

(B::BraidGroup)(word::Farey.Monoids.AbstractMonoidElement; generators=nothing) =
    B(farey_letters(word; generators))
(S::SymmetricGroup)(word::Farey.Monoids.AbstractMonoidElement; generators=nothing) =
    _word_to_permutation(S, farey_letters(word; generators))

function _word_to_permutation(S::SymmetricGroup, word::AbstractVector{<:Integer})
    result = one(S)
    for i in word
        i != 0 && 1 <= abs(i) < S.n ||
            throw(ArgumentError("symmetric-group generators are ±1,…,±$(S.n-1)"))
        result *= S(abs(i))
    end
    result
end

(S::SymmetricGroup)(word::AbstractVector{<:Integer}) = _word_to_permutation(S, word)
