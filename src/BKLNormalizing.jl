function _bkl_initial_decomposition(D::BKLStructure, word::AbstractVector{<:Integer})
    n = nstrands(D)
    delta = 0
    factors = Permutation[]
    for i in word
        i != 0 && 1 <= abs(i) < n ||
            throw(ArgumentError("Artin generators on $n strands are ±1,…,±$(n-1)"))
        atom = transposition(n, abs(i))
        if i > 0
            push!(factors, atom)
        else
            # σ⁻¹ = (σ⁻¹δ)δ⁻¹. Move the new δ⁻¹ left through H
            # using Hδ⁻¹ = δ⁻¹τ⁻¹(H).
            push!(factors, bkl_right_complement(D, atom))
            factors = [bkl_flip(D, factor, -1) for factor in factors]
            delta -= 1
        end
    end
    delta, factors
end

"""Normalize a sequence of dual simples by Proposition 2.6 local slidings."""
function _bkl_normalize_simples(D::BKLStructure, delta::Integer,
                                input::Vector{Permutation})
    factors = Permutation[]
    d = Int(delta)
    for atom in input
        isone(atom) && continue
        push!(factors, atom)
        for i in length(factors)-1:-1:1
            slid = bkl_local_sliding(D, factors[i], factors[i+1])
            factors[i] = slid.left
            factors[i+1] = slid.right
        end

        # Proposition 2.6 allows only a leading δ or a trailing identity after
        # this backward pass. Leading δ factors are already at the left, so no
        # additional flip is required when absorbing them.
        while !isempty(factors) && first(factors) == D.delta
            popfirst!(factors)
            d += 1
        end
        !isempty(factors) && isone(first(factors)) && popfirst!(factors)
        !isempty(factors) && isone(last(factors)) && pop!(factors)

        any(==(D.delta), factors) &&
            throw(ErrorException("non-leading δ remained after dual local sliding"))
        any(isone, factors) &&
            throw(ErrorException("identity remained after dual local sliding"))
    end
    BKLNormalForm(D.group, d, factors)
end

function _bkl_from_artin(D::BKLStructure, word::AbstractVector{<:Integer})
    # Use the established Artin normalizer to guarantee a canonical signed
    # input word, then substitute simple inverses in the dual structure.
    x = D.group(word)
    delta, factors = _bkl_initial_decomposition(D, letters(x))
    _bkl_normalize_simples(D, delta, factors)
end

bkl_normal_form(D::BKLStructure, x::Braid) = begin
    nstrands(x) == nstrands(D) || throw(ArgumentError("braid has the wrong strand count"))
    _bkl_from_artin(D, letters(x))
end
bkl_normal_form(D::BKLStructure, word::AbstractVector{<:Integer}) =
    _bkl_from_artin(D, word)
bkl_normal_form(D::BKLStructure, atoms::AbstractVector{BKLAtom}) =
    _bkl_from_artin(D, artin_letters(atoms))
bkl_normal_form(D::BKLStructure, x::BKLNormalForm) = begin
    nstrands(x) == nstrands(D) || throw(ArgumentError("braid has the wrong strand count"))
    x
end

bkl_normal_form(D::BKLStructure, word::Farey.Monoids.AbstractMonoidElement;
                generators=nothing) = bkl_normal_form(D, D.group(word; generators))
bkl_normal_form(B::BraidGroup, word; kwargs...) = bkl_normal_form(BKLStructure(B), word; kwargs...)
bkl_normal_form(n::Integer, word; kwargs...) = bkl_normal_form(BKLStructure(n), word; kwargs...)
bkl_normal_form(x::BKLNormalForm) = x

Base.inv(x::BKLNormalForm) = bkl_normal_form(BraidGroup(nstrands(x)), inv(Braid(x)))
Base.one(x::BKLNormalForm) = BKLNormalForm(x.group, 0, Permutation[])
Base.one(D::BKLStructure) = BKLNormalForm(D.group, 0, Permutation[])

function Base.:*(a::BKLNormalForm, b::BKLNormalForm)
    nstrands(a) == nstrands(b) || throw(DimensionMismatch("braids have different strand counts"))
    bkl_normal_form(a.group, Braid(a) * Braid(b))
end

function Base.:^(x::BKLNormalForm, n::Integer)
    n < 0 && return inv(x)^(-n)
    Base.power_by_squaring(x, n)
end
