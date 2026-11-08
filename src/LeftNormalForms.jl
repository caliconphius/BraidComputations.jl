"""
    LeftNormalForm

The standard Garside left normal form `Δ^k x₁ … xᵣ`: all factors are proper
standard braids and every adjacent pair is left weighted.
"""
struct LeftNormalForm
    group::BraidGroup
    delta_power::Int
    factors::Vector{Permutation}
    function LeftNormalForm(group::BraidGroup, delta_power::Integer,
                            factors::Vector{Permutation})
        all(p -> nstrands(p) == group.n, factors) ||
            throw(ArgumentError("factor has the wrong strand count"))
        new(group, Int(delta_power), copy(factors))
    end
end

nstrands(x::LeftNormalForm) = x.group.n
group(x::LeftNormalForm) = x.group
delta_power(x::LeftNormalForm) = x.delta_power
infimum(x::LeftNormalForm) = x.delta_power
supremum(x::LeftNormalForm) = x.delta_power + length(x.factors)
canonical_length(x::LeftNormalForm) = length(x.factors)
standard_factors(x::LeftNormalForm) = StandardBraid.(x.factors)
normal_form(x::LeftNormalForm) = (delta=delta_power(x), factors=standard_factors(x))
Base.copy(x::LeftNormalForm) = LeftNormalForm(x.group, x.delta_power, copy(x.factors))
Base.one(x::LeftNormalForm) = left_normal_form(one(x.group))
Base.isone(x::LeftNormalForm) = isempty(x.factors) && iszero(x.delta_power)

function Base.:(==)(a::LeftNormalForm, b::LeftNormalForm)
    nstrands(a) == nstrands(b) && a.delta_power == b.delta_power && a.factors == b.factors
end
Base.isequal(a::LeftNormalForm, b::LeftNormalForm) = a == b
Base.hash(x::LeftNormalForm, h::UInt) =
    hash((x.group.n, x.delta_power, Tuple(x.factors)), h)

function letters(x::LeftNormalForm)
    result = Int[]
    Δletters = letters(longest_permutation(nstrands(x)))
    if x.delta_power > 0
        for _ in 1:x.delta_power
            append!(result, Δletters)
        end
    else
        for _ in 1:-x.delta_power
            append!(result, -reverse(Δletters))
        end
    end
    for factor in x.factors
        append!(result, letters(factor))
    end
    result
end

permutation(x::LeftNormalForm) = permutation(Braid(x))

function Base.show(io::IO, x::LeftNormalForm)
    isone(x) && return print(io, "1")
    pieces = String[]
    d = x.delta_power
    d != 0 && push!(pieces, d == 1 ? "Δ" : "Δ^($d)")
    append!(pieces, [join(("σ$(i)" for i in letters(factor)), "·") for factor in x.factors])
    print(io, join(pieces, " · "))
end

"""The anti-automorphism reversing a word while preserving generator signs."""
Base.reverse(x::StandardBraid) = StandardBraid(inv(x.permutation))
Base.reverse(x::Braid) = x.group(reverse(letters(x)))
Base.reverse(x::LeftNormalForm) = left_normal_form(x.group(reverse(letters(x))))

"""
    left_normal_form(x)

Convert a braid to its left normal form. Reversal sends right-greedy form to
left-greedy form; the second conversion therefore requires only reversing the
factors and replacing each standard permutation by its inverse.
"""
function left_normal_form(x::Braid)
    reversed = reverse(x)
    reversed.delta_power == x.delta_power ||
        throw(ErrorException("word reversal changed the Δ exponent"))
    factors = [inv(factor) for factor in reverse(reversed.factors)]
    # reverse(Δ^k H) = reverse(H)Δ^k = Δ^k flip^k(reverse(H)).
    factors = [flip_power(factor, x.delta_power) for factor in factors]
    LeftNormalForm(x.group, x.delta_power, factors)
end

left_normal_form(x::LeftNormalForm) = x
left_normal_form(B::BraidGroup, word::AbstractVector{<:Integer}) =
    left_normal_form(B(word))
LeftNormalForm(x::Braid) = left_normal_form(x)

Braid(x::LeftNormalForm) = x.group(letters(x))
(B::BraidGroup)(x::LeftNormalForm) = begin
    nstrands(x) == B.n || throw(ArgumentError("braid has the wrong strand count"))
    B(letters(x))
end

Base.inv(x::LeftNormalForm) = left_normal_form(inv(Braid(x)))
is_left_weighted(a::StandardBraid, b::StandardBraid) =
    isone(left_gcd_standard(right_complement(a), b))

"""Check proper factors and every adjacent left-weighted pair."""
function is_left_normal(x::LeftNormalForm)
    Δperm = longest_permutation(nstrands(x))
    all(f -> !isone(f) && f != Δperm, x.factors) || return false
    all(i -> is_left_weighted(StandardBraid(x.factors[i]),
                              StandardBraid(x.factors[i+1])), 1:length(x.factors)-1)
end
