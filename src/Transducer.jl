"""The maximal-simple-tail transducer from book §9.2."""
struct BraidTransducer
    group::BraidGroup
end

transducer(B::BraidGroup) = BraidTransducer(B)
transducer(n::Integer) = transducer(BraidGroup(n))
nstrands(T::BraidTransducer) = T.group.n

"""
    transducer_step(T, state, i)

One positive transition. Returns `(output, state)` with standard-braid
entries satisfying

    output * new_state = old_state * σᵢ,

with `new_state` the maximal simple tail after reading `σᵢ`. The transition
uses the weak-order meet from book §9.2.
"""
function transducer_step(T::BraidTransducer, state::StandardBraid, i::Integer)
    B = T.group
    _check_letter(B, i)
    i > 0 || throw(ArgumentError("the maximal-tail transducer accepts positive letters only"))

    a = state.permutation
    τ = transposition(B.n, i)
    b = weak_order_meet(a, complement(τ))
    output = StandardBraid(a * inv(b))
    new_state = StandardBraid(b * τ)
    (output=output, state=new_state)
end

"""
    transduce(T, word; initial=identity)

Run the transducer on a positive word. The output is a vector of standard
factors; concatenating their words gives the transducer's positive output.
Returns `(output, state)`. It computes on demand rather than allocating the
`n! × (n-1)` transition table.
"""
function transduce(T::BraidTransducer, word::AbstractVector{<:Integer};
                   initial=StandardBraid(identity_permutation(T.group.n)))
    _check_positive_word(T.group, word)
    nstrands(initial) == T.group.n ||
        throw(ArgumentError("initial state has the wrong strand count"))
    output = Permutation[]
    state = initial
    for i in word
        step = transducer_step(T, state, i)
        isone(step.output) || push!(output, step.output.permutation)
        state = step.state
    end
    (output=output, state=state)
end

transduce(B::BraidGroup, word::AbstractVector{<:Integer}; kwargs...) =
    transduce(transducer(B), word; kwargs...)

function maximal_simple_tail(B::BraidGroup, word::AbstractVector{<:Integer})
    transduce(B, word).state
end

function _flatten_standard_words(factors)
    result = Int[]
    for factor in factors
        append!(result, letters(factor))
    end
    result
end

"""
    right_greedy_factors(B, positive_word)

Apply Theorem 9.2.2 recursively to find the right-greedy standard factors of
a positive braid. Factors are returned in their written left-to-right order.
"""
function _right_greedy_permutations(B::BraidGroup, word::AbstractVector{<:Integer})
    _check_positive_word(B, word)
    remaining = Int.(word)
    factors = Permutation[]
    while !isempty(remaining)
        previous_length = length(remaining)
        result = transduce(B, remaining)
        state = result.state
        isone(state) && throw(ErrorException("nonempty positive word has trivial maximal tail"))
        push!(factors, state.permutation)
        remaining = _flatten_standard_words(result.output)
        isempty(remaining) || length(remaining) < previous_length ||
            throw(ErrorException("transducer recursion failed to make progress"))
    end
    reverse!(factors)
    factors
end

right_greedy_factors(B::BraidGroup, word::AbstractVector{<:Integer}) =
    StandardBraid.(_right_greedy_permutations(B, word))

function _signed_replacement(B::BraidGroup, word::AbstractVector{<:Integer})
    _check_signed_word(B, word)
    positive = Int[]
    delta_adjustment = 0
    Δword = letters(longest_permutation(B.n))
    for i in word
        if i > 0
            push!(positive, i)
        else
            # Δ² is central and Δ²σᵢ⁻¹ = Δ(Δσᵢ⁻¹) is positive.
            # Thus Δ²w becomes positive after one such replacement per
            # inverse letter, exactly as in book §9.2.5.
            append!(positive, Δword)
            complement_i = StandardBraid(complement(transposition(B.n, abs(i))))
            append!(positive, letters(complement_i))
            delta_adjustment -= 2
        end
    end
    delta_adjustment, positive
end

"""Normalize a signed Artin-generator word into Δ-prefix right-greedy form."""
function _word_to_braid(B::BraidGroup, word::AbstractVector{<:Integer})
    adjustment, positive = _signed_replacement(B, word)
    factors = _right_greedy_permutations(B, positive)
    _delta_normal(B, adjustment, factors)
end

(B::BraidGroup)(word::AbstractVector{<:Integer}) = _word_to_braid(B, word)
(B::BraidGroup)(i::Integer) = B([i])

function right_greedy_form(B::BraidGroup, word::AbstractVector{<:Integer})
    B(word)
end

function right_greedy_factors(x::Braid)
    standard_factors(x)
end

function Base.:*(a::Braid, b::Braid)
    nstrands(a) == nstrands(b) || throw(DimensionMismatch("braids have different strand counts"))
    a.group(_product_word(a, b))
end
