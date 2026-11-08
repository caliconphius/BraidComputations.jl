"""The right complement `∂(s) = s⁻¹Δ`, again a standard braid."""
right_complement(s::StandardBraid) =
    StandardBraid(inv(s.permutation) * longest_permutation(nstrands(s)))

"""The left complement `∂⁻¹(s) = Δs⁻¹`, again a standard braid."""
left_complement(s::StandardBraid) =
    StandardBraid(longest_permutation(nstrands(s)) * inv(s.permutation))

"""
    left_gcd_standard(a, b)

Greatest common positive prefix of two standard braids. Reversal maps prefix
order to the suffix order computed by `weak_order_meet`.
"""
function left_gcd_standard(a::StandardBraid, b::StandardBraid)
    nstrands(a) == nstrands(b) || throw(DimensionMismatch("braids have different strand counts"))
    reverse(StandardBraid(weak_order_meet(inv(a.permutation), inv(b.permutation))))
end

"""
    local_sliding(a, b)

Slide `s = ∂(a) ∧ b` from `b` to `a`, returning `(as, s⁻¹b)`. This is the
local operation in Gebhardt–González-Meneses, *The cyclic sliding operation in
Garside groups*, §2.2. The returned pair is left weighted.
"""
function local_sliding(a::StandardBraid, b::StandardBraid)
    s = left_gcd_standard(right_complement(a), b)
    left = StandardBraid(a.permutation * s.permutation)
    right = StandardBraid(inv(s.permutation) * b.permutation)
    (left=left, right=right)
end

function flip_power(p::Permutation, exponent::Integer)
    result = p
    for _ in 1:abs(exponent)
        result = flip(result)
    end
    # flip is an involution, so positive and negative powers coincide.
    result
end

"""The simple initial factor `ι(x) = τ^(-inf(x))(x₁)` (identity if `ℓ(x)=0`)."""
function initial_factor(x::Braid)
    nf = left_normal_form(x)
    isempty(nf.factors) && return StandardBraid(identity_permutation(nstrands(x)))
    StandardBraid(flip_power(first(nf.factors), -nf.delta_power))
end

"""The simple final factor `φ(x)` (or `Δ` if `ℓ(x)=0`)."""
function final_factor(x::Braid)
    nf = left_normal_form(x)
    isempty(nf.factors) ? standard_delta(x.group) : StandardBraid(last(nf.factors))
end

"""
    preferred_prefix(x)

Compute `ι(x) ∧ ι(x⁻¹) = ι(x) ∧ ∂(φ(x))`, following Definition 3.2 of the
cyclic-sliding paper. This avoids constructing a second complete inverse form.
"""
function preferred_prefix(x::Braid)
    canonical_length(x) == 0 && return StandardBraid(identity_permutation(nstrands(x)))
    left_gcd_standard(initial_factor(x), right_complement(final_factor(x)))
end

is_rigid(x::Braid) = isone(preferred_prefix(x))

"""The conjugating simple braid used by cyclic sliding."""
sliding_conjugator(x::Braid) = BraidGroup(nstrands(x))(preferred_prefix(x))

"""Cyclic sliding `s(x) = p(x)⁻¹ x p(x)`."""
function cyclic_sliding(x::Braid)
    conjugator = sliding_conjugator(x)
    inv(conjugator) * x * conjugator
end

function cyclic_sliding_with_conjugator(x::Braid)
    conjugator = sliding_conjugator(x)
    (element=inv(conjugator)*x*conjugator, conjugator=conjugator)
end

"""Elrifai–Morton cycling by the initial factor; trivial if `ℓ(x)=0`."""
function cycling(x::Braid)
    canonical_length(x) == 0 && return x
    conjugator = BraidGroup(nstrands(x))(initial_factor(x))
    inv(conjugator) * x * conjugator
end

"""Elrifai–Morton decycling by the inverse final factor; trivial if `ℓ(x)=0`."""
function decycling(x::Braid)
    canonical_length(x) == 0 && return x
    final = BraidGroup(nstrands(x))(final_factor(x))
    final * x * inv(final)
end

"""
    transport(α, x)

Transport `α` at `x` under cyclic sliding:
`p(x)⁻¹ α p(α⁻¹xα)`. The returned braid conjugates `s(x)` to `s(α⁻¹xα)`.
"""
function transport(α::Braid, x::Braid)
    nstrands(α) == nstrands(x) || throw(DimensionMismatch("braids have different strand counts"))
    p = sliding_conjugator(x)
    conjugated = inv(α) * x * α
    inv(p) * α * sliding_conjugator(conjugated)
end

"""
    iterate_cyclic_sliding(x)

Iterate until the first repetition. Returns the representative, a verified
conjugator from `x`, the index at which the circuit starts, and one circuit in
cyclic order. This is only the first step of the full sliding-circuit search;
it does not enumerate `SC(x)`.
"""
function iterate_cyclic_sliding(x::Braid)
    seen = Dict{Braid,Int}()
    elements = Braid[]
    conjugators = Braid[]
    current = x
    conjugator = one(x.group)
    while !haskey(seen, current)
        seen[current] = length(elements) + 1
        push!(elements, current)
        push!(conjugators, conjugator)
        step = cyclic_sliding_with_conjugator(current)
        current = step.element
        conjugator *= step.conjugator
    end
    circuit_start = seen[current]
    return (element=current,
            conjugator=conjugator,
            circuit_start=circuit_start,
            steps=length(elements),
            circuit=elements[circuit_start:end])
end
