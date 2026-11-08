"""True if repeated cyclic sliding returns to `x` in at least one step."""
function has_sliding_circuit(x::Braid)
    start = x
    current = cyclic_sliding(x)
    seen = Set{Braid}([start])
    while !(current in seen)
        push!(seen, current)
        current = cyclic_sliding(current)
    end
    current == start
end

"""The least positive period of `x` under cyclic sliding, or zero if none."""
function sliding_circuit_period(x::Braid)
    current = x
    seen = Dict{Braid,Int}()
    steps = 0
    while !haskey(seen, current)
        seen[current] = steps
        current = cyclic_sliding(current)
        steps += 1
    end
    current == x ? steps : 0
end

function _standard_braid_list(n::Integer)
    [StandardBraid(Permutation(values)) for values in _all_permutation_values(n)]
end

"""
    sliding_circuits(x)

Enumerate all of `SC(x)` by breadth-first search over conjugation by every
standard braid/symmetric permutation braid. Returns:

- `representative`: a vertex from `x`;
- `elements`: all vertices of the sliding-circuits graph `SCG(x)`;
- `conjugators`: a dictionary with `v^conjugators[v] == x`;
- `edges`: triples `(v, s, v^s)` with `s == 1` omitted.

This exhaustive implementation is intended for small braid indices. Every
standard positive element is considered, whereas Gebhardt–González-Meneses
Algorithm 2 constructs only indecomposable arrows. Restricting to periodic
vertices keeps the search finite, and every stored conjugator is independently
verified before returning.
"""
function sliding_circuits(x::Braid)
    seed = iterate_cyclic_sliding(x)
    representative = seed.element

    elements = Braid[representative]
    # The seed conjugator maps x to the representative. CRS needs the inverse
    # direction: γ^conjugators[γ] = x.
    conjugators = Dict{Braid,Braid}(representative => inv(seed.conjugator))
    edges = Tuple{Braid,Braid,Braid}[]
    simples = _standard_braid_list(nstrands(x))

    index = 1
    while index <= length(elements)
        v = elements[index]
        for simple in simples
            isone(simple) && continue
            s = BraidGroup(nstrands(x))(simple)
            w = inv(s) * v * s  # v^s
            has_sliding_circuit(w) || continue
            push!(edges, (v, s, w))
            haskey(conjugators, w) && continue
            push!(elements, w)
            conjugators[w] = inv(s) * conjugators[v]
        end
        index += 1
    end

    for (index, γ) in enumerate(elements)
        sliding_circuit_period(γ) > 0 ||
            throw(ErrorException("nonperiodic vertex emitted during SC enumeration"))
        inv(conjugators[γ]) * γ * conjugators[γ] == x ||
            throw(ErrorException("sliding-circuit conjugator verification failed"))
    end
    return (representative=representative,
            elements=elements,
            conjugators=conjugators,
            edges=edges)
end

"""Whether `x` and `y` are conjugate, returning a witness `x^c == y`."""
function conjugacy_witness(x::Braid, y::Braid)
    nstrands(x) == nstrands(y) || return nothing
    sx = sliding_circuits(x)
    sy = sliding_circuits(y)
    shared = findfirst(in(sy.elements), sx.elements)
    shared === nothing && return nothing
    γ = sx.elements[shared]
    c = inv(sx.conjugators[γ]) * sy.conjugators[γ]
    inv(c) * x * c == y ||
        throw(ErrorException("internal error constructing a conjugacy witness"))
    c
end

is_conjugate(x::Braid, y::Braid) = !isnothing(conjugacy_witness(x, y))
