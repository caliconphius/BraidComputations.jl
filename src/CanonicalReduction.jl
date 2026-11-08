"""
    delete_strands(x, keep)

Delete strands whose initial labels are absent from `keep`, tracking crossings
by their current positions. This implements the strand-deletion step in §5.2
of the CRS paper, including signs and reindexing as strands are removed.
"""
function delete_strands(x::Braid, keep::AbstractVector{<:Integer})
    n = nstrands(x)
    labels = unique(Int.(keep))
    all(i -> 1 <= i <= n, labels) ||
        throw(ArgumentError("strand labels must belong to 1:n"))
    active = falses(n)
    active[labels] .= true
    result = Int[]
    for i in letters(x)
        k = abs(i)
        if active[k] && active[k+1]
            rank = count(active[1:k-1]) + 1
            push!(result, i < 0 ? -rank : rank)
        end
        active[k], active[k+1] = active[k+1], active[k]
    end
    result
end

is_full_twist(x::Braid) = x == standard_delta(x.group) * standard_delta(x.group)
is_central(x::Braid) = canonical_length(x) == 0

"""A concrete finite-order/finite-order-modulo-center component test."""
function is_periodic(x::Braid)
    r = nstrands(x)
    r == 1 && return true
    return is_central(x^(r-1)) || is_central(x^r)
end

function _fixed_curve_exponent(P::StandardParabolic, α::Braid)
    z0 = central_element(P)
    current = z0
    m = 0
    while true
        m += 1
        current = inv(α) * current * α
        current == z0 && return m
        m > nstrands(α) &&
            throw(ArgumentError("the standard curve does not have a finite orbit of the expected length"))
    end
end

function _parent_curve(C::StandardParabolic, M::ParabolicSystem)
    parents = filter(P -> C != P && issubset(C, P), M.curves)
    isempty(parents) && return nothing
    first(sort(parents; by=puncture_count))
end

"""Build the component `α_{M,external}` from §5.2 of the CRS paper."""
function braid_component(α::Braid, M::ParabolicSystem,
                         external::Union{Nothing,StandardParabolic})
    n = nstrands(α)
    all(P -> nstrands(P) == n && P in M, external === nothing ? StandardParabolic[] : [external]) ||
        throw(ArgumentError("external boundary must belong to the reduction system"))

    m = external === nothing ? 1 : _fixed_curve_exponent(external, α)
    outer_punctures = external === nothing ? collect(1:n) : punctures(external)
    inner = external === nothing ? M.curves : filter(P -> P != external && issubset(P, external), M.curves)

    # Collapse each remaining inner curve to one arbitrary original puncture;
    # all the other punctures it encloses disappear from this component.
    collapsed = Int[]
    representatives = Int[]
    for P in inner
        append!(collapsed, punctures(P))
        push!(representatives, first(punctures(P)))
    end
    keep = [x for x in outer_punctures if !(x in collapsed)]
    append!(keep, representatives)
    unique!(keep)
    sort!(keep)

    word = delete_strands(α^m, keep)
    projected = BraidGroup(maximum((length(keep), 2)))(word)
    (braid=projected,
     exponent=m,
     punctures=keep,
     strands=length(keep))
end

"""The output `(standard_system, conjugator)` from CRS Algorithm 3."""
struct CanonicalReductionSystem
    standard::ParabolicSystem
    conjugator::Braid
end

standard_system(C::CanonicalReductionSystem) = C.standard
base_conjugator(C::CanonicalReductionSystem) = C.conjugator
Base.isempty(C::CanonicalReductionSystem) = isempty(C.standard)
Base.length(C::CanonicalReductionSystem) = length(C.standard)

function Base.show(io::IO, C::CanonicalReductionSystem)
    print(io, "CanonicalReductionSystem(", C.standard, ", conjugator=", C.conjugator, ")")
end

function _coorbit(C::StandardParabolic, α::Braid)
    orbit = parabolic_orbit(C, α)
    orbit === nothing &&
        throw(ErrorException("a standard reduction curve's orbit left the standard simplex"))
    orbit
end

"""
    canonical_reduction_system(x)

Compute the braid canonical reduction system using the braid-specific
Algorithm 3 of Cumplido–González-Meneses–Perego (2025). A maximizing member
of `SC(x)` supplies a standard system; nonessential curves are discarded one
orbit at a time using restriction and periodicity tests.

The full sliding-circuit search currently uses exhaustive standard conjugacy
and is intended for small braid indices. The returned `(system,c)` satisfies
`CRS(x) = system^c`, with `c` mapping the chosen circuit braid to `x`.
"""
function canonical_reduction_system(x::Braid)
    circuits = sliding_circuits(x)
    best = nothing
    best_simplices = ParabolicSystem[]
    best_dimension = -2

    for γ in circuits.elements
        simplices = maximal_standard_reduction_simplices(γ)
        d = maximum(dimension, simplices)
        if d > best_dimension
            best = γ
            best_simplices = simplices
            best_dimension = d
        end
    end

    best === nothing && throw(ErrorException("no sliding-circuit representative found"))
    M = reduce(intersect, best_simplices)
    initial_curves = copy(M.curves)

    for C in initial_curves
        C in M || continue
        parent = _parent_curve(C, M)
        reduced = ParabolicSystem(filter(P -> P != C, M.curves))
        component = braid_component(best, reduced, parent)
        if component.strands == 1 || is_periodic(component.braid)
            for Q in _coorbit(C, best)
                M = ParabolicSystem(filter(P -> P != Q, M.curves))
            end
        end
    end

    CanonicalReductionSystem(M, circuits.conjugators[best])
end

crs(x::Braid) = canonical_reduction_system(x)

function nielsen_thurston_type(x::Braid)
    (isone(x) || is_periodic(x)) && return :periodic
    isempty(canonical_reduction_system(x)) ? :pseudo_anosov : :reducible
end

is_reducible(x::Braid) = !isempty(canonical_reduction_system(x))
is_pseudo_anosov(x::Braid) = nielsen_thurston_type(x) == :pseudo_anosov

"""Partition of punctures for the outer canonical reduction system, up to conjugacy."""
function outer_crs_partition(system::ParabolicSystem)
    outers = filter(C -> !any(P -> C != P && issubset(C, P), system.curves), system.curves)
    covered = reduce(vcat, punctures.(outers); init=Int[])
    sizes = puncture_count.(outers)
    n = nstrands(system)
    append!(sizes, ones(Int, n-length(unique(covered))))
    sort!(sizes; rev=true)
end

outer_crs_partition(x::Braid) = begin
    system = canonical_reduction_system(x).standard
    isempty(system) ? ones(Int, nstrands(x)) : outer_crs_partition(system)
end
outer_crs_partition(C::CanonicalReductionSystem) = outer_crs_partition(C.standard)

function act_on_system(C::CanonicalReductionSystem)
    result = act_on_system(C.standard, C.conjugator)
    result === nothing &&
        throw(ErrorException("CRS standardizer could not return the system to standard form"))
    result
end
