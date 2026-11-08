using BraidComputations

B = BraidGroup(4)
s = B([1, 2, 3])

# Published example from Gebhardt–González-Meneses, Proposition 5.1.
nf = left_normal_form(s^2)
@assert letters.(nf.factors) == [[1, 2, 3, 1, 2], [3]]
@assert is_left_normal(nf)
@assert Braid(nf) == s^2

@show nf
@show preferred_prefix(s)
@assert letters(preferred_prefix(s)) == [1, 2]
@assert cyclic_sliding(s) == B([3, 1, 2])

# Local sliding sends a simple pair (a,b) to a left-weighted pair (as,s⁻¹b).
a = StandardBraid(Permutation([2, 3, 1, 4]))
b = StandardBraid(Permutation([3, 1, 2, 4]))
slid = local_sliding(a, b)
@show slid
@assert is_left_weighted(slid.left, slid.right)
@assert B(slid.left) * B(slid.right) == B(a) * B(b)

@show cycling(s) decycling(s)

step = cyclic_sliding_with_conjugator(s)
@assert inv(step.conjugator) * s * step.conjugator == step.element

circuit = iterate_cyclic_sliding(s)
@show circuit.circuit_start circuit.steps circuit.circuit
@assert inv(circuit.conjugator) * s * circuit.conjugator == circuit.element
@assert cyclic_sliding(last(circuit.circuit)) == first(circuit.circuit)

α = B([1, -2, 1])
β = transport(α, s)
@assert inv(β) * cyclic_sliding(s) * β ==
        cyclic_sliding(inv(α) * s * α)

println("Finished Garside operations example.")
