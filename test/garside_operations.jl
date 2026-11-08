@testset "preferred prefixes, sliding, and transport" begin
    B = BraidGroup(4)
    x = B([1, 2, 3])
    α = B([1, -2, 1])

    @test preferred_prefix(one(B)) == StandardBraid(identity_permutation(4))
    @test preferred_prefix(Δ(B)) == StandardBraid(identity_permutation(4))
    @test final_factor(one(B)) == standard_delta(B)
    @test initial_factor(Δ(B)^3) == StandardBraid(identity_permutation(4))
    @test is_rigid(Δ(B)^3)
    @test cyclic_sliding(Δ(B)^3) == Δ(B)^3

    for word in ([1, 2, 3], [1, 2, 1, 3, -2, 1], [-1, 2, -3, 2, 1],
                 [1, 1, 2, 3, 1, 2, 1], [3, 2, 1])
        x = B(word)
        ix = inv(x)
        @test preferred_prefix(x) == preferred_prefix(ix)
        @test inv(cyclic_sliding(x)) == cyclic_sliding(ix)
        @test B(final_factor(x)) * B(initial_factor(ix)) == Δ(B)
        @test canonical_length(cyclic_sliding(x)) <= canonical_length(x)
        @test infimum(cyclic_sliding(x)) >= infimum(x)
        @test supremum(cyclic_sliding(x)) <= supremum(x)

        c = B(initial_factor(x))
        d = B(final_factor(x))
        @test cycling(x) == inv(c) * x * c
        @test decycling(x) == d * x * inv(d)
        @test canonical_length(cycling(x)) <= canonical_length(x)
        @test canonical_length(decycling(x)) <= canonical_length(x)
        D = Δ(B)
        τx = inv(D) * x * D
        @test cycling(τx) == inv(D) * cycling(x) * D
        @test cyclic_sliding(τx) == inv(D) * cyclic_sliding(x) * D

        β = transport(α, x)
        expected = cyclic_sliding(inv(α) * x * α)
        @test inv(β) * cyclic_sliding(x) * β == expected

        result = iterate_cyclic_sliding(x)
        @test inv(result.conjugator) * x * result.conjugator == result.element
        @test result.circuit_start <= result.steps
        circuit = result.circuit
        @test !isempty(circuit)
        for i in 1:length(circuit)-1
            @test cyclic_sliding(circuit[i]) == circuit[i+1]
        end
        @test cyclic_sliding(last(circuit)) == first(circuit)
        @test all(y -> nstrands(y) == 4, circuit)
    end

    @test_throws DimensionMismatch transport(BraidGroup(3)(1), x)
end
