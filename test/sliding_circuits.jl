@testset "complete sliding circuits and conjugacy" begin
    @test sliding_circuit_period(BraidGroup(4)([1, 2, 3])) == 0
    @test sliding_circuit_period(BraidGroup(4)([3, 1, 2])) == 1

    # Cyclic-sliding paper Proposition 5.1: |SC(δ)|=2^(n-2)-2, giving two
    # vertices for the reversed four-braid σ3σ2σ1.
    @test length(sliding_circuits(BraidGroup(3)([2, 1])).elements) == 2
    @test length(sliding_circuits(BraidGroup(4)([3, 2, 1])).elements) == 2

    for (n, word) in ((2, Int[]), (3, [1, 2, 1, 2]), (4, [1, 2, 1, 2]))
        x = BraidGroup(n)(word)
        sc = sliding_circuits(x)
        @test !isempty(sc.elements)
        for γ in sc.elements
            @test sliding_circuit_period(γ) > 0
            @test γ^sc.conjugators[γ] == x
        end
        for (v, s, w) in sc.edges
            @test v^s == w
            @test sliding_circuit_period(w) > 0
        end
    end

    for n in 2:4
        B = BraidGroup(n)
        x = B(-1)
        g = n == 2 ? B([1, 1, -1]) : B([2, 1, -1, 2])
        z = g * x * inv(g)
        c = conjugacy_witness(x, z)
        @test !isnothing(c)
        @test x^c == z
        @test is_conjugate(x, z)
        @test !is_conjugate(B(1), B([1, 1]))
    end
end
