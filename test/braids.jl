@testset "standard braids and signed Garside arithmetic" begin
    B3 = BraidGroup(3)
    B4 = BraidGroup(4)

    @test standard_braid_from_word(B3, [1, 2, 1]) == standard_delta(B3)
    @test standard_braid_from_word(B3, Int[]) == one(StandardBraid(identity_permutation(3)))
    @test_throws ArgumentError standard_braid_from_word(B3, [1, 1])
    @test_throws ArgumentError standard_braid_from_word(B3, [-1])

    @test B3([1, 2, 1]) == B3([2, 1, 2])
    @test B4([1, 3]) == B4([3, 1])
    @test B3([1]) != B3([2])
    @test B3([1, 1]) != one(B3)
    @test one(B3) == B3(Int[])

    # Flip comes from Δb = flip(b)Δ.
    for i in 1:3
        @test Δ(B4) * B4(i) == B4(4-i) * Δ(B4)
        @test Δ(B4)^2 * B4(i) == B4(i) * Δ(B4)^2
    end

    # An explicit inverse example: Δσ₂⁻¹ is the complementary simple braid.
    @test B3([-2]) == inv(Δ(B3)) * standard_braid(B3, complement(transposition(3, 2)))
    @test delta_power(B3([-1])) == -1
    @test canonical_length(B3([-1])) == 1

    rng = MersenneTwister(0xbaced)
    for n in 2:4, trial in 1:15
        B = BraidGroup(n)
        S = SymmetricGroup(n)
        random_letters = () -> begin
            result = Int[]
            while length(result) < 12
                i = rand(rng, -n+1:n-1)
                i == 0 || push!(result, i)
            end
            result
        end
        x = B(random_letters())
        y = B(random_letters())
        z = B(random_letters())

        @test x * inv(x) == one(B)
        @test inv(x) * x == one(B)
        @test (x * y) * inv(y) == x
        @test (x * y) * z == x * (y * z)
        @test x^0 == one(B)
        @test x^1 == x
        @test x^(-1) == inv(x)
        @test x^2 == x*x
        @test x^y == inv(y)*x*y
        @test B(letters(x)) == x
        @test permutation(x) == S(letters(x))

        # Moving trailing Δ factors through flip must not affect the canonical
        # output or the displayed generator word's induced permutation.
        nf = normal_form(x)
        @test nf.delta == delta_power(x)
        @test nf.factors == standard_factors(x)
        @test all(f -> !isone(f) && f != standard_delta(B), nf.factors)
    end

    @test_throws DimensionMismatch B3(1) * B4(1)
    @test_throws ArgumentError B3(0)
    @test_throws ArgumentError B3(3)
end
