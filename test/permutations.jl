@testset "symmetric groups and reduced words" begin
    S = SymmetricGroup(4)
    @test S([1, 2, 1]) == Permutation([3, 2, 1, 4])
    @test S([1, 2, 1]) * S([1, 2, 1]) == one(S)
    @test S([-1]) == S(1)
    @test length(generators(S)) == 3

    for values in all_permutations(4)
        p = Permutation(values)
        word = letters(p)
        @test S(word) == p
        @test length(word) == inversion_number(p)
        @test p * inv(p) == one(p)
        @test inv(p) * p == one(p)
        @test p^2 == p*p
        @test permutation_from_inversions(inversion_set(p), 4) == p
        @test flip(flip(p)) == p
        @test Permutation(images(p)) == p
        @test hash(Permutation(images(p))) == hash(p)
    end

    # Meet in the book's weak order: inversion sets give an independent
    # greatest-lower-bound test over every pair and every possible lower bound.
    permutations = Permutation.(all_permutations(4))
    for a in permutations, b in permutations
        m = weak_order_meet(a, b)
        rm = inversion_set(m)
        @test rm ⊆ inversion_set(a)
        @test rm ⊆ inversion_set(b)
        for c in permutations
            rc = inversion_set(c)
            if rc ⊆ inversion_set(a) && rc ⊆ inversion_set(b)
                @test rc ⊆ rm
            end
        end
    end

    @test flip(transposition(4, 1)) == transposition(4, 3)
    @test flip(transposition(4, 2)) == transposition(4, 2)
    @test is_simple_product(transposition(3, 1), transposition(3, 2))
    @test !is_simple_product(transposition(3, 1), transposition(3, 1))
    @test_throws ArgumentError Permutation([1, 1, 2])
    @test_throws ArgumentError S([4])
end
