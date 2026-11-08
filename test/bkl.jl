@testset "BKL simple factors and native generators" begin
    @test bkl_atom(3, 1) == bkl_atom(1, 3)
    @test artin_letters(bkl_atom(3, 1)) == [2, 1, -2]
    @test artin_letters(bkl_atom(4, 1)) == [3, 2, 1, -2, -3]
    @test_throws ArgumentError bkl_atom(1, 1)
    @test_throws ArgumentError bkl_atom(1, 2, 0)

    catalan = Dict(2=>2, 3=>5, 4=>14, 5=>42, 6=>132)
    for n in 2:6
        D = BKLStructure(n)
        simples = bkl_simple_permutations(D)
        @test length(simples) == catalan[n]
        @test length(unique(simples)) == catalan[n]
        @test D.delta in simples
        @test identity_permutation(n) in simples
        @test all(is_bkl_simple, simples)
        for p in simples
            atoms = bkl_letters(p)
            @test permutation(BKLSimpleBraid(p)) == p
            @test SymmetricGroup(n)(artin_letters(atoms)) == p
            @test absolute_prefix(p, D.delta)
        end
    end

    D = BKLStructure(4)
    # The original BKL paper lists these 14 canonical factors in B4. Indices
    # are sorted only in storage; a31 and a13 denote the same band generator.
    source_factors = [
        Int[], [1,2], [2,3], [1,3], [3,4], [2,4], [1,4],
        [2,3,1,2], [3,4,2,3], [3,4,1,3], [3,4,2,1],
        [2,4,1,2], [1,4,2,3], [3,4,2,3,1,2],
    ]
    band(word) = BKLAtom[bkl_atom(word[k], word[k+1]) for k in 1:2:length(word)-1]
    source_permutations = [SymmetricGroup(4)(artin_letters(band(word))) for word in source_factors]
    @test all(p -> p in D.simples, source_permutations)
    @test Set(source_permutations) == Set(D.simples)

    for ts in bkl_generators(D)
        @test ts.i < ts.j
        @test absolute_length(SymmetricGroup(4)(artin_letters(ts))) == 1
    end

    @test Braid(δ(D)) == D.group([3, 2, 1])
    @test artin_letters(δ(D)) == [3, 2, 1]
    @test bkl_flip(D, transposition(4, 1)) == transposition(4, 2)
    @test bkl_flip(D, transposition(4, 3)) == transposition(4, 1, 4)
    @test Braid(bkl_normal_form(D, Int[])) == one(D.group)

    for a in D.simples
        @test SymmetricGroup(4)(artin_letters(bkl_letters(bkl_right_complement(D, a)))) ==
              inv(a) * D.delta
    end
end

@testset "BKL absolute gcd and local sliding" begin
    D = BKLStructure(4)
    simples = D.simples
    for a in simples
        complement = bkl_right_complement(D, a)
        @test absolute_prefix(complement, D.delta)
        @test a * complement == D.delta
    end

    for a in simples, b in simples
        g = bkl_left_gcd(D, a, b)
        @test absolute_prefix(g, a)
        @test absolute_prefix(g, b)
        for q in simples
            if absolute_prefix(q, a) && absolute_prefix(q, b)
                @test absolute_prefix(q, g)
            end
        end

        slid = bkl_local_sliding(D, a, b)
        @test slid.left * slid.right == a * b
        @test is_bkl_simple(slid.left)
        @test is_bkl_simple(slid.right)
        @test is_bkl_left_weighted(D, slid.left, slid.right)
    end
end
