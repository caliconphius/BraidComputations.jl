@testset "BKL normal forms, products, and Farey input" begin
    # Both orientations of the ordinary braid relation have one dual form.
    D4 = BKLStructure(4)
    @test bkl_normal_form(D4, [1, 2, 1]) == bkl_normal_form(D4, [2, 1, 2])
    @test bkl_normal_form(D4, [1, 3]) == bkl_normal_form(D4, [3, 1])
    @test bkl_normal_form(D4, [3, 2, 1]) == δ(D4)
    @test infimum(δ(D4)) == 1
    @test canonical_length(δ(D4)) == 0

    # δ^n is the classical full twist Δ².
    full_twist = Braid(δ(D4)^4)
    @test full_twist == Δ(D4.group)^2

    rng = MersenneTwister(0xb17)
    for n in 2:5, trial in 1:5
        D = BKLStructure(n)
        B = D.group
        make_word() = begin
            word = Int[]
            while length(word) < 10
                i = rand(rng, -n+1:n-1)
                i == 0 || push!(word, i)
            end
            word
        end
        word = make_word()
        x = B(word)
        nf = bkl_normal_form(D, word)

        @test is_bkl_left_normal(D, nf)
        @test Braid(nf) == x
        @test B(nf) == x
        @test permutation(nf) == SymmetricGroup(n)(word)
        @test bkl_normal_form(D, artin_letters(nf)) == nf
        @test bkl_normal_form(D, bkl_letters(nf)) == nf
        @test bkl_normal_form(D, x) == nf
        @test bkl_normal_form(D, nf) === nf

        @test nf * inv(nf) == one(D)
        @test inv(nf) * nf == one(D)
        @test infimum(inv(nf)) == -supremum(nf)
        @test supremum(inv(nf)) == -infimum(nf)
        @test canonical_length(inv(nf)) == canonical_length(nf)
        @test nf^0 == one(D)
        @test nf^1 == nf
        @test nf^2 == nf * nf

        native = BKLAtom[]
        for _ in 1:8
            i, j = sort([rand(rng, 1:n), rand(rng, 1:n)])
            i == j && continue
            push!(native, bkl_atom(i, j, rand(rng, Bool) ? 1 : -1))
        end
        native_nf = bkl_normal_form(D, native)
        @test Braid(native_nf) == B(artin_letters(native))
        @test bkl_normal_form(D, bkl_letters(native_nf)) == native_nf
    end

    f, g, h = FreeGen("f g h")
    D = BKLStructure(4)
    w = f * g * inv(f) * h
    nf = bkl_normal_form(D, w; generators=(f, g, h))
    @test Braid(nf) == D.group(w; generators=(f, g, h))
    @test bkl_normal_form(D, farey_word(letters(D.group(w; generators=(f,g,h)));
                                         generators=(f,g,h));
                      generators=(f,g,h)) == nf

    @test_throws ArgumentError BKLStructure(2).group(2)
end
