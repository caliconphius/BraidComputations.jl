@testset "Farey.Monoids word boundary" begin
    f, g, h = FreeGen("f g h")
    B = BraidGroup(4)
    S = SymmetricGroup(4)

    w = f^2 * g * inv(f)
    @test farey_letters(w) == [6, 6, 7, -6]
    @test farey_letters(w; generators=(f, g, h)) == [1, 1, 2, -1]
    @test B(w; generators=(f, g, h)) == B([1, 1, 2, -1])
    @test S(w; generators=(f, g, h)) == S([1, 1, 2, -1])

    x = B(w * h * inv(g); generators=(f, g, h))
    @test B(farey_word(x)) == x
    @test B(farey_word(x; generators=(f, g, h)); generators=(f, g, h)) == x
    @test farey_letters(farey_word(x)) == farey_letters(farey_word(letters(x)))
    @test S(farey_word(permutation(x))) == permutation(x)

    empty_word = farey_word(Int[])
    @test isempty(farey_letters(empty_word))
    @test B(empty_word) == one(B)
    @test S(empty_word) == one(S)

    @test letters(standard_braid(B, transposition(4, 2))) == [2]
    @test B(farey_word([1, 2, 1])) == B(farey_word([2, 1, 2]))
end
