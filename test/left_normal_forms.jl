@testset "left normal form" begin
    B = BraidGroup(4)
    s = B([1, 2, 3])

    # Gebhardt–González-Meneses, cyclic sliding paper, Proposition 5.1.
    square_form = left_normal_form(s^2)
    @test letters.(square_form.factors) == [[1, 2, 3, 1, 2], [3]]
    @test is_left_normal(square_form)
    @test Braid(square_form) == s^2
    @test letters(preferred_prefix(s)) == [1, 2]
    @test cyclic_sliding(s) == B([3, 1, 2])

    for word in (Int[], [1], [1, 2, 1], [1, 2, 1, 2, -1, 3],
                 [-1, -2, 1, 1, 2], [3, 2, 1, -3, -2, -1])
        x = B(word)
        nf = left_normal_form(x)
        @test is_left_normal(nf)
        @test Braid(nf) == x
        @test B(nf) == x
        @test B(letters(nf)) == x
        @test left_normal_form(nf) === nf
        @test infimum(nf) == infimum(x)
        @test supremum(nf) == supremum(x)
        @test canonical_length(nf) == canonical_length(x)
        @test permutation(nf) == permutation(x)
        @test infimum(inv(nf)) == -supremum(nf)
        @test supremum(inv(nf)) == -infimum(nf)
        @test canonical_length(inv(nf)) == canonical_length(nf)
        @test reverse(reverse(x)) == x
        @test B(farey_word(nf)) == x  # Farey performs free cancellation on output.
    end

    # Every adjacent pair in a result is left weighted, including identity and
    # Δ endpoints of local sliding.
    for a_values in all_permutations(4), b_values in all_permutations(4)
        a = StandardBraid(Permutation(a_values))
        b = StandardBraid(Permutation(b_values))
        result = local_sliding(a, b)
        @test result.left.permutation * result.right.permutation == a.permutation * b.permutation
        @test is_left_weighted(result.left, result.right)
        @test result.left == StandardBraid(result.left.permutation)
        @test result.right == StandardBraid(result.right.permutation)
    end

    for values in all_permutations(4)
        a = StandardBraid(Permutation(values))
        @test B(a) * B(right_complement(a)) == Δ(B)
        @test B(left_complement(a)) * B(a) == Δ(B)
    end
end
