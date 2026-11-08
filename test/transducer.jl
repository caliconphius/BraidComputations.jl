@testset "Chapter 9 transducer" begin
    B = BraidGroup(4)
    T = transducer(B)
    Δperm = longest_permutation(4)

    # Exhaust all states for n=4. Check both the permutation identity and the
    # positive-length identity O(a,t) M(a,t) = at.
    for values in all_permutations(4), i in 1:3
        state = StandardBraid(Permutation(values))
        step = transducer_step(T, state, i)
        @test step.output.permutation * step.state.permutation ==
              state.permutation * transposition(4, i)
        @test inversion_number(step.output) + inversion_number(step.state) ==
              inversion_number(state) + 1
        @test B([letters(state.permutation); i]) ==
              B([letters(step.output); letters(step.state)])
    end

    # Extended output identity from a nonidentity initial state.
    initial = StandardBraid(Permutation([2, 3, 1, 4]))
    input = [2, 1, 3, 2, 2, 1]
    result = transduce(T, input; initial)
    output_letters = reduce(vcat, letters.(result.output); init=Int[])
    @test B([letters(initial); input]) ==
          B([output_letters; letters(result.state)])
    @test maximal_simple_tail(B, input) == transduce(B, input).state

    @test standard_factors(Δ(B)) == StandardBraid[]
    @test B([1, 2, 3, 1, 2, 1]) == Δ(B)
    @test B([1, 2, 1]) == B([2, 1, 2])
    @test letters(standard_delta(B)) == [1, 2, 3, 1, 2, 1]

    for word in (Int[], [1], [1, 1], [1, 2, 3, 1, 2, 1], [1, 2, 1, 3, 2])
        x = right_greedy_form(B, word)
        @test x isa Braid
        @test B(letters(x)) == x
        @test all(p -> !isone(p) && p != Δperm, x.factors)
        @test infimum(x) + canonical_length(x) == supremum(x)
    end

    @test_throws ArgumentError transducer_step(T, StandardBraid(one(SymmetricGroup(4))), -1)
    @test_throws ArgumentError right_greedy_factors(B, [-1])
end
