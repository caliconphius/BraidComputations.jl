@testset "parabolic subgroups and standard reduction simplices" begin
    B = BraidGroup(4)
    P = StandardParabolic(1, 1, 4)
    @test puncture_count(P) == 2
    @test central_element(P) == B(1)
    Q = StandardParabolic(1, 2, 4)
    @test central_element(Q) == B([1, 2, 1, 1, 2, 1])
    @test act_on_parabolic(P, Δ(B)) == StandardParabolic(3, 3, 4)
    @test act_on_parabolic(P, B(1)) == P
    @test act_on_parabolic(Q, B(1)) == Q

    # Deleting strands tracks positions, not merely original labels.
    @test delete_strands(B([1, 2, 1, -3, 2]), [2, 3]) == [-1, 1, 1]
    @test delete_strands(B(1), [1]) == Int[]

    for x in (one(B), B(1), B([1, 2, 3]), B([1, 2, 1, 1, 2, 1]))
        for S in maximal_standard_reduction_simplices(x)
            for P in S, Q in S
                @test adjacent(P, Q)
            end
            @test !any(T -> length(S) < length(T) && S.curves ⊆ T.curves,
                       maximal_standard_reduction_simplices(x))
        end
    end
end

@testset "braid canonical reduction systems" begin
    @test outer_crs_partition(BraidGroup(4)(Int[])) == [1, 1, 1, 1]
    @test nielsen_thurston_type(BraidGroup(4)(Int[])) == :periodic
    @test is_central(BraidGroup(4)(Int[]))

    for n in 3:4
        B = BraidGroup(n)
        D = canonical_reduction_system(B(1))
        @test D.standard == ParabolicSystem([StandardParabolic(1, 1, n)])
        @test act_on_system(D.standard, D.conjugator) == D.standard
        @test outer_crs_partition(D) == [2; ones(Int, n-2)]
        @test nielsen_thurston_type(B(1)) == :reducible

        full_twist = Δ(B)^2
        @test isempty(canonical_reduction_system(full_twist))
        @test nielsen_thurston_type(full_twist) == :periodic
    end

    B = BraidGroup(4)
    @test isempty(canonical_reduction_system(B([3, 2, 1])))
    @test nielsen_thurston_type(B([3, 2, 1])) == :periodic

    twist3 = B([1, 2, 1, 1, 2, 1])
    reduction3 = canonical_reduction_system(twist3)
    @test reduction3.standard == ParabolicSystem([StandardParabolic(1, 2, 4)])
    @test outer_crs_partition(reduction3) == [3, 1]

    split = B(1) * B(3)
    @test canonical_reduction_system(split).standard ==
          ParabolicSystem([StandardParabolic(1, 1, 4), StandardParabolic(3, 3, 4)])
    @test outer_crs_partition(split) == [2, 2]

    # A canonical length-one pseudo-Anosov three-braid used in small examples.
    pa = BraidGroup(3)([1, -2])
    @test isempty(canonical_reduction_system(pa))
    @test nielsen_thurston_type(pa) == :pseudo_anosov
    @test is_pseudo_anosov(pa)

    # Conjugation covariance for a standard reducing curve.
    x = B([1, 2, 1, 1, 2, 1])
    c = B([2, 3, -2, 1])
    y = x^c
    result = canonical_reduction_system(y)
    @test !isempty(result)
    @test !isempty(result.standard)
    @test outer_crs_partition(result) == outer_crs_partition(x)
    @test nielsen_thurston_type(y) == :reducible

    # Periodicity preserves (possibly empty) outer partitions.
    @test outer_crs_partition(y^3) == outer_crs_partition(y)
end
