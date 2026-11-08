module BraidComputations

import Farey
import Farey: farey_word

include("Permutations.jl")
include("Braids.jl")
include("Transducer.jl")
include("LeftNormalForms.jl")
include("GarsideOperations.jl")
include("SlidingCircuits.jl")
include("Parabolic.jl")
include("CanonicalReduction.jl")
include("Visualize.jl")
include("BKL.jl")
include("BKLNormalizing.jl")
include("Words.jl")

export SymmetricGroup, Permutation, BraidGroup, Braid, StandardBraid
export BraidTransducer, transducer, transduce, transducer_step, maximal_simple_tail
export standard_braid, standard_braid_from_word, standard_delta, standard_factors
export right_greedy_form, right_greedy_factors, normal_form
export LeftNormalForm, left_normal_form, is_left_normal, is_left_weighted
export right_complement, left_complement, left_gcd_standard, local_sliding
export initial_factor, final_factor, preferred_prefix, is_rigid
export sliding_conjugator, cyclic_sliding, cyclic_sliding_with_conjugator
export cycling, decycling, transport, iterate_cyclic_sliding
export has_sliding_circuit, sliding_circuit_period, sliding_circuits
export conjugacy_witness, is_conjugate
export StandardParabolic, ParabolicSystem, standard_parabolics
export rank, puncture_count, punctures, central_element, z
export flip_parabolic, act_on_parabolic, act_on_system, adjacent
export parabolic_orbit, standard_reduction_simplices, maximal_standard_reduction_simplices
export CanonicalReductionSystem, canonical_reduction_system, crs
export delete_strands, is_full_twist, is_central, is_periodic, braid_component
export outer_crs_partition, nielsen_thurston_type, is_reducible, is_pseudo_anosov
export standard_system, base_conjugator
export braid_svg, crs_svg, outer_crs_svg
export save_braid_svg, save_crs_svg, save_outer_crs_svg
export BKLAtom, BKLStructure, BKLSimpleBraid, BKLNormalForm
export bkl_atom, bkl_normal_form, bkl_simple, bkl_generators
export bkl_simple_permutations, is_bkl_simple, is_bkl_left_normal, is_bkl_left_weighted
export bkl_left_gcd, bkl_right_complement, bkl_local_sliding, bkl_flip
export bkl_letters, bkl_factors, bkl_delta, δ, artin_letters
export cycles, absolute_length, absolute_prefix, ncycle_permutation
export delta, Δ, generators, nstrands, images, identity_permutation, transposition
export permutation, permutation_from_inversions, inversion_set, inversion_number
export weak_order_meet, complement, flip, is_simple_product, longest_permutation
export letters, farey_word, farey_letters, group, delta_power
export infimum, supremum, canonical_length

end # module BraidComputations
