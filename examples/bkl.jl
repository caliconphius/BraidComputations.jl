using BraidComputations
using Farey.Monoids

D = BKLStructure(4)
B = D.group
S = SymmetricGroup(4)

# BKL equation (5): δ is descending in the Artin generators.
@assert artin_letters(δ(D)) == [3, 2, 1]
@assert length(bkl_simple_permutations(D)) == 14  # Catalan(4)

x = bkl_normal_form(D, [1, 2, 1, 2, -1, 3])
@show x
@show normal_form(x)
@show bkl_letters(x)
@show artin_letters(x)
@show permutation(x)

@assert is_bkl_left_normal(D, x)
@assert Braid(x) == B([1, 2, 1, 2, -1, 3])
@assert bkl_normal_form(D, bkl_letters(x)) == x
@assert x * inv(x) == one(D)

# Native band generators. a31 and a13 are accepted as the same generator;
# storage/display sorts their two indices.
native = bkl_normal_form(D, [bkl_atom(3, 1), bkl_atom(4, 3, -1)])
@show native
@show artin_letters(native)
@assert Braid(native) == B(artin_letters(native))

a = D.simples[4]
b = D.simples[8]
@show bkl_left_gcd(D, a, b)
@show bkl_local_sliding(D, a, b)

f, g, h = FreeGen("f g h")
w = f * g * inv(f) * h
from_farey = bkl_normal_form(D, w; generators=(f, g, h))
@assert Braid(from_farey) == B(w; generators=(f, g, h))

println("Finished BKL normal-form example.")
