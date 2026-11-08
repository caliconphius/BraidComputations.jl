using BraidComputations
using Farey.Monoids

# f, g, h will represent σ₁, σ₂, σ₃ respectively.
f, g, h = FreeGen("f g h")
B = BraidGroup(4)
S = SymmetricGroup(4)

x = B(f * g * f; generators=(f, g, h))
y = B(g * f * g; generators=(f, g, h))

@assert x == y                       # the length-three braid relation
@assert B(f * h; generators=(f, g, h)) ==
        B(h * f; generators=(f, g, h)) # distant generators commute

product = x * inv(B(h; generators=(f, g, h)))
@show product
@show normal_form(product)
@show letters(product)               # signed Artin-generator word
@show permutation(product)            # induced permutation
@show farey_word(product; generators=(f, g, h))

@assert B(farey_word(product; generators=(f, g, h));
          generators=(f, g, h)) == product
@assert product * inv(product) == one(B)

simple = standard_braid(B, Permutation([4, 3, 2, 1]))
@assert simple == standard_delta(B)
@assert B(simple) == Δ(B)

T = transducer(B)
step = transducer_step(T, StandardBraid(identity_permutation(4)), 1)
@show step.output step.state

# Symmetric-group words use the same generator mapping, but τᵢ = τᵢ⁻¹.
w = f * g * inv(f)
@show S(w; generators=(f, g, h))
@show letters(S(w; generators=(f, g, h))) # a reduced Coxeter word

println("Finished BraidComputations getting-started example.")
