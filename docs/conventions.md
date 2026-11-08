# Conventions and the implemented normal form

## Generators and permutation products

On `n` strands, signed integer `i` represents `σ_i`; `-i` represents
`σ_i^{-1}`. Input vectors and Farey words are read from left to right.

A `Permutation([a₁,…,aₙ])` uses one-line notation: position `i` has image
`aᵢ`. Products compose as functions:

```text
(p * q)(i) = p(q(i)).
```

Consequently, `permutation(x*y) == permutation(x)*permutation(y)`, with the
transposition `τᵢ = (i,i+1)` as the image of `σᵢ`.

## Standard braids

A standard braid (also called simple, non-repeating, or a permutation braid)
is stored by its permutation. Its length is the permutation's inversion
number. `letters` returns a reduced adjacent-transposition word by removing a
right descent at each step. In this way, no preferred positive word is needed
to define equality of standard braids.

## Maximal-tail transducer

For a simple state `a` and a positive generator `tᵢ`, book §9.2 gives
`b = a ∧ complement(tᵢ)` in the weak order. The package returns

```text
output   = a b⁻¹,
newstate = b tᵢ,
```

so `output * newstate = a * tᵢ`. Both entries are standard braids. The meet
uses the inversion-set closure immediately preceding Proposition 9.1.8, not a
plain intersection of inversion sets.

Applying this recursively to the transducer output gives the right-greedy
factors of a positive braid. No `n! × (n-1)` table is allocated unless caching
is added later.

## Signed Δ-prefix right-greedy form

The current `Braid` representation is

```text
Δᵏ x₁ … xᵣ,
```

where `x₁ … xᵣ` is right-greedy and every `xᵢ` is a nontrivial standard braid
other than `Δ`. This is a signed Garside normal form, not yet the left normal
form used by the later sliding-circuit work.

To normalize a word containing `m` inverse generators, multiply on the left by
`Δ^{2m}`. Centrality of `Δ²` and positivity of `Δ(Δσᵢ⁻¹)` allow every inverse
letter to be replaced by a positive segment. After right-greedy normalization,
trailing `Δ` factors are moved to the prefix using

```text
HΔ = Δ flip(H),       flip(σᵢ) = σₙ₋ᵢ.
```

`infimum`, `supremum`, and `canonical_length` refer to this normal form.

Products currently concatenate the normalized signed words and rerun the
normalizer. A direct cascade multiplication algorithm can be introduced later
without changing this external representation.

## Left normal form and local operations

The standard left normal form is also written `Δᵏ x₁ … xᵣ`, but now each pair
`(xᵢ,xᵢ₊₁)` is left weighted:

```text
right_complement(xᵢ) left-gcd xᵢ₊₁ = 1,
right_complement(s) = s⁻¹Δ.
```

Conversion uses the word-reversing anti-automorphism. Reversal interchanges
right- and left-greedy forms. Because `reverse(ΔᵏH) = reverse(H)Δᵏ`, moving
`Δᵏ` back to the left also applies `flipᵏ` to the reversed factors. `Braid`
and `LeftNormalForm` are distinct types so these conventions cannot be
confused accidentally.

`local_sliding(a,b)` uses `s = right_complement(a) left-gcd b` and returns
`(as,s⁻¹b)`. `initial_factor`, `final_factor`, `preferred_prefix`, and
`cyclic_sliding` follow Definitions 2.4, 3.1, and 3.2 of the cyclic-sliding
paper. In particular, `preferred_prefix(x)` is
`initial_factor(x) left-gcd right_complement(final_factor(x))`, and its cyclic
sliding is `p⁻¹xp`. `transport(α,x)` returns
`p(x)⁻¹ α p(α⁻¹xα)`.

`iterate_cyclic_sliding` stops at the first repetition and returns a verified
conjugator. It does not enumerate the full set `SC(x)`.

## Birman–Ko–Lee normal forms

The BKL implementation follows the original band-generator convention in
Birman–Ko–Lee, *A new approach to the word and conjugacy problems in the
braid groups*, equations (4) and (5):

```text
a_{ts} = (σ_{t-1}…σ_{s+1}) σ_s (σ_{s+1}^{-1}…σ_{t-1}^{-1}),
δ = σ_{n-1}…σ_1.
```

`BKLAtom(i,j)` sorts `i,j`, so either index order is accepted and both are
displayed using `i<j`. This matches the unordered-pair convention used by
Bessis–Digne–Michel while preserving the original band-generator formula.

The simple set consists of permutations `p` below the δ permutation in
absolute order. With `c(p)` the number of cycles, including fixed points,
`absolute_length(p) = n - c(p)`, and the condition is

```text
absolute_length(p) + absolute_length(inv(p) * δperm) = n-1.
```

A `BKLStructure` checks all `n!` permutations once and retains these simple
factors. There are Catalan(n) of them; in particular, the implementation
reproduces all 14 canonical factors printed by BKL in `B₄`. This establishes
a clear small-`n` implementation, not an optimized large-`n` algorithm.
Non-crossing-partition ranking or Kreweras-complement methods can replace the
enumeration later.

`bkl_left_gcd` computes the greatest common absolute-order prefix. Dual local
sliding and left normal forms use `s⁻¹δ` as the right complement and
`τ^k(p) = δ^(-k)pδ^k` as rotation. `bkl_letters` provides native signed band
generator output, while `artin_letters` provides the input to the existing
classical braid normalizer. Every BKL normal form is checked in tests against
classical `Braid` equality.

## Farey boundary

Farey generator IDs come from the generator names: for example, `FreeGen("f")`
has ID 6. Supply a correspondence explicitly when using named generators:

```julia
x = B(f*g*inv(f); generators=(f,g,h))
w = farey_word(x; generators=(f,g,h))
```

With `FreeGen(1), FreeGen(2), …`, or when raw numeric IDs are intended, the
`generators` keyword can be omitted. `farey_word` extends Farey's existing
function of that name for the new braid and permutation inputs; the original
rational/Farey-word methods remain unchanged.
