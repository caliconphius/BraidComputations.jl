# Proposed architecture and decisions

The long-term package name and some API names remain provisional. The first
implementation milestone uses classical braid groups, `Farey.Monoids`, an
on-demand transducer, signed Δ-prefix right-greedy storage, and a checked
left-normal-form conversion. Left normal form supports local and cyclic
sliding while braid arithmetic continues to use the existing right-greedy core.

## 1. Scope and layering

Implement classical Artin braid groups `B_n` first. Organize the Garside
operations behind small internal interfaces so that spherical Artin–Tits or
other Garside groups can be considered later, but do not build a fully generic
Garside framework before the braid implementation is tested.

Proposed layers:

1. **Words and adapters** — signed generator words; conversion to and from
   `Farey.Monoids` words.
2. **Permutations and simple braids** — the `n!` simple elements, Coxeter
   length, complementation, products that remain simple, and `Δ`.
3. **Thurston transducer** — maximal-simple-tail states and outputs;
   right/left-greedy forms; mixed canonical form.
4. **Garside arithmetic** — left normal form `Δ^p x_1 … x_r`, infimum,
   supremum, canonical length, inverses, products, lattice operations, and
   preferred prefixes.
5. **Conjugacy** — cycling/decycling where useful, cyclic sliding, sliding
   circuits, conjugators/transport, and minimal simple conjugating elements.
6. **Reduction systems** — standard irreducible parabolic subgroups,
   multicurves, Algorithm 1, strand deletion, periodicity tests, and
   braid-specific Algorithm 3.
7. **Optional automatic-group layer** — explicit automata, multiplier
   acceptors/transducers, and the older conjugacy construction in book §9.4.

**Decision A:** Is braid-only the right first milestone, with generic
spherical Artin–Tits support deferred? The CRS paper's general Algorithms 1
and 2 require a substantially broader abstraction than its braid Algorithm 3.

## 2. Farey boundary, not Farey internal storage

Use `Farey.Monoids` (`FreeGen`, `MonoidWord`, `eachgen`) as the supported word
input/output API, consistent with Farey's current cleanup direction. Do not
build on the legacy `FreeGroupElement`/KnuthBendix representation.

Farey words are free words. They are a syntactic interchange format, not a
model of braid equality: multiplication in Farey performs free cancellation,
whereas braid multiplication here must normalize modulo the braid relations.
Consequently, conversion should occur at an explicit boundary such as

```julia
braid(B, signed_letters)
braid(B, farey_word)
farey_word(b)
```

The initial adapter can extract signed integer generator letters from
`eachgen`. If that requires relying on unstable field details of `FreeGen`,
add a small public signed-letter iterator to Farey instead of duplicating
its internals here.

**Decision B:** Should the package depend directly on `Farey.Monoids` as
proposed, and may we make a small upstream signed-letter API addition if
needed? Alternatively, should Farey support be a package extension behind a
plain `Vector{Int}` core API?

## 3. Internal representation

Recommended representation:

- An immutable braid-group/context object records `n` and any transducer data.
- A simple braid is represented by its permutation, not by one chosen word.
- The current general braids have a signed Δ-prefix right-greedy
  representation consisting of an infimum and an immutable sequence of
  proper simple factors; a distinct left-normal-form representation or
  conversion will be required for cyclic sliding.
- Raw signed words remain a separate input representation.
- Constructors validate strand counts and generator ranges; multiplication
  returns a normalized braid.

Canonical storage makes hashing and duplicate elimination in sliding-circuit
and reduction-system searches reliable. It avoids equality depending on the
original word. The cost is that every arithmetic operation normalizes; a lazy
representation can be considered later if profiling justifies it.

A minimal internal permutation type avoids another dependency. Any wrapping
of input vectors must copy them to preserve immutability and hash stability.

**Earlier Decision C:** The initial implementation uses right-greedy normal
form as its stored representation. Before cyclic sliding, decide whether to
convert the storage to left normal form or retain both forms with checked
conversions. Separate raw/cached braid objects remain flexible but complicate
equality and conjugacy search.

## 4. Transducer construction and caching

The book's transducer has simple braids/permutations as states. A transition
returns both the next state and the simple output left behind. A first
implementation can compute transitions on demand; a context-local cache can
be added without global mutable state.

Building the full table costs roughly `n! * (n-1)` transitions. This is useful
for experiments at small `n`, but should not be an unconditional startup cost
for large `n`. Normal forms can share the same transition code whether or not
an explicit automaton object is requested.

**Decision D:** Should explicit transducer/automaton objects and their accepted
languages be public deliverables, or is the immediate goal only transducer-
based word normalization? Also, are eager tables wanted for small `n`, or
should all caching be lazy?

## 5. Conventions to fix before coding

Use one documented convention across permutations, braid multiplication,
conjugation, curve/parabolic-subgroup actions, and displayed words. Proposed
choices:

- `σ_i`, for `1 ≤ i < n`, denotes the standard positive Artin generator;
  signed integer `-i` denotes its inverse in input/output words.
- Braid products concatenate in the written order, with permutation
  composition tested explicitly against that order.
- `x^c` means `inv(c) * x * c`, matching the existing Farey convention.
- For the CRS paper, `P^α = inv(α) * P * α`; its output `(C, c)` satisfies
  `C` standard for the sliding-circuit representative `γ` and `γ^c = α`.
- Chapter 9's maximal-tail transducer works on positive words. Inverse
  letters are handled by powers of `Δ`, not by incorrectly extending its
  positive transition table to signed input.
- Keep left/right normal form and left/right divisibility terminology
  explicit; the book's right-greedy form must not be silently called the CRS
  paper's left normal form without a checked conversion.

**Decision E:** Confirm the conjugation and action conventions above,
particularly whether geometric curve action or algebraic parabolic-subgroup
conjugation should be the primary public interface.

## 6. Conjugacy and CRS API

Return witnesses as well as booleans/normal forms. Useful result objects
include:

- normal-form data with conversion back to a generator word;
- a cyclic-sliding result carrying the representative, conjugator, and number
  of steps;
- a sliding-circuit collection carrying a conjugator for every element;
- standard reduction simplices as sets of parabolic subgroups;
- a CRS result `(C, c)` in exactly the paper's orientation.

Do not implement full sliding circuits as repeated cyclic sliding only: the
latter reaches one circuit but does not enumerate `SC(α)`. Algorithm 3 also
requires transport of conjugators and the minimal-simple-conjugator search.

**Decision F:** Is the target the complete braid CRS Algorithm 3, including
full `SC(α)` enumeration, or should the first release stop at reaching a
sliding circuit and computing standard reduction simplices?

## 7. Dependencies and validation

Keep runtime dependencies to Farey and the standard library initially.
Chevie already appears in local experiments and could serve as an optional
development-time oracle, but should not become a runtime dependency without
discussion.

Validation should combine:

- small-`n` exhaustive permutation/simple-braid checks;
- braid-relation and `Δ` conjugation identities;
- the transducer identity `output * final_state = state * input`;
- reconstruction and uniqueness of canonical forms;
- conjugation tests that verify every returned conjugating element;
- independent expected CRS examples: periodic/pseudo-Anosov empty systems and
  explicit reducible braids with nonempty systems;
- optional cross-checks against Chevie, clearly separated from hand-checked
  tests so an oracle is not treated as the definition of correctness.

**Decision G:** May Chevie be used in a separate test/development environment
for differential testing, or should validation be entirely independent of it?
