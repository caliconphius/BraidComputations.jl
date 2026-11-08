# Development plan

Current status: the Farey boundary, symmetric-group arithmetic, standard
braids, on-demand transducer, right-greedy factors, signed Δ-prefix
right-greedy arithmetic, left normal form, local sliding, preferred prefixes,
cycling/decycling, cyclic sliding, transport, and iteration to one circuit are
implemented and tested. A small-`n` BKL normal-form implementation using
absolute-order permutations is also complete. Full `SC(x)` enumeration is
still pending, and BKL uses an `n!` construction rather than an optimized
non-crossing-partition implementation. Each
subsequent stage should land with documentation,
small explicit examples, regression tests, and source references. Do not start
CRS implementation by stubbing its difficult conjugacy subroutines.

## Stage 0 — settle contracts and sources

- Resolve Decisions A–G in `architecture.md`.
- Fix product, permutation, inverse, conjugation, and parabolic-action
  conventions in executable examples.
- Work through the required Gebhardt–González-Meneses algorithms and the
  reduction-system prerequisites now linked in `source-map.md`.
- Add an executable Farey adapter contract test using `FreeGen`, `MonoidWord`,
  and `eachgen`.

Exit criterion: conventions and public scope are agreed, and every algorithm
in the first implementation milestone has a complete local or linked source.

## Stage 1 — words, permutations, and simple braids (initial implementation complete)

- Implement signed-word validation and Farey conversion in both directions.
- Implement permutations, inversion length, descent information,
  complementation, and permutation-to-simple-braid conversion.
- Implement `Δ`, positive-word construction, and simple-product tests.
- Exhaustively check the simple-braid/permutation correspondence for small
  `n`, together with braid and commutation relations.

Exit criterion: braid relations, permutation composition, inversion length,
and simple products pass independently constructed small-`n` tests.

## Stage 2 — book transducer and canonical forms (right-greedy core complete)

- Implement the maximal-simple-tail transition and simple output from book
  §9.2, recording left/right divisibility conventions beside the formulas.
- Apply the transducer recursively to obtain right-greedy forms.
- Implement checked conversion to left-greedy/left normal form.
- Handle arbitrary signed words using `Δ` powers; implement inverse,
  multiplication, infimum, supremum, and canonical length.
- Add mixed canonical form and, if selected, explicit acceptor/transducer
  objects from §9.3.

Exit criterion: transducer invariants, canonical-form reconstruction,
uniqueness, and word-problem tests pass for random and adversarial signed
words at several strand counts.

## Stage 3 — Garside conjugacy toolkit (basic operations complete)

- Implement preferred prefixes and cyclic sliding from the CRS paper's
  Definitions 4–6 and the full Gebhardt–González-Meneses sources.
- Return both representatives and verified conjugators.
- Implement iteration to a sliding circuit with cycle detection.
- Implement minimal simple conjugating elements and full sliding-circuit
  enumeration with transport data.
- Decide separately whether book §9.4's older positive-conjugacy search is a
  required compatibility feature or a later optional addition.

Exit criterion: exhaustive small cases and targeted conjugate pairs verify
representatives, circuit membership, conjugators, and conjugacy decisions.

## Stage 4 — standard parabolic subgroups and Algorithm 1

- Represent standard parabolic subgroups by generator subsets; specialize
  irreducible braid subgroups to intervals.
- Implement their standard positive elements `z_P`, adjacency/commutation,
  conjugation, and standardness tests from Lemmas 27–28.
- Implement orbits, reduction simplices, adjacency graphs, and complete-
  subgraph enumeration without silently assuming polynomial worst-case time.
- Return maximal simplices as well as all simplices if both are useful.

Exit criterion: tests cover nested versus adjacent subgroups, rejected
nonstandard orbits, disconnected adjacency graphs, and maximal-simplex
selection.

## Stage 5 — braid canonical reduction systems

- For each sliding-circuit representative, find the maximal dimension of a
  standard reduction simplex.
- Select a maximizing representative and intersect its maximal standard
  reduction simplices.
- Implement the strand-deletion/restriction construction of §5.2 and the
  periodicity test used by Algorithm 3.
- Remove the appropriate braid orbit of nonessential curves and return the
  standard system together with the conjugator to the input braid.
- Add classification helpers only after their relationship to CRS output is
  tested.

Exit criterion: hand-checked periodic, reducible, and pseudo-Anosov examples;
conjugacy covariance; power/central-element properties where applicable; and
independent small-`n` comparisons all pass.

## Stage 6 — optional expansions

- General spherical Artin–Tits groups and general CRS Algorithms 1–2.
- Centralizer calculations from Franco–González-Meneses and
  González-Meneses–Wiest, if these are part of the desired project.
- Dual Garside structures.
- Full automatic-structure objects and multiplier automata.
- Performance work, cached transition tables, serialization, and benchmarks
  only after correctness tests are established.

## Testing and documentation policy

- Keep small hand-derived expected values alongside property tests.
- Use randomized tests with fixed seeds and shrinking examples where useful.
- Verify witnesses mathematically instead of merely checking their types.
- Mark Chevie comparisons as optional differential tests, not as proofs.
- Cite theorem/algorithm numbers at implementation sites, with a note when a
  formula has been converted between left and right conventions.
- Do not advertise general Artin–Tits or CRS support until all delegated
  subroutines are actually implemented and tested.
