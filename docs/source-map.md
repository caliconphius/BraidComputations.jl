# Source map

## Located local sources

### General braid reference

C. Kassel and V. Turaev, *Braid Groups*, Graduate Texts in Mathematics 247
(2008). Moved from the phone's Books library to the resources library during
project setup and verified as the 349-page book.

[Local PDF](</home/tcsiz/Downloads/resources/Kassel-Turaev(2008)BraidGroups.pdf>)

Use this for general braid-group background and a second account of Garside
constructions. The book and paper algorithms still need independent
left/right-convention checks before implementation.

### Word processing and transducers

D. B. A. Epstein, J. W. Cannon, D. F. Holt, S. V. F. Levy, M. S. Paterson and
W. P. Thurston, *Word Processing in Groups* (1992), Chapter 9.

[Local PDF](</home/tcsiz/Downloads/resources/Epstein-Cannon-Holt-Levy-Paterson-Thurston(1992)WordProcessingInGroups.pdf>)

- §9.1: braid presentation, permutations, inversion sets, simple braids,
  lattice operations and `Δ`.
- §9.2: maximal-tail automaton and transducer, right-greedy canonical form,
  left-greedy dual, and regularity.
- Theorem 9.2.2: output/final-state identity and maximal-tail property.
- §9.3: automatic structure, mixed canonical form and lattice consequences.
- §9.4: older transducer-aided conjugacy approach; separate from modern
  cyclic-sliding enumeration.

### Canonical reduction systems

M. Cumplido, J. González-Meneses and D. Perego, *Canonical Reduction Systems
in Artin–Tits groups of spherical type*, arXiv:2510.00713v2 (2025).

[Local PDF](</home/tcsiz/Downloads/resources/Cumplido-GonzalezMeneses-Perego(2025)CanonicalReductionSystemsInArtinTitsGroups-v2.pdf>) ·
[arXiv abstract](https://arxiv.org/abs/2510.00713)

- §2.3, Definitions 4–6: left normal form, preferred prefix, cyclic sliding
  and sliding circuits.
- Lemmas 27–28: standardness checks used for parabolic-subgroup actions.
- Algorithm 1: standard reduction simplices, including orbit and adjacency
  graph calculations.
- Algorithms 2–3: general Artin–Tits and improved braid CRS procedures.
- §5.2 and Algorithm 3: strand deletion, periodicity test, full sliding-
  circuit loop, maximal-standard-simplex intersection and output `(C, c)`.

Algorithm 3 explicitly calls for the whole sliding-circuit set and a
conjugating element for each representative. Reaching a single circuit by
iterated cyclic sliding is therefore only a preliminary subroutine.

### Conjugacy, efficient operations and reduction prerequisites

The implementation-critical references missing from the initial local
inventory have now been obtained from arXiv or the authors' public pages.
Version suffixes below identify the downloaded preprint versions; the years
follow their journal publications where applicable.

- V. Gebhardt and J. González-Meneses, *Solving the conjugacy problem in
  Garside groups by cyclic sliding* (2010). This self-contained paper gives
  Algorithms 1–3 for reaching a circuit, computing outgoing graph arrows and
  solving conjugacy decision/search, including transports, pullbacks and
  implementation-cost analysis. Definition 1.14 gives the implemented
  transport formula, and Algorithm 1 guides first-repetition iteration. [Local PDF](</home/tcsiz/Downloads/resources/Gebhardt-GonzalezMeneses(2010)SolvingConjugacyProblemInGarsideGroupsByCyclicSliding-v1.pdf>) · [arXiv:0809.0948](https://arxiv.org/abs/0809.0948)
- V. Gebhardt and J. González-Meneses, *The cyclic sliding operation in
  Garside groups* (2010). §2.2 supplies left normal form and local sliding;
  §3.1 supplies cycling, decycling, preferred prefixes and cyclic sliding;
  Proposition 5.1 supplies explicit braid test examples. [Local PDF](</home/tcsiz/Downloads/resources/Gebhardt-GonzalezMeneses(2010)CyclicSlidingOperationInGarsideGroups-v2.pdf>) · [arXiv:0808.1430](https://arxiv.org/abs/0808.1430)
- E.-K. Lee and S.-J. Lee, *A Garside-theoretic approach to the reducibility
  problem in braid groups* (2008). Standardizers and reduction calculations.
  [Local PDF](</home/tcsiz/Downloads/resources/Lee-Lee(2008)GarsideTheoreticApproachToReducibilityProblemInBraidGroups-v4.pdf>) · [arXiv:math/0506188](https://arxiv.org/abs/math/0506188)
- M. Cumplido, *On the minimal positive standardizer of a parabolic subgroup
  of an Artin–Tits group* (2019). Especially relevant to the general
  Artin–Tits algorithm. [Local PDF](</home/tcsiz/Downloads/resources/Cumplido(2019)MinimalPositiveStandardizerOfParabolicSubgroup-v2.pdf>) · [arXiv:1708.09310](https://arxiv.org/abs/1708.09310)
- M. Cumplido, V. Gebhardt, J. González-Meneses and B. Wiest, *On parabolic
  subgroups of Artin–Tits groups of spherical type* (2019). Parabolic
  subgroups and the complex used by the CRS algorithm. [Local PDF](</home/tcsiz/Downloads/resources/Cumplido-Gebhardt-GonzalezMeneses-Wiest(2019)ParabolicSubgroupsOfArtinTitsGroups-v3.pdf>) · [arXiv:1712.06727](https://arxiv.org/abs/1712.06727)
- J. González-Meneses, *Basic results on braid groups* (2011). Survey of
  Garside solutions to word and conjugacy problems and centralizer
  calculations. [Local PDF](</home/tcsiz/Downloads/resources/GonzalezMeneses(2011)BasicResultsOnBraidGroups-v1.pdf>) · [arXiv:1010.0321](https://arxiv.org/abs/1010.0321)
- J. González-Meneses and B. Wiest, *Reducible braids and Garside theory*
  (2011). Stabilized sliding circuits and reducibility/Nielsen–Thurston type.
  [Local PDF](</home/tcsiz/Downloads/resources/GonzalezMeneses-Wiest(2011)ReducibleBraidsAndGarsideTheory-v1.pdf>) · [arXiv:1008.0238](https://arxiv.org/abs/1008.0238)

No further paper is currently required for the planned braid normalization,
cyclic-sliding, parabolic-subgroup and braid-CRS milestones. Their detailed
statements still need to be read and cross-checked during implementation.

### Optional centralizer references

Open copies were also downloaded for the possible centralizer extension:

- N. Franco and J. González-Meneses, *Computation of centralizers in braid
  groups and Garside groups* (2003). [Local PDF](</home/tcsiz/Downloads/resources/Franco-GonzalezMeneses(2003)ComputationOfCentralizersInBraidGroupsAndGarsideGroups-v2.pdf>) · [arXiv:math/0201243](https://arxiv.org/abs/math/0201243)
- J. González-Meneses and B. Wiest, *On the structure of the centralizer of a
  braid* (2004). [Local PDF](</home/tcsiz/Downloads/resources/GonzalezMeneses-Wiest(2004)StructureOfCentralizerOfBraid-v2.pdf>) · [arXiv:math/0305156](https://arxiv.org/abs/math/0305156)

These possible targets still need user confirmation before being added to the
implementation scope.

### Birman–Ko–Lee/dual normal forms

- J. S. Birman, K. H. Ko and S. J. Lee, *A new approach to the word and
  conjugacy problems in the braid groups* (1998), arXiv:math/9712211v2.
  Equation (4) defines the band generators in the Artin alphabet and
  equation (5) specifies δ in descending order; §3 supplies canonical factors
  and left-weighted normal forms, including the 14 examples in `B₄`, and §4
  gives the normalization algorithm. [Local PDF](</home/tcsiz/Downloads/resources/Birman-Ko-Lee(1998)NewApproachToWordAndConjugacyProblemsInBraidGroups-v2.pdf>) · [arXiv](https://arxiv.org/abs/math/9712211)
- D. Bessis, F. Digne and J. Michel, *Springer theory in braid groups and the
  Birman–Ko–Lee monoid* (2002). §1 develops the local-monoid construction;
  Remark 3.2 fixes clockwise puncture numbering and the correspondences
  between its generators, Artin generators and BKL generators. [Local PDF](</home/tcsiz/Downloads/resources/Bessis-Digne-Michel(2002)SpringerTheoryInBraidGroupsAndBirmanKoLeeMonoid.pdf>) · [MSP](https://msp.org/pjm/2002/205-2/pjm-v205-n2-p02-p.pdf)
- D. Bessis, *The dual braid monoid*. §4.1 identifies type-A dual braids with
  the BKL monoid and describes non-crossing partitions and the Coxeter-cycle
  correspondence. [Local PDF](</home/tcsiz/Downloads/resources/Bessis(2002)DualBraidMonoid-v2.pdf>)

## Local implementation references

- Farey package: `/home/tcsiz/Projects/phd-code/Kleinian.jl/packages/Farey`.
- Preferred word layer: `src/MonoidWords.jl`; public loading boundary:
  `src/Monoids.jl`.
- Farey cleanup direction: `docs/cleanup-tracker.md` explicitly prefers
  `MonoidWords.jl` over `FreeGroup.jl`.
- Existing Chevie experiment: `working/braid_testing.jl`; useful as a source
  of possible differential tests, not as an implementation dependency.

## Source discipline

A citation in the CRS paper is not a substitute for reading a delegated
algorithm. Record exact theorem/algorithm locations as sources are obtained.
Keep source statements, implementation conventions, and independently added
optimizations distinguishable in code comments and documentation.
