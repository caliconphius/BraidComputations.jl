# BraidComputations.jl

**Working title; initial transducer and normal-form implementation.**

A Julia package for classical braid-group computations based on:

- Chapter 9 of Epstein et al., *Word Processing in Groups*: simple braids,
  the maximal-tail transducer, greedy and mixed canonical forms, and the
  associated automatic structures.
- Garside-theoretic conjugacy calculations, especially cyclic sliding and
  sliding circuits.
- Cumplido, González-Meneses and Perego (2025), *Canonical Reduction Systems
  in Artin–Tits groups of spherical type*, concentrating first on its
  braid-specific Algorithm 3.

Implemented now:

- symmetric-group products, inverses, reduced Coxeter words, inversion sets,
  weak-order meets, and `Δ`/flip permutations;
- standard/simple braids represented by permutations;
- the Chapter 9 maximal-simple-tail transducer and right-greedy factors;
- signed **Δ-prefix right-greedy normal forms**, braid products, inverses,
  powers, equality, and induced permutations;
- standard **left normal forms**, left-weighted checks, and local sliding;
- initial/final factors, preferred prefixes, Elrifai–Morton cycling and
  decycling, cyclic sliding, transport, and iteration to a repeated circuit;
- full sliding-circuit enumeration by exhaustive standard conjugacy and
  conjugacy witnesses for small braid indices;
- standard irreducible parabolic subgroups, reduction simplices, braid
  canonical reduction systems, Nielsen–Thurston type, and outer partitions;
- dependency-free braid and standard-CRS SVG drawings, in horizontal or
  vertical orientation;
- Birman–Ko–Lee/dual simple braids, absolute-order gcds, dual local sliding,
  and left normal forms for Artin words, Farey words, and signed native band
  generators;
- conversion to and from `Farey.Monoids` words, including explicit named
  generator mappings.

Not yet implemented: mixed canonical form, the optimized
Gebhardt–González-Meneses transport/pullback circuit search, and CRS support
for general spherical Artin–Tits groups. `iterate_cyclic_sliding` reaches one
circuit; `sliding_circuits` supplies full `SC(x)` by breadth-first simple
conjugacy, with conjugators satisfying `γ^c = α`. CRS and conjugacy routines
are intended for small braid indices, not for large-scale production use.

## Quick start

```julia
using BraidComputations
using Farey.Monoids

f, g, h = FreeGen("f g h")  # σ₁, σ₂, σ₃
B = BraidGroup(4)

x = B(f*g*f; generators=(f,g,h))
y = B(g*f*g; generators=(f,g,h))
@assert x == y

z = x * inv(B(h; generators=(f,g,h)))
@show z
@show normal_form(z)
@show permutation(z)
@show farey_word(z; generators=(f,g,h))
```

Run the examples with:

```sh
./julia-project examples/getting_started.jl
./julia-project examples/garside_operations.jl
./julia-project examples/bkl.jl
./julia-project examples/huxford_schillewaert.jl
./julia-project examples/visualize.jl
```

See [`docs/conventions.md`](docs/conventions.md) for the exact product,
permutation, transducer, Δ-prefix normal-form, and Farey-mapping conventions.
The stored right-greedy form is deliberately not called left normal form.
See [`examples/huxford_schillewaert.jl`](examples/huxford_schillewaert.jl)
for explicit calculations from the paper of that name, including the bounded
free-rank-two equation verification.

To generate SVG files without String or plotting dependencies:

```julia
x = BraidGroup(4)([1, 2, 1, 1, 2, 1, 3])
save_braid_svg("braid.svg", x; orientation=:left_to_right)
save_crs_svg("crs.svg", x)
save_outer_crs_svg("outer.svg", x)
```

## Development setup

From this directory, with the current local checkout of Farey:

```sh
./julia-project -e 'using Pkg; Pkg.develop(path="../Kleinian.jl/packages/Farey"); Pkg.instantiate(); Pkg.precompile()'
./julia-project -e 'using Pkg; Pkg.test()'
```

Start a project REPL with `./julia-project`, or run a script with
`./julia-project path/to/script.jl`. The launcher uses the project-local
`.depot` for caches and newly installed packages, with your usual Julia depot
as a fallback; it does not modify your global startup file.

The initial project targets Julia 1.10 or later, matching Farey's declared
compatibility.

## Sources and planning

Local and online sources are recorded in [`docs/source-map.md`](docs/source-map.md).
See [`docs/architecture.md`](docs/architecture.md) for the broader design and
[`docs/development-plan.md`](docs/development-plan.md) for subsequent stages.
