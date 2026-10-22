# Cosmochlore

**COSMO**chlo**R**e (**Co**ntinuous **S**hape/**S**ymmetry **M**easurements & **O**ctahedral distortion, in **R**ust) is a fast, pure-Rust command line tool that quantifies how closely a molecular structure (or any set of 3D points) matches an idealised reference geometry.

**Documentation: https://yluro.github.io/cosmochlore-docs**. It covers installation, usage, every command's options and examples, the reference shapes and point groups, and the theory behind each measure.

| Command | Measures                                                                                                  |
|---------|-----------------------------------------------------------------------------------------------------------|
| `cshm`  | Continuous Shape Measures (CShM) against 92 built-in reference polyhedra or user-defined YAML shapes.    |
| `csom`  | Continuous Symmetry Operation Measures (CSoM) against 47 point groups.                                    |
| `odis`  | OctaDist-style octahedral distortion parameters of six-coordinate centres.                                |
| `gidx`  | Geometry indices τ4, τ4′, τ5, τ6, τ6′, τ8 and τ8′ of four-, five-, six- and eight-coordinate centres.     |

## Quick start

Precompiled binaries are on the [releases page](https://github.com/Yluro/cosmochlore/releases). To build from source with a stable [Rust](https://rustup.rs/) toolchain:

```bash
git clone https://github.com/Yluro/cosmochlore
cd cosmochlore
cargo build --release          # binary at target/release/cosmochlore
```

```bash
cosmochlore <COMMAND> <NAME.xyz> [OPTIONS]
cosmochlore cshm tests/FeHS.xyz
cosmochlore help <COMMAND>
```

## How it works

Every command reads one `.xyz` file. The first atom is the centre unless `-c/--center` picks another one or `-n/--nc` says there is none. Results are printed as a table, and CSV/XYZ/MOL2 files are written next to the input only when a flag asks for them.

- **`cshm`**: the problem and reference shapes are centred and normalised, and the optimal rotation comes from the SVD of their correlation matrix (Kabsch). The vertex pairing is found by an exact branch-and-bound search over permutations, with the centre pinned to the reference centre and partial assignments pruned by a lower bound on their final S. The cost grows factorially, so it is impractical beyond about 13 vertices.
- **`csom`**: the structure is rotated to the orientation that minimises the mean deviation over the group's operations. Seeds are axis directions on a Fibonacci hemisphere, each also spun through the group's in-plane period when it has one. Every seed is refined with Nelder-Mead (`argmin`) over a rotation vector. Each operation's deviation uses a Hungarian atom matching within element groups, with the centre pinned. The result is the best value found, an upper bound on the true minimum.
- **`odis`**: closed-form OctaDist parameters from the six centre-ligand vectors. `--full` also runs `cshm` (`OC-6`, `TPR-6`) and `csom` (eight octahedral subgroups).
- **`gidx`**: closed-form angle indices, picked by the ligand count.

Every failure path returns an error. The program never prints a number when something went wrong.

## Code structure

```
src/
  main.rs            runs the subcommand, hands its result to out; prints any Error to stderr
  cli.rs             clap definitions of every subcommand and flag
  error.rs           crate-wide Error wrapping each module's own error type (thiserror)
  xyz.rs             .xyz parser -> Structure (optional centre + ligands)
  yaml.rs            parser for user-defined reference shapes (cshm --ref)
  geometry.rs        centring, normalising, rotations, angles, planes, volumes
  out/
    mod.rs           say! (--silent), banner, output file naming
    table.rs         Table / Row / Width: terminal tables
    csv.rs           CSV writer (stopgap: rows listed by hand)
    cshm.rs ...      one file per subcommand: its tables, CSVs, XYZ and MOL2 files
  cshm/
    mod.rs           run, calc_cshm
    permutations.rs  branch-and-bound search for the best vertex pairing
    bounds.rs        lower bounds that prune the search
    linalg.rs        correlation matrix, SVD rotation, S value
    shape_lookup.rs  built-in shapes by vertex count and index, user-shape checks
    automorphism.rs  shape automorphisms (test-only for now)
    types.rs         CShMResult (ideal polyhedron with --ideal)
  csom/
    mod.rs           run, calc_csom
    prepare.rs       centring modes, element labels
    seeding.rs       Fibonacci-hemisphere axis seeds and in-plane spins
    optimize.rs      Nelder-Mead refinement of each seed
    deviation.rs     per-operation and per-group deviations
    assignment.rs    Hungarian algorithm, matching within label groups
    types.rs         settings, results, CsomError
  odis/              mod.rs (run, OdisResult, OdisRun), calc.rs (parameters)
  gidx/              mod.rs (run, GidxResult), calc/tau{4,5,6,8}.rs
  data/
    standard_shapes.rs   the 92 built-in shapes, GENERATED from shapes/*.yaml
    shapes/              shapes_<N>vertex.yaml, source of standard_shapes.rs
    pgs.rs               operation matrices of the 47 point groups (principal axis on z)
    symops.rs            single named operations (unused for now)
    elements.rs          covalent radii (MOL2 bonds)
tests/               .xyz/.yaml/.csv fixtures used by the unit tests
*.py                 dev scripts: shape-table generator, point-group checker,
                     fixture generators, csom visualisers
PLAN.md              open work
SHIPPED.md           invariants and baselines of shipped work
CSM_PLAN.md          brief for the planned csm module
```

Unit tests live next to the code (`tests.rs` or `#[cfg(test)] mod tests`). There is no `lib.rs`, so there are no integration tests. Contributor and agent rules are in [AGENTS.md](AGENTS.md).

## License

**The code, binaries and sample tests are provided as is, with no warranty of any kind.** Cosmochlore is licensed under the GNU General Public Licence v3.0 (GPL-3.0); see [LICENSE.md](LICENSE.md). It is an independent implementation, not affiliated with the authors of the programs below.

- The branch and bound, the partial assignment method and the model shapes of `cshm` are based on [`cosymlib`](https://github.com/GrupEstructuraElectronicaSimetria/cosymlib) (`shp.f90`), Copyright (c) 2021 Pere Alemany, Efrem Bernuz, Abel Carreras and Miquel Llunell, MIT Licence.
- The Fibonacci seeding and point-group tables of `csom` are based on [`csom_main.py`](https://github.com/VRMNielsen/Continous-Symmetry-Operation-Measure-Program), Copyright (c) 2024 Villads R. M. Nielsen, Boris Le Guennic and Thomas Just Sørensen, GPL-3.0.
- Parts of `odis` are based on [`OctaDist`](https://github.com/OctaDist/OctaDist), Copyright (c) 2019-2026 Rangsiman Ketkaew et al., GPL-3.0.

Dependencies and acknowledgements are listed on the docs site's [license page](https://yluro.github.io/cosmochlore-docs/license), and the literature on its [references page](https://yluro.github.io/cosmochlore-docs/references).
