# Cosmochlore Improvement Plan

Revision 4, reviewed at `1ea345a`. Status legend: `FIXED` `PARTLY` `MOOT` `OPEN`
`NON-ISSUE` `NEW`. Update statuses in place as items land; don't append new sections.

Totals: 32/52 fixed, 3 partly, 3 moot, 3 non-issue, 11 open (1 of them new: C17).
Steps 1–3 of the order of work landed after this review (B10, D7, E2, E5, E6, E7, E8 — see those entries); the
`cargo fmt` pass in `e163603` moved every line reference below, so grep before trusting one.
All 56 tests pass. Working tree clean apart from untracked run output in `tests/` and `.idea/`
(see E7). Since revision 3 (`4123ede`): C6, D4 fixed; C14 half done; E4 moot; E5 regressed;
csom whole-process on `FeHS.xyz -p Oh D4h` went 3.86 s -> 0.22 s with identical results.

## A · Correctness

- A1 `FIXED` — test path was Windows-only. `src/odis/calc.rs:290`
- A2 `FIXED` — vacuous assertion missing `.abs()`. `src/odis/calc.rs:298`
- A3 `FIXED` — ligand-count error text mismatched count. `src/odis/mod.rs:29`
- A4 `OPEN` — `odis --full` writes the cshm and csom CSVs even without `--table`; the odis CSV
  itself is also written on `--full` alone. `src/odis/mod.rs:43,57,73`
- A5 `FIXED` — symmetry ops mislabeled ("i", "sigma_h"). `src/data/symops.rs`
- A6 `NON-ISSUE` — `.unwrap()`s in table printer can't fire (shape list provably non-empty). `src/out.rs`
- A7 `NON-ISSUE` — 2-vertex shapes unreachable (centering makes them identical). `src/xyz.rs:95`
- A8 `FIXED` — CSV fields now RFC 4180 escaped. `src/out.rs:9-17`
- A9 `FIXED` — central atom pinned in permutation search (8.7x speedup). `src/cshm/permutations.rs:39-43`

## B · Consistency

- B1 `FIXED` — `ReferenceShape` stores `Vector3` not `[f64;3]`. `src/data/standard_shapes.rs` (since D9)
- B2 `FIXED` — unified `Structure::atoms()` / `ReferenceShape::points()`. `src/xyz.rs:17-24`
- B3 `FIXED` — `samples`/`iterations` plain `usize`, `&str` not `&String`. `src/cli.rs`
- B4 `FIXED` — one `thiserror`-based `error::Error`, not 5 hand-rolled enums. `src/error.rs`
- B5 `FIXED` — `File::create().expect()` now propagates via `?`. `src/out.rs:92`, `src/csom/optimize.rs:57`
- B6 `FIXED` — C3/C5 trig literals now named f64 consts (COS_72, SIN_72, SQRT_3_DIV_2, ...). `src/data/pgs.rs:13-17`
- B7 `FIXED` — shape generator import path fixed. `generate_builtin_shapes_rs.py:130`
- B8 `OPEN` — `rotation_matrix()` transposes, `rotation_matrix_from_vector()` doesn't. `src/geometry.rs:97-127`
- B9 `PARTLY` — `csom --ignore` is a flag, `cshm` always ignores labels. `src/cli.rs:107-108`
- B10 `FIXED` — `csom::types::OptimiserSettings { seeds, iterations, tolerance }` with `Default` =
  `(20, 200, 1e-6)` is the one place the values live: `cli.rs` reads its `default_value_t`s from
  it (`CsomArgs::search_settings()` builds one from the flags), `calc_csom`/`search_best_axis`
  take it, and `odis --full` passes `OptimiserSettings::default()`. odis' eight deviations are
  unchanged to 3 decimals (rotation matrices move in the 4th) and `odis --full` runs ~1.6x faster
  (600 -> 375 ms on FeHS). Was: odis hardcoded `(20, 1000, 1e-8)` vs the CLI's `(20, 200, 1e-6)`.

## C · Performance (csom inner loop; now ~2·10^5 calls/point-group search after C6)

- C1 `FIXED` — labels double-stripped per cost eval; now uses pre-stripped `CsomStructure::labels`. `src/csom/deviation.rs`
- C2 `FIXED` — element grouping now precomputed once in `CsomStructure::groups`, not rebuilt per op call. `src/csom/prepare.rs:46`, `src/csom/assignment.rs:131`
- C3 `OPEN` — the z/y/x sort in `best_permutation_multiple_atoms` still runs per op call. It only
  exists because `group_by_label` collects a `HashMap` into a `Vec` in random order, so without
  it the float sum in `sds_dev` would differ in the last bit across runs. Cheap now: build the
  groups deterministically once (`BTreeMap`, or sort by first index) and drop the sort; the
  scalar path (`operation_deviation_scalar`) then doesn't need `final_a`/`final_b` at all — sum
  the squared distances straight off `pairs`. `src/csom/assignment.rs:224-238`, `src/csom/deviation.rs:103`
- C4 `NON-ISSUE` — denominator is only always `n` when centroid-centered; with `--center`/`--vector` (`has_centre`) the origin isn't the centroid, so it genuinely varies. Recomputing it is O(n) against the O(n³) Hungarian assignment in the same loop, so branching on `has_centre` to shortcut it isn't worth the complexity. `src/csom/deviation.rs:11-25`
- C5 `OPEN` — all 20 seeds fully refined instead of scoring first, refining best 2-3. Still ~10x
  on the seed loop, on top of what C6 already bought. `src/csom/optimize.rs:120`
- C6 `FIXED` — `--tol` (default 1e-6) wired to `NelderMead::with_sd_tolerance`; default
  `--iterations` 1000 -> 200; starting simplex edges 0.001 -> 0.1. Whole-process csom on
  `FeHS.xyz -p Oh D4h`: 3.86 s -> 0.22 s, same values (Oh 5.380, D4h 5.611). `src/csom/optimize.rs:80-88`, `src/cli.rs:115-123`
- C15 `FIXED` — `point_group_dev` (the hot-loop entry point) now sums per-operation deviations
  directly via `operation_deviation_scalar`; no `name.to_string()` or `Vec<CsomOperation>` per
  call. `point_group_operation_deviations`/`operation_deviation` (full breakdown) kept for the
  one-off `--full`/`--operated` report path. `src/csom/deviation.rs`
- C16 `FIXED` — corollary of C15: the one `point_group_operation_deviations` call in `calc_csom`
  (gated by `with_operations`) is now the only full computation, not a redundant second one. `src/csom/mod.rs:102-108`
- C17 `NEW` — the point group is re-resolved by name inside the hot loop. `point_group_dev`
  calls `get_pointgroup(pg)` — a 47-arm `&str` match — on every cost evaluation, then
  `to_matrix3` re-converts every `[[f64;3];3]` on every call (47 conversions per Oh eval).
  `refine_axis_from_seed` also builds a whole `HashMap<&str, Vec<Matrix3>>` via
  `get_pointgroup_map` just to test `.is_none()`, once per seed. Resolve once in `calc_csom`
  to a `Vec<Matrix3>` (or `&'static [SymmetryOperation]`) and hand that to
  `OrientationProblem`. `src/csom/deviation.rs:82,101`, `src/csom/optimize.rs:76`
- C8 `OPEN` — shapes/seeds/point-groups loops are serial; `rayon` would help. Matters more now
  that omitting `--pg` runs all 47 groups. `src/cshm/mod.rs`, `src/csom/optimize.rs:120`, `src/csom/mod.rs:98`
- C9 `MOOT` — automorphism dedup removed outright (superseded, not fixed)
- C10 `FIXED` — closed-form 3x3 eigenvalues (Smith 1961). `src/cshm/linalg.rs:36`
- C11 `FIXED` — norms/suffix sums precomputed. `src/cshm/bounds.rs:25-37`
- C12 `OPEN` — `builtin_shapes()` rebuilds full 90-shape HashMap on every lookup. `src/data/standard_shapes.rs`, `src/cshm/shape_lookup.rs`
- C13 `MOOT` — superseded by C9
- C14 `PARTLY` — `license = "GPL-3.0-only"` fixed (`734c0a9`); still no `[profile.release]`
  (`lto = "fat"`, `codegen-units = 1`, `panic = "abort"`). `Cargo.toml`

## D · Simplification

- D1 `OPEN` — `center_and_normalise` duplicates `center_by_centroid`+`normalise`; doc comment still ends at "Returns the". `src/geometry.rs:10-48`
- D2 `OPEN` — dead `points` array built only to construct `vectors`. `src/odis/calc.rs:28-46`
- D3 `FIXED` — `scale`/`centroid` now read (feed operated-coordinate writers). `src/csom/types.rs:24-28`
- D4 `FIXED` — `todo!()` gone: omitting `--pg` now analyses all 47 supported point groups
  (`bb0a88e`). Note this is a brute-force stand-in for measure 14, not detection. `src/csom/mod.rs:31`
- D5 `FIXED` — args parsed before banner prints. `src/main.rs:22-23`
- D6 `PARTLY` — `print_crab()` is now live behind the hidden `cshm --crab` flag; residue is the
  stale `//if args.crab {print_crab()}` comment in `src/main.rs:38` and `find_best_permutation`'s
  4th return (rotation matrix) still discarded by every production caller. `src/cshm/mod.rs:39`
- D7 `FIXED` — clippy 31 -> 11 -> 10 -> **0** warnings, `--all-targets`, no allows. The last four:
  `calc_csom` takes `OptimiserSettings` (B10); `branch` takes `&SearchTables` (the two point sets +
  the precomputed `hi`/norm/suffix tables) and `&mut SearchState` (partial permutation + best so
  far) instead of 11 arguments; `find_automorphisms` reuses `SearchTables::new(reference,
  reference)`; `symops::get_operation` returns `Option<&[SymmetryOperation]>`. cshm output is
  byte-identical before/after on FeHS and La03 and the 11-vertex B&B timing is unchanged
  (404 ms min both). Earlier rounds: `writeln!(file)` ×3, `.clone()`-on-`Copy` ×4, dead code (E6).
- D8 `FIXED` — `csom` module split: `types.rs`, `prepare.rs` (was `io.rs`), `assignment.rs`
  (Hungarian + label grouping, out of `dev.rs`), `deviation.rs` (was `dev.rs`), `optimize.rs`,
  `mod.rs` left with `csom_main` + `calc_csom`. Rename only, no logic change. `src/csom/`
- D9 `FIXED` — `shapes.rs` -> `cshm/shape_lookup.rs` (built-in lookup by vertex count/index,
  `check_vertex_count`, `ShapeLookupError`): reference shapes are a cshm-only concept and `csom`
  already keeps its own modules. `ReferenceShape` itself now lives next to its table in
  `data/standard_shapes.rs`, emitted by the generator, which removes the `shapes` <-> `data`
  import cycle. `yaml.rs` stays top-level for future non-shape inputs (user point groups).
  The test-only `structure_from_shape` helper moved to `cshm/test_utils.rs` with the other test
  scaffolding. Pure move: all cshm/csom/odis outputs byte-identical. `src/cshm/shape_lookup.rs`, `generate_builtin_shapes_rs.py`

## E · Build, tests, CI

- E1 `OPEN` — no `lib.rs`; no integration tests possible, all tests are `#[cfg(test)]` inline
- E2 `FIXED` — `.github/workflows/ci.yml` (`f9cccad`, `ac5035d`): Test, Clippy and Rustfmt jobs on
  push to master and on PRs, `RUSTFLAGS=-D warnings`, `--locked`; clippy runs with plain
  `-D warnings` since step 3 (the two temporary `-A` flags are gone). One-time `cargo fmt` in
  `e163603` (blame-ignored); generated tables
  in `data/pgs.rs` and `data/standard_shapes.rs` carry `#[rustfmt::skip]`, and the shapes
  generator emits it. `rustfmt.toml` pins `newline_style = "Native"` for the CRLF checkout.
- E3 `OPEN` — `bnb_is_faster_than_bf` asserts on wall-clock time. `src/cshm/tests.rs:310`
- E4 `MOOT` — `tests/FeCl6_table.csv` (and `FeCl6_ideal.xyz`) deleted in `7354442`. The SHAPE
  2.1 values are pinned inline in `src/cshm/tests.rs` (water 0.035/33.342, Eu7 2.109/2.630/7.288),
  so nothing was lost.
- E5 `FIXED` — every README row checked against `cli.rs`: `--iterations` 200, `--pg` optional (both
  places), `-e`/`--tol` row added, csom algorithm step 2 completed (and step 1 no longer claims seeds
  are scored before refinement — they are all refined, see C5). Also found and fixed: the odis
  `--full` sentence was truncated too; the odis example lacked the `Theta`/`Volume` rows and the
  `D3h` group; the cshm example used older coordinates than `tests/FeHS.xyz` while the other two
  examples used the fixture. All three example blocks are now verbatim output of the release binary
  on `tests/FeHS.xyz` + `tests/ebcT-6.yaml` (see E8 for why `ebcT-6` changed from 14.335 to 14.277).
  Was: `--iterations` said 1000, `--pg` "Currently required", no `--tol` row, step 2 cut off.
- E6 `FIXED` in `f9cccad` — `automorphism.rs` and `structure_from_shape` are `#[cfg(test)]`,
  `twist_top_face` deleted; `writeln!(file)` ×3 and four `.clone()`-on-Copy went with it.
  Was: 3 dead-code warnings in the binary: `find_automorphisms`, `automorphism_branch`
  (only caller: `cshm/tests.rs:108`), `structure_from_shape` (only callers: `odis/calc.rs` tests).
  `twist_top_face` is inside `#[cfg(test)]` but called from nowhere — dead even for tests.
  Blocks `-D warnings` in CI. `src/cshm/automorphism.rs`, `src/shapes.rs:36`, `src/odis/calc.rs:321`
- E7 `FIXED` — `.gitignore` now covers `.idea/` and every output name `out.rs` writes
  (`*_cshm_table.csv`, `*_csom_table.csv`, `*_odis_table.csv`, `*_details.csv`, `*_operated.xyz`,
  `*_merged.mol2`, `*_ideal.xyz`). `tests/FeHS_ideal.xyz`, a tracked run output, was removed with step 3 (was: added
  in `623eb10`, nothing reads it, differs from a fresh run) — `git rm` it or keep it deliberately.
  Was: 19 untracked `*_operated.xyz`/`*_merged.mol2`/`*_ideal.xyz` files in `tests/` plus `.idea/`.
- E8 `FIXED` — `tests/ebcT-6.yaml` had its `centre` and the fourth tetrahedron vertex swapped:
  `(0.5, 0.5, 0.5)` (the centroid of the `(1,1,0)/(1,0,1)/(0,1,1)/(0,0,0)` tetrahedron) was listed
  as a vertex and `(0,0,0)` as the centre. Before A9 the unpinned search silently found the swap,
  which is where the README's 14.335 came from; with the centre pinned the mis-specified shape
  scored 35.8. Swapped back in the fixture and in the README's yaml example; the README's inline
  structure now reproduces 14.335 exactly and `tests/FeHS.xyz` gives 14.277. No test used the file.

## Part II — Measures worth adding

Tier 1 (cheap, high value, reuses existing machinery):
1. `OPEN` Minimal-distortion paths / shape maps (OC-6↔TPR-6 interpolation) — Alvarez 2005.
   Fitting a structure *to* a path (the path coordinate plus the deviation from the path) is a
   one-variable dynamic shape: do it as 19c.
2. `SHIPPED` Θ face-twist + octahedral volume — Ketkaew 2021
3. `SHIPPED` τ4/τ4'/τ5/τ6/τ6' geometry indices for CN=4,5,6 — Addison 1984, Yang 2007, Okuniewski 2015,
   Stoeckli-Evans 2025. The `gidx` command (`src/gidx/`, `mod.rs` + `calc.rs`, laid out like `odis`): τ4 and τ4'
   for four ligands, τ5 for five, τ6 for six, plus the angles they are built from (two largest for CN 4/5; for
   CN 6 the three trans angles, chosen as the disjoint pairing of the ligands with the greatest angle sum, since the
   paper leaves the choice open). Pinned on the ideal built-in shapes and on `tests/ML4.xyz`, `tests/ML5.xyz` and
   `tests/FeHS.xyz`; the README lists the values for every 4-, 5- and 6-vertex ideal shape. Uses the exact tetrahedral
   angle arccos(−1/3) rather than the papers' rounded 109.5°/141°, so an ideal tetrahedron gives exactly 1. The
   paper's 1.00 for a pentagonal pyramid sums the five adjacent base angles (5 × 72°) rather than three trans angles;
   gidx uses three trans angles for every geometry (ideal shape: 0.900). The paper's own Cd centres (τ6 0.27, τ5 0.31)
   are reproduced from its SI: `tests/CdONCl4.xyz`, `tests/CdONCl3.xyz`.
   τ6' is cosmochlore's own complement to τ6 (not in the paper): (5×90° − sum of the five smallest angles)/90°, so 0 for
   OC, 0.456 for the equal-edge prism, 1 for a pentagonal pyramid, 5/3 for a hexagon. Five angles rather than three
   because the five base angles of a real pyramid sum to ~360° (three smallest read 1.31/1.18 on VAPSAG/CIJBII,
   five read 1.044/0.995); (1+τ6')/2 is the paper's pyramid value. Literature sites (CSD, one Cd each) in
   `tests/cd6/`: IQATAY Cd1-4, ULESAH, VAPSAG, CIJBII.
4. `OPEN` Classic distortion params (⟨λ⟩, σ², Baur D, ECoN) — Robinson 1971 et al.
5. `OPEN` Bond-valence sum — Brown & Altermatt 1985
6. `OPEN` Gyration-tensor descriptors (asphericity, κ²) — reuses `linalg.rs`
7. `OPEN` Planarity/pyramidalisation via covariance SVD
20. `OPEN` Polynator by-products — Link & Niewa 2023 (see 19), manual §3.5, §3.6, §5. Each is O(n) on
    results cosmochlore already has:
    - δ = 10·√S (`--delta`) as an extra column for `cshm` and for 19. δ grows linearly with small
      distortions while S grows quadratically, and 2 decimals of δ always separate "exact" from
      "slightly distorted". Display only: never change what is minimised or stored.
    - Shape-free values from the eigenvectors of Σ q qᵀ (centroid-centred): δ_linear (distance to the
      best line), δ_planar (to the best plane) and δ_spherical (spread of |q_i| about its mean). The
      same eigen-decomposition as items 6 and 7, so build all three together.
    - Per-vertex deviations of the fitted ideal shape (Polynator's `.aso` file): |q_i − v_i| in Å, its
      mean and SD, split into radial (|q_i| − |v_i|) and angular parts. `CShMResult.xyz` already
      holds v_i. Item 17 step 1 tests the same vector against the ADPs.
    - Convex-hull volume and surface area of the structure and of the fitted ideal shape for any CN
      (`odis` only has the octahedral volume). Quickhull is enough for n ≤ 60.
    - Ranking with a parameter penalty (Polynator's "variable tax": δ + 1.5 per free variable, used
      only to order the list, never reported). Only matters once 19 mixes models with different dof.

Tier 2 (moderate effort, distinctive):
8. `OPEN` Bailar/Ray-Dutt twist angles — Avdeef & Fackler 1975
9. `OPEN` Zabrodsky-Avnir CSM (folding/unfolding algorithm) — JACS 1992/93
10. `OPEN` Continuous chirality measure — Zabrodsky & Avnir 1995
11. `OPEN` Structure-to-structure CShM (`--ref` accepts `.xyz`)
12. `OPEN` Element-weighted CShM (closes B9)
19. `OPEN` Dynamic (deformable) shape measures — Link & Niewa, *J. Appl. Cryst.* 2023, 56, 1855
    (doi 10.1107/S1600576723008476); Polynator 1.7 manual (`poly_doc_1_7.pdf`) and the 1.7.1 Python
    source. The reference may deform within stated constraints instead of being rigid: OC-6 becomes
    "any D4h tetragonal bipyramid", "any D3 twisted prism" or "a point on the Bailar path".
    - **What Polynator does**: a model is a stack of *belts*, rings of n vertices perpendicular to a
      model axis. Each belt has a height h, a radius w, an azimuth φ, and the same three modulated
      by cos/sin(2πf(p+o)/n) over the vertex index p (`~h`, `~w`, `~φ`, which turn a square into a
      disphenoid, rhombus or rectangle). There is also a scale `sc`, and "bundles": parameters given
      as string formulas of 1–3 variables (Bailar twist, Berry pseudorotation, pyritohedron). It has
      241 models in 1.7.1 (211 in the paper), 139 of them with ≥ 2 variables. Its metric is
      δ = 100·√(Σ|a−v|²/Σ|a−c|²) = 10·√S. The vertex assignment is heuristic: 92 scan directions,
      a belt-cost estimate, dihedral ordering inside each belt, pairwise-swap repair, or an
      assignment inherited from a fitted parent model. The paper only vouches for it below δ ≈ 30
      (S ≈ 9). The fit is block coordinate descent: Kabsch, then φ, sc, w, h and the bundles, each by
      1-D step-halving, for at most 10 cycles. See `ModelVertexAssigner` and `ModelFit`,
      `polynator_main.py:1703-2232`. Rigid models match `cshm` exactly when the centre sits at the
      centroid. On `tests/FeHS.xyz`: OC-6 2.08926, TPR-6 11.07946, HP-6 33.20982.
    - **Key observation**: every Polynator model named after a point group ("[4/mmm]", "[-3m]",
      "[32]") is the set of *all* configurations with that symmetry and that orbit structure. Build it
      from `data/pgs.rs` instead of belts:
      - Split the parent's points (centre included) into orbits of H. Add the identity, which
        `pgs.rs` leaves out of its tables.
      - For an orbit with representative r, the free coordinates are the fixed subspace of r's
        stabiliser: the range of (1/|Stab|)·Σ_{g∈Stab} g, of rank 0–3.
      - Each orbit point is g·basis·θ_orbit.

      The model is **linear**, v(θ) = Bθ, with scale included. Twists are linear too, because
      (w cos φ, w sin φ) is just (x, y). So for a fixed permutation and rotation, θ is one
      least-squares solve, and for a fixed θ the rotation is Kabsch. Polynator's h, w and φ come back
      out of θ for the report. After re-centring, B loses one rank for polar groups (C_n, C_nv): all
      z-parameters can shift together along the axis. Solve with SVD or a pseudo-inverse.
    - **Reference values**: on `tests/FeHS.xyz` (centre included, exact permutations), every
      symmetry family reproduces Polynator 1.7.1 to 1e-5 in S. These are the regression fixtures,
      stored in two files:
      - `tests/dshm_reference.csv`, written by `generate_dshm_fixtures.py`. That script is a numpy
        reference implementation that reads `pgs.rs` and the YAML shapes. Per family it gives the
        frame, dof, anchors, the exhaustive S, and the planned search's S from the parent seed
        alone and from all seeds.
      - `tests/polynator_reference.csv`, written by `generate_polynator_fixtures.py`, which runs a
        local Polynator copy headless (`--polynator DIR`; not redistributed). Its header records
        the Polynator version and the sha256 of `polynator_main.py`.

      The Rust tests read these CSVs instead of hard-coding numbers:

      | Family: parent, H (frame) | dof | S, reference | Polynator model | S, Polynator |
      |---|---|---|---|---|
      | OC-6, Oh | 1 | 2.08926 | octahedron[platonic] | 2.08926 |
      | OC-6, D4h | 2 | 2.08550 | tetragonal_bipyramid[4/mmm] | 2.08550 |
      | OC-6, D3d (C3 ‖ [111]) | 2 | 2.06334 | trigonal_antiprism[-3m] | 2.06334 |
      | OC-6, D2h (C2′ through vertices) | 3 | 2.08434 | rhombic_bipyramid[mmm] | 2.08434 |
      | OC-6, D2h (C2′ between vertices) | 3 | 2.05755 | rectangular_bipyramid[mmm] | 2.05755 |
      | OC-6, D2d (C2′ through vertices) | 2 | 2.08550 | — (collapses onto D4h) | — |
      | OC-6, D2d (C2′ between vertices) | 3 | 0.36563 | didigonal_scalenohedron[-42m] | 0.36563 |
      | OC-6 or TPR-6, D3 | 3 | 0.80736 | twisted_trigonal_prism[32] | 0.80736 |
      | TPR-6, D3h | 2 | 10.46473 | trigonal_prism[-6m2] | 10.46473 |
      | OC-6, C4v, free centre | 4 | 1.99555 | — | — |
      | OC-6, C4v, `--fix-center` | 3 | 2.05807 | tetragonal_heterobipyramid[4mm] | 2.05807 |
      | OC-6, C3v, free centre | 4 | 1.99874 | — | — |
      | OC-6, C3v, `--fix-center` | 3 | 2.05403 | trigonal_antifrustum[3m] | 2.05403 |

      The two D2d rows show that the *embedding* of H in the parent's group defines the family, not
      the label H. With the C2′ axes through the equatorial vertices, D2d collapses onto D4h. With
      them between the vertices, it adds puckering. FeHS is a D2d-distorted octahedron: 0.366,
      against 2.089 for OC-6 and 2.086 for D4h. Polynator keeps such pairs as separate models (its
      cuboctahedron tree has 42m (1)/(2) and mmm (1)/(2)).
    - **Centre atom** (decided): the centre is an orbit of size 1. It has 0 dof when H fixes a single
      point, and 1 dof (along the axis) for C_n and C_nv.
      - **Default, free**: this gives the CSM-consistent "nearest H-symmetric structure, metal
        included".
      - **`--fix-center`**: pins the centre's model point to the centroid of the ligand model points,
        B_centre = (1/N_lig)·Σ B_ligand. This is still linear, removes the centre's own columns, and
        is Polynator's convention: it reproduces the polar rows above exactly.

      Both modes give the same result whenever H fixes a single point (all rows except C_n/C_nv).
      Report the mode in the output header and the CSV. `-n` (no centre) needs no special case.
    - **Permutation search**. This is Polynator's weak point, and cosmochlore's B&B is exact for any
      rigid reference:
      1. *Seeds*: the parent's θ0, plus the θ of every built-in rigid shape of the same vertex count
         that lies inside the family (its *anchors*). The parent alone is **not enough**. TPR-6 under
         C3v (σv through the vertices), seeded from TPR-6, converged to 10.462 even with all 12
         automorphisms. The family also contains OC-6 (bottom triangle with w → −w), and seeding
         from OC-6's θ reaches the brute-force minimum, 1.99874 (the `s_parent_seed` and
         `s_all_seeds` columns of `tests/dshm_reference.csv`). Compute anchor θs offline: fit each
         built-in shape into each family over all permutations (CN ≤ 8; for larger CN, the runtime
         pipeline from every seed found so far) and keep those with S < 1e-8. Store them in the
         generated data, and pin them with a test. `generate_dshm_fixtures.py` already does this
         for CN 6 (its `anchors` column; e.g. OC-6/D3 has HP-6, OC-6 and TPR-6).
      2. *Per seed*: B&B (`find_best_permutation`) of the problem against v(θ_seed) gives P0 in the
         family's vertex order. Refit P0 composed with each automorphism of v(θ_seed), and drop
         duplicates modulo H: about |Aut|/|H| fits (3 for OC-6 → D4h, 4 for → D3d, 6 for → D2d).
         The automorphisms come from `cshm::automorphism::find_automorphisms`, which is currently
         `#[cfg(test)]` and must be un-gated.
      3. *Fixed-permutation fit*: alternate Kabsch (improper allowed, as in `cshm`) and the θ solve,
         with B re-centred so the translation stays optimal. This takes 4–5 iterations to 1e-14 on
         FeHS.
      4. *Re-assign*: run the B&B with the fitted v(θ) as a rigid reference. If P changes and S
         drops, refit and repeat. Both steps only ever lower S, so this terminates.
      5. *Guarantees*, each pinned by a test:
         - S_dyn ≤ S_rigid(anchor) for every anchor.
         - S_dyn(child) ≤ S_dyn(parent) for nested families. Warm-start the child with
           θ = B_child⁺·v_parent.
         - The pipeline equals brute force over all permutations on the CN ≤ 7 fixtures.

      From CN 12 up, every B&B call in steps 2 and 4 follows item 21's `auto` rule: belt matching
      plus the capped B&B certificate, instead of the plain `find_best_permutation`.
    - **Path families** are nonlinear, with one variable t:
      - Bailar, OC-6 ↔ TPR-6 along D3. This is Polynator's bailar_twist[dynamic]: vertices on the
        unit sphere, twist t, h² = (¼ + ½cos t)/(1¼ + ½cos t).
      - Berry, TBPY-5 ↔ SPY-5.
      - The minimal-distortion paths of item 1.

      The fit is a golden-section/Brent search over t, with the exact rigid CShM at each step (29
      steps). On FeHS (the Bailar row of `tests/dshm_reference.csv`): t = 0.82862 rad, identical to
      Polynator, and the same h, w and twist.
      S is 1.20321 against Polynator's 1.19024. The cause: for these "pseudopolyhedron" models,
      Polynator shifts the ligands by (centre − centroid)/N but leaves the centre where it is
      (`vec_to_correct_centering`, `polynator_main.py:1404`). That is not a rigid translation, so
      keep cosmochlore's convention and document the gap. Two-variable paths use Nelder–Mead
      (`argmin`).
    - **Data**: `src/data/shapes/dynamic_<n>vertex.yaml`, one entry per family. For example
      `OC-6/D2d: {parent: OC-6, group: D2d, frame: {z: [0,0,1], x: [1,0,0]}, name: Didigonal
      scalenohedron}`. Here `frame` gives the directions, in the parent's coordinates, of the
      `pgs.rs` table's z and x axes. `pgs.rs` puts the D2d C2′ axes on the diagonals, so this is the
      puckering row; `x: [1,1,0]` gives the collapsed one. For OC-6/D3d the frame is
      `z: [1,1,1], x: [1,1,-2]`. Orbits and B are built at load time (|H|·N work). A frame whose
      operations do not map the parent's points onto themselves to 1e-9 is an `Error`, never a
      silent fit. User families via `--ref` come later.
    - **First library**:
      - CN 4: T-4 → D2d, C3v; SP-4 → D2h ×2, D2d.
      - CN 5: TBPY-5 → D3h, C3v; SPY-5 → C4v.
      - CN 6: OC-6 → D4h, D3d, D2h ×2, D2d, C4v, C3v, D3, C2v; TPR-6 → D3h, C3v.
      - CN 7: PBPY-7 → D5h; COC-7 → C3v; CTPR-7 → C2v.
      - CN 8: CU-8 → D4h, D2d; SAPR-8 → D4d; TDD-8 → D2d; BTPR-8 → C2v. In Polynator's Bi[VO4] case,
        the dynamic TDD has under a quarter of the δ of the dynamic BTPR or SAPR.
      - CN 12: COC-12 → Td (elpasolite), Th, D4h, D3d, D2d ×2, D3. D3 and one D2d are the best fits
        across Polynator's perovskite survey (CaTiO3, Table 3).
      - Paths: Bailar and Berry.
    - **Output**: family rows indented under their parent (a symmetry tree), with columns Family,
      H, dof and S, plus δ behind `--delta` (item 20). Fitted parameters in Å and degrees, derived
      from θ:
      - per orbit, the height along the family axis, the radius and the azimuth;
      - h/w ratios;
      - the twist between stacked polygons (its sign is only meaningful with proper rotations, see
        item 18's chirality note);
      - the centre displacement.

      Fitted coordinates go to an `.xyz` writer, as with `cshm --ideal`.
    - **CLI** (decided): a new **`dshm`** subcommand ("Dynamic Shape Measures"), so its settings stay
      out of `cshm`. Reuse `cshm`'s flag names where the meaning is the same. Draft:
      - `<NAME>`, `-n/--nc`, `-c/--center POS`, as in `cshm`.
      - `-s/--sh <IDX>...`: only the families of these built-in parents (same indices as `cshm`).
        `--family <SYMBOL>...` selects single families, e.g. `OC-6/D2d`. The default is every family
        for the vertex count, paths included.
      - `--fix-center`: centre pinned to the ligand centroid (see above). The default is free.
      - `--delta`: add a δ = 10·√S column.
      - `-t/--table`: `<name>_dshm_table.csv`, with the family, H, dof, centre mode, S and every
        fitted parameter.
      - `-i/--ideal`: one `.xyz` of fitted coordinates per family. Add the new output name patterns
        to `.gitignore` (E7).
      - `-r/--ref <YAML>...`: user families (later).
      - Search knobs, with "it is recommended not to change" help text as `csom` has:
        `-T/--tolerance` (fit stop, ΔS < 1e-12) and `--iterations` (fit cap). `--exhaustive` runs
        brute force over all permutations for verification (an error above CN 8).

      Keep one `DshmSettings` with a `Default` impl as the single source of these defaults, and have
      `cli.rs` read its `default_value_t`s from it (as B10 did for `OptimiserSettings`). Errors go
      through a `DshmError` in `error::Error`; never print a number from a fit that failed.
    - **Phases**:
      - 19a, linear-family core in `src/dshm/`: orbit builder, fixed-P fit with both centre modes,
        and the seed/coset/re-assign driver, with the OC-6 and TPR-6 families and the fixtures
        above. Shared `cshm` pieces (`find_best_permutation`, `find_automorphisms`, `linalg`) are
        reused, not copied. `cshm`, `csom` and `odis` output must stay byte-identical.
      - 19b, the YAML library for CN 4–8, plus the `dshm` subcommand and its CSV/xyz writers.
      - 19c, path families, together with item 1.
      - 19d, auto-derived trees: enumerate the non-conjugate subgroup embeddings of Aut(parent)
        instead of hand-written frames. Polynator's hand-made trees (paper Fig. 1 and Fig. 10) are
        exactly this enumeration.
      - 19e, constrained families (equal edges, planar faces) via a nonlinear θ and Nelder–Mead,
        only if someone asks for them.
    - **Shared machinery**:
      - The fixed-P θ solve is the folding/unfolding projection of `CSM_PLAN.md` §2.2, so 19a is the
        averaging kernel that Phase 7 there needs for multi-generator groups.
      - The orbit bases are the projectors item 13 needs.
      - Item 17's per-vertex thermal test applies to v(θ) unchanged.
    - **Not adopted from Polynator**:
      - CIF input with space-group expansion, graph tracing and Voronoi ligand selection: Timeo and
        Olex2 own the crystal side.
      - Its model tables: the 1.7.1 source carries no licence statement, so re-derive every family
        from `pgs.rs`, which 19d automates. Do not copy them.
21. `OPEN` Belt matching for large CN in `cshm` — Polynator's assignment heuristic (Link & Niewa 2023
    §2.1 steps 4–5; `ModelVertexAssigner`, `get_ordered_belts` and `find_beneficial_permutations`
    in `polynator_main.py`), made safe by a capped B&B certificate. **No Hungarian or ICP
    re-assignment in `cshm`** (decided). The Hungarian matching stays where it is: in `csom` and
    in `CSM_PLAN.md` Phase 6, the approximate CSM mode (also decided). Do not change either.
    - **Why**: the exact B&B stops being usable above CN ≈ 20. Measurements on the release build,
      centre included. Each structure is a built-in shape scaled to 2.2 Å, with Gaussian noise of σ =
      2 % or 8 % of the radius, and its atoms shuffled:

      | Structure (noise) vs reference | S | B&B | B&B nodes | B&B started from the belt result |
      |---|---|---|---|---|
      | CN ≤ 12, own shape | 0.06–1.6 | 10–20 ms | | |
      | IC-12 or COC-12, all 13 CN-12 shapes | | 2.6 s | | |
      | DD-20 (8 %) vs DD-20 | 1.419 | 0.62 s | 626 150 | 0.01 s, 2 726 nodes |
      | TCU-24 (2 %) vs TCU-24 | 0.093 | 15.9 s | 12 545 698 | 0.01 s, 1 306 nodes |
      | TOC-24 (2 %) vs TOC-24 | 0.098 | 30.0 s | 24 613 755 | 0.01 s, 1 180 nodes |
      | TCU-24 (8 %) vs TCU-24 | 1.513 | 1.74 s | 1 915 356 | 0.02 s, 4 882 nodes |
      | TCU-24 (8 %) vs TOC-24 | 8.172 | 164 s | 117 493 907 | 131 s, 73 391 644 nodes |
      | TCOC-48, TIC-60 (2 % and 8 %), own shape | | > 60 s, unfinished | | |

      Less noise is *slower* at CN 24, because many near-equivalent branches survive the bound.
      Shapes the structure does not resemble cost the most: the bound is weak at high S.
      The structures are committed as `tests/belt/<SHAPE>_<NN>.xyz` (NN = noise in %), written by
      `generate_belt_fixtures.py` with a fixed seed. The node counts came from a temporary
      instrumented build, so 21b's node counter reproduces them, not the script.
    - **Algorithm** (steps 1–5 are implemented in numpy in `generate_belt_fixtures.py`, as the
      reference for the Rust port):
      1. *Reference belts*, once per shape and cacheable. Candidate axes: vertex directions, sums
         of vertex pairs (C2 axes through edge midpoints) and normals of vertex triples (centres
         of planar rings; the C5 of PPR-10 is only found this way). Take the axis with the fewest
         height layers (tolerance 1e-3 of the radius) whose centroids all lie on it. Sort each
         layer by azimuth from one common zero: the first vertex of the largest layer. This gives
         OC-6 3/3, PPR-10 5/5, IC-12 1/5/5/1, COC-12 4/4/4, TCOC-48 8×6 and TIC-60
         5/5/10/10/10/10/5/5. The centre is not in any layer and stays pinned.
      2. *Layer assignment*: for each scan direction u (200 Fibonacci directions plus ± the three
         principal axes), sort the atoms by u·q̂ and fill the layers from the top down. Drop
         assignments already seen.
      3. *Cheap ranking* (Polynator's q1). The layer-plane normal is the smallest eigenvector of
         Σ(q − c_layer)(q − c_layer)ᵀ. The cost is the out-of-plane spread, plus the radius spread
         within each layer, plus n_layer·|c_layer × normal|². Keep the 8 best assignments.
      4. *In-layer pairing*. Take the axis from the layer centroids (the largest eigenvector, or the
         line through two centroids) and measure azimuths about it. Put each atom of the largest
         layer on the reference zero in turn, in both handednesses (`cshm` allows improper
         rotations). In every layer, sort by relative azimuth and take the cyclic shift that
         minimises Σ(1 − cos Δφ). Score every resulting permutation *exactly* (one 3×3 SVD), so
         Polynator's angular estimate q2 is not needed.
      5. *Repair*: Polynator's pairwise swaps (swap two atoms whenever it lowers S, repeat until
         none does) on the 3 best permutations.
      6. *Certificate*: start the existing B&B with the belt permutation as its incumbent (`best_s`,
         `best_perm`, `best_rot_matrix`), under a node budget. If the search completes, the value is
         exact. If the budget runs out, report the belt value as **not proven**. It is always a
         valid upper bound, since it is a real permutation's S.
    - **Reference results** (`tests/belt_reference.csv`: per structure and built-in shape of its
      CN, the layers, `s_belt`, cosmochlore's exact `s_exact` with an `exact_status` of ok, timeout
      or skipped, and the generating permutation's `s_oracle`):
      - **Structure vs its own shape**: the belt matching found the optimum every time, from CN 6 to
        60 (FeHS, PPR-10, JCPPR-11, IC-12, COC-12, DD-20, TCU-24, TOC-24). For TCOC-48 and TIC-60
        it equals the generating permutation, which bounds the optimum. Python time was 0.1–0.5 s
        at CN 20–24 and 2–4 s at CN 48–60, mostly the per-shape axis search over ~N³/6 triples,
        which is cacheable.
      - **Certificate**: started from the belt result, the B&B proves the optimum in 0.01–0.02 s,
        with 230–21 000× fewer nodes (table above). For TCU-24 vs TOC-24 the belt value was
        already optimal, but the proof still took 131 s. That is why the budget exists.
      - **Misses**: 34 of 104 pairs at CN ≤ 12 missed, never a structure against its own shape.
        They are of two kinds:
        - low-symmetry references whose layers are single vertices, which leaves no rings to
          exploit (JMBIC-10, JATDI-10, JSPC-10, JASPC-11);
        - shapes the structure does not resemble (S ≳ 6; misses from +0.07 to +6).

        This matches the paper's reliability limit (δ < 30, S < 9). Sorting the reference layers
        by unit-vector height, as the problem side does, left the count at 34. Real structures: all
        five FeHS shapes were found exactly. The La03 misses (JAPPR-11 +0.53, JASPC-11 +2.0, HP-11
        +0.75) are of these two kinds.
    - **CLI** (decided): `cshm -m/--match <auto|exact|belt>`, default `auto`. Model it on `csom`'s
      `-m/--mode`: a `#[derive(Debug, Clone, clap::ValueEnum)] pub enum MatchingMode { Auto,
      Exact, Belt }` (next to `CShMResult` in `src/cshm/types.rs`), declared in `CshmArgs` as
      `#[arg(short = 'm', long = "match", value_enum, default_value = "auto")]`, as `CenteringMode`
      is in `src/csom/prepare.rs` and `CsomArgs`. `-m` is still free in `cshm`.
      - `auto`, the default: `exact` when CN < 12, `belt` when CN ≥ 12. CN is the number of
        ligand vertices, `structure.ligands.len()` in `cshm_main`, with the centre not counted.
      - `exact`: today's B&B at any CN, unchanged. Its help text should say that it can run for
        minutes or hours above CN ≈ 20.
      - `belt`: steps 1–6 at any CN. Below 12 it is mainly for tests and comparisons.

      At CN = 12, `auto` switches to `belt`. In every CN-12 measurement above, the plain B&B
      finished in well under the budget, so the certificate completes and the values equal
      today's. The permutation can still be a different, symmetry-equivalent one, which reorders
      the atoms in the `--ideal` file. Byte-identical output is therefore guaranteed for
      `--match exact` always, and for `auto` only when CN < 12.

      `--nodes N` sets the certificate budget. Default 10⁷: the B&B ran at ~0.7 M nodes/s, so that
      is about 14 s per shape. Agents use this value and flag it in the PR; the maintainer may
      change it.

      Keep the mode, the budget and the belt knobs (200 scan directions, keep 8, repair 3) in one
      `CshmSettings` with a `Default` impl, with `cli.rs` reading its `default_value_t`s from it
      (B10 pattern). Output:
      - a header line naming the mode used (e.g. `Matching: auto (belt, CN 24)`);
      - unproven values printed with a marker (`8.172*`) and a footnote line under the table;
      - an `exact` column in the CSV.

      Never print an unproven value without its marker (README promise).
    - **Code**:
      - `src/cshm/belts.rs`: `ReferenceBelts { axis, layers, azimuths }` built once per
        `ReferenceShape`, and `belt_match(problem, reference, has_centre, &settings) -> (s, perm,
        rot)`.
      - `permutations.rs`: `find_best_permutation` takes an optional incumbent and node budget
        and returns whether it finished. The search logic is otherwise untouched.
      - `CShMResult` gains `exact: bool`.
      - Scan directions reuse `csom::seeding::fibonacci_hemisphere`, after it moves to a neutral
        module (`CSM_PLAN.md` P0.3).
    - **Tests** (no wall-clock asserts, E3):
      - Reference belts give the layer sizes in the `layers` column of `tests/belt_reference.csv`.
      - On every `tests/belt/*.xyz` structure against its own shape, belt S equals the certified
        exact S (with `exact` true), and matches `s_exact` (3 decimals) or, above CN 24, does not
        exceed `s_oracle`. No `rand` dependency is needed: the structures are committed.
      - Every row of `tests/belt_reference.csv` with `exact_status = ok` satisfies belt S ≥
        `s_exact` − 5e-4, i.e. the belt value is always an upper bound.
      - The Rust port misses no more often than the reference. Count rows where S > `s_exact` +
        5e-4, with the certificate off, and compare with the same count for `s_belt`. It need not
        match `s_belt` row by row, since scan directions and tie-breaking may differ.
      - A completed certificate equals plain B&B.
      - Node counts are asserted for the certified CN 20 and 24 cases.
      - CLI: `-m` defaults to `auto`; `auto` and `--match exact` give byte-identical output on the
        CN < 12 fixtures (`tests/FeHS.xyz`, `tests/La03.xyz`, plus `-n` runs); a budget of 1 node
        forces the "not proven" path, and the marker and the CSV `exact = false` appear.
    - **Phases** (one PR each; CI gates as in E2; `csom` and `odis` output byte-identical
      throughout):
      - 21a: `src/cshm/belts.rs` — `ReferenceBelts` and `belt_match` (steps 1–5) as a library
        function with its tests. No CLI change.
      - 21b: certificate (step 6). `find_best_permutation` takes an optional incumbent and node
        budget, and `CShMResult` gains `exact`. With no incumbent and no budget, results are
        unchanged.
      - 21c: `-m/--match`, `--nodes`, `CshmSettings`, the output marker, the CSV column, a README
        section ("The `--match` modes"), and new output patterns in `.gitignore` if any (E7).
    - **Not adopted from Polynator**: hand-written belts per model (derived automatically here);
      the q2 angular estimate (exact scoring is cheap); strategies 2 and 3 for prolate and oblate
      models, which built-in shapes above CN 12 do not need.

Tier 3 (larger, highest scientific return):
13. `OPEN` Symmetry-adapted distortion decomposition (irrep projection using `data/pgs.rs`). The
    orbit projectors of 19a are the same machinery, and the drop in S from parent to child family
    measures the distortion that the child's symmetry allows.
14. `OPEN` Automatic point-group detection — Nielsen 2024. D4's `todo!()` is gone, but the
    all-47-groups default is 47x the work of one group; detection is what makes it cheap.
15. `OPEN` Steric descriptors (%Vbur, cone angle) — Cavallo/Tolman
17. `OPEN` ADP-aware measures (thermal significance) — Munguba … Simas, *ACS Omega* 2025, 10,
    47189 (doi 10.1021/acsomega.5c05878). Timeo can pull each atom's ADPs from Olex2.
    - **What the paper does**: it is a *library-level* criterion, not per-structure propagation.
      Uiso = tr(U_cart)/3 (eq 6); radius r = √Uiso at 50 % probability coverage (c = 1.5382,
      eq 5 — inferred from their E ≈ 0.12; neither paper nor SI states the formula, the SI only
      gives the threshold 0.1541), divided by the M–L bond length. Pooled over
      ~42k CSD complexes, fitted log-normal (CN-6 μ = −2.128, σ = 0.251; CN-7 μ = −2.208,
      σ = 0.269), cutoff r_max = E + SD ≈ 85th percentile: **0.154 (CN-6), 0.145 (CN-7)**.
      Two shapes (unit sphere, RMSD-aligned with Marques et al.'s algorithm) are *thermally
      indistinguishable* if every matched vertex pair is ≤ r_max apart. A structure is accepted as
      shape P if RMSD(structure, P) < r_max (their HABLII DAC-6 0.121, Eu HECU-7 0.055).
    - **Input**: `Atom` gains `adp: Option<Matrix3<f64>>`. Parse 10-column xyz
      (`label x y z U11 U22 U33 U23 U13 U12`, SHELX/CIF order) or 5-column `Uiso`; ADPs on all atoms
      or on none; Cholesky check rejects non-positive-definite tensors. Timeo exports **U_cart
      in the xyz frame** (cctbx `adp_utils.u_cif_as_u_cart`), so cosmochlore never needs the cell.
      Plain 4-column files keep byte-identical output.
    - **Step 1, paper-faithful per-vertex test** (cheap, no sampling): `CShMResult.xyz` already holds
      the ideal shape fitted onto the structure in Å, so compare each ligand's
      |q_i − ideal_{perm(i)}| with its own 1.5382·√Uiso_i (or a Mahalanobis distance in U_i for the
      anisotropic version). Report per shape: max ratio and "within thermal ellipsoids yes/no".
      This uses the structure's own ADPs instead of the CSD-wide r_max.
    - **Step 2, `--rmax` fallback without ADPs**: same test against the paper's constants (scaled by
      mean bond length) for CN-6/7; plus an RMSD column alongside CShM, since the paper's criterion
      is RMSD, not CShM.
    - **Step 3, Monte Carlo propagation**: generic `thermal::propagate(atoms, n, seed, f)` sampling
      x_i = μ_i + chol(U_i)·z. CShM with the permutation fixed from the mean-structure search
      (per sample only Kabsch, `linalg.rs`), optional full re-search to report permutation-flip
      rate. Outputs mean/σ/P05/P95 and a *thermal floor* (the fitted ideal shape perturbed by the same
      ADPs). The same hook serves csom (fixed best axis) and odis (ζ, Δ, Σ, Θ). Analytic delta-method
      propagation is not enough on its own: ∂S/∂q = 0 at a perfect match. Caveats: CShM is biased up
      under noise (report S(mean) and the sampled mean separately); independent sampling ignores
      rigid-body correlation (TLS/rigid-bond is a later refinement). Flags `--samples`, `--seed`;
      new deps `rand`, `rand_distr` (and `rayon`, see C8).
    - **Output**: σ/P95/floor/ellipsoid columns in the table and CSV only when ADPs are present.
18. `OPEN` Thermally distinguishable shape library (same paper; SI `ao5c05878_si_001.zip`).
    - **Source**: `ii - TDPSs_SHAPE_ref_files/*.ref` (SHAPE format: symbol, name, point group, N
      vertices, centre `0 0 0` last). 12 files: DAC-6 (C2v); DPAC-7 (Cs), DTT-7 (C3v), FPSS-7 (C2v),
      HECU-7 (C2v), HEOB-7 (Cs), TSPPY-7 (Cs), TrBCSPY-7 (C2v), Δ/Λ-SHEAPR-7 (C2), Δ/Λ-THTB-7 (C2).
      Prefer these over compendium Tables S12/S13: those are 4 decimals and have typos (S12 prints
      OC-6 with 0.1 instead of 1.0).
    - **Full TDPS sets** (from the RMSD CSV headers): CN-6 = OC, TPR, PPY, DAC, HP; CN-7 = COC, CTPR,
      DTT, ETPY, HECU, HEOB, FPSS, HPY, TSPPY, PBPY, SHEAPR Δ/Λ, THTB Δ/Λ, DPAC, TrBCSPY, HP.
      cosmochlore already has OC/TPR/PPY/HP-6 (+ JPPY-6) and COC/CTPR/HPY/PBPY/HP-7 (+ JPBPY-7,
      JETPY-7). Add DAC-6 and the ten new CN-7 files. The paper's ETPY-7 has no `.ref`: check
      whether it is SHAPE's JETPY-7 before adding anything.
    - **Precision**: DAC-6 is given to only 4 decimals, and HECU-7 carries ~6e-5 off-axis residue,
      so the vertices are not exactly symmetric. Import verbatim first so the SHAPE values
      reproduce, then decide whether to symmetrise to full f64 like B6 (and note it in the yaml).
    - **Chirality**: cosmochlore's CShM sums all three singular values with no det(R) sign fix, so it
      allows improper rotations. That matches SHAPE ("does not distinguish between chiral
      enantiomers"), so Δ and Λ entries would score identically. Either ship one entry per pair, or
      add a proper-rotations-only mode (flip the smallest singular value when det(UVᵀ) < 0; the B&B
      nuclear-norm bound stays a valid upper bound) and ship both.
    - **Regression fixtures** (in the SI `.dat` files, centre first as cosmochlore expects):
      HABLII Th + 6 C → DAC-6 CShM 1.351, RMSD 0.1209; HOYGEL Eu + I + 6 N → HECU-7 CShM
      0.607, RMSD 0.0543. HOYGEL's `.dat` has the title `HABLII_0` by copy-paste; the geometry is Eu.
    - **Validation data**: `iv - TDPSs_RMSD_Spreadsheets` has per-shape RMSDs for 40,845 CN-6 and
      2,841 CN-7 CSD metal sites. It gives refcodes only, no coordinates, so it is only useful with
      CSD access (e.g. through Timeo). The CN-7 header repeats `lambda-SHEAPR-7,alpha-THTB-7`
      twice where it should be Δ/Λ, and every row has a trailing comma: map columns by position.
    - **Optional**: a `shapes --clusters` check that runs the pairwise r_max test over the built-in
      library, so thermally indistinguishable references (the K_n clusters of Figures 7–8;
      pairwise distances in compendium Tables S18/S19) are flagged in the output.

Tier 0:
16. `HALF DONE` Trajectory input + auto coordination-sphere extraction. `--center` done; multi-frame `.xyz` and `--cn <n>` cutoff remain.
    Polynator's two cheap rules cover `--cn`: *fixed number* (the n nearest atoms beyond d_min) and
    *gap* (sort the centre–ligand distances; a new sphere starts where two consecutive distances
    differ by more than 0.2·d_shortest). Its Voronoi solid-angle rule (O'Keeffe 1979, ≥ 20°) is
    heavier and not needed for molecular `.xyz` input.

## Order of work

1. ~~E6 dead-code cleanup → wire up CI (E2)~~ done (`f9cccad`…`ac5035d`)
2. ~~E5 README drift and E7 `.gitignore`~~ done (plus E8, the `ebcT-6` fixture fix it uncovered)
3. ~~B10 + rest of D7~~ done: `OptimiserSettings`, `SearchTables`/`SearchState`, `symops` type; clippy
   gate is now unconditional
4. C17, then C3 — resolve the point group once; make grouping deterministic and drop the
   per-call sort (both contained, both measurable with the 0.22 s baseline)
5. C5 (score seeds, refine best 2-3) — needs the CI gate from step 1 first
6. A4, D6 residue, D1, D2, B8, B9, C14 profile — the small residue
7. Measure 1 (paths/shape maps)
8. Rest of measure 16, then measure 14 (which retires the all-47 default), then 13
9. Measures 17–18 (ADPs, TDPS library): 17 step 1 + 18 first, they need no sampling
10. Measure 19 (dynamic shapes): 19a with the FeHS fixtures, then 19b and item 20. 19 only needs the
    existing `cshm` B&B, so it can move ahead of steps 4–9. Do 19c together with step 7's measure 1,
    and 19d last.
11. Measure 21 (21a → 21b → 21c), before any `dshm` family at CN ≥ 12. It touches only `cshm`, and
    `--match exact` keeps today's output.

Repo: [Yluro/cosmochlore](https://github.com/Yluro/cosmochlore). Line refs point at
`1ea345a` and drift as code changes — verify with grep before trusting a line number.
