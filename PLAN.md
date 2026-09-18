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
1. `OPEN` Minimal-distortion paths / shape maps (OC-6↔TPR-6 interpolation) — Alvarez 2005
2. `SHIPPED` Θ face-twist + octahedral volume — Ketkaew 2021
3. `OPEN` τ4/τ4'/τ5 geometry indices for CN=4,5 — Addison 1984 et al.
4. `OPEN` Classic distortion params (⟨λ⟩, σ², Baur D, ECoN) — Robinson 1971 et al.
5. `OPEN` Bond-valence sum — Brown & Altermatt 1985
6. `OPEN` Gyration-tensor descriptors (asphericity, κ²) — reuses `linalg.rs`
7. `OPEN` Planarity/pyramidalisation via covariance SVD

Tier 2 (moderate effort, distinctive):
8. `OPEN` Bailar/Ray-Dutt twist angles — Avdeef & Fackler 1975
9. `OPEN` Zabrodsky-Avnir CSM (folding/unfolding algorithm) — JACS 1992/93
10. `OPEN` Continuous chirality measure — Zabrodsky & Avnir 1995
11. `OPEN` Structure-to-structure CShM (`--ref` accepts `.xyz`)
12. `OPEN` Element-weighted CShM (closes B9)

Tier 3 (larger, highest scientific return):
13. `OPEN` Symmetry-adapted distortion decomposition (irrep projection using `data/pgs.rs`)
14. `OPEN` Automatic point-group detection — Nielsen 2024. D4's `todo!()` is gone, but the
    all-47-groups default is 47x the work of one group; detection is what makes it cheap.
15. `OPEN` Steric descriptors (%Vbur, cone angle) — Cavallo/Tolman

Tier 0:
16. `HALF DONE` Trajectory input + auto coordination-sphere extraction. `--center` done; multi-frame `.xyz` and `--cn <n>` cutoff remain.

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

Repo: [Yluro/cosmochlore](https://github.com/Yluro/cosmochlore). Line refs point at
`1ea345a` and drift as code changes — verify with grep before trusting a line number.
