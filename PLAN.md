# Cosmochlore Improvement Plan

Revision 3, reviewed at `4123ede`. Status legend: `FIXED` `PARTLY` `MOOT` `OPEN`
`NON-ISSUE` `NEW`. Update statuses in place as items land; don't append new sections.

Totals: 21/45 fixed, 3 partly, 2 moot, 3 non-issue, 16 open. All 56 tests pass.
Uncommitted in working tree: C1 fix (`src/csom/dev.rs`). B5/B6 already committed
(`39dd1ee`, `a2bfb6f`, `2ec9784`, `0246b87` on top of `4123ede`).

## A · Correctness

- A1 `FIXED` — test path was Windows-only. `src/odis/calc.rs:295`
- A2 `FIXED` — vacuous assertion missing `.abs()`. `src/odis/calc.rs:301`
- A3 `FIXED` — ligand-count error text mismatched count. `src/odis/mod.rs:29`
- A4 `OPEN` — `--full` writes CSVs even without `--table`. `src/odis/mod.rs:57,73`
- A5 `FIXED` — symmetry ops mislabeled ("i", "sigma_h"). `src/data/symops.rs:28-41`
- A6 `NON-ISSUE` — `.unwrap()`s in table printer can't fire (shape list provably non-empty). `src/out.rs:53,66`
- A7 `NON-ISSUE` — 2-vertex shapes unreachable (centering makes them identical). `src/xyz.rs:95`
- A8 `FIXED` — CSV fields now RFC 4180 escaped. `src/out.rs:9-17`
- A9 `FIXED` — central atom pinned in permutation search (8.7x speedup). `src/cshm/permutations.rs:39-43`

## B · Consistency

- B1 `FIXED` — `ReferenceShape` stores `Vector3` not `[f64;3]`. `src/shapes.rs:6-19`
- B2 `FIXED` — unified `Structure::atoms()` / `ReferenceShape::points()`. `src/xyz.rs:17-24`
- B3 `FIXED` — `samples`/`iterations` plain `usize`, `&str` not `&String`. `src/cli.rs:106-111`
- B4 `FIXED` — one `thiserror`-based `error::Error`, not 5 hand-rolled enums. `src/error.rs`
- B5 `FIXED` — `File::create().expect()` now propagates via `?`. `src/out.rs:92`, `src/csom/optimize.rs:56`
- B6 `FIXED` — C3/C5 trig literals (0.309, 0.951, 0.866, ...) now named f64 consts (COS_72, SIN_72, SQRT_3_DIV_2, ...). `src/data/pgs.rs:13-17`
- B7 `FIXED` — shape generator import path fixed. `generate_builtin_shapes_rs.py:130`
- B8 `OPEN` — `rotation_matrix()` transposes, `rotation_matrix_from_vector()` doesn't. `src/geometry.rs:97-115`
- B9 `PARTLY` — `csom --ignore` is a flag, `cshm` always ignores labels. `src/cli.rs:114-115`

## C · Performance (csom inner loop; ~10^6 calls/point-group search)

- C1 `FIXED` — labels double-stripped per cost eval; now uses pre-stripped `CsomStructure::labels`. `src/csom/dev.rs:261`
- C2 `FIXED` — element grouping now precomputed once in `CsomStructure::groups`, not rebuilt per op call. `src/csom/dev.rs`, `src/csom/io.rs`
- C3 `OPEN` — sort exists only to make that HashMap iteration deterministic. `src/csom/dev.rs:243`
- C4 `NON-ISSUE` — denominator is only always `n` when centroid-centered; with `--center`/`--vector` (`has_centre`) the origin isn't the centroid, so it genuinely varies. Recomputing it is O(n) against the O(n³) Hungarian assignment in the same loop, so branching on `has_centre` to shortcut it isn't worth the complexity. `src/csom/dev.rs:10-24`
- C5 `OPEN` — all 20 seeds fully refined instead of scoring first, refining best 2-3. `src/csom/optimize.rs:110`
- C6 `OPEN` — no convergence tolerance; simplex near-degenerate (0.001 edges). `src/csom/optimize.rs:77-87`
- C15 `FIXED` — `point_group_dev` (the hot-loop entry point) now sums per-operation deviations
  directly via `operation_deviation_scalar`, instead of going through
  `point_group_operation_deviations`. No more `name.to_string()` allocation or `Vec<CsomOperation>`
  collection per call; `point_group_operation_deviations`/`operation_deviation` (full breakdown,
  with name/image/pairing) is kept for the one-off `--full`/`--operated` report path. `src/csom/deviation.rs`
- C8 `OPEN` — shapes/seeds/point-groups loops are serial; `rayon` would help. `src/cshm/mod.rs`, `src/csom/optimize.rs:110`
- C9 `MOOT` — automorphism dedup removed outright (superseded, not fixed)
- C10 `FIXED` — closed-form 3x3 eigenvalues (Smith 1961). `src/cshm/linalg.rs:36`
- C11 `FIXED` — norms/suffix sums precomputed. `src/cshm/bounds.rs:25-37`
- C12 `OPEN` — `builtin_shapes()` rebuilds full 90-shape HashMap on every lookup. `src/data/standard_shapes.rs:7`
- C13 `MOOT` — superseded by C9
- C14 `OPEN` — no `[profile.release]` tuning; `license = "GPL-3"` not valid SPDX (want `GPL-3.0-only`). `Cargo.toml`
- C16 `FIXED` — resolved as a corollary of C15: `point_group_dev` no longer calls
  `point_group_operation_deviations` at all during the Nelder-Mead search, so there's no longer
  a "first, discarded" full computation for the winning axis to duplicate. The one call to
  `point_group_operation_deviations` in `calc_csom` (`src/csom/mod.rs:102-105`), gated by
  `with_operations`, is now the *only* full computation, not a redundant second one.

## D · Simplification

- D1 `OPEN` — `center_and_normalise` duplicates `center_by_centroid`+`normalise`. `src/geometry.rs:10-48`
- D2 `OPEN` — dead `points` array built only to construct `vectors`. `src/odis/calc.rs:28-45`
- D3 `FIXED` — `scale`/`centroid` now read (feed operated-coordinate writers). `src/csom/mod.rs:52-56`
- D4 `OPEN` — omitting `--pg` hits `todo!()` instead of a clap error. `src/csom/mod.rs:69-72`
- D5 `FIXED` — args parsed before banner prints. `src/main.rs:22-23`
- D6 `PARTLY` — dead `print_crab()`; `find_best_permutation`'s 4th return unused in prod. `src/main.rs:38`
- D7 `FIXED` — clippy 31→11 warnings (default target)
- D8 `FIXED` — `csom` module layout doesn't match its own pipeline; names mislead. `io.rs` isn't
  I/O (it's structure prep: centering/normalising/grouping). `dev.rs` fuses three unrelated
  things: deviation math (`sds_dev`, `point_group_dev`), a generic Hungarian-matching algorithm
  (`hungarian`, `best_permutation`) that isn't CSOM-specific, and the label-grouping assignment
  glue (`assign_group`, `best_permutation_multiple_atoms`). `mod.rs` mixes CLI entry
  (`csom_main`), orchestration (`calc_csom`), and domain types (`CsomOperationResult`, `CsomResult`,
  `CsomError`) in one file. Split done (rename only, no logic change): `types.rs` (the three
  structs/enum, out of `mod.rs`), `prepare.rs` (renamed `io.rs`), `assignment.rs` (split out of
  `dev.rs`: `hungarian`/`best_permutation`/`best_permutation_multiple_atoms`/`group_by_label`),
  `deviation.rs` (renamed/trimmed `dev.rs`: `sds_dev`/`point_group_dev`/
  `point_group_operation_deviations`), `optimize.rs` unchanged, `mod.rs` left with just
  `csom_main` + `calc_csom` + module wiring. Updated imports across `src/csom/*`, `src/out.rs`,
  `src/cli.rs`, and `src/odis/mod.rs`. All 56 tests pass; no new clippy warnings.
  `src/csom/mod.rs`, `src/csom/prepare.rs`, `src/csom/assignment.rs`, `src/csom/deviation.rs`

## E · Build, tests, CI

- E1 `OPEN` — no `lib.rs`; no integration tests possible, all tests are `#[cfg(test)]` inline
- E2 `OPEN` — CI only runs `cargo build` on tag push, no test/clippy/fmt job. **Highest leverage item.** `.github/workflows/release.yml`
- E3 `OPEN` — `bnb_is_faster_than_bf` asserts on wall-clock time. `src/cshm/tests.rs:310`
- E4 `OPEN` — `tests/FeCl6_table.csv` golden file unreferenced by any test
- E6 `OPEN` — 4 dead-code warnings: `find_automorphisms`, `automorphism_branch`, `structure_from_shape` (test-only), `twist_top_face` (fully dead). Blocks `-D warnings` in CI.
- E5 `FIXED` — README brought in line with code (`--iterations`, csom algorithm, odis caveat)

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
14. `OPEN` Automatic point-group detection (unblocks D4) — Nielsen 2024
15. `OPEN` Steric descriptors (%Vbur, cone angle) — Cavallo/Tolman

Tier 0:
16. `HALF DONE` Trajectory input + auto coordination-sphere extraction. `--center` done; multi-frame `.xyz` and `--cn <n>` cutoff remain.

## Order of work

1. E6 dead-code cleanup → wire up CI (E2): test/clippy/fmt on push+PR
2. Rest of D7 (11 remaining clippy warnings)
3. C2-C4, C15 (csom inner loop, contained/no-risk)
4. C5, C6 (seed selection, convergence — needs CI gate from step 1 first)
5. A4, D4 (remaining correctness residue)
6. B8, B9 (remaining consistency)
7. Measure 1 (paths/shape maps)
8. Rest of measure 16, then measure 13

Repo: [Yluro/cosmochlore](https://github.com/Yluro/cosmochlore). Line refs point at
`4123ede` and drift as code changes — verify with grep before trusting a line number.
