# Shipped changes — regression notes

Everything here has landed. It is kept for agents: before changing code these items touch,
check that the invariants below still hold, and re-run the measurements they cite. Open work
lives in [PLAN.md](PLAN.md). Item IDs match PLAN.md and old commit messages.

Line references were taken at `1ea345a`; the `cargo fmt` pass in `e163603` and later commits
moved them, so grep before trusting one.

## Baselines

- csom whole-process on `tests/FeHS.xyz -p Oh D4h`: **0.22 s** (was 3.86 s before C6), Oh 5.380, D4h 5.611.
- `odis --full` on FeHS: ~375 ms (was 600 ms before B10).
- cshm 11-vertex B&B timing: 404 ms min (unchanged by D7).
- Clippy: **0 warnings** with `--all-targets -D warnings`, no `allow`s added to silence it.
- The example blocks on the docs site (`content/Modules/*.md` in `Yluro/cosmochlore-docs`, moved
  there from the README) are verbatim release-binary output on the `tests/` fixtures
  (`FeHS.xyz` + `ebcT-6.yaml`, ...). Regenerate them if any printed value changes.

## Correctness

- **A1** — a test path was Windows-only; tests must stay path-portable. `src/odis/calc.rs`
- **A2** — a vacuous assertion was missing `.abs()`. `src/odis/calc.rs`
- **A3** — ligand-count error text now matches the count. `src/odis/mod.rs`
- **A5** — symmetry ops "i" and "sigma_h" were mislabeled. `src/data/symops.rs`
- **A8** — CSV fields are written unquoted. Instead, the parsers reject the characters that
  would need quoting (`xyz::FORBIDDEN_CHARS`: comma, double quote): `.xyz` atom labels
  (`XyzParseError::ForbiddenLabel`) and the YAML `symbol`, `symm` and `name`
  (`YamlParseError::ForbiddenCharacter`). Every other text field is built in (shapes, point
  groups, operations) and has none. Replaced the RFC 4180 `csv_field` escaping. A new text
  field from user input must go through the same check. `src/xyz.rs`, `src/yaml.rs`
- **A9** — the central atom is pinned in the permutation search (8.7× speedup). Any change to
  `src/cshm/permutations.rs` must keep the centre pinned; see E8 for what unpinning hid.
- **E8** — `tests/ebcT-6.yaml` had its `centre` and fourth vertex swapped. Fixed in the fixture and
  the README yaml example. Expected: README inline structure → 14.335, `tests/FeHS.xyz` → 14.277.
  A mis-specified (swapped) shape scores 35.8 with the centre pinned.

## Consistency

- **B1** — `ReferenceShape` stores `Vector3`, not `[f64;3]`. `src/data/standard_shapes.rs`
- **B2** — one accessor style: `Structure::atoms()` / `ReferenceShape::points()`. `src/xyz.rs`
- **B3** — `samples`/`iterations` are plain `usize`; `&str` not `&String`. `src/cli.rs`
- **B4** — a single `thiserror`-based `error::Error` replaces 5 hand-rolled enums. `src/error.rs`
- **B5** — `File::create` errors propagate via `?`, no `.expect()`. `src/out.rs`, `src/csom/optimize.rs`
- **B6** — C3/C5 trig literals are named f64 consts (`COS_72`, `SIN_72`, `SQRT_3_DIV_2`, …). `src/data/pgs.rs`
- **B7** — shape generator import path fixed. `generate_builtin_shapes_rs.py`
- **B10** — `csom::types::OptimiserSettings { seeds, iterations, tolerance }`, `Default` =
  `(20, 200, 1e-6)`, is the single source of these defaults. `cli.rs` reads its
  `default_value_t`s from it (`CsomArgs::search_settings()`), `calc_csom`/`search_best_axis` take
  it, and `odis --full` passes `OptimiserSettings::default()`. odis' eight deviations are unchanged
  to 3 decimals vs the old `(20, 1000, 1e-8)`; rotation matrices move in the 4th. New settings
  structs (`DshmSettings`, `CshmSettings`) should follow this pattern.

## Performance

- **C1** — cost evals use pre-stripped `CsomStructure::labels`. `src/csom/deviation.rs`
- **C2** — element grouping precomputed once in `CsomStructure::groups`. `src/csom/prepare.rs`, `src/csom/assignment.rs`
- **C6** — `--tol` (default 1e-6) wired to `NelderMead::with_sd_tolerance`; default `--iterations`
  1000 → 200; starting simplex edge 0.001 → 0.1. Gave the 3.86 s → 0.22 s baseline with identical
  values. `src/csom/optimize.rs`, `src/cli.rs`
- **C10** — closed-form 3×3 eigenvalues (Smith 1961). `src/cshm/linalg.rs`
- **C11** — norms and suffix sums precomputed for the B&B bound. `src/cshm/bounds.rs`
- **C15** — `point_group_dev` (hot loop) sums per-operation deviations via
  `operation_deviation_scalar`: no `String` or `Vec<CsomOperation>` per call. The full breakdown
  (`point_group_operation_deviations` / `operation_deviation`) is only for `--full`/`--operated`.
- **C17** — the point group is resolved once per `search_best_axis` call by `pointgroup_matrices`
  (a `Vec<Matrix3>`). `point_group_dev` and `refine_axis_from_seed` take `&[Matrix3<f64>]`, and
  `point_group_dev` returns `f64`, so no name lookup, `to_matrix3` or `expect` is left in the cost
  loop. `get_pointgroup_map` was removed. Output is byte-identical (`csom tests/FeHS.xyz -p Oh D4h
  C2v --full`, `tests/La03.xyz -p Oh Ci`). Wall time is unchanged within noise (~0.60 s release on
  the dev machine before and after), so the remaining cost is in the assignment (C3) and the seed
  loop (C5). `src/csom/deviation.rs`, `src/csom/optimize.rs`
- **C3** — label groups come from a `BTreeMap` (label order) instead of a `HashMap`, and the
  per-call z/y/x sort of the matched pairs is gone. `best_permutation_multiple_atoms` returns only
  the `pairing`; `sds_dev(reference, problem, pairing)` scores `reference[i]` against
  `problem[pairing[i]]` with every sum in atom order, so the score doesn't depend on the group
  order (`group_order_does_not_change_the_score` pins this to the bit). `best_permutation` returns
  just the assignment. Any new scorer must sum in atom order through the pairing, never in group
  or match order. Printed values are unchanged on the csom/odis/cshm fixtures (FeHS, FeCl6, La03,
  CdONCl3, ML5, sample, ideal_water, ideal_benzene; all 47 groups, `--full`, `--operated`,
  `--ignore`). The one byte that moved is the Ci rotation in `La03_csom_table.csv`: the inversion
  commutes with every rotation, so the Ci (and E) axis is arbitrary and follows last-bit noise.
  Wall time ~15 % lower on the dev machine (FeHS `-p Oh D4h` 545 → 460 ms, La03 all 47 groups
  8.96 → 7.83 s). `src/csom/assignment.rs`, `src/csom/deviation.rs`
- **C16** — the one `point_group_operation_deviations` call in `calc_csom` (gated by
  `with_operations`) is the only full computation. `src/csom/mod.rs`
- **C14 (part)** — `license = "GPL-3.0-only"` in `Cargo.toml` (`734c0a9`). Release profile still open.

## Simplification

- **D3** — `scale`/`centroid` in `src/csom/types.rs` feed the operated-coordinate writers.
- **D4** — omitting `--pg` analyses all 47 supported point groups (`bb0a88e`). Brute-force stand-in
  for automatic detection (measure 14), which should retire it.
- **D5** — args are parsed before the banner prints. `src/main.rs`
- **D6 (part)** — `print_crab()` is live behind the hidden `cshm --crab` flag.
- **D7** — clippy at 0. `branch` takes `&SearchTables` (both point sets + `hi`/norm/suffix tables)
  and `&mut SearchState` (partial permutation + best so far) instead of 11 arguments;
  `find_automorphisms` reuses `SearchTables::new(reference, reference)`;
  `symops::get_operation` returns `Option<&[SymmetryOperation]>`. cshm output byte-identical on
  FeHS and La03.
- **D8** — `csom` module split: `types.rs`, `prepare.rs` (was `io.rs`), `assignment.rs`
  (Hungarian + label grouping), `deviation.rs` (was `dev.rs`), `optimize.rs`, `mod.rs` with
  `csom_main` + `calc_csom`. Pure rename.
- **D9** — `shapes.rs` → `cshm/shape_lookup.rs` (lookup by vertex count/index,
  `check_vertex_count`, `ShapeLookupError`). `ReferenceShape` lives in
  `data/standard_shapes.rs`, emitted by `generate_builtin_shapes_rs.py`; this removed the
  `shapes` ↔ `data` import cycle. `yaml.rs` stays top-level. Test helper `structure_from_shape`
  is in `cshm/test_utils.rs`. Pure move, all outputs byte-identical.
- **`--silent` (`52188d9`)** — top-level `-s` (before the subcommand) silences banner, tables,
  progress lines and "Program finished"; errors stay on stderr; CSV/xyz files are still written.
- **D10** — `--silent` is a plain parameter: `main` passes `cli.silent` to `dispatch`, the `out`
  entry points and the banner/finished lines, and from there to every printer and writer
  (`create_output`, `Table::print`, `write_csv`, the xyz/mol2 writers). `say!(silent,
  ...)` is `println!` behind that flag. The `SILENT` `AtomicBool` and `set_silent()` are gone.
  `tests/output.rs` passed unchanged: same stdout, and empty stdout with the same files under `-s`.
- **D11** — output owned by `out`. `cshm::run`, `csom::run`, `odis::run` and `gidx::run` return
  `Vec<CShMResult>`, `CsomRun`, `OdisRun` and `GidxResult`; `main::dispatch` hands them with the
  args to `out::cshm` / `csom` / `odis` / `gidx`. No module imports `crate::out`. `src/out/` is
  one file per command (`cshm.rs`, `csom.rs`, `odis.rs`, `gidx.rs`, each `pub fn <command>`
  re-exported by `mod.rs`, with that command's xyz/mol2 writers), plus `table.rs`, `csv.rs` and
  `mod.rs` (`say!`, banner, `output_stem`, `create_output`, `table_path`). `odis --full` reuses
  `cshm::print` / `write_table_csv` and `csom::print` / `write_table_csv`.
  - Tables (decided): `Table::new(&[Width::Fixed(n) | Width::Dynamic])`, then `Row::Text`,
    `Row::Rule`, `Row::DoubleRule` and `Row::from([cells])` in order, `print(silent)`. Dynamic is
    the longest cell, header included, + 2; rules span the table; trailing spaces trimmed. Each
    writer lists its rows by hand with fixed widths as magic numbers; right-aligned values are
    padded by the caller (`{:>12.4}`). A row with the wrong cell count is a bug:
    `debug_assert` only. Two earlier attempts (`Field` / `Col` describing table and CSV at
    once, and the D10 `Table` builder / `CsvWriter`) were dropped as too abstract.
  - CSVs are not derived from tables: headers, decimals and layout differ (odis/gidx CSVs are one
    row of keys and one of values). `csv::write_csv(path, what, header, rows, silent)` is a
    stopgap until the CSV writer is designed.
  - `CShMResult` holds `ideal: Option<Structure>`, filled only with `--ideal`
    (`calc_cshm(.., with_ideal)`; `odis --full` passes false). `perm`, `xyz` and `CshmRun` are
    gone. `<file>_ideal.xyz` now lists atoms in input order (centre first, as
    `Structure::atoms()`), not reference-vertex order: same lines, reordered.
  - Every other CSV, xyz and mol2 file is byte-identical to before (release binary, the
    `tests/output.rs` cases plus La03 `-n -t`). stdout changed on purpose: no trailing padding,
    a blank line before every "Input file:" (gidx had none), and the rules span the table, so
    the odis/gidx rules are 37 wide and the cshm/csom ones shift by one or two. The docs
    examples need regenerating.
  - `odis --full` now computes everything before printing, so a csom failure no longer leaves
    tables on stdout ahead of the error.
  - `tests/output.rs`: golden stdout and files for 10 cases in `tests/expected/<case>/`
    (`stdout.txt`, `files.txt`, file contents), plus the same run under `-s` (empty stdout, same
    files). csom files (`_csom_table.csv`, `_details.csv`, `_operated.xyz`, `_merged.mol2`) are
    checked by name only, since libm may move a last digit across platforms. `.gitignore`
    re-includes `tests/expected/**` and `.gitattributes` pins it to LF. Goldens are not
    `export-ignore`d, so `cargo test` works from a source archive. ~20 s in debug: keep CN ≥ 12
    cases out (La03 `-n` alone took minutes).

## Build, tests, CI

- **E2** — `.github/workflows/ci.yml`: Test, Clippy and Rustfmt jobs on push to master and PRs,
  `RUSTFLAGS=-D warnings`, `--locked`. One-time `cargo fmt` in `e163603` (blame-ignored).
  Generated tables in `data/pgs.rs` and `data/standard_shapes.rs` carry `#[rustfmt::skip]`, and
  the shapes generator emits it. `rustfmt.toml` pins `newline_style = "Native"` (CRLF checkout).
- **E5** — README rows match `cli.rs` (`--iterations` 200, `--pg` optional, `-e`/`--tol`). Keep
  the docs options tables (`content/Modules/*.md`) in sync with `cli.rs` on every CLI change.
- **E6** (`f9cccad`) — `automorphism.rs` and `structure_from_shape` are `#[cfg(test)]`;
  `twist_top_face` deleted. Measure 19 will need `find_automorphisms` un-gated.
- **E7** — `.gitignore` covers `.idea/` and every output name `out.rs` writes
  (`*_cshm_table.csv`, `*_csom_table.csv`, `*_odis_table.csv`, `*_details.csv`, `*_operated.xyz`,
  `*_merged.mol2`, `*_ideal.xyz`). Add new output patterns whenever a writer is added.

## Measures

### 2 · Θ face-twist + octahedral volume (Ketkaew 2021)

Part of `odis`; printed as the `Theta`/`Volume` rows in the `odis` docs example.

### 3 · Geometry indices — `gidx` (τ4, τ4', τ5, τ6, τ6', τ8, τ8')

`src/gidx/`: `mod.rs`, `tests.rs`, `calc/{tau4,tau5,tau6,tau8}.rs`. Sources: Addison 1984,
Yang 2007, Okuniewski 2015, Stoeckli-Evans 2025, Turnbull/Wetmore/Gerken 2021.

Fixtures and pinned values:
- Ideal built-in shapes; the docs page `content/theory/geometry-indices.md` lists every 4-, 5-, 6-
  and 8-vertex ideal shape's values.
- `tests/ML4.xyz`, `tests/ML5.xyz`, `tests/FeHS.xyz`.
- `tests/CdONCl4.xyz` (τ6 0.27), `tests/CdONCl3.xyz` (τ5 0.31) — the paper's Cd centres from its SI.
- `tests/cd6/`: CSD Cd sites IQATAY Cd1-4, ULESAH, VAPSAG, CIJBII.
- `tests/tau8/`: reproduces SI Table S10 for 8 structures within 0.01.

Decisions to preserve:
- Exact tetrahedral angle arccos(−1/3), not 109.5°/141°: an ideal tetrahedron gives exactly 1.
- CN 4/5 use the two largest angles. CN 6 uses three trans angles, chosen as the disjoint pairing
  of the ligands with the greatest angle sum. A pentagonal pyramid gives 0.900 (the paper's 1.00
  sums the five adjacent base angles instead).
- τ6' is cosmochlore's own: (5×90° − sum of the five smallest angles)/90°. OC 0, equal-edge prism
  0.456, pentagonal pyramid 1, hexagon 5/3. Five angles, not three (VAPSAG/CIJBII read
  1.044/0.995 with five, 1.31/1.18 with three).
- τ8 (from the two square-face diagonals) and τ8' (dihedral, A and B planes) follow SI 2.1.1/2.2.1;
  face 2 is the complement of face 1. Kepert's exact values are used: 129.8208/97.0714 (root of
  12w³−7w²−2w+1, w = cos²tA) and arccos(√2−1) = 65.5302. `cube_delta` = τ8' − τ8.
- Known unexplained gaps vs the SI: the dpb chelate (1.10 vs 1.22) and one dpe state.

### 16 (part) · `--center`

Arbitrary centre position via `--center`. Multi-frame xyz and `--cn` are still open.

## Investigated — leave alone

- **A6** — the `.unwrap()`s in the table printer can't fire (shape list provably non-empty).
- **A7** — 2-vertex shapes are unreachable (centering makes them identical). `src/xyz.rs`
- **C4** — the csom deviation denominator is only constant when centroid-centred; with
  `--center`/`--vector` it varies. Recomputing is O(n) next to an O(n³) Hungarian; don't branch on it.
- **C5 (tried, failed)** — score all 20 seeds, refine only the best 2–3. Tried on a branch that was
  deleted long ago. Seeds that score badly at the start often refine into the global minimum, so
  the shortcut missed minima, and the speed-up was small. Don't retry it. Any seed-loop speed-up
  must keep every global minimum the full 20-seed refinement finds. `src/csom/optimize.rs`
- **C9 / C13** — automorphism dedup was removed outright.
- **E4** — `tests/FeCl6_table.csv` and `FeCl6_ideal.xyz` deleted (`7354442`). The SHAPE 2.1 values
  are pinned inline in `src/cshm/tests.rs` (water 0.035/33.342, Eu7 2.109/2.630/7.288).
