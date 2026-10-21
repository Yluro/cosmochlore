# AGENTS.md

Rules for AI agents and contributors working on cosmochlore. Read [README.md](README.md) first: it
says what the program does and where the code lives.

## Before changing code

- Open work is in [PLAN.md](PLAN.md). Shipped work, with the invariants and baselines it set, is in
  [SHIPPED.md](SHIPPED.md). Read the SHIPPED.md entries for any code you touch and re-run the
  measurements they cite. [CSM_PLAN.md](CSM_PLAN.md) is the self-contained brief for the planned
  `csm` module.
- When a PLAN.md item lands, delete it from PLAN.md, add its regression notes to SHIPPED.md, and
  name its ID (`C17`, `D10`, ...) in the commit message.
- Line references in the plan files drift. Grep before trusting one.
- The default branch is `master`.

## Build, test, lint

CI (`.github/workflows/ci.yml`) runs these with `RUSTFLAGS=-D warnings`, so any rustc warning
(dead code, unused import) fails the build:

```bash
cargo test --locked
cargo clippy --all-targets --locked -- -D warnings
cargo fmt --check
```

- Clippy stays at 0 warnings. Fix the code instead of adding `#[allow]`s.
- Timings and the example output on the docs site come from the release binary
  (`cargo build --release`), not a debug build.
- `rustfmt.toml` sets `newline_style = "Native"`, so Windows checkouts are CRLF. Don't convert
  line endings. The generated tables in `src/data/` carry `#[rustfmt::skip]`.

## Generated and checked data

- `src/data/standard_shapes.rs` is generated. Edit `src/data/shapes/*.yaml`, or the generator for
  the `ReferenceShape` type, then regenerate:

  ```bash
  py generate_builtin_shapes_rs.py src/data/shapes src/data/standard_shapes.rs
  ```

- After editing `src/data/pgs.rs`, run `py check_pointgroup_tables.py`. Matrices are written to
  full f64 precision, the identity is implicit, and the principal axis lies along `z`, which the
  `csom` seeding assumes.
- The Python scripts need Python 3 and numpy. The visualisers also need matplotlib, and
  `visualize_csom_landscape.py` needs scipy. Each script's docstring gives its usage.

## Code rules

- Never print a number when something went wrong. Every failure returns the module's error type
  (`thiserror`), wired into `crate::Error` in `src/error.rs`. Don't `unwrap`/`expect` on anything
  that comes from user input.
- Only `src/out/` prints to stdout or writes files, and stdout only through `say!(silent, ...)`
  with the `silent` flag passed down from `main`. Computation modules return results and never
  import `crate::out`.
- The output of `cshm`, `csom`, `odis` and `gidx` must stay byte-identical unless the change is
  meant to alter it. `tests/output.rs` checks it against `tests/expected/`; after an intended
  change, rerun it with `UPDATE_EXPECTED=1` and review the diff. Baselines are in SHIPPED.md
  (for example, `csom tests/FeHS.xyz -p Oh D4h` gives Oh 5.380 and D4h 5.611 in about 0.22 s).
- `cshm` keeps the centre pinned in the permutation search (SHIPPED.md A9).
- Tests must be path-portable and must not assert on wall-clock time. Count search nodes instead.
- Don't add a dependency without giving the reason in the PR.
- Output files are named `<input stem>_<kind>` and written next to the input. A new output kind
  needs a matching pattern in `.gitignore`.
- Dev-only files (helper scripts, plans) are marked `export-ignore` in `.gitattributes` so they stay
  out of source archives. Add new ones there.
- Keep doc comments and comments short: one or two lines saying what the item does and any
  non-obvious contract. History, rationale and plan IDs go in the commit message or SHIPPED.md,
  not the code.

## CLI changes and documentation

User documentation lives in a separate repository,
[Yluro/cosmochlore-docs](https://github.com/Yluro/cosmochlore-docs) (a Quartz site, pages under
`content/`), published at https://yluro.github.io/cosmochlore-docs. Keep the README barebones and
don't add user documentation to it.

When a change touches a flag, a default, the printed output or a computed value, update the docs
in the same piece of work, or say in the PR that they need updating:

- `content/Modules/<command>.md`: the options table, which must match `src/cli.rs`, and the
  examples, which are verbatim release-binary output. Regenerate them when a printed value changes.
- `content/getting-started/usage.md`: the command list and the output files table.
- `content/reference/`: built-in shapes, point groups and the YAML shape format.
- `content/theory/`: definitions and algorithms.

## Adding a subcommand

1. An args struct and a `Command` variant in `src/cli.rs`, and a `dispatch` arm in `src/main.rs`
   that hands the module's result to `out`.
2. A module directory whose `mod.rs` holds `run(&Args) -> Result<XRun, Error>`, which computes and
   prints nothing, with its tests in `tests.rs` behind `#[cfg(test)]`.
3. An error enum, wired into `crate::Error`.
4. A `src/out/<command>.rs` with `pub fn <command>(&run, &args, silent)`, re-exported from
   `src/out/mod.rs`. Build terminal tables row by row with `Table` / `Row` from
   `src/out/table.rs` (fixed widths are fine), and CSVs with `src/out/csv.rs`, listed by hand.
5. `.gitignore` patterns for its output files, a docs page under `content/Modules/`, a row in the
   README command table, and a case in `tests/output.rs`.

## Releases

Pushing a `v*` tag runs `.github/workflows/release.yml`. It builds the Linux, Windows and both
macOS binaries and attaches them to a GitHub release, with a changelog of the commits since the
previous tag. Bump `version` in `Cargo.toml` first. Don't tag or release unless asked.
