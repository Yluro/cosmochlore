# Cosmochlore Plan

Open work only. Finished items, with the invariants and baselines they set, are in
[SHIPPED.md](SHIPPED.md) — read it before touching code an item there changed. When an item
lands, delete it here and add its regression notes there.

Line references drift; grep before trusting one.

## Next up

1. **Small residue** — A4, D6, D1, D2, B8, B9, C14, D12.
2. **Measure 1** — minimal-distortion paths / shape maps.
3. **Measure 16** (rest), then **14** (retires the all-47 default), then **13**.
4. **Measures 17–18** — ADPs and the TDPS library. 17 step 1 and 18 first; they need no sampling.
5. **Measure 19** — dynamic shapes: 19a with the FeHS fixtures, then 19b and item 20. Only needs
   the existing `cshm` B&B, so it can move ahead of steps 1–4. Do 19c with measure 1; 19d last.
6. **Measure 21** — 21a → 21b → 21c, before any `dshm` family at CN ≥ 12. Touches only `cshm`;
   `--match exact` keeps today's output.

## Open fixes

| ID | Fix | Where |
|---|---|---|
| A4 | `odis --full` writes the cshm, csom and odis CSVs without `--table` | `src/out/odis.rs` |
| B8 | `rotation_matrix()` transposes, `rotation_matrix_from_vector()` doesn't | `src/geometry.rs` |
| B9 | `csom --ignore` is a flag, `cshm` always ignores labels (measure 12 closes it) | `src/cli.rs` |
| C8 | Shape / seed / point-group loops are serial; use `rayon` | `src/cshm/mod.rs`, `src/csom/` |
| C12 | `builtin_shapes()` rebuilds the 90-shape `HashMap` on every lookup | `src/data/standard_shapes.rs` |
| C14 | No `[profile.release]` (`lto = "fat"`, `codegen-units = 1`, `panic = "abort"`) | `Cargo.toml` |
| D1 | `center_and_normalise` duplicates `center_by_centroid` + `normalise`; doc comment truncated | `src/geometry.rs` |
| D2 | Dead `points` array built only to construct `vectors` | `src/odis/calc.rs` |
| D6 | `find_best_permutation`'s rotation return unused | `src/cshm/mod.rs` |
| D12 | Merge `csom`'s refinement flags (`-s`, `-i`, `-T`) into one (details below) | `src/cli.rs` |
| E1 | No `lib.rs`, so no integration tests | — |
| E3 | `bnb_is_faster_than_bf` asserts on wall-clock time | `src/cshm/tests.rs` |

### D12 — one flag for the csom refinement knobs

`csom` has three search flags, each documented as "recommended not to change": `-s/--seeds`,
`-i/--iterations` and `-T/--tolerance`. `CsomArgs::settings()` folds them into
`CsomSettings::search` (`OptimiserSettings`, B10). Merge them into a single flag. Suggested: `-r/--refine <SEEDS> <ITER>
<TOL>`, with clap `num_args = 3` and `value_names` as `--vector` does, defaulting to
`OptimiserSettings::default()` (20 200 1e-6). The maintainer may pick another name or format.

- It frees `-s` (today it means `--silent` before the subcommand and `--seeds` after it), `-i`
  (so measure 22 can be `-i/--ideal`, as in `cshm`) and `-T`.
- `settings()` builds `search` from the one flag; `OptimiserSettings::default()` stays the
  single source of the defaults.
- Output with the defaults is unchanged (`tests/output.rs`). Docs: the `csom` options table.

## Measures

Shipped: 2 (Θ + octahedral volume), 3 (`gidx`), part of 16 (`--center`) — see SHIPPED.md.

| # | Measure | Tier |
|---|---|---|
| 1 | Minimal-distortion paths / shape maps | 1 |
| 4 | Classic distortion params (⟨λ⟩, σ², Baur D, ECoN) | 1 |
| 5 | Bond-valence sum | 1 |
| 6 | Gyration-tensor descriptors | 1 |
| 7 | Planarity / pyramidalisation | 1 |
| 20 | Polynator by-products (δ, shape-free values, per-vertex deviations, hull volume) | 1 |
| 22 | Nearest symmetric structure for `csom` (`--ideal`) | 1 |
| 8 | Bailar / Ray-Dutt twist angles | 2 |
| 9 | Zabrodsky-Avnir CSM | 2 |
| 10 | Continuous chirality measure | 2 |
| 11 | Structure-to-structure CShM | 2 |
| 12 | Element-weighted CShM | 2 |
| 19 | Dynamic (deformable) shape measures — `dshm` | 2 |
| 21 | Belt matching for large CN in `cshm` | 2 |
| 13 | Symmetry-adapted distortion decomposition | 3 |
| 14 | Automatic point-group detection | 3 |
| 15 | Steric descriptors (%Vbur, cone angle) | 3 |
| 17 | ADP-aware measures | 3 |
| 18 | Thermally distinguishable shape library | 3 |
| 16 | Trajectory input + `--cn` coordination-sphere extraction | 0 |

## Measure specs

Agent briefs. Decisions marked "(decided)" are settled; don't reopen them.

### Tier 1 — cheap, high value, reuses existing machinery

- **1 ·** Minimal-distortion paths / shape maps (OC-6↔TPR-6 interpolation) — Alvarez 2005.
   Fitting a structure *to* a path (the path coordinate plus the deviation from the path) is a
   one-variable dynamic shape: do it as 19c.
- **4 ·** Classic distortion params (⟨λ⟩, σ², Baur D, ECoN) — Robinson 1971 et al.
- **5 ·** Bond-valence sum — Brown & Altermatt 1985
- **6 ·** Gyration-tensor descriptors (asphericity, κ²) — reuses `linalg.rs`
- **7 ·** Planarity/pyramidalisation via covariance SVD
- **20 ·** Polynator by-products — Link & Niewa 2023 (see 19), manual §3.5, §3.6, §5. Each is O(n) on
    results cosmochlore already has:
    - δ = 10·√S (`--delta`) as an extra column for `cshm` and for 19. δ grows linearly with small
      distortions while S grows quadratically, and 2 decimals of δ always separate "exact" from
      "slightly distorted". Display only: never change what is minimised or stored.
    - Shape-free values from the eigenvectors of Σ q qᵀ (centroid-centred): δ_linear (distance to the
      best line), δ_planar (to the best plane) and δ_spherical (spread of |q_i| about its mean). The
      same eigen-decomposition as items 6 and 7, so build all three together.
    - Per-vertex deviations of the fitted ideal shape (Polynator's `.aso` file): |q_i − v_i| in Å, its
      mean and SD, split into radial (|q_i| − |v_i|) and angular parts. `CShMResult::ideal` (from
      `calc_cshm(.., with_ideal)`) holds v_i, atoms in input order. Item 17 step 1 tests the same vector against the ADPs.
    - Convex-hull volume and surface area of the structure and of the fitted ideal shape for any CN
      (`odis` only has the octahedral volume). Quickhull is enough for n ≤ 60.
    - Ranking with a parameter penalty (Polynator's "variable tax": δ + 1.5 per free variable, used
      only to order the list, never reported). Only matters once 19 mixes models with different dof.
- **22 ·** Nearest symmetric structure for `csom` (`--ideal`): the folding/unfolding average of the
    paired images (Salomon & Avnir; `CSM_PLAN.md` §2.2), the "average Q̄ of all images ÔRQ" of the
    docs' CSoM theory page (Averaging). Every input already exists when operations are computed:
    `CsomOperationResult` holds `image` (refined, normalised frame) and `pairing`, and
    `CsomResult::operated` holds the same structures in the input frame: `[0]` is the input (E),
    then one per operation. `operated[g]` is g·q_i under atom i's own label, *not* yet paired, so
    Q̄ is not the plain atom-by-atom mean of `operated`: atom i averages `operated[g]` atom
    π_g(i) over g. The map back to the input frame is affine and the same for every g, so that
    average equals Rᵀ Q̄ / scale + centroid of the normalised-frame Q̄ below.
    - **Math**: work in the refined, normalised frame, with points q_i. For each table operation g,
      π_g is `pairing`, so the operation is scored on Σ_i |q_i − g q_{π_g(i)}|². Let G be {E} plus
      the table. `pgs.rs` leaves E out of every table except `E`'s own, and π_E = id. Don't count E
      twice. The nearest structure is Q̄_i = (1/|G|) Σ_{g∈G} g q_{π_g(i)}. If the π_g form a group
      action, (T_g Q)_i = g q_{π_g(i)} is an orthogonal representation, and (1/|G|) Σ T_g is the
      orthogonal projector onto the G-symmetric configurations. Q̄ is then exactly G-symmetric, and
      d² := 100·Σ|q_i − Q̄_i|²/D = CSoM·(|G|−1)/(2|G|). Here D is the `sds_dev` denominator, the same
      for every operation, and CSoM is the code's average over the |G|−1 non-identity operations.
      The factor is 1/4 for C2, Ci and Cs and 47/96 for Oh. For `E`, Q̄ = Q and d² = 0. Checked in
      Python on a noisy C3v structure: CSoM 1.045089, d² 0.435454, Q̄ symmetric to 5e-16.
    - **Consistency**: the group-action condition is π_{gh}(i) = π_h(π_g(i)) for all g, h ∈ G and
      all i. The swapped order π_g(π_h(i)) is wrong (checked numerically on C3v). Each operation
      runs its own Hungarian matching, so nothing enforces this. It holds for near-symmetric
      structures but can fail for strongly distorted ones or with `--ignore`. When it fails, the
      average is not symmetric, and CSoM·(|G|−1)/(2|G|) is only a lower bound on the distance to
      the nearest symmetric structure.
      - Check: build the Cayley table once per group. For each g, h, find the k with
        max|M_k − M_g M_h| < 1e-8 (the product may be E). The tables are exact f64 and closed
        (`check_pointgroup_tables.py`), so a missing product is a table bug: return a `CsomError`
        and don't fall back. Cost: |G|² products (14 400 for Ih), then O(|G|²·N) comparisons.
      - On failure: warn, naming the group and the first failing pair (g, h). Write no file for
        that group and carry on with the others. Never write a non-symmetric "ideal". For any
        input, the true nearest structure is `csm --ideal` (`CSM_PLAN.md` P3.1).
    - **Where**: new `src/csom/ideal.rs` with
      `nearest_structure(points, operations, pg) -> Result<Nearest, CsomError>`, computing Q̄ and d²
      in the refined, normalised frame. Suggested shape:
      `Nearest = Symmetric { structure: Structure, d2: f64 } | Inconsistent { g, h }`, where
      `structure` is Q̄ mapped back to the input frame with `Structure::with_coords`, as `calc_csom`
      builds `operated`, so it carries the input's labels like `CShMResult::ideal`. Call it in
      `calc_csom` where `operations` are built. Store `nearest: Option<Nearest>` in `CsomResult`,
      `None` when operations aren't computed. `odis` leaves `CsomSettings::operations` false, so
      `odis` is untouched.
    - **CLI**: `--ideal`, long-only until D12 frees `-i` (today `--iterations`); then
      `-i/--ideal`, as in `cshm`. In `CsomArgs::settings()`, set `operations: details || operated ||
      ideal`.
    - **Output**: one `<stem>_<PG>_ideal.xyz` per group, like `--operated`. `.gitignore`'s
      `*_ideal.xyz` already covers it, and the name can't collide with `cshm`'s `<stem>_ideal.xyz`.
      Atoms in `Structure::atoms()` order with their own labels. Comment line:
      `<PG> ideal csom = <S> d2 = <d²>`, 3 decimals as in `--operated`. The pinned centre atom has
      π_g(0) = 0 and sits at the origin, so it is written at its input position. With `--ignore`,
      Q̄_i can average atoms of different elements. Still write it under atom i's label, and say so
      in the README.
    - **Naming**: don't call d² the CSM in output or the README (`CSM_PLAN.md` §0 rule 1). When the
      pairings are consistent and the axis search found the global minimum, d² equals the
      Zabrodsky–Avnir CSM for the same centre. Add that cross-check test once `csm` ships.
    - **README**: in the `csom` section, add the flag row and output name. State the measure with
      the normalisation the code uses: S_G(Q) = min_R (1/(|G|−1)) Σ_{Ô∈G, Ô≠E} S_Ô(RQ), and S_E = 0.
      Then state d² = S_G·(|G|−1)/(2|G|). Equivalently, d² is exactly half of the average over all
      |G| operations with E included.
    - **Tests** (`src/csom/tests.rs`):
      1. `tests/ideal_water.xyz` (C2v, Cs, C2) and `tests/ideal_benzene.xyz` (`-n`; D6h, C6, C2v):
         Q̄ = Q to 1e-10 and d² < 1e-10.
      2. Distorted inputs (`FeCl6.xyz`, `La03.xyz`) for C2, Ci, Cs, C3v, D4h and Oh, wherever the
         pairings are consistent. Assert d² = CSoM·(|G|−1)/(2|G|) to 1e-9 relative. Also assert
         `point_group_dev(Q̄)` < 1e-10 at the identity rotation, since Q̄ is already in the refined
         frame.
      3. `E`: Q̄ = Q exactly and d² = 0. This guards against counting E twice.
      4. Inconsistent: find a seeded random structure that fails the check (e.g. 8 random points
         against Oh) and pin it as a fixture. Assert `Inconsistent` and that no file is written.
      5. Composition order: a non-abelian case (C3v) where the swapped rule π_g∘π_h gives a
         different answer, so a check written in the wrong order fails the test.
      6. Output of `csom` without `--ideal`, and of `cshm` and `odis`, stays byte-identical
         (as in D8/D9).


### Tier 2 — moderate effort, distinctive
- **8 ·** Bailar/Ray-Dutt twist angles — Avdeef & Fackler 1975
- **9 ·** Zabrodsky-Avnir CSM (folding/unfolding algorithm) — JACS 1992/93
- **10 ·** Continuous chirality measure — Zabrodsky & Avnir 1995
- **11 ·** Structure-to-structure CShM (`--ref` accepts `.xyz`)
- **12 ·** Element-weighted CShM (closes B9)
- **19 ·** Dynamic (deformable) shape measures — Link & Niewa, *J. Appl. Cryst.* 2023, 56, 1855
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
- **21 ·** Belt matching for large CN in `cshm` — Polynator's assignment heuristic (Link & Niewa 2023
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
      is in `src/cli.rs` and `CsomArgs`. `-m` is still free in `cshm`.
      - `auto`, the default: `exact` when CN < 12, `belt` when CN ≥ 12. CN is the number of
        ligand vertices, `structure.ligands.len()` in `cshm_main`, with the centre not counted.
      - `exact`: today's B&B at any CN, unchanged. Its help text should say that it can run for
        minutes or hours above CN ≈ 20.
      - `belt`: steps 1–6 at any CN. Below 12 it is mainly for tests and comparisons.

      At CN = 12, `auto` switches to `belt`. In every CN-12 measurement above, the plain B&B
      finished in well under the budget, so the certificate completes and the values equal
      today's. The permutation can still be a different, symmetry-equivalent one, which moves
      atoms to other vertices in the `--ideal` file. Byte-identical output is therefore guaranteed for
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
      - 21c: `-m/--match`, `--nodes`, `CshmSettings`, the output marker, the CSV column, a docs
        section in `content/Modules/cshm.md` ("The `--match` modes"), and new output patterns in `.gitignore` if any (E7).
    - **Not adopted from Polynator**: hand-written belts per model (derived automatically here);
      the q2 angular estimate (exact scoring is cheap); strategies 2 and 3 for prolate and oblate
      models, which built-in shapes above CN 12 do not need.


### Tier 3 — larger, highest scientific return
- **13 ·** Symmetry-adapted distortion decomposition (irrep projection using `data/pgs.rs`). The
    orbit projectors of 19a are the same machinery, and the drop in S from parent to child family
    measures the distortion that the child's symmetry allows.
- **14 ·** Automatic point-group detection — Nielsen 2024. D4's `todo!()` is gone, but the
    all-47-groups default is 47x the work of one group; detection is what makes it cheap.
- **15 ·** Steric descriptors (%Vbur, cone angle) — Cavallo/Tolman
- **17 ·** ADP-aware measures (thermal significance) — Munguba … Simas, *ACS Omega* 2025, 10,
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
    - **Step 1, paper-faithful per-vertex test** (cheap, no sampling): `CShMResult::ideal` holds
      the ideal shape fitted onto the structure in Å, atoms in input order, so compare each
      ligand's |q_i − ideal_i| with its own 1.5382·√Uiso_i (or a Mahalanobis distance in U_i for the
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
- **18 ·** Thermally distinguishable shape library (same paper; SI `ao5c05878_si_001.zip`).
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


### Tier 0 — input handling
- **16 ·** Trajectory input + auto coordination-sphere extraction. `--center` done; multi-frame `.xyz` and `--cn <n>` cutoff remain.
    Polynator's two cheap rules cover `--cn`: *fixed number* (the n nearest atoms beyond d_min) and
    *gap* (sort the centre–ligand distances; a new sphere starts where two consecutive distances
    differ by more than 0.2·d_shortest). Its Voronoi solid-angle rule (O'Keeffe 1979, ≥ 20°) is
    heavier and not needed for molecular `.xyz` input.

Repo: [Yluro/cosmochlore](https://github.com/Yluro/cosmochlore).
