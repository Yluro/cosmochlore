"""
Generate the fixtures for PLAN.md item 21 (`cshm -m/--match`, Polynator-style belt matching):

  tests/belt/<SHAPE>_<NN>.xyz   distorted copies of built-in shapes from CN 10 to 60 (centre M at
                                the origin, ligands X), with Gaussian noise of NN % of the radius
                                and the ligands shuffled. Fixed seed: the files are reproducible
                                and are the ones the PLAN.md timings were measured on.
  tests/belt_reference.csv      for every structure (those files plus tests/FeHS.xyz and
                                tests/La03.xyz) against every built-in shape of its CN:
                                  layers       the reference belts (layer sizes along the axis)
                                  s_belt       this script's reference belt matching
                                  s_exact      cosmochlore's exact branch-and-bound (cshm -s),
                                               3 decimals as printed; exact_status tells whether
                                               it ran (ok), hit --timeout, or was skipped above
                                               --max-exact-cn
                                  s_oracle     the S of the generating permutation, an upper
                                               bound on the optimum (own shape only)

The belt matching here is a numpy reference implementation of PLAN.md item 21 steps 1-5 (no
certificate, no Hungarian). It reproduces the miss statistics quoted there.

Usage:
    py generate_belt_fixtures.py [--cosmochlore target/release/cosmochlore.exe]
        [--timeout 300] [--max-exact-cn 24] [--out tests/belt_reference.csv]

Build the binary first (`cargo build --release`). The exact runs above CN 20 take minutes.
"""

import argparse
import csv
import itertools
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
SHAPES_DIR = ROOT / "src" / "data" / "shapes"
BELT_DIR = ROOT / "tests" / "belt"
SEED = 7
NOISES = (0.02, 0.08)
RADIUS = 2.2
CASES = ((10, "PPR-10"), (11, "JCPPR-11"), (12, "IC-12"), (12, "COC-12"), (20, "DD-20"),
         (24, "TCU-24"), (24, "TOC-24"), (48, "TCOC-48"), (60, "TIC-60"))
EXTRA = ("FeHS.xyz", "La03.xyz")


# ----------------------------------------------------------------------------- input


def builtin_shapes(n_vertices):
    """[(symbol, points)] in YAML order (= the cshm -s index), centre first."""
    text = (SHAPES_DIR / f"shapes_{n_vertices}vertex.yaml").read_text()
    shapes = []
    for block in re.split(r"\n(?=\S)", text):
        head = re.match(r"(\S+):\s*$", block.splitlines()[0])
        if not head:
            continue
        vertices = re.search(r"vertices:\s*\n(.*?)\n\s*center", block, re.S)
        centre = re.search(r"center:\s*\n\s*-\s*\[([^\]]*)\]", block)
        rows = re.findall(r"\[([^\]]*)\]", vertices.group(1))
        points = [[float(x) for x in centre.group(1).split(",")]]
        points += [[float(x) for x in r.split(",")] for r in rows]
        shapes.append((head.group(1), np.array(points)))
    return shapes


def read_xyz(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0])
    return np.array([[float(x) for x in l.split()[1:4]] for l in lines[2:2 + n]])


def write_structures():
    """Distorted structures; returns {file stem: generating permutation (centre first)}."""
    BELT_DIR.mkdir(exist_ok=True)
    rng = np.random.default_rng(SEED)
    perms = {}
    for n, symbol in CASES:
        vertices = dict(builtin_shapes(n))[symbol][1:]
        vertices = vertices / np.sqrt((vertices ** 2).sum(1).mean())
        for noise in NOISES:
            x = vertices * RADIUS + rng.normal(0, noise * RADIUS, vertices.shape)
            order = rng.permutation(len(x))
            x = x[order]
            stem = f"{symbol}_{round(noise * 100):02d}"
            lines = [str(len(x) + 1), f"{symbol} noise {noise} seed {SEED}", "M 0.0 0.0 0.0"]
            lines += [f"X {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}" for p in x]
            (BELT_DIR / f"{stem}.xyz").write_text("\n".join(lines) + "\n")
            perms[stem] = np.concatenate([[0], order + 1])
    return perms


# ----------------------------------------------------------------------------- belt matching


def rigid_s(q, v):
    """CShM for the fixed correspondence q[i] <-> v[i]; improper rotations allowed, as cshm."""
    qc, vc = q - q.mean(0), v - v.mean(0)
    a = np.linalg.svd(vc.T @ qc, compute_uv=False).sum()
    return 100 * (1 - a * a / ((qc ** 2).sum() * (vc ** 2).sum()))


def fibonacci_sphere(n):
    i = np.arange(n) + 0.5
    phi, z = np.pi * (1 + 5 ** 0.5) * i, 1 - 2 * i / n
    r = np.sqrt(1 - z * z)
    return np.column_stack([r * np.cos(phi), r * np.sin(phi), z])


def perpendicular_frame(a):
    t = np.eye(3)[np.argmin(np.abs(a))]
    x = np.cross(a, t)
    x /= np.linalg.norm(x)
    return x, np.cross(a, x)


def height_layers(u, a, tol=1e-3):
    h = u @ a
    order = np.argsort(-h, kind="stable")
    layers, current = [], [order[0]]
    for j in order[1:]:
        if abs(h[j] - h[current[-1]]) < tol:
            current.append(j)
        else:
            layers.append(current)
            current = [j]
    layers.append(current)
    return layers


def reference_belts(vertices):
    """Step 1: the axis with the fewest height layers whose centroids all lie on it. Candidates
    are vertex directions, vertex-pair sums (C2 axes) and vertex-triple plane normals (ring
    centres). Returns the layers (vertex indices, azimuth-sorted from one common zero), their
    azimuths, and the index of the largest layer."""
    u = vertices / np.sqrt((vertices ** 2).sum(1).mean())
    n = len(u)
    i, j = np.triu_indices(n, 1)
    t = np.array(list(itertools.combinations(range(n), 3)))
    raw = np.vstack([u, u[i] + u[j], np.cross(u[t[:, 1]] - u[t[:, 0]], u[t[:, 2]] - u[t[:, 0]])])
    raw = raw[np.linalg.norm(raw, axis=1) > 1e-6]
    raw /= np.linalg.norm(raw, axis=1, keepdims=True)
    raw *= np.where(raw @ np.array([0.3, 0.5, 0.81]) < 0, -1, 1)[:, None]  # +a and -a: one axis
    best = None
    for a in np.unique(np.round(raw, 6), axis=0):
        a = a / np.linalg.norm(a)
        layers = height_layers(u, a)
        if best is not None and len(layers) >= len(best[1]):
            continue
        centroids = [u[l].mean(0) for l in layers]
        if max(np.linalg.norm(c - (c @ a) * a) for c in centroids) < 1e-3:
            best = (a, layers)
    if best is None:
        best = (np.array([0.0, 0.0, 1.0]), height_layers(u, np.array([0.0, 0.0, 1.0])))
    a, layers = best
    x, y = perpendicular_frame(a)
    azimuth = np.arctan2(u @ y, u @ x)
    big = max(range(len(layers)), key=lambda b: len(layers[b]))
    relative = (azimuth - azimuth[layers[big][0]]) % (2 * np.pi)
    layers = [sorted(l, key=lambda k: relative[k]) for l in layers]
    return layers, [relative[l] for l in layers], big


def layer_cost(groups):
    """Step 3, Polynator's cost estimate: layers should be flat, round and on one axis."""
    centroids = [g.mean(0) for g in groups]
    rel = np.vstack([g - c for g, c in zip(groups, centroids)])
    normal = np.linalg.eigh(rel.T @ rel)[1][:, 0]
    cost = ((rel @ normal) ** 2).sum()
    for g, c in zip(groups, centroids):
        if len(g) > 1:
            widths = np.linalg.norm(np.cross(g, normal), axis=1)
            cost += ((widths - widths.mean()) ** 2).sum()
        cost += len(g) * np.linalg.norm(np.cross(c, normal)) ** 2
    return cost


def layer_axis(groups, fallback):
    centroids = np.array([g.mean(0) for g in groups])
    if len(centroids) >= 3:
        c = centroids - centroids.mean(0)
        a = np.linalg.eigh(c.T @ c)[1][:, 2]
    elif len(centroids) == 2 and np.linalg.norm(centroids[0] - centroids[1]) > 1e-6:
        a = centroids[0] - centroids[1]
    else:
        a = fallback
    a = a / np.linalg.norm(a)
    return -a if a @ centroids[0] < 0 else a


def swap_repair(q, v, perm, first):
    """Step 5, Polynator's pairwise swaps: swap two atoms whenever it lowers S, until none does."""
    best, improved = rigid_s(q, v[perm]), True
    while improved:
        improved = False
        for i, j in itertools.combinations(range(first, len(perm)), 2):
            trial = perm.copy()
            trial[i], trial[j] = trial[j], trial[i]
            s = rigid_s(q, v[trial])
            if s < best - 1e-12:
                best, perm, improved = s, trial, True
    return best


def belt_match(q, v, directions=200, keep=8, repair=3):
    """Steps 1-5 with the centre (index 0 on both sides) pinned."""
    layers, ref_azimuths, big = reference_belts(v[1:])
    layers = [[k + 1 for k in l] for l in layers]
    sizes = [len(l) for l in layers]
    qc = q - q.mean(0)
    ligands = np.arange(1, len(q))
    unit = qc[ligands] / np.linalg.norm(qc[ligands], axis=1, keepdims=True)
    principal = np.linalg.eigh(qc[ligands].T @ qc[ligands])[1].T
    seen, candidates = set(), []
    for d in np.vstack([fibonacci_sphere(directions), principal, -principal]):
        order = ligands[np.argsort(-(unit @ d), kind="stable")]
        groups = tuple(tuple(sorted(order[sum(sizes[:b]):sum(sizes[:b + 1])]))
                       for b in range(len(sizes)))
        if groups not in seen:
            seen.add(groups)
            candidates.append((layer_cost([qc[list(g)] for g in groups]), groups, d))
    candidates.sort(key=lambda c: c[0])
    scored = []
    for _, groups, d in candidates[:keep]:
        a = layer_axis([qc[list(g)] for g in groups], d)
        x, y = perpendicular_frame(a)
        azimuth = {i: np.arctan2(qc[i] @ y, qc[i] @ x) for g in groups for i in g}
        for hand in (1, -1):
            for i0 in groups[big]:
                perm = np.zeros(len(q), int)
                for g, l, ref in zip(groups, layers, ref_azimuths):
                    rel = {i: (hand * (azimuth[i] - azimuth[i0])) % (2 * np.pi) for i in g}
                    ordered = sorted(g, key=lambda i: rel[i])
                    m = len(g)
                    shift = min(range(m), key=lambda s: sum(
                        1 - np.cos(rel[ordered[(k + s) % m]] - ref[k]) for k in range(m)))
                    for k, vertex in enumerate(l):
                        perm[ordered[(k + shift) % m]] = vertex
                scored.append((rigid_s(q, v[perm]), perm))
    scored.sort(key=lambda c: c[0])
    s = min(swap_repair(q, v, p.copy(), 1) for _, p in scored[:repair])
    return s, "/".join(map(str, sizes))


# ----------------------------------------------------------------------------- exact cshm


def exact_cshm(binary, path, index, symbol, timeout):
    try:
        out = subprocess.run([str(binary), "cshm", str(path), "-s", str(index)],
                             capture_output=True, text=True, timeout=timeout).stdout
    except subprocess.TimeoutExpired:
        return "", "timeout"
    for line in out.splitlines():
        m = re.match(rf"\s*{re.escape(symbol)}\s.*\s(\d+\.\d+)\s*$", line)
        if m:
            return m.group(1), "ok"
    raise SystemExit(f"error: no {symbol} row in cshm output for {path}:\n{out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    default_bin = ROOT / "target" / "release" / ("cosmochlore.exe" if sys.platform == "win32"
                                                 else "cosmochlore")
    ap.add_argument("--cosmochlore", default=str(default_bin))
    ap.add_argument("--timeout", type=float, default=300.0, help="seconds per exact run")
    ap.add_argument("--max-exact-cn", type=int, default=24)
    ap.add_argument("--out", default=str(ROOT / "tests" / "belt_reference.csv"))
    args = ap.parse_args()
    binary = Path(args.cosmochlore)
    if not binary.is_file():
        sys.exit(f"error: {binary} not found; run `cargo build --release` or pass --cosmochlore")

    perms = write_structures()
    structures = [(BELT_DIR / f"{stem}.xyz", stem) for stem in perms]
    structures += [(ROOT / "tests" / name, None) for name in EXTRA]
    rows = []
    for path, stem in structures:
        q = read_xyz(path)
        cn = len(q) - 1
        noise = f"0.{stem[-2:]}" if stem else ""
        for index, (symbol, v) in enumerate(builtin_shapes(cn)):
            own = stem is not None and stem.rsplit("_", 1)[0] == symbol
            if cn > 24 and not own:
                continue
            s_belt, layers = belt_match(q, v)
            if cn <= args.max_exact_cn:
                s_exact, status = exact_cshm(binary, path, index, symbol, args.timeout)
            else:
                s_exact, status = "", "skipped"
            s_oracle = f"{rigid_s(q, v[perms[stem]]):.6f}" if own else ""
            rows.append({"structure": path.relative_to(ROOT).as_posix(), "cn": cn,
                         "noise": noise, "shape": symbol, "layers": layers,
                         "s_belt": f"{s_belt:.6f}", "s_exact": s_exact, "exact_status": status,
                         "s_oracle": s_oracle})
            print(f"{path.name} {symbol}: belt {s_belt:.3f} exact {s_exact or status}",
                  file=sys.stderr)

    with open(args.out, "w", newline="") as fh:
        fh.write("# Generated by generate_belt_fixtures.py; see PLAN.md item 21. Do not edit.\n")
        fh.write(f"# Exact runs: timeout {args.timeout:g} s, CN <= {args.max_exact_cn}.\n")
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} rows to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
