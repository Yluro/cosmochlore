"""
Generate tests/dshm_reference.csv: reference values for PLAN.md item 19 (`dshm`, dynamic shape
measures).

This is a small numpy reference implementation of the point-group orbit families described in
PLAN.md item 19, kept so the Rust implementation has independent numbers to match. For every
structure (default: tests/FeHS.xyz, centre = first atom) and every family it reports:

  dof            effective free parameters (rank of the re-centred model matrix B)
  anchors        built-in shapes of the same vertex count that lie inside the family (S < 1e-8)
  s_exhaustive   min over every ligand permutation (centre pinned) of the fixed-permutation fit:
                 the gold value, feasible for CN <= 8
  s_parent_seed  the planned search seeded from the parent shape only: exact rigid matching
                 against the seed, every automorphism of the seed, then re-assignment
  s_all_seeds    the same search seeded from the parent and every anchor (what dshm does)
  path_t         for path families, the optimal path coordinate (radians)
  polynator      the Polynator 1.7.1 model that covers the same shapes, if any

Families come from the src/data/pgs.rs tables and the built-in YAML shapes, rotated into the
parent's frame by `frame` (the directions of the table's z and x axes in parent coordinates).
Every operation must map the parent's points onto themselves or the script stops. Rotations may
be improper, as in cshm. The fixed-centre mode pins the centre to the ligand centroid
(Polynator's convention, `dshm --fix-center`).

The rigid matching here is exhaustive, standing in for cshm's branch-and-bound; both are exact.

Usage:
    py generate_dshm_fixtures.py [--out tests/dshm_reference.csv] [structure.xyz ...]
"""

import argparse
import csv
import itertools
import math
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
PGS_RS = ROOT / "src" / "data" / "pgs.rs"
SHAPES_DIR = ROOT / "src" / "data" / "shapes"

# (parent, group, frame z, frame x, centre mode, Polynator 1.7.1 model or "", note)
FAMILIES = [
    ("OC-6", "Oh", (0, 0, 1), (1, 0, 0), "free", "octahedron[platonic]", ""),
    ("OC-6", "D4h", (0, 0, 1), (1, 0, 0), "free", "tetragonal_bipyramid[4/mmm]", ""),
    ("OC-6", "D3d", (1, 1, 1), (1, 1, -2), "free", "trigonal_antiprism[-3m]", ""),
    ("OC-6", "D2h", (0, 0, 1), (1, 0, 0), "free", "rhombic_bipyramid[mmm]",
     "C2' through vertices"),
    ("OC-6", "D2h", (0, 0, 1), (1, 1, 0), "free", "rectangular_bipyramid[mmm]",
     "C2' between vertices"),
    ("OC-6", "D2d", (0, 0, 1), (1, 1, 0), "free", "",
     "C2' through vertices: collapses onto D4h"),
    ("OC-6", "D2d", (0, 0, 1), (1, 0, 0), "free", "didigonal_scalenohedron[-42m]",
     "C2' between vertices: puckering"),
    ("OC-6", "C4v", (0, 0, 1), (1, 0, 0), "free", "", ""),
    ("OC-6", "C4v", (0, 0, 1), (1, 0, 0), "fixed", "tetragonal_heterobipyramid[4mm]", ""),
    ("OC-6", "C3v", (1, 1, 1), (1, 1, -2), "free", "", ""),
    ("OC-6", "C3v", (1, 1, 1), (1, 1, -2), "fixed", "trigonal_antifrustum[3m]", ""),
    ("OC-6", "D3", (1, 1, 1), (1, 1, -2), "free", "twisted_trigonal_prism[32]", ""),
    ("OC-6", "C2v", (0, 0, 1), (1, 0, 0), "free", "", ""),
    ("TPR-6", "D3h", (0, 0, 1), (0, -1, 0), "free", "trigonal_prism[-6m2]", ""),
    ("TPR-6", "D3", (0, 0, 1), (0, -1, 0), "free", "twisted_trigonal_prism[32]", ""),
    ("TPR-6", "C3v", (0, 0, 1), (1, 0, 0), "free", "",
     "sigma_v through vertices; contains OC-6, so the parent seed alone is not enough"),
]


# ----------------------------------------------------------------------------- input


def read_xyz(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0])
    return np.array([[float(x) for x in l.split()[1:4]] for l in lines[2:2 + n]])


def builtin_shapes(n_vertices):
    """{symbol: points} for one shapes_<n>vertex.yaml, centre first (as cshm uses them)."""
    text = (SHAPES_DIR / f"shapes_{n_vertices}vertex.yaml").read_text()
    shapes = {}
    for block in re.split(r"\n(?=\S)", text):
        head = re.match(r"(\S+):\s*$", block.splitlines()[0])
        if not head:
            continue
        vertices = re.search(r"vertices:\s*\n(.*?)\n\s*center", block, re.S)
        centre = re.search(r"center:\s*\n\s*-\s*\[([^\]]*)\]", block)
        rows = re.findall(r"\[([^\]]*)\]", vertices.group(1))
        points = [[float(x) for x in centre.group(1).split(",")]]
        points += [[float(x) for x in r.split(",")] for r in rows]
        shapes[head.group(1)] = np.array(points)
    return shapes


def pgs_tables():
    """{group name: [3x3 matrices incl. identity]} parsed from src/data/pgs.rs."""
    text = PGS_RS.read_text()
    consts = {"FRAC_1_SQRT_2": 1 / math.sqrt(2)}
    for name, expr in re.findall(r"^const (\w+): f64 = (.+?);", text, re.M):
        consts[name] = eval(expr, {}, consts)   # Rust and Python share 0.866_025 literals
    names = dict(re.findall(r'"(\w+)" => Some\((POINTGROUP_\w+)\)', text))
    tables = {}
    for group, static in names.items():
        body = re.search(rf"pub static {static}: .*?= &\[(.*?)\n\];", text, re.S).group(1)
        mats = [np.eye(3)]
        for rows in re.findall(r"\[\[(.*?)\], \[(.*?)\], \[(.*?)\]\]", body):
            mats.append(np.array([[eval(x, {}, consts) for x in r.split(",")] for r in rows]))
        tables[group] = mats
    return tables


# ----------------------------------------------------------------------------- families


def framed(group_ops, z, x):
    z = np.asarray(z, float) / np.linalg.norm(z)
    x = np.asarray(x, float)
    x = x - (x @ z) * z
    x /= np.linalg.norm(x)
    frame = np.column_stack([x, np.cross(z, x), z])
    return [frame @ g @ frame.T for g in group_ops]


def index_of(points, p, tol=1e-6):
    hits = np.flatnonzero(np.linalg.norm(points - p, axis=1) < tol)
    return int(hits[0]) if len(hits) else None


def orbit_basis(parent, ops):
    """B (n x 3 x k) with v(theta) = B @ theta for every point of the family, and theta of the
    parent. Each orbit of the group contributes the fixed subspace of its representative's
    stabiliser."""
    n = len(parent)
    for g in ops:
        if any(index_of(parent, g @ p) is None for p in parent):
            raise SystemExit("frame error: an operation does not map the parent onto itself")
    seen, blocks = [False] * n, []
    for i in range(n):
        if seen[i]:
            continue
        stab = [g for g in ops if np.linalg.norm(g @ parent[i] - parent[i]) < 1e-7]
        u, s, _ = np.linalg.svd(sum(stab) / len(stab))
        basis = u[:, s > 0.5]
        members = {}
        for g in ops:
            members.setdefault(index_of(parent, g @ parent[i]), g @ basis)
        for j in members:
            seen[j] = True
        blocks.append(members)
    k = sum(next(iter(b.values())).shape[1] for b in blocks)
    B, col = np.zeros((n, 3, k)), 0
    for b in blocks:
        d = next(iter(b.values())).shape[1]
        for j, bj in b.items():
            B[j, :, col:col + d] = bj
        col += d
    theta = np.linalg.lstsq(B.reshape(-1, k), parent.reshape(-1), rcond=None)[0]
    assert np.allclose(B @ theta, parent, atol=1e-8), "parent is not a member of its family"
    return B, theta


def pin_centre_to_ligand_centroid(B):
    """--fix-center: centre point = centroid of the ligand points (still linear in theta)."""
    B = B.copy()
    B[0] = B[1:].mean(0)
    keep = np.linalg.norm(B.reshape(-1, B.shape[2]), axis=0) > 1e-12
    return B[:, :, keep]


def effective_dof(B):
    A = (B - B.mean(0)).reshape(-1, B.shape[2])
    s = np.linalg.svd(A, compute_uv=False)
    return int((s > 1e-9 * s[0]).sum())


# ----------------------------------------------------------------------------- measures


def rotation(q, v):
    """Orthogonal R (improper allowed, as cshm) minimising |q - v R^T|; both centred."""
    u, _, vt = np.linalg.svd(v.T @ q)
    return (u @ vt).T


def rigid_s(q, v):
    """CShM for the fixed correspondence q[i] <-> v[i]."""
    qc, vc = q - q.mean(0), v - v.mean(0)
    a = np.linalg.svd(vc.T @ qc, compute_uv=False).sum()
    return 100 * (1 - a * a / ((qc ** 2).sum() * (vc ** 2).sum()))


def fit(q, B, theta, perm, tol=1e-13, max_iter=500):
    """Fixed-permutation dynamic fit: alternate rotation and the linear least-squares solve."""
    qc = q - q.mean(0)
    qq = (qc ** 2).sum()
    Bp = B[perm] - B[perm].mean(0)            # re-centring keeps the translation optimal
    A = Bp.reshape(-1, Bp.shape[2])
    prev = math.inf
    for _ in range(max_iter):
        r = rotation(qc, Bp @ theta)
        theta = np.linalg.lstsq(A, (qc @ r).reshape(-1), rcond=None)[0]
        s = 100 * ((qc - (Bp @ theta) @ r.T) ** 2).sum() / qq
        if prev - s < tol:
            break
        prev = s
    return s, theta


def ligand_perms(n):
    for p in itertools.permutations(range(1, n)):
        yield [0, *p]


def best_rigid_perm(q, v):
    return min(ligand_perms(len(q)), key=lambda p: rigid_s(q, v[p]))


def automorphisms(v):
    return [p for p in ligand_perms(len(v)) if rigid_s(v, v[p]) < 1e-8]


def exhaustive(q, B, theta0):
    best = (math.inf, None)
    for perm in ligand_perms(len(q)):
        s, theta = fit(q, B, theta0.copy(), perm)
        if s < best[0]:
            best = (s, theta)
    return best


def planned_search(q, B, seeds):
    """dshm's search: per seed, exact rigid matching, all seed automorphisms, re-assignment."""
    best = (math.inf, None, None)
    for theta_seed in seeds:
        v = B @ theta_seed
        p0 = best_rigid_perm(q, v)
        for a in automorphisms(v):
            perm = [a[j] for j in p0]
            s, theta = fit(q, B, theta_seed.copy(), perm)
            if s < best[0]:
                best = (s, theta, perm)
    s, theta, perm = best
    while True:
        new = best_rigid_perm(q, B @ theta)
        if new == perm:
            return s
        s2, theta2 = fit(q, B, theta, new)
        if s2 >= s - 1e-12:
            return s
        s, theta, perm = s2, theta2, new


def bailar(t):
    """Polynator's bailar_twist[dynamic]: unit sphere, relative twist t; t = 0 is the trigonal
    prism, t = pi/3 the octahedron."""
    c = math.cos(t)
    h, w = math.sqrt((0.25 + 0.5 * c) / (1.25 + 0.5 * c)), 1 / math.sqrt(1.25 + 0.5 * c)
    pts = [[0.0, 0.0, 0.0]]
    for sign in (1, -1):
        for k in range(3):
            a = sign * t / 2 + 2 * math.pi * k / 3
            pts.append([w * math.cos(a), w * math.sin(a), sign * h])
    return np.array(pts)


def path_fit(q, path, lo, hi, tol=1e-7):
    """Golden-section search over the path coordinate, exact rigid CShM at every point."""
    f = lambda t: rigid_s(q, path(t)[best_rigid_perm(q, path(t))])
    g = (math.sqrt(5) - 1) / 2
    a, b = lo, hi
    while b - a > tol:
        c, d = b - g * (b - a), a + g * (b - a)
        if f(c) < f(d):
            b = d
        else:
            a = c
    t = (a + b) / 2
    return t, f(t)


# ----------------------------------------------------------------------------- main


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("structures", nargs="*", default=[str(ROOT / "tests" / "FeHS.xyz")])
    ap.add_argument("--out", default=str(ROOT / "tests" / "dshm_reference.csv"))
    args = ap.parse_args()

    tables = pgs_tables()
    rows = []
    for structure in args.structures:
        q = read_xyz(structure)
        n = len(q) - 1
        shapes = builtin_shapes(n)
        name = Path(structure).name
        for parent, group, z, x, centre, polynator, note in FAMILIES:
            if parent not in shapes:
                continue
            B, theta_parent = orbit_basis(shapes[parent], framed(tables[group], z, x))
            if centre == "fixed":
                B = pin_centre_to_ligand_centroid(B)
                theta_parent = np.linalg.lstsq(B.reshape(-1, B.shape[2]),
                                               shapes[parent].reshape(-1), rcond=None)[0]
            anchors, seeds = [], [theta_parent]
            for symbol, pts in shapes.items():
                s, theta = exhaustive(pts, B, theta_parent)
                if s < 1e-8:
                    anchors.append(symbol)
                    if symbol != parent:
                        seeds.append(theta)
            s_all, _ = exhaustive(q, B, theta_parent)
            rows.append({
                "structure": name, "family": f"{parent}/{group}", "parent": parent,
                "group": group, "frame_z": " ".join(map(str, z)),
                "frame_x": " ".join(map(str, x)), "centre": centre, "dof": effective_dof(B),
                "anchors": " ".join(anchors), "s_exhaustive": f"{s_all:.6f}",
                "s_parent_seed": f"{planned_search(q, B, [theta_parent]):.6f}",
                "s_all_seeds": f"{planned_search(q, B, seeds):.6f}", "path_t": "",
                "polynator": polynator, "note": note,
            })
            print(f"{name} {parent}/{group} {centre}: {rows[-1]['s_exhaustive']}", file=sys.stderr)
        if n == 6:
            t, s = path_fit(q, bailar, 0.0, math.pi / 3)
            rows.append({
                "structure": name, "family": "Bailar(TPR-6~OC-6)", "parent": "OC-6", "group": "D3",
                "frame_z": "", "frame_x": "", "centre": "free", "dof": 2,
                "anchors": "TPR-6 OC-6", "s_exhaustive": f"{s:.6f}", "s_parent_seed": "",
                "s_all_seeds": "", "path_t": f"{t:.6f}", "polynator": "bailar_twist[dynamic]",
                "note": "Polynator moves the ligands but not the centre; its S differs",
            })

    with open(args.out, "w", newline="") as fh:
        fh.write("# Generated by generate_dshm_fixtures.py; see PLAN.md item 19. Do not edit.\n")
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} rows to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
