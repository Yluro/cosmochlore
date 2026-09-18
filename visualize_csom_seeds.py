"""
Visualise how the csom symmetry-axis search seeds one point group versus another.

Reproduces src/csom/seeding.rs in numpy, reading the point-group tables straight from
src/data/pgs.rs, and draws one panel per group. Every group shares the same `--seeds`
candidate axes (black dots: the Fibonacci hemisphere, each carried onto z); what changes
from group to group is the in-plane spin about that axis:

  * a group whose operations are all about z (C2, Cs, S4, C6h, ...) has no in-plane
    period, so every axis is one seed and the spin is whatever the alignment produced;
  * a group with vertical planes or perpendicular C2 axes repeats after a spin of
    `period` degrees (90 for C2v and Oh, 60 for D3h, 30 for D6h, 72 for Ih, ...), and each
    axis is spun through that period in ceil(period / lattice spacing) steps.

The coloured dots are, for each seed, the direction of the structure that the seed sends to
+x: the k = 0 spin is the plain alignment, k = 1, 2, ... the extra spins, joined per axis by an
arc along the great circle perpendicular to it. With --vectors the panels instead show the
seeds as rotation vectors, the parameters Nelder-Mead refines: the k = 0 seeds form a disk in
the v_z = 0 plane and each extra spin lifts a copy of it along v_z.

Usage:
    py visualize_csom_seeds.py                                # C2 D6h D3h D4h C2v Oh
    py visualize_csom_seeds.py --groups Cs D2d Td Ih --seeds 40
    py visualize_csom_seeds.py --vectors --save seeds.png
"""

import argparse
import math
import re
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

PGS_RS = Path(__file__).resolve().parent / "src" / "data" / "pgs.rs"

# The named constants pgs.rs uses inside its matrices.
CONSTANTS = {
    "SQRT_3_DIV_2": math.sqrt(3) / 2,
    "COS_72": math.cos(math.radians(72)),
    "SIN_72": math.sin(math.radians(72)),
    "COS_144": -math.cos(math.radians(144)),  # pgs.rs stores the magnitude, 0.809
    "SIN_144": math.sin(math.radians(144)),
    "FRAC_1_SQRT_2": 1 / math.sqrt(2),
}

Z = np.array([0.0, 0.0, 1.0])
X = np.array([1.0, 0.0, 0.0])


# --- Point-group tables, read from pgs.rs ----------------------------------------------------

def read_point_group(name):
    """The operation matrices of `name` as pgs.rs::get_pointgroup returns them (E excluded)."""
    source = PGS_RS.read_text(encoding="utf-8")
    static = f"POINTGROUP_{name.upper()}"
    match = re.search(rf"pub static {static}: [^=]+= &\[(.*?)^\];", source, re.S | re.M)
    if match is None:
        raise SystemExit(f"{name}: no {static} table in {PGS_RS}")

    def value(token):
        token = token.strip()
        sign = -1.0 if token.startswith("-") else 1.0
        token = token.lstrip("-")
        return sign * (CONSTANTS[token] if token in CONSTANTS else float(token))

    matrices = []
    for row_block in re.finditer(r'\("[^"]*",\s*\[(\[.*?\]),\s*(\[.*?\]),\s*(\[.*?\])\]\)', match.group(1), re.S):
        rows = [[value(t) for t in row.strip("[]").split(",")] for row in row_block.groups()]
        matrices.append(np.array(rows))
    return matrices


# --- The same maths as seeding.rs -----------------------------------------------------------

def rodrigues(v):
    """geometry::rotation_matrix_from_vector."""
    angle = np.linalg.norm(v)
    if angle < 1e-10:
        return np.eye(3)
    k = v / angle
    kx = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + math.sin(angle) * kx + (1 - math.cos(angle)) * kx @ kx


def log_map(r):
    """nalgebra's Rotation3::scaled_axis."""
    angle = math.acos(max(-1.0, min(1.0, (np.trace(r) - 1) / 2)))
    if angle < 1e-10:
        return np.zeros(3)
    axis = np.array([r[2, 1] - r[1, 2], r[0, 2] - r[2, 0], r[1, 0] - r[0, 1]])
    return axis / (2 * math.sin(angle)) * angle


def rotation_between(a, b):
    """nalgebra's Rotation3::rotation_between."""
    axis = np.cross(a, b)
    sin = np.linalg.norm(axis)
    if sin < 1e-12:
        return np.eye(3)
    return rodrigues(axis / sin * math.atan2(sin, a @ b))


def spin_about_z(angle):
    return rodrigues(np.array([0.0, 0.0, angle]))


def fibonacci_hemisphere(n):
    """seeding::fibonacci_hemisphere."""
    golden_angle = math.pi * (3 - math.sqrt(5))
    points = []
    for i in range(max(n, 1)):
        z = 1 - i / n if n else 1.0
        r = math.sqrt(1 - z * z)
        theta = golden_angle * i
        points.append([r * math.cos(theta), r * math.sin(theta), z])
    return np.array(points)


def in_plane_period(matrices):
    """seeding::in_plane_period: None for a group invariant under any spin about z, else the
    smallest spin that maps the operation set onto itself."""
    eps = 1e-9

    def axial(m):
        return (abs(m[0, 2]) < eps and abs(m[1, 2]) < eps and abs(m[2, 0]) < eps and abs(m[2, 1]) < eps
                and abs(m[0, 0] - m[1, 1]) < eps and abs(m[0, 1] + m[1, 0]) < eps)

    if all(axial(m) for m in matrices):
        return None

    n = 1
    for m in matrices:
        if axial(m) and abs(m[2, 2] - 1) < eps:
            angle = abs(math.atan2(m[1, 0], m[0, 0]))
            if angle > eps:
                n = max(n, round(2 * math.pi / angle))

    def maps_group_onto_itself(angle):
        spin = spin_about_z(angle)
        return all(any(np.linalg.norm(spin @ m @ spin.T - other) < 1e-2 for other in matrices)
                   for m in matrices)

    for q in range(2 * n, 0, -1):
        if maps_group_onto_itself(2 * math.pi / q):
            return 2 * math.pi / q
    raise AssertionError("q = 1 always works")


def seed_rotation_vectors(n_directions, matrices):
    """seeding::seed_rotation_vectors. Returns (seeds as (axis index, k, rotation vector),
    period, spins per axis, step)."""
    spacing = math.sqrt(2 * math.pi / max(n_directions, 1))
    period = in_plane_period(matrices)
    if period is None:
        n_in_plane, step = 1, 0.0
    else:
        n_in_plane = max(1, math.ceil(period / spacing))
        step = period / n_in_plane

    seeds = []
    for i, d in enumerate(fibonacci_hemisphere(n_directions)):
        align = rotation_between(d, Z)
        for k in range(n_in_plane):
            seeds.append((i, k, log_map(spin_about_z(k * step) @ align)))
    return seeds, period, n_in_plane, step


# --- Plotting -------------------------------------------------------------------------------

def draw_sphere(ax, alpha=0.08):
    u = np.linspace(0, 2 * np.pi, 36)
    v = np.linspace(0, np.pi, 18)
    ax.plot_wireframe(np.outer(np.cos(u), np.sin(v)), np.outer(np.sin(u), np.sin(v)),
                      np.outer(np.ones_like(u), np.cos(v)), color="grey", linewidth=0.4, alpha=alpha)


def draw_arrow(ax, v, **kwargs):
    ax.quiver(0, 0, 0, v[0], v[1], v[2], arrow_length_ratio=0.12, **kwargs)


def style_axes(ax, labels=("x", "y", "z"), lim=1.05, elev=28, azim=-60):
    ax.view_init(elev=elev, azim=azim)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-lim, lim)
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])
    ax.set_zlabel(labels[2])
    ax.set_box_aspect((1, 1, 1))
    ax.tick_params(labelsize=7)


def panel_title(group, n_ops, period, n_in_plane, step, n_seeds):
    if period is None:
        return f"{group} ({n_ops} ops): all about z, no in-plane period\n1 seed per axis -> {n_seeds} seeds"
    spins = "1 spin" if n_in_plane == 1 else f"{n_in_plane} spins"
    return (f"{group} ({n_ops} ops): period {math.degrees(period):g} deg\n"
            f"{spins} per axis, step {math.degrees(step):g} deg -> {n_seeds} seeds")


def draw_directions(ax, directions, seeds, period, n_in_plane, step, colours):
    """The direction of the structure each seed sends to +x, around its candidate axis."""
    draw_sphere(ax)
    draw_arrow(ax, Z, color="k", linewidth=0.8)
    draw_arrow(ax, X, color="k", linewidth=0.8)
    ax.scatter(*directions.T, c="k", s=12)

    to_x = {(i, k): rodrigues(v).T @ X for i, k, v in seeds}
    for k in range(n_in_plane):
        pts = np.array([to_x[i, k] for i in range(len(directions))])
        ax.scatter(*pts.T, c=["silver"] if period is None else [colours[k]], s=22)

    if n_in_plane > 1:
        t = np.linspace(0, 1, 12)
        for i, d in enumerate(directions):
            # From the k = 0 image, spinning about the axis d retraces the other spins.
            arc = np.array([rodrigues(d * (-s * (n_in_plane - 1) * step)) @ to_x[i, 0] for s in t])
            ax.plot(*arc.T, color="grey", linewidth=0.6)
    style_axes(ax)


def draw_vectors(ax, seeds, period, n_in_plane, colours):
    """The seeds as rotation vectors."""
    draw_sphere(ax, alpha=0.05)
    vectors = np.array([v for _, _, v in seeds])
    ks = np.array([k for _, k, _ in seeds])
    for k in range(n_in_plane):
        ax.scatter(*vectors[ks == k].T, c=["silver"] if period is None else [colours[k]], s=22)
    ax.scatter([0], [0], [0], c="red", marker="*", s=90)
    style_axes(ax, labels=("v_x (rad)", "v_y (rad)", "v_z (rad)"), lim=1.6, elev=8)


def main():
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument("--groups", nargs="+", default=["C2", "D6h", "D3h", "D4h", "C2v", "Oh"],
                        help="point groups to draw, one panel each (default: C2 D6h D3h D4h C2v Oh)")
    parser.add_argument("--seeds", type=int, default=20, help="number of candidate axes (csom --seeds)")
    parser.add_argument("--vectors", action="store_true",
                        help="show the seeds as rotation vectors instead of as directions on the sphere")
    parser.add_argument("--save", metavar="FILE", help="write the figure to FILE instead of showing it")
    args = parser.parse_args()

    n = args.seeds
    directions = fibonacci_hemisphere(n)
    spacing = math.sqrt(2 * math.pi / max(n, 1))
    print(f"--seeds {n}: {len(directions)} candidate axes, lattice spacing "
          f"sqrt(2pi/{n}) = {math.degrees(spacing):.1f} deg (shared by every group)\n")
    print(f"{'group':<6} {'ops':>4} {'period':>8} {'spins':>6} {'step':>7} {'seeds':>6}")

    per_group = []
    for group in args.groups:
        matrices = read_point_group(group)
        seeds, period, n_in_plane, step = seed_rotation_vectors(n, matrices)
        per_group.append((group, matrices, seeds, period, n_in_plane, step))
        print(f"{group:<6} {len(matrices):>4} {'-' if period is None else f'{math.degrees(period):.1f}':>8} "
              f"{n_in_plane:>6} {'-' if period is None else f'{math.degrees(step):.1f}':>7} {len(seeds):>6}")

    max_spins = max(g[4] for g in per_group)
    colours = plt.cm.viridis(np.linspace(0.15, 0.9, max(max_spins, 2)))

    cols = min(3, len(per_group))
    rows = math.ceil(len(per_group) / cols)
    fig = plt.figure(figsize=(5 * cols, 5.4 * rows + 1.2))
    if args.vectors:
        view = (f"seeds as rotation vectors: the {len(directions)} k = 0 alignments form a disk in the "
                "v_z = 0 plane, each extra spin lifts a copy of it")
    else:
        view = f"{len(directions)} candidate axes (black) shared by every group; coloured: direction sent to +x per seed"
    fig.suptitle(f"csom seeds per point group, --seeds {n}, spun in-plane as each group requires\n[{view}]",
                 fontsize=11)

    for idx, (group, matrices, seeds, period, n_in_plane, step) in enumerate(per_group):
        ax = fig.add_subplot(rows, cols, idx + 1, projection="3d")
        if args.vectors:
            draw_vectors(ax, seeds, period, n_in_plane, colours)
        else:
            draw_directions(ax, directions, seeds, period, n_in_plane, step, colours)
        ax.set_title(panel_title(group, len(matrices), period, n_in_plane, step, len(seeds)), fontsize=9)

    # One legend for the whole figure.
    handles = [plt.Line2D([], [], marker="o", color="k", linestyle="", markersize=5,
                          label="candidate axis (sent to z)" if not args.vectors else "identity: input orientation as-is"),
               plt.Line2D([], [], marker="o", color="silver", linestyle="", markersize=5,
                          label="no in-plane period: spin left as aligned")]
    if args.vectors:
        handles[0] = plt.Line2D([], [], marker="*", color="red", linestyle="", markersize=9,
                                label="identity: input orientation as-is")
    for k in range(max_spins):
        handles.append(plt.Line2D([], [], marker="o", color=colours[k], linestyle="", markersize=5,
                                  label=f"k = {k}: {'alignment only' if k == 0 else f'extra spin {k} x step'}"))
    fig.legend(handles=handles, loc="lower center", ncol=len(handles), fontsize=8, frameon=False)

    fig.subplots_adjust(left=0.02, right=0.98, top=0.88, bottom=0.07, wspace=0.05, hspace=0.3)
    if args.save:
        fig.savefig(args.save, dpi=130)
        print(f"\nwrote {args.save}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
