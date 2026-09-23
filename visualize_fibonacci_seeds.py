"""
Visualise the first N Fibonacci-hemisphere samples that `csom --seeds N` uses as candidate axis
directions (seeding::fibonacci_hemisphere), numbered in the order they are generated.

The structure (default tests/FeHS.xyz) is drawn in the same 3D space as balls and sticks: centred
on its first atom and scaled so the mean distance to the other atoms is 1, so the ligands sit on
the unit sphere next to the seeds. Atoms closer than --bond angstrom are joined by a line.

Usage:
    py visualize_fibonacci_seeds.py                   # 20 seeds (the csom default), FeHS
    py visualize_fibonacci_seeds.py tests/FeHS.xyz --seeds 50
    py visualize_fibonacci_seeds.py --seeds 50 --save seeds.png
"""

import argparse
import math
import re
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

DEFAULT_XYZ = Path(__file__).resolve().parent / "tests" / "FeHS.xyz"
ELEMENT_COLOURS = {"Fe": "#d95f02", "N": "#3050f8", "O": "#e7298a", "C": "#666666", "Cl": "#66a61e"}


def fibonacci_hemisphere(n):
    """seeding::fibonacci_hemisphere: n unit vectors on the upper hemisphere, z from 1 to 1/n."""
    if n == 0:
        return np.array([[0.0, 0.0, 1.0]])
    golden_angle = math.pi * (3 - math.sqrt(5))
    i = np.arange(n)
    z = 1 - i / n
    r = np.sqrt(1 - z * z)
    theta = golden_angle * i
    return np.column_stack([r * np.cos(theta), r * np.sin(theta), z])


def read_xyz(path):
    """Element symbols and coordinates (angstrom) of the atoms in `path`, in file order."""
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    elements, coords = [], []
    for line in lines[2:2 + int(lines[0])]:
        label, *xyz = line.split()
        elements.append(re.match(r"[A-Z][a-z]?", label).group())
        coords.append([float(c) for c in xyz[:3]])
    return elements, np.array(coords)


def draw_molecule(ax, elements, coords, bond_cutoff):
    """Ball-and-stick model, centred on the first atom and scaled to a mean radius of 1."""
    bonds = [(i, j) for i in range(len(coords)) for j in range(i + 1, len(coords))
             if np.linalg.norm(coords[i] - coords[j]) < bond_cutoff]
    centred = coords - coords[0]
    radii = np.linalg.norm(centred[1:], axis=1)
    scaled = centred / radii.mean() if len(radii) else centred

    for i, j in bonds:
        ax.plot(*scaled[[i, j]].T, color="dimgrey", linewidth=2.5)
    for element in dict.fromkeys(elements):
        mask = np.array([e == element for e in elements])
        ax.scatter(*scaled[mask].T, s=250 if element == elements[0] else 150,
                   facecolors=ELEMENT_COLOURS.get(element, "#999999"), edgecolors="k",
                   alpha=0.85, depthshade=False, label=element)


def main():
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument("xyz", nargs="?", default=DEFAULT_XYZ, help="structure (default tests/FeHS.xyz)")
    parser.add_argument("-s", "--seeds", type=int, default=20, help="number of samples (default 20)")
    parser.add_argument("--no-labels", action="store_true", help="do not number the samples")
    parser.add_argument("--bond", type=float, default=2.5, help="bond cutoff in angstrom (default 2.5)")
    parser.add_argument("--save", help="write the figure to this file instead of showing it")
    args = parser.parse_args()

    points = fibonacci_hemisphere(args.seeds)
    order = np.arange(len(points))

    fig = plt.figure(figsize=(13, 6.5))
    fig.suptitle(f"First {len(points)} Fibonacci-hemisphere seeds, {Path(args.xyz).name}")

    # 3D view on a wireframe hemisphere.
    ax3 = fig.add_subplot(1, 2, 1, projection="3d")
    u, v = np.mgrid[0:2 * np.pi:40j, 0:np.pi / 2:12j]
    ax3.plot_wireframe(np.cos(u) * np.sin(v), np.sin(u) * np.sin(v), np.cos(v),
                       color="lightgrey", linewidth=0.5)
    draw_molecule(ax3, *read_xyz(args.xyz), args.bond)
    sc = ax3.scatter(*points.T, c=order, cmap="cool", s=40, depthshade=False)
    if not args.no_labels:
        for k, p in enumerate(points):
            ax3.text(*(p * 1.06), str(k), fontsize=8)
    ax3.set_box_aspect((1, 1, 1))
    ax3.set_xlim(-1.1, 1.1)
    ax3.set_ylim(-1.1, 1.1)
    ax3.set_zlim(-1.1, 1.1)
    ax3.legend(loc="upper left", fontsize=8)
    ax3.set_xlabel("x")
    ax3.set_ylabel("y")
    ax3.set_zlabel("z")

    # Top-down view (looking along -z), unit circle = equator.
    ax2 = fig.add_subplot(1, 2, 2)
    ax2.add_patch(plt.Circle((0, 0), 1, fill=False, color="grey", linewidth=0.8))
    ax2.plot(points[:, 0], points[:, 1], color="lightgrey", linewidth=0.6, zorder=1)
    ax2.scatter(points[:, 0], points[:, 1], c=order, cmap="cool", s=40, zorder=2)
    if not args.no_labels:
        for k, p in enumerate(points):
            ax2.annotate(str(k), p[:2], xytext=(4, 4), textcoords="offset points", fontsize=8)
    ax2.set_aspect("equal")
    ax2.set_xlim(-1.1, 1.1)
    ax2.set_ylim(-1.1, 1.1)
    ax2.set_xlabel("x")
    ax2.set_ylabel("y")
    ax2.set_title("top view (z up)")

    fig.colorbar(sc, ax=[ax3, ax2], label="sample index", shrink=0.8)

    if args.save:
        fig.savefig(args.save, dpi=150, bbox_inches="tight")
        print(f"saved {args.save}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
