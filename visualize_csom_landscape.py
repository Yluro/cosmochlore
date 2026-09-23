"""
Visualise the csom deviation landscape of a structure over the direction of its symmetry axis.

For every direction n(theta, phi) on the unit sphere -- theta the polar angle from +z, phi the
azimuth from +x, the usual spherical coordinates, in degrees -- the structure is scored against
itself rotated by `--angle` degrees about n (default 180, a C2) through the same steps as
src/csom: centre on the first atom (or on the centroid with --nc), normalise, find the best
atom-to-atom assignment within each element (Hungarian, centre pinned) and take
deviation::sds_dev of the pair. The result is drawn twice: as a surface over the (theta, phi)
rectangle, and painted onto the unit sphere in the structure's own frame with the ligand
directions overlaid, so each valley can be matched to the bonds it sits on. The local minima of
that surface are the structure's approximate C_n axes: the targets a csom seed has to start
within reach of. With --seeds N the N Fibonacci-hemisphere candidate axes of `csom --seeds N`
are overlaid on the sphere.

With --pg the single rotation is replaced by the point group's whole operation set (read from
src/data/pgs.rs, as visualize_csom_seeds.py does): n is carried onto z by the smallest rotation
that does so, spun by --spin degrees about z, every operation of the group is applied and the
average sds_dev is taken -- exactly deviation::point_group_dev at the orientation a k = 0 seed
starts from (see seeding.rs). --angle A is the same thing for the one-operation group {C(z, A)}.

A rotation by A about n and by A about -n score identically (one is the inverse of the other,
and the assignment is a bijection), so with --angle the sphere is antipodally symmetric and
every axis shows up twice, at (theta, phi) and (180 - theta, phi + 180).

Usage:
    py visualize_csom_landscape.py                                # tests/FeHS.xyz, C2 (180 deg)
    py visualize_csom_landscape.py tests/FeHS.xyz --angle 90      # the C4 axes instead
    py visualize_csom_landscape.py tests/FeHS.xyz --pg Oh --step 3
    py visualize_csom_landscape.py --seeds 20 --save landscape.png
"""

import argparse
import math
import time
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import art3d
from scipy.ndimage import minimum_filter
from scipy.optimize import linear_sum_assignment

from visualize_csom_seeds import fibonacci_hemisphere, read_point_group, rodrigues

Z = np.array([0.0, 0.0, 1.0])
DEFAULT_XYZ = Path(__file__).resolve().parent / "tests" / "FeHS.xyz"

# Fixed element colours for the ligand markers (fallback: tab10 in order of appearance).
ELEMENT_COLOURS = {"Fe": "#d95f02", "N": "#1b9e77", "O": "#e7298a", "C": "#666666", "Cl": "#66a61e"}


# --- The structure, prepared as csom::prepare does ------------------------------------------

def read_xyz(path):
    """Labels and coordinates of the atoms in `path`, in file order (xyz::parse_xyz)."""
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    labels, coords = [], []
    for line in lines[2:]:
        parts = line.split()
        if not parts:
            continue
        if len(parts) != 4:
            raise SystemExit(f"{path}: expected 'label x y z', got {line!r}")
        labels.append(parts[0])
        coords.append([float(v) for v in parts[1:]])
    if len(labels) != int(lines[0].strip()):
        raise SystemExit(f"{path}: header says {lines[0].strip()} atoms, found {len(labels)}")
    return labels, np.array(coords)


def strip_label(label):
    """prepare::strip_label: the element symbol, two letters when the second is lowercase."""
    return label[:2] if len(label) > 1 and label[1].islower() else label[:1]


def prepare(coords, has_centre):
    """prepare::prepare_csom_structure with --mode auto: centre on the first atom when the
    structure has one, else on the centroid, then normalise so that sum |p|^2 = n."""
    points = coords - (coords[0] if has_centre else coords.mean(axis=0))
    return points * math.sqrt(len(points) / np.sum(points ** 2))


def assignment_groups(labels, has_centre, ignore_labels):
    """The index groups the assignment is restricted to (assignment::group_by_label, or every
    atom at once with --ignore), minus the centre when it is pinned to itself."""
    groups = {}
    for i, label in enumerate(labels):
        if has_centre and i == 0:
            continue
        groups.setdefault("" if ignore_labels else strip_label(label), []).append(i)
    return [np.array(g) for g in groups.values()]


# --- The same maths as csom::deviation ------------------------------------------------------

def sds_dev(reference, problem):
    """deviation::sds_dev: 100 * sum |q - p|^2 / sum |q - centroid(q)|^2."""
    centred = problem - problem.mean(axis=0)
    return 100.0 * np.sum((problem - reference) ** 2) / np.sum(centred ** 2)


def operation_deviation(rotated, image, groups, pinned):
    """deviation::operation_deviation_scalar: sds_dev of `rotated` against its image under the
    best assignment within each group (assignment::best_permutation_multiple_atoms), the centre
    paired with itself when pinned."""
    a_rows = [rotated[:1]] if pinned else []
    b_rows = [image[:1]] if pinned else []
    for idx in groups:
        # Minimising the squared distance <=> minimising -dot, as sq_dist_cost does.
        rows, cols = linear_sum_assignment(-rotated[idx] @ image[idx].T)
        a_rows.append(rotated[idx][rows])
        b_rows.append(image[idx][cols])
    return sds_dev(np.concatenate(a_rows), np.concatenate(b_rows))


def direction(theta, phi):
    """Unit vector at polar angle `theta` and azimuth `phi` (radians)."""
    return np.array([math.sin(theta) * math.cos(phi), math.sin(theta) * math.sin(phi), math.cos(theta)])


def align_to_z(d):
    """The smallest rotation carrying unit vector `d` onto +z (nalgebra's rotation_between, as
    seeding.rs uses it), with the antiparallel case d = -z, which rotation_between leaves
    undefined, taken as a half turn about x."""
    axis = np.cross(d, Z)
    s = np.linalg.norm(axis)
    c = d @ Z
    if s < 1e-12:
        return np.eye(3) if c > 0 else np.diag([1.0, -1.0, -1.0])
    return rodrigues(axis / s * math.atan2(s, c))


def landscape(points, groups, pinned, ops, spin, thetas, phis):
    """deviation::point_group_dev over every direction n(theta, phi): the structure is rotated
    so that n lands on z (then spun by `spin` about z), each of the `ops` matrices is applied,
    and the per-operation deviations are averaged."""
    spin_matrix = rodrigues(np.array([0.0, 0.0, spin]))
    values = np.empty((len(thetas), len(phis)))
    for i, theta in enumerate(thetas):
        for j, phi in enumerate(phis):
            rotated = points @ (spin_matrix @ align_to_z(direction(theta, phi))).T
            values[i, j] = np.mean([operation_deviation(rotated, rotated @ op.T, groups, pinned)
                                    for op in ops])
    return values


def local_minima(values, thetas, phis, merge_deg):
    """Grid cells no higher than their eight neighbours (phi periodic, theta clamped at the
    poles), keeping only the lowest of any set that point within `merge_deg` of each other,
    which also folds each all-equal pole row into one entry. Returns (i, j) sorted by value."""
    is_min = values == minimum_filter(values, size=3, mode=("nearest", "wrap"))
    candidates = sorted(zip(*np.nonzero(is_min)), key=lambda ij: values[ij])
    kept, kept_dirs = [], []
    for i, j in candidates:
        d = direction(thetas[i], phis[j])
        if all(math.degrees(math.acos(min(1.0, d @ k))) > merge_deg for k in kept_dirs):
            kept.append((i, j))
            kept_dirs.append(d)
    return kept


# --- Plotting -------------------------------------------------------------------------------

def element_colour(element, fallback):
    return ELEMENT_COLOURS.get(element) or fallback.setdefault(element, plt.cm.tab10(len(fallback) % 10))


def draw_surface(ax, thetas_deg, phis_deg, values, minima, norm):
    """The landscape as a surface over the (theta, phi) rectangle, closed at phi = 360."""
    closed = np.concatenate([values, values[:, :1]], axis=1)
    th, ph = np.meshgrid(thetas_deg, np.append(phis_deg, 360.0), indexing="ij")
    ax.computed_zorder = False  # the marker collection must not hide behind the surface
    ax.plot_surface(th, ph, closed, cmap="viridis", norm=norm, rcount=closed.shape[0],
                    ccount=closed.shape[1], linewidth=0, antialiased=False, zorder=1)
    for rank, (i, j) in enumerate(minima, 1):
        ax.scatter(thetas_deg[i], phis_deg[j], values[i, j], color="red", s=30, zorder=2,
                   depthshade=False)
        ax.text(thetas_deg[i], phis_deg[j], values[i, j], f" {rank}", fontsize=8, zorder=3)
    ax.set_xlabel("theta (deg, from +z)")
    ax.set_ylabel("phi (deg, from +x)")
    ax.set_zlabel("sds_dev")
    ax.set_xlim(0, 180)
    ax.set_ylim(0, 360)
    ax.set_xticks(range(0, 181, 45))
    ax.set_yticks(range(0, 361, 90))
    ax.view_init(elev=35, azim=-50)
    ax.tick_params(labelsize=8)


def view_direction(ax):
    """Unit vector from the origin towards the eye, as Axes3D.get_proj places it."""
    elev, azim = math.radians(ax.elev), math.radians(ax.azim)
    return np.array([math.cos(elev) * math.cos(azim), math.cos(elev) * math.sin(azim), math.sin(elev)])


class NearSideText(art3d.Text3D):
    """A label on the sphere that is only drawn while its position faces the viewer, so it
    disappears together with the marker it names when that rolls round the back."""

    def draw(self, renderer):
        if np.dot(self.get_position_3d(), view_direction(self.axes)) > 0:
            super().draw(renderer)


def sphere_marker(ax, d, label=None, label_colour="k", **kwargs):
    """One marker just outside the unit sphere at direction `d`, as its own artist so the
    depth sort can hide it behind the sphere, with an optional near-side-only label."""
    ax.scatter(*(1.04 * d), depthshade=False, **kwargs)
    if label is not None:
        # Anchored at the marker and offset on screen (a radial offset vanishes in projection
        # for markers facing the viewer). zdir=None keeps the text upright; the high zorder
        # keeps it above the depth-sorted collections (it is only drawn while facing the viewer
        # anyway).
        ax.add_artist(NearSideText(*(1.04 * d), f" {label}", zdir=None, fontsize=7, color=label_colour, ha="left",
                                   va="bottom", clip_on=False, transform=ax.transData, zorder=1000))


def draw_sphere(ax, thetas, phis, values, minima, norm, ligands, seeds):
    """The landscape painted onto the unit sphere in the structure's frame. The sphere is
    depth-sorted at its centre and every marker sits just outside it, so matplotlib draws the
    markers facing the viewer on top and hides the ones on the far side behind the surface."""
    th, ph = np.meshgrid(thetas, np.append(phis, 2 * math.pi), indexing="ij")
    closed = np.concatenate([values, values[:, :1]], axis=1)
    sphere = ax.plot_surface(np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th),
                             facecolors=cm.viridis(norm(closed)), rcount=closed.shape[0],
                             ccount=closed.shape[1], shade=False, linewidth=0, antialiased=False)
    sphere.set_sort_zpos(0.0)

    fallback = {}
    for label, d in ligands:
        sphere_marker(ax, d, label, color=element_colour(strip_label(label), fallback), marker="^", s=45,
                      edgecolors="k", linewidths=0.5)
    for d in seeds:
        sphere_marker(ax, d, color="k", s=10)
    for rank, (i, j) in enumerate(minima, 1):
        sphere_marker(ax, direction(thetas[i], phis[j]), str(rank), label_colour="red", color="red", s=30)

    lim = 1.15
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-lim, lim)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=28, azim=-60)
    ax.tick_params(labelsize=7)
    return fallback


def main():
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument("xyz", nargs="?", default=str(DEFAULT_XYZ), help="structure to analyse (default: tests/FeHS.xyz)")
    what = parser.add_mutually_exclusive_group()
    what.add_argument("--angle", type=float, default=180.0,
                      help="rotate the structure by this many degrees about n(theta, phi) (default: 180, a C2)")
    what.add_argument("--pg", metavar="GROUP",
                      help="instead, average over every operation of this point group with n as its principal axis")
    parser.add_argument("--spin", type=float, default=0.0,
                        help="with --pg: in-plane spin about the axis, in degrees (default: 0, the k = 0 seed)")
    parser.add_argument("--nc", action="store_true", help="the structure has no central atom (csom --nc)")
    parser.add_argument("--ignore", action="store_true", help="ignore atom labels in the assignment (csom --ignore)")
    parser.add_argument("--step", type=float, default=2.0, help="grid spacing in degrees for both angles (default: 2)")
    parser.add_argument("--seeds", type=int, default=0, metavar="N",
                        help="overlay the N Fibonacci-hemisphere candidate axes of csom --seeds N on the sphere")
    parser.add_argument("--save", metavar="FILE", help="write the figure to FILE instead of showing it")
    args = parser.parse_args()

    labels, coords = read_xyz(args.xyz)
    has_centre = not args.nc
    points = prepare(coords, has_centre)
    groups = assignment_groups(labels, has_centre, args.ignore)

    if args.pg:
        ops = read_point_group(args.pg)
        what = f"{args.pg} ({len(ops)} ops, spin {args.spin:g} deg) with n(theta, phi) as principal axis"
    else:
        ops = [rodrigues(np.array([0.0, 0.0, math.radians(args.angle)]))]
        what = f"rotation by {args.angle:g} deg about n(theta, phi)"

    thetas = np.radians(np.arange(0.0, 180.0 + 1e-9, args.step))
    phis = np.radians(np.arange(0.0, 360.0, args.step))
    print(f"{Path(args.xyz).name}: {len(labels)} atoms, "
          f"{'centred on ' + labels[0] if has_centre else 'centred on the centroid'}, "
          f"{'labels ignored' if args.ignore else f'{len(groups)} assignment group(s)'}")
    print(f"sds_dev of the structure against itself under {what}")
    print(f"grid: {len(thetas)} x {len(phis)} = {len(thetas) * len(phis)} directions at {args.step:g} deg, "
          f"{len(ops)} operation(s) each")

    start = time.perf_counter()
    values = landscape(points, groups, has_centre, ops, math.radians(args.spin), thetas, phis)
    print(f"evaluated in {time.perf_counter() - start:.1f} s: min {values.min():.4f}, "
          f"max {values.max():.4f}, mean {values.mean():.4f}\n")

    # Ligand directions from the centre, to name the nearest bond of every minimum.
    ligands = [(label, p / np.linalg.norm(p)) for label, p in zip(labels, points) if np.linalg.norm(p) > 1e-9]

    minima = local_minima(values, thetas, phis, merge_deg=3 * args.step)
    print(f"{len(minima)} local minima on the grid (antipodes counted separately):")
    print(f"{'#':>3} {'theta':>7} {'phi':>7} {'sds_dev':>9}   direction (x, y, z)        nearest ligand")
    for rank, (i, j) in enumerate(minima, 1):
        d = direction(thetas[i], phis[j])
        label, ld = min(ligands, key=lambda ligand: -ligand[1] @ d)
        off = math.degrees(math.acos(min(1.0, ld @ d)))
        print(f"{rank:>3} {math.degrees(thetas[i]):>7.1f} {math.degrees(phis[j]):>7.1f} {values[i, j]:>9.4f}   "
              f"({d[0]:+.3f}, {d[1]:+.3f}, {d[2]:+.3f})   {label} ({off:.1f} deg away)")

    seeds = fibonacci_hemisphere(args.seeds) if args.seeds else []
    norm = plt.Normalize(values.min(), values.max())
    thetas_deg, phis_deg = np.degrees(thetas), np.degrees(phis)

    fig = plt.figure(figsize=(15, 7))
    fig.suptitle(f"{Path(args.xyz).name}: sds_dev against itself under {what}\n"
                 f"{len(minima)} grid local minima (red, numbered by value), global minimum "
                 f"{values.min():.3f} at theta = {thetas_deg[minima[0][0]]:g}, phi = {phis_deg[minima[0][1]]:g}",
                 fontsize=11)

    ax_surface = fig.add_subplot(1, 2, 1, projection="3d")
    draw_surface(ax_surface, thetas_deg, phis_deg, values, minima, norm)
    ax_surface.set_title("as a surface over (theta, phi)", fontsize=10)

    ax_sphere = fig.add_subplot(1, 2, 2, projection="3d")
    fallback = draw_sphere(ax_sphere, thetas, phis, values, minima, norm, ligands, seeds)
    ax_sphere.set_title("on the unit sphere, in the structure's frame", fontsize=10)
    fig.colorbar(cm.ScalarMappable(norm=norm, cmap="viridis"), ax=ax_sphere, shrink=0.55, pad=0.1, label="sds_dev")

    handles = [plt.Line2D([], [], marker="o", color="red", linestyle="", markersize=6, label="local minimum")]
    for element in dict.fromkeys(strip_label(label) for label, _ in ligands):
        handles.append(plt.Line2D([], [], marker="^", color=element_colour(element, fallback), markeredgecolor="k",
                                  linestyle="", markersize=7, label=f"{element} ligand direction"))
    if args.seeds:
        handles.append(plt.Line2D([], [], marker="o", color="k", linestyle="", markersize=4,
                                  label=f"csom --seeds {args.seeds} candidate axis"))
    fig.legend(handles=handles, loc="lower center", ncol=len(handles), fontsize=8, frameon=False)

    fig.subplots_adjust(left=0.03, right=0.97, top=0.86, bottom=0.08, wspace=0.08)
    if args.save:
        fig.savefig(args.save, dpi=130)
        print(f"\nwrote {args.save}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
