"""
Check that the point-group tables in src/data/pgs.rs are well formed.

For every point group it checks that:

  1. the group is registered consistently: its name is in POINTGROUP_NAMES,
     get_pointgroup maps it to a POINTGROUP_* table, and every table is reachable;
  2. every matrix is orthogonal to full f64 precision (so no rounded literals)
     and has determinant +1 or -1;
  3. no matrix is listed twice and the identity is left out, since it is implicit
     (E's own table is the exception: it holds just the identity);
  4. the operations plus the identity are closed under multiplication, i.e. they
     form a group;
  5. the group has the order its name implies (C6v: 12, Oh: 48, Ih: 120, ...);
  6. it is the group its name says: its elements, classified by determinant and
     rotation angle, match the same group built from textbook generators. The
     classification does not depend on orientation, so this assumes no axis
     convention;
  7. its principal axis lies along z, as the csom axis search assumes.

Labels are only names in the --full and -o output, so they are checked as
warnings: every label should use the file's notation (2C6, 2C6^2, 3C2', C2(z),
i, 12S10^3, 3sigma_v, sigma(xy), ...), agree with its matrix (C_n^k is a proper
rotation by k*360/n degrees, S_n^k an improper rotation, sigma a reflection, i
the inversion), and follow the conjugacy classes: every operation of a class
carries the same label, no two classes share one, and the prefix (the 3 in
3C2') is the size of the class.

Usage:
    py check_pointgroup_tables.py [path/to/pgs.rs]

Exits with status 1 if any check fails; label warnings alone do not fail it.
"""

import ast
import math
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np

DEFAULT_PATH = Path(__file__).resolve().parent / "src" / "data" / "pgs.rs"

# Tight enough that a literal rounded to even 6 decimals fails.
ORTHOGONALITY_TOL = 1e-12
# Two matrices closer than this (largest entry difference) are the same operation.
MATCH_TOL = 1e-9

# std::f64::consts the tables may use.
STD_CONSTS = {"FRAC_1_SQRT_2": 1 / math.sqrt(2), "SQRT_2": math.sqrt(2)}

ROW_RE = re.compile(
    r'\(\s*"([^"]*)"\s*,\s*\[\s*\[([^\]]*)\]\s*,\s*\[([^\]]*)\]\s*,\s*\[([^\]]*)\]\s*\]\s*\)'
)
TABLE_RE = re.compile(r"pub static (POINTGROUP_\w+): [^=]+= &\[(.*?)\n\];", re.S)
ARM_RE = re.compile(r'"([^"]+)" => Some\((POINTGROUP_\w+)\)')
CONST_RE = re.compile(r"const (\w+): f64 = (.+?);")
NAMES_RE = re.compile(r"pub const POINTGROUP_NAMES: &\[&str\] = &\[(.*?)\];", re.S)

# The label notation: an optional class size, then E, i, C_n^k with primes or an axis for
# classes that would otherwise share a name (3C2', C2(z)), S_n^k, or sigma with a subscript or
# a plane (3sigma_v, sigma(xy)).
BODY = (r"(?:E|i|C\d+(?:\^\d+)?(?:'{1,2}|\([xyz]\))?|S\d+(?:\^\d+)?"
        r"|sigma(?:_[hvd])?(?:\((?:xy|xz|yz)\))?)")
STANDARD_LABEL_RE = re.compile(r"(?:[2-9]|[1-9]\d+)?" + BODY)


# --- Parsing ----------------------------------------------------------------------------------

def evaluate(expr, consts):
    """Evaluates a Rust float expression of literals and constants: + - * / and brackets."""
    def ev(node):
        if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
            return float(node.value)
        if isinstance(node, ast.Name) and node.id in consts:
            return consts[node.id]
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
            return -ev(node.operand) if isinstance(node.op, ast.USub) else ev(node.operand)
        if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Add, ast.Sub, ast.Mult, ast.Div)):
            a, b = ev(node.left), ev(node.right)
            return {ast.Add: a + b, ast.Sub: a - b, ast.Mult: a * b, ast.Div: a / b}[type(node.op)]
        raise ValueError(f"cannot evaluate '{expr.strip()}'")
    try:
        return ev(ast.parse(expr.strip(), mode="eval").body)
    except SyntaxError:
        raise ValueError(f"cannot evaluate '{expr.strip()}'") from None


def parse(text):
    """Returns (POINTGROUP_NAMES, {name: static}, {static: [(label, matrix)]}, parse errors)."""
    errors = []
    consts = dict(STD_CONSTS)
    for name, expr in CONST_RE.findall(text):  # in file order, so later ones can use earlier ones
        try:
            consts[name] = evaluate(expr, consts)
        except (ValueError, ZeroDivisionError) as e:
            errors.append(f"constant {name}: {e}")

    def number(token):
        return evaluate(token, consts)

    names_match = NAMES_RE.search(text)
    names = re.findall(r'"([^"]+)"', names_match.group(1)) if names_match else []
    arms = dict(ARM_RE.findall(text))

    tables = {}
    for static, body in TABLE_RE.findall(text):
        rows = []
        for label, *matrix_rows in ROW_RE.findall(body):
            try:
                rows.append((label, np.array([[number(t) for t in r.split(",")] for r in matrix_rows])))
            except (ValueError, KeyError) as e:
                errors.append(f"{static}: cannot read the matrix of '{label}': {e}")
        listed = body.count('("')
        if listed != len(rows):
            errors.append(f"{static}: read {len(rows)} of {listed} rows")
        tables[static] = rows
    return names, arms, tables, errors


# --- Reference groups -------------------------------------------------------------------------

def rz(angle):
    c, s = math.cos(angle), math.sin(angle)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def rotation(axis, angle):
    a = np.asarray(axis, float) / np.linalg.norm(axis)
    k = np.array([[0.0, -a[2], a[1]], [a[2], 0.0, -a[0]], [-a[1], a[0], 0.0]])
    return np.eye(3) + math.sin(angle) * k + (1 - math.cos(angle)) * k @ k


SIGMA_H = np.diag([1.0, 1.0, -1.0])
SIGMA_V = np.diag([1.0, -1.0, 1.0])
C2_X = np.diag([1.0, -1.0, -1.0])
INVERSION = -np.eye(3)
C3_BODY_DIAGONAL = np.array([[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
SIGMA_D = np.array([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
C5_OFF_AXIS = rotation([0.0, 2.0, 1.0], 2 * math.pi / 5)  # a second C5 axis of the icosahedron


def reference(name):
    """(generators, order, operation that puts the principal axis along z) for a group name,
    or None for a name this script does not know."""
    special = {
        "E": ([], 1, None),
        "Ci": ([INVERSION], 2, None),
        "Cs": ([SIGMA_H], 2, SIGMA_H),
        "T": ([rz(math.pi), C3_BODY_DIAGONAL], 12, rz(math.pi)),
        "Td": ([rz(math.pi), C3_BODY_DIAGONAL, SIGMA_D], 24, rz(math.pi)),
        "Th": ([rz(math.pi), C3_BODY_DIAGONAL, INVERSION], 24, rz(math.pi)),
        "O": ([rz(math.pi / 2), C3_BODY_DIAGONAL], 24, rz(math.pi / 2)),
        "Oh": ([rz(math.pi / 2), C3_BODY_DIAGONAL, INVERSION], 48, rz(math.pi / 2)),
        "I": ([rz(2 * math.pi / 5), C5_OFF_AXIS], 60, rz(2 * math.pi / 5)),
        "Ih": ([rz(2 * math.pi / 5), C5_OFF_AXIS, INVERSION], 120, rz(2 * math.pi / 5)),
    }
    if name in special:
        return special[name]

    m = re.fullmatch(r"([CDS])(\d+)([vhd]?)", name)
    if not m:
        return None
    kind, n, sub = m.group(1), int(m.group(2)), m.group(3)
    cn = rz(2 * math.pi / n)
    if kind == "C":
        return {"": ([cn], n, cn), "v": ([cn, SIGMA_V], 2 * n, cn), "h": ([cn, SIGMA_H], 2 * n, cn)}.get(sub)
    if kind == "D":
        s2n = SIGMA_H @ rz(math.pi / n)
        return {"": ([cn, C2_X], 2 * n, cn), "h": ([cn, C2_X, SIGMA_H], 4 * n, cn),
                "d": ([s2n, C2_X], 4 * n, s2n)}.get(sub)
    if kind == "S" and sub == "" and n % 2 == 0:
        return [SIGMA_H @ cn], n, SIGMA_H @ cn
    return None


def closure(generators):
    elements = [np.eye(3)]
    frontier = [np.eye(3)]
    while frontier:
        new = []
        for e in frontier:
            for g in generators:
                p = g @ e
                if not any(np.abs(p - x).max() < MATCH_TOL for x in elements):
                    elements.append(p)
                    new.append(p)
        frontier = new
    return elements


def classify(m):
    """(determinant, rotation angle in degrees of the proper part): invariant under a change of
    axes, so two tables of the same group agree on it whatever their orientation."""
    det = round(np.linalg.det(m))
    angle = math.degrees(math.acos(max(-1.0, min(1.0, (np.trace(det * m) - 1) / 2))))
    return det, round(angle, 3)


def fold(angle):
    angle %= 360
    return round(360 - angle if angle > 180 else angle, 3)


def label_class(body):
    """The (determinant, angle) a label promises, or None if it cannot be read."""
    if body == "E":
        return 1, 0.0
    if body == "i":
        return -1, 0.0
    if re.fullmatch(r"sigma(?:_[hvd])?(?:\((?:xy|xz|yz)\))?", body):
        return -1, 180.0
    m = re.fullmatch(r"C(\d+)(?:\^(\d+))?(?:'{1,2}|\([xyz]\))?", body)
    if m:
        n, k = int(m.group(1)), int(m.group(2) or 1)
        return 1, fold(360 * k / n)
    m = re.fullmatch(r"S(\d+)(?:\^(\d+))?", body)
    if m:
        n, k = int(m.group(1)), int(m.group(2) or 1)
        return (-1, fold(180 + 360 * k / n)) if k % 2 else (1, fold(360 * k / n))
    return None


def conjugacy_classes(stack):
    """The conjugacy classes of a group given as an (n, 3, 3) array, as lists of indices."""
    classes, seen = [], set()
    for i, x in enumerate(stack):
        if i not in seen:
            conjugates = np.einsum("gij,jk,glk->gil", stack, x, stack)  # g x g^-1
            members = sorted({int(np.abs(stack - c).max(axis=(1, 2)).argmin()) for c in conjugates})
            seen.update(members)
            classes.append(members)
    return classes


# --- Checks -----------------------------------------------------------------------------------

def check_group(name, rows):
    """Returns (errors, warnings) for one group."""
    errors, warnings = [], []
    if not rows:
        return [f"{name}: the table is empty"], warnings

    off = [(np.abs(m @ m.T - np.eye(3)).max(), label) for label, m in rows]
    bad = [(err, label) for err, label in off if err > ORTHOGONALITY_TOL]
    if bad:
        worst, label = max(bad)
        errors.append(f"{len(bad)} of {len(rows)} matrices are not orthogonal, worst '{label}' "
                      f"off by {worst:.1e}: rounded literals?")
    for label, m in rows:
        if abs(abs(np.linalg.det(m)) - 1) > ORTHOGONALITY_TOL and not bad:
            errors.append(f"'{label}' has determinant {np.linalg.det(m):.6f}")
    if errors:
        return errors, warnings  # the group checks below assume exact matrices

    is_identity = [np.abs(m - np.eye(3)).max() < MATCH_TOL for _, m in rows]
    if name == "E":
        if len(rows) != 1 or not all(is_identity):
            errors.append("E's table should hold exactly one row, the identity")
        elements = [m for _, m in rows]
        element_labels = [label for label, _ in rows]
    else:
        if any(is_identity):
            errors.append(f"lists the identity ('{rows[is_identity.index(True)][0]}'), which is implicit")
        kept = [row for row, identity in zip(rows, is_identity) if not identity]
        elements = [np.eye(3)] + [m for _, m in kept]
        element_labels = [None] + [label for label, _ in kept]

    stack = np.array(elements)
    distance = np.abs(stack[:, None] - stack[None]).max(axis=(2, 3))
    np.fill_diagonal(distance, np.inf)
    for a, b in zip(*np.nonzero(np.triu(distance < MATCH_TOL))):
        errors.append(f"element {a} and element {b} are the same operation")

    missing = 0
    for a in stack:
        products = np.einsum("ij,bjk->bik", a, stack)
        nearest = np.abs(products[:, None] - stack[None]).max(axis=(2, 3)).min(axis=1)
        missing += int((nearest > MATCH_TOL).sum())
    if missing:
        errors.append(f"not closed under multiplication: {missing} products fall outside the table")

    ref = reference(name)
    if ref is None:
        errors.append("unknown group name: no reference to compare against")
        return errors, warnings
    generators, order, principal = ref
    if len(elements) != order:
        errors.append(f"has {len(elements)} elements, {name} has {order}")
    expected = Counter(classify(m) for m in closure(generators))
    found = Counter(classify(m) for m in elements)
    if found != expected:
        diff = (found - expected) + (expected - found)
        errors.append(f"is not {name}: its elements differ in {sorted(diff.elements())} (det, angle)")
    if principal is not None and not any(np.abs(principal - m).max() < MATCH_TOL for m in elements):
        errors.append("its principal axis is not along z")

    # Labels (warnings only).
    for label, m in rows:
        promised = label_class(re.fullmatch(r"\d*(.*)", label).group(1))
        if promised is None:
            warnings.append(f"'{label}' cannot be read as an operation")
            continue
        if not STANDARD_LABEL_RE.fullmatch(label):
            warnings.append(f"'{label}' is not in the file's usual notation")
        if promised != classify(m):
            det, angle = classify(m)
            kind = "a rotation" if det == 1 else "an improper operation"
            warnings.append(f"'{label}' is {kind} by {angle:g} deg, which the label does not describe")

    classes_of_label = Counter()
    for members in conjugacy_classes(stack):
        labels = sorted({element_labels[i] for i in members} - {None})
        if not labels:
            continue  # the implicit identity
        if len(labels) > 1:
            warnings.append(f"one class of {len(members)} is labelled {', '.join(map(repr, labels))}")
            continue
        label = labels[0]
        classes_of_label[label] += 1
        prefix = re.match(r"\d*", label).group()
        if (int(prefix) if prefix else 1) != len(members):
            warnings.append(f"'{label}' labels a class of {len(members)}, its prefix says {prefix or 1}")
        # A rotation about z without its inverse in the class must turn counterclockwise, the
        # sense the file's header documents (C_n^k: 2*pi*k/n about +z, S_n^k: sigma_h^k times it).
        power = re.fullmatch(r"\d*([CS])(\d+)(?:\^(\d+))?", label)
        if power and len(members) == 1:
            kind, n, k = power.group(1), int(power.group(2)), int(power.group(3) or 1)
            expected = rz(2 * math.pi * k / n)
            if kind == "S" and k % 2:
                expected = SIGMA_H @ expected
            if np.abs(stack[members[0]] - expected).max() > MATCH_TOL:
                warnings.append(f"'{label}' turns the wrong way: C_n^k is counterclockwise about +z")
    for label, count in classes_of_label.items():
        if count > 1:
            warnings.append(f"'{label}' labels {count} different classes")
    return errors, sorted(set(warnings), key=warnings.index)


def main():
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_PATH
    names, arms, tables, errors = parse(path.read_text(encoding="utf-8"))
    print(f"checking {path}: {len(names)} point groups\n")

    for name in names:
        if name not in arms:
            errors.append(f"{name}: in POINTGROUP_NAMES but not in get_pointgroup")
    for name, static in arms.items():
        if name not in names:
            errors.append(f"{name}: in get_pointgroup but not in POINTGROUP_NAMES")
        if static not in tables:
            errors.append(f"{name}: get_pointgroup points at {static}, which does not exist")
    for static in tables:
        if static not in arms.values():
            errors.append(f"{static}: not reachable from get_pointgroup")
    for e in errors:
        print(f"error: {e}")
    if errors:
        print()

    all_warnings = []
    for name in names:
        static = arms.get(name)
        if static not in tables:
            continue
        group_errors, warnings = check_group(name, tables[static])
        note = f"  ({len(warnings)} label warning{'s' if len(warnings) != 1 else ''})" if warnings else ""
        print(f"{name:<4} {len(tables[static]) + (name != 'E'):>4} elements  "
              f"{'FAIL' if group_errors else 'ok'}{note}")
        for e in group_errors:
            print(f"       error: {e}")
        errors += [f"{name}: {e}" for e in group_errors]
        all_warnings += [f"{name}: {w}" for w in warnings]

    if all_warnings:
        print("\nlabel warnings (names in the --full and -o output only):")
        for w in all_warnings:
            print(f"  {w}")
    print(f"\n{len(errors)} error(s), {len(all_warnings)} label warning(s)")
    sys.exit(1 if errors else 0)


if __name__ == "__main__":
    main()
