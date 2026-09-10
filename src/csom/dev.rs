use std::collections::HashMap;
use nalgebra::Vector3;
use crate::csom::{CsomError, CsomOperation};
use crate::csom::io::strip_all_labels;
use crate::data::pgs::{get_pointgroup, to_matrix3};

/// Shape deviation for two given shapes. It returns squared-distance-sum deviation (0 to 100)
///
/// Assumes centered and normalized points.
/// Assumes correct point-to-point correspondence.
pub(crate) fn sds_dev(reference: &[Vector3<f64>], problem: &[Vector3<f64>]) -> f64 {
    let n = reference.len() as f64;
    let problem_centroid = problem.iter().sum::<Vector3<f64>>() / n;

    let denominator: f64 = problem.iter()
        .map(|p| (p - problem_centroid).norm_squared())
        .sum();

    let numerator: f64 = reference.iter().zip(problem)
        .map(
            |(p , q)| (q - p).norm_squared()
        ).sum();

    100.0 * numerator / denominator
}

/// Returns the cost matrix for two given sets of points A and B.
///
/// Minimizing squared distance <=> maximizing dot product <=> minimizing -dot.
fn sq_dist_cost(a: &[Vector3<f64>], b: &[Vector3<f64>]) -> Vec<Vec<f64>> {
    let cost: Vec<Vec<f64>> = a
        .iter()
        .map(|ai| b.iter().map(|bj| -ai.dot(bj)).collect())
        .collect();
    cost
}


/// Jonker-Volgenant / Kuhn-Munkres algorithm with dual potentials. O(n^3). Assumes a square cost matrix.
/// Returns, for each row i, the column it's assigned to (0-indexed).
pub fn hungarian(cost: &[Vec<f64>]) -> Vec<usize> {
    let n = cost.len();
    const INF: f64 = f64::INFINITY;

    // Dual potentials for rows (u) and columns (v). Kept so that
    // reduced_cost = cost[i][j] - u[i] - v[j] is always >= 0.
    let mut u = vec![0.0; n + 1];
    let mut v = vec![0.0; n + 1];

    // p[j] = which row is currently matched to column j (1-indexed, 0 = dummy/none).
    let mut p = vec![0usize; n + 1];
    // way[j] = predecessor column on the augmenting path that reached column j.
    let mut way = vec![0usize; n + 1];

    // Process one row at a time, each time growing an augmenting path
    // from a dummy column (0) until a free column is found.
    for i in 1..=n {
        p[0] = i; // temporarily "assign" row i to the dummy column
        let mut j0 = 0usize; // current column in the path search
        let mut minv = vec![INF; n + 1]; // best reduced cost found so far to reach col j
        let mut used = vec![false; n + 1]; // columns already visited this round

        loop {
            used[j0] = true;
            let i0 = p[j0]; // row currently sitting at column j0
            let mut delta = INF; // smallest reduced cost among unvisited columns
            let mut j1 = 0usize; // column achieving that smallest cost

            // Relax distances to every unvisited column via row i0
            // (Dijkstra-style expansion, like a shortest-path search).
            for j in 1..=n {
                if !used[j] {
                    let cur = cost[i0 - 1][j - 1] - u[i0] - v[j];
                    if cur < minv[j] {
                        minv[j] = cur;
                        way[j] = j0; // remember how we got here, for backtracking later
                    }
                    if minv[j] < delta {
                        delta = minv[j];
                        j1 = j;
                    }
                }
            }

            // Update potentials: tighten for visited columns, relax for the rest,
            // keeping all reduced costs valid (non-negative) after moving by delta.
            for j in 0..=n {
                if used[j] {
                    u[p[j]] += delta;
                    v[j] -= delta;
                } else {
                    minv[j] -= delta;
                }
            }

            j0 = j1; // step onto the closest next column
            if p[j0] == 0 {
                // Reached a column with no row assigned yet -> augmenting path complete.
                break;
            }
        }

        // Backtrack along the path, flipping assignments so each column
        // on the path now belongs to the row that should own it.
        loop {
            let j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
            if j0 == 0 {
                break; // back at the dummy column, done for this row
            }
        }
    }

    // Convert 1-indexed p[col] = row into 0-indexed result[row] = col.
    let mut result = vec![0usize; n];
    for j in 1..=n {
        if p[j] != 0 {
            result[p[j] - 1] = j - 1;
        }
    }
    result
}

/// Finds the best one-to-one matching between point sets A and B (equal length)
/// that minimizes total squared distance, and reorders B accordingly.
///
/// Returns (B reordered to best match A, permutation indices into B)
pub fn best_permutation(a: &[Vector3<f64>], b: &[Vector3<f64>]) -> (Vec<Vector3<f64>>, Vec<usize>) {
    let cost = sq_dist_cost(a, b);

    let assignment = hungarian(&cost);
    let b_permuted: Vec<Vector3<f64>> = assignment.iter().map(|&j| b[j]).collect();

    (b_permuted, assignment)
}


/// Splits an array of points A and labels L given the different labels of L
///
/// ["Cl", "Cl2", "O"] -> ["Cl", "Cl"], ["O"]
fn split_by_atoms(labels: &[String]) -> HashMap<String, Vec<usize>> {
    let mut result: HashMap<String, Vec<usize>> = HashMap::new();

    for (i, label) in labels.iter().enumerate() {
        result
            .entry(label.clone())// Get the entry for the element
            .or_default() // If element is not present, insert an empty Vec<>
            .push(i); // Push the Vector3 to the first element.
    }
    result
}



/// Finds the best permutation of b to a using the Hungarian algorithm.
///
/// When `has_centre` is true, index `0` is pulled out of the group and matched to
/// itself directly.
///
/// Returns (a_idx, b_idx, a_point, b_point), where each index is that point's own index
/// into the underlying `a`/`b`/`labels` arrays. `b_idx` is the atom `b_point` belongs to,
/// which is *not* `a_idx` whenever the assignment moved that atom somewhere else.
fn assign_group(
    indices: &[usize],
    a: &[Vector3<f64>],
    b: &[Vector3<f64>],
    has_centre: bool,
    pairs: &mut Vec<(usize, usize, Vector3<f64>, Vector3<f64>)>,
) {
    let (pinned, rest): (Vec<usize>, Vec<usize>) = if has_centre {
        indices.iter().copied().partition(|&i| i == 0)
    } else {
        (Vec::new(), indices.to_vec())
    };

    for i in pinned {
        pairs.push((i, i, a[i], b[i]));
    }

    if rest.is_empty() {
        return;
    }

    let a_subset: Vec<Vector3<f64>> = rest.iter().map(|&i| a[i]).collect();
    let b_subset: Vec<Vector3<f64>> = rest.iter().map(|&i| b[i]).collect();

    // `assignment[k]` indexes into the subset, so `rest[assignment[k]]` is the atom whose
    // b-point ended up paired with `rest[k]`.
    let (perm_b, assignment) = best_permutation(&a_subset, &b_subset);

    pairs.extend(
        rest.iter()
            .zip(a_subset)
            .zip(perm_b)
            .zip(&assignment)
            .map(|(((&i, pa), pb), &k)| (i, rest[k], pa, pb)),
    );
}

/// Finds the best one-to-one matching between the subsets A and B by atom type. (equal length)
/// that minimizes total squared distance, and reorders B accordingly.
///
/// If `ignore_labels` is true, atom labels are not used to restrict the search. If `has_centre` is true, atom
/// index `0` is pinned to itself instead of entering the Hungarian assignment.
///
/// Returns (A and B reordered to best match, pairing), where `pairing[i]` is the atom whose
/// B-point got matched to `a[i]`. With labels honoured that atom always shares atom `i`'s
/// element; with `ignore_labels` it need not, so a caller reporting the matched points must
/// take each one's identity from `pairing[i]` rather than from `i`.
pub fn best_permutation_multiple_atoms(
    a: &[Vector3<f64>],
    b: &[Vector3<f64>],
    labels: &[String],
    ignore_labels: bool,
    has_centre: bool,
) -> (Vec<Vector3<f64>>, Vec<Vector3<f64>>, Vec<usize>)
{
    debug_assert_eq!(labels.len(), b.len());
    debug_assert_eq!(labels.len(), a.len());

    let mut pairs: Vec<(usize, usize, Vector3<f64>, Vector3<f64>)> = Vec::new();

    if ignore_labels {
        let all_indices: Vec<usize> = (0..labels.len()).collect();
        assign_group(&all_indices, a, b, has_centre, &mut pairs);
    } else {
        for (_, indices) in split_by_atoms(labels) {
            assign_group(&indices, a, b, has_centre, &mut pairs);
        }
    }

    debug_assert_eq!(pairs.len(), labels.len());

    let mut pairing = vec![0usize; labels.len()];
    for &(a_idx, b_idx, _, _) in &pairs {
        pairing[a_idx] = b_idx;
    }

    // Sort pairs by ascending order of z- y- x- values so output from hashmap is deterministic.
    pairs.sort_by(
        |(_, _, a1, _), (_, _, a2, _)| {
            a1.z.total_cmp(&a2.z)
                .then(a1.y.total_cmp(&a2.y))
                .then(a1.x.total_cmp(&a2.x))
        }
    );

    let mut final_a = Vec::with_capacity(pairs.len());
    let mut final_b = Vec::with_capacity(pairs.len());
    for (_, _, pa, pb) in pairs {
        final_a.push(pa);
        final_b.push(pb);
    }

    (final_a, final_b, pairing)
}


/// Deviation of `points` from every individual symmetry operation of `pg`.
pub(crate) fn point_group_operation_deviations(
    points: &[Vector3<f64>],
    labels: &[String],
    pg: &str,
    ignore_labels: bool,
    has_centre: bool,
) -> Result<Vec<CsomOperation>, CsomError> {
    let ops = get_pointgroup(pg).ok_or_else(|| CsomError::WrongSpaceGroup { pg: pg.to_string() })?;
    let stripped = strip_all_labels(labels);

    Ok(ops.iter().map(|(name, matrix)| {
        let sym_op = to_matrix3(*matrix);
        let image: Vec<Vector3<f64>> = points.iter().map(|p| sym_op * p).collect();

        // The assignment only decides which image point each atom is *scored against*; it
        // never changes whose image a point is, so `image` stays in the input's atom order
        // and the pairing is reported separately.
        let (a, b, pairing) = best_permutation_multiple_atoms(points, &image, &stripped, ignore_labels, has_centre);

        CsomOperation {
            name: name.to_string(),
            matrix: sym_op,
            deviation: sds_dev(&a, &b),
            image,
            pairing,
        }
    }).collect())
}

/// Average CSM deviation of `points` from every symmetry operation of point group `pg`.
pub fn point_group_dev(
    points: &[Vector3<f64>],
    labels: &[String],
    pg: &str,
    ignore_labels: bool,
    has_centre: bool,
) -> Result<f64, CsomError> {
    let devs = point_group_operation_deviations(points, labels, pg, ignore_labels, has_centre)?;

    Ok(devs.iter().map(|op| op.deviation).sum::<f64>() / devs.len() as f64)
}

