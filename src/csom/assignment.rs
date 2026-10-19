use nalgebra::Vector3;
use std::collections::BTreeMap;

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
        let mut minv = vec![f64::INFINITY; n + 1]; // best reduced cost found so far to reach col j
        let mut used = vec![false; n + 1]; // columns already visited this round

        loop {
            used[j0] = true;
            let i0 = p[j0]; // row currently sitting at column j0
            let mut delta = f64::INFINITY; // smallest reduced cost among unvisited columns
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

/// Min squared-distance matching of A to B: `a[i]` pairs with `b[result[i]]`.
pub fn best_permutation(a: &[Vector3<f64>], b: &[Vector3<f64>]) -> Vec<usize> {
    hungarian(&sq_dist_cost(a, b))
}

/// Splits an array of points A and labels L given the different labels of L
///
/// ["Cl", "Cl2", "O"] -> ["Cl", "Cl"], ["O",]
///
/// `BTreeMap` keeps the group order fixed across runs.
fn split_by_atoms(labels: &[String]) -> BTreeMap<String, Vec<usize>> {
    let mut result: BTreeMap<String, Vec<usize>> = BTreeMap::new();

    for (i, label) in labels.iter().enumerate() {
        result
            .entry(label.clone()) // Get the entry for the element
            .or_default() // If element is not present, insert an empty Vec<>
            .push(i); // Push the Vector3 to the first element.
    }
    result
}

/// Precomputes the index groups `best_permutation_multiple_atoms` restricts its search to.
pub(crate) fn group_by_label(labels: &[String]) -> Vec<Vec<usize>> {
    split_by_atoms(labels).into_values().collect()
}

/// Hungarian matching of the atoms in `indices`, written into `pairing`.
/// With `has_centre`, atom `0` is pinned to itself.
fn assign_group(
    indices: &[usize],
    a: &[Vector3<f64>],
    b: &[Vector3<f64>],
    has_centre: bool,
    pairing: &mut [usize],
) {
    let (pinned, rest): (Vec<usize>, Vec<usize>) = if has_centre {
        indices.iter().copied().partition(|&i| i == 0)
    } else {
        (Vec::new(), indices.to_vec())
    };

    for i in pinned {
        pairing[i] = i;
    }

    if rest.is_empty() {
        return;
    }

    let a_subset: Vec<Vector3<f64>> = rest.iter().map(|&i| a[i]).collect();
    let b_subset: Vec<Vector3<f64>> = rest.iter().map(|&i| b[i]).collect();

    // `assignment` indexes into the subset; map it back through `rest`.
    let assignment = best_permutation(&a_subset, &b_subset);

    for (&i, &k) in rest.iter().zip(&assignment) {
        pairing[i] = rest[k];
    }
}

/// Min squared-distance matching of A to B within each label group (from [`group_by_label`]).
///
/// Returns `pairing`: `a[i]` pairs with `b[pairing[i]]`. With `ignore_labels` the groups are
/// ignored, so `pairing[i]` may be another element. With `has_centre`, atom `0` is pinned.
pub fn best_permutation_multiple_atoms(
    a: &[Vector3<f64>],
    b: &[Vector3<f64>],
    groups: &[Vec<usize>],
    ignore_labels: bool,
    has_centre: bool,
) -> Vec<usize> {
    debug_assert_eq!(a.len(), b.len());

    let mut pairing = vec![0usize; a.len()];

    if ignore_labels {
        let all_indices: Vec<usize> = (0..a.len()).collect();
        assign_group(&all_indices, a, b, has_centre, &mut pairing);
    } else {
        for indices in groups {
            assign_group(indices, a, b, has_centre, &mut pairing);
        }
    }

    pairing
}
