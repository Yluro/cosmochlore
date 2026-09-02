use nalgebra::{Vector3};



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

/// Solves the min-cost assignment problem via the Jonker-Volgenant / Kuhn-Munkres
/// algorithm with dual potentials. O(n^3). Assumes a square cost matrix.
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
pub fn best_permutation(a: &[Vector3<f64>], b: &[Vector3<f64>]) -> (Vec<Vector3<f64>>, Vec<usize>) {
    // Since the matrix is square (bijection), the |a_i|^2 + |b_j|^2 terms
    // of |a_i - b_j|^2 are constant over any assignment and can be dropped.
    // Minimizing squared distance <=> maximizing dot product <=> minimizing -dot.
    let cost: Vec<Vec<f64>> = a
        .iter()
        .map(|ai| b.iter().map(|bj| -ai.dot(bj)).collect())
        .collect();

    let assignment = hungarian(&cost);
    let b_permuted: Vec<Vector3<f64>> = assignment.iter().map(|&j| b[j]).collect();

    (b_permuted, assignment)
}