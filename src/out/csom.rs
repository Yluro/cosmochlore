//! `csom` output: the point-group table, `--table` and `--full` CSVs, and the `--operated` xyz and
//! merged mol2.

use super::csv::write_csv;
use super::table::{Row, Table, Width};
use super::{create_output, output_stem, table_path};
use crate::cli::CsomArgs;
use crate::csom::prepare::strip_all_labels;
use crate::csom::types::{CsomResult, CsomRun};
use crate::data::elements::covalent_radius;
use nalgebra::{Matrix3, Vector3};
use std::io::{Result, Write};

pub fn csom(run: &CsomRun, args: &CsomArgs, silent: bool) -> Result<()> {
    print(&run.results, &args.name, silent);
    if args.table {
        write_table_csv(&run.results, &args.name, silent)?;
    }
    if args.full {
        write_details(run, &args.name, silent)?;
    }
    if args.operated {
        write_operated(run, &args.name, silent)?;
        write_merged_mol2(run, &args.name, silent)?;
    }
    Ok(())
}

pub(super) fn print(results: &[CsomResult], file_name: &str, silent: bool) {
    let mut table = Table::new(&[Width::Fixed(12), Width::Fixed(10)]);
    table.row(Row::Text(format!("\nInput file: {file_name}")));
    table.row(Row::DoubleRule);
    table.row(Row::from(["Point group", "CSoM"]));
    table.row(Row::Rule);
    for r in results {
        table.row(Row::from([
            r.point_group.clone(),
            format!("{:.3}", r.deviation),
        ]));
    }
    table.row(Row::Rule);
    table.print(silent);
}

/// `<file>_csom_table.csv`: point group, deviation and refined rotation per group.
pub(super) fn write_table_csv(results: &[CsomResult], file_name: &str, silent: bool) -> Result<()> {
    let rows: Vec<Vec<String>> = results
        .iter()
        .map(|r| {
            let deviation = format!("{:.3}", r.deviation);
            vec![
                r.point_group.clone(),
                deviation,
                format_matrix3(&r.rotation),
            ]
        })
        .collect();
    let header = ["PointGroup", "Dev", "Rotation Matrix"];
    let path = table_path(file_name, "csom");
    write_csv(&path, "output table", &header, &rows, silent)
}

/// `<file>_<point group>_details.csv`: name, matrix, deviation and atom pairing per operation.
fn write_details(run: &CsomRun, file_name: &str, silent: bool) -> Result<()> {
    let header = ["name", "op_matrix", "dev", "pairing"];
    for result in &run.results {
        let rows: Vec<Vec<String>> = result
            .operations
            .iter()
            .map(|op| {
                vec![
                    op.name.clone(),
                    format_matrix3(&op.matrix),
                    format!("{:.3}", op.deviation),
                    format_pairing(&op.pairing, &run.labels),
                ]
            })
            .collect();
        let path = format!(
            "{}_{}_details.csv",
            output_stem(file_name),
            result.point_group
        );
        write_csv(&path, "operation details", &header, &rows, silent)?;
    }
    Ok(())
}

fn format_matrix3(m: &Matrix3<f64>) -> String {
    let row = |i: usize| format!("{:.4} {:.4} {:.4}", m[(i, 0)], m[(i, 1)], m[(i, 2)]);
    format!("[{}; {}; {}]", row(0), row(1), row(2))
}

/// Formats an operation's atom pairing as `Source>Target` entries in the input's atom order:
/// `A>B` means the operation carries atom `A` onto the site of atom `B`.
fn format_pairing(pairing: &[usize], labels: &[String]) -> String {
    // `pairing[target] = source`, so invert it to walk the atoms in input order.
    let mut target_of = vec![0usize; pairing.len()];
    for (target, &source) in pairing.iter().enumerate() {
        target_of[source] = target;
    }
    let pairs: Vec<String> = (0..pairing.len())
        .map(|source| format!("{}>{}", labels[source], labels[target_of[source]]))
        .collect();
    pairs.join(" ")
}

/// The structure carried by each operation of `result`, back in the input frame, after an
/// "E" block holding the input itself: `(operation name, deviation, coordinates in file order)`.
fn operated_blocks<'a>(
    result: &'a CsomResult,
    original_coords: &[Vector3<f64>],
) -> Vec<(&'a str, f64, Vec<Vector3<f64>>)> {
    // `rotation` is orthogonal, so its transpose undoes the alignment.
    let rotation_inv = result.rotation.transpose();
    let mut blocks = vec![("E", 0.0, original_coords.to_vec())];
    blocks.extend(result.operations.iter().map(|op| {
        let points = op
            .image
            .iter()
            .map(|p| (rotation_inv * p) / result.scale + result.centroid)
            .collect();
        (op.name.as_str(), op.deviation, points)
    }));
    blocks
}

/// `<file>_<point group>_operated.xyz`: the input and its image under every operation.
fn write_operated(run: &CsomRun, file_name: &str, silent: bool) -> Result<()> {
    let labels = &run.labels;
    for result in &run.results {
        let path = format!(
            "{}_{}_operated.xyz",
            output_stem(file_name),
            result.point_group
        );
        let mut file = create_output(&path, "operated coordinates", silent)?;

        // The tables in data/pgs.rs omit E, so the "E" block is the input with dev 0. Every atom
        // keeps its own label: an image landing on another atom's site shows up as such.
        for (name, deviation, points) in operated_blocks(result, &run.coords) {
            writeln!(file, "{}", labels.len())?;
            writeln!(
                file,
                "{} {} dev = {:.3}",
                result.point_group, name, deviation
            )?;
            for (label, point) in labels.iter().zip(&points) {
                writeln!(
                    file,
                    "{}  {:.6}  {:.6}  {:.6}",
                    label, point.x, point.y, point.z
                )?;
            }
            writeln!(file)?;
        }
    }
    Ok(())
}

/// Two atoms are bonded in the merged .mol2 when they are closer than this factor times the
/// sum of their covalent radii (the same criterion as RDKit's connectivity perception).
const BOND_TOLERANCE: f64 = 1.3;

/// Bonds of one block of the merged .mol2: every pair within bond distance. times the sum of their
/// covalent radii. Every ligand is bonded to the central atom regardless of the distance.
fn perceive_bonds(
    elements: &[String],
    coords: &[Vector3<f64>],
    has_centre_atom: bool,
) -> Vec<(usize, usize)> {
    let radii: Vec<f64> = elements.iter().map(|e| covalent_radius(e)).collect();

    let mut bonds = Vec::new();
    for i in 0..coords.len() {
        for j in (i + 1)..coords.len() {
            let cutoff = BOND_TOLERANCE * (radii[i] + radii[j]);
            if (has_centre_atom && i == 0) || (coords[i] - coords[j]).norm() <= cutoff {
                bonds.push((i, j));
            }
        }
    }
    bonds
}

/// `<file>_<point group>_merged.mol2`: the blocks of [`write_operated`] overlaid in one
/// structure, bonds from [`perceive_bonds`].
fn write_merged_mol2(run: &CsomRun, file_name: &str, silent: bool) -> Result<()> {
    let stem = output_stem(file_name);
    let labels = &run.labels;
    let elements = strip_all_labels(labels);
    let n_per_block = labels.len();
    let bonds = perceive_bonds(&elements, &run.coords, run.has_centre);

    for result in &run.results {
        let path = format!("{}_{}_merged.mol2", stem, result.point_group);
        let mut file = create_output(&path, "merged mol2", silent)?;
        let blocks = operated_blocks(result, &run.coords);

        writeln!(file, "@<TRIPOS>MOLECULE")?;
        writeln!(file, "{}_{}_merged", stem, result.point_group)?;
        writeln!(
            file,
            "{} {} {} 0 0",
            n_per_block * blocks.len(),
            bonds.len() * blocks.len(),
            blocks.len()
        )?;
        writeln!(file, "SMALL")?;
        writeln!(file, "NO_CHARGES")?;
        writeln!(file)?;

        writeln!(file, "@<TRIPOS>ATOM")?;
        let mut atom_id = 0usize;
        for (block_idx, (name, _, points)) in blocks.iter().enumerate() {
            let subst_id = block_idx + 1;
            for ((label, element), point) in labels.iter().zip(&elements).zip(points) {
                atom_id += 1;
                writeln!(
                    file,
                    "{} {}.{} {:.6} {:.6} {:.6} {} {} {} 0.0000",
                    atom_id, label, name, point.x, point.y, point.z, element, subst_id, name
                )?;
            }
        }

        if !bonds.is_empty() {
            writeln!(file, "@<TRIPOS>BOND")?;
            let mut bond_id = 0usize;
            for block in 0..blocks.len() {
                let base = block * n_per_block;
                for &(i, j) in &bonds {
                    bond_id += 1;
                    writeln!(file, "{} {} {} 1", bond_id, base + i + 1, base + j + 1)?;
                }
            }
        }

        writeln!(file, "@<TRIPOS>SUBSTRUCTURE")?;
        for (block, (name, _, _)) in blocks.iter().enumerate() {
            writeln!(file, "{} {} {}", block + 1, name, block * n_per_block + 1)?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn labels(symbols: &[&str]) -> Vec<String> {
        symbols.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn perceive_bonds_joins_a_ring_without_a_centre() {
        // A regular hexagon of carbons with 1.39 A edges, like a benzene ring.
        let coords: Vec<Vector3<f64>> = (0..6)
            .map(|k| {
                let angle = k as f64 * std::f64::consts::FRAC_PI_3;
                Vector3::new(1.39 * angle.cos(), 1.39 * angle.sin(), 0.0)
            })
            .collect();

        let bonds = perceive_bonds(&labels(&["C"; 6]), &coords, false);

        // Neighbours only: the 1,3 (2.41 A) and 1,4 (2.78 A) pairs are well past the cutoff.
        assert_eq!(bonds, vec![(0, 1), (0, 5), (1, 2), (2, 3), (3, 4), (4, 5)]);
    }

    #[test]
    fn perceive_bonds_always_joins_the_centre_to_every_ligand() {
        // An octahedron with 3.0 A bonds: past the covalent cutoff for Fe-N (2.9 A), so
        // only the centre rule can draw the polyhedron.
        let mut coords = vec![Vector3::zeros()];
        for axis in [Vector3::x(), Vector3::y(), Vector3::z()] {
            coords.push(3.0 * axis);
            coords.push(-3.0 * axis);
        }
        let elements = labels(&["Fe", "N", "N", "N", "N", "N", "N"]);

        let bonds = perceive_bonds(&elements, &coords, true);
        assert_eq!(bonds, vec![(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6)]);

        // Without a declared centre the same coordinates are simply seven loose atoms.
        assert!(perceive_bonds(&elements, &coords, false).is_empty());
    }

    #[test]
    fn perceive_bonds_uses_the_fallback_radius_for_unknown_labels() {
        // 2.5 A is a bond for a metal-sized dummy atom next to a nitrogen ((1.5 + 0.71) * 1.3
        // = 2.87 A) but not between two nitrogens (1.85 A).
        let coords = vec![Vector3::zeros(), Vector3::new(2.5, 0.0, 0.0)];

        assert_eq!(
            perceive_bonds(&labels(&["Xx", "N"]), &coords, false),
            vec![(0, 1)]
        );
        assert!(perceive_bonds(&labels(&["N", "N"]), &coords, false).is_empty());
    }
}
