use crate::cshm::types::CShMResult;
use crate::csom::prepare::strip_all_labels;
use crate::csom::types::CsomResult;
use crate::data::elements::covalent_radius;
use crate::odis::OdisResult;
use nalgebra::{Matrix3, Vector3};
use std::fs::File;
use std::io::Write;

/// Escapes a CSV field per RFC 4180: wraps it in double quotes, doubling any embedded quotes,
/// whenever it contains a comma, a quote or a newline.
fn csv_field(s: &str) -> String {
    if s.contains(',') || s.contains('"') || s.contains('\n') || s.contains('\r') {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_string()
    }
}

fn format_matrix3(m: &Matrix3<f64>) -> String {
    format!(
        "[{:.4} {:.4} {:.4}; {:.4} {:.4} {:.4}; {:.4} {:.4} {:.4}]",
        m[(0, 0)],
        m[(0, 1)],
        m[(0, 2)],
        m[(1, 0)],
        m[(1, 1)],
        m[(1, 2)],
        m[(2, 0)],
        m[(2, 1)],
        m[(2, 2)],
    )
}

pub fn welcome_msg() {
    let msg: &str = {
        r"
  _  ______   _____ __  __  ____   _____ _    _ _      ____  _____
 | |/ / __ \ / ____|  \/  |/ __ \ / ____| |  | | |    / __ \|  __ \
 | ' / |  | | (___ | \  / | |  | | |    | |__| | |   | |  | | |__) |
 |  /  |  | |\___ \| |\/| | |  | | |    |  __  | |   | |  | |  _  /
 | . \ |__| |____) | |  | | |__| | |____| |  | | |___| |__| | | \ \
 |_|\_\____/|_____/|_|  |_|\____/ \_____|_|  |_|______\____/|_|  \_\
"
    };
    println!("{}", msg);
    println!("{}", env!("CARGO_PKG_DESCRIPTION"));
    println!("Version: {}", env!("CARGO_PKG_VERSION"));
    println!("Authors: {}", env!("CARGO_PKG_AUTHORS"));
    println!("Repository: {}", env!("CARGO_PKG_REPOSITORY"));
}

pub fn print_cshm_table(results: &[CShMResult], file: &str) {
    println!("\nInput file: {}", file);

    let name_width = results.iter().map(|r| r.name.len()).max().unwrap() + 2;
    let symbol_width = results.iter().map(|r| r.symbol.len()).max().unwrap() + 2;
    let symm_width = "Symmetry".len() + 2;
    let total_width = symbol_width + name_width + symm_width + 7 + 4;

    println!("{}", "=".repeat(total_width));
    println!(
        " {:<sw$} {:<nw$} {:<syw$} {:<7}",
        "Symbol",
        "Shape",
        "Symmetry",
        "CShM",
        sw = symbol_width,
        nw = name_width,
        syw = symm_width
    );
    println!("{}", "-".repeat(total_width));

    for result in results {
        println!(
            " {:<sw$} {:<nw$} {:<syw$} {:<7.3}",
            result.symbol,
            result.name,
            result.symm,
            result.cshm,
            sw = symbol_width,
            nw = name_width,
            syw = symm_width
        );
    }

    let min_s = results
        .iter()
        .map(|r| r.cshm)
        .min_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();
    println!("{}", "-".repeat(total_width));
    if min_s > 10.0 {
        println!(
            "Only extremely distorted geometries were found for this shape. Make sure the .xyz file is correct."
        )
    }
}

/// Writes the cshm results table (symbol, name, symmetry, s-value)
/// to a `<file>_cshm_table.csv` file.
pub fn write_cshm_csv(results: &[CShMResult], file_name: &str) -> Result<(), std::io::Error> {
    let out_name = file_name
        .strip_suffix(".xyz")
        .unwrap_or(file_name)
        .to_owned()
        + "_cshm_table.csv";
    let mut file = File::create(&out_name)?;

    println!("Writing output table to {}...", out_name);

    writeln!(file, "Symbol,Name,Symmetry,CShM").expect("Unable to write to file.");
    for r in results {
        writeln!(
            file,
            "{},{},{},{:.3}",
            csv_field(&r.symbol),
            csv_field(&r.name),
            csv_field(&r.symm),
            r.cshm
        )?;
    }

    Ok(())
}

/// Writes all ideal reference shapes scaled and rotated to align to the problem structure
/// to a `<file>_ideal.xyz` file.
pub fn write_cshm_reconstructed_xyz(
    file_name: &str,
    results: &[CShMResult],
    labels: &[String],
) -> Result<(), std::io::Error> {
    let out_name = file_name
        .strip_suffix(".xyz")
        .unwrap_or(file_name)
        .to_owned()
        + "_ideal.xyz";
    let mut file = File::create(&out_name)?;

    println!(
        "Writing idealised polyhedra coordinates to table to {}...",
        out_name
    );

    for result in results {
        // Write the preamble to each xyz block.
        writeln!(file, "{}", labels.len())?;
        writeln!(
            file,
            "{} {} CShM = {:.3}",
            result.symbol, result.symm, result.cshm
        )?;

        let mut inverse_perm = vec![0usize; result.perm.len()];
        for (problem_idx, &ref_idx) in result.perm.iter().enumerate() {
            inverse_perm[ref_idx] = problem_idx;
        }

        for (ref_idx, point) in result.xyz.iter().enumerate() {
            let problem_idx = inverse_perm[ref_idx];
            writeln!(
                file,
                "{}  {:.6}  {:.6}  {:.6}",
                labels[problem_idx], point.x, point.y, point.z
            )?;
        }
        writeln!(file)?;
    }

    Ok(())
}

pub fn print_crab() {
    println!(
        r"
     /\
    ( /   @ @    ()
     \\ __| |__  /
      \/   v   \/
     /-|       |-\
    / /-\     /-\ \
     / /-`---'-\ \
      /         \ "
    )
}

pub fn print_odis_table(result: &OdisResult, file: &str) {
    println!("\nInput file: {}", file);
    println!("{}", "=".repeat(34));
    println!(" Octahedral distortion parameters");
    println!("{}", "-".repeat(34));
    println!(
        "{:<16}{:>12.4}  {:<12}",
        " Mean d(M-X)", result.d_mean, "Ang"
    );
    println!("{:<16}{:>12.4}  {:<12}", " Zeta", result.zeta, "Ang");
    println!("{:<16}{:>12.6}  {:<12}", " Delta", result.delta, "");
    println!("{:<16}{:>12.2}  {:<12}", " Sigma", result.sigma, "deg");
    println!("{:<16}{:>12.2}  {:<12}", " Theta", result.theta, "deg");
    println!("{:<16}{:>12.4}  {:<12}", " Volume", result.vol, "Ang^3");
    println!("{}", "-".repeat(34));
    println!("{:<16}{:>12.2}  {:<12}", " Tau", result.tau, "deg");
    println!("{:<16}{:>12.2}  {:<12}", " Mu", result.mu, "Ang");
    println!("{}", "=".repeat(34));
}

pub fn print_csom_table(results: &[CsomResult], file: &str) {
    println!("\nInput file: {}", file);
    println!("{}", "=".repeat(20));
    println!(" {:<12} {:<10}", "Point group", "CSoM");
    println!("{}", "-".repeat(20));
    for result in results {
        println!(" {:<12} {:<10.3}", result.point_group, result.deviation);
    }
    println!("{}", "-".repeat(20));
}

pub fn write_odis_csv(result: OdisResult, file_name: &str) -> Result<(), std::io::Error> {
    let out_name = file_name
        .strip_suffix(".xyz")
        .unwrap_or(file_name)
        .to_owned()
        + "_odis_table.csv";
    let mut file = File::create(&out_name)?;

    println!("Writing output table to {}...", out_name);

    writeln!(file, "d_mean,zeta,delta,sigma,theta,vol,tau,mu")?;
    writeln!(
        file,
        "{:.4},{:.4},{:.6},{:.2},{:.2},{:.4},{:.4},{:.4}",
        result.d_mean,
        result.zeta,
        result.delta,
        result.sigma,
        result.theta,
        result.vol,
        result.tau,
        result.mu
    )?;

    Ok(())
}

/// Writes the csom summary table (point group, deviation, refined rotation matrix) to a
///  `<file>_csom_table.csv` file.
pub fn write_csom_csv(results: &[CsomResult], file_name: &str) -> Result<(), std::io::Error> {
    let out_name = file_name
        .strip_suffix(".xyz")
        .unwrap_or(file_name)
        .to_owned()
        + "_csom_table.csv";
    let mut file = File::create(&out_name)?;

    println!("Writing output table to {}...", out_name);

    writeln!(file, "PointGroup,Dev,Rotation Matrix")?;
    for r in results {
        writeln!(
            file,
            "{},{:.3},{}",
            csv_field(&r.point_group),
            r.deviation,
            format_matrix3(&r.rotation)
        )?;
    }

    Ok(())
}

/// Formats an operation's atom pairing as `Source>Target` entries in the input's atom order:
/// `A>B` means the operation carries atom `A` onto the site of atom `B`.
fn format_pairing(pairing: &[usize], labels: &[String]) -> String {
    // `pairing[target] = source`, so invert it to walk the atoms in input order.
    let mut target_of = vec![0usize; pairing.len()];
    for (target, &source) in pairing.iter().enumerate() {
        target_of[source] = target;
    }

    target_of
        .iter()
        .enumerate()
        .map(|(source, &target)| format!("{}>{}", labels[source], labels[target]))
        .collect::<Vec<String>>()
        .join(" ")
}

/// Writes the individual symmetry operation deviations (name, matrix, s-value, atom pairing)
/// per point group analysed to a  `<file>_<point group>_details.csv` file.
pub fn write_csom_details_csv(
    results: &[CsomResult],
    file_name: &str,
    labels: &[String],
) -> Result<(), std::io::Error> {
    let stem = file_name.strip_suffix(".xyz").unwrap_or(file_name);

    for result in results {
        let out_name = format!("{}_{}_details.csv", stem, result.point_group);
        let mut file = File::create(&out_name)?;

        println!("Writing operation details to {}...", out_name);

        writeln!(file, "name,op_matrix,dev,pairing")?;
        for op in &result.operations {
            writeln!(
                file,
                "{},{},{:.3},{}",
                csv_field(&op.name),
                format_matrix3(&op.matrix),
                op.deviation,
                csv_field(&format_pairing(&op.pairing, labels)),
            )?;
        }
    }

    Ok(())
}

/// Writes, for each point group, a single `<file>_<point group>_operated.xyz` with all
/// reconstructed operated structures by each symmetry element.
pub fn write_csom_operated_xyz(
    results: &[CsomResult],
    file_name: &str,
    labels: &[String],
    original_coords: &[Vector3<f64>],
) -> Result<(), std::io::Error> {
    let stem = file_name.strip_suffix(".xyz").unwrap_or(file_name);

    for result in results {
        let out_name = format!("{}_{}_operated.xyz", stem, result.point_group);
        let mut file = File::create(&out_name)?;

        println!("Writing operated coordinates to {}...", out_name);

        // The identity isn't stored among `result.operations` (the point-group tables in
        // data/pgs.rs omit E -- it trivially gives zero deviation for any structure), so write
        // the untouched original structure as its own "E" block first.
        writeln!(file, "{}", labels.len())?;
        writeln!(file, "{} E dev = 0.000", result.point_group)?;
        for (label, point) in labels.iter().zip(original_coords) {
            writeln!(
                file,
                "{}  {:.6}  {:.6}  {:.6}",
                label, point.x, point.y, point.z
            )?;
        }
        writeln!(file)?;

        // `rotation` is orthogonal, so its transpose is its inverse -- undoes the
        // CSOM-alignment rotation without an explicit matrix inversion.
        let rotation_inv = result.rotation.transpose();

        // `op.image[i]` is where the operation sends atom `i`, so every atom keeps its own
        // label: an atom whose image lands on a *different* atom's site shows up as such.
        for op in &result.operations {
            writeln!(file, "{}", labels.len())?;
            writeln!(
                file,
                "{} {} dev = {:.3}",
                result.point_group, op.name, op.deviation
            )?;
            for (label, point) in labels.iter().zip(&op.image) {
                let real_point = (rotation_inv * point) / result.scale + result.centroid;
                writeln!(
                    file,
                    "{}  {:.6}  {:.6}  {:.6}",
                    label, real_point.x, real_point.y, real_point.z
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

/// Bonds of one block of the merged .mol2, as 0-based atom index pairs `(i, j)` with `i < j`:
/// every pair within `BOND_TOLERANCE` times the sum of their covalent radii, plus, when the
/// structure has a centre atom (index 0), the centre to every ligand whatever its distance,
/// so the coordination polyhedron is always drawn.
///
/// Each block is a rigid image of the original structure, so one list serves all of them.
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

/// Writes, for each point group, a single `<file>_<point group>_merged.mol2` with all
/// the reconstructed operated structures. Bonds come from [`perceive_bonds`].
pub fn write_csom_merged_mol2(
    results: &[CsomResult],
    file_name: &str,
    labels: &[String],
    original_coords: &[Vector3<f64>],
    has_centre_atom: bool,
) -> Result<(), std::io::Error> {
    let stem = file_name.strip_suffix(".xyz").unwrap_or(file_name);
    let elements = strip_all_labels(labels);
    let n_per_block = labels.len();
    let bonds = perceive_bonds(&elements, original_coords, has_centre_atom);

    for result in results {
        let out_name = format!("{}_{}_merged.mol2", stem, result.point_group);
        let mut file = File::create(&out_name)?;

        println!("Writing merged mol2 to {}...", out_name);

        // Each block (the original "E" structure, then one per symmetry operation) as
        // (substructure name, that block's atom coordinates recovered to the original frame).
        let rotation_inv = result.rotation.transpose();
        let mut blocks: Vec<(&str, Vec<Vector3<f64>>)> = vec![("E", original_coords.to_vec())];
        blocks.extend(result.operations.iter().map(|op| {
            let real_points = op
                .image
                .iter()
                .map(|p| (rotation_inv * p) / result.scale + result.centroid)
                .collect();
            (op.name.as_str(), real_points)
        }));

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
        for (block_idx, (name, points)) in blocks.iter().enumerate() {
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
        for (block, (name, _)) in blocks.iter().enumerate() {
            writeln!(file, "{} {} {}", block + 1, name, block * n_per_block + 1)?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn csv_field_leaves_plain_text_untouched() {
        assert_eq!(csv_field("vOC-5"), "vOC-5");
    }

    #[test]
    fn csv_field_quotes_a_comma() {
        assert_eq!(csv_field("My Shape, Custom"), "\"My Shape, Custom\"");
    }

    #[test]
    fn csv_field_doubles_embedded_quotes() {
        assert_eq!(csv_field("6\" wide"), "\"6\"\" wide\"");
    }

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
