// Symmetry operations of every supported point group, principal axis along z. The identity
// is left out of each table (except E's own), since every operation set implicitly has it.
//
// Each matrix is stored row-major and acts on column vectors (image = M * p):
// [[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]]
//
// Operations are labelled by conjugacy class, as in standard character tables: the prefix is the
// size of the class (3C2'), and C_n is the counterclockwise rotation by 360/n degrees about +z.
//
// Every matrix must be exact to full f64 precision. Run check_pointgroup_tables.py after
// editing this file: it checks that every table is the group its name says.

use nalgebra::Matrix3;
use std::collections::HashMap;
use std::f64::consts::FRAC_1_SQRT_2;

// std::f64::consts has no sqrt(3)/2 or 36°/72° trig constants, so these are
// spelled out to full f64 precision.
const SQRT_3_DIV_2: f64 = 0.866_025_403_784_438_6;
const COS_36: f64 = 0.809_016_994_374_947_4;
const SIN_36: f64 = 0.587_785_252_292_473_1;
const COS_72: f64 = 0.309_016_994_374_947_4;
const SIN_72: f64 = 0.951_056_516_295_153_6;

// The seven-fold groups (C7, D7h) are built from multiples of pi/7.
const COS_PI_7: f64 = 0.900_968_867_902_419_1;
const COS_2PI_7: f64 = 0.623_489_801_858_733_6;
const COS_3PI_7: f64 = 0.222_520_933_956_314_45;
const SIN_PI_7: f64 = 0.433_883_739_117_558_1;
const SIN_2PI_7: f64 = 0.781_831_482_468_029_8;
const SIN_3PI_7: f64 = 0.974_927_912_181_823_6;

// The icosahedral groups (I, Ih) have z along a C5 axis and a second C5 axis in the yz plane,
// tilted arccos(1/sqrt(5)) from z, so their entries combine 1/sqrt(5) with the 36° and 72°
// constants.
const FRAC_1_SQRT_5: f64 = 0.447_213_595_499_957_9;
const FRAC_2_SQRT_5: f64 = 2.0 * FRAC_1_SQRT_5;
const FRAC_3_SQRT_20: f64 = 1.5 * FRAC_1_SQRT_5;
const HALF_PLUS_FRAC_1_SQRT_5: f64 = 0.5 + FRAC_1_SQRT_5;
const HALF_MINUS_FRAC_1_SQRT_5: f64 = 0.5 - FRAC_1_SQRT_5;
const COS_36_DIV_SQRT_5: f64 = COS_36 * FRAC_1_SQRT_5;
const COS_72_DIV_SQRT_5: f64 = COS_72 * FRAC_1_SQRT_5;
const SIN_36_DIV_SQRT_5: f64 = SIN_36 * FRAC_1_SQRT_5;
const SIN_72_DIV_SQRT_5: f64 = SIN_72 * FRAC_1_SQRT_5;
const TWO_COS_36_DIV_SQRT_5: f64 = 2.0 * COS_36_DIV_SQRT_5;
const TWO_COS_72_DIV_SQRT_5: f64 = 2.0 * COS_72_DIV_SQRT_5;
const TWO_SIN_36_DIV_SQRT_5: f64 = 2.0 * SIN_36_DIV_SQRT_5;
const TWO_SIN_72_DIV_SQRT_5: f64 = 2.0 * SIN_72_DIV_SQRT_5;
const HALF_PLUS_COS_36_DIV_SQRT_5: f64 = 0.5 + COS_36_DIV_SQRT_5;
const HALF_PLUS_COS_72_DIV_SQRT_5: f64 = 0.5 + COS_72_DIV_SQRT_5;
const SIN_72_PLUS_SIN_36_DIV_SQRT_5: f64 = (SIN_72 + SIN_36) * FRAC_1_SQRT_5;
const SIN_72_MINUS_SIN_36_DIV_SQRT_5: f64 = (SIN_72 - SIN_36) * FRAC_1_SQRT_5;

/// A single symmetry operation: its Schoenflies name paired with its row-major matrix.
pub type SymmetryOperation = (&'static str, [[f64; 3]; 3]);

/// Look up a point group's raw symmetry-operation list.
pub fn get_pointgroup(name: &str) -> Option<&'static [SymmetryOperation]> {
    match name {
        "C2" => Some(POINTGROUP_C2),
        "C2h" => Some(POINTGROUP_C2H),
        "C2v" => Some(POINTGROUP_C2V),
        "C3" => Some(POINTGROUP_C3),
        "C3h" => Some(POINTGROUP_C3H),
        "C3v" => Some(POINTGROUP_C3V),
        "C4" => Some(POINTGROUP_C4),
        "C4h" => Some(POINTGROUP_C4H),
        "C4v" => Some(POINTGROUP_C4V),
        "C5" => Some(POINTGROUP_C5),
        "C5h" => Some(POINTGROUP_C5H),
        "C5v" => Some(POINTGROUP_C5V),
        "C6" => Some(POINTGROUP_C6),
        "C6h" => Some(POINTGROUP_C6H),
        "C6v" => Some(POINTGROUP_C6V),
        "C7" => Some(POINTGROUP_C7),
        "C8" => Some(POINTGROUP_C8),
        "Ci" => Some(POINTGROUP_CI),
        "Cs" => Some(POINTGROUP_CS),
        "D2" => Some(POINTGROUP_D2),
        "D2d" => Some(POINTGROUP_D2D),
        "D2h" => Some(POINTGROUP_D2H),
        "D3" => Some(POINTGROUP_D3),
        "D3d" => Some(POINTGROUP_D3D),
        "D3h" => Some(POINTGROUP_D3H),
        "D4" => Some(POINTGROUP_D4),
        "D4d" => Some(POINTGROUP_D4D),
        "D4h" => Some(POINTGROUP_D4H),
        "D5" => Some(POINTGROUP_D5),
        "D5d" => Some(POINTGROUP_D5D),
        "D5h" => Some(POINTGROUP_D5H),
        "D6" => Some(POINTGROUP_D6),
        "D6h" => Some(POINTGROUP_D6H),
        "D7h" => Some(POINTGROUP_D7H),
        "D8h" => Some(POINTGROUP_D8H),
        "E" => Some(POINTGROUP_E),
        "I" => Some(POINTGROUP_I),
        "Ih" => Some(POINTGROUP_IH),
        "O" => Some(POINTGROUP_O),
        "Oh" => Some(POINTGROUP_OH),
        "S10" => Some(POINTGROUP_S10),
        "S4" => Some(POINTGROUP_S4),
        "S6" => Some(POINTGROUP_S6),
        "S8" => Some(POINTGROUP_S8),
        "T" => Some(POINTGROUP_T),
        "Td" => Some(POINTGROUP_TD),
        "Th" => Some(POINTGROUP_TH),
        _ => None,
    }
}

/// Look up a point group as a HashMap<operation name, Vec<Matrix3<f64>>>.
/// An operation name (e.g. "3C2") can repeat with a different matrix
/// each time; all matrices sharing a name are collected under that key.
pub fn get_pointgroup_map(name: &str) -> Option<HashMap<&'static str, Vec<Matrix3<f64>>>> {
    get_pointgroup(name).map(|ops| {
        let mut map: HashMap<&'static str, Vec<Matrix3<f64>>> = HashMap::new();
        for (op_name, matrix) in ops {
            map.entry(*op_name).or_default().push(to_matrix3(*matrix));
        }
        map
    })
}

#[rustfmt::skip]
pub static POINTGROUP_C2: &[(&str, [[f64; 3]; 3])] = &[
    ("C2", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_C2H: &[(&str, [[f64; 3]; 3])] = &[
    ("C2", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("i", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("sigma_h", [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_C2V: &[(&str, [[f64; 3]; 3])] = &[
    ("C2", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("sigma_v(xz)", [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("sigma_v(yz)", [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_C3: &[(&str, [[f64; 3]; 3])] = &[
    ("C3^2", [[-0.5, SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("C3", [[-0.5, -SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, 1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_C3H: &[(&str, [[f64; 3]; 3])] = &[
    ("C3^2", [[-0.5, SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("C3", [[-0.5, -SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("S3^5", [[-0.5, SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, -1.0]]),
    ("S3", [[-0.5, -SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, -1.0]]),
    ("sigma_h", [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_C3V: &[(&str, [[f64; 3]; 3])] = &[
    ("2C3", [[-0.5, SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("2C3", [[-0.5, -SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("3sigma_v", [[-0.5, SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, 0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("3sigma_v", [[-0.5, -SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, 0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("3sigma_v", [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_C4: &[(&str, [[f64; 3]; 3])] = &[
    ("C4^3", [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
    ("C2", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("C4", [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_C4H: &[(&str, [[f64; 3]; 3])] = &[
    ("C4^3", [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
    ("C4", [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
    ("C2", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("i", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("S4^3", [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]),
    ("S4", [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]),
    ("sigma_h", [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_C4V: &[(&str, [[f64; 3]; 3])] = &[
    ("2C4", [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
    ("2C4", [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
    ("C2", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("2sigma_v", [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("2sigma_v", [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("2sigma_d", [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
    ("2sigma_d", [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_C5: &[(&str, [[f64; 3]; 3])] = &[
    ("C5^4", [[COS_72, SIN_72, 0.0], [-SIN_72, COS_72, 0.0], [0.0, 0.0, 1.0]]),
    ("C5", [[COS_72, -SIN_72, 0.0], [SIN_72, COS_72, 0.0], [0.0, 0.0, 1.0]]),
    ("C5^3", [[-COS_36, SIN_36, 0.0], [-SIN_36, -COS_36, 0.0], [0.0, 0.0, 1.0]]),
    ("C5^2", [[-COS_36, -SIN_36, 0.0], [SIN_36, -COS_36, 0.0], [0.0, 0.0, 1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_C5H: &[(&str, [[f64; 3]; 3])] = &[
    ("C5^4", [[COS_72, SIN_72, 0.0], [-SIN_72, COS_72, 0.0], [0.0, 0.0, 1.0]]),
    ("C5^3", [[-COS_36, SIN_36, 0.0], [-SIN_36, -COS_36, 0.0], [0.0, 0.0, 1.0]]),
    ("C5^2", [[-COS_36, -SIN_36, 0.0], [SIN_36, -COS_36, 0.0], [0.0, 0.0, 1.0]]),
    ("C5", [[COS_72, -SIN_72, 0.0], [SIN_72, COS_72, 0.0], [0.0, 0.0, 1.0]]),
    ("sigma_h", [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("S5^9", [[COS_72, SIN_72, 0.0], [-SIN_72, COS_72, 0.0], [0.0, 0.0, -1.0]]),
    ("S5^7", [[-COS_36, -SIN_36, 0.0], [SIN_36, -COS_36, 0.0], [0.0, 0.0, -1.0]]),
    ("S5^3", [[-COS_36, SIN_36, 0.0], [-SIN_36, -COS_36, 0.0], [0.0, 0.0, -1.0]]),
    ("S5", [[COS_72, -SIN_72, 0.0], [SIN_72, COS_72, 0.0], [0.0, 0.0, -1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_C5V: &[(&str, [[f64; 3]; 3])] = &[
    ("2C5", [[COS_72, SIN_72, 0.0], [-SIN_72, COS_72, 0.0], [0.0, 0.0, 1.0]]),
    ("2C5", [[COS_72, -SIN_72, 0.0], [SIN_72, COS_72, 0.0], [0.0, 0.0, 1.0]]),
    ("2C5^2", [[-COS_36, SIN_36, 0.0], [-SIN_36, -COS_36, 0.0], [0.0, 0.0, 1.0]]),
    ("2C5^2", [[-COS_36, -SIN_36, 0.0], [SIN_36, -COS_36, 0.0], [0.0, 0.0, 1.0]]),
    ("5sigma_v", [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("5sigma_v", [[COS_36, -SIN_36, 0.0], [-SIN_36, -COS_36, 0.0], [0.0, 0.0, 1.0]]),
    ("5sigma_v", [[-COS_72, SIN_72, 0.0], [SIN_72, COS_72, 0.0], [0.0, 0.0, 1.0]]),
    ("5sigma_v", [[-COS_72, -SIN_72, 0.0], [-SIN_72, COS_72, 0.0], [0.0, 0.0, 1.0]]),
    ("5sigma_v", [[COS_36, SIN_36, 0.0], [SIN_36, -COS_36, 0.0], [0.0, 0.0, 1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_C6: &[(&str, [[f64; 3]; 3])] = &[
    ("C6^5", [[0.5, SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, 0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("C6", [[0.5, -SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, 0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("C3^2", [[-0.5, SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("C3", [[-0.5, -SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("C2", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_C6H: &[(&str, [[f64; 3]; 3])] = &[
    ("C6", [[0.5, -SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, 0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("C6^5", [[0.5, SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, 0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("C3", [[-0.5, -SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("C3^2", [[-0.5, SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("C2", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("i", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("S3", [[-0.5, -SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, -1.0]]),
    ("S3^5", [[-0.5, SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, -1.0]]),
    ("S6", [[0.5, -SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, 0.5, 0.0], [0.0, 0.0, -1.0]]),
    ("S6^5", [[0.5, SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, 0.5, 0.0], [0.0, 0.0, -1.0]]),
    ("sigma_h", [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_C6V: &[(&str, [[f64; 3]; 3])] = &[
    ("2C6", [[0.5, SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, 0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("2C6", [[0.5, -SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, 0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("2C3", [[-0.5, SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("2C3", [[-0.5, -SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("C2", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("3sigma_d", [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("3sigma_d", [[0.5, -SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("3sigma_d", [[0.5, SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("3sigma_v", [[-0.5, -SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, 0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("3sigma_v", [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("3sigma_v", [[-0.5, SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, 0.5, 0.0], [0.0, 0.0, 1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_C7: &[(&str, [[f64; 3]; 3])] = &[
    ("C7^6", [[COS_2PI_7, SIN_2PI_7, 0.0], [-SIN_2PI_7, COS_2PI_7, 0.0], [0.0, 0.0, 1.0]]),
    ("C7^5", [[-COS_3PI_7, SIN_3PI_7, 0.0], [-SIN_3PI_7, -COS_3PI_7, 0.0], [0.0, 0.0, 1.0]]),
    ("C7^4", [[-COS_PI_7, SIN_PI_7, 0.0], [-SIN_PI_7, -COS_PI_7, 0.0], [0.0, 0.0, 1.0]]),
    ("C7^3", [[-COS_PI_7, -SIN_PI_7, 0.0], [SIN_PI_7, -COS_PI_7, 0.0], [0.0, 0.0, 1.0]]),
    ("C7^2", [[-COS_3PI_7, -SIN_3PI_7, 0.0], [SIN_3PI_7, -COS_3PI_7, 0.0], [0.0, 0.0, 1.0]]),
    ("C7", [[COS_2PI_7, -SIN_2PI_7, 0.0], [SIN_2PI_7, COS_2PI_7, 0.0], [0.0, 0.0, 1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_C8: &[(&str, [[f64; 3]; 3])] = &[
    ("C8^7", [[FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0], [-FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0], [0.0, 0.0, 1.0]]),
    ("C4^3", [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
    ("C8^5", [[-FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0], [-FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0], [0.0, 0.0, 1.0]]),
    ("C2", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("C8^3", [[-FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0], [FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0], [0.0, 0.0, 1.0]]),
    ("C4", [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
    ("C8", [[FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0], [FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0], [0.0, 0.0, 1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_CI: &[(&str, [[f64; 3]; 3])] = &[
    ("i", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_CS: &[(&str, [[f64; 3]; 3])] = &[
    ("sigma_h", [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_D2: &[(&str, [[f64; 3]; 3])] = &[
    ("C2(z)", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("C2(y)", [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("C2(x)", [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_D2D: &[(&str, [[f64; 3]; 3])] = &[
    ("2S4", [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]),
    ("2S4", [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]),
    ("C2", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("2C2'", [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]),
    ("2C2'", [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]),
    ("2sigma_d", [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("2sigma_d", [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_D2H: &[(&str, [[f64; 3]; 3])] = &[
    ("C2(x)", [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("C2(y)", [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("C2(z)", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("i", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("sigma(yz)", [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("sigma(xz)", [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("sigma(xy)", [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_D3: &[(&str, [[f64; 3]; 3])] = &[
    ("2C3", [[-0.5, SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("2C3", [[-0.5, -SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("3C2", [[0.5, -SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, -1.0]]),
    ("3C2", [[0.5, SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, -1.0]]),
    ("3C2", [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_D3D: &[(&str, [[f64; 3]; 3])] = &[
    ("2C3", [[-0.5, SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("2C3", [[-0.5, -SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("3C2", [[0.5, -SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, -1.0]]),
    ("3C2", [[0.5, SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, -1.0]]),
    ("3C2", [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("i", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("2S6", [[0.5, SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, 0.5, 0.0], [0.0, 0.0, -1.0]]),
    ("2S6", [[0.5, -SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, 0.5, 0.0], [0.0, 0.0, -1.0]]),
    ("3sigma_d", [[-0.5, SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, 0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("3sigma_d", [[-0.5, -SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, 0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("3sigma_d", [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_D3H: &[(&str, [[f64; 3]; 3])] = &[
    ("2C3", [[-0.5, SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("2C3", [[-0.5, -SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("3C2", [[0.5, -SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, -1.0]]),
    ("3C2", [[0.5, SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, -1.0]]),
    ("3C2", [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("sigma_h", [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("2S3", [[-0.5, SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, -1.0]]),
    ("2S3", [[-0.5, -SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, -1.0]]),
    ("3sigma_v", [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("3sigma_v", [[0.5, -SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("3sigma_v", [[0.5, SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, 1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_D4: &[(&str, [[f64; 3]; 3])] = &[
    ("2C4", [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
    ("2C4", [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
    ("C2", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("2C2'", [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("2C2'", [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("2C2''", [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]),
    ("2C2''", [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_D4D: &[(&str, [[f64; 3]; 3])] = &[
    ("2S8", [[FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0], [-FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0], [0.0, 0.0, -1.0]]),
    ("2S8", [[FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0], [FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0], [0.0, 0.0, -1.0]]),
    ("2C4", [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
    ("2C4", [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
    ("2S8^3", [[-FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0], [-FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0], [0.0, 0.0, -1.0]]),
    ("2S8^3", [[-FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0], [FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0], [0.0, 0.0, -1.0]]),
    ("C2", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("4C2'", [[FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0], [FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0], [0.0, 0.0, -1.0]]),
    ("4C2'", [[-FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0], [FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0], [0.0, 0.0, -1.0]]),
    ("4C2'", [[-FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0], [-FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0], [0.0, 0.0, -1.0]]),
    ("4C2'", [[FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0], [-FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0], [0.0, 0.0, -1.0]]),
    ("4sigma_d", [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("4sigma_d", [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("4sigma_d", [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
    ("4sigma_d", [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_D4H: &[(&str, [[f64; 3]; 3])] = &[
    ("2C4", [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
    ("2C4", [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
    ("C2", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("2C2'", [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("2C2'", [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("2C2''", [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]),
    ("2C2''", [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]),
    ("i", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("2S4", [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]),
    ("2S4", [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]),
    ("sigma_h", [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("2sigma_v", [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("2sigma_v", [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("2sigma_d", [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
    ("2sigma_d", [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_D5: &[(&str, [[f64; 3]; 3])] = &[
    ("2C5", [[COS_72, SIN_72, 0.0], [-SIN_72, COS_72, 0.0], [0.0, 0.0, 1.0]]),
    ("2C5", [[COS_72, -SIN_72, 0.0], [SIN_72, COS_72, 0.0], [0.0, 0.0, 1.0]]),
    ("2C5^2", [[-COS_36, SIN_36, 0.0], [-SIN_36, -COS_36, 0.0], [0.0, 0.0, 1.0]]),
    ("2C5^2", [[-COS_36, -SIN_36, 0.0], [SIN_36, -COS_36, 0.0], [0.0, 0.0, 1.0]]),
    ("5C2", [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("5C2", [[COS_36, SIN_36, 0.0], [SIN_36, -COS_36, 0.0], [0.0, 0.0, -1.0]]),
    ("5C2", [[-COS_72, -SIN_72, 0.0], [-SIN_72, COS_72, 0.0], [0.0, 0.0, -1.0]]),
    ("5C2", [[-COS_72, SIN_72, 0.0], [SIN_72, COS_72, 0.0], [0.0, 0.0, -1.0]]),
    ("5C2", [[COS_36, -SIN_36, 0.0], [-SIN_36, -COS_36, 0.0], [0.0, 0.0, -1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_D5D: &[(&str, [[f64; 3]; 3])] = &[
    ("2C5", [[COS_72, SIN_72, 0.0], [-SIN_72, COS_72, 0.0], [0.0, 0.0, 1.0]]),
    ("2C5", [[COS_72, -SIN_72, 0.0], [SIN_72, COS_72, 0.0], [0.0, 0.0, 1.0]]),
    ("2C5^2", [[-COS_36, SIN_36, 0.0], [-SIN_36, -COS_36, 0.0], [0.0, 0.0, 1.0]]),
    ("2C5^2", [[-COS_36, -SIN_36, 0.0], [SIN_36, -COS_36, 0.0], [0.0, 0.0, 1.0]]),
    ("5C2", [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("5C2", [[COS_36, SIN_36, 0.0], [SIN_36, -COS_36, 0.0], [0.0, 0.0, -1.0]]),
    ("5C2", [[-COS_72, -SIN_72, 0.0], [-SIN_72, COS_72, 0.0], [0.0, 0.0, -1.0]]),
    ("5C2", [[-COS_72, SIN_72, 0.0], [SIN_72, COS_72, 0.0], [0.0, 0.0, -1.0]]),
    ("5C2", [[COS_36, -SIN_36, 0.0], [-SIN_36, -COS_36, 0.0], [0.0, 0.0, -1.0]]),
    ("i", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("2S10^3", [[-COS_72, SIN_72, 0.0], [-SIN_72, -COS_72, 0.0], [0.0, 0.0, -1.0]]),
    ("2S10^3", [[-COS_72, -SIN_72, 0.0], [SIN_72, -COS_72, 0.0], [0.0, 0.0, -1.0]]),
    ("2S10", [[COS_36, -SIN_36, 0.0], [SIN_36, COS_36, 0.0], [0.0, 0.0, -1.0]]),
    ("2S10", [[COS_36, SIN_36, 0.0], [-SIN_36, COS_36, 0.0], [0.0, 0.0, -1.0]]),
    ("5sigma_d", [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("5sigma_d", [[-COS_36, -SIN_36, 0.0], [-SIN_36, COS_36, 0.0], [0.0, 0.0, 1.0]]),
    ("5sigma_d", [[COS_72, SIN_72, 0.0], [SIN_72, -COS_72, 0.0], [0.0, 0.0, 1.0]]),
    ("5sigma_d", [[COS_72, -SIN_72, 0.0], [-SIN_72, -COS_72, 0.0], [0.0, 0.0, 1.0]]),
    ("5sigma_d", [[-COS_36, SIN_36, 0.0], [SIN_36, COS_36, 0.0], [0.0, 0.0, 1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_D5H: &[(&str, [[f64; 3]; 3])] = &[
    ("2C5", [[COS_72, SIN_72, 0.0], [-SIN_72, COS_72, 0.0], [0.0, 0.0, 1.0]]),
    ("2C5", [[COS_72, -SIN_72, 0.0], [SIN_72, COS_72, 0.0], [0.0, 0.0, 1.0]]),
    ("2C5^2", [[-COS_36, SIN_36, 0.0], [-SIN_36, -COS_36, 0.0], [0.0, 0.0, 1.0]]),
    ("2C5^2", [[-COS_36, -SIN_36, 0.0], [SIN_36, -COS_36, 0.0], [0.0, 0.0, 1.0]]),
    ("5C2", [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("5C2", [[-COS_36, SIN_36, 0.0], [SIN_36, COS_36, 0.0], [0.0, 0.0, -1.0]]),
    ("5C2", [[COS_72, -SIN_72, 0.0], [-SIN_72, -COS_72, 0.0], [0.0, 0.0, -1.0]]),
    ("5C2", [[COS_72, SIN_72, 0.0], [SIN_72, -COS_72, 0.0], [0.0, 0.0, -1.0]]),
    ("5C2", [[-COS_36, -SIN_36, 0.0], [-SIN_36, COS_36, 0.0], [0.0, 0.0, -1.0]]),
    ("sigma_h", [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("2S5", [[COS_72, SIN_72, 0.0], [-SIN_72, COS_72, 0.0], [0.0, 0.0, -1.0]]),
    ("2S5", [[COS_72, -SIN_72, 0.0], [SIN_72, COS_72, 0.0], [0.0, 0.0, -1.0]]),
    ("2S5^3", [[-COS_36, SIN_36, 0.0], [-SIN_36, -COS_36, 0.0], [0.0, 0.0, -1.0]]),
    ("2S5^3", [[-COS_36, -SIN_36, 0.0], [SIN_36, -COS_36, 0.0], [0.0, 0.0, -1.0]]),
    ("5sigma_v", [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("5sigma_v", [[-COS_36, -SIN_36, 0.0], [-SIN_36, COS_36, 0.0], [0.0, 0.0, 1.0]]),
    ("5sigma_v", [[COS_72, SIN_72, 0.0], [SIN_72, -COS_72, 0.0], [0.0, 0.0, 1.0]]),
    ("5sigma_v", [[COS_72, -SIN_72, 0.0], [-SIN_72, -COS_72, 0.0], [0.0, 0.0, 1.0]]),
    ("5sigma_v", [[-COS_36, SIN_36, 0.0], [SIN_36, COS_36, 0.0], [0.0, 0.0, 1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_D6: &[(&str, [[f64; 3]; 3])] = &[
    ("2C6", [[0.5, -SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, 0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("2C6", [[0.5, SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, 0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("2C3", [[-0.5, -SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("2C3", [[-0.5, SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("C2", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("3C2'", [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("3C2'", [[-0.5, SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, 0.5, 0.0], [0.0, 0.0, -1.0]]),
    ("3C2'", [[-0.5, -SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, 0.5, 0.0], [0.0, 0.0, -1.0]]),
    ("3C2''", [[0.5, SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, -1.0]]),
    ("3C2''", [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("3C2''", [[0.5, -SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, -1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_D6H: &[(&str, [[f64; 3]; 3])] = &[
    ("2C6", [[0.5, -SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, 0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("2C6", [[0.5, SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, 0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("2C3", [[-0.5, -SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("2C3", [[-0.5, SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("C2", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("3C2'", [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("3C2'", [[-0.5, SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, 0.5, 0.0], [0.0, 0.0, -1.0]]),
    ("3C2'", [[-0.5, -SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, 0.5, 0.0], [0.0, 0.0, -1.0]]),
    ("3C2''", [[0.5, SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, -1.0]]),
    ("3C2''", [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("3C2''", [[0.5, -SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, -1.0]]),
    ("i", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("2S3", [[-0.5, -SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, -1.0]]),
    ("2S3", [[-0.5, SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, -1.0]]),
    ("2S6", [[0.5, -SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, 0.5, 0.0], [0.0, 0.0, -1.0]]),
    ("2S6", [[0.5, SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, 0.5, 0.0], [0.0, 0.0, -1.0]]),
    ("sigma_h", [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("3sigma_v", [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("3sigma_v", [[-0.5, SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, 0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("3sigma_v", [[-0.5, -SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, 0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("3sigma_d", [[0.5, SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("3sigma_d", [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("3sigma_d", [[0.5, -SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, 1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_D7H: &[(&str, [[f64; 3]; 3])] = &[
    ("2C7", [[COS_2PI_7, SIN_2PI_7, 0.0], [-SIN_2PI_7, COS_2PI_7, 0.0], [0.0, 0.0, 1.0]]),
    ("2C7", [[COS_2PI_7, -SIN_2PI_7, 0.0], [SIN_2PI_7, COS_2PI_7, 0.0], [0.0, 0.0, 1.0]]),
    ("2C7^2", [[-COS_3PI_7, SIN_3PI_7, 0.0], [-SIN_3PI_7, -COS_3PI_7, 0.0], [0.0, 0.0, 1.0]]),
    ("2C7^2", [[-COS_3PI_7, -SIN_3PI_7, 0.0], [SIN_3PI_7, -COS_3PI_7, 0.0], [0.0, 0.0, 1.0]]),
    ("2C7^3", [[-COS_PI_7, SIN_PI_7, 0.0], [-SIN_PI_7, -COS_PI_7, 0.0], [0.0, 0.0, 1.0]]),
    ("2C7^3", [[-COS_PI_7, -SIN_PI_7, 0.0], [SIN_PI_7, -COS_PI_7, 0.0], [0.0, 0.0, 1.0]]),
    ("7C2", [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("7C2", [[-COS_3PI_7, SIN_3PI_7, 0.0], [SIN_3PI_7, COS_3PI_7, 0.0], [0.0, 0.0, -1.0]]),
    ("7C2", [[-COS_PI_7, -SIN_PI_7, 0.0], [-SIN_PI_7, COS_PI_7, 0.0], [0.0, 0.0, -1.0]]),
    ("7C2", [[COS_2PI_7, -SIN_2PI_7, 0.0], [-SIN_2PI_7, -COS_2PI_7, 0.0], [0.0, 0.0, -1.0]]),
    ("7C2", [[COS_2PI_7, SIN_2PI_7, 0.0], [SIN_2PI_7, -COS_2PI_7, 0.0], [0.0, 0.0, -1.0]]),
    ("7C2", [[-COS_PI_7, SIN_PI_7, 0.0], [SIN_PI_7, COS_PI_7, 0.0], [0.0, 0.0, -1.0]]),
    ("7C2", [[-COS_3PI_7, -SIN_3PI_7, 0.0], [-SIN_3PI_7, COS_3PI_7, 0.0], [0.0, 0.0, -1.0]]),
    ("sigma_h", [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("2S7", [[COS_2PI_7, SIN_2PI_7, 0.0], [-SIN_2PI_7, COS_2PI_7, 0.0], [0.0, 0.0, -1.0]]),
    ("2S7", [[COS_2PI_7, -SIN_2PI_7, 0.0], [SIN_2PI_7, COS_2PI_7, 0.0], [0.0, 0.0, -1.0]]),
    ("2S7^3", [[-COS_PI_7, SIN_PI_7, 0.0], [-SIN_PI_7, -COS_PI_7, 0.0], [0.0, 0.0, -1.0]]),
    ("2S7^3", [[-COS_PI_7, -SIN_PI_7, 0.0], [SIN_PI_7, -COS_PI_7, 0.0], [0.0, 0.0, -1.0]]),
    ("2S7^5", [[-COS_3PI_7, -SIN_3PI_7, 0.0], [SIN_3PI_7, -COS_3PI_7, 0.0], [0.0, 0.0, -1.0]]),
    ("2S7^5", [[-COS_3PI_7, SIN_3PI_7, 0.0], [-SIN_3PI_7, -COS_3PI_7, 0.0], [0.0, 0.0, -1.0]]),
    ("7sigma_v", [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("7sigma_v", [[-COS_3PI_7, -SIN_3PI_7, 0.0], [-SIN_3PI_7, COS_3PI_7, 0.0], [0.0, 0.0, 1.0]]),
    ("7sigma_v", [[-COS_PI_7, SIN_PI_7, 0.0], [SIN_PI_7, COS_PI_7, 0.0], [0.0, 0.0, 1.0]]),
    ("7sigma_v", [[COS_2PI_7, SIN_2PI_7, 0.0], [SIN_2PI_7, -COS_2PI_7, 0.0], [0.0, 0.0, 1.0]]),
    ("7sigma_v", [[COS_2PI_7, -SIN_2PI_7, 0.0], [-SIN_2PI_7, -COS_2PI_7, 0.0], [0.0, 0.0, 1.0]]),
    ("7sigma_v", [[-COS_PI_7, -SIN_PI_7, 0.0], [-SIN_PI_7, COS_PI_7, 0.0], [0.0, 0.0, 1.0]]),
    ("7sigma_v", [[-COS_3PI_7, SIN_3PI_7, 0.0], [SIN_3PI_7, COS_3PI_7, 0.0], [0.0, 0.0, 1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_D8H: &[(&str, [[f64; 3]; 3])] = &[
    ("2C8", [[FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0], [-FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0], [0.0, 0.0, 1.0]]),
    ("2C8", [[FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0], [FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0], [0.0, 0.0, 1.0]]),
    ("2C8^3", [[-FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0], [-FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0], [0.0, 0.0, 1.0]]),
    ("2C8^3", [[-FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0], [FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0], [0.0, 0.0, 1.0]]),
    ("2C4", [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
    ("2C4", [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
    ("C2", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("4C2'", [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("4C2'", [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]),
    ("4C2'", [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("4C2'", [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]),
    ("4C2''", [[-FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0], [FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0], [0.0, 0.0, -1.0]]),
    ("4C2''", [[FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0], [FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0], [0.0, 0.0, -1.0]]),
    ("4C2''", [[FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0], [-FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0], [0.0, 0.0, -1.0]]),
    ("4C2''", [[-FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0], [-FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0], [0.0, 0.0, -1.0]]),
    ("i", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("2S8", [[FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0], [-FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0], [0.0, 0.0, -1.0]]),
    ("2S8", [[FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0], [FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0], [0.0, 0.0, -1.0]]),
    ("2S8^3", [[-FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0], [-FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0], [0.0, 0.0, -1.0]]),
    ("2S8^3", [[-FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0], [FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0], [0.0, 0.0, -1.0]]),
    ("2S4", [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]),
    ("2S4", [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]),
    ("sigma_h", [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("4sigma_v", [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("4sigma_v", [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
    ("4sigma_v", [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("4sigma_v", [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
    ("4sigma_d", [[FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0], [-FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0], [0.0, 0.0, 1.0]]),
    ("4sigma_d", [[-FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0], [-FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0], [0.0, 0.0, 1.0]]),
    ("4sigma_d", [[-FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0], [FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0], [0.0, 0.0, 1.0]]),
    ("4sigma_d", [[FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0], [FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0], [0.0, 0.0, 1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_E: &[(&str, [[f64; 3]; 3])] = &[
    ("E", [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_I: &[(&str, [[f64; 3]; 3])] = &[
    ("12C5", [[COS_72, SIN_72, 0.0], [-SIN_72, COS_72, 0.0], [0.0, 0.0, 1.0]]),
    ("12C5", [[COS_72, -SIN_72_DIV_SQRT_5, -TWO_SIN_72_DIV_SQRT_5], [SIN_72_DIV_SQRT_5, HALF_PLUS_COS_36_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5], [TWO_SIN_72_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("12C5", [[COS_72, SIN_72_DIV_SQRT_5, TWO_SIN_72_DIV_SQRT_5], [-SIN_72_DIV_SQRT_5, HALF_PLUS_COS_36_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5], [-TWO_SIN_72_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("12C5", [[COS_72, -SIN_72, 0.0], [SIN_72, COS_72, 0.0], [0.0, 0.0, 1.0]]),
    ("12C5", [[0.5, SIN_72_PLUS_SIN_36_DIV_SQRT_5, -TWO_SIN_36_DIV_SQRT_5], [-SIN_72_MINUS_SIN_36_DIV_SQRT_5, FRAC_3_SQRT_20, TWO_COS_36_DIV_SQRT_5], [TWO_SIN_72_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("12C5", [[0.5, -SIN_72_PLUS_SIN_36_DIV_SQRT_5, TWO_SIN_36_DIV_SQRT_5], [SIN_72_MINUS_SIN_36_DIV_SQRT_5, FRAC_3_SQRT_20, TWO_COS_36_DIV_SQRT_5], [-TWO_SIN_72_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("12C5", [[0.5, SIN_72_MINUS_SIN_36_DIV_SQRT_5, -TWO_SIN_72_DIV_SQRT_5], [-SIN_72_PLUS_SIN_36_DIV_SQRT_5, FRAC_3_SQRT_20, -TWO_COS_72_DIV_SQRT_5], [TWO_SIN_36_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("12C5", [[0.5, -SIN_72_MINUS_SIN_36_DIV_SQRT_5, TWO_SIN_72_DIV_SQRT_5], [SIN_72_PLUS_SIN_36_DIV_SQRT_5, FRAC_3_SQRT_20, -TWO_COS_72_DIV_SQRT_5], [-TWO_SIN_36_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("12C5", [[COS_36, SIN_36_DIV_SQRT_5, TWO_SIN_36_DIV_SQRT_5], [-SIN_36, COS_36_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5], [0.0, -FRAC_2_SQRT_5, FRAC_1_SQRT_5]]),
    ("12C5", [[COS_36, -SIN_36_DIV_SQRT_5, -TWO_SIN_36_DIV_SQRT_5], [SIN_36, COS_36_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5], [0.0, -FRAC_2_SQRT_5, FRAC_1_SQRT_5]]),
    ("12C5", [[COS_36, SIN_36, 0.0], [-SIN_36_DIV_SQRT_5, COS_36_DIV_SQRT_5, -FRAC_2_SQRT_5], [-TWO_SIN_36_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("12C5", [[COS_36, -SIN_36, 0.0], [SIN_36_DIV_SQRT_5, COS_36_DIV_SQRT_5, -FRAC_2_SQRT_5], [TWO_SIN_36_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("12C5^2", [[-COS_36, SIN_36, 0.0], [-SIN_36, -COS_36, 0.0], [0.0, 0.0, 1.0]]),
    ("12C5^2", [[-COS_36, -SIN_36_DIV_SQRT_5, -TWO_SIN_36_DIV_SQRT_5], [SIN_36_DIV_SQRT_5, HALF_PLUS_COS_72_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5], [TWO_SIN_36_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("12C5^2", [[-COS_36, SIN_36_DIV_SQRT_5, TWO_SIN_36_DIV_SQRT_5], [-SIN_36_DIV_SQRT_5, HALF_PLUS_COS_72_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5], [-TWO_SIN_36_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("12C5^2", [[-COS_36, -SIN_36, 0.0], [SIN_36, -COS_36, 0.0], [0.0, 0.0, 1.0]]),
    ("12C5^2", [[-COS_72, SIN_72, 0.0], [SIN_72_DIV_SQRT_5, COS_72_DIV_SQRT_5, FRAC_2_SQRT_5], [TWO_SIN_72_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("12C5^2", [[-COS_72, -SIN_72, 0.0], [-SIN_72_DIV_SQRT_5, COS_72_DIV_SQRT_5, FRAC_2_SQRT_5], [-TWO_SIN_72_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("12C5^2", [[-COS_72, -SIN_72_DIV_SQRT_5, -TWO_SIN_72_DIV_SQRT_5], [-SIN_72, COS_72_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5], [0.0, FRAC_2_SQRT_5, -FRAC_1_SQRT_5]]),
    ("12C5^2", [[-COS_72, SIN_72_DIV_SQRT_5, TWO_SIN_72_DIV_SQRT_5], [SIN_72, COS_72_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5], [0.0, FRAC_2_SQRT_5, -FRAC_1_SQRT_5]]),
    ("12C5^2", [[0.5, -SIN_72_MINUS_SIN_36_DIV_SQRT_5, TWO_SIN_72_DIV_SQRT_5], [-SIN_72_PLUS_SIN_36_DIV_SQRT_5, -FRAC_3_SQRT_20, TWO_COS_72_DIV_SQRT_5], [TWO_SIN_36_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("12C5^2", [[0.5, SIN_72_MINUS_SIN_36_DIV_SQRT_5, -TWO_SIN_72_DIV_SQRT_5], [SIN_72_PLUS_SIN_36_DIV_SQRT_5, -FRAC_3_SQRT_20, TWO_COS_72_DIV_SQRT_5], [-TWO_SIN_36_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("12C5^2", [[0.5, SIN_72_PLUS_SIN_36_DIV_SQRT_5, -TWO_SIN_36_DIV_SQRT_5], [SIN_72_MINUS_SIN_36_DIV_SQRT_5, -FRAC_3_SQRT_20, -TWO_COS_36_DIV_SQRT_5], [-TWO_SIN_72_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("12C5^2", [[0.5, -SIN_72_PLUS_SIN_36_DIV_SQRT_5, TWO_SIN_36_DIV_SQRT_5], [-SIN_72_MINUS_SIN_36_DIV_SQRT_5, -FRAC_3_SQRT_20, -TWO_COS_36_DIV_SQRT_5], [TWO_SIN_72_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("20C3", [[0.0, TWO_SIN_72_DIV_SQRT_5, TWO_SIN_36_DIV_SQRT_5], [-TWO_SIN_36_DIV_SQRT_5, -FRAC_1_SQRT_5, TWO_COS_36_DIV_SQRT_5], [TWO_SIN_72_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("20C3", [[0.0, -TWO_SIN_36_DIV_SQRT_5, TWO_SIN_72_DIV_SQRT_5], [TWO_SIN_72_DIV_SQRT_5, -FRAC_1_SQRT_5, -TWO_COS_72_DIV_SQRT_5], [TWO_SIN_36_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("20C3", [[COS_36, SIN_36_DIV_SQRT_5, TWO_SIN_36_DIV_SQRT_5], [SIN_36, -COS_36_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5], [0.0, FRAC_2_SQRT_5, -FRAC_1_SQRT_5]]),
    ("20C3", [[COS_36, SIN_36, 0.0], [SIN_36_DIV_SQRT_5, -COS_36_DIV_SQRT_5, FRAC_2_SQRT_5], [TWO_SIN_36_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("20C3", [[0.0, TWO_SIN_72_DIV_SQRT_5, TWO_SIN_36_DIV_SQRT_5], [TWO_SIN_36_DIV_SQRT_5, FRAC_1_SQRT_5, -TWO_COS_36_DIV_SQRT_5], [-TWO_SIN_72_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("20C3", [[0.0, TWO_SIN_36_DIV_SQRT_5, -TWO_SIN_72_DIV_SQRT_5], [TWO_SIN_72_DIV_SQRT_5, FRAC_1_SQRT_5, TWO_COS_72_DIV_SQRT_5], [TWO_SIN_36_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("20C3", [[COS_36, -SIN_36, 0.0], [-SIN_36_DIV_SQRT_5, -COS_36_DIV_SQRT_5, FRAC_2_SQRT_5], [-TWO_SIN_36_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("20C3", [[COS_36, -SIN_36_DIV_SQRT_5, -TWO_SIN_36_DIV_SQRT_5], [-SIN_36, -COS_36_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5], [0.0, FRAC_2_SQRT_5, -FRAC_1_SQRT_5]]),
    ("20C3", [[-COS_72, SIN_72_DIV_SQRT_5, TWO_SIN_72_DIV_SQRT_5], [-SIN_72, -COS_72_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5], [0.0, -FRAC_2_SQRT_5, FRAC_1_SQRT_5]]),
    ("20C3", [[-COS_72, -SIN_72, 0.0], [SIN_72_DIV_SQRT_5, -COS_72_DIV_SQRT_5, -FRAC_2_SQRT_5], [TWO_SIN_72_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("20C3", [[-COS_72, SIN_72, 0.0], [-SIN_72_DIV_SQRT_5, -COS_72_DIV_SQRT_5, -FRAC_2_SQRT_5], [-TWO_SIN_72_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("20C3", [[-COS_72, -SIN_72_DIV_SQRT_5, -TWO_SIN_72_DIV_SQRT_5], [SIN_72, -COS_72_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5], [0.0, -FRAC_2_SQRT_5, FRAC_1_SQRT_5]]),
    ("20C3", [[0.0, TWO_SIN_36_DIV_SQRT_5, -TWO_SIN_72_DIV_SQRT_5], [-TWO_SIN_72_DIV_SQRT_5, -FRAC_1_SQRT_5, -TWO_COS_72_DIV_SQRT_5], [-TWO_SIN_36_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("20C3", [[0.0, -TWO_SIN_72_DIV_SQRT_5, -TWO_SIN_36_DIV_SQRT_5], [TWO_SIN_36_DIV_SQRT_5, -FRAC_1_SQRT_5, TWO_COS_36_DIV_SQRT_5], [-TWO_SIN_72_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("20C3", [[-0.5, -SIN_72_MINUS_SIN_36_DIV_SQRT_5, TWO_SIN_72_DIV_SQRT_5], [SIN_72_MINUS_SIN_36_DIV_SQRT_5, HALF_PLUS_FRAC_1_SQRT_5, TWO_COS_72_DIV_SQRT_5], [-TWO_SIN_72_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("20C3", [[-0.5, SIN_72_MINUS_SIN_36_DIV_SQRT_5, -TWO_SIN_72_DIV_SQRT_5], [-SIN_72_MINUS_SIN_36_DIV_SQRT_5, HALF_PLUS_FRAC_1_SQRT_5, TWO_COS_72_DIV_SQRT_5], [TWO_SIN_72_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("20C3", [[0.0, -TWO_SIN_36_DIV_SQRT_5, TWO_SIN_72_DIV_SQRT_5], [-TWO_SIN_72_DIV_SQRT_5, FRAC_1_SQRT_5, TWO_COS_72_DIV_SQRT_5], [-TWO_SIN_36_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("20C3", [[0.0, -TWO_SIN_72_DIV_SQRT_5, -TWO_SIN_36_DIV_SQRT_5], [-TWO_SIN_36_DIV_SQRT_5, FRAC_1_SQRT_5, -TWO_COS_36_DIV_SQRT_5], [TWO_SIN_72_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("20C3", [[-0.5, SIN_72_PLUS_SIN_36_DIV_SQRT_5, -TWO_SIN_36_DIV_SQRT_5], [-SIN_72_PLUS_SIN_36_DIV_SQRT_5, HALF_MINUS_FRAC_1_SQRT_5, TWO_COS_36_DIV_SQRT_5], [TWO_SIN_36_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("20C3", [[-0.5, -SIN_72_PLUS_SIN_36_DIV_SQRT_5, TWO_SIN_36_DIV_SQRT_5], [SIN_72_PLUS_SIN_36_DIV_SQRT_5, HALF_MINUS_FRAC_1_SQRT_5, TWO_COS_36_DIV_SQRT_5], [-TWO_SIN_36_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("15C2", [[COS_72, -SIN_72_DIV_SQRT_5, -TWO_SIN_72_DIV_SQRT_5], [-SIN_72_DIV_SQRT_5, -HALF_PLUS_COS_36_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5], [-TWO_SIN_72_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("15C2", [[-1.0, 0.0, 0.0], [0.0, FRAC_1_SQRT_5, FRAC_2_SQRT_5], [0.0, FRAC_2_SQRT_5, -FRAC_1_SQRT_5]]),
    ("15C2", [[COS_72, SIN_72_DIV_SQRT_5, TWO_SIN_72_DIV_SQRT_5], [SIN_72_DIV_SQRT_5, -HALF_PLUS_COS_36_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5], [TWO_SIN_72_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("15C2", [[COS_72, -SIN_72, 0.0], [-SIN_72, -COS_72, 0.0], [0.0, 0.0, -1.0]]),
    ("15C2", [[COS_72, SIN_72, 0.0], [SIN_72, -COS_72, 0.0], [0.0, 0.0, -1.0]]),
    ("15C2", [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("15C2", [[-COS_36, -SIN_36_DIV_SQRT_5, -TWO_SIN_36_DIV_SQRT_5], [-SIN_36_DIV_SQRT_5, -HALF_PLUS_COS_72_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5], [-TWO_SIN_36_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("15C2", [[-COS_36, SIN_36_DIV_SQRT_5, TWO_SIN_36_DIV_SQRT_5], [SIN_36_DIV_SQRT_5, -HALF_PLUS_COS_72_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5], [TWO_SIN_36_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("15C2", [[-0.5, -SIN_72_MINUS_SIN_36_DIV_SQRT_5, TWO_SIN_72_DIV_SQRT_5], [-SIN_72_MINUS_SIN_36_DIV_SQRT_5, -HALF_PLUS_FRAC_1_SQRT_5, -TWO_COS_72_DIV_SQRT_5], [TWO_SIN_72_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("15C2", [[-1.0, 0.0, 0.0], [0.0, -FRAC_1_SQRT_5, -FRAC_2_SQRT_5], [0.0, -FRAC_2_SQRT_5, FRAC_1_SQRT_5]]),
    ("15C2", [[-0.5, SIN_72_MINUS_SIN_36_DIV_SQRT_5, -TWO_SIN_72_DIV_SQRT_5], [SIN_72_MINUS_SIN_36_DIV_SQRT_5, -HALF_PLUS_FRAC_1_SQRT_5, -TWO_COS_72_DIV_SQRT_5], [-TWO_SIN_72_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("15C2", [[-0.5, SIN_72_PLUS_SIN_36_DIV_SQRT_5, -TWO_SIN_36_DIV_SQRT_5], [SIN_72_PLUS_SIN_36_DIV_SQRT_5, -HALF_MINUS_FRAC_1_SQRT_5, -TWO_COS_36_DIV_SQRT_5], [-TWO_SIN_36_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("15C2", [[-0.5, -SIN_72_PLUS_SIN_36_DIV_SQRT_5, TWO_SIN_36_DIV_SQRT_5], [-SIN_72_PLUS_SIN_36_DIV_SQRT_5, -HALF_MINUS_FRAC_1_SQRT_5, -TWO_COS_36_DIV_SQRT_5], [TWO_SIN_36_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("15C2", [[-COS_36, -SIN_36, 0.0], [-SIN_36, COS_36, 0.0], [0.0, 0.0, -1.0]]),
    ("15C2", [[-COS_36, SIN_36, 0.0], [SIN_36, COS_36, 0.0], [0.0, 0.0, -1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_IH: &[(&str, [[f64; 3]; 3])] = &[
    ("12C5", [[COS_72, SIN_72, 0.0], [-SIN_72, COS_72, 0.0], [0.0, 0.0, 1.0]]),
    ("12C5", [[COS_72, -SIN_72_DIV_SQRT_5, -TWO_SIN_72_DIV_SQRT_5], [SIN_72_DIV_SQRT_5, HALF_PLUS_COS_36_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5], [TWO_SIN_72_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("12C5", [[COS_72, SIN_72_DIV_SQRT_5, TWO_SIN_72_DIV_SQRT_5], [-SIN_72_DIV_SQRT_5, HALF_PLUS_COS_36_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5], [-TWO_SIN_72_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("12C5", [[COS_72, -SIN_72, 0.0], [SIN_72, COS_72, 0.0], [0.0, 0.0, 1.0]]),
    ("12C5", [[0.5, SIN_72_PLUS_SIN_36_DIV_SQRT_5, -TWO_SIN_36_DIV_SQRT_5], [-SIN_72_MINUS_SIN_36_DIV_SQRT_5, FRAC_3_SQRT_20, TWO_COS_36_DIV_SQRT_5], [TWO_SIN_72_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("12C5", [[0.5, -SIN_72_PLUS_SIN_36_DIV_SQRT_5, TWO_SIN_36_DIV_SQRT_5], [SIN_72_MINUS_SIN_36_DIV_SQRT_5, FRAC_3_SQRT_20, TWO_COS_36_DIV_SQRT_5], [-TWO_SIN_72_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("12C5", [[0.5, SIN_72_MINUS_SIN_36_DIV_SQRT_5, -TWO_SIN_72_DIV_SQRT_5], [-SIN_72_PLUS_SIN_36_DIV_SQRT_5, FRAC_3_SQRT_20, -TWO_COS_72_DIV_SQRT_5], [TWO_SIN_36_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("12C5", [[0.5, -SIN_72_MINUS_SIN_36_DIV_SQRT_5, TWO_SIN_72_DIV_SQRT_5], [SIN_72_PLUS_SIN_36_DIV_SQRT_5, FRAC_3_SQRT_20, -TWO_COS_72_DIV_SQRT_5], [-TWO_SIN_36_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("12C5", [[COS_36, SIN_36_DIV_SQRT_5, TWO_SIN_36_DIV_SQRT_5], [-SIN_36, COS_36_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5], [0.0, -FRAC_2_SQRT_5, FRAC_1_SQRT_5]]),
    ("12C5", [[COS_36, -SIN_36_DIV_SQRT_5, -TWO_SIN_36_DIV_SQRT_5], [SIN_36, COS_36_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5], [0.0, -FRAC_2_SQRT_5, FRAC_1_SQRT_5]]),
    ("12C5", [[COS_36, SIN_36, 0.0], [-SIN_36_DIV_SQRT_5, COS_36_DIV_SQRT_5, -FRAC_2_SQRT_5], [-TWO_SIN_36_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("12C5", [[COS_36, -SIN_36, 0.0], [SIN_36_DIV_SQRT_5, COS_36_DIV_SQRT_5, -FRAC_2_SQRT_5], [TWO_SIN_36_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("12C5^2", [[-COS_36, SIN_36, 0.0], [-SIN_36, -COS_36, 0.0], [0.0, 0.0, 1.0]]),
    ("12C5^2", [[-COS_36, -SIN_36_DIV_SQRT_5, -TWO_SIN_36_DIV_SQRT_5], [SIN_36_DIV_SQRT_5, HALF_PLUS_COS_72_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5], [TWO_SIN_36_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("12C5^2", [[-COS_36, SIN_36_DIV_SQRT_5, TWO_SIN_36_DIV_SQRT_5], [-SIN_36_DIV_SQRT_5, HALF_PLUS_COS_72_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5], [-TWO_SIN_36_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("12C5^2", [[-COS_36, -SIN_36, 0.0], [SIN_36, -COS_36, 0.0], [0.0, 0.0, 1.0]]),
    ("12C5^2", [[-COS_72, SIN_72, 0.0], [SIN_72_DIV_SQRT_5, COS_72_DIV_SQRT_5, FRAC_2_SQRT_5], [TWO_SIN_72_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("12C5^2", [[-COS_72, -SIN_72, 0.0], [-SIN_72_DIV_SQRT_5, COS_72_DIV_SQRT_5, FRAC_2_SQRT_5], [-TWO_SIN_72_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("12C5^2", [[-COS_72, -SIN_72_DIV_SQRT_5, -TWO_SIN_72_DIV_SQRT_5], [-SIN_72, COS_72_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5], [0.0, FRAC_2_SQRT_5, -FRAC_1_SQRT_5]]),
    ("12C5^2", [[-COS_72, SIN_72_DIV_SQRT_5, TWO_SIN_72_DIV_SQRT_5], [SIN_72, COS_72_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5], [0.0, FRAC_2_SQRT_5, -FRAC_1_SQRT_5]]),
    ("12C5^2", [[0.5, -SIN_72_MINUS_SIN_36_DIV_SQRT_5, TWO_SIN_72_DIV_SQRT_5], [-SIN_72_PLUS_SIN_36_DIV_SQRT_5, -FRAC_3_SQRT_20, TWO_COS_72_DIV_SQRT_5], [TWO_SIN_36_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("12C5^2", [[0.5, SIN_72_MINUS_SIN_36_DIV_SQRT_5, -TWO_SIN_72_DIV_SQRT_5], [SIN_72_PLUS_SIN_36_DIV_SQRT_5, -FRAC_3_SQRT_20, TWO_COS_72_DIV_SQRT_5], [-TWO_SIN_36_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("12C5^2", [[0.5, SIN_72_PLUS_SIN_36_DIV_SQRT_5, -TWO_SIN_36_DIV_SQRT_5], [SIN_72_MINUS_SIN_36_DIV_SQRT_5, -FRAC_3_SQRT_20, -TWO_COS_36_DIV_SQRT_5], [-TWO_SIN_72_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("12C5^2", [[0.5, -SIN_72_PLUS_SIN_36_DIV_SQRT_5, TWO_SIN_36_DIV_SQRT_5], [-SIN_72_MINUS_SIN_36_DIV_SQRT_5, -FRAC_3_SQRT_20, -TWO_COS_36_DIV_SQRT_5], [TWO_SIN_72_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("20C3", [[0.0, TWO_SIN_72_DIV_SQRT_5, TWO_SIN_36_DIV_SQRT_5], [-TWO_SIN_36_DIV_SQRT_5, -FRAC_1_SQRT_5, TWO_COS_36_DIV_SQRT_5], [TWO_SIN_72_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("20C3", [[0.0, -TWO_SIN_36_DIV_SQRT_5, TWO_SIN_72_DIV_SQRT_5], [TWO_SIN_72_DIV_SQRT_5, -FRAC_1_SQRT_5, -TWO_COS_72_DIV_SQRT_5], [TWO_SIN_36_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("20C3", [[COS_36, SIN_36_DIV_SQRT_5, TWO_SIN_36_DIV_SQRT_5], [SIN_36, -COS_36_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5], [0.0, FRAC_2_SQRT_5, -FRAC_1_SQRT_5]]),
    ("20C3", [[COS_36, SIN_36, 0.0], [SIN_36_DIV_SQRT_5, -COS_36_DIV_SQRT_5, FRAC_2_SQRT_5], [TWO_SIN_36_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("20C3", [[0.0, TWO_SIN_72_DIV_SQRT_5, TWO_SIN_36_DIV_SQRT_5], [TWO_SIN_36_DIV_SQRT_5, FRAC_1_SQRT_5, -TWO_COS_36_DIV_SQRT_5], [-TWO_SIN_72_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("20C3", [[0.0, TWO_SIN_36_DIV_SQRT_5, -TWO_SIN_72_DIV_SQRT_5], [TWO_SIN_72_DIV_SQRT_5, FRAC_1_SQRT_5, TWO_COS_72_DIV_SQRT_5], [TWO_SIN_36_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("20C3", [[COS_36, -SIN_36, 0.0], [-SIN_36_DIV_SQRT_5, -COS_36_DIV_SQRT_5, FRAC_2_SQRT_5], [-TWO_SIN_36_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("20C3", [[COS_36, -SIN_36_DIV_SQRT_5, -TWO_SIN_36_DIV_SQRT_5], [-SIN_36, -COS_36_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5], [0.0, FRAC_2_SQRT_5, -FRAC_1_SQRT_5]]),
    ("20C3", [[-COS_72, SIN_72_DIV_SQRT_5, TWO_SIN_72_DIV_SQRT_5], [-SIN_72, -COS_72_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5], [0.0, -FRAC_2_SQRT_5, FRAC_1_SQRT_5]]),
    ("20C3", [[-COS_72, -SIN_72, 0.0], [SIN_72_DIV_SQRT_5, -COS_72_DIV_SQRT_5, -FRAC_2_SQRT_5], [TWO_SIN_72_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("20C3", [[-COS_72, SIN_72, 0.0], [-SIN_72_DIV_SQRT_5, -COS_72_DIV_SQRT_5, -FRAC_2_SQRT_5], [-TWO_SIN_72_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("20C3", [[-COS_72, -SIN_72_DIV_SQRT_5, -TWO_SIN_72_DIV_SQRT_5], [SIN_72, -COS_72_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5], [0.0, -FRAC_2_SQRT_5, FRAC_1_SQRT_5]]),
    ("20C3", [[0.0, TWO_SIN_36_DIV_SQRT_5, -TWO_SIN_72_DIV_SQRT_5], [-TWO_SIN_72_DIV_SQRT_5, -FRAC_1_SQRT_5, -TWO_COS_72_DIV_SQRT_5], [-TWO_SIN_36_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("20C3", [[0.0, -TWO_SIN_72_DIV_SQRT_5, -TWO_SIN_36_DIV_SQRT_5], [TWO_SIN_36_DIV_SQRT_5, -FRAC_1_SQRT_5, TWO_COS_36_DIV_SQRT_5], [-TWO_SIN_72_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("20C3", [[-0.5, -SIN_72_MINUS_SIN_36_DIV_SQRT_5, TWO_SIN_72_DIV_SQRT_5], [SIN_72_MINUS_SIN_36_DIV_SQRT_5, HALF_PLUS_FRAC_1_SQRT_5, TWO_COS_72_DIV_SQRT_5], [-TWO_SIN_72_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("20C3", [[-0.5, SIN_72_MINUS_SIN_36_DIV_SQRT_5, -TWO_SIN_72_DIV_SQRT_5], [-SIN_72_MINUS_SIN_36_DIV_SQRT_5, HALF_PLUS_FRAC_1_SQRT_5, TWO_COS_72_DIV_SQRT_5], [TWO_SIN_72_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("20C3", [[0.0, -TWO_SIN_36_DIV_SQRT_5, TWO_SIN_72_DIV_SQRT_5], [-TWO_SIN_72_DIV_SQRT_5, FRAC_1_SQRT_5, TWO_COS_72_DIV_SQRT_5], [-TWO_SIN_36_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("20C3", [[0.0, -TWO_SIN_72_DIV_SQRT_5, -TWO_SIN_36_DIV_SQRT_5], [-TWO_SIN_36_DIV_SQRT_5, FRAC_1_SQRT_5, -TWO_COS_36_DIV_SQRT_5], [TWO_SIN_72_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("20C3", [[-0.5, SIN_72_PLUS_SIN_36_DIV_SQRT_5, -TWO_SIN_36_DIV_SQRT_5], [-SIN_72_PLUS_SIN_36_DIV_SQRT_5, HALF_MINUS_FRAC_1_SQRT_5, TWO_COS_36_DIV_SQRT_5], [TWO_SIN_36_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("20C3", [[-0.5, -SIN_72_PLUS_SIN_36_DIV_SQRT_5, TWO_SIN_36_DIV_SQRT_5], [SIN_72_PLUS_SIN_36_DIV_SQRT_5, HALF_MINUS_FRAC_1_SQRT_5, TWO_COS_36_DIV_SQRT_5], [-TWO_SIN_36_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("15C2", [[COS_72, -SIN_72_DIV_SQRT_5, -TWO_SIN_72_DIV_SQRT_5], [-SIN_72_DIV_SQRT_5, -HALF_PLUS_COS_36_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5], [-TWO_SIN_72_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("15C2", [[-1.0, 0.0, 0.0], [0.0, FRAC_1_SQRT_5, FRAC_2_SQRT_5], [0.0, FRAC_2_SQRT_5, -FRAC_1_SQRT_5]]),
    ("15C2", [[COS_72, SIN_72_DIV_SQRT_5, TWO_SIN_72_DIV_SQRT_5], [SIN_72_DIV_SQRT_5, -HALF_PLUS_COS_36_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5], [TWO_SIN_72_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("15C2", [[COS_72, -SIN_72, 0.0], [-SIN_72, -COS_72, 0.0], [0.0, 0.0, -1.0]]),
    ("15C2", [[COS_72, SIN_72, 0.0], [SIN_72, -COS_72, 0.0], [0.0, 0.0, -1.0]]),
    ("15C2", [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("15C2", [[-COS_36, -SIN_36_DIV_SQRT_5, -TWO_SIN_36_DIV_SQRT_5], [-SIN_36_DIV_SQRT_5, -HALF_PLUS_COS_72_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5], [-TWO_SIN_36_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("15C2", [[-COS_36, SIN_36_DIV_SQRT_5, TWO_SIN_36_DIV_SQRT_5], [SIN_36_DIV_SQRT_5, -HALF_PLUS_COS_72_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5], [TWO_SIN_36_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("15C2", [[-0.5, -SIN_72_MINUS_SIN_36_DIV_SQRT_5, TWO_SIN_72_DIV_SQRT_5], [-SIN_72_MINUS_SIN_36_DIV_SQRT_5, -HALF_PLUS_FRAC_1_SQRT_5, -TWO_COS_72_DIV_SQRT_5], [TWO_SIN_72_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("15C2", [[-1.0, 0.0, 0.0], [0.0, -FRAC_1_SQRT_5, -FRAC_2_SQRT_5], [0.0, -FRAC_2_SQRT_5, FRAC_1_SQRT_5]]),
    ("15C2", [[-0.5, SIN_72_MINUS_SIN_36_DIV_SQRT_5, -TWO_SIN_72_DIV_SQRT_5], [SIN_72_MINUS_SIN_36_DIV_SQRT_5, -HALF_PLUS_FRAC_1_SQRT_5, -TWO_COS_72_DIV_SQRT_5], [-TWO_SIN_72_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("15C2", [[-0.5, SIN_72_PLUS_SIN_36_DIV_SQRT_5, -TWO_SIN_36_DIV_SQRT_5], [SIN_72_PLUS_SIN_36_DIV_SQRT_5, -HALF_MINUS_FRAC_1_SQRT_5, -TWO_COS_36_DIV_SQRT_5], [-TWO_SIN_36_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("15C2", [[-0.5, -SIN_72_PLUS_SIN_36_DIV_SQRT_5, TWO_SIN_36_DIV_SQRT_5], [-SIN_72_PLUS_SIN_36_DIV_SQRT_5, -HALF_MINUS_FRAC_1_SQRT_5, -TWO_COS_36_DIV_SQRT_5], [TWO_SIN_36_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("15C2", [[-COS_36, -SIN_36, 0.0], [-SIN_36, COS_36, 0.0], [0.0, 0.0, -1.0]]),
    ("15C2", [[-COS_36, SIN_36, 0.0], [SIN_36, COS_36, 0.0], [0.0, 0.0, -1.0]]),
    ("i", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("12S10", [[COS_36, SIN_36, 0.0], [-SIN_36, COS_36, 0.0], [0.0, 0.0, -1.0]]),
    ("12S10", [[COS_36, -SIN_36_DIV_SQRT_5, -TWO_SIN_36_DIV_SQRT_5], [SIN_36_DIV_SQRT_5, -HALF_PLUS_COS_72_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5], [TWO_SIN_36_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("12S10", [[COS_36, SIN_36_DIV_SQRT_5, TWO_SIN_36_DIV_SQRT_5], [-SIN_36_DIV_SQRT_5, -HALF_PLUS_COS_72_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5], [-TWO_SIN_36_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("12S10", [[COS_36, -SIN_36, 0.0], [SIN_36, COS_36, 0.0], [0.0, 0.0, -1.0]]),
    ("12S10", [[COS_72, -SIN_72_DIV_SQRT_5, -TWO_SIN_72_DIV_SQRT_5], [-SIN_72, -COS_72_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5], [0.0, -FRAC_2_SQRT_5, FRAC_1_SQRT_5]]),
    ("12S10", [[COS_72, SIN_72_DIV_SQRT_5, TWO_SIN_72_DIV_SQRT_5], [SIN_72, -COS_72_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5], [0.0, -FRAC_2_SQRT_5, FRAC_1_SQRT_5]]),
    ("12S10", [[COS_72, SIN_72, 0.0], [SIN_72_DIV_SQRT_5, -COS_72_DIV_SQRT_5, -FRAC_2_SQRT_5], [TWO_SIN_72_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("12S10", [[COS_72, -SIN_72, 0.0], [-SIN_72_DIV_SQRT_5, -COS_72_DIV_SQRT_5, -FRAC_2_SQRT_5], [-TWO_SIN_72_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("12S10", [[-0.5, SIN_72_PLUS_SIN_36_DIV_SQRT_5, -TWO_SIN_36_DIV_SQRT_5], [SIN_72_MINUS_SIN_36_DIV_SQRT_5, FRAC_3_SQRT_20, TWO_COS_36_DIV_SQRT_5], [-TWO_SIN_72_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("12S10", [[-0.5, -SIN_72_PLUS_SIN_36_DIV_SQRT_5, TWO_SIN_36_DIV_SQRT_5], [-SIN_72_MINUS_SIN_36_DIV_SQRT_5, FRAC_3_SQRT_20, TWO_COS_36_DIV_SQRT_5], [TWO_SIN_72_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("12S10", [[-0.5, -SIN_72_MINUS_SIN_36_DIV_SQRT_5, TWO_SIN_72_DIV_SQRT_5], [-SIN_72_PLUS_SIN_36_DIV_SQRT_5, FRAC_3_SQRT_20, -TWO_COS_72_DIV_SQRT_5], [TWO_SIN_36_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("12S10", [[-0.5, SIN_72_MINUS_SIN_36_DIV_SQRT_5, -TWO_SIN_72_DIV_SQRT_5], [SIN_72_PLUS_SIN_36_DIV_SQRT_5, FRAC_3_SQRT_20, -TWO_COS_72_DIV_SQRT_5], [-TWO_SIN_36_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("12S10^3", [[-COS_72, -SIN_72, 0.0], [SIN_72, -COS_72, 0.0], [0.0, 0.0, -1.0]]),
    ("12S10^3", [[-COS_72, SIN_72_DIV_SQRT_5, TWO_SIN_72_DIV_SQRT_5], [-SIN_72_DIV_SQRT_5, -HALF_PLUS_COS_36_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5], [-TWO_SIN_72_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("12S10^3", [[-COS_72, -SIN_72_DIV_SQRT_5, -TWO_SIN_72_DIV_SQRT_5], [SIN_72_DIV_SQRT_5, -HALF_PLUS_COS_36_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5], [TWO_SIN_72_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("12S10^3", [[-COS_72, SIN_72, 0.0], [-SIN_72, -COS_72, 0.0], [0.0, 0.0, -1.0]]),
    ("12S10^3", [[-0.5, -SIN_72_PLUS_SIN_36_DIV_SQRT_5, TWO_SIN_36_DIV_SQRT_5], [SIN_72_MINUS_SIN_36_DIV_SQRT_5, -FRAC_3_SQRT_20, -TWO_COS_36_DIV_SQRT_5], [-TWO_SIN_72_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("12S10^3", [[-0.5, SIN_72_PLUS_SIN_36_DIV_SQRT_5, -TWO_SIN_36_DIV_SQRT_5], [-SIN_72_MINUS_SIN_36_DIV_SQRT_5, -FRAC_3_SQRT_20, -TWO_COS_36_DIV_SQRT_5], [TWO_SIN_72_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("12S10^3", [[-0.5, -SIN_72_MINUS_SIN_36_DIV_SQRT_5, TWO_SIN_72_DIV_SQRT_5], [SIN_72_PLUS_SIN_36_DIV_SQRT_5, -FRAC_3_SQRT_20, TWO_COS_72_DIV_SQRT_5], [-TWO_SIN_36_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("12S10^3", [[-0.5, SIN_72_MINUS_SIN_36_DIV_SQRT_5, -TWO_SIN_72_DIV_SQRT_5], [-SIN_72_PLUS_SIN_36_DIV_SQRT_5, -FRAC_3_SQRT_20, TWO_COS_72_DIV_SQRT_5], [TWO_SIN_36_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("12S10^3", [[-COS_36, -SIN_36_DIV_SQRT_5, -TWO_SIN_36_DIV_SQRT_5], [SIN_36, -COS_36_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5], [0.0, FRAC_2_SQRT_5, -FRAC_1_SQRT_5]]),
    ("12S10^3", [[-COS_36, SIN_36_DIV_SQRT_5, TWO_SIN_36_DIV_SQRT_5], [-SIN_36, -COS_36_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5], [0.0, FRAC_2_SQRT_5, -FRAC_1_SQRT_5]]),
    ("12S10^3", [[-COS_36, -SIN_36, 0.0], [SIN_36_DIV_SQRT_5, -COS_36_DIV_SQRT_5, FRAC_2_SQRT_5], [TWO_SIN_36_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("12S10^3", [[-COS_36, SIN_36, 0.0], [-SIN_36_DIV_SQRT_5, -COS_36_DIV_SQRT_5, FRAC_2_SQRT_5], [-TWO_SIN_36_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("20S6", [[0.0, TWO_SIN_36_DIV_SQRT_5, -TWO_SIN_72_DIV_SQRT_5], [-TWO_SIN_72_DIV_SQRT_5, FRAC_1_SQRT_5, TWO_COS_72_DIV_SQRT_5], [-TWO_SIN_36_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("20S6", [[0.0, -TWO_SIN_72_DIV_SQRT_5, -TWO_SIN_36_DIV_SQRT_5], [TWO_SIN_36_DIV_SQRT_5, FRAC_1_SQRT_5, -TWO_COS_36_DIV_SQRT_5], [-TWO_SIN_72_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("20S6", [[-COS_36, -SIN_36, 0.0], [-SIN_36_DIV_SQRT_5, COS_36_DIV_SQRT_5, -FRAC_2_SQRT_5], [-TWO_SIN_36_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("20S6", [[-COS_36, -SIN_36_DIV_SQRT_5, -TWO_SIN_36_DIV_SQRT_5], [-SIN_36, COS_36_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5], [0.0, -FRAC_2_SQRT_5, FRAC_1_SQRT_5]]),
    ("20S6", [[0.0, -TWO_SIN_36_DIV_SQRT_5, TWO_SIN_72_DIV_SQRT_5], [-TWO_SIN_72_DIV_SQRT_5, -FRAC_1_SQRT_5, -TWO_COS_72_DIV_SQRT_5], [-TWO_SIN_36_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("20S6", [[0.0, -TWO_SIN_72_DIV_SQRT_5, -TWO_SIN_36_DIV_SQRT_5], [-TWO_SIN_36_DIV_SQRT_5, -FRAC_1_SQRT_5, TWO_COS_36_DIV_SQRT_5], [TWO_SIN_72_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("20S6", [[-COS_36, SIN_36_DIV_SQRT_5, TWO_SIN_36_DIV_SQRT_5], [SIN_36, COS_36_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5], [0.0, -FRAC_2_SQRT_5, FRAC_1_SQRT_5]]),
    ("20S6", [[-COS_36, SIN_36, 0.0], [SIN_36_DIV_SQRT_5, COS_36_DIV_SQRT_5, -FRAC_2_SQRT_5], [TWO_SIN_36_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("20S6", [[COS_72, SIN_72, 0.0], [-SIN_72_DIV_SQRT_5, COS_72_DIV_SQRT_5, FRAC_2_SQRT_5], [-TWO_SIN_72_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("20S6", [[COS_72, -SIN_72_DIV_SQRT_5, -TWO_SIN_72_DIV_SQRT_5], [SIN_72, COS_72_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5], [0.0, FRAC_2_SQRT_5, -FRAC_1_SQRT_5]]),
    ("20S6", [[COS_72, SIN_72_DIV_SQRT_5, TWO_SIN_72_DIV_SQRT_5], [-SIN_72, COS_72_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5], [0.0, FRAC_2_SQRT_5, -FRAC_1_SQRT_5]]),
    ("20S6", [[COS_72, -SIN_72, 0.0], [SIN_72_DIV_SQRT_5, COS_72_DIV_SQRT_5, FRAC_2_SQRT_5], [TWO_SIN_72_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("20S6", [[0.0, TWO_SIN_72_DIV_SQRT_5, TWO_SIN_36_DIV_SQRT_5], [-TWO_SIN_36_DIV_SQRT_5, FRAC_1_SQRT_5, -TWO_COS_36_DIV_SQRT_5], [TWO_SIN_72_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("20S6", [[0.0, -TWO_SIN_36_DIV_SQRT_5, TWO_SIN_72_DIV_SQRT_5], [TWO_SIN_72_DIV_SQRT_5, FRAC_1_SQRT_5, TWO_COS_72_DIV_SQRT_5], [TWO_SIN_36_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("20S6", [[0.5, -SIN_72_MINUS_SIN_36_DIV_SQRT_5, TWO_SIN_72_DIV_SQRT_5], [SIN_72_MINUS_SIN_36_DIV_SQRT_5, -HALF_PLUS_FRAC_1_SQRT_5, -TWO_COS_72_DIV_SQRT_5], [-TWO_SIN_72_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("20S6", [[0.5, SIN_72_MINUS_SIN_36_DIV_SQRT_5, -TWO_SIN_72_DIV_SQRT_5], [-SIN_72_MINUS_SIN_36_DIV_SQRT_5, -HALF_PLUS_FRAC_1_SQRT_5, -TWO_COS_72_DIV_SQRT_5], [TWO_SIN_72_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("20S6", [[0.0, TWO_SIN_72_DIV_SQRT_5, TWO_SIN_36_DIV_SQRT_5], [TWO_SIN_36_DIV_SQRT_5, -FRAC_1_SQRT_5, TWO_COS_36_DIV_SQRT_5], [-TWO_SIN_72_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("20S6", [[0.0, TWO_SIN_36_DIV_SQRT_5, -TWO_SIN_72_DIV_SQRT_5], [TWO_SIN_72_DIV_SQRT_5, -FRAC_1_SQRT_5, -TWO_COS_72_DIV_SQRT_5], [TWO_SIN_36_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("20S6", [[0.5, SIN_72_PLUS_SIN_36_DIV_SQRT_5, -TWO_SIN_36_DIV_SQRT_5], [-SIN_72_PLUS_SIN_36_DIV_SQRT_5, -HALF_MINUS_FRAC_1_SQRT_5, -TWO_COS_36_DIV_SQRT_5], [TWO_SIN_36_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("20S6", [[0.5, -SIN_72_PLUS_SIN_36_DIV_SQRT_5, TWO_SIN_36_DIV_SQRT_5], [SIN_72_PLUS_SIN_36_DIV_SQRT_5, -HALF_MINUS_FRAC_1_SQRT_5, -TWO_COS_36_DIV_SQRT_5], [-TWO_SIN_36_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("15sigma", [[-COS_72, SIN_72_DIV_SQRT_5, TWO_SIN_72_DIV_SQRT_5], [SIN_72_DIV_SQRT_5, HALF_PLUS_COS_36_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5], [TWO_SIN_72_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("15sigma", [[1.0, 0.0, 0.0], [0.0, -FRAC_1_SQRT_5, -FRAC_2_SQRT_5], [0.0, -FRAC_2_SQRT_5, FRAC_1_SQRT_5]]),
    ("15sigma", [[-COS_72, -SIN_72_DIV_SQRT_5, -TWO_SIN_72_DIV_SQRT_5], [-SIN_72_DIV_SQRT_5, HALF_PLUS_COS_36_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5], [-TWO_SIN_72_DIV_SQRT_5, -TWO_COS_72_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("15sigma", [[-COS_72, SIN_72, 0.0], [SIN_72, COS_72, 0.0], [0.0, 0.0, 1.0]]),
    ("15sigma", [[-COS_72, -SIN_72, 0.0], [-SIN_72, COS_72, 0.0], [0.0, 0.0, 1.0]]),
    ("15sigma", [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("15sigma", [[COS_36, SIN_36_DIV_SQRT_5, TWO_SIN_36_DIV_SQRT_5], [SIN_36_DIV_SQRT_5, HALF_PLUS_COS_72_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5], [TWO_SIN_36_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("15sigma", [[COS_36, -SIN_36_DIV_SQRT_5, -TWO_SIN_36_DIV_SQRT_5], [-SIN_36_DIV_SQRT_5, HALF_PLUS_COS_72_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5], [-TWO_SIN_36_DIV_SQRT_5, -TWO_COS_36_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("15sigma", [[0.5, SIN_72_MINUS_SIN_36_DIV_SQRT_5, -TWO_SIN_72_DIV_SQRT_5], [SIN_72_MINUS_SIN_36_DIV_SQRT_5, HALF_PLUS_FRAC_1_SQRT_5, TWO_COS_72_DIV_SQRT_5], [-TWO_SIN_72_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("15sigma", [[1.0, 0.0, 0.0], [0.0, FRAC_1_SQRT_5, FRAC_2_SQRT_5], [0.0, FRAC_2_SQRT_5, -FRAC_1_SQRT_5]]),
    ("15sigma", [[0.5, -SIN_72_MINUS_SIN_36_DIV_SQRT_5, TWO_SIN_72_DIV_SQRT_5], [-SIN_72_MINUS_SIN_36_DIV_SQRT_5, HALF_PLUS_FRAC_1_SQRT_5, TWO_COS_72_DIV_SQRT_5], [TWO_SIN_72_DIV_SQRT_5, TWO_COS_72_DIV_SQRT_5, -FRAC_1_SQRT_5]]),
    ("15sigma", [[0.5, -SIN_72_PLUS_SIN_36_DIV_SQRT_5, TWO_SIN_36_DIV_SQRT_5], [-SIN_72_PLUS_SIN_36_DIV_SQRT_5, HALF_MINUS_FRAC_1_SQRT_5, TWO_COS_36_DIV_SQRT_5], [TWO_SIN_36_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("15sigma", [[0.5, SIN_72_PLUS_SIN_36_DIV_SQRT_5, -TWO_SIN_36_DIV_SQRT_5], [SIN_72_PLUS_SIN_36_DIV_SQRT_5, HALF_MINUS_FRAC_1_SQRT_5, TWO_COS_36_DIV_SQRT_5], [-TWO_SIN_36_DIV_SQRT_5, TWO_COS_36_DIV_SQRT_5, FRAC_1_SQRT_5]]),
    ("15sigma", [[COS_36, SIN_36, 0.0], [SIN_36, -COS_36, 0.0], [0.0, 0.0, 1.0]]),
    ("15sigma", [[COS_36, -SIN_36, 0.0], [-SIN_36, -COS_36, 0.0], [0.0, 0.0, 1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_O: &[(&str, [[f64; 3]; 3])] = &[
    ("8C3", [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]]),
    ("8C3", [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
    ("8C3", [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]]),
    ("8C3", [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]]),
    ("8C3", [[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]]),
    ("8C3", [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]]),
    ("8C3", [[0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]]),
    ("8C3", [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
    ("6C2", [[0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]]),
    ("6C2", [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [-1.0, 0.0, 0.0]]),
    ("6C2", [[-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]]),
    ("6C2", [[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]]),
    ("6C2", [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]),
    ("6C2", [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]),
    ("6C4", [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]]),
    ("6C4", [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]]),
    ("6C4", [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]]),
    ("6C4", [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]]),
    ("6C4", [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
    ("6C4", [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
    ("3C2", [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("3C2", [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("3C2", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_OH: &[(&str, [[f64; 3]; 3])] = &[
    ("8C3", [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]]),
    ("8C3", [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
    ("8C3", [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]]),
    ("8C3", [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]]),
    ("8C3", [[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]]),
    ("8C3", [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]]),
    ("8C3", [[0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]]),
    ("8C3", [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
    ("6C2", [[0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]]),
    ("6C2", [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [-1.0, 0.0, 0.0]]),
    ("6C2", [[-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]]),
    ("6C2", [[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]]),
    ("6C2", [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]),
    ("6C2", [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]),
    ("6C4", [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]]),
    ("6C4", [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]]),
    ("6C4", [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]]),
    ("6C4", [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]]),
    ("6C4", [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
    ("6C4", [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
    ("3C2", [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("3C2", [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("3C2", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("i", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("6S4", [[-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]]),
    ("6S4", [[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]]),
    ("6S4", [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]]),
    ("6S4", [[0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [-1.0, 0.0, 0.0]]),
    ("6S4", [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]),
    ("6S4", [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]),
    ("8S6", [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]]),
    ("8S6", [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]]),
    ("8S6", [[0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]]),
    ("8S6", [[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
    ("8S6", [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]]),
    ("8S6", [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
    ("8S6", [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]]),
    ("8S6", [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]]),
    ("3sigma_h", [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("3sigma_h", [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("3sigma_h", [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("6sigma_d", [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
    ("6sigma_d", [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]]),
    ("6sigma_d", [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]]),
    ("6sigma_d", [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
    ("6sigma_d", [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]]),
    ("6sigma_d", [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_S10: &[(&str, [[f64; 3]; 3])] = &[
    ("C5^4", [[COS_72, SIN_72, 0.0], [-SIN_72, COS_72, 0.0], [0.0, 0.0, 1.0]]),
    ("C5", [[COS_72, -SIN_72, 0.0], [SIN_72, COS_72, 0.0], [0.0, 0.0, 1.0]]),
    ("C5^3", [[-COS_36, SIN_36, 0.0], [-SIN_36, -COS_36, 0.0], [0.0, 0.0, 1.0]]),
    ("C5^2", [[-COS_36, -SIN_36, 0.0], [SIN_36, -COS_36, 0.0], [0.0, 0.0, 1.0]]),
    ("i", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("S10^7", [[-COS_72, SIN_72, 0.0], [-SIN_72, -COS_72, 0.0], [0.0, 0.0, -1.0]]),
    ("S10^3", [[-COS_72, -SIN_72, 0.0], [SIN_72, -COS_72, 0.0], [0.0, 0.0, -1.0]]),
    ("S10", [[COS_36, -SIN_36, 0.0], [SIN_36, COS_36, 0.0], [0.0, 0.0, -1.0]]),
    ("S10^9", [[COS_36, SIN_36, 0.0], [-SIN_36, COS_36, 0.0], [0.0, 0.0, -1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_S4: &[(&str, [[f64; 3]; 3])] = &[
    ("C2", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("S4^3", [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]),
    ("S4", [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_S6: &[(&str, [[f64; 3]; 3])] = &[
    ("C3^2", [[-0.5, SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("C3", [[-0.5, -SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, -0.5, 0.0], [0.0, 0.0, 1.0]]),
    ("i", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("S6^5", [[0.5, SQRT_3_DIV_2, 0.0], [-SQRT_3_DIV_2, 0.5, 0.0], [0.0, 0.0, -1.0]]),
    ("S6", [[0.5, -SQRT_3_DIV_2, 0.0], [SQRT_3_DIV_2, 0.5, 0.0], [0.0, 0.0, -1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_S8: &[(&str, [[f64; 3]; 3])] = &[
    ("S8^7", [[FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0], [-FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0], [0.0, 0.0, -1.0]]),
    ("S8^5", [[-FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0], [-FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0], [0.0, 0.0, -1.0]]),
    ("S8^3", [[-FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0], [FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0], [0.0, 0.0, -1.0]]),
    ("S8", [[FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0], [FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0], [0.0, 0.0, -1.0]]),
    ("C4^3", [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
    ("C2", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("C4", [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_T: &[(&str, [[f64; 3]; 3])] = &[
    ("4C3^2", [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]]),
    ("4C3", [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
    ("4C3", [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]]),
    ("4C3^2", [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]]),
    ("4C3", [[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]]),
    ("4C3^2", [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]]),
    ("4C3^2", [[0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]]),
    ("4C3", [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
    ("3C2", [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("3C2", [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("3C2", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_TD: &[(&str, [[f64; 3]; 3])] = &[
    ("8C3", [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]]),
    ("8C3", [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
    ("8C3", [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]]),
    ("8C3", [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]]),
    ("8C3", [[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]]),
    ("8C3", [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]]),
    ("8C3", [[0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]]),
    ("8C3", [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
    ("3C2", [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("3C2", [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("3C2", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("6S4", [[-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]]),
    ("6S4", [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]]),
    ("6S4", [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]),
    ("6S4", [[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]]),
    ("6S4", [[0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [-1.0, 0.0, 0.0]]),
    ("6S4", [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]),
    ("6sigma_d", [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
    ("6sigma_d", [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]]),
    ("6sigma_d", [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]]),
    ("6sigma_d", [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
    ("6sigma_d", [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]]),
    ("6sigma_d", [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]]),
];

#[rustfmt::skip]
pub static POINTGROUP_TH: &[(&str, [[f64; 3]; 3])] = &[
    ("4C3^2", [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]]),
    ("4C3", [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
    ("4C3", [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]]),
    ("4C3^2", [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]]),
    ("4C3", [[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]]),
    ("4C3^2", [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]]),
    ("4C3^2", [[0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]]),
    ("4C3", [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
    ("3C2", [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("3C2", [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("3C2", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("i", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]),
    ("4S6^5", [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]]),
    ("4S6", [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]]),
    ("4S6", [[0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]]),
    ("4S6^5", [[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
    ("4S6", [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]]),
    ("4S6^5", [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
    ("4S6^5", [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]]),
    ("4S6", [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]]),
    ("3sigma_h", [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("3sigma_h", [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
    ("3sigma_h", [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]),
];

#[rustfmt::skip]
pub(crate) fn to_matrix3(m: [[f64; 3]; 3]) -> Matrix3<f64> {
    Matrix3::new(
        m[0][0], m[0][1], m[0][2],
        m[1][0], m[1][1], m[1][2],
        m[2][0], m[2][1], m[2][2],
    )
}

#[rustfmt::skip]
pub const POINTGROUP_NAMES: &[&str] = &[
    "C2",
    "C2h",
    "C2v",
    "C3",
    "C3h",
    "C3v",
    "C4",
    "C4h",
    "C4v",
    "C5",
    "C5h",
    "C5v",
    "C6",
    "C6h",
    "C6v",
    "C7",
    "C8",
    "Ci",
    "Cs",
    "D2",
    "D2d",
    "D2h",
    "D3",
    "D3d",
    "D3h",
    "D4",
    "D4d",
    "D4h",
    "D5",
    "D5d",
    "D5h",
    "D6",
    "D6h",
    "D7h",
    "D8h",
    "E",
    "I",
    "Ih",
    "O",
    "Oh",
    "S10",
    "S4",
    "S6",
    "S8",
    "T",
    "Td",
    "Th",
];
