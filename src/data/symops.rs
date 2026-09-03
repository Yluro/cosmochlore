/// Look up the `operation` symmetry-operation set by name.
pub fn get_operation(name: &str) -> Option<&'static [(&'static str, [[f64; 3]; 3])]> {
    match name {
        "C2" => Some(OPERATION_C2),
        "C3" => Some(OPERATION_C3),
        "C3_2" => Some(OPERATION_C3_2),
        "S4" => Some(OPERATION_S4),
        "i" => Some(OPERATION_I),
        "sigma_h" => Some(OPERATION_SIGMA_H),
        "sigma_v" => Some(OPERATION_SIGMA_V),
        _ => None,
    }
}

pub static OPERATION_C2: &[(&str, [[f64; 3]; 3])] = &[
    ("C2", [[-1.0, 0.0, 0.0], [-0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
];

pub static OPERATION_C3: &[(&str, [[f64; 3]; 3])] = &[
    ("C3", [[-0.5, 0.866, 0.0], [-0.866, -0.5, 0.0], [0.0, 0.0, 1.0]]),
];

pub static OPERATION_C3_2: &[(&str, [[f64; 3]; 3])] = &[
    ("C3_2", [[-0.5, -0.866, 0.0], [0.866, -0.5, 0.0], [0.0, 0.0, 1.0]]),
];

pub static OPERATION_I: &[(&str, [[f64; 3]; 3])] = &[
    ("sigma_h", [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]),
];

pub static OPERATION_S4: &[(&str, [[f64; 3]; 3])] = &[
    ("S4", [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]),
];

pub static OPERATION_SIGMA_H: &[(&str, [[f64; 3]; 3])] = &[
    ("sigma_v", [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]),
];

pub static OPERATION_SIGMA_V: &[(&str, [[f64; 3]; 3])] = &[
    ("sigma_v", [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
];
