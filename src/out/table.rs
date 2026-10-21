//! Terminal tables: columns of fixed or fitted width, filled row by row and printed once.

/// How wide a column is. Every cell is printed after one space.
#[derive(Clone, Copy)]
pub enum Width {
    /// Exactly this many characters.
    Fixed(usize),
    /// The longest cell of the column, header included, plus two.
    Dynamic,
}

pub enum Row {
    /// Printed as is, outside the columns.
    Text(String),
    /// `-` across the table.
    Rule,
    /// `=` across the table.
    DoubleRule,
    /// One cell per column, left to right.
    Cells(Vec<String>),
}
// Makes a row from a slice of Strings
impl<S: Into<String>, const N: usize> From<[S; N]> for Row {
    fn from(cells: [S; N]) -> Self {
        Row::Cells(cells.into_iter().map(Into::into).collect())
    }
}

pub struct Table {
    widths: Vec<Width>,
    rows: Vec<Row>,
}

impl Table {
    pub fn new(widths: &[Width]) -> Self {
        Table {
            widths: widths.to_vec(),
            rows: Vec::new(),
        }
    }

    pub fn row(&mut self, row: Row) {
        if let Row::Cells(cells) = &row {
            debug_assert_eq!(cells.len(), self.widths.len(), "one cell per column");
        }
        self.rows.push(row);
    }

    /// The width of every column, fitting the `Dynamic` ones to their cells.
    fn calc_widths(&self) -> Vec<usize> {
        let cells = || {
            self.rows.iter().filter_map(|row| match row {
                Row::Cells(cells) => Some(cells),
                _ => None,
            })
        };
        let fit = |i: usize| cells().map(|c| c[i].len()).max().unwrap_or(0) + 2;
        (0..self.widths.len())
            .map(|i| match self.widths[i] {
                Width::Fixed(w) => w,
                Width::Dynamic => fit(i),
            })
            .collect()
    }

    /// The table as text, one line per row, without trailing spaces.
    pub fn render(&self) -> String {
        let widths = self.calc_widths();
        let rule_width = widths.iter().sum::<usize>() + widths.len();
        let lines: Vec<String> = self
            .rows
            .iter()
            .map(|row| match row {
                Row::Text(text) => text.clone(),
                Row::Rule => "-".repeat(rule_width),
                Row::DoubleRule => "=".repeat(rule_width),
                Row::Cells(cells) => {
                    let line: String = cells
                        .iter()
                        .zip(&widths)
                        .map(|(cell, w)| format!(" {cell:<w$}"))
                        .collect();
                    line.trim_end().to_string()
                }
            })
            .collect();
        lines.join("\n")
    }

    pub fn print(&self, silent: bool) {
        say!(silent, "{}", self.render());
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn renders_fixed_and_dynamic_columns_between_rules() {
        let mut table = Table::new(&[Width::Dynamic, Width::Fixed(8)]);
        table.row(Row::Text("Title".to_string()));
        table.row(Row::DoubleRule);
        table.row(Row::from(["Symbol", "S"]));
        table.row(Row::Rule);
        table.row(Row::from(["OC-6".to_string(), format!("{:>7.3}", 0.0)]));
        table.row(Row::from(["TPR-6", " 16.500"]));
        table.row(Row::Rule);

        // "Symbol" (6) + 2, then 8, plus one space per column.
        let expected: [&str; 7] = [
            "Title",
            &"=".repeat(18),
            " Symbol   S",
            &"-".repeat(18),
            " OC-6       0.000",
            " TPR-6     16.500",
            &"-".repeat(18),
        ];
        assert_eq!(table.render(), expected.join("\n"));
    }
}
