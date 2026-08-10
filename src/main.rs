mod xyz;

fn main() {
    let result = xyz::parse_xyz("sample.xyz", true);

    match result {
        Ok(structure) => println!("{:?}", structure),
        Err(e) => eprintln!("error: {}", e),
    }
}
