// From vergen documentation
extern crate vergen;
use std::process::Command;

use vergen::{ConstantsFlags, generate_cargo_keys};

fn main() {
    let flags = ConstantsFlags::all();
    generate_cargo_keys(flags).expect("Must be able to generate cargo keys");
    // Based on 
    // https://vallentin.io/2019/06/06/versioning
    // https://stackoverflow.com/questions/43753491/include-git-commit-hash-as-string-into-rust-program
    // https://unix.stackexchange.com/questions/155046/determine-if-git-working-directory-is-clean-from-a-script
    let out = Command::new("git")
        .arg("status")
        .arg("--porcelain")
        .output()
        .expect("Can run git command successfully");

    if out.stdout.is_empty() {
        println!("cargo:rustc-env=WD_IS_CLEAN=true");
    } else {
        println!("cargo:rustc-env=WD_IS_CLEAN=false");
    }
}
