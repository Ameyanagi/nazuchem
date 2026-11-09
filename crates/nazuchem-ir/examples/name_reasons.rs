// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | crates/nazuchem-ir/examples/name_reasons.rs | binding vocabulary checks

fn main() {
    println!(
        "{}",
        serde_json::to_string(nazuchem_ir::NAME_REASONS).expect("reason vocabulary serializes")
    );
}
