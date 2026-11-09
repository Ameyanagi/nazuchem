// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | crates/nazuchem/tests/ir.rs | Phase B step 8a

#![cfg(feature = "n2s")]

use nazuchem::{ErrorKind, Molecule, n2s, name_to_molecule, name_to_smiles};

#[test]
fn facade_returns_ethanol_ir_with_parent_locants_and_json_roundtrip() {
    let molecule = name_to_molecule("ethanol").unwrap();
    molecule.validate().unwrap();
    assert_eq!(molecule.atoms.len(), 3);
    assert_eq!(molecule.bonds.len(), 2);
    assert_eq!(
        molecule.displayed_locants(),
        [Some("1".into()), Some("2".into()), None]
    );
    assert_eq!(
        molecule.semantic_cxsmiles().unwrap(),
        name_to_smiles("ethanol").unwrap()
    );
    let restored = Molecule::from_json(&molecule.to_json().unwrap()).unwrap();
    assert_eq!(restored, molecule);
    assert_eq!(
        restored.to_cxsmiles().unwrap(),
        molecule.to_cxsmiles().unwrap()
    );
}

#[test]
fn strict_ir_facade_preserves_warning_and_parse_error_classification() {
    let parser = n2s::Parser::new().unwrap();
    for (name, status, kind) in [
        ("(+)-lactic acid", n2s::Status::Warning, ErrorKind::Warning),
        (
            "not a chemical name",
            n2s::Status::Failure,
            ErrorKind::Parse,
        ),
    ] {
        let parsed = parser.parse(name, &n2s::ParseOptions::strict());
        assert_eq!(parsed.status, status);
        let error = name_to_molecule(name).unwrap_err();
        assert_eq!(error.kind(), kind);
        assert_eq!(error.to_string(), parsed.message);
    }
}

#[test]
fn axial_ir_rejection_is_typed_and_classified_as_facade_serialization_error() {
    let name = "(1M)-[1,1'-binaphthalene]-2,2'-diol";
    let parsed = n2s::Parser::new()
        .unwrap()
        .parse(name, &n2s::ParseOptions::strict());
    assert_eq!(parsed.status, n2s::Status::Success, "{}", parsed.message);
    let unsupported = parsed.structure.unwrap().to_ir().unwrap_err();
    assert!(matches!(unsupported, n2s::IrError::UnsupportedAxialStereo));
    let error = name_to_molecule(name).unwrap_err();
    assert_eq!(error.kind(), ErrorKind::Serialize);
    assert_eq!(error.to_string(), unsupported.to_string());
}
