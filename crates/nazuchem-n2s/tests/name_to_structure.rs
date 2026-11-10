// SPDX-License-Identifier: MIT
// Origin: adapted
// Source: OPSIN | opsin-core/src/main/java/uk/ac/cam/ch/wwmm/opsin/NameToStructure.java | 42e87ac6efb1bd9f578986f5ea7ee1faa8e93c91
// Source: OPSIN | opsin-core/src/test/java/uk/ac/cam/ch/wwmm/opsin/NameToStructureTest.java | 42e87ac6efb1bd9f578986f5ea7ee1faa8e93c91
// Source: ReShiki | crates/opsin/tests/name_to_structure.rs | 5c3d90c2852227ad644ac07ad91300f032c697f5
// Copyright 2017 Daniel Lowe
// Changes: Rust graph/tree construction and regression harness around source-derived cases/assertions; local regressions retained.

//! Upstream NameToStructureTest additions and complete-result regressions.
use nazuchem_n2s::{ParseOptions, Parser, Status};

#[test]
fn unicode_superscript_bridge_locants_round_trip_to_the_ascii_structure() {
    let parser = Parser::new().unwrap();
    let ascii = parser.parse("tricyclo[3.2.1.0^{2,4}]octane", &ParseOptions::strict());
    let formatted = parser.parse("tricyclo[3.2.1.0²,⁴]octane", &ParseOptions::strict());
    assert_eq!(ascii.status, Status::Success, "{}", ascii.message);
    assert_eq!(formatted.status, Status::Success, "{}", formatted.message);
    assert_eq!(
        formatted.structure.unwrap().semantic_cxsmiles().unwrap(),
        ascii.structure.unwrap().semantic_cxsmiles().unwrap()
    );
}

#[test]
fn isotope_specification_between_fusion_components_fails_cleanly() {
    let result = Parser::new()
        .unwrap()
        .parse("pyrido[2,3-b](513C)pyrazine", &ParseOptions::strict());
    assert_eq!(result.status, Status::Failure);
    assert!(
        result
            .message
            .starts_with("Unable to combine ring components"),
        "Expected a diagnosable failure but was: {}",
        result.message
    );
}

#[test]
fn unlocanted_hydro_ignores_suffix_bonds_in_parent_valency() {
    let result = Parser::new()
        .unwrap()
        .parse("hexahydropyridin-2-one", &ParseOptions::strict());
    assert_eq!(result.status, Status::Success, "{}", result.message);
    assert_eq!(
        result.structure.unwrap().semantic_cxsmiles().unwrap(),
        "N1C(CCCC1)=O"
    );
}

#[test]
fn peroxycyanate_esters_preserve_both_oxygens() {
    let parser = Parser::new().unwrap();
    for (name, expected) in [
        ("methyl peroxycyanate", "COOC#N"),
        ("phenyl peroxycyanate", "C1(=CC=CC=C1)OOC#N"),
    ] {
        let result = parser.parse(name, &ParseOptions::strict());
        assert_eq!(result.status, Status::Success, "{name}: {}", result.message);
        assert_eq!(
            result.structure.unwrap().semantic_cxsmiles().unwrap(),
            expected,
            "{name}"
        );
    }
}

// Issue #33. P-72.6.2: oxido is –O⁻; P-72.6.1: –S(=O)O⁻ is sulfinato. On an S/P
// that already has =O/=S, oxido stays anionic instead of becoming a second oxo.
#[test]
fn oxido_on_sulfinyl_like_groups_keeps_its_charge() {
    let parser = Parser::new().unwrap();
    for (name, expected) in [
        (
            "[methyl(oxidosulfinyl)amino]benzene",
            "CN(S(=O)[O-])C1=CC=CC=C1",
        ),
        // The structure-to-name spelling of the same anion.
        (
            "[methyl(sulfinato)amino]benzene",
            "CN(S(=O)[O-])C1=CC=CC=C1",
        ),
        ("oxidosulfinylbenzene", "[O-]S(=O)C1=CC=CC=C1"),
        ("oxidophosphinylbenzene", "[O-]P(=O)C1=CC=CC=C1"),
        ("oxidosulfinothioylbenzene", "[O-]S(=S)C1=CC=CC=C1"),
        // Without an existing =O/=S the upstream oxido-to-oxo reading is kept.
        ("1-oxidothiolane", "O=S1CCCC1"),
        ("dimethyl(oxido)sulfur", "CS(=O)C"),
    ] {
        let result = parser.parse(name, &ParseOptions::strict());
        assert_eq!(result.status, Status::Success, "{name}: {}", result.message);
        assert_eq!(
            result.structure.unwrap().semantic_cxsmiles().unwrap(),
            expected,
            "{name}"
        );
    }
}
