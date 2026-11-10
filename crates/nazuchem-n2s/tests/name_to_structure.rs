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

#[test]
fn cited_indicated_hydrogen_rejects_hydro_prefixes_that_strand_an_atom() {
    // P-14.7.1, P-31.2.2: once indicated hydrogen is cited, hydro prefixes that
    // leave a ring atom without a double-bond partner describe no structure.
    let parser = Parser::new().unwrap();
    for name in [
        "3a,4,7,7a-tetrahydro-2H-isoindole",
        "3a,4,7,7a-tetrahydro-2H-isoindol-2-ium",
        "4,5-dihydro-2H-imidazole",
        "4,6-dihydro-1H-indene",
        "hexahydro-1H,6H-pyrrolo[3,4-b]pyridine",
        // The spiro junction replaces the cited atom.
        "spiro[4,5-dihydro-2H-imidazole-2,1'-cyclopentane]",
        "spiro[cyclopentane-1,2'-4,5-dihydro-2H-imidazole]",
    ] {
        let result = parser.parse(name, &ParseOptions::strict());
        assert_eq!(result.status, Status::Failure, "{name}");
        assert!(
            result.message.contains("without a double bond"),
            "{name}: {}",
            result.message
        );
    }
    // Valid combinations; P-58.2.1.1 lets indicated hydrogen be omitted.
    for (name, expected) in [
        (
            "2,3,3a,4,7,7a-hexahydro-1H-isoindol-2-ium",
            "C1[NH2+]CC2CC=CCC12",
        ),
        ("3a,4,7,7a-tetrahydro-1H-isoindole", "C1N=CC2CC=CCC12"),
        ("1,3-dihydroisoindole", "C1NCC2=CC=CC=C12"),
        ("1,3-dihydro-2H-isoindole", "C1NCC2=CC=CC=C12"),
        ("octahydro-1H-pyrrolo[3,4-b]pyridine", "N1C2C(CCC1)CNC2"),
        ("pyrimidine-2,4(1H,3H)-dione", "N1C(NC(C=C1)=O)=O"),
        // Indicated hydrogen cited in one spiro component only.
        (
            "spiro[1,3-dihydroisoindole-1,2'-2,3-dihydro-1H-indene]",
            "C1C2(CC3=CC=CC=C13)NCC1=CC=CC=C12",
        ),
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
