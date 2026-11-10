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

/// Issue #26: chalcogen peroxy names place each element where the name says.
#[test]
fn chalcogen_peroxy_names_keep_their_element_order() {
    let parser = Parser::new().unwrap();
    for (name, expected) in [
        // P-65.1.5.3: the acid locant names the atom on C, then the acidic atom.
        ("ethane(thioperoxoic) OS-acid", "C(C)(=O)OS"),
        ("benzenecarbo(thioperoxoic) SO-acid", "C1(=CC=CC=C1)C(=O)SO"),
        (
            "ethaneseleno(tellurothioperoxoic) STe-acid",
            "C(C)(S[TeH])=[Se]",
        ),
        ("prop-2-ene(dithioperoxoic) acid", "C(C=C)(=O)SS"),
        ("methanesulfono(thioperoxoic) OS-acid", "CS(=O)(=O)OS"),
        // P-63.4.2.1, Table 6.1.
        ("methane-SO-thioperoxol", "CSO"),
        ("propane-1-OS-thioperoxol", "C(CC)OS"),
        ("methane-SSe-selenothioperoxol", "CS[SeH]"),
        ("methane-TeSe-selenotelluroperoxol", "C[Te][SeH]"),
        // P-66.6.3, Table 6.4.
        ("ethanetellanal", "C(C)=[Te]"),
        // P-65.7.4: the element locant names the bridge atom bonded to that acid.
        ("acetic thioperoxyanhydride", "C(C)(=O)SOC(C)=O"),
        (
            "S-acetic O-propanoic thioperoxyanhydride",
            "C(CC)(=O)OSC(C)=O",
        ),
        (
            "O-acetic S-propanoic thioperoxyanhydride",
            "C(CC)(=O)SOC(C)=O",
        ),
        // P-63.3.2, method (2).
        ("O-ethyl S-methyl thioperoxide", "CSOCC"),
        ("S-cyclohexyl O-methyl thioperoxide", "COSC1CCCCC1"),
        ("Se-ethyl S-methyl selenothioperoxide", "CS[Se]CC"),
        ("dimethyl selenoperoxide", "CO[Se]C"),
        // A single leading locant still belongs to its substituent.
        (
            "acetic S-methylcysteine anhydride",
            "CSC[C@H](N)C(=O)OC(C)=O",
        ),
        (
            "O-methylserine S-methylcysteine anhydride",
            "CSC[C@H](N)C(=O)OC([C@@H](N)COC)=O",
        ),
        // Tellanal keeps its aldehydic hydrogen like ethanal.
        ("perchloroethanetellanal", "ClC(C=[Te])(Cl)Cl"),
    ] {
        let result = parser.parse(name, &ParseOptions::strict());
        assert_eq!(result.status, Status::Success, "{name}: {}", result.message);
        assert_eq!(
            result.structure.unwrap().semantic_cxsmiles().unwrap(),
            expected,
            "{name}"
        );
    }
    for name in [
        "S-acetic S-propanoic thioperoxyanhydride",
        "S-ethyl S-methyl thioperoxide",
        "O-ethyl O-methyl ether",
        // Without locants the bridge could face either way.
        "acetic propanoic thioperoxyanhydride",
        "ethyl methyl thioperoxide",
        "2-methylbutanedioic thioperoxyanhydride",
    ] {
        let result = parser.parse(name, &ParseOptions::strict());
        assert_eq!(result.status, Status::Failure, "{name}");
    }
}
