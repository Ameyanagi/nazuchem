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

/// P-73.6 method (1) and the 2020 corrections to P-73.6 and P-74.1.3: a
/// substituted cationic nitrogen prefix is the aminium cation name plus
/// 'yl' or 'ylidene', e.g. (N,N-dimethylmethanaminiumyl)acetate.
#[test]
fn aminiumyl_prefixes_keep_nitrogen_charge_hydrogens_and_attachment_order() {
    let parser = Parser::new().unwrap();
    for (name, expected) in [
        (
            "(N,N-dimethylmethanaminiumyl)acetate",
            "C[N+](C)(C)CC(=O)[O-]",
        ),
        ("(methanaminiumyl)acetate", "C[NH2+]CC(=O)[O-]"),
        ("(N-methylmethanaminiumyl)acetate", "C[NH+](C)CC(=O)[O-]"),
        (
            "2-(N,N-dimethylmethanaminiumyl)propan-2-ide",
            "C[N+](C)(C)[C-](C)C",
        ),
        (
            "(N-ethyl-N-methylpropan-1-aminiumyl)acetate",
            "C(C)[N+](CCC)(C)CC(=O)[O-]",
        ),
        (
            "(N-methylmethanaminiumylidene)acetate",
            "C[N+](C)=CC(=O)[O-]",
        ),
        ("(methanaminiumylidene)acetate", "C[NH+]=CC(=O)[O-]"),
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

// P-72.2.2.1: "-ide" removes a hydron from the named atom. Where a mancude
// ring has no hydrogen there, the name needs added hydrogen (P-72.2.2.1.1,
// pyridin-1(2H)-ide), so the bare names are rejected instead of hydrogenated.
#[test]
fn ide_on_a_mancude_ring_atom_without_hydrogen_is_rejected() {
    let parser = Parser::new().unwrap();
    for name in [
        "pyridin-1-ide",
        "pyrimidin-1-ide",
        "quinolin-1-ide",
        "3H-pyrrol-1-ide",
        "2H-pyrrol-1-ide",
        "pyrazine-1,4-diide",
        "perchloropyridin-1-ide",
    ] {
        let result = parser.parse(name, &ParseOptions::strict());
        assert_eq!(result.status, Status::Failure, "{name}: {result:?}");
        assert!(
            result.message.contains("no hydrogen")
                // P-14.7.1: the cited indicated hydrogen check (#34) can
                // reject 3H-/2H-pyrrol-1-ide first.
                || result.message.contains("without a double bond"),
            "{name}: {}",
            result.message
        );
    }
    for (name, expected) in [
        ("pyridin-1(2H)-ide", "[N-]1CC=CC=C1"),
        ("1,2-dihydropyridin-1-ide", "[N-]1CC=CC=C1"),
        ("pyridin-2-ide", "N1=[C-]C=CC=C1"),
        ("1H-pyrrol-1-ide", "[N-]1C=CC=C1"),
        ("pyrrolide", "[N-]1C=CC=C1"),
        ("2,3-dihydro-1H-pyrrol-1-ide", "[N-]1CCC=C1"),
        ("indol-1-ide", "[N-]1C=CC2=CC=CC=C12"),
        ("2H-pyrazol-2-ide", "N=1[N-]C=CC1"),
        // P-14.7.1: the cited indicated hydrogen, not pyrazole's default 1H,
        // places the hydrogen that the unlocanted "-ide" removes.
        ("2H-pyrazolide", "N=1[N-]C=CC1"),
        ("pyrazolide", "[N-]1N=CC=C1"),
        ("uracil-1-ide", "[N-]1C(=O)NC(=O)C=C1"),
        ("uracil-3-ide", "N1C(=O)[N-]C(=O)C=C1"),
        ("perchloropyrrol-1-ide", "ClC=1[N-]C(=C(C1Cl)Cl)Cl"),
        // Omitted indicated hydrogen may sit on the anionic atom (9H-purine).
        ("purin-9-ide", "N1=CN=C2[N-]C=NC2=C1"),
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

// P-92.5 and P-92.6 (rules 4 and 5): ligands that differ only in the
// configuration of other stereogenic units are ranked at a deeper sphere. C4
// here depends on the r/s centres C3 and C5, and the ring carbons of a
// cyclohexanehexol each depend on the others (compare P-92.5.2.2, example 3).
#[test]
fn centres_ranked_by_other_stereodescriptors_are_assigned() {
    let parser = Parser::new().unwrap();
    for (name, expected) in [
        (
            "(2R,3r,4R,5s,6R)-3,5-bis[(1S)-1-fluoroethyl]-2,6-difluoroheptane-3,4,5-triol",
            "F[C@@H](C)[C@@]([C@@H](C)F)([C@@H]([C@]([C@@H](C)F)(O)[C@H](C)F)O)O",
        ),
        (
            "(1R,2R,3S,4S,5R,6S)-cyclohexane-1,2,3,4,5,6-hexol",
            "[C@@H]1([C@H]([C@H]([C@H]([C@@H]([C@@H]1O)O)O)O)O)O",
        ),
        (
            "(1R,2S,3r,4R,5S,6r)-cyclohexane-1,2,3,4,5,6-hexol",
            "[C@@H]1([C@H]([C@@H]([C@H]([C@H]([C@@H]1O)O)O)O)O)O",
        ),
        (
            "(1R,2R,3R,4R,5S,6S)-1,2,3,4,5,6-hexachlorocyclohexane",
            "Cl[C@H]1[C@@H]([C@H]([C@@H]([C@@H]([C@H]1Cl)Cl)Cl)Cl)Cl",
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
    // No configuration of the six ring carbons gives these labels.
    let result = parser.parse(
        "(1R,2R,3R,4R,5R,6R)-cyclohexane-1,2,3,4,5,6-hexol",
        &ParseOptions::strict(),
    );
    assert_eq!(result.status, Status::Failure, "{result:?}");
    // C3 is ranked by C2 and C4 and is also fixed by its endo descriptor
    // (P-93.5.2.2.1); only the label that agrees with endo is accepted.
    let endo = |c3: &str| format!("(1S,2R,3{c3},4R,5R)-bicyclo[3.2.1]octane-2,3-endo,4-triol");
    let result = parser.parse(&endo("S"), &ParseOptions::strict());
    assert_eq!(result.status, Status::Success, "{}", result.message);
    let result = parser.parse(&endo("R"), &ParseOptions::strict());
    assert_eq!(result.status, Status::Failure, "{result:?}");
    // In warning mode a failed centre keeps the configuration it had.
    let options = ParseOptions {
        warn_rather_than_fail_on_uninterpretable_stereochemistry: true,
        ..ParseOptions::strict()
    };
    let warning = parser.parse("(3R)-ribitol", &options);
    assert_eq!(warning.status, Status::Warning, "{}", warning.message);
    let ribitol = parser.parse("ribitol", &ParseOptions::strict());
    assert_eq!(
        warning.structure.unwrap().semantic_cxsmiles().unwrap(),
        ribitol.structure.unwrap().semantic_cxsmiles().unwrap()
    );
    // Bridgeheads with three equal ring branches are not stereogenic.
    for name in ["(1R,4S)-bicyclo[2.2.2]octane", "(1Xi)-bicyclo[2.2.2]octane"] {
        let result = parser.parse(name, &ParseOptions::strict());
        assert_eq!(result.status, Status::Failure, "{name}: {result:?}");
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
        // P-73.1.1.2: adding hydrons does not supply the missing hydro prefixes.
        "4,5-dihydro-2H-imidazole-1,3-diium",
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
        // P-14.7.2: added hydrogen is not the parent's indicated hydrogen.
        ("4,5-dihydroimidazol-2(3H)-one", "N1C(NCC1)=O"),
        ("4,5-dihydro-1H-imidazol-2(3H)-one", "N1C(NCC1)=O"),
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
        // An =O made from one oxido must not keep the next oxido anionic.
        ("1,1-dioxidothiolane", "O=S1(CCCC1)=O"),
        ("1,1-dioxidothiomorpholine", "O=S1(CCNCC1)=O"),
        ("1-oxido-1-sulfidothiolane", "O=S1(CCCC1)=S"),
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
