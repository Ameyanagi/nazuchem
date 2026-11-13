// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | crates/nazuchem/tests/phosphorus_acid_prefixes.rs | issue #160

#![cfg(all(feature = "n2s", feature = "s2n"))]

//! Issue #160: prefixes for phosphorus, arsenic and antimony acid groups follow
//! the IUPAC 2013 recommendations, are read back as the structure they came
//! from, and are reported as preferred only when they are.

use nazuchem::{
    Budget, NameReason, NameStatus, name_to_smiles, s2n::canonical_smiles, smiles_to_name,
};

fn assert_round_trip(cases: &[(&str, &str, NameStatus)]) {
    let budget = Budget::default();
    for &(smiles, expected, status) in cases {
        let name = smiles_to_name(smiles).unwrap();
        assert_eq!(name.text, expected, "{smiles}");
        assert_eq!(name.status, status, "{smiles}");
        let parsed = name_to_smiles(&name.text).unwrap_or_else(|e| panic!("{expected}: {e:?}"));
        assert_eq!(
            canonical_smiles(&parsed, &budget),
            canonical_smiles(smiles, &budget),
            "{smiles} -> {expected} -> {parsed}"
        );
    }
}

#[test]
fn acid_groups_use_the_retained_and_acyl_prefixes() {
    let preferred = NameStatus::Preferred;
    assert_round_trip(&[
        // P-67.1.4.1.1.1: phosphono, arsono and stibono are retained;
        // P-67.1.5.1: phosphonoacetic acid.
        ("OC(=O)CP(=O)(O)O", "phosphonoacetic acid", preferred),
        ("OC(=O)C[As](=O)(O)O", "arsonoacetic acid", preferred),
        ("OC(=O)C[Sb](=O)(O)O", "stibonoacetic acid", preferred),
        (
            "OC(=O)C(P(=O)(O)O)P(=O)(O)O",
            "diphosphonoacetic acid",
            preferred,
        ),
        (
            "OS(=O)(=O)CP(=O)(O)O",
            "phosphonomethanesulfonic acid",
            preferred,
        ),
        // P-67.1.4.1.1.2: PH2(O)- is phosphinoyl and AsH2(O)- is arsinoyl;
        // oxophosphanyl is -P=O (P-67.1.4.1.1.6).
        ("OC(=O)C[PH2]=O", "phosphinoylacetic acid", preferred),
        ("OC(=O)C[AsH2]=O", "arsinoylacetic acid", preferred),
        ("OC(=O)CP=O", "oxophosphanylacetic acid", preferred),
        // P-67.1.4.1.1.3: methylphosphonoyl, ethylstibinoyl and
        // methyl(phenyl)arsinoyl are substituted acyl groups.
        (
            "OC(=O)C[PH](C)=O",
            "(methylphosphinoyl)acetic acid",
            preferred,
        ),
        (
            "OC(=O)CP(C)(C)=O",
            "(dimethylphosphinoyl)acetic acid",
            preferred,
        ),
        (
            "O=C(O)c1ccc(cc1)P(C)(C)=O",
            "4-(dimethylphosphinoyl)benzoic acid",
            preferred,
        ),
        (
            "OC(=O)C[SbH](CC)=O",
            "(ethylstibinoyl)acetic acid",
            preferred,
        ),
        (
            "OC(=O)C[As](C)(c1ccccc1)=O",
            "[methyl(phenyl)arsinoyl]acetic acid",
            preferred,
        ),
        // P-67.1.4.1.1.5: hydroxy, alkoxy and sulfanyl are concatenated with
        // phosphoryl, as in hydroxy(sulfanyl)phosphoryl;
        // P-67.1.4.1.3: (dimethoxyphosphoryl)sulfanyl.
        (
            "OC(=O)CP(=O)(OC)OC",
            "(dimethoxyphosphoryl)acetic acid",
            preferred,
        ),
        (
            "OC(=O)CP(=O)(O)S",
            "[hydroxy(sulfanyl)phosphoryl]acetic acid",
            preferred,
        ),
    ]);
}

#[test]
fn prefixes_that_are_not_the_preferred_acyl_form_say_so() {
    let acyl = NameStatus::NotPreferred(NameReason::AcylPrefix);
    assert_round_trip(&[
        // P-67.1.4.1.1.5: the preferred prefixes are
        // hydroxy(methylphosphonoyl) and its methoxy and phenyl analogues.
        (
            "OC(=O)CP(C)(=O)O",
            "[hydroxy(methyl)phosphoryl]acetic acid",
            acyl,
        ),
        (
            "OC(=O)CP(C)(=O)OC",
            "[methoxy(methyl)phosphoryl]acetic acid",
            acyl,
        ),
        (
            "OC(=O)CP(=O)(O)c1ccccc1",
            "[hydroxy(phenyl)phosphoryl]acetic acid",
            acyl,
        ),
        // P-67.1.4.1.1.4: the preferred prefixes derive from the
        // phosphorochloridic and phosphoramidic acids, as in
        // phosphorodichloridoyl, not dichlorophosphoryl.
        (
            "OC(=O)CP(=O)(O)Cl",
            "[chloro(hydroxy)phosphoryl]acetic acid",
            acyl,
        ),
        (
            "OC(=O)CP(=O)(O)N",
            "[amino(hydroxy)phosphoryl]acetic acid",
            acyl,
        ),
        // P-67.1.4.1.1.5: hydrogen of phosphinoyl is not substituted in a
        // preferred prefix; hydroxy is concatenated with phosphonoyl, as in
        // hydroxyazonoyl and selanylphosphonoyl.
        ("OC(=O)C[PH](=O)O", "(hydroxyphosphinoyl)acetic acid", acyl),
        // P-68.3.2.3.1: with no principal characteristic group the
        // phosphorus is the parent, a lambda5-phosphanone.
        ("C[PH2]=O", "phosphinoylmethane", acyl),
        ("C[PH](C)=O", "(methylphosphinoyl)methane", acyl),
    ]);
}

#[test]
fn unsubstituted_stibinoyl_is_declined() {
    // P-67.1.4.1.1.2 names SbH2(O)- stibinoyl, which the parser reads as
    // -Sb=O; no name is better than one for another structure.
    assert!(smiles_to_name("OC(=O)C[SbH2]=O").is_err());
}

#[test]
fn blue_book_prefix_names_are_read() {
    let budget = Budget::default();
    for (name, smiles) in [
        // P-67.1.5.1.
        ("phosphonoacetic acid", "OC(=O)CP(=O)(O)O"),
        (
            "(4-arsonobutyl)phosphonic acid",
            "OP(=O)(O)CCCC[As](=O)(O)O",
        ),
        // P-67.1.4.1.3.
        ("(phosphonooxy)acetic acid", "OC(=O)COP(=O)(O)O"),
        // P-67.1.4.1.1.5.
        (
            "4,4'-(hydroxyphosphoryl)dibenzoic acid",
            "OC(=O)c1ccc(cc1)P(=O)(O)c1ccc(cc1)C(=O)O",
        ),
        // P-67.1.4.1.1.6: -P(O)2 is dioxo-lambda5-phosphanyl.
        ("(dioxo-lambda^5-phosphanyl)acetic acid", "OC(=O)CP(=O)=O"),
    ] {
        let parsed = name_to_smiles(name).unwrap_or_else(|e| panic!("{name}: {e:?}"));
        assert_eq!(
            canonical_smiles(&parsed, &budget),
            canonical_smiles(smiles, &budget),
            "{name} -> {parsed}"
        );
    }
}
