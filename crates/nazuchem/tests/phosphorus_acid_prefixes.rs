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
fn acid_derived_prefixes_preserve_hydrogen_and_add_ligands() {
    let preferred = NameStatus::Preferred;
    assert_round_trip(&[
        // P-67.1.4.1.1.5: the methyl group substitutes H within the
        // phosphonoyl group, and the chalcogen group fills a free valence.
        (
            "OC(=O)CP(C)(=O)O",
            "[hydroxy(methylphosphonoyl)]acetic acid",
            preferred,
        ),
        (
            "OC(=O)CP(C)(=O)OC",
            "[methoxy(methylphosphonoyl)]acetic acid",
            preferred,
        ),
        (
            "OC(=O)CP(=O)(O)c1ccccc1",
            "[hydroxy(phenylphosphonoyl)]acetic acid",
            preferred,
        ),
        (
            "OC(=O)C[PH](=O)O",
            "(hydroxyphosphonoyl)acetic acid",
            preferred,
        ),
        (
            "OC(=O)C[PH](=O)[SeH]",
            "(selanylphosphonoyl)acetic acid",
            preferred,
        ),
        (
            "OC(=O)C[As](C)(=O)O",
            "[hydroxy(methylarsonoyl)]acetic acid",
            preferred,
        ),
        (
            "OC(=O)C[As](C)(=O)OC",
            "[methoxy(methylarsonoyl)]acetic acid",
            preferred,
        ),
        (
            "OC(=O)C[SbH](=O)O",
            "(hydroxystibonoyl)acetic acid",
            preferred,
        ),
        // P-67.1.4.1.1.4: halide and amide replacement changes the acyl
        // stem. Hydroxy then fills one of the free valences.
        (
            "OC(=O)CP(=O)(O)Cl",
            "(hydroxyphosphorochloridoyl)acetic acid",
            preferred,
        ),
        (
            "OC(=O)CP(=O)(Cl)Cl",
            "phosphorodichloridoylacetic acid",
            preferred,
        ),
        (
            "OC(=O)CP(=O)(O)N",
            "(hydroxyphosphoramidoyl)acetic acid",
            preferred,
        ),
        (
            "OC(=O)CP(=S)(O)O",
            "(dihydroxyphosphorothioyl)acetic acid",
            preferred,
        ),
        (
            "OC(=O)CP(=[Se])(OC)OC",
            "(dimethoxyphosphoroselenoyl)acetic acid",
            preferred,
        ),
        (
            "OC(=O)CCNP(=S)(O)S",
            "3-{[hydroxy(sulfanyl)phosphorothioyl]amino}propanoic acid",
            preferred,
        ),
        // P-67.1.4.1.1.2 and P-67.1.4.1.1.6 distinguish these Sb prefixes.
        ("OC(=O)C[SbH2]=O", "stibinoylacetic acid", preferred),
        ("OC(=O)C[Sb]=O", "oxostibanylacetic acid", preferred),
        (
            "OC(=O)CP(=O)=O",
            "(dioxo-lambda^5-phosphanyl)acetic acid",
            preferred,
        ),
        // P-67.1.5.1: the phosphorus acid is senior to the arsenic acid.
        (
            "OP(=O)(O)CCCC[As](=O)(O)O",
            "(4-arsonobutyl)phosphonic acid",
            preferred,
        ),
    ]);
}

#[test]
fn phosphorus_heterones_with_hydrogen_use_lambda_five() {
    let preferred = NameStatus::Preferred;
    assert_round_trip(&[
        // P-64.4.1 and P-68.3.2.3.1: the oxygen is a suffix on the
        // phosphorus parent, which retains any remaining hydrogens.
        ("[PH3]=O", "lambda5-phosphanone", preferred),
        ("C[PH2]=O", "methyl-lambda5-phosphanone", preferred),
        ("C[PH](C)=O", "dimethyl-lambda5-phosphanone", preferred),
        (
            "OC[PH](C)=O",
            "(hydroxymethyl)(methyl)-lambda5-phosphanone",
            preferred,
        ),
    ]);
}

#[test]
fn the_heterone_shortcut_preserves_charge() {
    for smiles in ["[PH3+][O-]", "C[PH2+][O-]", "C[PH+](C)[O-]"] {
        assert!(smiles_to_name(smiles).is_err(), "{smiles}");
    }

    assert_round_trip(&[(
        "[Cl-].[PH3+]C",
        "phosphaniumylmethane chloride",
        NameStatus::Unverified(NameReason::CationParent),
    )]);
}

#[test]
fn isotopic_phosphorus_heterones_need_locants() {
    for smiles in ["C[PH]([2H])=O", "C[PH2]=[18O]"] {
        assert!(smiles_to_name(smiles).is_err(), "{smiles}");
    }

    assert_round_trip(&[(
        "C[PH]([13CH3])=O",
        "(methylphosphinoyl)(1-13C)methane",
        NameStatus::NotPreferred(NameReason::AcylPrefix),
    )]);
}

#[test]
fn competing_ketones_use_the_ordinary_parent_selection() {
    let acyl = NameStatus::NotPreferred(NameReason::AcylPrefix);
    assert_round_trip(&[
        // P-64.7.2: the dione can express both carbonyls as suffixes.
        (
            "C[PH](=O)CC(=O)CC(=O)C",
            "1-(methylphosphinoyl)pentane-2,4-dione",
            acyl,
        ),
        (
            "C[PH](=O)CC(=O)C",
            "1-(methylphosphinoyl)propan-2-one",
            acyl,
        ),
    ]);
}

#[test]
fn prefixes_that_are_not_the_preferred_acyl_form_say_so() {
    let acyl = NameStatus::NotPreferred(NameReason::AcylPrefix);
    assert_round_trip(&[
        // P-64.7.1, P-68.3.2.3.1: heterones are senior to hydroxy compounds
        // and amines, so without a more senior class the phosphorus is the
        // parent, a lambda5-phosphanone.
        ("OCCP(C)(C)=O", "2-(dimethylphosphinoyl)ethan-1-ol", acyl),
        ("NCCP(C)(C)=O", "2-(dimethylphosphinoyl)ethan-1-amine", acyl),
        ("OCC[As](C)(C)=O", "2-(dimethylarsinoyl)ethan-1-ol", acyl),
    ]);
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
