// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | crates/nazuchem/tests/phosphorus_acid_anions.rs | issue #150

#![cfg(all(feature = "n2s", feature = "s2n"))]

//! Issue #150: anions, acid salts and partial esters of the organyl-substituted
//! phosphorus, arsenic and antimony oxoacids are named from the acid after the
//! IUPAC 2013 recommendations and read back as the structure they came from.

use nazuchem::{Budget, name_to_smiles, s2n::canonical_smiles, smiles_to_name};

fn assert_names_round_trip(cases: &[(&str, &str)]) {
    let budget = Budget::default();
    for &(smiles, expected) in cases {
        let name = smiles_to_name(smiles).unwrap().text;
        assert_eq!(name, expected, "{smiles}");
        let parsed = name_to_smiles(&name).unwrap_or_else(|e| panic!("{name}: {e:?}"));
        assert_eq!(
            canonical_smiles(&parsed, &budget),
            canonical_smiles(smiles, &budget),
            "{smiles} -> {name} -> {parsed}"
        );
    }
}

#[test]
fn organyl_pnictogen_acid_anions_are_named_from_the_acid() {
    assert_names_round_trip(&[
        // P-67.1.1.1 (phosphinic, phosphinous, phosphonic, phosphonous acid),
        // P-67.1.1.2 (organyl groups substitute the hydrogen on phosphorus),
        // P-72.2.2.2.1.1 (`ic acid` to `ate`, `ous acid` to `ite`, as in
        // dibenzylphosphinite).
        ("CP[O-]", "methylphosphinite"),
        ("C[PH](=O)[O-]", "methylphosphinate"),
        ("CP(=O)[O-]", "methylphosphinate"),
        ("CP(C)[O-]", "dimethylphosphinite"),
        ("CP(C)(=O)[O-]", "dimethylphosphinate"),
        ("CCP(C)(=O)[O-]", "ethyl(methyl)phosphinate"),
        ("[O-]P(C)c1ccccc1", "methyl(phenyl)phosphinite"),
        ("c1ccccc1[PH](=O)[O-]", "phenylphosphinate"),
        ("CP([O-])[O-]", "methylphosphonite"),
        ("CP(=O)([O-])[O-]", "methylphosphonate"),
        ("c1ccccc1P(=O)([O-])[O-]", "phenylphosphonate"),
        // P-72.2.2.2.1.2: hydrogen phenylphosphonate.
        ("CP(=O)(O)[O-]", "hydrogen methylphosphonate"),
        ("CP(O)[O-]", "hydrogen methylphosphonite"),
        // P-67.1.3.2: ester groups precede the anion name.
        ("CP(=O)(OC)[O-]", "methyl methylphosphonate"),
        ("CP(OC)[O-]", "methyl methylphosphonite"),
        // P-67.1.1.1 arsinic, arsinous, arsonic, stibinic and stibinous acid;
        // P-67.1.3.1: potassium dimethylarsinate.
        ("C[AsH](=O)[O-]", "methylarsinate"),
        ("C[As](C)(=O)[O-]", "dimethylarsinate"),
        ("C[AsH][O-]", "methylarsinite"),
        ("C[As](=O)([O-])[O-]", "methylarsonate"),
        ("C[SbH](=O)[O-]", "methylstibinate"),
        ("C[Sb](C)[O-]", "dimethylstibinite"),
        // P-67.1.3.1: cations, then `hydrogen`, then the anion; disodium
        // methylphosphonite and potassium dimethylarsinate are its examples.
        ("[Na+].C[PH](=O)[O-]", "sodium methylphosphinate"),
        ("[Na+].[Na+].CP([O-])[O-]", "disodium methylphosphonite"),
        ("[K+].C[As](C)(=O)[O-]", "potassium dimethylarsinate"),
        ("[Na+].CP(=O)(O)[O-]", "sodium hydrogen methylphosphonate"),
        ("[Na+].CP(=O)(OC)[O-]", "sodium methyl methylphosphonate"),
    ]);
}

#[test]
fn organyl_arsenic_and_antimony_acids_are_named_like_the_phosphorus_acids() {
    assert_names_round_trip(&[
        // P-67.1.1.2: diphenylarsinous acid, phenylstibonous acid; the
        // pentavalent acids follow methyl(phenyl)arsinic acid and
        // ethylstibinic acid (P-67.1.4.1.1.3).
        ("C[As](=O)(O)O", "methylarsonic acid"),
        ("C[As](C)(=O)O", "dimethylarsinic acid"),
        ("C[AsH](=O)O", "methylarsinic acid"),
        ("C[Sb](=O)(O)O", "methylstibonic acid"),
        // P-67.1.3.2: diethyl methylphosphonate.
        ("C[As](=O)(OC)OC", "dimethyl methylarsonate"),
        // The phosphorus acids are unchanged.
        ("CP(=O)(O)O", "methylphosphonic acid"),
        ("C[PH](=O)O", "methylphosphinic acid"),
        ("CP(=O)(OC)O", "methyl hydrogen methylphosphonate"),
    ]);
}

#[test]
fn pnictogen_prefixes_cite_their_oxygen_ligands() {
    assert_names_round_trip(&[
        // P-16.5.1.1: the compound prefix is enclosed, so that oxido stays
        // on phosphorus.
        ("[O-][PH](=O)CC(=O)[O-]", "[oxido(oxo)phosphanyl]acetate"),
        ("[O-][AsH](=O)CC(=O)[O-]", "[oxido(oxo)arsanyl]acetate"),
        // P-67.1.4.1.1.5 (concatenation with the acyl prefix, as in
        // dimethoxyphosphoroselenoyl) with oxido from P-72.6.2.
        (
            "[O-]P(=O)(OC)CC(=O)[O-]",
            "[methoxy(oxido)phosphoryl]acetate",
        ),
        (
            "[O-][As](=O)(OC)CC(=O)[O-]",
            "[methoxy(oxido)arsoryl]acetate",
        ),
        // P-67.1.4.1.1.2 arsoryl; P-67.1.4.1.1.6 hydroxyarsanyl and
        // dimethoxyphosphanyl.
        ("O[As](C)(=O)CC(=O)O", "[hydroxy(methyl)arsoryl]acetic acid"),
        (
            "CO[As](C)(=O)CC(=O)O",
            "[methoxy(methyl)arsoryl]acetic acid",
        ),
        ("O[AsH]CC(=O)O", "(hydroxyarsanyl)acetic acid"),
        ("O[As](C)CC(=O)O", "[hydroxy(methyl)arsanyl]acetic acid"),
        ("CO[As](C)CC(=O)O", "[methoxy(methyl)arsanyl]acetic acid"),
        ("O[Sb](C)CC(=O)O", "[hydroxy(methyl)stibanyl]acetic acid"),
        // P-72.6.2 oxido on the substitutive prefixes of P-67.1.4.1.1.6.
        ("[O-][AsH]CC(=O)[O-]", "(oxidoarsanyl)acetate"),
        ("[O-][As](C)CC(=O)[O-]", "[methyl(oxido)arsanyl]acetate"),
        ("[O-][SbH]CC(=O)[O-]", "(oxidostibanyl)acetate"),
    ]);
}

#[test]
fn phosphorus_anions_without_a_readable_name_are_declined() {
    // `oxidophosphanyl` is read as the phosphane oxide P=O, and the thio
    // analogues need the phosphinothioate names of P-67.1.2.4.1.1, which
    // are not generated; none of these may be named as an acyl group that
    // drops the charge.
    for smiles in [
        "[O-]PCC(=O)[O-]",
        "[O-]P([O-])CC(=O)[O-]",
        "[O-]P(O)CC(=O)[O-]",
        "[O-]P(C)CC(=O)[O-]",
        "CP(=S)[O-]",
        "CP(C)(=S)[O-]",
    ] {
        assert!(smiles_to_name(smiles).is_err(), "{smiles}");
    }
}
