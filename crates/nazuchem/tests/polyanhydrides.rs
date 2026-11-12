// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | crates/nazuchem/tests/polyanhydrides.rs | issue #110

#![cfg(all(feature = "n2s", feature = "s2n"))]

//! Issue #110: a dibasic acid whose two acid groups each form an anhydride
//! with a monobasic acid is named as a dianhydride (IUPAC 2013 P-65.7.6),
//! and the name reads back as the structure it was generated from.

use nazuchem::{Budget, NameStatus, name_to_smiles, s2n::canonical_smiles, smiles_to_name};

fn same_structure(name: &str, smiles: &str) {
    let parsed = name_to_smiles(name).unwrap_or_else(|e| panic!("{name}: {e:?}"));
    let budget = Budget::default();
    assert_eq!(
        canonical_smiles(&parsed, &budget).unwrap(),
        canonical_smiles(smiles, &budget).unwrap(),
        "{name}"
    );
}

#[test]
fn dianhydrides_follow_the_blue_book_and_round_trip() {
    for (smiles, expected, preferred) in [
        // P-65.7.6.1 examples: the acids in their order of occurrence,
        // beginning with the end acid that is first in alphanumerical order.
        (
            "CC(=O)OC(=O)CCC(=O)OC(C)=O",
            "diacetic butanedioic dianhydride",
            true,
        ),
        (
            "CCC(=O)OC(=O)CCC(=O)OC(=O)CCC",
            "butanoic butanedioic propanoic dianhydride",
            true,
        ),
        // The four heteroatoms on the O-S-C-S-O path keep the conservative
        // unverified label of the skeletal replacement check (P-51.4.1).
        (
            "CC(=O)OS(=O)(=O)CS(=O)(=O)OC(=O)CC",
            "acetic methanedisulfonic propanoic dianhydride",
            false,
        ),
        // P-65.7.6.1 writes "diacetic 2-methylbutanedioic dianhydride"; the
        // acid namer omits the locant, so the name stays unverified.
        (
            "CC(=O)OC(=O)C(C)CC(=O)OC(C)=O",
            "diacetic methylbutanedioic dianhydride",
            false,
        ),
        // P-65.7.6.1 with other dibasic acids.
        (
            "CC(=O)OC(=O)CCC(=O)OC(=O)CC",
            "acetic butanedioic propanoic dianhydride",
            true,
        ),
        (
            "CC(=O)OC(=O)c1ccc(cc1)C(=O)OC(C)=O",
            "diacetic benzene-1,4-dicarboxylic dianhydride",
            true,
        ),
        (
            "CC(=O)OS(=O)(=O)CCS(=O)(=O)OC(C)=O",
            "diacetic ethane-1,2-disulfonic dianhydride",
            false,
        ),
        // P-16.3.5: a substituted acid is multiplied with bis (the acid of
        // P-65.7.1 "chloroacetic anhydride").
        (
            "ClCC(=O)OC(=O)CCC(=O)OC(=O)CCl",
            "bis(chloroacetic) butanedioic dianhydride",
            true,
        ),
        // P-16.3.4 (c), (d), (e): di(pentanal), di(decanoic acid) and
        // di(cyclohexane-1-carboxylic acid) enclose the multiplied acid.
        (
            "CCCCCCCCCC(=O)OC(=O)CCC(=O)OC(=O)CCCCCCCCC",
            "di(decanoic) butanedioic dianhydride",
            true,
        ),
        (
            "CCCCC(=O)OC(=O)CCC(=O)OC(=O)CCCC",
            "di(pentanoic) butanedioic dianhydride",
            true,
        ),
        (
            "C=CCC(=O)OC(=O)CCC(=O)OC(=O)CC=C",
            "di(but-3-enoic) butanedioic dianhydride",
            true,
        ),
        // P-65.7.6.4.1 examples: identical chalcogen linkages.
        (
            "CC(=O)SC(=O)CCC(=O)SC(C)=O",
            "diacetic butanedioic bis(thioanhydride)",
            false,
        ),
        (
            "CC(=O)SS(=O)(=O)CS(=O)(=O)SC(C)=O",
            "diacetic methanedisulfonic bis(thioanhydride)",
            false,
        ),
        // P-65.7.6.4.3 example: a thioacid end.
        (
            "CC(=S)OC(=O)CCC(=O)OC(C)=O",
            "acetic butanedioic ethanethioic dianhydride",
            true,
        ),
        // P-65.7.6.4.2 example: unlike linkages keep the oxygen anhydride
        // and cite the other substitutively.
        (
            "CC(=O)OC(=O)CCCC(=O)SC(C)=O",
            "acetic 5-(acetylsulfanyl)-5-oxopentanoic anhydride",
            false,
        ),
        // P-65.7.2 (benzoic ethanethioic anhydride): a thioacid end on an
        // oxygen bridge, here a sulfonothioic acid.
        (
            "CC(=O)OS(=O)(=S)C",
            "acetic methanesulfonothioic anhydride",
            true,
        ),
    ] {
        let name = smiles_to_name(smiles).unwrap();
        assert_eq!(name.text, expected, "{smiles}");
        assert_eq!(name.status == NameStatus::Preferred, preferred, "{smiles}");
        same_structure(&name.text, smiles);
    }
}

#[test]
fn blue_book_dianhydride_names_are_read() {
    for (name, smiles) in [
        // P-65.7.6.1 example, with the locant the generated name omits.
        (
            "diacetic 2-methylbutanedioic dianhydride",
            "CC(=O)OC(=O)C(C)CC(=O)OC(C)=O",
        ),
        // P-65.7.6.4.3 example.
        (
            "diacetic butanebis(thioic) dianhydride",
            "CC(=O)OC(=S)CCC(=S)OC(C)=O",
        ),
    ] {
        same_structure(name, smiles);
    }
    // The dibasic acid stands between the two acids it links (P-65.7.6.1).
    assert!(name_to_smiles("acetic propanoic butanedioic dianhydride").is_err());
    // P-65.7.6.3: unlike ends on unlike acid groups need locants.
    assert!(name_to_smiles("acetic 2-methylbutanedioic propanoic dianhydride").is_err());
}

#[test]
fn unlike_ends_on_a_chiral_dibasic_acid_keep_distinct_names() {
    // P-65.7.6.3, P-91.2.1.1: without end locants the dianhydride name
    // would be the same for both enantiomers, so it is not used.
    let names = [
        "CC(=O)OC(=O)[C@H](O)[C@H](O)C(=O)OC(=O)CC",
        "CC(=O)OC(=O)[C@@H](O)[C@@H](O)C(=O)OC(=O)CC",
    ]
    .map(|smiles| {
        let name = smiles_to_name(smiles).unwrap();
        assert_ne!(name.status, NameStatus::Preferred, "{smiles}");
        same_structure(&name.text, smiles);
        name.text
    });
    assert_ne!(names[0], names[1]);
}
