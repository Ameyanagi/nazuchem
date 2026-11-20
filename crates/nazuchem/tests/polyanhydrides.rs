// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | crates/nazuchem/tests/polyanhydrides.rs | issue #110
// Source: Nazuchem | crates/nazuchem/tests/polyanhydrides.rs | issue #182

#![cfg(all(feature = "n2s", feature = "s2n"))]

//! Issue #110: anhydrides with several linkages (IUPAC 2013 P-65.7.6),
//! including the connection of unlike ends to specific acid groups.

use nazuchem::{
    Budget, NameReason, NameStatus, Stereo, StereoGroupKind, n2s, name_to_molecule, name_to_smiles,
    s2n::canonical_smiles, smiles_to_name,
};

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
        (
            "CC(=O)OC(=O)C(C)CC(=O)OC(C)=O",
            "diacetic 2-methylbutanedioic dianhydride",
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
            "CCCCCCCCCCCC(=O)OC(=O)CCC(=O)OC(=O)CCCCCCCCCCC",
            "di(dodecanoic) butanedioic dianhydride",
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
        // P-65.7.6.4.3 example: a dibasic thioacid.
        (
            "CC(=O)OC(=S)CCC(=S)OC(C)=O",
            "diacetic butanebis(thioic) dianhydride",
            true,
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
        // P-65.7.6.1 example.
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
    // P-65.7.6.2: unlike ends on a substituted dibasic acid need locants.
    assert!(name_to_smiles("acetic 2-methylbutanedioic propanoic dianhydride").is_err());
}

#[test]
fn unlike_ends_on_a_stereogenic_dibasic_acid_keep_distinct_names() {
    // P-65.7.6.2, P-91.2.1.1: without end locants the dianhydride name
    // would be the same for both structures, so it is not used.
    for pair in [
        [
            "CC(=O)OC(=O)[C@H](O)[C@H](O)C(=O)OC(=O)CC",
            "CC(=O)OC(=O)[C@@H](O)[C@@H](O)C(=O)OC(=O)CC",
        ],
        [
            "CC(=O)OC(=O)/C=C/CC/C=C\\C(=O)OC(=O)CC",
            "CCC(=O)OC(=O)/C=C/CC/C=C\\C(=O)OC(C)=O",
        ],
    ] {
        let names = pair.map(|smiles| {
            let name = smiles_to_name(smiles).unwrap();
            assert_ne!(name.status, NameStatus::Preferred, "{smiles}");
            same_structure(&name.text, smiles);
            name.text
        });
        assert_ne!(names[0], names[1]);
    }
}

#[test]
fn linear_polyanhydrides_use_substitutive_end_acids() {
    let examples = [
        (
            "CC(=O)OC(=O)CCC(=O)OC(=O)CCC(=O)OC(=O)CC",
            "acetic butanedioic 4-oxo-4-(propanoyloxy)butanoic dianhydride",
        ),
        (
            "CC(=O)OC(=O)CC(=O)OC(=O)CCC(=O)OC(C)=O",
            "acetic 3-(acetyloxy)-3-oxopropanoic butanedioic dianhydride",
        ),
        (
            "CC(=O)OC(=O)CCC(=O)OC(=O)CC(C)C(=O)OC(=O)CC",
            "4-[4-(acetyloxy)-4-oxobutanoic] 2-methylbutanedioic 1-propanoic dianhydride",
        ),
    ];
    for (smiles, expected) in examples {
        let name = smiles_to_name(smiles).unwrap();
        assert_eq!(name.text, expected, "{smiles}");
        same_structure(&name.text, smiles);
    }

    // The alternative linear description (P-65.7.6.2 method (2)) and the
    // example spelling of the first substitutive end denote the same graph.
    same_structure("acetic dibutanedioic propanoic trianhydride", examples[0].0);
    same_structure(
        "acetic butanedioic 4-(propanoyloxy)-4-oxobutanoic dianhydride",
        examples[0].0,
    );
    same_structure(
        "acetic butanedioic propanedioic acetic trianhydride",
        examples[1].0,
    );
}

#[test]
fn unlike_ends_are_located_on_the_dibasic_acid() {
    for (smiles, expected) in [
        (
            "CC(=O)OC(=O)C(C)CC(=O)OC(=O)CC",
            "1-acetic 2-methylbutanedioic 4-propanoic dianhydride",
        ),
        (
            "CC(=O)OC(=O)CC(C)C(=O)OC(=O)CC",
            "4-acetic 2-methylbutanedioic 1-propanoic dianhydride",
        ),
    ] {
        let name = smiles_to_name(smiles).unwrap();
        assert_eq!(name.text, expected, "{smiles}");
        same_structure(expected, smiles);
    }
    for invalid in [
        "1-acetic 2-methylbutanedioic 1-propanoic dianhydride",
        "1-acetic 2-methylbutanedioic 2-propanoic dianhydride",
        "1-acetic 2-methylbutanedioic propanoic dianhydride",
    ] {
        assert!(name_to_smiles(invalid).is_err(), "{invalid}");
    }
}

#[test]
fn polybasic_acid_names_use_locants_and_multipliers() {
    for (name, smiles) in [
        (
            "3,6-diacetic 2-propanoic naphthalene-2,3,6-tricarboxylic trianhydride",
            "CCC(=O)OC(=O)c1cc2ccc(C(=O)OC(C)=O)cc2cc1C(=O)OC(C)=O",
        ),
        (
            "triacetic benzene-1,3,5-tricarboxylic trianhydride",
            "CC(=O)OC(=O)c1cc(C(=O)OC(C)=O)cc(C(=O)OC(C)=O)c1",
        ),
        (
            "3,6-diacetic 2-propanoic naphthalene-2,3,6-tricarboxylic tris(thioanhydride)",
            "CCC(=O)SC(=O)c1cc2ccc(C(=O)SC(C)=O)cc2cc1C(=O)SC(C)=O",
        ),
    ] {
        let actual = smiles_to_name(smiles).unwrap();
        assert_eq!(actual.text, name, "{smiles}");
        same_structure(name, smiles);
    }
    assert!(
        name_to_smiles("2,3-diacetic 2-propanoic naphthalene-2,3,6-tricarboxylic trianhydride")
            .is_err()
    );
}

#[test]
fn locanted_polyanhydrides_with_unsaturated_ends_round_trip() {
    for smiles in [
        "CC(=O)OC(=O)C(C)CC(=O)OC(=O)C=CC",
        "CC(=O)OC(=O)C(C)CC(=O)OC(=O)/C=C/C",
        "CC(=O)OC(=O)c1cc(C(=O)OC(C)=O)cc(C(=O)OC(=O)C=CC)c1",
        "CC(=O)OC(=O)c1cc(C(=O)OC(C)=O)cc(C(=O)OC(=O)/C=C/C)c1",
    ] {
        let name = smiles_to_name(smiles).unwrap();
        same_structure(&name.text, smiles);
    }
}

#[test]
fn cyclic_dianhydride_keeps_the_acid_prefix_locants() {
    let smiles = "O=C1OC(=O)c2cc3c(cc12)C(=O)OC3=O";
    for name in [
        "1,2,4,5-benzenetetracarboxylic acid 1,2:4,5-dianhydride",
        "1,2,4,5-benzenetetracarboxylic 1,2:4,5-dianhydride",
    ] {
        same_structure(name, smiles);
    }
}

#[test]
fn symmetric_polybasic_acid_gives_lowest_locants_to_the_first_end_acid() {
    let expected = "1,3-diacetic 5-propanoic benzene-1,3,5-tricarboxylic trianhydride";
    for smiles in [
        "CC(=O)OC(=O)c1cc(C(=O)OC(C)=O)cc(C(=O)OC(=O)CC)c1",
        "CCC(=O)OC(=O)c1cc(C(=O)OC(C)=O)cc(C(=O)OC(C)=O)c1",
        "CC(=O)OC(=O)c1cc(C(=O)OC(=O)CC)cc(C(=O)OC(C)=O)c1",
    ] {
        let name = smiles_to_name(smiles).unwrap();
        assert_eq!(name.text, expected, "{smiles}");
        same_structure(&name.text, smiles);
    }
}

#[test]
fn enantiomeric_polyanhydride_ends_give_r_the_lower_locant() {
    let expected = "1-acetic 3-[(2R)-2-methylbutanoic] 5-[(2S)-2-methylbutanoic] benzene-1,3,5-tricarboxylic trianhydride";
    for smiles in [
        "CC(=O)OC(=O)c1cc(C(=O)OC(=O)[C@H](C)CC)cc(C(=O)OC(=O)[C@@H](C)CC)c1",
        "CC(=O)OC(=O)c1cc(C(=O)OC(=O)[C@@H](C)CC)cc(C(=O)OC(=O)[C@H](C)CC)c1",
    ] {
        let name = smiles_to_name(smiles).unwrap();
        assert_eq!(name.text, expected, "{smiles}");
        same_structure(&name.text, smiles);
    }
}

#[test]
fn identical_chiral_ends_are_multiplied_before_citing_the_enantiomer() {
    let expected = "1,3-bis[(2R)-2-methylbutanoic] 5-[(2S)-2-methylbutanoic] benzene-1,3,5-tricarboxylic trianhydride";
    for smiles in [
        "CC[C@@H](C)C(=O)OC(=O)c1cc(C(=O)OC(=O)[C@@H](C)CC)cc(C(=O)OC(=O)[C@H](C)CC)c1",
        "CC[C@@H](C)C(=O)OC(=O)c1cc(C(=O)OC(=O)[C@H](C)CC)cc(C(=O)OC(=O)[C@@H](C)CC)c1",
    ] {
        let name = smiles_to_name(smiles).unwrap();
        assert_eq!(name.text, expected, "{smiles}");
        same_structure(&name.text, smiles);
    }
}

#[test]
fn implicit_anhydride_ends_cannot_discard_locants() {
    for invalid in [
        "8-acetic benzene-1,2,4-tricarboxylic dianhydride",
        "9-acetic butanedioic dianhydride",
        "acetic 9-butanedioic dianhydride",
    ] {
        assert!(name_to_smiles(invalid).is_err(), "{invalid}");
    }

    let succinic = "CC(=O)OC(=O)CCC(=O)OC(C)=O";
    same_structure("acetic butanedioic dianhydride", succinic);
    same_structure("acetic 1,4-butanedioic dianhydride", succinic);
    same_structure(
        "1,4-diacetic 2-methylbutanedioic dianhydride",
        "CC(=O)OC(=O)C(C)CC(=O)OC(C)=O",
    );
}

#[test]
fn locanted_compound_acid_ends_are_enclosed() {
    for (smiles, expected) in [
        (
            "O=C(OC(=O)C)c1cc(C(=O)OC(=O)C)cc(C(=O)OC(=O)C(F)(F)F)c1",
            "1,3-diacetic 5-(trifluoroacetic) benzene-1,3,5-tricarboxylic trianhydride",
        ),
        (
            "CC(=O)OC(=O)c1cc(C(=O)OC(C)=O)cc(C(=O)OC(=O)CCl)c1",
            "1,3-diacetic 5-(chloroacetic) benzene-1,3,5-tricarboxylic trianhydride",
        ),
        (
            "CC(=O)OC(=O)c1cc(C(=O)OC(=O)CCl)cc(C(=O)OC(=O)C(F)(F)F)c1",
            "1-acetic 3-(chloroacetic) 5-(trifluoroacetic) benzene-1,3,5-tricarboxylic trianhydride",
        ),
        (
            "CC(=O)OC(=O)c1cc(C(=O)OC(=O)C(F)(F)F)cc(C(=O)OC(=O)C(F)(F)F)c1",
            "1-acetic 3,5-bis(trifluoroacetic) benzene-1,3,5-tricarboxylic trianhydride",
        ),
    ] {
        let name = smiles_to_name(smiles).unwrap();
        assert_eq!(name.text, expected, "{smiles}");
        same_structure(&name.text, smiles);
    }
}

#[test]
fn locanted_carbamic_ends_use_substitutive_names() {
    for (smiles, expected) in [
        (
            "O=C(OC(=O)C)C(C)CC(=O)OC(=O)N",
            "acetic 4-(carbamoyloxy)-2-methyl-4-oxobutanoic anhydride",
        ),
        (
            "O=C(OC(=O)C)C(C)CC(=O)OC(=O)NC",
            "acetic 2-methyl-4-[(methylcarbamoyl)oxy]-4-oxobutanoic anhydride",
        ),
        (
            "CC(=O)OC(=O)c1cc(C(=O)OC(C)=O)cc(C(=O)OC(=O)N)c1",
            "acetic 3-[(acetyloxy)carbonyl]-5-[(carbamoyloxy)carbonyl]benzoic anhydride",
        ),
    ] {
        let name = smiles_to_name(smiles).unwrap();
        assert_eq!(name.text, expected, "{smiles}");
        assert_ne!(name.status, NameStatus::Preferred, "{smiles}");
        same_structure(&name.text, smiles);
    }
}

#[test]
fn anhydrides_inherit_the_preference_of_each_cited_acid() {
    // P-65.1.8.1 recommends functional replacement for this structure.
    let acid = smiles_to_name("[N-]=[N+]=NC(=O)O").unwrap();
    assert_eq!(acid.text, "azidoformic acid");
    assert!(matches!(acid.status, NameStatus::NotPreferred(_)));
    for (smiles, expected) in [
        ("[N-]=[N+]=NC(=O)OC(=O)N=[N+]=[N-]", "azidoformic anhydride"),
        ("[N-]=[N+]=NC(=O)OC(C)=O", "acetic azidoformic anhydride"),
        ("CC(=O)OC(=O)N=[N+]=[N-]", "acetic azidoformic anhydride"),
        (
            "[N-]=[N+]=NC(=O)OC(=O)CCC(=O)OC(C)=O",
            "acetic butanedioic azidoformic dianhydride",
        ),
    ] {
        let name = smiles_to_name(smiles).unwrap();
        assert_eq!(name.text, expected, "{smiles}");
        assert_eq!(name.status, acid.status, "{smiles}");
        same_structure(&name.text, smiles);
    }

    // P-65.1.8.2 explicitly allows this substitution on formic acid.
    for (smiles, expected) in [
        ("[O-][N+](=O)C(=O)O", "nitroformic acid"),
        ("[O-][N+](=O)C(=O)OC(C)=O", "acetic nitroformic anhydride"),
    ] {
        let name = smiles_to_name(smiles).unwrap();
        assert_eq!(name.text, expected, "{smiles}");
        assert_eq!(name.status, NameStatus::Preferred, "{smiles}");
        same_structure(&name.text, smiles);
    }
}

#[test]
fn anhydride_preference_reasons_are_independent_of_atom_order() {
    for (pair, expected) in [
        (
            [
                "[N-]=[N+]=NC(=O)OC(=O)CCN(C(C)=O)C(C)=O",
                "CC(=O)N(C(C)=O)CCC(=O)OC(=O)N=[N+]=[N-]",
            ],
            "azidoformic 3-(diacetylamino)propanoic anhydride",
        ),
        (
            [
                "[N-]=[N+]=NC(=O)OC(=O)CCC(=O)OC(=O)CCN(C(C)=O)C(C)=O",
                "CC(=O)N(C(C)=O)CCC(=O)OC(=O)CCC(=O)OC(=O)N=[N+]=[N-]",
            ],
            "azidoformic butanedioic 3-(diacetylamino)propanoic dianhydride",
        ),
        (
            [
                "[N-]=[N+]=NC(=O)OC(=O)CCC(=O)OC(=O)CCC(=O)OC(=O)CCN(C(C)=O)C(C)=O",
                "CC(=O)N(C(C)=O)CCC(=O)OC(=O)CCC(=O)OC(=O)CCC(=O)OC(=O)N=[N+]=[N-]",
            ],
            "azidoformic butanedioic 4-{[3-(diacetylamino)propanoyl]oxy}-4-oxobutanoic dianhydride",
        ),
        (
            [
                "CC(=O)OC(=O)c1cc(C(=O)OC(C)=O)cc(C(=O)OC(=O)CCN(C(C)=O)C(C)=O)c1",
                "CC(=O)N(C(C)=O)CCC(=O)OC(=O)c1cc(C(=O)OC(C)=O)cc(C(=O)OC(C)=O)c1",
            ],
            "1,3-diacetic 5-[3-(diacetylamino)propanoic] benzene-1,3,5-tricarboxylic trianhydride",
        ),
    ] {
        for smiles in pair {
            let name = smiles_to_name(smiles).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(
                name.status,
                NameStatus::NotPreferred(NameReason::AminoPrefix),
                "{smiles}"
            );
            same_structure(&name.text, smiles);
        }
    }
}

#[test]
fn unlocanted_dianhydrides_can_have_free_acid_groups_on_the_ends() {
    for (name, smiles) in [
        (
            "acetic butanedioic butanedioic dianhydride",
            "CC(=O)OC(=O)CCC(=O)OC(=O)CCC(=O)O",
        ),
        (
            "acetic propanedioic butanedioic dianhydride",
            "CC(=O)OC(=O)CC(=O)OC(=O)CCC(=O)O",
        ),
        (
            "butanedioic propanedioic butanedioic dianhydride",
            "OC(=O)CCC(=O)OC(=O)CC(=O)OC(=O)CCC(=O)O",
        ),
        (
            "dibutanedioic propanedioic dianhydride",
            "OC(=O)CCC(=O)OC(=O)CC(=O)OC(=O)CCC(=O)O",
        ),
    ] {
        same_structure(name, smiles);
    }
    for invalid in [
        "1-acetic butanedioic 4-butanedioic dianhydride",
        "acetic 2-methylbutanedioic butanedioic dianhydride",
    ] {
        assert!(name_to_smiles(invalid).is_err(), "{invalid}");
    }
}

#[test]
fn unlocanted_dianhydrides_distinguish_enantiotopic_and_homotopic_sites() {
    let meso = "(2R,3S)-2,3-dihydroxybutanedioic";
    for name in [
        format!("acetic {meso} butanedioic dianhydride"),
        format!("acetic butanedioic {meso} dianhydride"),
        format!("acetic {meso} propanoic dianhydride"),
        format!("acetic {meso} butanedioic propanoic trianhydride"),
    ] {
        let error = name_to_smiles(&name).unwrap_err();
        assert!(error.to_string().contains("locants"), "{name}: {error}");
    }

    let locanted = [
        format!("1-acetic {meso} 4-propanoic dianhydride"),
        format!("4-acetic {meso} 1-propanoic dianhydride"),
    ];
    let expected = [
        "CC(=O)OC(=O)[C@H](O)[C@H](O)C(=O)OC(=O)CC",
        "CC(=O)OC(=O)[C@@H](O)[C@@H](O)C(=O)OC(=O)CC",
    ];
    for (name, smiles) in locanted.iter().zip(expected) {
        same_structure(name, smiles);
    }

    let chiral = "(2R,3R)-2,3-dihydroxybutanedioic";
    for (name, smiles) in [
        (
            format!("acetic {chiral} butanedioic dianhydride"),
            "CC(=O)OC(=O)[C@H](O)[C@@H](O)C(=O)OC(=O)CCC(=O)O",
        ),
        (
            format!("acetic butanedioic {chiral} dianhydride"),
            "CC(=O)OC(=O)CCC(=O)OC(=O)[C@H](O)[C@@H](O)C(=O)O",
        ),
        (
            format!("acetic {chiral} propanoic dianhydride"),
            "CC(=O)OC(=O)[C@H](O)[C@@H](O)C(=O)OC(=O)CC",
        ),
        (
            format!("diacetic {meso} dianhydride"),
            "CC(=O)OC(=O)[C@H](O)[C@H](O)C(=O)OC(C)=O",
        ),
    ] {
        same_structure(&name, smiles);
    }
}

#[test]
fn linear_polyanhydride_linkage_sites_include_double_bond_geometry() {
    let error =
        name_to_smiles("acetic (2E,6Z)-octa-2,6-dienedioic butanedioic dianhydride").unwrap_err();
    assert!(error.to_string().contains("locants"), "{error}");
    same_structure(
        "acetic (2E,6E)-octa-2,6-dienedioic butanedioic propanoic trianhydride",
        "CC(=O)OC(=O)/C=C/CC/C=C/C(=O)OC(=O)CCC(=O)OC(=O)CC",
    );
}

#[test]
fn equivalent_anhydride_sites_preserve_enhanced_stereo() {
    for (descriptor, kind) in [
        ("2R*,3R*", StereoGroupKind::Or(1)),
        ("2R*,3S*", StereoGroupKind::Or(1)),
        ("2RS,3RS", StereoGroupKind::And(1)),
        ("2RS,3SR", StereoGroupKind::And(1)),
    ] {
        let acid = format!("({descriptor})-2,3-dihydroxybutanedioic");
        for name in [
            format!("acetic {acid} butanedioic dianhydride"),
            format!("acetic butanedioic {acid} dianhydride"),
        ] {
            let molecule =
                name_to_molecule(&name).unwrap_or_else(|error| panic!("{name}: {error}"));
            assert_eq!(molecule.stereo_groups.len(), 1, "{name}");
            assert_eq!(molecule.stereo_groups[0].kind, kind, "{name}");
            assert_eq!(molecule.stereo_groups[0].stereo.len(), 2, "{name}");
        }
    }
}

#[test]
fn equivalent_anhydride_sites_preserve_an_axis() {
    let acid = "(1M)-[1,1'-binaphthalene]-2,2'-dicarboxylic";
    let parser = n2s::Parser::new().unwrap();
    for name in [
        format!("acetic {acid} butanedioic dianhydride"),
        format!("acetic butanedioic {acid} dianhydride"),
    ] {
        let result = parser.parse(&name, &n2s::ParseOptions::strict());
        assert_eq!(
            result.status,
            n2s::Status::Success,
            "{name}: {}",
            result.message
        );
        let structure = result.structure.unwrap();
        let (_, _, axes) = structure.semantic_cxsmiles_with_locants_and_axes().unwrap();
        assert_eq!(axes.len(), 1, "{name}");
        assert_eq!(axes[0].descriptor, "M", "{name}");
        let molecule = name_to_molecule(&name).unwrap();
        assert_eq!(molecule, structure.to_ir().unwrap(), "{name}");
        assert!(
            matches!(&molecule.stereo[..], [Stereo::Atropisomer { .. }]),
            "{name}"
        );
    }
}

#[test]
fn cumulene_geometry_distinguishes_dianhydride_linkage_sites() {
    let parser = n2s::Parser::new().unwrap();
    let central = "undeca-2,3,4,7,8,9-hexaenedioic";
    let ambiguous = format!("acetic (2E,7Z)-{central} propanoic dianhydride");
    let result = parser.parse(&ambiguous, &n2s::ParseOptions::strict());
    assert_eq!(result.status, n2s::Status::Failure);
    assert!(result.message.contains("locants"), "{}", result.message);

    for (name, expected) in [
        (
            format!("1-acetic (2E,7Z)-{central} 11-propanoic dianhydride"),
            r"CC(=O)OC(=O)/C=C=C=C/C/C=C=C=C\C(=O)OC(=O)CC",
        ),
        (
            format!("11-acetic (2E,7Z)-{central} 1-propanoic dianhydride"),
            r"CC(=O)OC(=O)/C=C=C=C\C/C=C=C=C/C(=O)OC(=O)CC",
        ),
        (
            format!("acetic (2E,7E)-{central} propanoic dianhydride"),
            "CC(=O)OC(=O)/C=C=C=C/C/C=C=C=C/C(=O)OC(=O)CC",
        ),
    ] {
        let result = parser.parse(&name, &n2s::ParseOptions::strict());
        assert_eq!(
            result.status,
            n2s::Status::Success,
            "{name}: {}",
            result.message
        );
        same_structure(&name, expected);
    }
}
