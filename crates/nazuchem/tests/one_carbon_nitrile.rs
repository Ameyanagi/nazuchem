// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | crates/nazuchem/tests/one_carbon_nitrile.rs | issue #158

#![cfg(all(feature = "n2s", feature = "s2n"))]

//! Issue #158: a nitrile carbon bonded to a sulfur acid group is not the
//! methane of a `cyanomethanesulfonic acid`, and the dianion of sulfoformic
//! acid keeps both charges. Names follow the IUPAC 2013 recommendations and
//! read back as the structure they came from.

use nazuchem::{Budget, NameStatus, name_to_smiles, s2n::canonical_smiles, smiles_to_name};

fn assert_round_trip(smiles: &str, name: &str) {
    let parsed = name_to_smiles(name).unwrap_or_else(|e| panic!("{name}: {e:?}"));
    let budget = Budget::default();
    assert_eq!(
        canonical_smiles(&parsed, &budget),
        canonical_smiles(smiles, &budget),
        "{smiles} -> {name} -> {parsed}"
    );
}

#[test]
fn cyanides_of_sulfuric_and_sulfurous_acid_use_the_cyanido_infix() {
    for (smiles, expected) in [
        // P-67.1.2.4.1.1: HO-SO2-CN sulfurocyanidic acid (PIN), sulfuric
        // acid with the infix cyanido (P-67.1.2.3.2).
        ("N#CS(O)(=O)=O", "sulfurocyanidic acid"),
        // P-67.1.3.1: `ic acid` becomes `ate`, `ous acid` becomes `ite`.
        ("N#CS(=O)(=O)[O-]", "sulfurocyanidate"),
        ("[Na+].N#CS(=O)(=O)[O-]", "sodium sulfurocyanidate"),
        // P-16.3.6 (a): `disulfuro...` would be read as a derivative of
        // disulfuric acid (P-67.2.1).
        (
            "[Ca+2].N#CS(=O)(=O)[O-].N#CS(=O)(=O)[O-]",
            "calcium bis(sulfurocyanidate)",
        ),
        // P-67.1.3.2: esters, as in methyl phenylphosphonochloridate.
        ("N#CS(=O)(=O)OC", "methyl sulfurocyanidate"),
        ("N#CS(=O)(=O)Oc1ccccc1", "phenyl sulfurocyanidate"),
        (
            "N#CS(=O)(=O)O[C@H](C)CC",
            "(2R)-butan-2-yl sulfurocyanidate",
        ),
        // P-67.1.2.1: sulfurous acid is modified by infixes in the same way,
        // as in N-methylsulfurimidous dichloride (P-67.1.2.5.1).
        ("N#CS(O)=O", "sulfurocyanidous acid"),
        ("N#CS(=O)[O-]", "sulfurocyanidite"),
        ("N#CS(=O)OCC", "ethyl sulfurocyanidite"),
    ] {
        let name = smiles_to_name(smiles).unwrap();
        assert_eq!(name.text, expected, "{smiles}");
        assert_eq!(name.status, NameStatus::Preferred, "{smiles}");
        assert_round_trip(smiles, &name.text);
    }
}

#[test]
fn a_nitrile_carbon_that_is_the_whole_parent_is_not_cited_as_cyano() {
    // P-66.5.1.1.4: cyano is the whole -CN group on a parent, so the nitrile
    // carbon cannot also be the methane of `cyanomethanesulfon...`. The
    // names these need (P-67.1.2.5 halides, P-67.1.2.6 amides, chalcogen
    // infixes, selenium and tellurium acids, anhydrides) are not generated,
    // and an ester whose alcohol part holds a principal group is not ranked
    // against it here.
    for smiles in [
        "N#CS(=O)(=O)Cl",
        "N#CS(N)(=O)=O",
        "N#CS(=O)(=O)NC",
        "N#CS(O)(=O)=S",
        "N#C[Se](O)(=O)=O",
        "N#CS(=O)(=O)OS(=O)(=O)C#N",
        "N#CS(=O)(=O)OC(C)=O",
        "N#CS(=O)(=O)OCCO",
        // A core in another charge or hydrogen state is not the plain
        // ester, acid or anion.
        "N#CS(=O)(=O)[OH2+]",
        "N#CS(=O)(=O)[OH+]C",
        "N#C[SH2](=O)O",
    ] {
        assert!(smiles_to_name(smiles).is_err(), "{smiles}");
    }
}

#[test]
fn neighbours_of_the_one_carbon_nitrile_keep_their_names() {
    for (smiles, expected) in [
        // P-66.5.1.1.4: cyano on a parent that has its own carbon.
        ("N#CCS(O)(=O)=O", "cyanomethanesulfonic acid"),
        ("OS(=O)(=O)C(C#N)S(O)(=O)=O", "cyanomethanedisulfonic acid"),
        // Unchanged and not labelled preferred: the PIN is carbonocyanidic
        // acid (P-65.2.1.4), which is not generated.
        ("N#CC(=O)O", "cyanoformic acid"),
        // P-67.1.4.4: cyanosulfonyl is the preferred prefix for -SO2-CN.
        ("N#CS(=O)(=O)OCC(=O)O", "(cyanosulfonyloxy)acetic acid"),
        ("N#CS(=O)(=O)CC(=O)O", "(cyanosulfonyl)acetic acid"),
        // P-41: a cation is senior to an ester, so a charged alcohol part
        // stays the parent.
        (
            "N#CS(=O)(=O)Oc1cc[n+](C)cc1",
            "4-(cyanosulfonyloxy)-1-methylpyridin-1-ium",
        ),
        ("[Ca+2].N#C[S-].N#C[S-]", "calcium dithiocyanate"),
    ] {
        let name = smiles_to_name(smiles).unwrap().text;
        assert_eq!(name, expected, "{smiles}");
        assert_round_trip(smiles, &name);
    }
}

#[test]
fn the_dianion_of_sulfoformic_acid_keeps_both_charges() {
    for (smiles, expected) in [
        // P-72.6.1, P-65.6.2.3.1: -SO2-O(-) under a more senior anion is
        // sulfonato; the free acid group stays sulfo.
        ("[O-]C(=O)S(=O)(=O)[O-]", "sulfonatoformate"),
        (
            "[Na+].[Na+].[O-]C(=O)S(=O)(=O)[O-]",
            "disodium sulfonatoformate",
        ),
        ("[O-]C(=O)S(=O)[O-]", "sulfinatoformate"),
        ("[O-]C(=O)S(=O)(=O)O", "sulfoformate"),
        ("OC(=O)S(=O)(=O)O", "sulfoformic acid"),
    ] {
        let name = smiles_to_name(smiles).unwrap();
        assert_eq!(name.text, expected, "{smiles}");
        // These names were not labelled preferred before this fix and
        // still are not; whether P-65.1.8.2 allows them was not checked.
        assert_ne!(name.status, NameStatus::Preferred, "{smiles}");
        assert_round_trip(smiles, &name.text);
    }
}
