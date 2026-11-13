// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | crates/nazuchem/tests/multiplied_ester_anions.rs | issue #170

#![cfg(all(feature = "n2s", feature = "s2n"))]

//! Issue #170: a multiplied ester anion of a salt, as in
//! `calcium bis(methyl hydrogen phosphate)`, reads back as the structure it
//! was generated from. Names follow the IUPAC 2013 recommendations.

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
fn a_generated_multiplied_ester_anion_reads_back() {
    for (smiles, expected) in [
        // P-67.1.3.2: organyl groups are separate words before the anion
        // name, then `hydrogen` for a remaining acid group, and the cation
        // comes first (sodium methyl hydrogen phosphate). P-16.3.6 (a), (e):
        // bis(hydrogen sulfate) and tris(dihydrogen phosphate) multiply an
        // anion of several words.
        (
            "[Ca+2].COP(=O)(O)[O-].COP(=O)(O)[O-]",
            "calcium bis(methyl hydrogen phosphate)",
        ),
        (
            "[Ca+2].COS(=O)[O-].COS(=O)[O-]",
            "calcium bis(methyl sulfite)",
        ),
        (
            "[Ca+2].CCOP(=O)(O)[O-].CCOP(=O)(O)[O-]",
            "calcium bis(ethyl hydrogen phosphate)",
        ),
        (
            "[Ca+2].CC(C)OP(=O)(O)[O-].CC(C)OP(=O)(O)[O-]",
            "calcium bis(propan-2-yl hydrogen phosphate)",
        ),
        (
            "[Ca+2].COP(O)[O-].COP(O)[O-]",
            "calcium bis(methyl hydrogen phosphite)",
        ),
        (
            "[Ca+2].COP(OC)[O-].COP(OC)[O-]",
            "calcium bis(dimethyl phosphite)",
        ),
        (
            "[Ca+2].COB(O)[O-].COB(O)[O-]",
            "calcium bis(methyl hydrogen borate)",
        ),
        (
            "[Ca+2].CO[As](=O)(O)[O-].CO[As](=O)(O)[O-]",
            "calcium bis(methyl hydrogen arsorate)",
        ),
        (
            "[Ca+2].CO[Se](=O)[O-].CO[Se](=O)[O-]",
            "calcium bis(methyl selenite)",
        ),
    ] {
        let name = smiles_to_name(smiles).unwrap();
        assert_eq!(name.text, expected, "{smiles}");
        assert_eq!(name.status, NameStatus::Preferred, "{smiles}");
        assert_round_trip(smiles, &name.text);
    }
}

#[test]
fn a_multiplied_ester_anion_is_read_as_its_words() {
    for (name, smiles) in [
        // P-67.1.3.2: a word before the anion name is an ester group,
        // HP(O)(O-CH3)2 dimethyl phosphonate, and a group joined to it is a
        // substituent of the acid, CH3-P(O)(O-CH2-CH3)2 diethyl
        // methylphosphonate.
        (
            "calcium bis(methyl phosphonate)",
            "[Ca+2].CO[PH](=O)[O-].CO[PH](=O)[O-]",
        ),
        (
            "calcium bis(phenyl phosphonate)",
            "[Ca+2].[O-][PH](=O)Oc1ccccc1.[O-][PH](=O)Oc1ccccc1",
        ),
        (
            "calcium bis(ethyl methylphosphonate)",
            "[Ca+2].CCOP(C)(=O)[O-].CCOP(C)(=O)[O-]",
        ),
        (
            "calcium bis(methyl thiosulfate)",
            "[Ca+2].COS(=O)(=S)[O-].COS(=O)(=S)[O-]",
        ),
        (
            "barium bis(phenyl hydrogen phosphate)",
            "[Ba+2].OP(=O)([O-])Oc1ccccc1.OP(=O)([O-])Oc1ccccc1",
        ),
        (
            "aluminium tris(methyl sulfite)",
            "[Al+3].COS(=O)[O-].COS(=O)[O-].COS(=O)[O-]",
        ),
        (
            "dicalcium bis(methyl hydrogen diphosphate)",
            "[Ca+2].[Ca+2].COP(=O)(O)OP(=O)([O-])[O-].COP(=O)(O)OP(=O)([O-])[O-]",
        ),
        // P-65.6.3.2.1: esters of organic acids are named the same way.
        (
            "calcium bis(methyl malonate)",
            "[Ca+2].COC(=O)CC(=O)[O-].COC(=O)CC(=O)[O-]",
        ),
        // Without the space the group is a substituent of the acid where
        // the acid has a hydrogen for it (P-67.1.3.1, disodium
        // methylphosphonite), as it is without the multiplier.
        (
            "calcium bis(methylphosphonate)",
            "[Ca+2].[Ca+2].CP(=O)([O-])[O-].CP(=O)([O-])[O-]",
        ),
        (
            "calcium bis(methylmalonate)",
            "[Ca+2].[Ca+2].CC(C(=O)[O-])C(=O)[O-].CC(C(=O)[O-])C(=O)[O-]",
        ),
        ("calcium methylsulfite", "[Ca+2].COS(=O)[O-].COS(=O)[O-]"),
        (
            "calcium bis(methylsulfite)",
            "[Ca+2].COS(=O)[O-].COS(=O)[O-]",
        ),
        // An ester that would leave no anion is not taken out of the
        // brackets; these names are read as before.
        (
            "calcium bis(methyl carbamate)",
            "[Ca+2].CNC(=O)[O-].CNC(=O)[O-]",
        ),
        (
            "calcium bis(4-methyl benzoate)",
            "[Ca+2].Cc1ccc(cc1)C(=O)[O-].Cc1ccc(cc1)C(=O)[O-]",
        ),
        (
            "calcium bis(hydrogen methylphosphonate)",
            "[Ca+2].CP(=O)(O)[O-].CP(=O)(O)[O-]",
        ),
        (
            "calcium bis(hydrogen sulfate)",
            "[Ca+2].OS(=O)(=O)[O-].OS(=O)(=O)[O-]",
        ),
        (
            "calcium bis(ethyl methyl phosphate)",
            "[Ca+2].CCOP(=O)(OC)[O-].CCOP(=O)(OC)[O-]",
        ),
        ("calcium diacetate", "[Ca+2].CC(=O)[O-].CC(=O)[O-]"),
    ] {
        assert_round_trip(smiles, name);
    }
}

#[test]
fn the_salt_names_of_the_golden_fixture_read_back() {
    let golden = include_str!("../../../tests/fixtures/golden/s2n_n2s.jsonl");
    let mut salts = 0;
    for line in golden.lines() {
        let record: serde_json::Value = serde_json::from_str(line).unwrap();
        let smiles = record["smiles"].as_str().unwrap();
        let Some(name) = record["s2n"]["Ok"]["text"].as_str() else {
            continue;
        };
        if smiles.contains('.') {
            assert_round_trip(smiles, name);
            salts += 1;
        }
    }
    assert!(salts > 0);
}
