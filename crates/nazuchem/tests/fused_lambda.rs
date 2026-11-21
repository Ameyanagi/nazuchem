// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | crates/nazuchem/tests/fused_lambda.rs | issue #6

#![cfg(all(feature = "n2s", feature = "s2n"))]

use nazuchem::{
    Budget, NameReason, NameStatus, name_to_smiles, s2n::canonical_smiles, smiles_to_name,
};

#[test]
fn fused_bonding_number_breaks_equivalent_sulfur_locants() {
    // P-25.6 gives the first example. Reverse the input traversal and
    // move the ordinary sulfur to check completed-system locant selection.
    for (smiles, expected) in [
        ("[SH]1=CC=CSC2=C1C=CC=C2", "1lambda^4,5-benzodithiepine"),
        ("S1C=CC=[SH]C2=C1C=CC=C2", "1lambda^4,5-benzodithiepine"),
        ("[SH]1=CSC=CC2=C1C=CC=C2", "1lambda^4,3-benzodithiepine"),
    ] {
        let budget = Budget::default();
        let parsed = name_to_smiles(expected).unwrap();
        assert_eq!(
            canonical_smiles(&parsed, &budget).unwrap(),
            canonical_smiles(smiles, &budget).unwrap()
        );
        let name = smiles_to_name(smiles).unwrap();
        assert_eq!(name.text, expected);
        assert_eq!(name.status, NameStatus::Preferred);
    }
}

#[test]
fn fused_lambda_fallbacks_preserve_localized_and_exocyclic_double_bonds() {
    for smiles in [
        "C1[SH]=[SH]C(C)=NN=C2C=1C=CC=C2",
        "[SH]1=CC=CS(=O)C2=C1C=CC=C2",
        "[SH]1=CS(=O)C=CC2=C1C=CC=C2",
        "[SH]1=CC=CS(=O)(=O)C2=C1C=CC=C2",
        "[SH]1=CC=C[Se](=O)C2=C1C=CC=C2",
    ] {
        let name = smiles_to_name(smiles).unwrap();
        assert_eq!(
            name.status,
            NameStatus::NotPreferred(NameReason::VonBaeyerFusionFallback)
        );
        let parsed = name_to_smiles(&name.text).unwrap();
        let budget = Budget::default();
        assert_eq!(
            canonical_smiles(&parsed, &budget).unwrap(),
            canonical_smiles(smiles, &budget).unwrap(),
            "{}",
            name.text
        );
    }
}

#[test]
fn mixed_lambda_phosphine_oxide_keeps_its_stereochemical_fallback() {
    let base = "2lambda^4-thia-6lambda^5-phosphabicyclo[5.4.0]undeca-1(11),2,4,7,9-pentaen-6-one";
    for (smiles, stereo) in [
        ("[SH]1=CC=C[P@H](=O)C2=C1C=CC=C2", "(6S)"),
        ("[SH]1=CC=C[P@@H](=O)C2=C1C=CC=C2", "(6R)"),
    ] {
        let name = smiles_to_name(smiles).unwrap();
        assert_eq!(name.text, format!("{stereo}-{base}"));
        assert_eq!(
            name.status,
            NameStatus::NotPreferred(NameReason::VonBaeyerFusionFallback)
        );
        let parsed = name_to_smiles(&name.text).unwrap();
        let budget = Budget::default();
        assert_eq!(
            canonical_smiles(&parsed, &budget).unwrap(),
            canonical_smiles(smiles, &budget).unwrap(),
            "{}",
            name.text
        );
    }
}
