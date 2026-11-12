// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | crates/nazuchem/tests/cumulene_ez.rs | issue #89

#![cfg(all(feature = "n2s", feature = "s2n"))]

//! Issue #89: E/Z of a cumulene with an odd number of double bonds
//! (P-93.4.2.3) is named, parsed and kept in the structure identity.

use nazuchem::{
    Budget, ErrorKind, name_to_molecule, name_to_smiles, s2n::canonical_smiles, smiles_to_name,
};

fn identity(smiles: &str) -> String {
    canonical_smiles(smiles, &Budget::default()).unwrap_or_else(|error| panic!("{smiles}: {error}"))
}

#[test]
fn odd_cumulene_ez_names_follow_the_blue_book_and_round_trip() {
    for (smiles, expected) in [
        // P-93.4.2.3: (2Z)-hexa-2,3,4-triene is the Blue Book example; the
        // descriptor takes the locant at which the cumulated system starts.
        (r"C/C=C=C=C\C", "(2Z)-hexa-2,3,4-triene"),
        ("C/C=C=C=C/C", "(2E)-hexa-2,3,4-triene"),
        // The same structures written from the hydrogen atoms.
        ("[H]/C(C)=C=C=C(C)/[H]", "(2E)-hexa-2,3,4-triene"),
        (r"[H]/C(C)=C=C=C/C", "(2Z)-hexa-2,3,4-triene"),
        // P-93.4.2.3, (2E,6E)-octa-2,3,4,6-tetraene: descriptors of a
        // cumulene and of a double bond are cited together, by locant.
        ("C/C=C/C=C=C=C/C", "(2E,6E)-octa-2,3,4,6-tetraene"),
        (r"C/C=C/C=C=C=C\C", "(2Z,6E)-octa-2,3,4,6-tetraene"),
        ("C/C=C=C=C=C=C/C", "(2E)-octa-2,3,4,5,6-pentaene"),
        // The CIP-senior ligand of each terminal atom sets E or Z: fluorine
        // and ethyl are cis when the two methyl groups are.
        ("F/C=C=C=C/Cl", "(1E)-1-chloro-4-fluorobuta-1,2,3-triene"),
        (
            r"C/C(F)=C=C=C(\C)CC",
            "(2Z)-2-fluoro-5-methylhepta-2,3,4-triene",
        ),
        (
            "C/C(F)=C=C=C(/C)CC",
            "(2E)-2-fluoro-5-methylhepta-2,3,4-triene",
        ),
        (r"OC(=O)/C=C=C=C\C", "(2Z)-hexa-2,3,4-trienoic acid"),
        // Issue #89: the two examples differ only at the fluorine end.
        (
            r"C/C=C/C(=C=C=C/F)/C=C\C",
            "(1Z,5E)-1-fluoro-4-[(1Z)-prop-1-en-1-yl]hepta-1,2,3,5-tetraene",
        ),
        (
            r"C/C=C/C(=C=C=C\F)/C=C\C",
            "(1E,5E)-1-fluoro-4-[(1Z)-prop-1-en-1-yl]hepta-1,2,3,5-tetraene",
        ),
    ] {
        let name = smiles_to_name(smiles).unwrap_or_else(|error| panic!("{smiles}: {error:?}"));
        assert_eq!(name.text, expected, "{smiles}");
        let parsed =
            name_to_smiles(expected).unwrap_or_else(|error| panic!("{expected}: {error:?}"));
        assert_eq!(
            identity(&parsed),
            identity(smiles),
            "{expected} -> {parsed}"
        );
        assert_eq!(smiles_to_name(&parsed).unwrap().text, expected, "{parsed}");
    }
}

#[test]
fn odd_cumulene_configuration_is_part_of_the_identity() {
    let [e, z, unspecified] =
        ["C/C=C=C=C/C", r"C/C=C=C=C\C", "CC=C=C=CC"].map(|smiles| identity(smiles));
    assert_ne!(e, z);
    assert_ne!(e, unspecified);
    assert_ne!(z, unspecified);
    assert_eq!(identity(r"C\C=C=C=C\C"), e);
    assert_eq!(
        smiles_to_name("CC=C=C=CC").unwrap().text,
        "hexa-2,3,4-triene"
    );
    // P-93.4.2.3: cis and trans are the general-nomenclature descriptors.
    for (name, expected) in [
        ("cis-hexa-2,3,4-triene", &z),
        ("trans-hexa-2,3,4-triene", &e),
        ("(Z)-hexa-2,3,4-triene", &z),
        ("(5E)-hexa-2,3,4-triene", &e),
    ] {
        assert_eq!(
            &identity(&name_to_smiles(name).unwrap()),
            expected,
            "{name}"
        );
    }
}

#[test]
fn unsupported_cumulene_configurations_are_refused() {
    // The locant of an inner atom does not cite the cumulene, and a cumulene
    // with an even number of double bonds has no E/Z (P-93.4.2.2).
    for name in [
        "(3E)-hexa-2,3,4-triene",
        "(2E)-penta-2,3-diene",
        "(1E)-buta-1,2,3-triene",
    ] {
        assert!(name_to_smiles(name).is_err(), "{name}");
    }
    // Slashes at both ends of a cumulene specify it (OpenSMILES 3.8.3), so
    // SMILES cannot leave this one unspecified between two E double bonds.
    assert_eq!(
        name_to_smiles("(2E,8E)-deca-2,4,5,6,8-pentaene")
            .unwrap_err()
            .kind(),
        ErrorKind::Serialize
    );
    assert_eq!(
        name_to_smiles("(2E)-octa-2,4,5,6-tetraene").unwrap(),
        r"C\C=C\C=C=C=CC"
    );
    // Two identical ligands on a terminal atom leave nothing to specify.
    assert!(name_to_smiles("(2E,8E)-4-[(1E)-prop-1-en-1-yl]deca-2,4,5,6,8-pentaene").is_ok());
    // An unknown configuration of a cumulene has no representation yet.
    assert!(name_to_smiles("(2Ξ)-hexa-2,3,4-triene").is_err());
    // The ligands of the central cumulene differ only in the E/Z of the two
    // cumulenyl branches, which the ranking cannot see: refused, not dropped.
    for smiles in [
        r"F/C=C=C=C(/C=C=C=C/C)\C=C=C=C/C",
        r"F/C=C=C=C(\C=C=C=C\C)/C=C=C=C\C",
    ] {
        assert!(smiles_to_name(smiles).is_err(), "{smiles}");
    }
    // IR v1 has no record for this configuration.
    let error = name_to_molecule("(2E)-hexa-2,3,4-triene").unwrap_err();
    assert_eq!(error.kind(), ErrorKind::Serialize);
    // A cumulene that leaves the parent chain or lies in a ring is not named.
    for smiles in [
        // Both ligands of one terminal atom on the same side.
        r"C/C(\F)=C=C=C/C",
        r"C/C=C=C=C1/CCCCC1C",
        r"C/N=C=C=C\C",
        r"C1CCCC/C=C=C=C\CCCC1",
    ] {
        assert!(smiles_to_name(smiles).is_err(), "{smiles}");
    }
}
