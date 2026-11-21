// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | crates/nazuchem/tests/chalcogen_anions.rs | issue #138

#![cfg(all(feature = "n2s", feature = "s2n"))]

//! Issue #138: peroxolate anions and anionic chalcogen prefixes are named
//! after the IUPAC 2013 recommendations and read back as the structure they
//! were generated from.

use nazuchem::{
    Budget, ErrorKind, NameReason, NameStatus, molecule_to_name, name_to_molecule, name_to_smiles,
    s2n::{canonical_smiles, verify_identity},
    smiles_to_name,
};

#[test]
fn chalcogen_anion_names_follow_the_blue_book_and_round_trip() {
    for (smiles, expected) in [
        // P-72.2.2.2.2: methaneperoxolate and ethane(SO-thioperoxolate) are
        // the Blue Book examples; `ate` is added to the Table 6.1 suffixes
        // (P-63.4.2.1), and a letter locant stays inside the brackets.
        ("CO[O-]", "methaneperoxolate"),
        ("CCS[O-]", "ethane(SO-thioperoxolate)"),
        ("CCO[S-]", "ethane(OS-thioperoxolate)"),
        ("CCS[S-]", "ethanedithioperoxolate"),
        ("CC[Se][Se-]", "ethanediselenoperoxolate"),
        ("CC(C)S[O-]", "propane-2-(SO-thioperoxolate)"),
        ("CCCS[S-]", "propane-1-dithioperoxolate"),
        ("[O-]Sc1ccccc1", "benzene(SO-thioperoxolate)"),
        // P-72.2.2.2.2: ethane-1,2-bis(peroxolate),
        // benzene-1,4-bis(dithioperoxolate), benzene-1,4-bis(OS-thioperoxolate).
        ("[O-]OCCO[O-]", "ethane-1,2-bis(peroxolate)"),
        ("[S-]Sc1ccc(S[S-])cc1", "benzene-1,4-bis(dithioperoxolate)"),
        ("[S-]Oc1ccc(O[S-])cc1", "benzene-1,4-bis(OS-thioperoxolate)"),
        // P-63.8.1: the cation, then the anion as a separate word.
        ("[Na+].[O-]SCC", "sodium ethane(SO-thioperoxolate)"),
        // P-41 (Table 4.1): an anion is senior to an acid, which becomes the
        // carboxy prefix as in 5-carboxypentanoate (P-72.2.2.2.1.2).
        ("[S-]SCC(=O)O", "carboxymethanedithioperoxolate"),
        ("[O-]SCC(=O)O", "carboxymethane(SO-thioperoxolate)"),
        ("[S-]SCCC(=O)O", "2-carboxyethane-1-dithioperoxolate"),
        // P-72.6.3: -S-S(-) is disulfanidyl, the parent anion name and `yl`;
        // dioxidanidyl is formed the same way from dioxidanide.
        ("[S-]SCC(=O)[O-]", "disulfanidylacetate"),
        ("[S-]Sc1ccc(cc1)C(=O)[O-]", "4-(disulfanidyl)benzoate"),
        ("[S-]SCC[S-]", "2-(disulfanidyl)ethane-1-thiolate"),
        ("[O-]OCC(=O)[O-]", "dioxidanidylacetate"),
        // P-73.6, P-74.1.3: in a zwitterion the cation is a prefix on the
        // anionic parent.
        ("[O-]OCC[NH3+]", "2-azaniumylethane-1-peroxolate"),
        ("[O-]SCC[NH3+]", "2-azaniumylethane-1-(SO-thioperoxolate)"),
        // P-72.6.3: azanidyl on sulfanyl is not a chalcogen anion prefix.
        ("[NH-]SCC(=O)[O-]", "(azanidylsulfanyl)acetate"),
        // P-65.1.4.1: ethanediperoxoic acid; P-72.2.2.2.1.1 for the anion.
        ("OOC(=O)C(=O)OO", "ethanediperoxoic acid"),
        ("OOC(=O)CCC(=O)OO", "butanediperoxoic acid"),
        ("COOC(=O)CCC(=O)OOC", "dimethyl butanediperoxoate"),
        ("[O-]OC(=O)CCC(=O)O[O-]", "butanediperoxoate"),
    ] {
        let name = smiles_to_name(smiles).unwrap().text;
        assert_eq!(name, expected, "{smiles}");
        let parsed = name_to_smiles(&name).unwrap_or_else(|e| panic!("{name}: {e:?}"));
        let budget = Budget::default();
        assert_eq!(
            canonical_smiles(&parsed, &budget),
            canonical_smiles(smiles, &budget),
            "{smiles} -> {name} -> {parsed}"
        );
    }
}

#[test]
fn thioate_anion_has_one_name_for_both_mesomeric_forms() {
    // P-72.2.2.2.1.1: CH3-CO-S(-) <-> CH3-CS-O(-) is ethanethioate; the
    // anion has no letter locant, so the name reads back as one of the two.
    let budget = Budget::default();
    for (s_form, o_form, expected) in [
        ("CC(=O)[S-]", "CC(=S)[O-]", "ethanethioate"),
        ("CCC(=O)[S-]", "CCC(=S)[O-]", "propanethioate"),
        ("O=C[S-]", "S=C[O-]", "methanethioate"),
        (
            "O=C([S-])c1ccccc1",
            "S=C([O-])c1ccccc1",
            "benzenecarbothioate",
        ),
        (
            "[S-]C(=O)CC(=S)[O-]",
            "[O-]C(=S)CC(=S)[O-]",
            "propanebis(thioate)",
        ),
        (
            "[Na+].CC(=O)[S-]",
            "[Na+].CC(=S)[O-]",
            "sodium ethanethioate",
        ),
    ] {
        let parsed = name_to_smiles(expected).unwrap();
        let exact_o = canonical_smiles(o_form, &budget).unwrap();
        assert_eq!(canonical_smiles(&parsed, &budget).unwrap(), exact_o);
        assert_ne!(canonical_smiles(s_form, &budget).unwrap(), exact_o);
        for smiles in [s_form, o_form] {
            let name = smiles_to_name(smiles).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, NameStatus::Preferred, "{smiles}");
            verify_identity(smiles, &parsed, &budget).unwrap();
            verify_identity(&parsed, smiles, &budget).unwrap();
        }
    }
}

#[test]
fn thiocarboxylate_round_trip_keeps_explicit_isotopes() {
    let budget = Budget::default();
    for (left, right) in [
        ("[1H]C(=O)[S-]", "[1H]C(=S)[O-]"),
        ("[2H]C(=O)[S-]", "[2H]C(=S)[O-]"),
        ("[3H]C(=O)[S-]", "[3H]C(=S)[O-]"),
        ("[2H][13C](=[18O])[34S-]", "[2H][13C](=[34S])[18O-]"),
        ("C[C@]([2H])(F)C(=O)[S-]", "C[C@]([2H])(F)C(=S)[O-]"),
    ] {
        assert_ne!(
            canonical_smiles(left, &budget).unwrap(),
            canonical_smiles(right, &budget).unwrap()
        );
        verify_identity(left, right, &budget).unwrap();
        verify_identity(right, left, &budget).unwrap();
    }
    let deuterated = "[2H]C(=O)[S-]";
    for other in ["C(=S)[O-]", "[1H]C(=S)[O-]", "[3H]C(=S)[O-]", "[2H]C(=S)O"] {
        assert!(
            verify_identity(deuterated, other, &budget).is_err(),
            "{other}"
        );
        assert!(
            verify_identity(other, deuterated, &budget).is_err(),
            "{other}"
        );
    }
}

#[test]
fn thiocarboxylate_round_trip_preserves_other_features() {
    let budget = Budget::default();
    for (left, right) in [
        ("CC(=[18O])[34S-]", "CC(=[34S])[18O-]"),
        ("C[C@H](F)C(=O)[S-]", "C[C@H](F)C(=S)[O-]"),
        ("C/C=C/C(=O)[S-]", "C/C=C/C(=S)[O-]"),
        ("[S-]C(=O)CC(=S)[O-]", "[O-]C(=S)CC(=O)[S-]"),
    ] {
        verify_identity(left, right, &budget).unwrap();
        verify_identity(right, left, &budget).unwrap();
    }
    for (left, right) in [
        ("CC(=O)S", "CC(=S)O"),
        ("CC(=O)SC", "CC(=S)OC"),
        ("CC(=O)[S-]", "CC(=S)O"),
        ("CC(=O)[S-]", "CC(=S)O[O-]"),
        ("CC(=[18O])[S-]", "CC(=S)[O-]"),
        ("C[C@H](F)C(=O)[S-]", "C[C@@H](F)C(=S)[O-]"),
        ("C[C@H](C(=O)[S-])C(=S)[O-]", "CC(C(=S)[O-])C(=S)[O-]"),
        ("C/C=C/C(=O)[S-]", "C/C=C\\C(=S)[O-]"),
        ("[O-]SCC(=O)[S-]", "[S-]OCC(=S)[O-]"),
    ] {
        assert!(
            verify_identity(left, right, &budget).is_err(),
            "{left} vs {right}"
        );
        assert!(
            verify_identity(right, left, &budget).is_err(),
            "{right} vs {left}"
        );
    }
}

#[test]
fn heteroatom_bound_thioates_preserve_the_location_of_the_charge() {
    let budget = Budget::default();
    for (sulfur_anion, sulfur_name, oxygen_anion, expected) in [
        (
            "NC(=O)[S-]",
            "amino(oxo)methanethiolate",
            "NC(=S)[O-]",
            "aminomethanethioate",
        ),
        (
            "N1(C(=O)[S-])CCCCC1",
            "oxo(piperidin-1-yl)methanethiolate",
            "N1(C(=S)[O-])CCCCC1",
            "piperidine-1-carbothioate",
        ),
    ] {
        let sulfur_molecule = name_to_molecule(sulfur_name).unwrap();
        let parsed = name_to_smiles(sulfur_name).unwrap();
        verify_identity(sulfur_anion, &parsed, &budget).unwrap();
        for result in [
            smiles_to_name(sulfur_anion),
            molecule_to_name(&sulfur_molecule),
        ] {
            let error = result.unwrap_err();
            assert_eq!(error.kind(), ErrorKind::Naming);
            assert!(
                error
                    .to_string()
                    .contains("heteroatom-bound thioate charge on sulfur")
            );
        }

        let oxygen_molecule = name_to_molecule(expected).unwrap();
        for result in [
            smiles_to_name(oxygen_anion),
            molecule_to_name(&oxygen_molecule),
        ] {
            let name = result.unwrap();
            assert_eq!(name.text, expected);
            assert_eq!(name.status, NameStatus::Unverified(NameReason::AnionParent));
            let parsed = name_to_smiles(&name.text).unwrap();
            assert_eq!(
                canonical_smiles(&parsed, &budget).unwrap(),
                canonical_smiles(oxygen_anion, &budget).unwrap()
            );
        }
        assert!(verify_identity(sulfur_anion, oxygen_anion, &budget).is_err());
        assert!(verify_identity(oxygen_anion, sulfur_anion, &budget).is_err());
    }
}

#[test]
fn anionic_chalcogen_prefixes_that_read_back_wrong_are_declined() {
    // `oxidosulfanyl` and `selenidosulfanyl` (P-72.6.2) are read as sulfinyl
    // and selenosulfinyl groups, so no name is better than a wrong one.
    for smiles in ["[O-]SCC(=O)[O-]", "[O-][Se]CC(=O)[O-]", "[Se-]SCC(=O)[O-]"] {
        assert!(smiles_to_name(smiles).is_err(), "{smiles}");
    }
}

#[test]
fn multiplied_peroxy_infix_beside_another_infix_stays_on_one_acid() {
    // P-65.2.1.2 (carbonodiperoxoic acid): the two singly bonded oxygens take
    // the peroxy groups and thio replaces the doubly bonded one.
    let budget = Budget::default();
    for (name, smiles) in [
        ("carbonodiperoxothioic acid", "OOC(=S)OO"),
        ("carbonodiperoxoic acid", "OOC(=O)OO"),
        ("phosphorotriperoxothioic acid", "OOP(=S)(OO)OO"),
    ] {
        let parsed = name_to_smiles(name).unwrap_or_else(|e| panic!("{name}: {e:?}"));
        assert_eq!(
            canonical_smiles(&parsed, &budget),
            canonical_smiles(smiles, &budget),
            "{name}"
        );
    }
}
