// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | crates/nazuchem/tests/chalcogen_anions.rs | issue #138

#![cfg(all(feature = "n2s", feature = "s2n"))]

//! Issue #138: peroxolate anions and anionic chalcogen prefixes are named
//! after the IUPAC 2013 recommendations and read back as the structure they
//! were generated from.

use nazuchem::{Budget, name_to_smiles, s2n::canonical_smiles, smiles_to_name};

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
    let forms = ["CC(=O)[S-]", "CC(=S)[O-]"].map(|s| canonical_smiles(s, &budget));
    for smiles in ["CC(=O)[S-]", "CC(=S)[O-]"] {
        assert_eq!(smiles_to_name(smiles).unwrap().text, "ethanethioate");
    }
    let parsed = name_to_smiles("ethanethioate").unwrap();
    assert!(forms.contains(&canonical_smiles(&parsed, &budget)));
}

#[test]
fn anionic_chalcogen_prefixes_that_read_back_wrong_are_declined() {
    // `oxidosulfanyl` and `selenidosulfanyl` (P-72.6.2) are read as sulfinyl
    // and selenosulfinyl groups, so no name is better than a wrong one.
    for smiles in ["[O-]SCC(=O)[O-]", "[O-][Se]CC(=O)[O-]", "[Se-]SCC(=O)[O-]"] {
        assert!(smiles_to_name(smiles).is_err(), "{smiles}");
    }
}
