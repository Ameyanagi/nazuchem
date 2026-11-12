// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | crates/nazuchem/tests/anionic_prefixes.rs | issue #145

#![cfg(all(feature = "n2s", feature = "s2n"))]

//! Issue #145: acid-anion prefixes under a more senior anion and anionic
//! oxygen on a chalcogen chain keep their charge, are named after the IUPAC
//! 2013 recommendations and read back as the structure they came from.

use nazuchem::{Budget, name_to_smiles, s2n::canonical_smiles, smiles_to_name};

#[test]
fn anionic_prefix_names_follow_the_blue_book_and_round_trip() {
    for (smiles, expected) in [
        // P-72.6.1, P-65.6.2.3.1: -SO2-O(-) is sulfonato; the free acid
        // stays sulfo, as in potassium 6-carboxyhexanoate.
        ("[O-]S(=O)(=O)CC(=O)[O-]", "sulfonatoacetate"),
        ("[O-]C(=O)CCS(=O)(=O)[O-]", "3-sulfonatopropanoate"),
        ("[O-]S(=O)(=O)c1ccc(cc1)C(=O)[O-]", "4-sulfonatobenzoate"),
        (
            "[O-]C(=O)CC(S(=O)(=O)[O-])CS(=O)(=O)[O-]",
            "3,4-disulfonatobutanoate",
        ),
        (
            "[O-]C(=O)CC(S(=O)(=O)[O-])CS(=O)(=O)O",
            "4-sulfo-3-sulfonatobutanoate",
        ),
        ("OS(=O)(=O)CC(=O)[O-]", "sulfoacetate"),
        // P-63.8.1 / P-65.6.2.1: cations, then the anion as a separate word.
        (
            "[Na+].[Na+].[O-]S(=O)(=O)CC(=O)[O-]",
            "disodium sulfonatoacetate",
        ),
        // P-72.6.1: sulfinato, phosphonato and arsonato.
        ("[O-]S(=O)CC(=O)[O-]", "sulfinatoacetate"),
        ("[O-]P(=O)([O-])CC(=O)[O-]", "phosphonatoacetate"),
        ("[O-][As](=O)([O-])CC(=O)[O-]", "arsonatoacetate"),
        // P-65.6.2.3.1: "and similarly for the corresponding selenium and
        // tellurium acids".
        ("[O-][Se](=O)(=O)CC(=O)[O-]", "selenonatoacetate"),
        ("[O-][Se](=O)CC(=O)[O-]", "seleninatoacetate"),
        // P-72.6.1 asks for a hydron removed from all hydroxy groups, so
        // -P(O)(OH)O(-) is not phosphonato: it is phosphoryl
        // (P-67.1.4.1.1.2) substituted as in hydroxy(sulfanyl)phosphoryl
        // (P-67.1.4.1.1.5), with oxido from P-72.6.2.
        (
            "OP(=O)([O-])CC(=O)[O-]",
            "[hydroxy(oxido)phosphoryl]acetate",
        ),
        (
            "OP(=O)([O-])c1ccc(cc1)C(=O)[O-]",
            "4-[hydroxy(oxido)phosphoryl]benzoate",
        ),
        (
            "O[As](=O)([O-])CC(=O)[O-]",
            "[hydroxy(oxido)arsoryl]acetate",
        ),
        // P-72.2.2.2.4: an anionic centre that no suffix expresses is named
        // from the anionic parent hydride, as in tert-butyltrioxidanide and
        // (ethylsulfanyl)oxidanide (P-72.2.2.2.2).
        ("[O-]SSC", "[(methylsulfanyl)sulfanyl]oxidanide"),
        ("[O-]SSc1ccccc1", "[(phenylsulfanyl)sulfanyl]oxidanide"),
        ("[O-]SSSC", "{[(methylsulfanyl)sulfanyl]sulfanyl}oxidanide"),
        ("[O-]SS", "(disulfanyl)oxidanide"),
        ("[O-][Se][Se]C", "[(methylselanyl)selanyl]oxidanide"),
        (
            "[Na+].[O-]SSC",
            "sodium [(methylsulfanyl)sulfanyl]oxidanide",
        ),
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
fn anionic_prefixes_without_a_readable_name_are_declined() {
    // P-72.6.1 gives no example for the chalcogen analogues of carboxylato,
    // and `carboxy` would drop the charge; two oxidanide centres on one chain
    // have no single anionic parent hydride.
    for smiles in [
        "[O-]C(=S)CC(=O)[O-]",
        "[S-]C(=S)CC(=O)[O-]",
        "[S-]C(=O)CC(=O)[O-]",
        "[O-]SS[O-]",
    ] {
        assert!(smiles_to_name(smiles).is_err(), "{smiles}");
    }
}

#[test]
fn selenium_and_tellurium_acid_anion_prefixes_are_read() {
    // P-65.6.2.3.1, P-72.6.1.
    let budget = Budget::default();
    for (name, smiles) in [
        ("selenonatoacetate", "[O-][Se](=O)(=O)CC(=O)[O-]"),
        ("4-seleninatobenzoate", "[O-][Se](=O)c1ccc(cc1)C(=O)[O-]"),
        ("telluronatoacetate", "[O-][Te](=O)(=O)CC(=O)[O-]"),
        ("tellurinatoacetate", "[O-][Te](=O)CC(=O)[O-]"),
    ] {
        let parsed = name_to_smiles(name).unwrap_or_else(|e| panic!("{name}: {e:?}"));
        assert_eq!(
            canonical_smiles(&parsed, &budget),
            canonical_smiles(smiles, &budget),
            "{name}"
        );
    }
}
