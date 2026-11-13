// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | crates/nazuchem/tests/anionic_prefixes.rs | issue #145

#![cfg(all(feature = "n2s", feature = "s2n"))]

//! Issue #145: acid-anion prefixes under a more senior anion and anionic
//! oxygen on a chalcogen chain keep their charge, are named after the IUPAC
//! 2013 recommendations and read back as the structure they came from.

use nazuchem::{Budget, NameStatus, name_to_smiles, s2n::canonical_smiles, smiles_to_name};

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
        (
            "O[Sb](=O)([O-])CC(=O)[O-]",
            "[hydroxy(oxido)stiboryl]acetate",
        ),
        // P-65.6.3.2.3: an ester cited as a prefix keeps its own name.
        ("COS(=O)(=O)CC(=O)[O-]", "(methoxysulfonyl)acetate"),
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
    // or a second charge beside the oxidanide have no single anionic parent
    // hydride.
    for smiles in [
        "[O-]C(=S)CC(=O)[O-]",
        "[S-]C(=S)CC(=O)[O-]",
        "[S-]C(=O)CC(=O)[O-]",
        "[O-]SS[O-]",
        "[O-]S(C)([O-])SC",
        "[Na+].[Na+].[O-]S(C)([O-])SC",
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

#[test]
fn aminide_with_an_acid_anion_cites_the_acid_anion_as_a_prefix() {
    // Issue #173. P-72.6: one anion is the parent anion and the other a
    // prefix. P-72.7(d): nitrogen precedes oxygen among anionic centres, so
    // the aminide (P-72.2.2.2.3) is the suffix. P-72.6.1: the acid anion is
    // carboxylato, sulfonato or phosphonato, and the carboxylato carbon is
    // not a chain atom. P-72.2.2.2.3: two aminide suffixes take "bis".
    for (smiles, expected, preferred) in [
        ("[NH-]CC(=O)[O-]", "1-carboxylatomethanaminide", true),
        ("[NH-]CCC(=O)[O-]", "2-carboxylatoethan-1-aminide", true),
        ("[NH-]C(C)C(=O)[O-]", "1-carboxylatoethan-1-aminide", true),
        (
            "C[N-]CC(=O)[O-]",
            "1-carboxylato-N-methylmethanaminide",
            true,
        ),
        (
            "[Na+].[Na+].[NH-]CC(=O)[O-]",
            "disodium 1-carboxylatomethanaminide",
            true,
        ),
        (
            "[NH-]c1ccc(cc1)C(=O)[O-]",
            "4-carboxylatobenzen-1-aminide",
            true,
        ),
        ("[NH-]CCS(=O)(=O)[O-]", "2-sulfonatoethan-1-aminide", true),
        ("[NH-]CC[O-]", "2-oxidoethan-1-aminide", true),
        ("[NH-]CC[NH-]", "ethane-1,2-bis(aminide)", true),
        (
            "[NH-]CC([NH-])C(=O)[O-]",
            "1-carboxylatoethane-1,2-bis(aminide)",
            true,
        ),
        // P-72.7(a): a parent with more anionic centres is senior, so these
        // aminide names are not labelled preferred.
        (
            "[NH-]CC(C(=O)[O-])C(=O)[O-]",
            "2,2-bis(carboxylato)ethan-1-aminide",
            false,
        ),
        (
            "[NH-]C(CC(=O)[O-])CC(=O)[O-]",
            "1,3-bis(carboxylato)propan-2-aminide",
            false,
        ),
        (
            "[NH-]CCP(=O)([O-])[O-]",
            "2-phosphonatoethan-1-aminide",
            false,
        ),
        // The carbanion could stay in a dicarboxylate parent (P-72.5.3).
        (
            "[NH-]C([CH-]C(=O)[O-])C(=O)[O-]",
            "1,2-bis(carboxylato)ethan-1-id-2-aminide",
            false,
        ),
        // P-61.5.1: nitro is a neutral prefix, not an anionic centre.
        (
            "[NH-]C([N+](=O)[O-])C(=O)[O-]",
            "1-carboxylato-1-nitromethanaminide",
            true,
        ),
        // An aminide on sulfur is no suffix; the carboxylate stays the parent.
        (
            "[NH-]S(=O)(=O)CC(=O)[O-]",
            "(azanidylsulfonyl)acetate",
            false,
        ),
    ] {
        let name = smiles_to_name(smiles).unwrap();
        assert_eq!(name.text, expected, "{smiles}");
        assert_eq!(
            name.status == NameStatus::Preferred,
            preferred,
            "{smiles}: {:?}",
            name.status
        );
        let parsed = name_to_smiles(&name.text).unwrap_or_else(|e| panic!("{expected}: {e:?}"));
        let budget = Budget::default();
        assert_eq!(
            canonical_smiles(&parsed, &budget),
            canonical_smiles(smiles, &budget),
            "{smiles} -> {expected} -> {parsed}"
        );
    }
}
