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
        // Issue #180, P-72.7(a): two sulfonate centres outrank the one
        // carboxylate, which is then carboxylato.
        (
            "[O-]C(=O)CC(S(=O)(=O)[O-])CS(=O)(=O)[O-]",
            "3-carboxylatopropane-1,2-disulfonate",
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
    // prefix. P-72.7(a) ties at one anionic centre each, and P-72.7(b)
    // prefers the parent with more `ide` centres (diide is senior to
    // idoate), so the aminide (P-72.2.2.2.3) is the suffix. P-72.6.1: the
    // acid anion is carboxylato, sulfonato or phosphonato, and the
    // carboxylato carbon is not a chain atom. P-72.2.2.2.3: two aminide
    // suffixes take "bis".
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
        // Issue #180. P-72.7(a): the parent with more anionic centres,
        // suffixes included, is senior, so two carboxylates outrank one
        // aminide, which is then azanidyl (P-72.6.3).
        (
            "[NH-]CC(C(=O)[O-])C(=O)[O-]",
            "(azanidylmethyl)propanedioate",
            true,
        ),
        (
            "[NH-]C(CC(=O)[O-])CC(=O)[O-]",
            "3-azanidylpentanedioate",
            true,
        ),
        // P-72.7(a), P-72.5.3: the carbanion stays in the dicarboxylate.
        (
            "[NH-]C([CH-]C(=O)[O-])C(=O)[O-]",
            "3-azanidylbutan-2-ide-1,4-dioate",
            true,
        ),
        // P-72.7(a) ties at two centres each; P-72.7(b) keeps the aminides.
        (
            "[NH-]CC([NH-])CC(C(=O)[O-])C(=O)[O-]",
            "4,4-bis(carboxylato)butane-1,2-bis(aminide)",
            true,
        ),
        // A phosphonate has two anionic centres but is no suffix on carbon
        // here, so this aminide name is not labelled preferred.
        (
            "[NH-]CCP(=O)([O-])[O-]",
            "2-phosphonatoethan-1-aminide",
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

#[test]
fn anionic_parent_has_the_most_anionic_centres() {
    // Issue #180. P-72.7(a): the parent with the maximum number of anionic
    // centres, anionic suffixes included, as in
    // 1-(borinan-1-uid-4-yl)ethane-1,2-bis(olate), not
    // 4-(1,2-dioxidoethyl)borinan-1-uide. The acid anion beside the senior
    // suffix is carboxylato or sulfonato (P-72.6.1).
    for (smiles, expected, preferred) in [
        (
            "[O-]CC([O-])C(=O)[O-]",
            "1-carboxylatoethane-1,2-bis(olate)",
            true,
        ),
        (
            "[O-]CC([O-])CC(=O)[O-]",
            "3-carboxylatopropane-1,2-bis(olate)",
            true,
        ),
        (
            "[O-]c1cc([O-])cc(c1)C(=O)[O-]",
            "5-carboxylatobenzene-1,3-bis(olate)",
            true,
        ),
        (
            "[S-]CC([S-])C(=O)[O-]",
            "1-carboxylatoethane-1,2-bis(thiolate)",
            true,
        ),
        (
            "[O-]S(=O)(=O)CC(S(=O)(=O)[O-])C(=O)[O-]",
            "1-carboxylatoethane-1,2-disulfonate",
            true,
        ),
        (
            "[O-]CC([O-])S(=O)(=O)[O-]",
            "1-sulfonatoethane-1,2-bis(olate)",
            true,
        ),
        // P-72.7(a) ties, so P-72.7(e) decides: carboxylate is senior to
        // olate, as in 3-oxidonaphthalene-2-carboxylate.
        ("[O-]CCC(=O)[O-]", "3-oxidopropanoate", true),
        ("[O-]c1ccccc1C(=O)[O-]", "2-oxidobenzoate", true),
        (
            "[O-]C(=O)C([O-])C([O-])C(=O)[O-]",
            "2,3-dioxidobutanedioate",
            true,
        ),
        ("[O-]CC([S-])C(=O)[O-]", "3-oxido-2-sulfidopropanoate", true),
        // The olates are on different parents, one centre each.
        (
            "[O-]c1ccc(cc1)CC(C(=O)[O-])Cc1ccc([O-])cc1",
            "3-(4-oxidophenyl)-2-[(4-oxidophenyl)methyl]propanoate",
            true,
        ),
        // Three olates cannot share one chain, so the name is not certified.
        (
            "[O-]CC(C[O-])(C[O-])CC(=O)[O-]",
            "2-(carboxylatomethyl)-2-(oxidomethyl)propane-1,3-bis(olate)",
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

#[test]
fn amide_nitrogen_anions_are_azanides() {
    // Issue #180. P-72.2.2.2.4: an anionic centre on an amide nitrogen is
    // named from the anionic parent hydride, acetylazanide; `amidide` and
    // `carboxamidide` are not recommended, and `acetamide` drops the charge.
    for (smiles, expected, preferred) in [
        ("[NH-]C(=O)C", "acetylazanide", true),
        ("[NH-]C=O", "formylazanide", true),
        ("[NH-]C(=O)c1ccccc1", "benzoylazanide", true),
        ("C[N-]C(=O)C", "acetyl(methyl)azanide", true),
        ("CC(=O)[N-]c1ccccc1", "acetyl(phenyl)azanide", true),
        ("[N-](C(=O)C)C(=O)C", "diacetylazanide", true),
        ("[NH-]C(=S)C", "(ethanethioyl)azanide", true),
        ("[Na+].[NH-]C(=O)C", "sodium acetylazanide", true),
        // A ligand with its own characteristic group is not audited.
        ("[NH-]C(=O)CCO", "(3-hydroxypropanoyl)azanide", false),
        ("[NH-]C(=O)OC", "(methoxycarbonyl)azanide", false),
        // P-72.6.3: beside a senior anion the nitrogen is azanidyl.
        ("[NH-]C(=O)CC(=O)[O-]", "3-azanidyl-3-oxopropanoate", true),
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
    // A neutral amide suffix would drop the charge beside a nitro group.
    assert!(smiles_to_name("[NH-]C(=O)c1ccc(cc1)[N+](=O)[O-]").is_err());
}

#[test]
fn peroxy_acid_anion_prefix_reads_back() {
    // Issue #180. P-72.6.3: -O-O(-) is dioxidanidyl, formed like
    // disulfanidyl, here on carbonyl. P-72.6.1 has no example of an `ato`
    // prefix for a peroxy acid anion, so the names are not labelled
    // preferred. P-72.7(b): the aminide is the suffix.
    for (smiles, expected) in [
        ("[NH-]CC(=O)O[O-]", "1-(dioxidanidylcarbonyl)methanaminide"),
        ("[O-]OC(=O)CC(=O)[O-]", "(dioxidanidylcarbonyl)acetate"),
        (
            "[O-]OC(=O)c1ccccc1C(=O)[O-]",
            "2-(dioxidanidylcarbonyl)benzoate",
        ),
    ] {
        let name = smiles_to_name(smiles).unwrap();
        assert_eq!(name.text, expected, "{smiles}");
        assert_ne!(name.status, NameStatus::Preferred, "{smiles}");
        let parsed = name_to_smiles(&name.text).unwrap_or_else(|e| panic!("{expected}: {e:?}"));
        let budget = Budget::default();
        assert_eq!(
            canonical_smiles(&parsed, &budget),
            canonical_smiles(smiles, &budget),
            "{smiles} -> {expected} -> {parsed}"
        );
    }
    // The chalcogen analogues have no validated prefix.
    assert!(smiles_to_name("[S-]OC(=O)CC(=O)[O-]").is_err());
}
