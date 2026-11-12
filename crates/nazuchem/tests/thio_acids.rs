// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | crates/nazuchem/tests/thio_acids.rs | issue #120

#![cfg(all(feature = "n2s", feature = "s2n"))]

//! Issue #120: thio diacids, thio acids cited as prefixes, multiplied
//! compound ester groups and -SSH prefixes are named after the IUPAC 2013
//! recommendations and read back as the structure they were generated from.

use nazuchem::{Budget, name_to_smiles, s2n::canonical_smiles, smiles_to_name};

#[test]
fn thio_acid_and_ester_names_follow_the_blue_book_and_round_trip() {
    for (smiles, expected) in [
        // P-16.3.5 (b), P-65.1.5.1 (butanebis(thioic acid)): bis encloses
        // the whole suffix, and the letter locant stays with each group.
        ("SC(=O)CCC(=O)S", "butanebis(thioic S-acid)"),
        ("OC(=S)CCC(=S)O", "butanebis(thioic O-acid)"),
        ("SC(=O)C(=O)S", "ethanebis(thioic S-acid)"),
        ("[SeH]C(=O)CCC(=O)[SeH]", "butanebis(selenoic Se-acid)"),
        ("SC(=[Se])CCC(=[Se])S", "butanebis(selenothioic S-acid)"),
        // P-16.3.5 (b): bis(dithioic acid), not didithioic acid.
        ("SC(=S)CCC(=S)S", "butanebis(dithioic acid)"),
        // P-65.1.5.2 (benzene-1,2-dicarbothioic acid): carbo forms keep di.
        ("SC(=O)c1ccc(cc1)C(=O)S", "benzene-1,4-dicarbothioic S-acid"),
        // P-65.2.1.6: sulfanylcarbonyl, hydroxycarbonothioyl and
        // dithiocarboxy; P-65.1.5.1 (3-(sulfanylcarbonyl)propanoic acid,
        // 4-(sulfanylcarbonyl)pyridine-2-carboxylic acid). On a chain the
        // PIN is the chain-end form, 4-oxo-4-sulfanylbutanoic acid, which
        // the namer does not generate for any acyl prefix yet.
        ("O=C(O)CCC(=O)S", "3-(sulfanylcarbonyl)propanoic acid"),
        ("O=C(O)CCC(=S)O", "3-(hydroxycarbonothioyl)propanoic acid"),
        ("O=C(O)CCC(=S)S", "3-(dithiocarboxy)propanoic acid"),
        ("O=C(O)c1ccc(C(=O)S)cc1", "4-(sulfanylcarbonyl)benzoic acid"),
        (
            "O=C(O)c1ccc(C(=S)O)cc1",
            "4-(hydroxycarbonothioyl)benzoic acid",
        ),
        // P-65.1.5.2: 4-(selanylcarbonyl)benzoic acid (PIN).
        (
            "O=C(O)c1ccc(C(=O)[SeH])cc1",
            "4-(selanylcarbonyl)benzoic acid",
        ),
        (
            "SC(=O)CCC(=S)O",
            "3-(hydroxycarbonothioyl)propanethioic S-acid",
        ),
        // P-16.3.5 (a) (bis(bromomethyl)): a compound organyl group takes bis.
        ("ClCOC(=O)C(=O)OCCl", "bis(chloromethyl) ethanedioate"),
        ("ClCOC(=O)OCCl", "bis(chloromethyl) carbonate"),
        (
            "ClCSC(=O)C(=O)SCCl",
            "S,S-bis(chloromethyl) ethanebis(thioate)",
        ),
        ("CCOC(=O)C(=O)OCC", "diethyl ethanedioate"),
        // P-63.4.2.2 (disulfanylacetic acid, 3,4-bis(disulfanyl)benzamide).
        ("SSCC(=O)O", "disulfanylacetic acid"),
        ("SSCCC(=O)O", "3-(disulfanyl)propanoic acid"),
        ("SSC(CO)SS", "2,2-bis(disulfanyl)ethan-1-ol"),
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
