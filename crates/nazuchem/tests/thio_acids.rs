// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | crates/nazuchem/tests/thio_acids.rs | issue #120

#![cfg(all(feature = "n2s", feature = "s2n"))]

//! Issue #120: thio diacids, thio acids cited as prefixes, multiplied
//! compound ester groups and -SSH prefixes are named after the IUPAC 2013
//! recommendations and read back as the structure they were generated from.

use nazuchem::{Budget, NameStatus, name_to_smiles, s2n::canonical_smiles, smiles_to_name};

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
        // dithiocarboxy; P-65.1.5.1 (4-(sulfanylcarbonyl)pyridine-2-carboxylic
        // acid). Ending a chain, a mixed acid is cited with simple prefixes;
        // see the next test.
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

/// Issue #138: a peroxy acid or a mixed chalcogen acid cited as a prefix
/// keeps its carbon in the chain it ends and names its own atoms.
#[test]
fn acid_prefixes_ending_a_chain_follow_the_blue_book_and_round_trip() {
    for (smiles, expected) in [
        // P-65.1.5.1: sulfanyl with oxo and hydroxy with sulfanylidene
        // (4-oxo-4-sulfanylbutanoic acid, 4-hydroxy-4-sulfanylidenebutanoic
        // acid, oxo(sulfanyl)acetic acid); P-65.1.5.2
        // (hydroxy(sulfanylidene)acetic acid).
        ("O=C(O)CCC(=O)S", "4-oxo-4-sulfanylbutanoic acid"),
        ("O=C(O)CCC(=S)O", "4-hydroxy-4-sulfanylidenebutanoic acid"),
        ("O=C(O)C(=O)S", "oxo(sulfanyl)acetic acid"),
        ("O=C(O)C(=S)O", "hydroxy(sulfanylidene)acetic acid"),
        ("O=C(O)CC(=O)[SeH]", "3-oxo-3-selanylpropanoic acid"),
        (
            "OC(=O)CCC(=S)[SeH]",
            "4-selanyl-4-sulfanylidenebutanoic acid",
        ),
        ("OC(=O)C=CC(=O)S", "4-oxo-4-sulfanylbut-2-enoic acid"),
        (
            "[Na+].[O-]C(=O)CCC(=O)S",
            "sodium 4-oxo-4-sulfanylbutanoate",
        ),
        (
            "SC(=O)CCC(=S)O",
            "4-hydroxy-4-sulfanylidenebutanethioic S-acid",
        ),
        // P-65.1.4.2 (2): hydroperoxy with oxo (6-hydroperoxy-6-oxohexanoic
        // acid).
        ("O=C(O)CCCCC(=O)OO", "6-hydroperoxy-6-oxohexanoic acid"),
        ("OC(=O)C(=O)OO", "hydroperoxy(oxo)acetic acid"),
        (
            "O=C(O)[C@H](O)[C@@H](O)C(=O)OO",
            "(2R,3R)-4-hydroperoxy-2,3-dihydroxy-4-oxobutanoic acid",
        ),
        // P-65.1.5.3: sulfanyloxy or hydroxysulfanyl with oxo or
        // sulfanylidene (4-sulfanylidene-4-(sulfanyloxy)butanoic acid,
        // (hydroxysulfanyl)oxoacetic acid).
        (
            "O=C(O)CCC(=S)OS",
            "4-sulfanylidene-4-(sulfanyloxy)butanoic acid",
        ),
        ("O=C(O)C(=O)SO", "(hydroxysulfanyl)(oxo)acetic acid"),
        ("O=C(O)CCC(=O)SO", "4-(hydroxysulfanyl)-4-oxobutanoic acid"),
        ("O=C(O)CCC(=O)OS", "4-oxo-4-(sulfanyloxy)butanoic acid"),
        (
            "O=C(O)CCC(=S)OO",
            "4-hydroperoxy-4-sulfanylidenebutanoic acid",
        ),
        // P-63.4.2.2, P-63.3.1: -SSH is disulfanyl, -SeSeH diselanyl.
        ("O=C(O)CCC(=O)SS", "4-(disulfanyl)-4-oxobutanoic acid"),
        ("OC(=O)CCC(=O)[Se][SeH]", "4-(diselanyl)-4-oxobutanoic acid"),
        ("[SeH][Se]CC(=O)O", "diselanylacetic acid"),
        // P-16.5.1.3.2: cited after another unlocanted prefix, a prefix
        // that begins with a multiplying term is enclosed.
        ("ClC(SS)C(=O)O", "chloro(disulfanyl)acetic acid"),
        ("NC([Se][SeH])C(=O)O", "amino(diselanyl)acetic acid"),
        ("ClC([Te][TeH])C(=O)O", "chloro(ditellanyl)acetic acid"),
        ("OC(=O)C(=O)SS", "disulfanyl(oxo)acetic acid"),
        // P-45.2.1: of two chains of equal length, the one through the acid
        // carbon has more substituents cited as prefixes.
        (
            "OC(=O)CC(C)C(=O)S",
            "3-methyl-4-oxo-4-sulfanylbutanoic acid",
        ),
        (
            "OC(=O)CC(C)C(=O)OO",
            "4-hydroperoxy-3-methyl-4-oxobutanoic acid",
        ),
        (
            "OC(=O)CC(CC(=O)S)CC(=O)S",
            "5-oxo-3-(2-oxo-2-sulfanylethyl)-5-sulfanylpentanoic acid",
        ),
        (
            "O=C(O)c1ccc(CC(C)C(=O)S)cc1",
            "4-(2-methyl-3-oxo-3-sulfanylpropyl)benzoic acid",
        ),
        // P-65.1.5.1, P-65.2.1.6: off the chain end the compound prefix
        // stays, on a nonterminal chain atom (a longer chain wins) and on a
        // ring (4-(sulfanylcarbonyl)pyridine-2-carboxylic acid).
        ("OC(=O)CC(CC)C(=O)S", "3-(sulfanylcarbonyl)pentanoic acid"),
        (
            "OC(=O)CC(CC(O)=O)C(=O)S",
            "3-(sulfanylcarbonyl)pentanedioic acid",
        ),
        (
            "OC(=O)c1cc(ccn1)C(=O)S",
            "4-(sulfanylcarbonyl)pyridine-2-carboxylic acid",
        ),
        // P-65.1.2.2.3, P-65.2.1.6: -COOH stays carboxy and -CS-SH
        // dithiocarboxy; neither needs the position of a chalcogen.
        (
            "O=C(O)CC(CC(=O)O)CC(=O)O",
            "3-(carboxymethyl)pentanedioic acid",
        ),
        ("OC(=O)CC(C)C(=S)S", "3-(dithiocarboxy)butanoic acid"),
    ] {
        let name = smiles_to_name(smiles).unwrap();
        assert_eq!(name.text, expected, "{smiles}");
        assert_eq!(name.status, NameStatus::Preferred, "{smiles}");
        let parsed = name_to_smiles(&name.text).unwrap_or_else(|e| panic!("{expected}: {e:?}"));
        let budget = Budget::default();
        assert_eq!(
            canonical_smiles(&parsed, &budget),
            canonical_smiles(smiles, &budget),
            "{smiles} -> {expected} -> {parsed}"
        );
    }
}
