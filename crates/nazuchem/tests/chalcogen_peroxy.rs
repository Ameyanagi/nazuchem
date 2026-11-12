// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | crates/nazuchem/tests/chalcogen_peroxy.rs | issue #94

#![cfg(all(feature = "n2s", feature = "s2n"))]

//! Issue #94: names of chalcogen peroxy compounds and of multiplied
//! chalcogen esters follow the IUPAC 2013 recommendations and read back as
//! the structure they were generated from.

use nazuchem::{Budget, name_to_smiles, s2n::canonical_smiles, smiles_to_name};

#[test]
fn chalcogen_peroxy_names_follow_the_blue_book_and_round_trip() {
    for (smiles, expected) in [
        // P-63.4.2.1, Table 6.1: methane-SO-thioperoxol and
        // propane-1-OS-thioperoxol are the Blue Book examples.
        ("CSO", "methane-SO-thioperoxol"),
        ("CO[SH]", "methane-OS-thioperoxol"),
        ("CCCOS", "propane-1-OS-thioperoxol"),
        ("OSc1ccccc1", "benzene-SO-thioperoxol"),
        ("CCCO[SeH]", "propane-1-OSe-selenoperoxol"),
        // P-72.2.2.2.2 (benzene-1,4-bis(OS-thioperoxolate)): a suffix with a
        // letter locant takes bis.
        ("OSCCSO", "ethane-1,2-bis(SO-thioperoxol)"),
        ("SOc1ccc(OS)cc1", "benzene-1,4-bis(OS-thioperoxol)"),
        // P-63.4.2.2: as a prefix the group stays hydroxysulfanyl or
        // sulfanyloxy (3-(sulfanyloxy)propanenitrile).
        ("OSCCO", "2-(hydroxysulfanyl)ethan-1-ol"),
        ("SOCCC#N", "3-(sulfanyloxy)propanenitrile"),
        // P-65.1.5.3: SO and OS are needed; two sulfur atoms leave no
        // choice, so (dithioperoxoic) takes no locant.
        ("CC(=O)SO", "ethane(thioperoxoic) SO-acid"),
        ("CC(=O)OS", "ethane(thioperoxoic) OS-acid"),
        ("CC(=O)SS", "ethane(dithioperoxoic) acid"),
        ("CS(=O)(=O)SO", "methanesulfono(thioperoxoic) SO-acid"),
        ("CS(=O)(=O)SS", "methanesulfono(dithioperoxoic) acid"),
        ("CS(=O)SS", "methanesulfino(dithioperoxoic) acid"),
        // P-65.2.1.2 (carbonobis(thioperoxoic) OS,SO-acid): bis before the
        // bracketed suffix and one locant for each group.
        ("OSC(=O)CCC(=O)SO", "butanebis(thioperoxoic) SO,SO-acid"),
        ("O=C(OS)CCC(=O)OS", "butanebis(thioperoxoic) OS,OS-acid"),
        ("SSC(=O)CCC(=O)SS", "butanebis(dithioperoxoic) acid"),
        (
            "OSC(=O)C1CCC(CC1)C(=O)SO",
            "cyclohexane-1,4-dicarbo(thioperoxoic) SO,SO-acid",
        ),
        (
            "OSS(=O)(=O)CCS(=O)(=O)SO",
            "ethane-1,2-disulfono(thioperoxoic) SO,SO-acid",
        ),
        // P-72.2.2.2.1.1 (ethane(OS-thioperoxoate)): the anion carries its
        // locant inside the brackets.
        ("CC(=O)S[O-]", "ethane(SO-thioperoxoate)"),
        ("CC(=O)O[S-]", "ethane(OS-thioperoxoate)"),
        ("CC(=O)S[S-]", "ethane(dithioperoxoate)"),
        ("CC(=O)S[Se-]", "ethane(SSe-selenothioperoxoate)"),
        ("CC(=O)S[O-].[Na+]", "sodium ethane(SO-thioperoxoate)"),
        ("[O-]SC(=O)CCC(=O)S[O-]", "butanebis(SO-thioperoxoate)"),
        ("CS(=O)(=O)S[O-]", "methanesulfono(SO-thioperoxoate)"),
        ("[O-]SC(=O)c1ccccc1", "benzenecarbo(SO-thioperoxoate)"),
        // P-65.6.3.3.7.1: SO-methyl benzene(carbothioperoxoate) and
        // S,S-dimethyl ethanebis(thioate); P-16.3.5 (b) for bis(thioate).
        ("CC(=O)SOC", "SO-methyl ethane(thioperoxoate)"),
        ("CSC(=O)C(=O)SC", "S,S-dimethyl ethanebis(thioate)"),
        ("COC(=S)CCC(=S)OC", "O,O-dimethyl butanebis(thioate)"),
        ("CSC(=S)CCC(=S)SC", "S,S-dimethyl butanebis(dithioate)"),
        (
            "CSC(=O)c1ccc(cc1)C(=O)SC",
            "S,S-dimethyl benzene-1,4-dicarbothioate",
        ),
        (
            "COSC(=O)CCC(=O)SOC",
            "SO,SO-dimethyl butanebis(thioperoxoate)",
        ),
        (
            "CSSC(=O)CCC(=O)SSC",
            "SS,SS-dimethyl butanebis(dithioperoxoate)",
        ),
        ("COOC(=S)CCC(=S)OOC", "OO,OO-dimethyl butanediperoxothioate"),
        (
            "COSC(=O)c1ccc(cc1)C(=O)SOC",
            "SO,SO-dimethyl benzene-1,4-dicarbo(thioperoxoate)",
        ),
        // P-65.6.3.3.7.2.1 (3-ethyl 1-S-methyl 1-thiodicarbonate): the
        // position of the acid group, then the letter locant.
        (
            "CCSC(=O)CC(Cl)C(=O)SC",
            "4-S-ethyl 1-S-methyl 2-chlorobutanebis(thioate)",
        ),
        (
            "CCOSC(=O)CC(Cl)C(=O)SOC",
            "4-SO-ethyl 1-SO-methyl 2-chlorobutanebis(thioperoxoate)",
        ),
        // P-65.1.5.3 (4-[(hydroxysulfanyl)carbonyl]cyclohexane-1-carboxylic
        // acid): a peroxy acid cited as a prefix names its own chalcogens.
        (
            "O=C(O)C1CCC(CC1)C(=O)SO",
            "4-[(hydroxysulfanyl)carbonyl]cyclohexane-1-carboxylic acid",
        ),
        (
            "O=C(O)C1CCC(CC1)C(=S)OS",
            "4-[(sulfanyloxy)carbonothioyl]cyclohexane-1-carboxylic acid",
        ),
        (
            "O=C(O)C1CCC(CC1)C(=O)SS",
            "4-(disulfanylcarbonyl)cyclohexane-1-carboxylic acid",
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
