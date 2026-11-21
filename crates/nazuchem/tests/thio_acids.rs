// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | crates/nazuchem/tests/thio_acids.rs | issue #120

#![cfg(all(feature = "n2s", feature = "s2n"))]

//! Issue #120: thio diacids, thio acids cited as prefixes, multiplied
//! compound ester groups and -SSH prefixes are named after the IUPAC 2013
//! recommendations and read back as the structure they were generated from.

use nazuchem::{
    Budget, NameReason, NameStatus, name_to_smiles, s2n::canonical_smiles, smiles_to_name,
};

/// Issue #27. A cyano group on the S-organyl part does not change the
/// thiocarboxylate modifier's connectivity (P-65.2.1.5, P-65.6.3.3.7.1).
#[test]
fn acyl_carbonotrithioate_preserves_the_s_organyl_nitrile() {
    let budget = Budget::default();
    for smiles in ["Cc1ccc(C(=O)SC(=S)SCC#N)cc1", "N#CCSC(=S)SC(=O)c1ccc(C)cc1"] {
        let name = smiles_to_name(smiles).unwrap();
        assert_eq!(
            name.text,
            "S-{[(cyanomethyl)sulfanyl]carbonothioyl} 4-methylbenzene-1-carbothioate"
        );
        assert_eq!(
            name.status,
            NameStatus::Unverified(NameReason::FunctionalClass)
        );
        let parsed = name_to_smiles(&name.text).unwrap();
        assert_eq!(
            canonical_smiles(&parsed, &budget).unwrap(),
            canonical_smiles(smiles, &budget).unwrap()
        );
    }
}

#[test]
fn sulfur_acid_carbonotrithioate_anhydrides_keep_their_acid_citations() {
    let budget = Budget::default();
    for (smiles, expected) in [
        (
            "CS(=O)(=O)SC(=S)SC#N",
            "methanesulfonic thiocyanatomethanethioic thioanhydride",
        ),
        (
            "CS(=O)SC(=S)SC#N",
            "methanesulfinic thiocyanatomethanethioic thioanhydride",
        ),
    ] {
        let name = smiles_to_name(smiles).unwrap();
        assert_eq!(name.text, expected);
        assert_eq!(
            name.status,
            NameStatus::Unverified(NameReason::UnauditedParent)
        );
        let parsed = name_to_smiles(&name.text).unwrap();
        assert_eq!(
            canonical_smiles(&parsed, &budget).unwrap(),
            canonical_smiles(smiles, &budget).unwrap()
        );
    }
}

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
        // P-43.1 Table 4.3: -CS-OH is senior to -CO-SH; the prefix follows
        // P-65.1.5.1.
        ("SC(=O)CCC(=S)O", "4-oxo-4-sulfanylbutanethioic O-acid"),
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
        (
            "OC(=O)CCC(=O)[Se]O",
            "4-(hydroxyselanyl)-4-oxobutanoic acid",
        ),
        // -Te-OH keeps the compound prefix: its simple prefix is not
        // generated yet.
        (
            "OC(=O)CCC(=O)[Te]O",
            "3-[(hydroxytellanyl)carbonyl]propanoic acid",
        ),
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

/// Issue #184: of two acid groups the suffix goes to the senior one of P-43.1
/// Table 4.3, a -CO-OOH prefix off a chain end is carbonoperoxoyl, and the
/// chalcogen and peroxy analogues of carbonic acid keep the carbonic stem.
#[test]
fn acid_suffix_seniority_and_carbonic_analogues_follow_the_blue_book() {
    for (smiles, expected) in [
        // Table 4.3: -COOH precedes every modified suffix. P-65.2.1.6:
        // hydroxycarbonothioyl keeps the sulfur of -CS-OH.
        (
            "OC(=S)c1ccc(C(=O)O)cc1",
            "4-(hydroxycarbonothioyl)benzoic acid",
        ),
        // P-65.1.4.2 (1), P-65.2.1.5: carbonoperoxoyl, not
        // hydroperoxycarbonyl (2-carbonoperoxoylbenzoic acid).
        ("OOC(=O)c1ccccc1C(=O)O", "2-carbonoperoxoylbenzoic acid"),
        ("OC(=O)CC(CC)C(=O)OO", "3-carbonoperoxoylpentanoic acid"),
        (
            "OC(=O)c1ccc(cc1)C1CC(C(=O)OO)C1",
            "4-(3-carbonoperoxoylcyclobutyl)benzoic acid",
        ),
        (
            "OC(=O)c1ccc(OC(=O)OO)cc1",
            "4-(carbonoperoxoyloxy)benzoic acid",
        ),
        // Table 4.3: -CO-OOH (3O) precedes the thio acids.
        (
            "OOC(=O)c1ccc(C(=O)S)cc1",
            "4-(sulfanylcarbonyl)benzene-1-carboperoxoic acid",
        ),
        ("SC(=O)CCC(=O)OO", "4-oxo-4-sulfanylbutaneperoxoic acid"),
        (
            "OOC(=O)c1ccc(cc1)C(=S)O",
            "4-(hydroxycarbonothioyl)benzene-1-carboperoxoic acid",
        ),
        // Table 4.3: -CS-OOH (2O,1S; OO) precedes -CO-SOH (2O,1S; OS; OH),
        // which precedes -CO-OSH (2O,1S; OS; SH).
        (
            "OOC(=S)c1ccc(C(=O)OS)cc1",
            "4-[(sulfanyloxy)carbonyl]benzene-1-carboperoxothioic acid",
        ),
        (
            "OSC(=O)CCC(=O)OS",
            "4-oxo-4-(sulfanyloxy)butane(thioperoxoic) SO-acid",
        ),
        // Table 4.3: every peroxy acid, here -CO-OSH, precedes -CS-OH.
        (
            "SOC(=O)CCC(=S)O",
            "4-hydroxy-4-sulfanylidenebutane(thioperoxoic) OS-acid",
        ),
        // Table 4.3: -CS-OH (1O,1S; OH) precedes -CO-SH (1O,1S; SH).
        (
            "OC(=S)c1ccc(C(=O)S)cc1",
            "4-(sulfanylcarbonyl)benzene-1-carbothioic O-acid",
        ),
        ("OC(=S)CCC(=O)S", "4-oxo-4-sulfanylbutanethioic O-acid"),
        // Table 4.3: -CO-SH precedes -CS-SH (2S); P-43.1 (a): S before Se.
        ("SC(=S)CCC(=O)S", "3-(dithiocarboxy)propanethioic S-acid"),
        ("[SeH]C(=O)CCC(=O)S", "4-oxo-4-selanylbutanethioic S-acid"),
        // Table 4.3: -CSe-SH (1S,1Se; SH) precedes -CS-SeH (1S,1Se; SeH).
        (
            "[SeH]C(=S)CCC(=[Se])S",
            "4-selanyl-4-sulfanylidenebutaneselenothioic S-acid",
        ),
        // P-65.2.1.2: carbonoperoxoic acid, carbonodiperoxoic acid,
        // carbonothioic S-acid, carbonodiselenoic Se,Se-acid,
        // carbonotrithioic acid. P-65.2.1.5: carbonothioic O,O-acid.
        ("O=C(O)OO", "carbonoperoxoic acid"),
        ("OOC(=O)OO", "carbonodiperoxoic acid"),
        ("O=C(O)OS", "carbono(thioperoxoic) OS-acid"),
        ("O=C(O)SO", "carbono(thioperoxoic) SO-acid"),
        ("OSC(=O)OS", "carbonobis(thioperoxoic) OS,SO-acid"),
        ("O=C(OS)OS", "carbonobis(thioperoxoic) OS,OS-acid"),
        ("O=C(SO)SO", "carbonobis(thioperoxoic) SO,SO-acid"),
        ("O=C(O)O[SeH]", "carbono(selenoperoxoic) OSe-acid"),
        (
            "O[Te]C(=O)O[TeH]",
            "carbonobis(telluroperoxoic) OTe,TeO-acid",
        ),
        ("O=C(O)S", "carbonothioic S-acid"),
        ("O=C(O)[SeH]", "carbonoselenoic Se-acid"),
        ("OC(=S)O", "carbonothioic O,O-acid"),
        ("O=C(S)S", "carbonodithioic S,S-acid"),
        ("O=C([SeH])[SeH]", "carbonodiselenoic Se,Se-acid"),
        ("SC(=S)S", "carbonotrithioic acid"),
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
    // The names on methane for other carbonic replacements are unverified.
    for smiles in ["O=C(S)OO", "OC(=S)S", "O=C(S)[SeH]"] {
        let name = smiles_to_name(smiles).unwrap();
        assert_ne!(name.status, NameStatus::Preferred, "{smiles} {}", name.text);
    }
}

#[test]
fn carbonic_mixed_peroxy_acid_locants_match_the_suffix() {
    for name in [
        "carbono(thioperoxoic) acid",
        "carbono(thioperoxoic) SeO-acid",
        "carbono(thioperoxoic) OS,SO-acid",
        "carbonobis(thioperoxoic) OS-acid",
        "carbonobis(thioperoxoic) OS,SO,SO-acid",
        "phosphorobis(thioperoxoic) OS,SO-acid",
    ] {
        assert!(name_to_smiles(name).is_err(), "{name}");
    }
}

#[test]
fn chain_end_acid_halides_and_esters_cite_separate_ligands() {
    for (smiles, expected) in [
        // P-65.5.4 (3): the halide and oxo at the end of the parent chain.
        ("O=C(O)CC(=O)Cl", "3-chloro-3-oxopropanoic acid"),
        ("O=C(O)CC(=O)C(=O)Br", "4-bromo-3,4-dioxobutanoic acid"),
        (
            "O=C(O)CC(=O)OC(=O)Br",
            "3-(carbonobromidoyloxy)-3-oxopropanoic acid",
        ),
        ("O=C(O)CC(=S)Cl", "3-chloro-3-sulfanylidenepropanoic acid"),
        // The same form can occur in a substituent chain.
        (
            "O=C(O)c1ccc(CC(=O)Cl)cc1",
            "4-(2-chloro-2-oxoethyl)benzoic acid",
        ),
        // P-65.6.3.2.3 and P-65.6.3.3.5: neutral partial ester and its salt.
        ("CCOC(=O)CCC(=O)O", "4-ethoxy-4-oxobutanoic acid"),
        (
            "CCOC(=O)CCC(=O)[O-].[Na+]",
            "sodium 4-ethoxy-4-oxobutanoate",
        ),
        ("O=C(O)CCC(=O)SCC", "4-(ethylsulfanyl)-4-oxobutanoic acid"),
        ("COOC(=O)CC(=O)O", "3-(methylperoxy)-3-oxopropanoic acid"),
        (
            "CC(C)OOC(=O)CC(=O)O",
            "3-oxo-3-[(propan-2-yl)peroxy]propanoic acid",
        ),
        // P-65.6.3.3.7.1: the ordinary ester takes the suffix and
        // the chain-end thiono ester cites its ligands separately.
        (
            "CCCCOC(=S)CCC(=O)OCC",
            "ethyl 4-butoxy-4-sulfanylidenebutanoate",
        ),
        ("CCOC(=O)c1ccccc1C(=O)O", "2-(ethoxycarbonyl)benzoic acid"),
        // P-65.6.3.3.4.2: the attachment carbonyl appears in propanoyl.
        (
            "CCOC(=O)CC(=O)Oc1ccc(OC(=O)CCC(=O)OC)cc1",
            "4-[(3-ethoxy-3-oxopropanoyl)oxy]phenyl methyl butanedioate",
        ),
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

#[test]
fn ester_organyl_locants_are_omitted_when_acid_sites_are_interchangeable() {
    for (smiles, expected) in [
        ("CCOC(=O)CCC(=O)OC", "ethyl methyl butanedioate"),
        ("CCSC(=O)CCC(=O)SC", "S-ethyl S-methyl butanebis(thioate)"),
        (
            "CCOSC(=O)CCC(=O)SOC",
            "SO-ethyl SO-methyl butanebis(thioperoxoate)",
        ),
        (
            "CCOC(=O)c1cc(C(=O)OCC)cc(C(=O)OC)c1",
            "diethyl methyl benzene-1,3,5-tricarboxylate",
        ),
        (
            "CCOC(=O)c1cc(C(=O)OC)c(C(=O)OC)cc1C(=O)OC",
            "ethyl trimethyl benzene-1,2,4,5-tetracarboxylate",
        ),
        // A chloro substituent breaks the exchange symmetry.
        (
            "CCSC(=O)CC(Cl)C(=O)SC",
            "4-S-ethyl 1-S-methyl 2-chlorobutanebis(thioate)",
        ),
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
    // Two ethyl groups may sit at several distinct relative positions.
    let smiles = "CCOC(=O)c1cc(C(=O)OCC)c(C(=O)OC)cc1C(=O)OC";
    let name = smiles_to_name(smiles).unwrap();
    assert!(
        name.text.starts_with(|c: char| c.is_ascii_digit()),
        "{name:?}"
    );
    let parsed = name_to_smiles(&name.text).unwrap();
    assert_eq!(
        canonical_smiles(&parsed, &Budget::default()),
        canonical_smiles(smiles, &Budget::default())
    );
}
