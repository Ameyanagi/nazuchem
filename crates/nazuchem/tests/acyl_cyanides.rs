// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | crates/nazuchem/tests/acyl_cyanides.rs | issue #163

#![cfg(all(feature = "n2s", feature = "s2n"))]

//! Issue #163: a cyanide or other pseudohalide on a sulfur, selenium or
//! tellurium acyl group is named like the halide, phosphoric acid takes the
//! cyanido infix, and a multiplied ester anion reads back. Names follow the
//! IUPAC 2013 recommendations and read back as the structure they came from.

use nazuchem::{Budget, NameStatus, name_to_smiles, s2n::canonical_smiles, smiles_to_name};

fn assert_round_trip(smiles: &str, name: &str) {
    let parsed = name_to_smiles(name).unwrap_or_else(|e| panic!("{name}: {e:?}"));
    let budget = Budget::default();
    assert_eq!(
        canonical_smiles(&parsed, &budget),
        canonical_smiles(smiles, &budget),
        "{smiles} -> {name} -> {parsed}"
    );
}

fn assert_preferred(smiles: &str, expected: &str) {
    let name = smiles_to_name(smiles).unwrap();
    assert_eq!(name.text, expected, "{smiles}");
    assert_eq!(name.status, NameStatus::Preferred, "{smiles}");
    assert_round_trip(smiles, &name.text);
}

#[test]
fn pseudohalides_of_sulfonic_and_sulfinic_acids_are_named_like_the_halides() {
    for (smiles, expected) in [
        // P-66.5.1.3.2: CH3-SO2-CN methanesulfonyl cyanide (PIN) and
        // C6H5-SeO-CN benzeneseleninyl cyanide (PIN).
        ("CS(=O)(=O)C#N", "methanesulfonyl cyanide"),
        ("N#C[Se](=O)c1ccccc1", "benzeneseleninyl cyanide"),
        // P-65.5.2.1: acyl pseudohalides of suffix acids (carboxylic,
        // sulfonic, sulfinic, selenonic, etc.) are named like the halides
        // of P-65.5.1, the acyl group (P-65.3.2.2.2) and then the class.
        // Acyl cyanides are senior to nitriles (P-66.5.1.3.1).
        ("N#CS(=O)(=O)c1ccccc1", "benzenesulfonyl cyanide"),
        (
            "Cc1ccc(cc1)S(=O)(=O)C#N",
            "4-methylbenzene-1-sulfonyl cyanide",
        ),
        ("CS(=O)C#N", "methanesulfinyl cyanide"),
        ("C[Se](=O)(=O)C#N", "methaneselenonyl cyanide"),
        ("C[Se](=O)C#N", "methaneseleninyl cyanide"),
        ("C[Te](=O)(=O)C#N", "methanetelluronyl cyanide"),
        ("CS(=O)(=S)C#N", "methanesulfonothioyl cyanide"),
        ("CS(=O)(=O)N=[N+]=[N-]", "methanesulfonyl azide"),
        (
            "Cc1ccc(cc1)S(=O)(=O)N=[N+]=[N-]",
            "4-methylbenzene-1-sulfonyl azide",
        ),
        ("CS(=O)(=O)N=C=O", "methanesulfonyl isocyanate"),
        ("CS(=O)(=O)N=C=S", "methanesulfonyl isothiocyanate"),
        ("CS(=O)(=O)N=C=[Se]", "methanesulfonyl isoselenocyanate"),
        // The cyanide carbon is not a parent carbon; another nitrile is
        // cited as cyano (P-66.5.1.1.4).
        ("N#CCS(=O)(=O)C#N", "cyanomethanesulfonyl cyanide"),
        (
            "N#CS(=O)(=O)c1ccc(cc1)C#N",
            "4-cyanobenzene-1-sulfonyl cyanide",
        ),
        // P-65.5.1: ethane-1,2-disulfonyl dibromide multiplies the class
        // name; oxalyl diisothiocyanate does for a pseudohalide (P-65.5.2.1).
        (
            "ClS(=O)(=O)CCS(=O)(=O)Cl",
            "ethane-1,2-disulfonyl dichloride",
        ),
        (
            "BrS(=O)(=O)CCS(=O)(=O)Br",
            "ethane-1,2-disulfonyl dibromide",
        ),
        (
            "N#CS(=O)(=O)CCS(=O)(=O)C#N",
            "ethane-1,2-disulfonyl dicyanide",
        ),
        (
            "O=C=NS(=O)(=O)c1ccc(cc1)S(=O)(=O)N=C=O",
            "benzene-1,4-disulfonyl diisocyanate",
        ),
        // Carbon acyl groups are unchanged.
        ("CC(=O)C#N", "acetyl cyanide"),
        ("CS(=O)(=O)Cl", "methanesulfonyl chloride"),
    ] {
        assert_preferred(smiles, expected);
    }
}

#[test]
fn unlike_sulfonyl_halides_and_pseudohalides_are_not_multiplied_together() {
    for (smiles, expected) in [
        // P-65.5.2.1: halides are senior to pseudohalides and azide to
        // cyanide; cyanosulfonyl is the prefix of P-65.5.4. These read as
        // `disulfonyl chloride` before, a name for another structure. The
        // form with both class names, ethane-1,2-disulfonyl bromide chloride
        // (P-65.5.1, butanedioyl bromide chloride), is not generated, for
        // carbonyl halides either.
        (
            "N#CS(=O)(=O)c1ccc(cc1)S(=O)(=O)Cl",
            "4-(cyanosulfonyl)benzene-1-sulfonyl chloride",
        ),
        (
            "N#CS(=O)(=O)CCS(=O)(=O)N=[N+]=[N-]",
            "2-(cyanosulfonyl)ethane-1-sulfonyl azide",
        ),
        (
            "ClS(=O)(=O)CCS(=O)(=O)Br",
            "2-(bromosulfonyl)ethane-1-sulfonyl chloride",
        ),
        // Without another carbon on the sulfur these are unchanged: the
        // PIN, sulfuryl dicyanide (P-67.1.2.5.1), is not generated.
        ("N#CS(=O)(=O)C#N", "(cyanosulfonyl)methanenitrile"),
        // P-65.3.2.2.2: the acyl group of a thio acid keeps its infix,
        // which the prefix dropped before.
        ("N#CS(=O)(=S)CC(=O)O", "(cyanosulfonothioyl)acetic acid"),
        ("ClS(=O)(=S)CC(=O)O", "(chlorosulfonothioyl)acetic acid"),
        ("ClS(=S)CC(=O)O", "(chlorosulfinothioyl)acetic acid"),
        // A more senior class keeps the suffix (P-41).
        ("N#CS(=O)(=O)CC(=O)O", "(cyanosulfonyl)acetic acid"),
        (
            "N#CS(=O)(=O)c1ccc(cc1)C(=O)Cl",
            "4-(cyanosulfonyl)benzoyl chloride",
        ),
    ] {
        let name = smiles_to_name(smiles).unwrap().text;
        assert_eq!(name, expected, "{smiles}");
        assert_round_trip(smiles, &name);
    }
}

#[test]
fn a_group_that_only_resembles_a_pseudohalide_is_not_given_its_class_name() {
    // A protonated cyanide or isocyanate and a triazane zwitterion are not
    // the -CN, -NCO and -N3 of P-65.5.2.1.
    for smiles in [
        "CS(=O)(=O)C#[NH+]",
        "CS(=O)(=O)[NH+]=C=O",
        "CS(=O)(=O)N[NH2+][NH-]",
    ] {
        if let Ok(name) = smiles_to_name(smiles) {
            assert!(!name.text.starts_with("methanesulfonyl "), "{smiles}");
        }
    }
}

#[test]
fn the_cyanide_of_phosphoric_acid_uses_the_cyanido_infix() {
    for (smiles, expected) in [
        // P-67.1.2.1: phosphoric acid is modified by infixes; cyanido is
        // -CN (P-67.1.2.3.2). The acid keeps its name while an -OH is left
        // (P-67.1.2.4), as sulfurocyanidic acid (P-67.1.2.4.1.1) and
        // phosphorodicyanatidic acid (P-67.1.2.4.1.3) do.
        ("N#CP(=O)(O)O", "phosphorocyanidic acid"),
        // P-67.1.3.1: `ic acid` becomes `ate`; an acid salt keeps
        // `hydrogen`, as in sodium hydrogen borocyanatidate.
        ("N#CP(=O)([O-])[O-]", "phosphorocyanidate"),
        ("N#CP(=O)(O)[O-]", "hydrogen phosphorocyanidate"),
        (
            "[Na+].[Na+].N#CP(=O)([O-])[O-]",
            "disodium phosphorocyanidate",
        ),
        (
            "[Na+].N#CP(=O)(O)[O-]",
            "sodium hydrogen phosphorocyanidate",
        ),
        // P-16.3.6 (a), (e): bis(hydrogen sulfate), tris(dihydrogen phosphate).
        (
            "[Ca+2].N#CP(=O)(O)[O-].N#CP(=O)(O)[O-]",
            "calcium bis(hydrogen phosphorocyanidate)",
        ),
        // P-67.1.3.2: trimethyl phosphate, ethyl methyl phenyl phosphate,
        // methyl dihydrogen phosphate, sodium methyl hydrogen phosphate.
        ("N#CP(=O)(OC)OC", "dimethyl phosphorocyanidate"),
        ("N#CP(=O)(OC)OCC", "ethyl methyl phosphorocyanidate"),
        ("N#CP(=O)(O)OC", "methyl hydrogen phosphorocyanidate"),
        (
            "N#CP(=O)(O)O[C@H](C)CC",
            "(2R)-butan-2-yl hydrogen phosphorocyanidate",
        ),
        ("N#CP(=O)(OC)[O-]", "methyl phosphorocyanidate"),
        ("[Na+].N#CP(=O)(OC)[O-]", "sodium methyl phosphorocyanidate"),
        (
            "[Ca+2].N#CP(=O)(OC)[O-].N#CP(=O)(OC)[O-]",
            "calcium bis(methyl phosphorocyanidate)",
        ),
    ] {
        assert_preferred(smiles, expected);
    }
}

#[test]
fn other_cyanophosphorus_compounds_keep_their_names() {
    // Not generated here: phosphonic and phosphinic cyanides, halides,
    // amides, chalcogen analogues, anhydrides and a configured phosphorus
    // centre. Their names are unchanged and still read back.
    for (smiles, expected) in [
        ("CP(=O)(O)C#N", "[hydroxy(methyl)phosphoryl]methanenitrile"),
        ("N#CP(=O)(O)Cl", "[chloro(hydroxy)phosphoryl]methanenitrile"),
        ("N#CP(=O)(O)N", "[amino(hydroxy)phosphoryl]methanenitrile"),
        (
            "N#CP(=O)(O)OC(=O)C",
            "[(acetyloxy)(hydroxy)phosphoryl]methanenitrile",
        ),
        (
            "N#CP(=O)(O)OCC(=O)O",
            "{[cyano(hydroxy)phosphoryl]oxy}acetic acid",
        ),
        (
            "N#CP1(=O)OCCO1",
            "2-oxo-1,3,2lambda^5-dioxaphospholane-2-carbonitrile",
        ),
    ] {
        let name = smiles_to_name(smiles).unwrap().text;
        assert_eq!(name, expected, "{smiles}");
        assert_round_trip(smiles, &name);
    }
}

#[test]
fn a_multiplied_ester_anion_reads_back() {
    for (smiles, expected) in [
        // P-16.3.6 (a), (e): bis(hydrogen sulfate) and tris(dihydrogen
        // phosphate) multiply an anion of several words; the ester anions
        // of P-67.1.3.2 (sodium methyl hydrogen phosphate) are multiplied
        // the same way.
        (
            "[Ca+2].COS(=O)(=O)[O-].COS(=O)(=O)[O-]",
            "calcium bis(methyl sulfate)",
        ),
        (
            "[Ca+2].CCOS(=O)(=O)[O-].CCOS(=O)(=O)[O-]",
            "calcium bis(ethyl sulfate)",
        ),
        (
            "[Ca+2].COC(=O)[O-].COC(=O)[O-]",
            "calcium bis(methyl carbonate)",
        ),
        (
            "[Ca+2].COP(=O)(OC)[O-].COP(=O)(OC)[O-]",
            "calcium bis(dimethyl phosphate)",
        ),
        (
            "[Ca+2].C[C@H](CC)OS(=O)(=O)[O-].C[C@H](CC)OS(=O)(=O)[O-]",
            "calcium bis[(2R)-butan-2-yl sulfate]",
        ),
        (
            "[Ca+2].[Ca+2].COS(=O)(=O)[O-].COS(=O)(=O)[O-].CCOS(=O)(=O)[O-].CCOS(=O)(=O)[O-]",
            "dicalcium bis(ethyl sulfate) bis(methyl sulfate)",
        ),
    ] {
        assert_preferred(smiles, expected);
    }
    for (name, smiles) in [
        (
            "iron tris(methyl sulfate)",
            "[Fe+3].COS(=O)(=O)[O-].COS(=O)(=O)[O-].COS(=O)(=O)[O-]",
        ),
        (
            "calcium bis(methyl sulfate) dihydrate",
            "[Ca+2].COS(=O)(=O)[O-].COS(=O)(=O)[O-].O.O",
        ),
        // Each anion keeps its own ester groups and hydrogen.
        (
            "dicalcium bis(ethyl methyl phosphate) bis(methyl hydrogen phosphate)",
            "[Ca+2].[Ca+2].CCOP(=O)(OC)[O-].CCOP(=O)(OC)[O-].COP(=O)(O)[O-].COP(=O)(O)[O-]",
        ),
        (
            "heptasodium heptakis(methyl sulfate)",
            "[Na+].[Na+].[Na+].[Na+].[Na+].[Na+].[Na+].COS(=O)(=O)[O-].COS(=O)(=O)[O-].COS(=O)(=O)[O-].COS(=O)(=O)[O-].COS(=O)(=O)[O-].COS(=O)(=O)[O-].COS(=O)(=O)[O-]",
        ),
        // Anions with only `hydrogen` words are read as before, as salts
        // and as the acid esters of P-16.3.6 (a).
        (
            "calcium bis(hydrogen sulfate)",
            "[Ca+2].OS(=O)(=O)[O-].OS(=O)(=O)[O-]",
        ),
        (
            "calcium bis(dihydrogen phosphate)",
            "[Ca+2].OP(=O)(O)[O-].OP(=O)(O)[O-]",
        ),
        (
            "ethane-1,2-diyl bis(hydrogen sulfate)",
            "OS(=O)(=O)OCCOS(=O)(=O)O",
        ),
        (
            "ethane-1,2-diyl bis(Hydrogen sulfate)",
            "OS(=O)(=O)OCCOS(=O)(=O)O",
        ),
        // A name that was read before is read the same way: a spaced
        // substituent of the acid, and the acid part of an ester after a
        // substituent word.
        (
            "calcium bis(4-methyl benzoate)",
            "[Ca+2].Cc1ccc(cc1)C(=O)[O-].Cc1ccc(cc1)C(=O)[O-]",
        ),
        (
            "ethane-1,2-diyl bis(4-methyl benzoate)",
            "Cc1ccc(cc1)C(=O)OCCOC(=O)c1ccc(C)cc1",
        ),
        (
            "CALCIUM BIS(METHYL SULFATE)",
            "[Ca+2].COS(=O)(=O)[O-].COS(=O)(=O)[O-]",
        ),
        (
            "D-fructofuranose 1,6-bis(phosphate)",
            "O=P(O)(O)OC[C@H]1OC(O)(COP(=O)(O)O)[C@@H](O)[C@@H]1O",
        ),
        ("butanebis(thioic acid)", "OC(=S)CCC(O)=S"),
    ] {
        assert_round_trip(smiles, name);
    }
}
