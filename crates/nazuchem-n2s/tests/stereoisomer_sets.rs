// SPDX-License-Identifier: MIT
// Origin: original
// Source: ReShiki | crates/opsin/tests/stereoisomer_sets.rs | 82d110bdf311502cdb3b3b4fc74827d613c2768a

use std::collections::BTreeSet;

use nazuchem_kernel::smiles::stereoisomers;
use nazuchem_n2s::{ParseOptions, Parser, Status};

// These independent representatives were checked with RDKit CIP assignments
// and conformers before being used as expectations.
const BICYCLO_A: &str = "O[C@@H]1[C@H](C)[C@H]2CC[C@@H]1C2";
const BICYCLO_B: &str = "O[C@@H]1[C@@H](C)[C@H]2CC[C@@H]1C2";
const BICYCLO_C: &str = "[C@@H]12[C@@H](C[C@@H](CC1)C2)O";
const BICYCLO_D: &str = "C[C@@H]1[C@H]2CC[C@H](C2)[C@@H]1[C@@H](C)CO";
const BICYCLO_E: &str = "C[C@H](F)[C@@H]1[C@H]2CC[C@H](C2)[C@@H]1O";

fn check_relative_mirror(name: &str, actual: &BTreeSet<String>) -> Result<(), String> {
    // Leading rel-/rac- applies to the whole compound:
    // https://iupac.qmul.ac.uk/BlueBook/P9.html#930102
    // https://iupac.qmul.ac.uk/BlueBook/P9.html#930103
    let relative = name.split_once('-').is_some_and(|(prefix, _)| {
        prefix.eq_ignore_ascii_case("rel") || prefix.eq_ignore_ascii_case("rac")
    });
    if !relative {
        return Ok(());
    }
    let first = actual.first().unwrap();
    let mirror_pair = stereoisomers(&format!("{first} |r|"), &Default::default()).unwrap();
    if actual != &mirror_pair {
        return Err(format!(
            "leading rel-/rac- requires one mirror pair: actual {actual:?}, mirror pair {mirror_pair:?}"
        ));
    }
    Ok(())
}

/// Parsing must succeed; only a stereoisomer-set mismatch is an `Err`, so a
/// known failure cannot hide a name that stopped parsing.
fn check(name: &str, intended: &str) -> Result<(), String> {
    let result = Parser::new().unwrap().parse(name, &ParseOptions::strict());
    assert_eq!(result.status, Status::Success, "{name}: {}", result.message);
    let actual = result.structure.unwrap().semantic_cxsmiles().unwrap();
    let actual_set = stereoisomers(&actual, &Default::default()).unwrap();
    let intended_set = stereoisomers(intended, &Default::default()).unwrap();
    if actual_set != intended_set {
        return Err(format!(
            "{actual}: actual {actual_set:?}, intended {intended_set:?}"
        ));
    }
    check_relative_mirror(name, &actual_set)
}

#[test]
fn names_match_intended_stereoisomer_sets() {
    let a_pair = format!("{BICYCLO_A} |r|");
    let b_pair = format!("{BICYCLO_B} |r|");
    let c_pair = format!("{BICYCLO_C} |r|");
    let d_pair = format!("{BICYCLO_D} |r|");
    let e_pair = format!("{BICYCLO_E} |r|");

    // Cyclohexane representatives were independently checked with
    // reshiki --cli name --smiles '<plain representative>' --json.
    for (name, intended, known_failure) in [
        (
            "(1R,2S,3R,4S)-3-methylbicyclo[2.2.1]heptan-2-ol",
            BICYCLO_A,
            None::<u32>,
        ),
        (
            "rel-(1R,2S,3R,4S)-3-methylbicyclo[2.2.1]heptan-2-ol",
            a_pair.as_str(),
            None,
        ),
        (
            "rac-(1R,2S,3R,4S)-3-methylbicyclo[2.2.1]heptan-2-ol",
            a_pair.as_str(),
            None,
        ),
        (
            "rel-(1R,2S,3S,4S)-3-methylbicyclo[2.2.1]heptan-2-ol",
            b_pair.as_str(),
            None,
        ),
        (
            "rac-(1R,2S,3S,4S)-3-methylbicyclo[2.2.1]heptan-2-ol",
            b_pair.as_str(),
            None,
        ),
        (
            "rac-(1R,2R,4S)-bicyclo[2.2.1]heptan-2-ol",
            c_pair.as_str(),
            None,
        ),
        (
            "rel-(1R,2S)-2-methylcyclohexan-1-ol",
            // CLI: (1R,2S).
            "O[C@H]1[C@@H](C)CCCC1 |r|",
            None,
        ),
        (
            "rel-1-[(1R,2S)-2-methylcyclohexyl]ethan-1-one",
            // CLI: ring (1R,2S).
            "CC(=O)[C@H]1[C@@H](C)CCCC1 |r|",
            None,
        ),
        (
            "(1R*,2S*)-2-methylcyclohexan-1-ol",
            "O[C@H]1[C@@H](C)CCCC1 |r|",
            None,
        ),
        (
            "(1RS,2SR)-2-methylcyclohexan-1-ol",
            "O[C@H]1[C@@H](C)CCCC1 |r|",
            None,
        ),
        // A ring cis/trans descriptor keeps the rel/rac group of the centre
        // an earlier R*/RS descriptor configured. CLI: cis is (1R,2S),
        // trans (1R,2R), cis acid (1S,2R).
        (
            "(1R*)-cis-2-methylcyclohexan-1-ol",
            "O[C@H]1[C@@H](C)CCCC1 |r|",
            None,
        ),
        (
            "rel-(1R)-cis-2-methylcyclohexan-1-ol",
            "O[C@H]1[C@@H](C)CCCC1 |r|",
            None,
        ),
        (
            "rac-(1R)-cis-2-methylcyclohexan-1-ol",
            "O[C@H]1[C@@H](C)CCCC1 |r|",
            None,
        ),
        (
            "(1RS)-cis-2-methylcyclohexan-1-ol",
            "O[C@H]1[C@@H](C)CCCC1 |r|",
            None,
        ),
        (
            "(1R*)-trans-2-methylcyclohexan-1-ol",
            "O[C@H]1[C@H](C)CCCC1 |r|",
            None,
        ),
        (
            "rac-(1S)-cis-2-phenylcyclopropane-1-carboxylic acid",
            "OC(=O)[C@H]1C[C@H]1c1ccccc1 |r|",
            None,
        ),
        (
            "(1R)-cis-2-methylcyclohexan-1-ol",
            "O[C@H]1[C@@H](C)CCCC1",
            None,
        ),
        // An agreeing 5alpha keeps C5 in the rel group.
        (
            "rel-(5R,8R,9S,10S,13R,14S,17R,20R)-5alpha-cholestane",
            "CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC[C@H]4CCCC[C@]4(C)[C@H]3CC[C@]12C |r|",
            None,
        ),
        // The anomer joins the DL group: alpha-D plus alpha-L, not beta-L.
        (
            "alpha-DL-glucopyranose",
            "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O |r|",
            None,
        ),
        (
            "beta-DL-glucopyranose",
            "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O |r|",
            None,
        ),
        // #42: rel/rac at the front of the name give the configuration of the
        // entire molecule (P-93.1.2.1, P-93.1.3), so they reach nested
        // descriptor sets, other words of an ester and stereoparent centres.
        (
            "rel-(2R)-2-[(1R,2S,3R,4S)-3-methylbicyclo[2.2.1]heptan-2-yl]propan-1-ol",
            d_pair.as_str(),
            None,
        ),
        (
            "rac-(2R)-2-[(1R,2S,3R,4S)-3-methylbicyclo[2.2.1]heptan-2-yl]propan-1-ol",
            d_pair.as_str(),
            None,
        ),
        (
            "rel-(1R,2S,3R,4S)-3-[(1S)-1-fluoroethyl]bicyclo[2.2.1]heptan-2-ol",
            e_pair.as_str(),
            None,
        ),
        (
            "rac-(1R,2S,3R,4S)-3-[(1S)-1-fluoroethyl]bicyclo[2.2.1]heptan-2-ol",
            e_pair.as_str(),
            None,
        ),
        (
            "rel-(1R)-1-[(1R,2S)-2-methylcyclohexyl]ethan-1-ol",
            // CLI: side chain (1R), ring (1R,2S).
            "C[C@@H](O)[C@H]1[C@@H](C)CCCC1 |r|",
            None,
        ),
        (
            "rac-(1R)-1-[(1R,2S)-2-methylcyclohexyl]ethan-1-ol",
            "C[C@@H](O)[C@H]1[C@@H](C)CCCC1 |r|",
            None,
        ),
        (
            "rac-(1R,2S)-2-[(1S)-1-fluoroethyl]cyclohexan-1-ol",
            // CLI: ring (1R,2S), fluoroethyl (1S).
            "O[C@H]1[C@@H]([C@@H](F)C)CCCC1 |r|",
            None,
        ),
        (
            "rel-(1R)-1-phenylethyl (2R)-2-fluoropropanoate",
            "F[C@@H](C(=O)O[C@H](C)c1ccccc1)C |r|",
            None,
        ),
        // The configurations implied by a stereoparent are relative too
        // (P-101.8.2, P-101.8.3); the unspecified C3 stays unspecified.
        (
            "rel-(5alpha)-androstane",
            "C[C@@]12CCC[C@H]1[C@@H]1CC[C@H]3CCCC[C@]3(C)[C@H]1CC2 |r|",
            None,
        ),
        (
            "rac-(5beta)-androstane",
            "C[C@@]12CCC[C@H]1[C@@H]1CC[C@@H]3CCCC[C@]3(C)[C@H]1CC2 |r|",
            None,
        ),
        (
            "rel-(5alpha)-androstan-3-ol",
            "C[C@@]12CCC[C@H]1[C@@H]1CC[C@H]3CC(CC[C@]3(C)[C@H]1CC2)O |r|",
            None,
        ),
        (
            "rel-5alpha-androstan-3-ol",
            "C[C@@]12CCC[C@H]1[C@@H]1CC[C@H]3CC(CC[C@]3(C)[C@H]1CC2)O |r|",
            None,
        ),
        (
            "rac-(3alpha,5alpha)-cholestan-3-yl acetate",
            "C(C)(=O)O[C@H]1C[C@@H]2CC[C@H]3[C@@H]4CC[C@H]([C@@H](CCCC(C)C)C)[C@]4(CC[C@@H]3[C@]2(CC1)C)C |r|",
            None,
        ),
        (
            "rel-(2R)-2-[(1R)-1-fluoroethyl]-5alpha-androstane",
            "F[C@H](C)[C@@H]1CC[C@@H]2CC[C@H]3[C@@H]4CCC[C@@]4(C)CC[C@@H]3[C@]2(C1)C |r|",
            None,
        ),
    ] {
        match (known_failure, check(name, intended)) {
            (None, Ok(())) => eprintln!("PASS {name}"),
            (Some(issue), Err(error)) => eprintln!("KNOWN #{issue} {name}: {error}"),
            (None, Err(error)) => panic!("{name}: {error}"),
            (Some(issue), Ok(())) => {
                panic!("{name} now passes: move this row to the passing cases and close #{issue}")
            }
        }
    }
}

#[test]
fn leading_rel_rac_stops_at_the_molecular_entity() {
    // #42: rel/rac are cited at the front of the name and describe one whole
    // molecular entity (P-93.1.2.1, P-93.1.3). Other ions keep their own
    // descriptors, and a qualifier that does not lead the name keeps the scope
    // of its own descriptor set.
    for (name, expected) in [
        (
            "rel-(1R)-1-[(1R,2S)-2-methylcyclohexyl]ethan-1-ol",
            "C[C@@H]1[C@@H](CCCC1)[C@@H](C)O |o1:1,2,7|",
        ),
        (
            "rel-(1R)-1-[(1R,2S)-2-methylcyclohexyl]ethan-1-ol hydrochloride",
            "Cl.C[C@@H]1[C@@H](CCCC1)[C@@H](C)O |o1:2,3,8|",
        ),
        (
            "rel-(2R)-butan-2-aminium (2R)-2-hydroxypropanoate",
            "O[C@@H](C(=O)[O-])C.C[C@H](CC)[NH3+] |o1:7|",
        ),
        // The anion is its own compound, so its name starts at the qualifier.
        (
            "sodium rel-(2R)-2-[(1R,2S)-2-methylcyclohexyl]propanoate",
            "C[C@@H]1[C@@H](CCCC1)[C@H](C(=O)[O-])C.[Na+] |o1:1,2,7|",
        ),
        // Each ion keeps the qualifier at the front of its own name.
        (
            "rel-(2R)-2-[(1R)-1-fluoroethyl]butan-2-aminium rac-(2R)-2-[(1R)-1-fluoroethyl]propanoate",
            "F[C@H](C)[C@@H](C(=O)[O-])C.F[C@H](C)[C@@](C)(CC)[NH3+] |&1:1,3,o1:9,11|",
        ),
        (
            "rac-(2R)-2-[(1R)-1-fluoroethyl]butan-2-aminium rel-(2R)-2-[(1R)-1-fluoroethyl]propanoate",
            "F[C@H](C)[C@@H](C(=O)[O-])C.F[C@H](C)[C@@](C)(CC)[NH3+] |o1:1,3,&1:9,11|",
        ),
        // A standalone rac- still turns an unanchored cage's relative group
        // into a racemic one.
        (
            "rac-1-[bicyclo[2.2.1]heptan-2-exo-yl]ethan-1-one",
            "[C@@H]12[C@@H](C[C@@H](CC1)C2)C(C)=O |r|",
        ),
        (
            "(1R)-1-phenylethyl rac-(2R)-2-fluoropropanoate",
            "F[C@@H](C(=O)O[C@H](C)C1=CC=CC=C1)C |&1:1|",
        ),
        (
            "N-[(1R)-1-phenylethyl]-(rel)-(2R)-2-fluoropropanamide",
            "C1(=CC=CC=C1)[C@@H](C)NC([C@@H](C)F)=O |o1:10|",
        ),
        // A standalone qualifier inside the first bracket is not at the front.
        (
            "N-[(rel)-(1R)-1-phenylethyl]-(2R)-2-fluoropropanamide",
            "C1(=CC=CC=C1)[C@@H](C)NC([C@@H](C)F)=O |o1:6|",
        ),
        (
            "N-[(rac)-(1R)-1-phenylethyl]-(2R)-2-fluoropropanamide",
            "C1(=CC=CC=C1)[C@@H](C)NC([C@@H](C)F)=O |&1:6|",
        ),
        // An unspecified centre stays unspecified beside configured ones.
        (
            "rel-(2R,3xi)-3,4-dichloropentan-2-ol",
            "ClC([C@@H](C)O)C(C)Cl |o1:2|",
        ),
        (
            "rel-1-[(1R,2S)-2-methylcyclohexyl]ethan-1-ol",
            "C[C@@H]1[C@@H](CCCC1)C(C)O |o1:1,2|",
        ),
        // An unbracketed descriptor is covered like a bracketed one.
        ("rac-R-butan-2-ol", "C[C@H](CC)O |r|"),
        ("rel-R-butan-2-ol", "C[C@H](CC)O |o1:1|"),
        // With no configured centre, rac- still describes the single centre.
        ("rac-butan-2-ol", "C[C@H](CC)O |r|"),
        (
            "rac-N-(butan-2-yl)-2-methylalanine",
            "C[C@H](CC)NC(C)(C(=O)O)C |r|",
        ),
    ] {
        let result = Parser::new().unwrap().parse(name, &ParseOptions::strict());
        assert_eq!(result.status, Status::Success, "{name}: {}", result.message);
        let actual = result.structure.unwrap().semantic_cxsmiles().unwrap();
        assert_eq!(actual, expected, "{name}");
    }
    // A stereo group cannot hold the axis of the same molecule.
    for name in [
        "rel-(1R)-1-fluoroethyl (1P)-2',5'-dimethoxy-6-nitro[1,1'-biphenyl]-2-carboxylate",
        "rac-(1R)-1-fluoroethyl (1P)-2',5'-dimethoxy-6-nitro[1,1'-biphenyl]-2-carboxylate",
    ] {
        let result = Parser::new().unwrap().parse(name, &ParseOptions::strict());
        assert_eq!(result.status, Status::Failure, "{name}");
        assert!(
            result.message.contains("axial"),
            "{name}: {}",
            result.message
        );
        let mut options = ParseOptions::strict();
        options.warn_rather_than_fail_on_uninterpretable_stereochemistry = true;
        let result = Parser::new().unwrap().parse(name, &options);
        assert_eq!(result.status, Status::Warning, "{name}");
    }
    // rac/rel need at least one defined centre: a leading qualifier over
    // only ξ centres is refused (P-93.1.3; #14, #36).
    for name in ["(rac)-N-acetyl-(ξ)-alanine", "(rel)-N-acetyl-(ξ)-alanine"] {
        let result = Parser::new().unwrap().parse(name, &ParseOptions::strict());
        assert_eq!(result.status, Status::Failure, "{name}: {}", result.message);
        assert!(
            result.message.contains("defined stereocentre"),
            "{name}: {}",
            result.message
        );
    }
    // Warning mode reports the conflict once, also when ξ clears a group
    // that a descriptor set had already applied.
    let mut options = ParseOptions::strict();
    options.warn_rather_than_fail_on_uninterpretable_stereochemistry = true;
    for name in ["(rac)-N-acetyl-(ξ)-alanine", "rac-(2R,2xi)-butan-2-ol"] {
        let result = Parser::new().unwrap().parse(name, &options);
        assert_eq!(result.status, Status::Warning, "{name}: {}", result.message);
        let conflicts = result
            .warnings
            .iter()
            .filter(|warning| warning.message.contains("defined stereocentre"))
            .count();
        assert_eq!(conflicts, 1, "{name}: {:?}", result.warnings);
    }
}

#[test]
fn literal_cxsmiles_audit_rejects_split_or_lost_mirror_pairs() {
    let a_pair = format!("{BICYCLO_A} |r|");
    let c_pair = format!("{BICYCLO_C} |r|");
    let d_pair = format!("{BICYCLO_D} |r|");
    let e_pair = format!("{BICYCLO_E} |r|");
    for (name, literal, intended, size) in [
        // #898 round 3: one mirror pair was split into independent groups.
        (
            "rel-(1R,3R)-3-methylbicyclo[2.2.1]heptan-2-exo-ol",
            "C[C@H]1[C@H]([C@@H]2CC[C@H]1C2)O |o1:1,o2:2,3,6|",
            a_pair.as_str(),
            4,
        ),
        // #898 round 1: racemic bridgehead configurations were lost.
        (
            "rac-(2-exo)-bicyclo[2.2.1]heptan-2-ol",
            BICYCLO_C,
            c_pair.as_str(),
            1,
        ),
        (
            "rac-(1R,2-exo)-bicyclo[2.2.1]heptan-2-ol",
            "[C@@H]12[C@@H](C[C@@H](CC1)C2)O |&1:0|",
            c_pair.as_str(),
            2,
        ),
        // #898 rounds 4/5: nested descriptor sets remained absolute.
        (
            "rel-(2R)-2-[(1R,2S,3R,4S)-3-methylbicyclo[2.2.1]heptan-2-yl]propan-1-ol",
            "C[C@H]1[C@H]([C@@H]2CC[C@H]1C2)[C@H](CO)C |o1:8|",
            d_pair.as_str(),
            2,
        ),
        (
            "rel-(1R,2S,3R,4S)-3-[(1S)-1-fluoroethyl]bicyclo[2.2.1]heptan-2-ol",
            "F[C@@H](C)[C@H]1[C@H]([C@@H]2CC[C@H]1C2)O |o1:3,4,5,8|",
            e_pair.as_str(),
            2,
        ),
    ] {
        let actual = stereoisomers(literal, &Default::default()).unwrap();
        assert_eq!(actual.len(), size, "{literal}");
        assert_ne!(
            actual,
            stereoisomers(intended, &Default::default()).unwrap(),
            "{literal}"
        );
        assert!(check_relative_mirror(name, &actual).is_err(), "{literal}");
    }
}
