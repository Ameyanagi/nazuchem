// SPDX-License-Identifier: MIT
// Origin: original
// Source: Nazuchem | crates/nazuchem-n2s/tests/ir_oracle.rs | Phase B step 8a

//! The IR adapter is an additional gate over all frozen name-parse records.
//!
//! Exact legacy semantic CXSMILES is preferred. Different ring numbers, atom
//! traversals or equivalent slash directions are accepted only after the
//! kernel's canonical isomeric identity agrees. Canonical identity does not
//! include labels, enhanced groups, polymers or position variation: these are
//! compared separately after mapping writer indices to stable IR identities.
//! Locants are checked in that same identity order. Unsupported canonical
//! inputs (for example wildcards) still require exact molecular SMILES.

use std::io::{BufRead, BufReader};

use flate2::read::GzDecoder;
use nazuchem_ir::{CxsmilesError, Molecule};
use nazuchem_kernel::{budget::Budget, identity::canonical_identity, smiles};
use nazuchem_n2s::{ParseOptions, ParseResult, Parser, Status, Structure};
use serde_json::{Value, json};

#[path = "support/ir_cxsmiles.rs"]
mod cx;
#[path = "support/ir_schema.rs"]
mod schema;
#[path = "support/ir_stereo.rs"]
mod stereo;

#[derive(Debug, Default, PartialEq, Eq)]
struct Counts {
    records: usize,
    success: usize,
    warning: usize,
    failure: usize,
    supported: usize,
    axial: usize,
    exact: usize,
    canonical: usize,
    not_expressible: usize,
    direct_tetrahedral: usize,
    direct_double_bond: usize,
    direct_atropisomer: usize,
}

impl Counts {
    fn record(&mut self, result: &ParseResult) {
        self.records += 1;
        match result.status {
            Status::Success => self.success += 1,
            Status::Warning => self.warning += 1,
            Status::Failure => self.failure += 1,
        }
    }

    fn assert_baseline(&self, records: usize, success: usize, warning: usize, failure: usize) {
        assert_eq!(
            (self.records, self.success, self.warning, self.failure),
            (records, success, warning, failure),
            "parse coverage/status counts changed: {self:?}"
        );
        assert_eq!(
            self.supported,
            success + warning,
            "every parsed structure must be checked: {self:?}"
        );
        assert_eq!(
            self.exact + self.canonical + self.not_expressible,
            self.supported,
            "every supported IR must reproduce molecular identity: {self:?}"
        );
    }
}

fn published_schema() -> schema::Validator {
    schema::Validator::new(
        serde_json::from_str(nazuchem_ir::JSON_SCHEMA).expect("published schema is JSON"),
    )
    .expect("published schema uses only the independently validated vocabulary")
}

fn status(result: &ParseResult) -> &'static str {
    match result.status {
        Status::Success => "SUCCESS",
        Status::Warning => "WARNING",
        Status::Failure => "FAILURE",
    }
}

fn check_parse(result: &ParseResult, expected: &Value) -> Result<(), String> {
    let warnings: Vec<_> = result
        .warnings
        .iter()
        .map(|warning| {
            json!({
                "kind": warning.kind.as_str(), "message": warning.message,
            })
        })
        .collect();
    let actual = json!({"status": status(result), "message": result.message, "warnings": warnings});
    let frozen = json!({"status": expected["status"], "message": expected["message"], "warnings": expected["warnings"]});
    if actual != frozen {
        return Err(format!(
            "parse result changed: expected {frozen}, actual {actual}"
        ));
    }
    if result.structure.is_some() != (expected["status"] != "FAILURE") {
        return Err("parse structure presence changed".into());
    }
    Ok(())
}

fn canonical(text: &str) -> Result<String, String> {
    let budget = Budget::default();
    let raw = smiles::parse(cx::bare(text)).map_err(|error| error.to_string())?;
    let raw_metadata = nazuchem_kernel::stereo::bond_stereo_from_directions(
        &raw.graph,
        &raw.metadata,
        &raw.directions,
    )?;
    let prepared = smiles::prepare(cx::bare(text), &budget).map_err(|error| error.to_string())?;
    // Use the kernel's own H-removal identity map, rather than assuming atom
    // counts or canonical equivalence can prove preservation of every unit.
    let removed = nazuchem_kernel::hydrogens::remove(&nazuchem_kernel::hydrogens::Input {
        unknown_atoms: vec![false; raw.graph.atoms.len()],
        graph: raw.graph,
        metadata: raw.metadata,
        directions: raw.directions,
        annotations: Default::default(),
    })
    .map_err(|error| error.to_string())?;
    for (id, atom) in raw_metadata.atoms.iter().enumerate() {
        if atom.chiral_tag != 0 {
            let kept = removed
                .kept_atoms
                .iter()
                .position(|&old| old == id)
                .and_then(|new| prepared.state.metadata.atoms.get(new));
            if kept.is_none_or(|atom| atom.chiral_tag == 0) {
                return Err(format!(
                    "kernel preparation dropped specified atom stereo {id}"
                ));
            }
        }
    }
    for (id, bond) in raw_metadata.bonds.iter().enumerate() {
        if matches!(bond.stereo, 2..=7) {
            let kept = removed
                .kept_bonds
                .iter()
                .position(|&old| old == id)
                .and_then(|new| prepared.state.metadata.bonds.get(new));
            if kept.is_none_or(|bond| !matches!(bond.stereo, 2..=7)) {
                return Err(format!(
                    "kernel preparation dropped specified bond stereo {id}"
                ));
            }
        }
    }
    let counts = |metadata: &nazuchem_kernel::ranking::Metadata| {
        (
            metadata.atoms.iter().filter(|a| a.chiral_tag != 0).count(),
            metadata
                .bonds
                .iter()
                .filter(|b| matches!(b.stereo, 2..=5))
                .count(),
            metadata
                .bonds
                .iter()
                .filter(|b| matches!(b.stereo, 6 | 7))
                .count(),
        )
    };
    if counts(&raw_metadata) != counts(&prepared.state.metadata) {
        return Err("kernel preparation dropped a specified stereo unit".into());
    }
    canonical_identity(&prepared.state, &budget).map_err(|error| format!("{error:?}"))
}

fn check_structure(
    structure: &Structure,
    frozen_cxsmiles: Option<&str>,
    frozen_locants: Option<&Value>,
    validator: &schema::Validator,
    counts: &mut Counts,
) -> Result<(), String> {
    let (legacy, legacy_locants, axes) = structure
        .semantic_cxsmiles_with_locants_and_axes()
        .map_err(|error| error.to_string())?;
    if let Some(expected) = frozen_cxsmiles
        && legacy != expected
    {
        return Err(format!(
            "legacy semantic output changed: expected {expected}, actual {legacy}"
        ));
    }
    if let Some(expected) = frozen_locants
        && json!(legacy_locants) != *expected
    {
        return Err(format!(
            "legacy locants changed: expected {expected}, actual {legacy_locants:?}"
        ));
    }
    let (molecule, source_order) = structure
        .to_ir_with_atom_mapping()
        .map_err(|error| format!("IR projection: {error}"))?;
    let direct = stereo::source(structure, &molecule, &source_order)?;
    counts.direct_tetrahedral += direct.tetrahedral;
    counts.direct_double_bond += direct.double_bond;
    counts.direct_atropisomer += direct.atropisomer;
    counts.axial += usize::from(!axes.is_empty());
    molecule
        .validate()
        .map_err(|error| format!("IR validation: {error}"))?;
    if molecule.atoms.len() < legacy_locants.len() {
        return Err("IR atom identities differ from legacy written atoms".into());
    }
    cx::legacy_radicals(&molecule, &legacy)?;
    if molecule.displayed_locants()[..legacy_locants.len()] != legacy_locants {
        return Err(format!(
            "IR displayed locants differ: expected {legacy_locants:?}, actual {:?}",
            molecule.displayed_locants()
        ));
    }
    let encoded = molecule
        .to_json()
        .map_err(|error| format!("IR JSON write: {error}"))?;
    let document: Value = serde_json::from_str(&encoded).map_err(|error| error.to_string())?;
    validator
        .validate(&document)
        .map_err(|error| format!("published schema: {error}"))?;
    let restored =
        Molecule::from_json(&encoded).map_err(|error| format!("IR JSON read: {error}"))?;
    if restored != molecule {
        return Err("JSON round trip changed the complete IR".into());
    }
    let projection = match molecule.to_cxsmiles() {
        Ok(text) => text,
        Err(CxsmilesError::StereoNotExpressible) => {
            counts.not_expressible += 1;
            counts.supported += 1;
            return Ok(());
        }
        Err(error) => return Err(format!("IR convenience writer: {error}")),
    };
    let (actual, semantic_order) = molecule
        .semantic_cxsmiles_with_atom_order_for_tests()
        .map_err(|error| format!("IR semantic writer: {error}"))?;
    if restored
        .semantic_cxsmiles_for_tests()
        .map_err(|error| error.to_string())?
        != actual
        || restored.displayed_locants() != molecule.displayed_locants()
        || restored.to_cxsmiles().map_err(|error| error.to_string())?
            != molecule.to_cxsmiles().map_err(|error| error.to_string())?
    {
        return Err("JSON round trip changed CXSMILES or locant projection".into());
    }
    let order = molecule
        .smiles_atom_order()
        .map_err(|error| error.to_string())?;
    let identity: Vec<u32> = (0..molecule.atoms.len() as u32).collect();
    let mut sorted_order = order.clone();
    sorted_order.sort_unstable();
    if sorted_order != identity {
        return Err(format!(
            "writer traversal is not an identity permutation: {order:?}"
        ));
    }
    cx::written_graph(&molecule, &projection, &order)?;
    stereo::written(&molecule, &projection, &order)?;
    let expected_layers = cx::semantic_layers(&legacy, &identity[..legacy_locants.len()])?;
    let actual_layers = cx::semantic_layers(&actual, &semantic_order)?;
    if expected_layers != actual_layers {
        return Err(format!(
            "standard semantic CX layers changed: expected {expected_layers:?}, actual {actual_layers:?}\nlegacy {legacy}\nIR {actual}"
        ));
    }
    if legacy == actual {
        counts.exact += 1;
    } else {
        let expected = canonical(&legacy).map_err(|error| format!("legacy canonical identity unavailable ({error}); exact output required: {legacy} != {actual}"))?;
        let actual_identity = canonical(&actual)
            .map_err(|error| format!("IR canonical identity: {error}: {actual}"))?;
        if expected != actual_identity {
            return Err(format!(
                "canonical identity changed: {expected} != {actual_identity}\nlegacy {legacy}\nIR {actual}"
            ));
        }
        counts.canonical += 1;
    }
    counts.supported += 1;
    Ok(())
}

fn finish(corpus: &str, counts: &Counts, differences: Vec<String>) {
    eprintln!("{corpus}: {counts:?}");
    assert!(
        differences.is_empty(),
        "{corpus} IR gate failed; first differences:\n{}",
        differences.join("\n")
    );
}

#[test]
fn every_structure_oracle_record_passes_through_ir() {
    let parser = Parser::new().expect("embedded OPSIN resources");
    let validator = published_schema();
    let reader = BufReader::new(GzDecoder::new(
        include_bytes!("fixtures/structure-golden.jsonl.gz").as_slice(),
    ));
    let (mut counts, mut differences) = (Counts::default(), Vec::new());
    for line in reader.lines() {
        let frozen: Value =
            serde_json::from_str(&line.expect("fixture decompression")).expect("fixture JSON");
        let name = frozen["input"].as_str().unwrap();
        let result = parser.parse(
            name,
            &ParseOptions {
                allow_radicals: frozen["allow_radicals"].as_bool().unwrap(),
                ..ParseOptions::strict()
            },
        );
        counts.record(&result);
        let checked = check_parse(&result, &frozen).and_then(|()| match &result.structure {
            Some(structure) => check_structure(
                structure,
                frozen["cxsmiles"].as_str(),
                None,
                &validator,
                &mut counts,
            ),
            None => Ok(()),
        });
        if let Err(error) = checked
            && differences.len() < 30
        {
            differences.push(format!("{}: {name}: {error}", frozen["family"]));
        }
    }
    finish("structure oracle", &counts, differences);
    counts.assert_baseline(2624, 2535, 70, 19);
    assert_eq!(
        counts.axial, 0,
        "frozen structure oracle contains no axial structures"
    );
}

#[test]
fn every_configuration_oracle_record_passes_through_ir() {
    let parser = Parser::new().expect("embedded OPSIN resources");
    let validator = published_schema();
    let reader = BufReader::new(GzDecoder::new(
        include_bytes!("fixtures/configuration-golden.jsonl.gz").as_slice(),
    ));
    let (mut counts, mut differences) = (Counts::default(), Vec::new());
    for line in reader.lines() {
        let frozen: Value =
            serde_json::from_str(&line.expect("fixture decompression")).expect("fixture JSON");
        let name = frozen["input"].as_str().unwrap();
        let options = &frozen["options"];
        let result = parser.parse(
            name,
            &ParseOptions {
                allow_radicals: options["allow_radicals"].as_bool().unwrap(),
                output_radicals_as_wildcard_atoms: options["output_radicals_as_wildcard_atoms"]
                    .as_bool()
                    .unwrap(),
                detailed_failure_analysis: options["detailed_failure_analysis"].as_bool().unwrap(),
                interpret_acids_without_the_word_acid:
                    options["interpret_acids_without_the_word_acid"]
                        .as_bool()
                        .unwrap(),
                warn_rather_than_fail_on_uninterpretable_stereochemistry:
                    options["warn_rather_than_fail_on_uninterpretable_stereochemistry"]
                        .as_bool()
                        .unwrap(),
            },
        );
        counts.record(&result);
        let checked = check_parse(&result, &frozen).and_then(|()| match &result.structure {
            Some(structure) => check_structure(
                structure,
                frozen["cxsmiles"].as_str(),
                None,
                &validator,
                &mut counts,
            ),
            None => Ok(()),
        });
        if let Err(error) = checked
            && differences.len() < 30
        {
            differences.push(format!("{name} mask={}: {error}", frozen["mask"]));
        }
    }
    finish("configuration oracle", &counts, differences);
    counts.assert_baseline(640, 304, 112, 224);
    assert_eq!(
        counts.axial, 0,
        "frozen configuration oracle contains no axial structures"
    );
}

#[test]
fn every_frozen_generated_name_parse_passes_through_ir() {
    let parser = Parser::new().expect("embedded OPSIN resources");
    let validator = published_schema();
    let (mut counts, mut differences, mut records) = (Counts::default(), Vec::new(), 0);
    for line in include_str!("../../../tests/fixtures/golden/s2n_n2s.jsonl").lines() {
        let frozen: Value = serde_json::from_str(line).expect("golden JSON");
        records += 1;
        for expected in frozen["n2s"].as_array().unwrap() {
            let name = expected["name"].as_str().unwrap();
            let result = parser.parse(name, &ParseOptions::strict());
            counts.record(&result);
            let checked = check_parse(&result, expected).and_then(|()| match &result.structure {
                Some(structure) => {
                    let expected = expected["structure"].get("Ok").ok_or("frozen successful parse has no serialized structure")?;
                    let (_, _, axes) = structure.semantic_cxsmiles_with_locants_and_axes().map_err(|error| error.to_string())?;
                    let actual_axes: Vec<_> = axes.iter().map(|axis| json!({"atoms": axis.atoms, "refs": axis.refs, "clockwise": axis.clockwise, "descriptor": axis.descriptor})).collect();
                    if json!(actual_axes) != expected["axes"] { return Err("legacy axes changed".into()); }
                    check_structure(structure, expected["cxsmiles"].as_str(), Some(&expected["locants"]), &validator, &mut counts)
                }
                None => Ok(()),
            });
            if let Err(error) = checked
                && differences.len() < 30
            {
                differences.push(format!("{name}: {error}"));
            }
        }
    }
    finish("generated-name golden", &counts, differences);
    assert_eq!(records, 1650);
    counts.assert_baseline(1651, 1610, 0, 41);
    assert_eq!(
        counts.axial, 0,
        "frozen generated names contain no axial structures"
    );
}

#[test]
fn canonical_fallback_is_isomeric_and_preserves_element_isotope_charge_and_radicals() {
    for (left, right) in [
        ("CCO", "OCC"),
        ("N[C@@H](C)C(=O)O", "C[C@H](N)C(=O)O"),
        ("F/C=C/F", r"F\C=C\F"),
    ] {
        assert_eq!(
            canonical(left).unwrap(),
            canonical(right).unwrap(),
            "{left} / {right}"
        );
    }
    for (left, right) in [
        ("CCO", "COC"),
        ("[13CH4]", "C"),
        ("[NH4+]", "N"),
        ("[CH3]", "C"),
        ("N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O"),
        ("F/C=C/F", r"F/C=C\F"),
    ] {
        assert_ne!(
            canonical(left).unwrap(),
            canonical(right).unwrap(),
            "{left} / {right}"
        );
    }
}

#[test]
fn semantic_cx_layers_are_checked_under_writer_atom_index_mapping() {
    let legacy = "*CCC |$R&#32;one;;;tail$,o1:1,2,m:0:1.2,Sg:n:1,2,3::ht|";
    let reordered = "CCC* |$tail;;;R&#32;one$,o7:2,1,m:3:2.1,Sg:n:2,1,0::ht|";
    let expected = cx::semantic_layers(legacy, &[0, 1, 2, 3]).unwrap();
    assert_eq!(
        expected,
        cx::semantic_layers(reordered, &[3, 2, 1, 0]).unwrap()
    );
    for changed in [
        legacy.replace("R&#32;one", "other"),
        legacy.replace("o1:", "&1:"),
        legacy.replace("m:0:1.2", "m:0:1.3"),
        legacy.replace("::ht", "::hh"),
        legacy.replace("Sg:n:1,2,3", "Sg:n:1,2"),
    ] {
        assert_ne!(
            expected,
            cx::semantic_layers(&changed, &[0, 1, 2, 3]).unwrap(),
            "{changed}"
        );
    }
    assert!(cx::semantic_layers("CC |future:1|", &[0, 1]).is_err());
}

#[test]
fn published_schema_rejects_invalid_json_independently_of_deserialization() {
    let validator = published_schema();
    let structure = Parser::new()
        .unwrap()
        .parse("ethanol", &ParseOptions::strict())
        .structure
        .unwrap();
    let document: Value =
        serde_json::from_str(&structure.to_ir().unwrap().to_json().unwrap()).unwrap();
    validator.validate(&document).unwrap();
    let mut mutations = Vec::new();
    let mut changed = document.clone();
    changed.as_object_mut().unwrap().remove("nazuchem_ir");
    mutations.push(changed);
    let mut changed = document.clone();
    changed["nazuchem_ir"] = json!(2);
    mutations.push(changed);
    let mut changed = document.clone();
    changed["future"] = json!(true);
    mutations.push(changed);
    let mut changed = document.clone();
    changed["atoms"][0]["element"] = json!(119);
    mutations.push(changed);
    let mut changed = document.clone();
    changed["atoms"][0]["element"] = json!("carbon");
    mutations.push(changed);
    let mut changed = document.clone();
    changed["atoms"][0]["charge"] = json!(-129);
    mutations.push(changed);
    let mut changed = document.clone();
    changed["atoms"][0]["future"] = json!(true);
    mutations.push(changed);
    let mut changed = document.clone();
    changed["bonds"][0]["atoms"] = json!([0]);
    mutations.push(changed);
    let mut changed = document.clone();
    changed["bonds"][0]["order"] = json!("aromatic");
    mutations.push(changed);
    let mut changed = document;
    changed["stereo"] = json!([{"kind":"future"}]);
    mutations.push(changed);
    for mutation in mutations {
        assert!(
            validator.validate(&mutation).is_err(),
            "invalid document accepted: {mutation}"
        );
    }
    assert!(schema::Validator::new(json!({"$defs":{"unused":{"unsupported":true}}})).is_err());
    let recursive = schema::Validator::new(json!({"$defs":{"node":{"type":"object","properties":{"next":{"$ref":"#/$defs/node"}},"additionalProperties":false}},"$ref":"#/$defs/node"})).unwrap();
    recursive.validate(&json!({"next":{"next":{}}})).unwrap();
    assert!(
        recursive
            .validate(&json!({"next":{"invalid":true}}))
            .is_err()
    );
}

#[test]
fn canonical_fallback_fails_closed_when_preparation_erases_stereo() {
    for text in ["[N@](C)(CC)CCC", "[N@@](C)(CC)CCC", "F/C(F)=C/F"] {
        assert!(
            canonical(text).is_err(),
            "cleanup must not conceal raw assertions: {text}"
        );
    }
}

#[test]
fn unspecified_middle_double_bond_cannot_be_exported_with_extra_stereo() {
    let structure = Parser::new()
        .unwrap()
        .parse("(2E,6E)-octa-2,4,6-triene", &ParseOptions::strict())
        .structure
        .unwrap();
    let molecule = structure.to_ir().unwrap();
    assert_eq!(molecule.stereo.len(), 2);
    assert!(matches!(
        molecule.to_cxsmiles(),
        Err(CxsmilesError::StereoNotExpressible)
    ));
}

fn cumulene_between_alkenes(
    length: usize,
    states: [u8; 2],
    unknown_chain: bool,
    reversed: bool,
) -> Molecule {
    use nazuchem_ir::{Atom, Bond, BondOrder, Relation, Stereo};
    let orders = [BondOrder::Single, BondOrder::Double, BondOrder::Single]
        .into_iter()
        .chain(std::iter::repeat_n(BondOrder::Double, length))
        .chain([BondOrder::Single, BondOrder::Double, BondOrder::Single])
        .collect::<Vec<_>>();
    let atoms = (0..=orders.len())
        .map(|i| {
            let valence = [i.checked_sub(1), (i < orders.len()).then_some(i)]
                .into_iter()
                .flatten()
                .map(|bond| {
                    if orders[bond] == BondOrder::Double {
                        2
                    } else {
                        1
                    }
                })
                .sum::<u8>();
            Atom {
                element: 6,
                hydrogens: 4 - valence,
                ..Atom::default()
            }
        })
        .collect();
    let bonds = orders
        .into_iter()
        .enumerate()
        .map(|(i, order)| Bond {
            atoms: if reversed {
                [(i + 1) as u32, i as u32]
            } else {
                [i as u32, (i + 1) as u32]
            },
            order,
        })
        .collect();
    let mut molecule = Molecule::new(atoms, bonds).unwrap();
    for (bond, state) in [(1, states[0]), (length + 4, states[1])]
        .into_iter()
        .chain((3..3 + length).filter_map(|bond| unknown_chain.then_some((bond, 3))))
    {
        let relation = match state {
            0 => continue,
            1 => Relation::Cis,
            2 => Relation::Trans,
            3 => Relation::Unknown,
            _ => unreachable!(),
        };
        let refs = [(bond - 1) as u32, (bond + 2) as u32];
        molecule.stereo.push(Stereo::DoubleBond {
            bond: bond as u32,
            refs: if reversed { [refs[1], refs[0]] } else { refs },
            relation,
        });
    }
    molecule
}

#[test]
fn specified_outer_alkenes_cannot_assert_absent_odd_cumulene_stereo() {
    // Explicit IR for (2E,8E)-deca-2,4,5,6,8-pentaene. Its central
    // three-double-bond cumulene has no extended cis/trans assertion.
    let molecule = cumulene_between_alkenes(3, [2, 2], false, false);
    assert_eq!(molecule.atoms.len(), 10);
    assert_eq!(molecule.stereo.len(), 2);
    molecule.validate().unwrap();
    for result in [
        molecule.to_cxsmiles(),
        molecule.semantic_cxsmiles_for_tests(),
    ] {
        assert!(matches!(result, Err(CxsmilesError::StereoNotExpressible)));
    }
    // This was the faulty output: ordinary per-bond parsing sees only the
    // two intended alkenes, but the independent full-chain audit rejects it.
    let bad = r"C\C=C\C=C=C=C\C=C\C";
    let error = stereo::written(&molecule, bad, &(0..10).collect::<Vec<_>>()).unwrap_err();
    assert!(error.contains("extended cis/trans"), "{error}");
}

#[test]
fn odd_and_even_cumulene_specification_patterns_preserve_stereo_when_exported() {
    let mut emitted = 0;
    let mut refused = 0;
    for length in 2..=7 {
        for first in 0..4 {
            for last in 0..4 {
                for unknown_chain in [false, true] {
                    for reversed in [false, true] {
                        let molecule = cumulene_between_alkenes(
                            length,
                            [first, last],
                            unknown_chain,
                            reversed,
                        );
                        let must_refuse =
                            length % 2 == 1 && matches!(first, 1 | 2) && matches!(last, 1 | 2);
                        match molecule.to_cxsmiles() {
                            Ok(text) => {
                                assert!(!must_refuse, "{length}: {text}");
                                stereo::written(
                                    &molecule,
                                    &text,
                                    &molecule.smiles_atom_order().unwrap(),
                                )
                                .unwrap_or_else(|error| panic!("{length}: {text}: {error}"));
                                emitted += 1;
                            }
                            Err(CxsmilesError::StereoNotExpressible) => {
                                assert!(must_refuse, "{length}: {first}/{last}");
                                refused += 1;
                            }
                            Err(error) => panic!("unexpected export error: {error}"),
                        }
                    }
                }
            }
        }
    }
    assert_eq!((emitted, refused), (336, 48));
    eprintln!("cumulene property cases: emitted={emitted}, StereoNotExpressible={refused}");
}

#[test]
fn unknown_tetrahedral_carrier_patterns_preserve_stereo_when_exported() {
    use nazuchem_ir::{Atom, Bond, BondOrder, Ref, Stereo, Winding};
    // Structural IR does not perform valence perception. Exhaust incident
    // single/double carrier choices even for unusual valence patterns, and
    // both record orders, to measure conservative carrier refusals explicitly.
    let mut emitted = 0;
    let mut refused = 0;
    for mask in 0..32 {
        for reversed in [false, true] {
            let atoms = [6, 6, 9, 17, 35, 53]
                .into_iter()
                .enumerate()
                .map(|(id, element)| Atom {
                    element,
                    hydrogens: u8::from(id < 2),
                    ..Atom::default()
                })
                .collect();
            let bonds = [[0, 1], [0, 2], [0, 3], [1, 4], [1, 5]]
                .into_iter()
                .enumerate()
                .map(|(id, atoms)| Bond {
                    atoms,
                    order: if mask & (1 << id) != 0 {
                        BondOrder::Single
                    } else {
                        BondOrder::Double
                    },
                })
                .collect();
            let mut molecule = Molecule::new(atoms, bonds).unwrap();
            for (center, others) in [(0, [1, 2, 3]), (1, [0, 4, 5])] {
                molecule.stereo.push(Stereo::Tetrahedral {
                    center,
                    neighbors: [
                        Ref::ImplicitHydrogen,
                        Ref::Atom(others[0]),
                        Ref::Atom(others[1]),
                        Ref::Atom(others[2]),
                    ],
                    winding: Winding::Unknown,
                });
            }
            if reversed {
                molecule.stereo.reverse();
            }
            match molecule.to_cxsmiles() {
                Ok(text) => {
                    stereo::written(&molecule, &text, &molecule.smiles_atom_order().unwrap())
                        .unwrap_or_else(|error| panic!("mask={mask}: {text}: {error}"));
                    emitted += 1;
                }
                Err(CxsmilesError::StereoNotExpressible) => refused += 1,
                Err(error) => panic!("mask={mask}: unexpected export error: {error}"),
            }
        }
    }
    assert_eq!(emitted + refused, 64);
    assert!(emitted > 0 && refused > 0);
    eprintln!("unknown-carrier property cases: emitted={emitted}, StereoNotExpressible={refused}");
}

#[test]
fn all_conjugated_chain_specification_patterns_preserve_stereo_when_exported() {
    use nazuchem_ir::{Atom, Bond, BondOrder, Relation, Stereo};
    // Exhaust all unspecified/cis/trans choices on one through five alternating
    // double bonds, plus unknown entries (which must not become assertions).
    let mut emitted = 0;
    let mut refused = 0;
    for n in 1..=5 {
        for mut pattern in 0..4usize.pow(n as u32) {
            let atoms = (0..2 * n + 2)
                .map(|i| Atom {
                    element: 6,
                    hydrogens: if i == 0 || i == 2 * n + 1 { 3 } else { 1 },
                    ..Atom::default()
                })
                .collect();
            let bonds = (0..2 * n + 1)
                .map(|i| Bond {
                    atoms: [i as u32, (i + 1) as u32],
                    order: if i % 2 == 0 {
                        BondOrder::Single
                    } else {
                        BondOrder::Double
                    },
                })
                .collect();
            let mut molecule = Molecule::new(atoms, bonds).unwrap();
            for i in 0..n {
                let choice = pattern % 4;
                pattern /= 4;
                if choice > 0 {
                    molecule.stereo.push(Stereo::DoubleBond {
                        bond: (2 * i + 1) as u32,
                        refs: [(2 * i) as u32, (2 * i + 3) as u32],
                        relation: match choice {
                            1 => Relation::Cis,
                            2 => Relation::Trans,
                            _ => Relation::Unknown,
                        },
                    });
                }
            }
            match molecule.to_cxsmiles() {
                Ok(text) => {
                    stereo::written(&molecule, &text, &molecule.smiles_atom_order().unwrap())
                        .unwrap();
                    emitted += 1;
                }
                Err(CxsmilesError::StereoNotExpressible) => refused += 1,
                Err(error) => panic!("unexpected export error: {error}"),
            }
        }
    }
    assert!(emitted > 0 && refused > 0);
    eprintln!("conjugated-chain property cases: emitted={emitted}, StereoNotExpressible={refused}");
}

#[test]
fn direct_source_comparison_rejects_inverted_configurations() {
    use nazuchem_ir::{AxialWinding, Relation, Stereo, Winding};
    let parser = Parser::new().unwrap();
    for name in [
        "(2R)-butan-2-ol",
        "(E)-but-2-ene",
        "(1M)-[1,1'-binaphthalene]-2,2'-diol",
    ] {
        let structure = parser
            .parse(name, &ParseOptions::strict())
            .structure
            .unwrap();
        let (mut molecule, order) = structure.to_ir_with_atom_mapping().unwrap();
        stereo::source(&structure, &molecule, &order).unwrap();
        match &mut molecule.stereo[0] {
            Stereo::Tetrahedral { winding, .. } => {
                *winding = if *winding == Winding::Clockwise {
                    Winding::Anticlockwise
                } else {
                    Winding::Clockwise
                }
            }
            Stereo::DoubleBond { relation, .. } => {
                *relation = if *relation == Relation::Cis {
                    Relation::Trans
                } else {
                    Relation::Cis
                }
            }
            Stereo::Atropisomer { winding, .. } => {
                *winding = if *winding == AxialWinding::Anticlockwise {
                    AxialWinding::Clockwise
                } else {
                    AxialWinding::Anticlockwise
                }
            }
            _ => panic!("new untested stereo variant"),
        }
        assert!(
            stereo::source(&structure, &molecule, &order).is_err(),
            "inversion escaped direct audit: {name}"
        );
    }
}

#[test]
fn branched_and_ring_stereo_survives_reversed_bond_and_neighbor_orders() {
    use nazuchem_ir::Stereo;
    let parser = Parser::new().unwrap();
    for name in [
        "(2R)-butan-2-ol",
        "(1R,2S)-1,2-dichlorocyclohexane",
        "(E)-3-methylhex-3-ene",
        "(E)-cyclodec-1-ene",
    ] {
        let result = parser.parse(name, &ParseOptions::strict());
        let original = result
            .structure
            .unwrap_or_else(|| panic!("{name}: {}", result.message))
            .to_ir()
            .unwrap();
        for reverse in [false, true] {
            let mut molecule = original.clone();
            if reverse {
                for bond in &mut molecule.bonds {
                    bond.atoms.swap(0, 1);
                }
                for unit in &mut molecule.stereo {
                    match unit {
                        Stereo::Tetrahedral { neighbors, .. } => neighbors.rotate_left(2), // even permutation
                        Stereo::DoubleBond { refs, .. } => refs.swap(0, 1),
                        _ => {}
                    }
                }
            }
            let text = molecule.to_cxsmiles().unwrap();
            stereo::written(&molecule, &text, &molecule.smiles_atom_order().unwrap())
                .unwrap_or_else(|error| panic!("{name}: {error}: {text}"));
        }
    }
}

#[test]
fn adjacent_unknown_centers_use_distinct_parseable_wavy_bonds() {
    use nazuchem_ir::{Stereo, Winding};
    let structure = Parser::new()
        .unwrap()
        .parse("(2R,3S)-2,3-dichlorobutane", &ParseOptions::strict())
        .structure
        .unwrap();
    let original = structure.to_ir().unwrap();
    let centers = original
        .stereo
        .iter()
        .filter_map(|unit| match unit {
            Stereo::Tetrahedral { center, .. } => Some(*center),
            _ => None,
        })
        .collect::<Vec<_>>();
    assert_eq!(centers.len(), 2);
    let shared = original
        .bonds
        .iter()
        .position(|bond| centers.iter().all(|atom| bond.atoms.contains(atom)))
        .unwrap();
    for unknown_mask in 0..4 {
        let mut molecule = original.clone();
        // Make the shared carrier the first choice at both centers. Selecting
        // it twice emits invalid duplicate CX w entries, even though each
        // center independently has a usable single bond.
        molecule.bonds.swap(0, shared);
        for (index, unit) in molecule.stereo.iter_mut().enumerate() {
            if unknown_mask & (1 << index) != 0 {
                let Stereo::Tetrahedral { winding, .. } = unit else {
                    unreachable!()
                };
                *winding = Winding::Unknown;
            }
        }
        let text = molecule.to_cxsmiles().unwrap();
        stereo::written(&molecule, &text, &molecule.smiles_atom_order().unwrap())
            .unwrap_or_else(|error| panic!("{error}: {text}"));
    }
}
