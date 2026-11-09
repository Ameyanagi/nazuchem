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
use nazuchem_ir::Molecule;
use nazuchem_kernel::{budget::Budget, identity::canonical_identity, smiles};
use nazuchem_n2s::{IrError, ParseOptions, ParseResult, Parser, Status, Structure};
use serde_json::{Value, json};

#[path = "support/ir_cxsmiles.rs"]
mod cx;
#[path = "support/ir_schema.rs"]
mod schema;

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
            self.supported + self.axial,
            success + warning,
            "every parsed structure must be checked: {self:?}"
        );
        assert_eq!(
            self.exact + self.canonical,
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
    let prepared = smiles::prepare(cx::bare(text), &budget).map_err(|error| error.to_string())?;
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
    if !axes.is_empty() {
        return match structure.to_ir() {
            Err(IrError::UnsupportedAxialStereo) => {
                counts.axial += 1;
                Ok(())
            }
            other => Err(format!(
                "axial stereo must return typed UnsupportedAxialStereo, got {other:?}"
            )),
        };
    }
    let molecule = structure
        .to_ir()
        .map_err(|error| format!("IR projection: {error}"))?;
    molecule
        .validate()
        .map_err(|error| format!("IR validation: {error}"))?;
    if molecule.atoms.len() != legacy_locants.len() {
        return Err("IR atom identities differ from legacy written atoms".into());
    }
    cx::legacy_radicals(&molecule, &legacy)?;
    if molecule.displayed_locants() != legacy_locants {
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
    let actual = molecule
        .semantic_cxsmiles()
        .map_err(|error| format!("IR semantic writer: {error}"))?;
    if restored
        .semantic_cxsmiles()
        .map_err(|error| error.to_string())?
        != actual
        || restored.displayed_locants() != legacy_locants
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
    cx::written_graph(&molecule, &actual, &order)?;
    let expected_layers = cx::semantic_layers(&legacy, &identity)?;
    let actual_layers = cx::semantic_layers(&actual, &order)?;
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
    counts.assert_baseline(1650, 1609, 0, 41);
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
