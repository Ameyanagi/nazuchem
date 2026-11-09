// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | crates/nazuchem-s2n/src/ir_proof_tests.rs | Phase B step 8b

//! Corpus proof that localized IR bonds retain naming and aromatic identity.
//!
//! Naming rejection does not bypass the aromatic graph comparison. Every
//! distinct source SMILES is accounted for, including inputs the test converter
//! cannot represent, and every parsed n2s structure enters `Structure::to_ir`.

use std::collections::BTreeMap;
use std::io::{BufRead, BufReader};

use flate2::read::GzDecoder;
use nazuchem_ir::IsotopeStyle;
use nazuchem_kernel::budget::Budget;
use nazuchem_n2s::{ParseOptions, ParseResult, Parser, Status, Structure};
use serde_json::Value;

use crate::{Name, graph_io, ir_input, molecule, name_molecule, name_smiles};

#[derive(Debug, Default, PartialEq, Eq)]
struct NameCounts {
    records: usize,
    identical: usize,
    differing: usize,
    both_rejected: usize,
    one_rejected: usize,
    rejection_causes: BTreeMap<String, usize>,
}

impl NameCounts {
    fn merge(&mut self, other: Self) {
        let Self {
            records,
            identical,
            differing,
            both_rejected,
            one_rejected,
            rejection_causes,
        } = other;
        self.records += records;
        self.identical += identical;
        self.differing += differing;
        self.both_rejected += both_rejected;
        self.one_rejected += one_rejected;
        for (cause, count) in rejection_causes {
            *self.rejection_causes.entry(cause).or_default() += count;
        }
    }

    fn compare(
        &mut self,
        label: &str,
        legacy: Result<Name, String>,
        adapted: Result<Name, String>,
        differences: &mut Vec<String>,
    ) {
        self.records += 1;
        match (legacy, adapted) {
            (Ok(legacy), Ok(adapted)) => {
                // Compare every serialized field, including status reasons,
                // alternative order, omitted units and recovery verification.
                let legacy = serde_json::to_value(legacy).expect("Name serializes");
                let adapted = serde_json::to_value(adapted).expect("Name serializes");
                if legacy == adapted {
                    self.identical += 1;
                } else {
                    self.differing += 1;
                    difference(
                        differences,
                        format!("{label}: complete Name differs\nSMILES {legacy}\nIR {adapted}"),
                    );
                }
            }
            (Err(_), Err(_)) => self.both_rejected += 1,
            (legacy, adapted) => {
                self.one_rejected += 1;
                let cause = match (&legacy, &adapted) {
                    (Err(error), Ok(_)) => format!("SMILES rejected: {error}"),
                    (Ok(_), Err(error)) => format!("IR rejected: {error}"),
                    _ => unreachable!("two-sided results handled above"),
                };
                *self.rejection_causes.entry(cause).or_default() += 1;
                difference(
                    differences,
                    format!(
                        "{label}: naming acceptance differs\nSMILES {legacy:?}\nIR {adapted:?}"
                    ),
                );
            }
        }
    }

    fn assert_complete(&self) {
        assert_eq!(
            self.records,
            self.identical + self.differing + self.both_rejected + self.one_rejected,
            "every naming comparison must be accounted for: {self:?}"
        );
        assert_eq!((self.differing, self.one_rejected), (0, 0), "{self:?}");
    }
}

#[derive(Debug, PartialEq, Eq)]
struct AromaticSnapshot {
    atoms: Vec<(usize, bool)>,
    // Localized pi orders are checked independently of aromatic flags. An
    // inter-ring bond between aromatic atoms must retain its own identity.
    bonds: Vec<(usize, usize, usize, usize, bool)>,
}

fn aromatic_snapshot(molecule: &molecule::Molecule) -> AromaticSnapshot {
    AromaticSnapshot {
        atoms: molecule
            .atoms
            .iter()
            .map(|(&id, atom)| (id, atom.is_aromatic))
            .collect(),
        bonds: molecule
            .bonds
            .iter()
            .map(|(&id, bond)| (id, bond.u, bond.v, bond.order, bond.is_aromatic))
            .collect(),
    }
}

#[derive(Debug, Default, PartialEq, Eq)]
struct SmilesCounts {
    names: NameCounts,
    rejection_messages_identical: usize,
    rejection_messages_differing: usize,
    converter_rejected: usize,
    perception_identical: usize,
    perception_differing: usize,
    perception_rejected: usize,
    graph_identical: usize,
    graph_differing: usize,
    graph_both_rejected: usize,
    graph_one_rejected: usize,
    aromatic_graphs: usize,
    aromatic_atoms: usize,
    aromatic_bonds: usize,
    localized_aromatic_endpoint_bonds: usize,
}

fn difference(differences: &mut Vec<String>, text: String) {
    if std::env::var_os("NAZUCHEM_IR_PROOF_TRACE").is_some() {
        eprintln!("IR proof difference: {text}");
    }
    // Repeated early acceptance failures must not hide successful names whose
    // status/alternatives/omitted units changed later in the same corpus.
    if differences.len() < 40 || (differences.len() < 80 && text.contains("complete Name differs"))
    {
        differences.push(text);
    }
}

fn finish(label: &str, counts: &impl std::fmt::Debug, differences: Vec<String>) {
    eprintln!("{label}: {counts:?}");
    assert!(
        differences.is_empty(),
        "{label} IR/localized-pi proof failed; first differences:\n{}",
        differences.join("\n")
    );
}

#[test]
fn every_distinct_smiles_retains_complete_names_and_aromatic_graphs() {
    let budget = Budget::default();
    let mut inputs = BTreeMap::<String, Vec<String>>::new();
    for (source, fixture, rows) in [
        (
            "bakeoff",
            include_str!("../../../tests/fixtures/s2n/bakeoff_corpus.tsv"),
            703,
        ),
        (
            "openclatura",
            include_str!("../../../tests/fixtures/s2n/openclatura_tests.tsv"),
            894,
        ),
        (
            "rules",
            include_str!("../../../tests/fixtures/s2n/rules_cases.tsv"),
            94,
        ),
    ] {
        assert_eq!(
            fixture.lines().next().unwrap().split('\t').next(),
            Some("smiles")
        );
        assert_eq!(
            fixture.lines().count() - 1,
            rows,
            "{source} coverage changed"
        );
        for (line, record) in fixture.lines().enumerate().skip(1) {
            assert!(
                !record.is_empty(),
                "{source}:{}: empty fixture row",
                line + 1
            );
            let smiles = record.split('\t').next().expect("SMILES column");
            assert!(!smiles.is_empty(), "{source}:{}: missing SMILES", line + 1);
            inputs
                .entry(smiles.into())
                .or_default()
                .push(format!("{source}:{}", line + 1));
        }
    }
    for (line, record) in include_str!("../../../tests/fixtures/golden/s2n_n2s.jsonl")
        .lines()
        .enumerate()
    {
        let record: Value = serde_json::from_str(record).expect("golden JSON");
        let smiles = record["smiles"].as_str().expect("golden source SMILES");
        inputs
            .entry(smiles.into())
            .or_default()
            .push(format!("golden:{}", line + 1));
    }

    let (mut counts, mut differences) = (SmilesCounts::default(), Vec::new());
    let mut conversion_errors = BTreeMap::<String, usize>::new();
    for (text, sources) in inputs {
        let label = format!("{text} ({})", sources.join(", "));
        let legacy_name = name_smiles(&text, &budget);
        match ir_input::from_smiles(&text, &budget) {
            Ok(ir) => {
                // This comparison precedes naming admission, so disconnected
                // salts and mixtures are covered even when both graph readers
                // reject them. Localized naming orders are checked below.
                let legacy_state = nazuchem_kernel::smiles::prepare(&text, &budget)
                    .map(|prepared| prepared.state)
                    .map_err(|error| error.to_string());
                let adapted_state = ir_input::prepare(&ir, &budget);
                match (legacy_state, adapted_state) {
                    (Ok(legacy), Ok(adapted)) => {
                        let snapshot = |state: &nazuchem_kernel::stereo::perception::State| {
                            (
                                state
                                    .graph
                                    .atoms
                                    .iter()
                                    .map(|atom| atom.aromatic)
                                    .collect::<Vec<_>>(),
                                state
                                    .graph
                                    .bonds
                                    .iter()
                                    .map(|bond| (bond.a, bond.b, bond.aromatic))
                                    .collect::<Vec<_>>(),
                            )
                        };
                        let before = snapshot(&legacy);
                        let after = snapshot(&adapted);
                        if before == after {
                            counts.perception_identical += 1;
                        } else {
                            counts.perception_differing += 1;
                            difference(
                                &mut differences,
                                format!(
                                    "{label}: pre-admission aromatic snapshot differs\nSMILES {before:?}\nIR {after:?}"
                                ),
                            );
                        }
                    }
                    (legacy, adapted) => {
                        counts.perception_rejected += 1;
                        difference(
                            &mut differences,
                            format!(
                                "{label}: source/IR aromatic preparation unavailable\nSMILES {legacy:?}\nIR {adapted:?}"
                            ),
                        );
                    }
                }
                // Graph comparison runs even if either naming route rejects.
                let legacy = graph_io::read_smiles(&text, &budget);
                let adapted = graph_io::read_molecule(&ir, &budget);
                match (&legacy, &adapted) {
                    (Ok(legacy), Ok(adapted)) => {
                        let before = aromatic_snapshot(legacy);
                        let after = aromatic_snapshot(adapted);
                        if before == after {
                            counts.graph_identical += 1;
                        } else {
                            counts.graph_differing += 1;
                            difference(
                                &mut differences,
                                format!(
                                    "{label}: aromatic/localized bond snapshot differs\nSMILES {before:?}\nIR {after:?}"
                                ),
                            );
                        }
                        let atoms = legacy
                            .atoms
                            .values()
                            .filter(|atom| atom.is_aromatic)
                            .count();
                        counts.aromatic_graphs += usize::from(atoms != 0);
                        counts.aromatic_atoms += atoms;
                        counts.aromatic_bonds += legacy
                            .bonds
                            .values()
                            .filter(|bond| bond.is_aromatic)
                            .count();
                        counts.localized_aromatic_endpoint_bonds += legacy
                            .bonds
                            .values()
                            .filter(|bond| {
                                !bond.is_aromatic
                                    && legacy.atoms[&bond.u].is_aromatic
                                    && legacy.atoms[&bond.v].is_aromatic
                            })
                            .count();
                    }
                    (Err(_), Err(_)) => counts.graph_both_rejected += 1,
                    _ => {
                        counts.graph_one_rejected += 1;
                        difference(
                            &mut differences,
                            format!(
                                "{label}: graph acceptance differs\nSMILES {legacy:?}\nIR {adapted:?}"
                            ),
                        );
                    }
                }
                let adapted_name = name_molecule(&ir, IsotopeStyle::Iupac, &budget);
                if let (Err(before), Err(after)) = (&legacy_name, &adapted_name) {
                    if before == after {
                        counts.rejection_messages_identical += 1;
                    } else {
                        counts.rejection_messages_differing += 1;
                        difference(
                            &mut differences,
                            format!(
                                "{label}: same-source rejection message differs\nSMILES {before}\nIR {after}"
                            ),
                        );
                    }
                }
                counts
                    .names
                    .compare(&label, legacy_name, adapted_name, &mut differences);
            }
            Err(error) => {
                // Invalid/unsupported source strings cannot be fabricated into
                // an IR. Require legacy rejection and expose the count/reason.
                counts.converter_rejected += 1;
                *conversion_errors.entry(error.clone()).or_default() += 1;
                counts.names.compare(
                    &label,
                    legacy_name,
                    Err(format!("source-to-IR: {error}")),
                    &mut differences,
                );
            }
        }
    }
    eprintln!("source-to-IR rejections: {conversion_errors:?}");
    finish("distinct SMILES", &counts, differences);
    counts.names.assert_complete();
    assert_eq!(counts.graph_differing, 0);
    assert_eq!(counts.graph_one_rejected, 0);
    assert_eq!(
        counts.names.records,
        counts.converter_rejected + counts.graph_identical + counts.graph_both_rejected
    );
    assert_eq!(counts.names.records, 1650);
    assert_eq!(
        (counts.names.identical, counts.names.both_rejected),
        (1616, 34)
    );
    assert_eq!(counts.converter_rejected, 0);
    assert_eq!(
        (
            counts.rejection_messages_identical,
            counts.rejection_messages_differing
        ),
        (34, 0)
    );
    assert_eq!(
        (
            counts.perception_identical,
            counts.perception_differing,
            counts.perception_rejected
        ),
        (1650, 0, 0)
    );
    assert_eq!(
        (counts.graph_identical, counts.graph_both_rejected),
        (1637, 13)
    );
    assert_eq!(
        (
            counts.aromatic_graphs,
            counts.aromatic_atoms,
            counts.aromatic_bonds,
            counts.localized_aromatic_endpoint_bonds
        ),
        (987, 11212, 11760, 149)
    );
}

#[derive(Debug, Default, PartialEq, Eq)]
struct N2sCounts {
    records: usize,
    success: usize,
    warning: usize,
    failure: usize,
    ir_unavailable: usize,
    legacy_serialization_unavailable: usize,
    names: NameCounts,
}

impl N2sCounts {
    fn merge(&mut self, other: Self) {
        let Self {
            records,
            success,
            warning,
            failure,
            ir_unavailable,
            legacy_serialization_unavailable,
            names,
        } = other;
        self.records += records;
        self.success += success;
        self.warning += warning;
        self.failure += failure;
        self.ir_unavailable += ir_unavailable;
        self.legacy_serialization_unavailable += legacy_serialization_unavailable;
        self.names.merge(names);
    }
}

fn compare_n2s(
    label: &str,
    result: ParseResult,
    expected_status: &str,
    budget: &Budget,
    counts: &mut N2sCounts,
    differences: &mut Vec<String>,
) {
    counts.records += 1;
    let actual = match result.status {
        Status::Success => {
            counts.success += 1;
            "SUCCESS"
        }
        Status::Warning => {
            counts.warning += 1;
            "WARNING"
        }
        Status::Failure => {
            counts.failure += 1;
            "FAILURE"
        }
    };
    if actual != expected_status {
        difference(
            differences,
            format!("{label}: parse status changed: {expected_status} -> {actual}"),
        );
    }
    match result.structure {
        Some(structure) => compare_structure(label, &structure, budget, counts, differences),
        None if actual != "FAILURE" => difference(
            differences,
            format!("{label}: successful parse omitted its structure"),
        ),
        None => {}
    }
}

fn compare_structure(
    label: &str,
    structure: &Structure,
    budget: &Budget,
    counts: &mut N2sCounts,
    differences: &mut Vec<String>,
) {
    // Opt-in tracing identifies unusually expensive corpus records without
    // changing the default gate's output or excluding those structures.
    if std::env::var_os("NAZUCHEM_IR_PROOF_TRACE").is_some() {
        eprintln!("IR naming proof record {}: {label}", counts.records);
    }
    let ir = match structure.to_ir() {
        Ok(ir) => ir,
        Err(error) => {
            counts.ir_unavailable += 1;
            difference(
                differences,
                format!("{label}: Structure::to_ir failed: {error}"),
            );
            return;
        }
    };
    let adapted = name_molecule(&ir, IsotopeStyle::Iupac, budget);
    match structure.semantic_cxsmiles() {
        Ok(text) => counts
            .names
            .compare(label, name_smiles(&text, budget), adapted, differences),
        Err(error) => {
            // Axes and other unsupported legacy exports are explicitly counted;
            // an unavailable baseline must never masquerade as a comparison.
            counts.legacy_serialization_unavailable += 1;
            difference(
                differences,
                format!(
                    "{label}: legacy semantic CXSMILES unavailable: {error}; IR result: {adapted:?}"
                ),
            );
        }
    }
}

fn assert_n2s_coverage(counts: &N2sCounts, baseline: (usize, usize, usize, usize)) {
    assert_eq!(
        (
            counts.records,
            counts.success,
            counts.warning,
            counts.failure
        ),
        baseline
    );
    assert_eq!(counts.ir_unavailable, 0);
    assert_eq!(counts.legacy_serialization_unavailable, 0);
    assert_eq!(counts.names.records, counts.success + counts.warning);
    counts.names.assert_complete();
}

#[test]
fn every_n2s_structure_oracle_structure_retains_complete_name_results() {
    let reader = BufReader::new(GzDecoder::new(
        include_bytes!("../../nazuchem-n2s/tests/fixtures/structure-golden.jsonl.gz").as_slice(),
    ));
    // Only the test harness uses threads. Every record still invokes to_ir and
    // both naming routes; separate parsers and thread-local fusion caches keep
    // requests independent. Four bounded shards reduce the full oracle gate's
    // wall time, and joining in shard order keeps diagnostics deterministic.
    let mut shards: [Vec<(usize, Value)>; 4] = std::array::from_fn(|_| Vec::new());
    for (row, line) in reader.lines().enumerate() {
        let record =
            serde_json::from_str(&line.expect("fixture decompression")).expect("fixture JSON");
        shards[row % 4].push((row + 1, record));
    }
    let results = std::thread::scope(|scope| {
        let workers: Vec<_> = shards
            .into_iter()
            .map(|rows| {
                scope.spawn(move || {
                    let parser = Parser::new().expect("embedded OPSIN resources");
                    let budget = Budget::default();
                    let (mut counts, mut differences) = (N2sCounts::default(), Vec::new());
                    for (row, record) in rows {
                        let name = record["input"].as_str().unwrap();
                        let result = parser.parse(
                            name,
                            &ParseOptions {
                                allow_radicals: record["allow_radicals"].as_bool().unwrap(),
                                ..ParseOptions::strict()
                            },
                        );
                        compare_n2s(
                            &format!("row {row}, {}: {name}", record["family"]),
                            result,
                            record["status"].as_str().unwrap(),
                            &budget,
                            &mut counts,
                            &mut differences,
                        );
                    }
                    (counts, differences)
                })
            })
            .collect();
        workers
            .into_iter()
            .map(|worker| worker.join().expect("structure corpus shard did not panic"))
            .collect::<Vec<_>>()
    });
    let (mut counts, mut differences) = (N2sCounts::default(), Vec::new());
    for (shard_counts, shard_differences) in results {
        counts.merge(shard_counts);
        for text in shard_differences {
            difference(&mut differences, text);
        }
    }
    finish("n2s structure oracle", &counts, differences);
    assert_n2s_coverage(&counts, (2624, 2535, 70, 19));
    assert_eq!(
        (counts.names.identical, counts.names.both_rejected),
        (2293, 312)
    );
}

#[test]
fn every_n2s_configuration_oracle_structure_retains_complete_name_results() {
    let parser = Parser::new().expect("embedded OPSIN resources");
    let budget = Budget::default();
    let reader = BufReader::new(GzDecoder::new(
        include_bytes!("../../nazuchem-n2s/tests/fixtures/configuration-golden.jsonl.gz")
            .as_slice(),
    ));
    let (mut counts, mut differences) = (N2sCounts::default(), Vec::new());
    for line in reader.lines() {
        let record: Value =
            serde_json::from_str(&line.expect("fixture decompression")).expect("fixture JSON");
        let name = record["input"].as_str().unwrap();
        let options = &record["options"];
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
        compare_n2s(
            &format!("{name} mask={}", record["mask"]),
            result,
            record["status"].as_str().unwrap(),
            &budget,
            &mut counts,
            &mut differences,
        );
    }
    finish("n2s configuration oracle", &counts, differences);
    assert_n2s_coverage(&counts, (640, 304, 112, 224));
    assert_eq!(
        (counts.names.identical, counts.names.both_rejected),
        (336, 80)
    );
}

#[test]
fn every_n2s_golden_structure_retains_complete_name_results() {
    let parser = Parser::new().expect("embedded OPSIN resources");
    let budget = Budget::default();
    let (mut counts, mut differences, mut records) = (N2sCounts::default(), Vec::new(), 0);
    for line in include_str!("../../../tests/fixtures/golden/s2n_n2s.jsonl").lines() {
        let record: Value = serde_json::from_str(line).expect("golden JSON");
        records += 1;
        for expected in record["n2s"].as_array().unwrap() {
            let name = expected["name"].as_str().unwrap();
            compare_n2s(
                name,
                parser.parse(name, &ParseOptions::strict()),
                expected["status"].as_str().unwrap(),
                &budget,
                &mut counts,
                &mut differences,
            );
        }
    }
    finish("n2s generated-name golden", &counts, differences);
    assert_eq!(records, 1650);
    assert_n2s_coverage(&counts, (1650, 1609, 0, 41));
    assert_eq!(
        (counts.names.identical, counts.names.both_rejected),
        (1608, 1)
    );
}
