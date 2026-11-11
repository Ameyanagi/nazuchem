// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | crates/nazuchem-s2n/src/ir_proof_tests.rs | Phase B step 8b

//! Corpus proof that localized IR bonds retain naming and aromatic identity.
//!
//! Naming rejection does not bypass the aromatic graph comparison. Every
//! distinct source SMILES is accounted for, including inputs the test converter
//! cannot represent, and every parsed n2s structure enters `Structure::to_ir`.

use std::cell::RefCell;
use std::collections::BTreeMap;
use std::io::{BufRead, BufReader};
use std::time::Instant;

use flate2::read::GzDecoder;
use nazuchem_ir::IsotopeStyle;
use nazuchem_kernel::budget::Budget;
use nazuchem_n2s::{ParseOptions, ParseResult, Parser, Status, Structure};
use serde_json::Value;

use crate::{Name, graph_io, ir_input, molecule, name_molecule, name_smiles};

#[derive(Debug, Default)]
struct RouteCounts {
    fallbacks: BTreeMap<String, usize>,
    monatomic_probes: usize,
}

thread_local! {
    // Only the compared IR call enables collection. Other tests sharing this
    // module, and other test-harness threads, cannot contaminate the counters.
    static ROUTE_COUNTS: RefCell<Option<RouteCounts>> = const { RefCell::new(None) };
}

pub(crate) fn record_fallback(path: &str) {
    ROUTE_COUNTS.with_borrow_mut(|counts| {
        if let Some(counts) = counts {
            *counts.fallbacks.entry(path.into()).or_default() += 1;
        }
    });
}

pub(crate) fn record_monatomic_probe() {
    ROUTE_COUNTS.with_borrow_mut(|counts| {
        if let Some(counts) = counts {
            counts.monatomic_probes += 1;
        }
    });
}

struct RouteScope;

impl Drop for RouteScope {
    fn drop(&mut self) {
        // Clear collection even when a comparison panics.
        ROUTE_COUNTS.with_borrow_mut(|counts| *counts = None);
    }
}

fn compared_name_molecule(
    molecule: &nazuchem_ir::Molecule,
    budget: &Budget,
    counts: &mut NameCounts,
) -> Result<Name, String> {
    ROUTE_COUNTS.with_borrow_mut(|current| {
        assert!(current.is_none(), "proof route collection cannot nest");
        *current = Some(RouteCounts::default());
    });
    let _scope = RouteScope;
    let result = name_molecule(molecule, IsotopeStyle::Iupac, budget);
    let routes = ROUTE_COUNTS.with_borrow_mut(|current| current.take().unwrap());
    counts.ir_calls += 1;
    counts.monatomic_probes += routes.monatomic_probes;
    counts.fallback_comparisons += usize::from(!routes.fallbacks.is_empty());
    for (path, count) in routes.fallbacks {
        *counts.fallback_paths.entry(path).or_default() += count;
    }
    result
}

#[derive(Debug, Default, PartialEq, Eq)]
struct NameCounts {
    records: usize,
    identical: usize,
    differing: usize,
    both_rejected: usize,
    one_rejected: usize,
    rejection_causes: BTreeMap<String, usize>,
    ir_calls: usize,
    fallback_comparisons: usize,
    fallback_paths: BTreeMap<String, usize>,
    monatomic_probes: usize,
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
            ir_calls,
            fallback_comparisons,
            fallback_paths,
            monatomic_probes,
        } = other;
        self.records += records;
        self.identical += identical;
        self.differing += differing;
        self.both_rejected += both_rejected;
        self.one_rejected += one_rejected;
        self.ir_calls += ir_calls;
        self.fallback_comparisons += fallback_comparisons;
        self.monatomic_probes += monatomic_probes;
        for (cause, count) in rejection_causes {
            *self.rejection_causes.entry(cause).or_default() += count;
        }
        for (path, count) in fallback_paths {
            *self.fallback_paths.entry(path).or_default() += count;
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
        assert_eq!(
            self.fallback_comparisons,
            self.fallback_paths.values().sum::<usize>(),
            "each compared IR call selects at most one fallback: {self:?}"
        );
        assert!(self.monatomic_probes <= self.ir_calls, "{self:?}");
    }

    fn assert_routes(&self, paths: &[(&str, usize)], monatomic_probes: usize) {
        assert_eq!(self.ir_calls, self.records, "{self:?}");
        assert_eq!(
            self.fallback_paths,
            paths
                .iter()
                .map(|&(path, count)| (path.into(), count))
                .collect(),
            "{self:?}"
        );
        assert_eq!(self.monatomic_probes, monatomic_probes, "{self:?}");
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

fn finish(label: &str, counts: &impl std::fmt::Debug, differences: Vec<String>, started: Instant) {
    eprintln!(
        "{label}: {counts:?}; elapsed_seconds={:.3}",
        started.elapsed().as_secs_f64()
    );
    assert!(
        differences.is_empty(),
        "{label} IR/localized-pi proof failed; first differences:\n{}",
        differences.join("\n")
    );
}

fn distinct_smiles_inputs() -> BTreeMap<String, Vec<String>> {
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
    inputs
}

#[test]
fn every_distinct_smiles_retains_complete_names_and_aromatic_graphs() {
    let started = Instant::now();
    let budget = Budget::default();
    let (mut counts, mut differences) = (SmilesCounts::default(), Vec::new());
    let mut conversion_errors = BTreeMap::<String, usize>::new();
    for (text, sources) in distinct_smiles_inputs() {
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
                let adapted_name = compared_name_molecule(&ir, &budget, &mut counts.names);
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
    finish("distinct SMILES", &counts, differences, started);
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
        (1620, 30)
    );
    counts
        .names
        .assert_routes(&[("monatomic cation", 1), ("salt", 11)], 1638);
    assert_eq!(counts.converter_rejected, 0);
    assert_eq!(
        (
            counts.rejection_messages_identical,
            counts.rejection_messages_differing
        ),
        (30, 0)
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

#[derive(Debug, Default)]
struct ExplicitHydrogenCounts {
    sources: usize,
    converter_rejected: usize,
    stereo_excluded: usize,
    no_expandable_hydrogens: usize,
    explicit_hydrogen_vertices: usize,
    radical_sources: usize,
    radical_stereo_excluded: usize,
    radical_no_expandable_hydrogens: usize,
    radical_compared: usize,
    radical_identical: usize,
    radical_both_rejected: usize,
    names: NameCounts,
}

impl ExplicitHydrogenCounts {
    fn merge(&mut self, other: Self) {
        self.sources += other.sources;
        self.converter_rejected += other.converter_rejected;
        self.stereo_excluded += other.stereo_excluded;
        self.no_expandable_hydrogens += other.no_expandable_hydrogens;
        self.explicit_hydrogen_vertices += other.explicit_hydrogen_vertices;
        self.radical_sources += other.radical_sources;
        self.radical_stereo_excluded += other.radical_stereo_excluded;
        self.radical_no_expandable_hydrogens += other.radical_no_expandable_hydrogens;
        self.radical_compared += other.radical_compared;
        self.radical_identical += other.radical_identical;
        self.radical_both_rejected += other.radical_both_rejected;
        self.names.merge(other.names);
    }
}

fn ordinary_hydrogen_variant(source: &nazuchem_ir::Molecule) -> nazuchem_ir::Molecule {
    use nazuchem_ir::{Atom, Bond, BondOrder, Molecule};

    assert!(source.stereo.is_empty() && source.stereo_groups.is_empty());
    assert!(source.position_variation.is_empty() && source.polymers.is_empty());
    let mut variant = Molecule::default();
    let mut indices = Vec::with_capacity(source.atoms.len());
    for atom in &source.atoms {
        let mut atom = atom.clone();
        // Materialize one ordinary H at every H-bearing heavy atom, including
        // radical centers. Interleave H before its parent so folding must retain
        // sparse caller identities; methyl becomes H(0)-C(1), with two count H.
        let hydrogen = (atom.element > 1 && atom.hydrogens > 0).then(|| {
            atom.hydrogens -= 1;
            let hydrogen = variant.atoms.len() as u32;
            variant.atoms.push(Atom {
                element: 1,
                ..Atom::default()
            });
            hydrogen
        });
        let index = variant.atoms.len() as u32;
        indices.push(index);
        variant.atoms.push(atom);
        if let Some(hydrogen) = hydrogen {
            variant.bonds.push(Bond {
                atoms: [hydrogen, index],
                order: BondOrder::Single,
            });
        }
    }
    variant.bonds.extend(source.bonds.iter().map(|bond| Bond {
        atoms: [
            indices[bond.atoms[0] as usize],
            indices[bond.atoms[1] as usize],
        ],
        order: bond.order,
    }));
    variant.parent = source
        .parent
        .iter()
        .map(|&index| indices[index as usize])
        .collect();
    variant.validate().expect("ordinary-H variant is valid IR");
    variant
}

fn compare_ordinary_hydrogen_variant(
    label: &str,
    text: &str,
    ir: &nazuchem_ir::Molecule,
    budget: &Budget,
    counts: &mut ExplicitHydrogenCounts,
    differences: &mut Vec<String>,
) {
    let radical = ir.atoms.iter().any(|atom| atom.radical_electrons > 0);
    counts.radical_sources += usize::from(radical);
    if !ir.stereo.is_empty() || !ir.stereo_groups.is_empty() {
        counts.stereo_excluded += 1;
        counts.radical_stereo_excluded += usize::from(radical);
        return;
    }
    if !ir
        .atoms
        .iter()
        .any(|atom| atom.element > 1 && atom.hydrogens > 0)
    {
        counts.no_expandable_hydrogens += 1;
        counts.radical_no_expandable_hydrogens += usize::from(radical);
        return;
    }
    let variant = ordinary_hydrogen_variant(ir);
    counts.explicit_hydrogen_vertices += variant.atoms.len() - ir.atoms.len();
    let legacy = name_smiles(text, budget);
    let adapted = compared_name_molecule(&variant, budget, &mut counts.names);
    let before = (counts.names.identical, counts.names.both_rejected);
    counts.names.compare(label, legacy, adapted, differences);
    if radical {
        counts.radical_compared += 1;
        counts.radical_identical += counts.names.identical - before.0;
        counts.radical_both_rejected += counts.names.both_rejected - before.1;
    }
}

#[test]
fn stereo_free_corpus_ordinary_hydrogen_variants_retain_complete_name_results() {
    let started = Instant::now();
    // The baseline helper reads prepared SMILES and therefore folds source H
    // vertices already. Explicitly materializing H independently here covers
    // representations absent from that helper, including all eligible radicals.
    // Eligibility is exhaustive and structural: no stereo entries/groups, and
    // at least one aggregate H on a non-H, non-wildcard atom. Rejected names stay
    // in the comparison; successful naming is never an eligibility condition.
    let mut shards: [Vec<(String, Vec<String>)>; 4] = std::array::from_fn(|_| Vec::new());
    for (row, input) in distinct_smiles_inputs().into_iter().enumerate() {
        shards[row % 4].push(input);
    }
    let results = std::thread::scope(|scope| {
        let workers: Vec<_> = shards
            .into_iter()
            .map(|rows| {
                scope.spawn(move || {
                    let budget = Budget::default();
                    let (mut counts, mut differences) =
                        (ExplicitHydrogenCounts::default(), Vec::new());
                    for (text, sources) in rows {
                        counts.sources += 1;
                        let label = format!("ordinary-H {text} ({})", sources.join(", "));
                        let ir = match ir_input::from_smiles(&text, &budget) {
                            Ok(ir) => ir,
                            Err(error) => {
                                counts.converter_rejected += 1;
                                difference(&mut differences, format!("{label}: {error}"));
                                continue;
                            }
                        };
                        compare_ordinary_hydrogen_variant(
                            &label,
                            &text,
                            &ir,
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
            .map(|worker| {
                worker
                    .join()
                    .expect("ordinary-H corpus shard did not panic")
            })
            .collect::<Vec<_>>()
    });
    let (mut counts, mut differences) = (ExplicitHydrogenCounts::default(), Vec::new());
    for (shard_counts, shard_differences) in results {
        counts.merge(shard_counts);
        for text in shard_differences {
            difference(&mut differences, text);
        }
    }
    finish("ordinary-H variants", &counts, differences, started);
    counts.names.assert_complete();
    assert_eq!(counts.sources, 1650);
    assert_eq!(counts.converter_rejected, 0);
    assert_eq!(
        (
            counts.stereo_excluded,
            counts.no_expandable_hydrogens,
            counts.explicit_hydrogen_vertices
        ),
        (266, 8, 10953)
    );
    assert_eq!(
        (
            counts.names.records,
            counts.names.identical,
            counts.names.both_rejected
        ),
        (1376, 1348, 28)
    );
    counts.names.assert_routes(&[("salt", 9)], 1366);
    assert_eq!(
        (
            counts.radical_sources,
            counts.radical_stereo_excluded,
            counts.radical_no_expandable_hydrogens,
            counts.radical_compared,
            counts.radical_identical,
            counts.radical_both_rejected
        ),
        (5, 0, 0, 5, 1, 4)
    );
    assert_eq!(
        counts.sources,
        counts.converter_rejected
            + counts.stereo_excluded
            + counts.no_expandable_hydrogens
            + counts.names.records
    );
    assert_eq!(
        counts.radical_sources,
        counts.radical_stereo_excluded
            + counts.radical_no_expandable_hydrogens
            + counts.radical_compared
    );
    assert_eq!(
        counts.radical_compared,
        counts.radical_identical + counts.radical_both_rejected
    );
}

#[test]
fn n2s_radical_corpus_ordinary_hydrogen_variants_retain_complete_name_results() {
    let started = Instant::now();
    let parser = Parser::new().expect("embedded OPSIN resources");
    let budget = Budget::default();
    let reader = BufReader::new(GzDecoder::new(
        include_bytes!("../../nazuchem-n2s/tests/fixtures/structure-golden.jsonl.gz").as_slice(),
    ));
    let (mut counts, mut differences) = (ExplicitHydrogenCounts::default(), Vec::new());
    let (mut success, mut warning) = (0, 0);
    for line in reader.lines() {
        let record: Value =
            serde_json::from_str(&line.expect("fixture decompression")).expect("fixture JSON");
        if record["family"] != "radicals" || record["allow_radicals"] != true {
            continue;
        }
        counts.sources += 1;
        let name = record["input"].as_str().unwrap();
        let label = format!("n2s radical ordinary-H: {name}");
        let result = parser.parse(
            name,
            &ParseOptions {
                allow_radicals: true,
                ..ParseOptions::strict()
            },
        );
        let status = match result.status {
            Status::Success => {
                success += 1;
                "SUCCESS"
            }
            Status::Warning => {
                warning += 1;
                "WARNING"
            }
            Status::Failure => "FAILURE",
        };
        assert_eq!(status, record["status"].as_str().unwrap(), "{label}");
        let structure = result.structure.expect("radical corpus structure");
        let ir = structure.to_ir().expect("radical corpus IR");
        let text = structure
            .semantic_cxsmiles()
            .expect("radical corpus semantic CXSMILES");
        // This source is independent of the SMILES test bridge: the expanded
        // graph comes from Structure::to_ir, and the baseline names the original
        // semantic CXSMILES, without an IR -> SMILES serialization step.
        compare_ordinary_hydrogen_variant(
            &label,
            &text,
            &ir,
            &budget,
            &mut counts,
            &mut differences,
        );
    }
    finish(
        "n2s radical ordinary-H variants",
        &counts,
        differences,
        started,
    );
    counts.names.assert_complete();
    assert_eq!((counts.sources, success, warning), (12, 11, 1));
    assert_eq!(counts.converter_rejected, 0);
    assert_eq!(counts.radical_sources, 12);
    assert_eq!(
        (
            counts.stereo_excluded,
            counts.no_expandable_hydrogens,
            counts.explicit_hydrogen_vertices
        ),
        (0, 1, 32)
    );
    assert_eq!(
        (
            counts.names.records,
            counts.names.identical,
            counts.names.both_rejected
        ),
        (11, 7, 4)
    );
    counts.names.assert_routes(&[("radical ion", 1)], 10);
    assert_eq!(
        counts.sources,
        counts.stereo_excluded + counts.no_expandable_hydrogens + counts.names.records
    );
    assert_eq!(
        counts.radical_compared,
        counts.radical_identical + counts.radical_both_rejected
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
        // Issue #36: IR states ξ as unknown configuration (P-91.2.1.2.1),
        // which SMILES cannot express. Compare the rest with the baseline.
        Ok(mut ir) => {
            ir.stereo.retain(|stereo| {
                !matches!(
                    stereo,
                    nazuchem_ir::Stereo::Tetrahedral {
                        winding: nazuchem_ir::Winding::Unknown,
                        ..
                    } | nazuchem_ir::Stereo::DoubleBond {
                        relation: nazuchem_ir::Relation::Unknown,
                        ..
                    }
                )
            });
            ir
        }
        Err(error) => {
            counts.ir_unavailable += 1;
            difference(
                differences,
                format!("{label}: Structure::to_ir failed: {error}"),
            );
            return;
        }
    };
    let adapted = compared_name_molecule(&ir, budget, &mut counts.names);
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
    let started = Instant::now();
    let reader = BufReader::new(GzDecoder::new(
        include_bytes!("../../nazuchem-n2s/tests/fixtures/structure-golden.jsonl.gz").as_slice(),
    ));
    // Only the test harness uses threads. Every record still invokes to_ir and
    // both naming routes; separate parsers and thread-local fusion caches keep
    // requests independent. Eight bounded shards reduce the full oracle gate's
    // wall time, and joining in shard order keeps diagnostics deterministic.
    const SHARDS: usize = 8;
    let mut shards: [Vec<(usize, Value)>; SHARDS] = std::array::from_fn(|_| Vec::new());
    for (row, line) in reader.lines().enumerate() {
        let record =
            serde_json::from_str(&line.expect("fixture decompression")).expect("fixture JSON");
        shards[row % SHARDS].push((row + 1, record));
    }
    let results = std::thread::scope(|scope| {
        let workers: Vec<_> = shards
            .into_iter()
            .map(|rows| {
                scope.spawn(move || {
                    let shard_started = Instant::now();
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
                    (counts, differences, shard_started.elapsed())
                })
            })
            .collect();
        workers
            .into_iter()
            .map(|worker| worker.join().expect("structure corpus shard did not panic"))
            .collect::<Vec<_>>()
    });
    let (mut counts, mut differences) = (N2sCounts::default(), Vec::new());
    for (shard, (shard_counts, shard_differences, elapsed)) in results.into_iter().enumerate() {
        eprintln!(
            "n2s structure oracle shard {}/{}: records={}, comparisons={}, elapsed_seconds={:.3}",
            shard + 1,
            SHARDS,
            shard_counts.records,
            shard_counts.names.records,
            elapsed.as_secs_f64(),
        );
        counts.merge(shard_counts);
        for text in shard_differences {
            difference(&mut differences, text);
        }
    }
    finish("n2s structure oracle", &counts, differences, started);
    assert_n2s_coverage(&counts, (2624, 2535, 70, 19));
    assert_eq!(
        (counts.names.identical, counts.names.both_rejected),
        (2361, 244)
    );
    counts
        .names
        .assert_routes(&[("radical ion", 1), ("salt", 116)], 2457);
}

#[test]
fn every_n2s_configuration_oracle_structure_retains_complete_name_results() {
    let started = Instant::now();
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
    finish("n2s configuration oracle", &counts, differences, started);
    assert_n2s_coverage(&counts, (640, 304, 112, 224));
    assert_eq!(
        (counts.names.identical, counts.names.both_rejected),
        (336, 80)
    );
    counts.names.assert_routes(&[], 336);
}

#[test]
fn every_n2s_golden_structure_retains_complete_name_results() {
    let started = Instant::now();
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
    finish("n2s generated-name golden", &counts, differences, started);
    assert_eq!(records, 1650);
    // Two generated names cite impossible hydro/indicated hydrogen (P-14.7.1, P-31.2.2; #34).
    assert_n2s_coverage(&counts, (1659, 1648, 0, 11));
    assert_eq!(
        (counts.names.identical, counts.names.both_rejected),
        (1647, 1)
    );
    counts
        .names
        .assert_routes(&[("monatomic cation", 1), ("salt", 11)], 1636);
}
