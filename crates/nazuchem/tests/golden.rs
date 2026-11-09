// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | crates/nazuchem/tests/golden.rs | standalone library golden gate

#![cfg(all(feature = "n2s", feature = "s2n"))]

use std::collections::HashSet;
use std::fmt::Write;
use std::path::Path;

use nazuchem::{Budget, n2s, s2n};
use serde_json::{Value, json};

const CORPORA: [&str; 3] = [
    include_str!("../../../tests/fixtures/s2n/bakeoff_corpus.tsv"),
    include_str!("../../../tests/fixtures/s2n/openclatura_tests.tsv"),
    include_str!("../../../tests/fixtures/s2n/rules_cases.tsv"),
];

fn corpus_smiles() -> Vec<&'static str> {
    let mut seen = HashSet::new();
    CORPORA
        .iter()
        .flat_map(|corpus| corpus.lines().skip(1))
        .filter(|line| !line.trim().is_empty() && !line.starts_with('#'))
        .map(|line| line.split('\t').next().expect("SMILES column").trim())
        .filter(|smiles| !smiles.is_empty() && seen.insert(*smiles))
        .collect()
}

fn strict_parse(parser: &n2s::Parser, name: &str) -> Value {
    let parsed = parser.parse(name, &n2s::ParseOptions::strict());
    let status = match parsed.status {
        n2s::Status::Success => "SUCCESS",
        n2s::Status::Warning => "WARNING",
        n2s::Status::Failure => "FAILURE",
    };
    let warnings: Vec<_> = parsed
        .warnings
        .iter()
        .map(|warning| json!({"kind": warning.kind.as_str(), "message": warning.message}))
        .collect();
    let structure = parsed.structure.as_ref().map(|structure| {
        structure
            .semantic_cxsmiles_with_locants_and_axes()
            .map(|(cxsmiles, locants, axes)| {
                let axes: Vec<_> = axes
                    .iter()
                    .map(|axis| {
                        json!({
                            "atoms": axis.atoms,
                            "refs": axis.refs,
                            "clockwise": axis.clockwise,
                            "descriptor": axis.descriptor,
                        })
                    })
                    .collect();
                json!({"cxsmiles": cxsmiles, "locants": locants, "axes": axes})
            })
            .map_err(|error| error.to_string())
    });
    json!({
        "name": name,
        "status": status,
        "message": parsed.message,
        "warnings": warnings,
        "structure": structure,
    })
}

#[test]
fn s2n_n2s_golden() {
    let budget = Budget::default();
    let parser =
        n2s::Parser::with_budget(budget.clone()).expect("initialize embedded OPSIN parser");
    let mut actual = String::new();
    for smiles in corpus_smiles() {
        let generated = s2n::name_smiles(smiles, &budget);
        let mut parses = Vec::new();
        if let Ok(name) = &generated {
            let mut seen = HashSet::new();
            for text in std::iter::once(name.text.as_str())
                .chain(name.alternatives.iter().map(|name| name.name.as_str()))
                .chain(name.verification.iter().map(|target| target.name.as_str()))
            {
                if seen.insert(text) {
                    parses.push(strict_parse(&parser, text));
                }
            }
        }
        let record = json!({"smiles": smiles, "s2n": generated, "n2s": parses});
        writeln!(actual, "{record}").expect("append golden record");
    }

    let path =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/golden/s2n_n2s.jsonl");
    if std::env::var("NAZUCHEM_UPDATE_GOLDEN").as_deref() == Ok("1") {
        std::fs::create_dir_all(path.parent().expect("golden directory"))
            .expect("create golden directory");
        std::fs::write(&path, &actual).expect("write golden fixture");
    }
    let expected = std::fs::read_to_string(&path)
        .expect("read golden fixture; regenerate with NAZUCHEM_UPDATE_GOLDEN=1");
    if actual != expected {
        let line = actual
            .lines()
            .zip(expected.lines())
            .position(|(actual, expected)| actual != expected)
            .unwrap_or_else(|| actual.lines().count().min(expected.lines().count()));
        panic!(
            "golden differs at record {} (actual {} records, expected {}); \
             review changes before regenerating with NAZUCHEM_UPDATE_GOLDEN=1",
            line + 1,
            actual.lines().count(),
            expected.lines().count(),
        );
    }
}
