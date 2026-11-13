// SPDX-License-Identifier: MIT
// Origin: original
// Source: ReShiki | crates/opsin/tests/graph_resource_oracle.rs | 0e6e9009ca2e8bfebe149c260fc5a2f484f7dfb9
// Source: OPSIN | opsin-core/src/main/java/uk/ac/cam/ch/wwmm/opsin/SMILESFragmentBuilder.java | 42e87ac6efb1bd9f578986f5ea7ee1faa8e93c91
// Source: OPSIN | opsin-core/src/main/java/uk/ac/cam/ch/wwmm/opsin/FragmentTools.java | 42e87ac6efb1bd9f578986f5ea7ee1faa8e93c91
// Source: ReShiki | crates/opsin/tests/fixtures/graph-resource-inputs.tsv | 5f559fd482e237307a589dfcb7ea333cee0a6e13
// Source: ReShiki | crates/opsin/tests/fixtures/graph-resource-golden.tsv | 5f559fd482e237307a589dfcb7ea333cee0a6e13
// Copyright 2017 Daniel Lowe

//! Frozen dev-only Java oracle, OPSIN main MIT Daniel Lowe/contributors.
//! Source 42e87ac6efb1bd9f578986f5ea7ee1faa8e93c91; development jar and
//! fixture SHA256 identities are recorded in provenance/manifest.json.
//! Captured with package-local SMILESFragmentBuilder.build(value,type,labels),
//! then FragmentTools.convertSpareValenciesToDoubleBonds. Tests invoke no Java.
use nazuchem_n2s::{
    fragment_tools,
    graph::{FragmentId, Graph},
    smiles::build_fragment,
};
use std::collections::BTreeSet;
fn snapshot(graph: &Graph, fragment: FragmentId) -> String {
    let fragment = graph.fragment(fragment);
    let spare: String = fragment
        .atoms
        .iter()
        .map(|a| {
            if graph.atom(*a).spare_valency {
                '1'
            } else {
                '0'
            }
        })
        .collect();
    let orders: String = fragment
        .bonds
        .iter()
        .map(|b| char::from(b'0' + graph.bond(*b).order))
        .collect();
    let cycles: String = fragment
        .atoms
        .iter()
        .map(|a| if graph.atom(*a).in_cycle { '1' } else { '0' })
        .collect();
    let indicated = fragment
        .indicated_hydrogens
        .iter()
        .map(|a| {
            fragment
                .atoms
                .iter()
                .position(|other| other == a)
                .unwrap()
                .to_string()
        })
        .collect::<Vec<_>>()
        .join(",");
    format!("{spare}/{orders}/{cycles}/{indicated}")
}
#[test]
fn every_distinct_literal_resource_has_source_cycles_indicated_h_and_kekulization() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("resources/xml");
    // Native simpleSubstituents.xml acid-anion prefixes (issue #83, P-72.6.1;
    // issue #145, P-65.6.2.3.1 for selenium and tellurium)
    // postdate the frozen Java oracle; they are acyclic, so they must build
    // without spare valencies or cycles and kekulization must leave them as is.
    const NATIVE: [&str; 16] = [
        "-S(=S)[O-]",
        "-S(=[Se])[O-]",
        "-S(=[Te])[O-]",
        "-S(=O)(=S)[O-]",
        "-S(=O)(=[Se])[O-]",
        "-S(=O)(=[Te])[O-]",
        "-S(=S)(=S)[O-]",
        "-S(=[Se])(=[Se])[O-]",
        "-S(=[Te])(=[Te])[O-]",
        "-S(=S)(=[Se])[O-]",
        "-S(=[Se])(=[Te])[O-]",
        "-S(=S)(=[Te])[O-]",
        "-[Se](=O)(=O)[O-]",
        "-[Se](=O)[O-]",
        "-[Te](=O)(=O)[O-]",
        "-[Te](=O)[O-]",
    ];
    let mut native = Vec::new();
    let mut inputs = BTreeSet::new();
    let mut literals = 0;
    for entry in std::fs::read_dir(root).unwrap() {
        let path = entry.unwrap().path();
        if path.extension().is_none_or(|x| x != "xml") {
            continue;
        }
        let xml = std::fs::read_to_string(&path).unwrap();
        let document = roxmltree::Document::parse_with_options(
            &xml,
            roxmltree::ParsingOptions {
                allow_dtd: true,
                ..Default::default()
            },
        )
        .unwrap();
        for list in document
            .descendants()
            .filter(|n| n.has_tag_name("tokenList") && n.attribute("tagname") == Some("group"))
        {
            for token in list.children().filter(|n| n.has_tag_name("token")) {
                let Some(value) = token.attribute("value").or_else(|| list.attribute("value"))
                else {
                    continue;
                };
                if NATIVE.contains(&value) {
                    native.push(value.to_owned());
                    continue;
                }
                inputs.insert((
                    value.to_owned(),
                    token
                        .attribute("type")
                        .or_else(|| list.attribute("type"))
                        .unwrap_or("")
                        .to_owned(),
                    token
                        .attribute("labels")
                        .or_else(|| list.attribute("labels"))
                        .unwrap_or("none")
                        .to_owned(),
                ));
                literals += 1;
            }
        }
    }
    native.sort();
    let mut expected_native = NATIVE;
    expected_native.sort();
    assert_eq!(native, expected_native);
    for value in NATIVE {
        let mut graph = Graph::default();
        let fragment = build_fragment(&mut graph, value, "substituent", "none").unwrap();
        let before = snapshot(&graph, fragment);
        let [spare, _, cycles, indicated] = before.split('/').collect::<Vec<_>>()[..] else {
            panic!("{value}: {before}");
        };
        assert!(
            !spare.contains('1') && !cycles.contains('1') && indicated.is_empty(),
            "{value}: {before}"
        );
        fragment_tools::convert_spare_valencies_to_double_bonds(&mut graph, fragment).unwrap();
        assert_eq!(snapshot(&graph, fragment), before, "{value}");
    }
    assert_eq!(literals, 3183);
    assert_eq!(inputs.len(), 2643);
    let actual_inputs: String = inputs
        .iter()
        .map(|(value, kind, labels)| format!("{value}\t{kind}\t{labels}\n"))
        .collect();
    assert_eq!(
        actual_inputs,
        include_str!("fixtures/graph-resource-inputs.tsv"),
        "Sorted resource triples differ from the frozen oracle inputs"
    );
    let frozen = include_str!("fixtures/graph-resource-golden.tsv");
    assert_eq!(frozen.lines().count(), inputs.len());
    for (index, ((value, kind, labels), expected)) in
        inputs.into_iter().zip(frozen.lines()).enumerate()
    {
        let p: Vec<_> = expected.split('\t').collect();
        assert_eq!(p.len(), 3, "Oracle column count at row {index}");
        assert_eq!(p[0], index.to_string(), "Oracle row index");
        let mut graph = Graph::default();
        let fragment = build_fragment(&mut graph, &value, &kind, &labels).unwrap();
        assert_eq!(
            snapshot(&graph, fragment),
            p[1],
            "before conversion: {value} ({kind}, {labels})"
        );
        let before = graph.clone();
        let result = fragment_tools::convert_spare_valencies_to_double_bonds(&mut graph, fragment);
        if let Some(message) = p[2].strip_prefix("ERROR:") {
            assert_eq!(result.unwrap_err().0, message, "conversion error: {value}");
            assert_eq!(graph, before, "failed conversion must be transactional");
        } else {
            result.unwrap_or_else(|error| panic!("{value}: {error}"));
            assert_eq!(
                snapshot(&graph, fragment),
                p[2],
                "after conversion: {value} ({kind}, {labels})"
            );
        }
    }
}
