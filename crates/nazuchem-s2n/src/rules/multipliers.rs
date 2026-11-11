// SPDX-License-Identifier: MIT AND (MIT OR Apache-2.0)
// Origin: mixed
// Source: openclatura | src/openclatura/rules/multipliers.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: ReShiki | src/naming/openclatura/rules/multipliers.rs | 5f559fd482e237307a589dfcb7ea333cee0a6e13
// Source: openclatura | src/openclatura/principal_suffixes.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | src/openclatura/rules/elision.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/data/parser_xml_resources/suffixRules.json | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/unit/test_analysis.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/data/parser_xml_resources/suffixes.json | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/data/parser_xml_resources/carboxylicAcids.json | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/unit/test_spiro_composed_parent_families.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: Nazuchem | crates/nazuchem-s2n/src/rules/multipliers.rs | issue #21
// Copyright (c) 2026 lamalab-org
// Copyright (c) 2026 Ameyanagi and ReShiki contributors
// Changes: Inherited openclatura numerical-term constants/generation plus original suffix-specific elision and compound-anion multiplication.
// Local contributions: MIT OR Apache-2.0; see NOTICE and PROVENANCE-nazuchem.toml.

//! Shared P-14.2 numerical terms (counts through 1000).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Multiplier {
    pub count: usize,
    pub basic: String,
    pub complex: String,
}
const UNITS: [&str; 10] = [
    "", "hen", "do", "tri", "tetra", "penta", "hexa", "hepta", "octa", "nona",
];
const TENS: [&str; 10] = [
    "",
    "deca",
    "icosa",
    "triaconta",
    "tetraconta",
    "pentaconta",
    "hexaconta",
    "heptaconta",
    "octaconta",
    "nonaconta",
];
const HUNDREDS: [&str; 10] = [
    "", "hecta", "dicta", "tricta", "tetracta", "pentacta", "hexacta", "heptacta", "octacta",
    "nonacta",
];
// Provenance: fn numerical_term | adapted
// Source: openclatura | src/openclatura/rules/multipliers.py::numerical_term
pub fn numerical_term(count: usize) -> String {
    assert!((1..=1000).contains(&count));
    if count == 1000 {
        return "kilia".into();
    }
    let remainder = count % 100;
    let under = if remainder == 11 {
        "undeca".into()
    } else {
        let u = UNITS[remainder % 10];
        let t = TENS[remainder / 10];
        let t = if !u.is_empty() && remainder / 10 == 2 && u.ends_with(['a', 'e', 'i', 'o', 'u']) {
            &t[1..]
        } else {
            t
        };
        format!("{u}{t}")
    };
    format!("{under}{}", HUNDREDS[count / 100])
}
// Provenance: fn get | adapted
// Source: openclatura | src/openclatura/rules/multipliers.py::get
pub fn get(count: usize) -> Multiplier {
    assert!((2..=1000).contains(&count));
    let basic = match count {
        2 => "di".into(),
        3 => "tri".into(),
        _ => numerical_term(count),
    };
    let complex = match count {
        2 => "bis".into(),
        3 => "tris".into(),
        _ => format!("{}akis", basic.trim_end_matches('a')),
    };
    Multiplier {
        count,
        basic,
        complex,
    }
}
// Provenance: fn basic | adapted
// Source: openclatura | src/openclatura/rules/multipliers.py::basic
pub fn basic(count: usize) -> String {
    get(count).basic
}
/// Render a multiplied suffix with the elision required by P-16.7.1(c).
/// Only suffixes beginning with `a` or `o` elide a numerical prefix's terminal
/// `a`; substituent prefixes and other suffix initials keep it.
/// https://iupac.qmul.ac.uk/BlueBook/P1.html#160701
// Provenance: fn suffix | adapted
// Source: openclatura | src/openclatura/principal_suffixes.py::render_principal_suffix
// Source: openclatura | src/openclatura/rules/elision.py::elide_terminal_a
pub fn suffix(count: usize, suffix: &str) -> String {
    if count <= 1 {
        return suffix.into();
    }
    // P-72.2.2.2.2: these compound anionic suffixes take bis/tris, etc.
    if matches!(suffix, "olate" | "thiolate" | "selenolate" | "tellurolate") {
        return format!("{}({suffix})", complex_(count));
    }
    let mut multiplier = basic(count);
    if suffix.starts_with(['a', 'o']) && multiplier.ends_with('a') {
        multiplier.pop();
    }
    format!("{multiplier}{suffix}")
}
// Provenance: fn complex_ | adapted
// Source: openclatura | src/openclatura/rules/multipliers.py::complex_
pub fn complex_(count: usize) -> String {
    get(count).complex
}
// Provenance: fn candidate_splits | adapted
// Source: openclatura | src/openclatura/rules/multipliers.py::candidate_splits
pub fn candidate_splits(name: &str) -> Vec<(usize, String)> {
    let mut found = Vec::new();
    for n in 2..=1000 {
        let m = get(n);
        for p in [m.basic, m.complex] {
            if let Some(rest) = name.strip_prefix(&p) {
                found.push((p.len(), n, rest.to_owned()));
            }
        }
    }
    found.sort_by_key(|(len, _, _)| std::cmp::Reverse(*len));
    found.into_iter().map(|(_, n, s)| (n, s)).collect()
}

#[cfg(test)]
// Provenance: mod tests | original
// Source: ReShiki | src/naming/openclatura/rules/multipliers.rs::tests
mod tests {
    use super::*;

    #[test]
    // Provenance: fn tests::suffix_multiplier_elides_a_only_before_a_or_o | adapted
    // Source: ReShiki | src/naming/openclatura/rules/multipliers.rs::tests::suffix_multiplier_elides_a_only_before_a_or_o
    // Source: openclatura | tests/data/parser_xml_resources/suffixRules.json::<fixture data>
    // Source: openclatura | tests/unit/test_analysis.py::test_suffix_unsaturation_and_charge_operations_emit_precise_tokens
    // Source: openclatura | tests/unit/test_analysis.py::test_repeated_grouped_locants_use_following_structural_scope
    // Source: openclatura | tests/unit/test_analysis.py::test_principal_suffixes_render_programmatically_from_template_positions
    // Source: openclatura | tests/data/parser_xml_resources/suffixes.json::<fixture data>
    // Source: openclatura | tests/data/parser_xml_resources/carboxylicAcids.json::<fixture data>
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn suffix_multiplier_elides_a_only_before_a_or_o() {
        for (ending, expected) in [
            ("amine", "tetramine"),
            ("one", "tetrone"),
            ("ol", "tetrol"),
            ("oic acid", "tetroic acid"),
            ("ene", "tetraene"),
            ("imine", "tetraimine"),
            ("uide", "tetrauide"),
            ("yl", "tetrayl"),
            ("carboxylic acid", "tetracarboxylic acid"),
        ] {
            assert_eq!(suffix(4, ending), expected, "{ending}");
        }
        assert_eq!(basic(4), "tetra");
    }

    #[test]
    // Provenance: fn tests::suffix_multiplier_handles_single_and_larger_counts | adapted
    // Source: ReShiki | src/naming/openclatura/rules/multipliers.rs::tests::suffix_multiplier_handles_single_and_larger_counts
    // Source: openclatura | tests/unit/test_spiro_composed_parent_families.py::test_side_suffixes_group_same_kind_without_losing_component_primes
    // Source: openclatura | tests/unit/test_analysis.py::test_name_assembly_binds_every_final_token_to_graph_metadata
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn suffix_multiplier_handles_single_and_larger_counts() {
        for (count, ending, expected) in [
            (1, "one", "one"),
            (2, "one", "dione"),
            (3, "one", "trione"),
            (6, "ol", "hexol"),
            (11, "one", "undecone"),
            (20, "one", "icosone"),
            (1000, "one", "kilione"),
        ] {
            assert_eq!(suffix(count, ending), expected, "{count} {ending}");
        }
    }

    #[test]
    // Provenance: fn tests::parser_reads_every_generated_multiplier | original
    // Source: Nazuchem | crates/nazuchem-s2n/src/rules/multipliers.rs::tests::parser_reads_every_generated_multiplier
    // Issue #21. P-14.2.1 basic and P-14.2.2 "kis" terms through 1000, and von Baeyer
    // ring counts (P-23.1.9) above the tabulated prefixes.
    // https://iupac.qmul.ac.uk/BlueBook/P1.html#140201
    // https://iupac.qmul.ac.uk/BlueBook/P1.html#140202
    // https://iupac.qmul.ac.uk/BlueBook/P2.html#230109
    fn parser_reads_every_generated_multiplier() {
        let parser = nazuchem_n2s::Parser::new().unwrap();
        let options = nazuchem_n2s::ParseOptions::default();
        let reads = |name: &str, tag: &str, count: usize| {
            parser.parse_trees(name, &options).is_ok_and(|trees| {
                trees.iter().any(|tree| {
                    tree.arena
                        .descendants_named(tree.root, tag)
                        .into_iter()
                        .any(|id| tree.arena[id].attribute("value") == Some(&count.to_string()))
                })
            })
        };
        for count in 2..=1000 {
            let name = format!("{}fluoroethane", basic(count));
            assert!(reads(&name, "multiplier", count), "{name}");
            let name = format!("{}(chloromethyl)ethane", complex_(count));
            assert!(reads(&name, "multiplier", count), "{name}");
            if count > 10 {
                let name = format!(
                    "{}[2.2.1]heptane",
                    crate::ring_renderer::von_baeyer_kind(count)
                );
                assert!(reads(&name, "multiplier", count), "{name}");
            }
        }
    }

    #[test]
    // Provenance: fn tests::generated_large_multiplier_names_round_trip | original
    // Source: Nazuchem | crates/nazuchem-s2n/src/rules/multipliers.rs::tests::generated_large_multiplier_names_round_trip
    // Issue #21: hexaconta = 60 and hecta = 100 (P-14.2.1), hentriacontakis = 31
    // (P-14.2.2) and icosacyclo = 20 rings (P-23.1.9).
    // https://iupac.qmul.ac.uk/BlueBook/P1.html#140201
    // https://iupac.qmul.ac.uk/BlueBook/P1.html#140202
    // https://iupac.qmul.ac.uk/BlueBook/P2.html#230109
    fn generated_large_multiplier_names_round_trip() {
        let perfluoro = |ring: usize| format!("FC1(F){}C1(F)F", "C(F)(F)".repeat(ring - 2));
        // A saturated ladder of 20 fused cyclobutanes.
        let top: String = (10..30).map(|i| format!("C%{i}")).collect();
        let bottom: String = (10..30).rev().map(|i| format!("C%{i}")).collect();
        for (smiles, prefix) in [
            (perfluoro(30), "hexacontafluorocyclotriacontane"),
            (perfluoro(50), "hectafluorocyclopentacontane"),
            (
                format!("C1(CCCl){}1", "C(CCCl)".repeat(30)),
                "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31-hentriacontakis(2-chloroethyl)cyclohentriacontane",
            ),
            (format!("{top}CC{bottom}"), "icosacyclo[20.20.0."),
        ] {
            let name = crate::name_smiles(&smiles, &Default::default()).unwrap();
            assert!(name.text.starts_with(prefix), "{}", name.text);
            let parsed = crate::test_support::parse_strict(&name.text).unwrap();
            assert_eq!(
                crate::test_support::canonical_smiles(&parsed).unwrap(),
                crate::test_support::canonical_smiles(&smiles).unwrap()
            );
        }
    }
}
