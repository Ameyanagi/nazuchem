// SPDX-License-Identifier: MIT AND (MIT OR Apache-2.0)
// Origin: mixed
// Source: OPSIN | opsin-inchi/src/test/resources/uk/ac/cam/ch/wwmm/opsin/phanes.txt | 42e87ac6efb1bd9f578986f5ea7ee1faa8e93c91
// Source: ReShiki | src/naming/openclatura/phane/tests.rs | 5f559fd482e237307a589dfcb7ea333cee0a6e13
// Copyright (c) 2026 lamalab-org
// Source: openclatura | tests/unit/test_analysis.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: Nazuchem | crates/nazuchem-s2n/src/phane/tests.rs | issue #4
// Copyright 2017 Daniel Lowe
// Copyright (c) 2026 Ameyanagi and ReShiki contributors
// Changes: Original Rust phane test harness with explicit OPSIN phane-oracle fixture selectors and assertions. Adapted openclatura control cases are marked below.
// Local contributions: MIT OR Apache-2.0; see NOTICE and PROVENANCE-nazuchem.toml.

use super::*;
use crate::{graph_io, test_support::verify_identity};
use nazuchem_ir::{NameReason, NameStatus};

const PARA: &str = "C1Cc2ccc(cc2)CCc2ccc1cc2";
const META: &str = "C1Cc2cccc(c2)CCc2cccc1c2";

// Provenance: fn round_trip | original
// Source: ReShiki | src/naming/openclatura/phane/tests.rs::round_trip
fn round_trip(smiles: &str, name: &str) -> Result<(), String> {
    let parsed = nazuchem_n2s::Parser::new()
        .map_err(|e| e.to_string())?
        .parse(name, &nazuchem_n2s::ParseOptions::strict());
    assert_eq!(
        parsed.status,
        nazuchem_n2s::Status::Success,
        "{name}: {}",
        parsed.message
    );
    assert!(parsed.warnings.is_empty(), "{name}: {:?}", parsed.warnings);
    let actual = parsed
        .structure
        .ok_or("Missing phane graph")?
        .semantic_cxsmiles()
        .map_err(|e| e.to_string())?;
    verify_identity(smiles, &actual)
}

#[test]
// Provenance: fn cyclic_phane_names_reconstruct_exactly | original
// Source: ReShiki | src/naming/openclatura/phane/tests.rs::cyclic_phane_names_reconstruct_exactly
fn cyclic_phane_names_reconstruct_exactly() -> Result<(), String> {
    for (smiles, expected) in [
        (PARA, "1,4(1,4)-dibenzenacyclohexaphane"),
        (
            "C1CC2=CC=C(C=C2)CCC2=CC=C1C=C2",
            "1,4(1,4)-dibenzenacyclohexaphane",
        ),
        (META, "1,4(1,3)-dibenzenacyclohexaphane"),
        (
            "C1CCc2cccc(c2)CCc2ccc1cc2",
            "1(1,3),4(1,4)-dibenzenacycloheptaphane",
        ),
        (
            "C1Cc2ncc(cc2)CCc2ccc1cc2",
            "1(2,5)-pyridina-4(1,4)-benzenacyclohexaphane",
        ),
        (
            "C1Cc2ccc(o2)CCc2ccc1cc2",
            "1(2,5)-furana-4(1,4)-benzenacyclohexaphane",
        ),
        (
            "C1Cc2ccc(s2)CCc2ccc1cc2",
            "1(2,5)-thiophena-4(1,4)-benzenacyclohexaphane",
        ),
        (
            "C1OC2=CC=C(C=C2)CCC=2C=CC1=NC2",
            "3-oxa-1(2,5)-pyridina-4(1,4)-benzenacyclohexaphane",
        ),
        (
            "C1Nc2ccc(cc2)CNc2ccc1cc2",
            "2,5-diaza-1,4(1,4)-dibenzenacyclohexaphane",
        ),
    ] {
        assert_eq!(
            name(&graph_io::read_smiles(smiles, &Default::default())?)
                .map(|p| p.0)
                .as_deref(),
            Some(expected),
            "{smiles}"
        );
        round_trip(smiles, expected)?;
        // Input order and aromatic/Kekule spelling cannot change the result.
        let canonical = crate::test_support::canonical_smiles(smiles)?;
        assert_eq!(
            name(&graph_io::read_smiles(&canonical, &Default::default())?)
                .map(|p| p.0)
                .as_deref(),
            Some(expected)
        );
    }
    Ok(())
}

#[test]
// Provenance: fn declines_out_of_scope_graphs | adapted
// Source: ReShiki | src/naming/openclatura/phane/tests.rs::declines_out_of_scope_graphs
// Source: openclatura | tests/unit/test_analysis.py::test_biphenyl_is_graph_selected_not_exact_name_rewritten
// Test provenance: upstream inputs/expected names with local assertions and extensions;
// exact matches are conservatively adapted; independent ancestry is not established.
fn declines_out_of_scope_graphs() -> Result<(), String> {
    for smiles in [
        "CC1Cc2ccc(cc2)CCc2ccc1cc2",
        "C1Cc2cc(C)cc(c2)CCc2cccc1c2",
        "C1=Cc2ccc(cc2)CCc2ccc1cc2",
        "C1CCCCCCCCCCC1",
        "C1CCCCCCc2ccc1cc2",
        "c1ccccc1-c1ccccc1",
        "[13CH2]1Cc2ccc(cc2)CCc2ccc1cc2",
        "C1[NH2+]c2ccc(cc2)CCc2ccc1cc2",
        "C[C@H]1Cc2ccc(cc2)CCc2ccc1cc2",
        "C1Cc2ccc3ccccc3c2CCc2ccc1cc2",
        "C1Cc2ccc(cc2)C2(CC2)Cc2ccc1cc2",
        "C1Cc2ccc(cc2)C2CC(C2)c2ccc1cc2",
        // Fewer than six nodes is not a cyclophane (P-52.2.5.1).
        "C1c2ccccc2Cc2ccccc21",
        "C1c2cccc(c2)Cc2cccc1c2",
        // Adjacent attachment atoms make a fused system (P-52.2.5.2.1).
        "C1Cc2ccccc2CCc2ccccc21",
        "C1Cc2occc2CCc2ccc1cc2",
    ] {
        assert!(
            name(&graph_io::read_smiles(smiles, &Default::default())?).is_none(),
            "{smiles}"
        );
    }
    Ok(())
}

#[test]
// Provenance: fn mixed_attachment_sets_follow_amplificant_seniority | original
// Source: ReShiki | src/naming/openclatura/phane/tests.rs::mixed_attachment_sets_follow_amplificant_seniority
fn mixed_attachment_sets_follow_amplificant_seniority() -> Result<(), String> {
    // Pyridine/benzene/pyridine/pyridine/benzene/pyridine, each separated by CH2.
    // PhI-3.2.2 compares the pyridine attachment sets before the benzene sets.
    let smiles = "C1c2cc(ccn2)Cc2cccc(c2)Cc2ccc(cn2)Cc2cc(ccn2)Cc2cccc(c2)Cc2cc1ccn2";
    let expected = "1(2,4),3,7(4,2),9(5,2)-tetrapyridina-5,11(1,3)-dibenzenacyclododecaphane";
    assert_eq!(
        name(&graph_io::read_smiles(smiles, &Default::default())?)
            .map(|p| p.0)
            .as_deref(),
        Some(expected)
    );
    round_trip(smiles, expected)
}

#[test]
// Provenance: fn ordinary_names_with_phane_alternatives_are_never_preferred | original
// Source: ReShiki | src/naming/openclatura/phane/tests.rs::ordinary_names_with_phane_alternatives_are_never_preferred
fn ordinary_names_with_phane_alternatives_are_never_preferred() -> Result<(), String> {
    // PARA has a von Baeyer name now that its two-atom secondary bridge is
    // numbered from the higher bridgehead (P-23.2.5.2).
    for smiles in [META, PARA] {
        let named = super::super::name_smiles(smiles, &Default::default())?;
        assert_eq!(named.status, NameStatus::NotPreferred(NameReason::Phane));
        assert_eq!(named.alternatives.len(), 1);
        assert!(!named.text.contains("phane"));
        round_trip(smiles, &named.text)?;
        for alternative in named.alternatives {
            assert_eq!(
                alternative.status,
                NameStatus::Unverified(NameReason::Phane)
            );
            round_trip(smiles, &alternative.name)?;
        }
    }
    Ok(())
}

#[test]
// Provenance: fn linear_phanes_number_from_the_end_with_the_lower_locants | original
// Source: Nazuchem | crates/nazuchem-s2n/src/phane/tests.rs::linear_phanes_number_from_the_end_with_the_lower_locants
fn linear_phanes_number_from_the_end_with_the_lower_locants() -> Result<(), String> {
    for (smiles, expected) in [
        // P-52.2.5.3 Example 2.
        (
            "c1ccccc1Oc2ccc(Oc3ccc(Oc4ccccc4)cc3)cc2",
            "2,4,6-trioxa-1,7(1),3,5(1,4)-tetrabenzenaheptaphane",
        ),
        // P-51.3.2.2.
        (
            "c1ccccc1Oc1cccc(c1)Oc1cccc(c1)Oc1ccccc1",
            "2,4,6-trioxa-1,7(1),3,5(1,3)-tetrabenzenaheptaphane",
        ),
        // P-52.2.5.3 Example 1.
        (
            "c1ccccc1-c1cccc(c1)-c1cccc(c1)-c1cc(cnc1)-c1cccc(c1)-c1cccc(c1)-c1ccccc1",
            "4(3,5)-pyridina-1,7(1),2,3,5,6(1,3)-hexabenzenaheptaphane",
        ),
        // P-26.4.1.3: the senior amplificant takes the lower superatom locant.
        (
            "c1ccccc1Cc1ccc(cn1)Cc1ccc(cc1)Cc1ccccc1",
            "3(2,5)-pyridina-1,7(1),5(1,4)-tribenzenaheptaphane",
        ),
        // P-26.4.1.4: the lower attachment locant is next to the lower
        // skeletal locant.
        (
            "c1ccccc1Cc1ccc(o1)Cc1ccc(cc1)Cc1ccccc1",
            "3(2,5)-furana-1,7(1),5(1,4)-tribenzenaheptaphane",
        ),
        // P-26.4.2.2: the lower attachment set takes the lower superatom locant.
        (
            "c1ccccc1Cc1ccc(cn1)Cc1cccc(n1)Cc1ccccc1",
            "3(2,5),5(2,6)-dipyridina-1,7(1)-dibenzenaheptaphane",
        ),
        // P-26.4.2.4, P-26.5.1: the attachment locants in name order are
        // settled before the heteroatom takes its locant.
        (
            "c1ccccc1Cc1ccc(Cc2ccc(Oc3ccccc3)cn2)cn1",
            "6-oxa-3,5(2,5)-dipyridina-1,7(1)-dibenzenaheptaphane",
        ),
        (
            "c1ccccc1-c1ccc(-c2ccc(Oc3cccc(-c4cccc(-c5ccccn5)n4)n3)cc2)cc1",
            "4-oxa-1(2),2,3(2,6)-tripyridina-5,6(1,4),7(1)-tribenzenaheptaphane",
        ),
        (
            "c1ccccc1-c2ccc(-c3ccc(-c4ccc(-c5ccc(-c6ccc(-c7ccccc7)cc6)cc5)cc4)cc3)cc2",
            "1,7(1),2,3,4,5,6(1,4)-heptabenzenaheptaphane",
        ),
    ] {
        assert_eq!(
            name(&graph_io::read_smiles(smiles, &Default::default())?),
            Some((expected.to_owned(), true)),
            "{smiles}"
        );
        round_trip(smiles, expected)?;
        let named = crate::name_smiles(smiles, &Default::default())?;
        assert_eq!(named.text, expected, "{smiles}");
        assert_eq!(named.status, NameStatus::Preferred, "{smiles}");
        let [alternative] = &named.alternatives[..] else {
            panic!("{smiles}: {:?}", named.alternatives);
        };
        assert_eq!(
            alternative.status,
            NameStatus::NotPreferred(NameReason::LinearPhane),
            "{smiles}"
        );
        verify_identity(
            smiles,
            &crate::test_support::parse_strict(&alternative.name)?,
        )?;
    }
    // P-52.2.5.1 (2): four rings and seven nodes. Substituted nodes, amine
    // nodes and skeletons past nine nodes are not built.
    for smiles in [
        "c1ccccc1Oc2ccc(Oc3ccccc3)cc2",
        "c1ccccc1-c1cc(cnc1)-c1cccc(c1)-c1cccc(c1)-c1ccccc1",
        "c1ccccc1Oc1ccc(cc1)C(C)(C)c1ccc(Oc2ccccc2)cc1",
        "Oc1ccc(Oc2ccc(-c3ccc(CCc4ccccc4)cc3)cc2)cc1",
        "c1ccccc1SCc1ccc(cc1)Nc1ccc(cc1)CCc1ccccc1",
        "c1cocc1SCCSc1cocc1SCCSc1c(SCCSc2ccoc2)coc1",
    ] {
        assert!(
            name(&graph_io::read_smiles(smiles, &Default::default())?).is_none(),
            "{smiles}"
        );
    }
    Ok(())
}
