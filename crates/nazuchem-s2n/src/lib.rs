// SPDX-License-Identifier: MIT AND (MIT OR Apache-2.0)
// Origin: mixed
// Source: Nazuchem | crates/nazuchem-s2n/src/lib.rs | Phase B step 8b
// Source: openclatura | src/openclatura/namer.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | src/openclatura/engine.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: ReShiki | src/naming/openclatura/mod.rs | 5f559fd482e237307a589dfcb7ea333cee0a6e13
// Source: openclatura | src/openclatura/component_namer.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | src/openclatura/locant_elision.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | src/openclatura/fusion/charges.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/unit/test_analysis.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/roundtrip/test_stereo.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/fixtures/diverse_corpus.csv | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/fixtures/diverse_corpus.golden.json | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/unit/test_chalcogen_functional_groups.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/unit/test_fused_bonds.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/roundtrip/test_replacement_parents.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/data/parser_xml_resources/simpleGroups.json | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/unit/test_polyoxide_chains.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/unit/test_audit_reconstruction.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/unit/test_all_chalcogen_acyl_regressions.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/roundtrip/test_oxyacids.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/unit/test_locant_elision.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/roundtrip/test_polycycles.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/roundtrip/test_role_models.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/roundtrip/test_principal_group_breadth.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/exploratory/test_group14_ether_ligand_parent_selection.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/data/parser_xml_resources/arylGroups.json | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/data/parser_xml_resources/nonCarboxylicAcids.json | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/unit/test_chalcogen_citation_compatibility.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/data/parser_xml_resources/functionalTerms.json | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/data/fusion_pin_cases.json | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/roundtrip/test_alkenes_alkynes.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/unit/test_chalcogen_structural_roles.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/data/parser_xml_resources/multipliers.json | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/unit/test_von_baeyer_candidate_audit.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/roundtrip/test_retained_fused_production.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/roundtrip/test_fused_parents.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/unit/test_retained_bridge_h_relocation.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/unit/test_rules_tables.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/integration/test_systematic_fusion.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/unit/test_fusion_wrappers.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/roundtrip/test_p13_fusion.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/roundtrip/test_fused_hydrocarbon_retained.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Copyright (c) 2026 lamalab-org
// Copyright (c) 2026 Ameyanagi and ReShiki contributors
// Changes: Port entry sequence surrounded by local naming routes, assessment, result types and alternatives.
// Local contributions: MIT OR Apache-2.0; see NOTICE and PROVENANCE-nazuchem.toml.

//! Plain-name pipeline that began as a port of openclatura 0.4.3, extended in ReShiki.
//! Attribution and license terms are in `licenses/openclatura/`.
//!
//! Supports chains, monocycles, retained and systematic fused parents, pure
//! spiro and von Baeyer ring systems, with independent reconstruction proofs.
//!
//! Callers supply a [`nazuchem_kernel::budget::Budget`]; `Budget::default()` is
//! unbounded. The library performs no I/O. Worst-case time on adversarial input
//! is unbounded, so hosts needing hard limits should supervise requests externally.
// Direct graph/path indexing preserves upstream invariants and panic behavior;
// JSON Value indexing preserves null for missing optional fields.
#![cfg_attr(not(test), allow(clippy::indexing_slicing))]
#![forbid(unsafe_code)]
#![cfg_attr(
    not(test),
    deny(
        clippy::unwrap_used,
        clippy::expect_used,
        clippy::panic,
        clippy::unreachable,
        clippy::todo,
        clippy::unimplemented
    )
)]

mod identity;
mod ir_input;
#[cfg(test)]
mod ir_proof_tests;
pub use identity::{canonical_smiles, verify_alternative_identity, verify_identity};

mod additive;
mod amine_oxides;
mod assembler;
mod assembly_charge;
mod assembly_parent;
mod assembly_parts;
mod assembly_prefixes;
mod assembly_spiro;
mod assembly_utils;
mod chains;
mod component_group_rules;
mod component_modifiers;
mod component_namer;
mod endo_exo_stereo;
mod formatting;
mod functional_prefixes;
mod fused_ion_templates;
mod fusion;
mod graph_io;
mod graph_kernel;
mod heteroatom_subgraphs;
mod heterospiro;
#[cfg(test)]
mod heterospiro_tests;
mod identity_validation;
mod ionic_naming;
mod isotope_naming;
mod radical_ions;
pub use radical_ions::recover as recover_radical_ion;
mod locant_elision;
mod locants;
mod molecule;
mod multiplicative;
mod namer;
mod nomenclature;
mod numbering;
mod parent_pipeline;
mod parent_selection;
mod perception;
mod phane;
mod polycycle_topology;
mod principal_groups;
mod principal_suffixes;
mod python_int_set;
mod radical_naming;
mod retained_derivative_hydrogen;
mod retained_fused_production;
mod retained_fused_templates;
mod retained_graph_model;
mod retained_name_policy;
mod retained_specs;
mod ring_assembly;
mod ring_parent;
mod ring_renderer;
mod rules;
mod small_ring_stereo;
mod special_cases;
mod species;
mod spelling;
mod spiro_subgraph;
mod spirobi;
#[cfg(test)]
mod spirobi_tests;
mod subgraph_tools;
mod unknown_stereo;
mod von_baeyer;

#[cfg(test)]
mod amine_oxide_tests;
#[cfg(test)]
mod anion_parent_tests;
#[cfg(test)]
mod anion_spelling_tests;
#[cfg(test)]
mod ion_spelling_tests;
#[cfg(test)]
mod salt_anion_spelling_tests;

/// Engine identifier for callers that pin naming output.
pub const PROFILE: &str = "openclatura 0.4.3 Rust port";

/// A generated name and how it relates to the preferred IUPAC name.
#[derive(Debug, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Name {
    /// The generated chemical name.
    pub text: String,
    /// Whether the name is preferred, with a reason when it is not.
    pub status: nazuchem_ir::NameStatus,
    /// Other generated names, each with its own preference status.
    pub alternatives: Vec<nazuchem_ir::AlternativeName>,
    /// Locants of unknown stereochemistry omitted from the main name.
    pub omitted_stereo: Vec<String>,
    /// Number of unknown stereochemical units omitted from the main name.
    pub omitted_units: usize,
    /// A recovered neutral target for verifying a descriptive radical-ion name.
    #[serde(default)]
    pub verification: Option<Verification>,
}

/// Recovered neutral target used to verify a descriptive radical-ion name.
#[derive(Debug, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Verification {
    /// The generated name of the recovered neutral structure.
    pub name: String,
    /// SMILES of the recovered neutral structure to compare with the parsed name.
    pub smiles: String,
}

/// Generate a name from molecular SMILES using IUPAC isotope notation.
///
/// Molecular perception uses the caller's budget. Returns an error when the
/// input is unsupported, cannot be named, or exceeds the supplied budget.
// Provenance: fn name_smiles | adapted
// Source: openclatura | src/openclatura/namer.py::name_smiles
// Source: openclatura | src/openclatura/engine.py::Namer.name_smiles
pub fn name_smiles(smiles: &str, budget: &nazuchem_kernel::budget::Budget) -> Result<Name, String> {
    name_smiles_with_style(smiles, nazuchem_ir::IsotopeStyle::Iupac, budget)
}

/// Generate a name from molecular SMILES using the requested isotope notation.
///
/// Molecular perception uses the caller's budget. Returns an error when the
/// input is unsupported, cannot be named, or exceeds the supplied budget.
// Provenance: fn name_smiles_with_style | original
// Source: ReShiki | src/naming/openclatura/mod.rs::name_smiles_with_style
pub fn name_smiles_with_style(
    smiles: &str,
    style: nazuchem_ir::IsotopeStyle,
    budget: &nazuchem_kernel::budget::Budget,
) -> Result<Name, String> {
    finish_name(name_routes(smiles, style, budget)?)
}

// Provenance: fn finish_name | original
// Source: Nazuchem | crates/nazuchem-s2n/src/lib.rs::finish_name
fn finish_name(mut name: Name) -> Result<Name, String> {
    // A radical-ion name is verified through its neutral parent only.
    if name.verification.is_some() {
        return Ok(name);
    }
    // P-93.5.7.1: the general name cites each M/P axis as Ra/Sa.
    let general: Vec<_> = std::iter::once(&name.text)
        .chain(
            name.alternatives
                .iter()
                .map(|alternative| &alternative.name),
        )
        .filter_map(|text| general_axial_name(text))
        .collect();
    for text in general {
        if !name
            .alternatives
            .iter()
            .any(|existing| existing.name == text)
        {
            name.alternatives.push(nazuchem_ir::AlternativeName {
                name: text,
                status: nazuchem_ir::NameStatus::NotPreferred(
                    nazuchem_ir::NameReason::AxialAlternative,
                ),
            });
        }
    }
    Ok(name)
}

/// Generate a name from a validated IR graph.
///
/// Uses the same chemical admission policy and naming rules as [`name_smiles`].
/// Unsupported annotations and stereo return an error rather than being lost.
/// Admission and ordinary component naming preserve surviving caller atom IDs.
/// Identity checks serialize canonical SMILES; the monatomic-cation probe and
/// the unknown-stereo baseline, salt, cation and radical-ion naming routes can
/// reparse it into detached graphs. The caller's IR is never changed.
/// Molecular perception uses the caller's budget.
// Provenance: fn name_molecule | original
// Source: Nazuchem | crates/nazuchem-s2n/src/lib.rs::name_molecule
pub fn name_molecule(
    molecule: &nazuchem_ir::Molecule,
    style: nazuchem_ir::IsotopeStyle,
    budget: &nazuchem_kernel::budget::Budget,
) -> Result<Name, String> {
    let state = ir_input::prepare(molecule, budget)?;
    let canonical = identity::canonical_state(&state, budget)?;
    let name = if unknown_stereo::has_unknown(&state) {
        #[cfg(test)]
        ir_proof_tests::record_fallback("unknown stereo");
        let normal = name_smiles_with_style(&canonical, style, budget)?;
        // The legacy unknown route admits its unspecified baseline first, then
        // requires one connected graph before rendering explicit descriptors.
        if canonical.contains('.') {
            return Err("Openclatura requires a single connected molecule".into());
        }
        let mol = graph_io::read_state(state, Some(molecule))?;
        unknown_stereo::name_graph(normal, &mol)?
    } else if canonical.contains('.') {
        #[cfg(test)]
        ir_proof_tests::record_fallback("salt");
        species::name_adduct(&canonical, style, budget)?
    } else if let Some(name) = {
        #[cfg(test)]
        ir_proof_tests::record_monatomic_probe();
        species::monatomic_cation(&canonical, false, budget)?
    } {
        #[cfg(test)]
        ir_proof_tests::record_fallback("monatomic cation");
        name
    } else if let Some((neutral, suffix)) = radical_ions::recover_state(&state, &canonical, budget)?
    {
        #[cfg(test)]
        ir_proof_tests::record_fallback("radical ion");
        name_radical_ion(neutral, suffix, style, budget)?
    } else {
        name_graph(graph_io::read_state(state, Some(molecule))?, style)?
    };
    finish_name(name)
}

/// The name with Ra/Sa for M/P in every locanted stereodescriptor group, such
/// as (1M,2R) -> (1Ra,2R); descriptor groups are innermost parentheses.
// Provenance: fn general_axial_name | original
// Source: ReShiki | src/naming/openclatura/mod.rs::general_axial_name
fn general_axial_name(text: &str) -> Option<String> {
    let mut out = String::with_capacity(text.len() + 4);
    let mut rest = text;
    while let Some(open) = rest.find('(') {
        out.push_str(&rest[..=open]);
        rest = &rest[open + 1..];
        // A stereodescriptor group is followed by a hyphen; a nuclide group
        // such as (31P) is not.
        let Some(close) = rest
            .find(['(', ')'])
            .filter(|&i| rest[i..].starts_with(")-"))
        else {
            continue;
        };
        let terms: Option<Vec<_>> = rest[..close]
            .split(',')
            .map(|term| {
                let descriptor = term.chars().last()?;
                let locant = &term[..term.len() - descriptor.len_utf8()];
                (locant.starts_with(|c: char| c.is_ascii_digit())
                    && "RSrsEZMPΞξ".contains(descriptor))
                .then_some((locant, descriptor))
            })
            .collect();
        if let Some(terms) = terms.filter(|terms| terms.iter().any(|(_, d)| matches!(d, 'M' | 'P')))
        {
            let terms: Vec<_> = terms
                .into_iter()
                .map(|(locant, descriptor)| match descriptor {
                    'M' => format!("{locant}Ra"),
                    'P' => format!("{locant}Sa"),
                    _ => format!("{locant}{descriptor}"),
                })
                .collect();
            out.push_str(&terms.join(","));
            rest = &rest[close..];
        }
    }
    out.push_str(rest);
    (out != text).then_some(out)
}

// Provenance: fn name_routes | adapted
// Source: openclatura | src/openclatura/namer.py::name_smiles
// Source: openclatura | src/openclatura/namer.py::name_component
fn name_routes(
    smiles: &str,
    style: nazuchem_ir::IsotopeStyle,
    budget: &nazuchem_kernel::budget::Budget,
) -> Result<Name, String> {
    if let Some(name) = unknown_stereo::name_smiles(smiles, style, budget)? {
        return Ok(name);
    }
    let canonical = canonical_smiles(smiles, budget)?;
    if canonical.contains('.') {
        return species::name_adduct(&canonical, style, budget);
    }
    if let Some(name) = species::monatomic_cation(&canonical, false, budget)? {
        return Ok(name);
    }
    if let Some((neutral, suffix)) = radical_ions::recover(&canonical, budget)? {
        return name_radical_ion(neutral, suffix, style, budget);
    }
    let mol = graph_io::read_smiles(smiles, budget)?;
    name_graph(mol, style)
}

// Provenance: fn name_radical_ion | original
// Source: Nazuchem | crates/nazuchem-s2n/src/lib.rs::name_radical_ion
// Source: ReShiki | src/naming/openclatura/mod.rs::name_routes
fn name_radical_ion(
    neutral: String,
    suffix: &str,
    style: nazuchem_ir::IsotopeStyle,
    budget: &nazuchem_kernel::budget::Budget,
) -> Result<Name, String> {
    let name = name_smiles_with_style(&neutral, style, budget)?;
    Ok(Name {
        text: format!("{} {suffix}", name.text),
        status: nazuchem_ir::NameStatus::NotPreferred(nazuchem_ir::NameReason::RadicalIon),
        alternatives: vec![],
        omitted_stereo: vec![],
        omitted_units: 0,
        verification: Some(Verification {
            name: name.text,
            smiles: neutral,
        }),
    })
}

// Provenance: fn name_graph | adapted
// Source: openclatura | src/openclatura/namer.py::name_component
fn name_graph(mol: molecule::Molecule, style: nazuchem_ir::IsotopeStyle) -> Result<Name, String> {
    if mol.atoms.values().any(|atom| atom.isotope.is_some()) {
        return isotope_naming::name_molecule(&mol, style);
    }
    // Match upstream require_validated_chalcogens before any shortcut can consume Po.
    if mol.atoms.values().any(|atom| atom.symbol == "Po") {
        return Err("openclatura: polonium nomenclature is not validated".into());
    }
    if let Some(name) = radical_naming::name_radical(&mol)? {
        return Ok(name);
    }
    validate_parent(&mol)?;
    let components = graph_io::get_connected_components(&mol);
    if components.len() != 1 {
        return Err("openclatura: expected one connected component".into());
    }
    let oxide = amine_oxides::name(&mol, &components[0])?;
    let phane = phane::name(&mol);
    let component =
        match component_namer::name_component_with_citation(&mol, &components[0], None, None) {
            Ok(component) => component,
            Err(error) => {
                if let Some(name) = oxide {
                    return Ok(name);
                }
                return phane.map_or(Err(error), |(text, _)| {
                    Ok(Name {
                        text,
                        status: nazuchem_ir::NameStatus::Unverified(nazuchem_ir::NameReason::Phane),
                        alternatives: vec![],
                        omitted_stereo: vec![],
                        omitted_units: 0,
                        verification: None,
                    })
                });
            }
        };
    let unmarked = mol.name_status()
        == nazuchem_ir::NameStatus::Unverified(nazuchem_ir::NameReason::UnauditedPath);
    assess_derivation(&mol, &component.derivation);
    // Preserve already certified fused oxide templates, including isoquinoline 2-oxide.
    if mol.name_status() != nazuchem_ir::NameStatus::Preferred
        && let Some(name) = oxide
    {
        return Ok(name);
    }
    let multiplicative =
        nazuchem_ir::NameStatus::NotPreferred(nazuchem_ir::NameReason::Multiplicative);
    let status = mol.name_status();
    // P-51.3.1, P-52.2.5.3: the multiplicative or linear phane name is
    // preferred; the substitutive name stays as an alternative.
    let preferred = if !unmarked {
        None
    } else if status == multiplicative {
        (!requires_skeletal_replacement(&mol, &component.derivation))
            .then(|| multiplicative::name(&mol, &components[0]))
            .flatten()
    } else if status == nazuchem_ir::NameStatus::NotPreferred(nazuchem_ir::NameReason::LinearPhane)
    {
        phane.clone().filter(|p| p.1).map(|p| p.0)
    } else {
        None
    };
    if let Some(name) = preferred {
        return Ok(Name {
            text: nest_enclosing_marks(&name),
            status: nazuchem_ir::NameStatus::Preferred,
            alternatives: vec![nazuchem_ir::AlternativeName {
                name: nest_enclosing_marks(&component.text),
                status,
            }],
            omitted_stereo: vec![],
            omitted_units: 0,
            verification: None,
        });
    }
    let text = component.text;
    if phane.as_ref().is_some_and(|p| !p.1) {
        mol.mark_name(nazuchem_ir::NameStatus::NotPreferred(
            nazuchem_ir::NameReason::Phane,
        ));
    }
    let status = mol.name_status();
    let mut alternatives = vec![];
    if mol.has_higher_order_fusion_alternative() || mol.has_benzo_fusion_alternative() {
        // FR-3.5 / P-25.3.5: retain the pre-regrouping citation for general use.
        // Cloning discards fusion plans cached under the preferred component policy.
        let higher_order = mol.clone();
        higher_order.set_benzoheterocycle_names(false);
        let alternative = component_namer::name_component(&higher_order, &components[0]);
        mol.set_benzoheterocycle_names(true);
        if let Ok(alternative) = alternative
            && alternative != text
        {
            alternatives.push(nazuchem_ir::AlternativeName {
                name: nest_enclosing_marks(&alternative),
                status: nazuchem_ir::NameStatus::NotPreferred(
                    if mol.has_benzo_fusion_alternative() {
                        nazuchem_ir::NameReason::BenzoFusion
                    } else {
                        nazuchem_ir::NameReason::HigherOrderFusion
                    },
                ),
            });
        }
    }
    if mol.has_relative_ring_alternative() {
        mol.set_relative_ring_names(true);
        if let Ok(relative) = component_namer::name_component(&mol, &components[0])
            && relative != text
        {
            alternatives.push(nazuchem_ir::AlternativeName {
                name: nest_enclosing_marks(&relative),
                status: nazuchem_ir::NameStatus::NotPreferred(
                    nazuchem_ir::NameReason::RelativeRingStereo,
                ),
            });
        }
        mol.set_relative_ring_names(false);
    }
    if mol.has_endo_exo_alternative() {
        mol.set_endo_exo_names(true);
        if let Ok(relative) = component_namer::name_component(&mol, &components[0])
            && relative != text
        {
            alternatives.push(nazuchem_ir::AlternativeName {
                name: nest_enclosing_marks(&relative),
                status: nazuchem_ir::NameStatus::NotPreferred(
                    nazuchem_ir::NameReason::EndoExoStereo,
                ),
            });
        }
    }
    alternatives.extend(phane.map(|(name, _)| nazuchem_ir::AlternativeName {
        name,
        status: nazuchem_ir::NameStatus::Unverified(nazuchem_ir::NameReason::Phane),
    }));
    Ok(Name {
        text: nest_enclosing_marks(&text),
        status,
        alternatives,
        omitted_stereo: vec![],
        omitted_units: 0,
        verification: None,
    })
}

/// Keep the isotope parent on the ordinary engine's existing validation path.
// Provenance: fn validate_parent | original
// Source: ReShiki | src/naming/openclatura/mod.rs::validate_parent
pub(crate) fn validate_parent(mol: &molecule::Molecule) -> Result<(), String> {
    ionic_naming::validate_heteroatom_oxo_charge_states(mol)?;
    // Upstream divergence: nitrogen stereocentres and ring stereocentres with a
    // lone pair (for example three single bonds on P, As, S+) read back through OPSIN as
    // another stereoisomer, or lose their descriptor.
    let cyclic = chains::get_cyclic_atoms(mol, &molecule::AtomSet::new());
    if mol.atoms.values().any(|atom| {
        atom.has_stereo()
            && (atom.symbol == "N"
                || !atom.is_carbon()
                    && cyclic.contains(&atom.idx)
                    && mol.degree(atom.idx) == 3
                    && has_nonbonding_pair(atom, 3)
                    && mol
                        .get_neighbors(atom.idx)
                        .iter()
                        .all(|&n| mol.bond_order(atom.idx, n) == 1))
    }) {
        return Err("openclatura: N and ring lone-pair stereocentres are not validated".into());
    }
    Ok(())
}

// Provenance: fn assess_derivation | original
// Source: ReShiki | src/naming/openclatura/mod.rs::assess_derivation
fn assess_derivation(mol: &molecule::Molecule, derivation: &component_namer::Derivation) {
    species::assess_anion_parent(mol, derivation);
    species::assess_cation_parent(mol, derivation);
    // P-51.3.2 / P-52.2.5.3: mandatory phane nomenclature takes precedence.
    if phane::has_linear_skeleton(mol, derivation.component_splits()) {
        mol.mark_name(nazuchem_ir::NameStatus::NotPreferred(
            nazuchem_ir::NameReason::LinearPhane,
        ));
    } else if has_ring_assembly(mol) || has_anionic_ring_assembly(mol) {
        mol.mark_name(nazuchem_ir::NameStatus::NotPreferred(
            nazuchem_ir::NameReason::RingAssembly,
        ));
    } else if has_multiplicative_pin(mol, derivation) {
        mol.mark_name(nazuchem_ir::NameStatus::NotPreferred(
            nazuchem_ir::NameReason::Multiplicative,
        ));
    }
    if requires_skeletal_replacement(mol, derivation) {
        mol.mark_name(nazuchem_ir::NameStatus::Unverified(
            nazuchem_ir::NameReason::UnauditedParent,
        ));
    }
    if matches!(
        derivation,
        component_namer::Derivation::Plan { .. }
            | component_namer::Derivation::Shortcut {
                certified: true,
                ..
            }
    ) {
        mol.certify();
    }
}

/// Rewrite the namer's enclosing parentheses in the IUPAC nesting order
/// `{[({[( )]})]}` (P-16.5.4), choosing each mark from what it encloses.
///
/// Square brackets (fusion, von Baeyer, spiro descriptors, locanted ring
/// assemblies), the braces of superscripts such as `1^{3,7}`, the parentheses
/// of a ring assembly's parent hydride and of indicated or added hydrogen are
/// part of the parent name and do not count (P-16.5.4.1.1,
/// P-16.5.4.1.2). Any other parenthesis, including stereodescriptors and
/// compound locants, counts as a level (P-16.5.4.1.3), and a mark that would
/// directly follow or precede one of the same kind moves up a level
/// (P-16.5.4.1.5).
///
/// Public so downstream regression harnesses, including ReShiki's, can compare
/// parentheses-only reference names with the engine's generated spelling.
// Provenance: fn nest_enclosing_marks | original
// Source: ReShiki | src/naming/openclatura/mod.rs::nest_enclosing_marks
pub fn nest_enclosing_marks(name: &str) -> String {
    // Per open group: its start, the deepest counted level inside it, and the
    // levels of counted groups that start right after it opens and end right
    // before it closes.
    struct Group {
        open: usize,
        inner: usize,
        first: Option<usize>,
        last: Option<(usize, usize)>,
    }
    let mut out = name.to_owned();
    let mut stack: Vec<Group> = Vec::new();
    for (i, byte) in name.bytes().enumerate() {
        match byte {
            b'(' | b'[' => stack.push(Group {
                open: i,
                inner: 0,
                first: None,
                last: None,
            }),
            close @ (b')' | b']') => {
                let Some(group) = stack.pop() else { continue };
                // Indicated or added hydrogen such as `(4aH,9H)` and the
                // multiplied parent of a ring assembly do not count
                // (P-16.5.4.1.1, P-16.5.4.1.2).
                let counted = close == b')'
                    && !["-bi", "-ter", "-quater", "-quinque", "-sexi"]
                        .iter()
                        .any(|multiplier| name[..group.open].ends_with(multiplier))
                    && !name[group.open + 1..i].split(',').all(|locant| {
                        locant.strip_suffix('H').is_some_and(|locant| {
                            let digits = locant
                                .trim_end_matches(|c: char| c == '\'' || c.is_ascii_lowercase());
                            !digits.is_empty() && digits.bytes().all(|b| b.is_ascii_digit())
                        })
                    });
                let mut level = group.inner + usize::from(counted);
                if counted {
                    let last = group.last.filter(|&(end, _)| end + 1 == i).map(|(_, l)| l);
                    while [group.first, last]
                        .into_iter()
                        .flatten()
                        .any(|l| l % 3 == level % 3)
                    {
                        level += 1;
                    }
                    let mark = (level - 1) % 3;
                    out.replace_range(group.open..=group.open, ["(", "[", "{"][mark]);
                    out.replace_range(i..=i, [")", "]", "}"][mark]);
                }
                if let Some(parent) = stack.last_mut() {
                    parent.inner = parent.inner.max(level);
                    if counted {
                        if group.open == parent.open + 1 {
                            parent.first = Some(level);
                        }
                        parent.last = Some((i, level));
                    }
                }
            }
            _ => {}
        }
    }
    out
}

/// P-51.4.1 / P-15.4.3.2: four heterounits on an acyclic path require
/// skeletal replacement, including in substituents (with imine N or
/// oxo-substituted P in the chain, P-51.4.1.1). Each skeletal heteroatom
/// counts as a unit, which can only overcount heterounits such as -SiH2-O-SiH2-;
/// terminal B, Al, Ga, In, Tl, Si, Ge, Sn, Pb, P, As, Sb and Bi count too (P-51.4.1.4). Only the cited
/// principal group's atoms are excluded.
// Provenance: fn requires_skeletal_replacement | original
// Source: ReShiki | src/naming/openclatura/mod.rs::requires_skeletal_replacement
fn requires_skeletal_replacement(
    mol: &molecule::Molecule,
    derivation: &component_namer::Derivation,
) -> bool {
    let principal_atoms = match derivation {
        component_namer::Derivation::Plan {
            principal_group, ..
        } => principal_group.as_ref().map(|(_, atoms)| atoms),
        component_namer::Derivation::Shortcut { .. } => None,
    };
    let cyclic = chains::get_cyclic_atoms(mol, &molecule::AtomSet::new());
    let skeletal: molecule::AtomSet = mol
        .atoms
        .values()
        .filter(|atom| {
            !atom.is_carbon()
                && atom.symbol != "H"
                && atom.charge == 0
                && !cyclic.contains(&atom.idx)
                && principal_atoms.is_none_or(|atoms| !atoms.contains(&atom.idx))
                && (mol
                    .get_neighbors(atom.idx)
                    .iter()
                    .filter(|n| mol.atoms[n].symbol != "H")
                    .count()
                    >= 2
                    || matches!(
                        atom.symbol.as_str(),
                        "B" | "Al"
                            | "Ga"
                            | "In"
                            | "Tl"
                            | "Si"
                            | "Ge"
                            | "Sn"
                            | "Pb"
                            | "P"
                            | "As"
                            | "Sb"
                            | "Bi"
                    ))
        })
        .map(|atom| atom.idx)
        .collect();
    skeletal.iter().any(|&start| {
        let mut pending = vec![(start, start, 1)];
        while let Some((previous, current, units)) = pending.pop() {
            if units >= 4 {
                return true;
            }
            for &next in mol.get_neighbors(current) {
                if next != previous && !cyclic.contains(&next) && mol.atoms[&next].symbol != "H" {
                    pending.push((current, next, units + usize::from(skeletal.contains(&next))));
                }
            }
        }
        false
    })
}

/// P-51.1.5 / P-51.3.1: multiply identical senior parents, including all
/// substituents, suffixes and their attachment positions. P-15.3 permits an
/// asymmetric central group, but its concatenated arms must be identical.
// Provenance: fn has_multiplicative_pin | adapted
// Source: openclatura | src/openclatura/component_namer.py::name_component
// Source: openclatura | src/openclatura/component_namer.py::select_component_parent
fn has_multiplicative_pin(
    mol: &molecule::Molecule,
    derivation: &component_namer::Derivation,
) -> bool {
    if matches!(derivation, component_namer::Derivation::Shortcut { .. }) {
        return false;
    }
    use molecule::AtomSet;
    let Some((parents, (paths, rings))) = multiplicative::parents(mol) else {
        return false;
    };
    let parents: Vec<_> = parents
        .into_iter()
        .map(|(p, complete)| (p.primary_path()[0], p.atom_set(), complete))
        .collect();
    if parents.len() < 2 {
        return false;
    }
    // Simple central groups use a chain, cyclic parent (including an assembly),
    // or heteroatom. Side substituents belong to the central group and need
    // not be symmetric; the branches from it towards the parents must be.
    let mut cores: Vec<AtomSet> = paths.iter().map(|p| p.iter().copied().collect()).collect();
    cores.extend(rings.iter().map(|r| r.atoms.clone()));
    cores.extend(
        ring_assembly::discover_ring_assemblies(mol, &rings)
            .into_iter()
            .map(|a| a.locant_maps[0].keys().copied().collect()),
    );
    cores.extend(
        mol.atoms
            .values()
            .filter(|a| !a.is_carbon() && a.symbol != "H")
            .map(|a| AtomSet::from([a.idx])),
    );
    for (i, (root_a, a, required_a)) in parents.iter().enumerate() {
        for (root_b, b, required_b) in &parents[i + 1..] {
            if !required_a.is_disjoint(required_b) {
                continue;
            }
            for core in &cores {
                if !core.is_disjoint(required_a) || !core.is_disjoint(required_b) {
                    continue;
                }
                let left = subgraph_tools::subgraph_component(mol, *root_a, core);
                if !required_a.is_subset(&left) || !left.is_disjoint(required_b) {
                    continue;
                }
                let right = subgraph_tools::subgraph_component(mol, *root_b, core);
                if !required_b.is_subset(&right) || !left.is_disjoint(&right) {
                    continue;
                }
                let attachments = |branch: &AtomSet| {
                    branch
                        .iter()
                        .flat_map(|&u| {
                            mol.get_neighbors(u)
                                .iter()
                                .filter(|v| core.contains(v))
                                .map(move |&v| (u, v))
                        })
                        .collect::<Vec<_>>()
                };
                let (left_edges, right_edges) = (attachments(&left), attachments(&right));
                let ([(u, x)], [(v, y)]) = (&left_edges[..], &right_edges[..]) else {
                    continue;
                };
                if multiplicative_bond(mol, *u, *x) == multiplicative_bond(mol, *v, *y)
                    && identical_multiplicative_branches(mol, &left, &right, a, b, *u, *v)
                {
                    return true;
                }
            }
        }
    }
    false
}

// Provenance: fn multiplicative_bond | original
// Source: ReShiki | src/naming/openclatura/mod.rs::multiplicative_bond
fn multiplicative_bond(
    mol: &molecule::Molecule,
    u: usize,
    v: usize,
) -> Option<(usize, Option<&str>)> {
    mol.get_bond(u, v)
        .map(|b| (if b.is_aromatic { 4 } else { b.order }, b.stereo.as_deref()))
}

/// A rooted graph proof preserves complete parent units and configuration
/// (P-51.3.1, P-45.6.2), rather than comparing unrooted hydrides or global ranks.
// Provenance: fn identical_multiplicative_branches | adapted
// Source: openclatura | src/openclatura/locant_elision.py::_isomorphic
// Source: openclatura | src/openclatura/locant_elision.py::_has_parent_set_automorphism_mapping
fn identical_multiplicative_branches(
    mol: &molecule::Molecule,
    left: &molecule::AtomSet,
    right: &molecule::AtomSet,
    parent_left: &molecule::AtomSet,
    parent_right: &molecule::AtomSet,
    root_left: usize,
    root_right: usize,
) -> bool {
    use std::collections::BTreeMap;
    if left.len() != right.len() {
        return false;
    }
    let ranks_left = chains::canonical_ranks(&mol.subgraph(left), None);
    let ranks_right = chains::canonical_ranks(&mol.subgraph(right), None);
    let atom_key = |u, atoms: &molecule::AtomSet, parent: &molecule::AtomSet, root| {
        let a = &mol.atoms[&u];
        (
            a.symbol.as_str(),
            a.charge,
            a.isotope,
            a.is_aromatic,
            a.total_h_count,
            a.radical_electrons,
            a.stereo.as_deref(),
            mol.get_neighbors(u)
                .iter()
                .filter(|n| atoms.contains(n))
                .count(),
            parent.contains(&u),
            u == root,
        )
    };
    let candidates: BTreeMap<_, Vec<_>> = left
        .iter()
        .map(|&u| {
            (
                u,
                right
                    .iter()
                    .copied()
                    .filter(|&v| {
                        ranks_left[&u] == ranks_right[&v]
                            && atom_key(u, left, parent_left, root_left)
                                == atom_key(v, right, parent_right, root_right)
                    })
                    .collect(),
            )
        })
        .collect();
    // Provenance: fn identical_multiplicative_branches::matches | adapted
    // Source: openclatura | src/openclatura/locant_elision.py::_isomorphic
    fn matches(
        mol: &molecule::Molecule,
        candidates: &BTreeMap<usize, Vec<usize>>,
        mapping: &mut BTreeMap<usize, usize>,
    ) -> bool {
        let Some((&u, options)) = candidates
            .iter()
            .filter(|(u, _)| !mapping.contains_key(u))
            .min_by_key(|(_, options)| {
                options
                    .iter()
                    .filter(|v| !mapping.values().any(|x| x == *v))
                    .count()
            })
        else {
            return true;
        };
        for &v in options {
            if mapping.values().any(|&x| x == v)
                || mapping.iter().any(|(&a, &b)| {
                    multiplicative_bond(mol, u, a) != multiplicative_bond(mol, v, b)
                })
            {
                continue;
            }
            mapping.insert(u, v);
            if matches(mol, candidates, mapping) {
                return true;
            }
            mapping.remove(&u);
        }
        false
    }
    matches(mol, &candidates, &mut BTreeMap::new())
}

/// Whether two identical ring systems are joined directly by a single or
/// double bond that the rendered name does not cite as a ring assembly
/// (P-28, P-52.2.7): a double-bond junction or units the assembly parent
/// cannot number. Systems compare as parent hydrides, so
/// substituents count as hydrogen; different systems keep the substitutive
/// PIN, cyclohexylbenzene (P-61.2.2).
// Provenance: fn has_ring_assembly | original
// Source: ReShiki | src/naming/openclatura/mod.rs::has_ring_assembly
fn has_ring_assembly(mol: &molecule::Molecule) -> bool {
    ring_assembly::assembly_junctions(mol)
        .iter()
        .any(|&junction| !mol.ring_assembly_is_named(junction))
}

/// Removing or adding a hydron does not change the parent hydride, so a
/// pyrrolide joined to a pyrrole still has a ring assembly parent (P-28,
/// P-72.2.2.1, P-73.1).
// Provenance: fn has_anionic_ring_assembly | adapted
// Source: openclatura | src/openclatura/fusion/charges.py::fusion_component_charge_parent
fn has_anionic_ring_assembly(mol: &molecule::Molecule) -> bool {
    if mol.atoms.values().all(|atom| atom.charge != -1) {
        return false;
    }
    let mut neutral = mol.clone();
    for atom in neutral.atoms.values_mut() {
        if atom.charge == -1 {
            atom.charge = 0;
            atom.total_h_count += 1;
            atom.explicit_h_count += 1;
        } else if atom.charge == 1 && atom.total_h_count > 0 {
            atom.charge = 0;
            atom.total_h_count -= 1;
            atom.explicit_h_count = atom.explicit_h_count.saturating_sub(1);
        }
    }
    ring_assembly::assembly_junctions(&neutral)
        .iter()
        .any(|&junction| !mol.ring_assembly_is_named(junction))
}

/// Count nonbonding electrons on the main-group elements that can carry
/// tetrahedral tags. Attached hydrogens also consume bonding electrons:
/// a three-neighbour Si–H centre has four bonds and no lone pair.
// Provenance: fn has_nonbonding_pair | original
// Source: ReShiki | src/naming/openclatura/mod.rs::has_nonbonding_pair
fn has_nonbonding_pair(atom: &molecule::Atom, bond_order_sum: usize) -> bool {
    let outer_electrons = match atom.symbol.as_str() {
        "C" | "Si" | "Ge" | "Sn" | "Pb" => 4,
        "N" | "P" | "As" | "Sb" | "Bi" => 5,
        "O" | "S" | "Se" | "Te" | "Po" => 6,
        "F" | "Cl" | "Br" | "I" | "At" => 7,
        _ => return false,
    };
    outer_electrons
        - atom.charge
        - (bond_order_sum + atom.total_h_count + atom.radical_electrons) as i32
        >= 2
}

#[cfg(test)]
mod test_support;

#[cfg(test)]
// Provenance: mod tests | original
// Source: ReShiki | src/naming/openclatura/mod.rs::tests
mod tests {
    use super::{name_smiles, nest_enclosing_marks};
    use std::io::Write;

    #[test]
    // Provenance: fn tests::naming_routes_use_the_supplied_reservation_policy | adapted
    // Source: ReShiki | src/naming/openclatura/mod.rs::tests::naming_routes_use_the_supplied_reservation_policy
    // Source: openclatura | tests/unit/test_analysis.py::test_public_api_golden_names
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn naming_routes_use_the_supplied_reservation_policy() {
        let budget =
            nazuchem_kernel::budget::Budget::new(usize::MAX, usize::MAX, None, usize::MAX, |_| {
                None
            });
        for smiles in ["CCO", "[13CH3]CO", "[Na+].[Cl-]", "C[N+](C)C"] {
            let error = name_smiles(smiles, &budget).unwrap_err();
            assert!(error.contains("resources are busy"), "{smiles}: {error}");
        }
    }

    /// Issue #10: the name follows the rule, not the input atom order. Each
    /// input is also rewritten from several root atoms.
    #[test]
    fn names_do_not_depend_on_atom_order() {
        use nazuchem_kernel::smiles;

        for (expected, inputs) in [
            // P-31.1.4.2 (1)-(3): fewest compound locants, then lowest locants. The
            // third input carries atom maps.
            (
                "21,22,23,24-tetraazapentacyclo[16.2.1.1^{3,6}.1^{8,11}.1^{13,16}]tetracosa-1(21),2,4,6(24),7,9,11,13,15,17,19-undecaene",
                &[
                    "C1=Cc2cc3ccc(cc4ccc(cc5nc(cc1n2)C=C5)[nH]4)[nH]3",
                    "C1=Cc2cc3nc(cc4ccc(cc5ccc(cc1n2)[nH]5)[nH]4)C=C3",
                    "[CH:24]1=[CH:23][c:22]2[cH:21][c:20]3[cH:19][cH:18][c:17]([cH:16][c:15]4[cH:14][cH:13][c:12]([cH:11][c:10]5[n:9][c:8]([cH:7][c:6]1[n:5]2)[CH:4]=[CH:3]5)[nH:2]4)[nH:1]3",
                ][..],
            ),
            ("bicyclo[4.2.0]oct-6-ene", &["C1CCC2=CCC2C1"][..]),
            (
                "bicyclo[6.5.1]tetradec-8-ene",
                &["C12CCCCCCC(=CCCCC1)C2"][..],
            ),
            (
                "bicyclo[12.2.2]octadeca-1(16),14,17-triene",
                &["C1=CC2=CC=C1CCCCCCCCCCCC2"][..],
            ),
            ("bicyclo[4.1.0]hepta-1,3,5-triene", &["C1C2=CC=CC=C12"][..]),
            (
                "tricyclo[9.3.1.1^{4,8}]hexadeca-1(15),11,13-triene",
                &["C1CC2CCCC(C2)CCc2cccc1c2"][..],
            ),
            (
                "tricyclo[9.3.1.1^{4,8}]hexadeca-1(15),4(16),5,7,11,13-hexaene",
                &["C1Cc2cccc(c2)CCc2cccc1c2"][..],
            ),
            (
                "tetracyclo[7.7.1.1^{3,7}.1^{11,15}]nonadeca-3,11(18)-diene",
                &["C12CC3=CCCC(C3)CC(C2)CC(=C4)CCCC4C1"][..],
            ),
            ("bicyclo[2.2.2]oct-1-ene", &["C1CC2=CCC1CC2"][..]),
            // P-44.1.1, P-52.2.8: one suffix class on a ring and a chain; the
            // ring is the parent.
            (
                "N-[4-(3-methylbutanamido)phenyl]-2-(phenylsulfanyl)pyridine-3-carboxamide",
                &[
                    "CC(C)CC(=O)Nc1ccc(NC(=O)c2cccnc2Sc2ccccc2)cc1",
                    "c1c(NC(=O)c2cccnc2Sc2ccccc2)ccc(NC(CC(C)C)=O)c1",
                ][..],
            ),
            (
                "4-chloro-N-[(2S)-1-{[2-(cyclohex-1-en-1-yl)ethyl]amino}-3-methyl-1-oxobutan-2-yl]benzamide",
                &[
                    "CC(C)[C@H](NC(=O)c1ccc(Cl)cc1)C(=O)NCCC1=CCCCC1",
                    "CC([C@@H](C(NCCC1=CCCCC1)=O)NC(c1ccc(Cl)cc1)=O)C",
                ][..],
            ),
            (
                "4-(2-carboxyethyl)benzoic acid",
                &["OC(=O)CCc1ccc(C(=O)O)cc1"][..],
            ),
            ("4-(3-oxopropyl)benzaldehyde", &["O=CCCc1ccc(C=O)cc1"][..]),
            ("N'-acetylbenzohydrazide", &["CC(=O)NNC(=O)c1ccccc1"][..]),
            // P-66.1.6.1.1.5: a carboxylic amide is senior to urea.
            (
                "2-(cyclopropyl{[4-(difluoromethoxy)-3-ethoxyphenyl]methyl}amino)-N-[(2-methylpropyl)carbamoyl]acetamide",
                &[
                    "CCOc1cc(CN(CC(=O)NC(=O)NCC(C)C)C2CC2)ccc1OC(F)F",
                    "N(C(NC(CN(Cc1cc(OCC)c(OC(F)F)cc1)C1CC1)=O)=O)CC(C)C",
                ][..],
            ),
            ("N-carbamoylbenzamide", &["NC(=O)NC(=O)c1ccccc1"][..]),
            (
                "N-carbamoyl-2-phenylacetamide",
                &["NC(=O)NC(=O)Cc1ccccc1"][..],
            ),
            (
                "N-[2-(carbamoylamino)ethyl]acetamide",
                &["CC(=O)NCCNC(N)=O"][..],
            ),
            (
                "N-[3-(carbamoylamino)propyl]formamide",
                &["O=CNCCCNC(N)=O"][..],
            ),
            ("N-(methylcarbamothioyl)acetamide", &["CC(=O)NC(=S)NC"][..]),
            // P-66.1.2: (R-CO)2NH is an N-acyl derivative of the senior amide.
            ("N-acetylbenzamide", &["CC(=O)NC(=O)c1ccccc1"][..]),
            (
                "N-acetyl-N-(3-chloropropanoyl)benzamide",
                &["CC(=O)N(C(=O)CCCl)C(=O)c1ccccc1"][..],
            ),
            ("N-acetylpropanamide", &["CC(=O)NC(=O)CC"][..]),
            (
                "N-acetyl-N-cyclopentylacetamide",
                &["CC(=O)N(C1CCCC1)C(C)=O"][..],
            ),
            // P-14.4 (f): lowest locants to the prefixes.
            (
                "1-(3-ethyl-4-methylhexyl)octahydro-2,5-methanopentalene",
                &[
                    "CCC(C)C(CC)CCC1C2CC3CC(C2)C1C3",
                    "C1C2CC3C(CCC(C(CC)C)CC)C1C(C2)C3",
                ][..],
            ),
            (
                "N^3-(6-chloro-2-methoxyacridin-9-yl)-N^1-{3-[(6-chloro-2-methoxyacridin-9-yl)amino]propyl}-N^1-(cyclopropylmethyl)propane-1,3-diamine",
                &[
                    "COc1ccc2nc3cc(Cl)ccc3c(NCCCN(CCCNc3c4ccc(Cl)cc4nc4ccc(OC)cc34)CC3CC3)c2c1",
                    "C(CCNc1c2c(nc3ccc(OC)cc31)cc(Cl)cc2)N(CCCNc1c2c(nc3ccc(OC)cc13)cc(Cl)cc2)CC1CC1",
                ][..],
            ),
            // P-14.4 (g): the prefix cited first takes the lower locant.
            (
                "N-(1'-methyl[1,4'-bipiperidin]-4-yl)-N'-[3-(1H-pyrazol-1-yl)phenyl]ethanediamide",
                &[
                    "CN1CCC(N2CCC(NC(=O)C(=O)Nc3cccc(-n4cccn4)c3)CC2)CC1",
                    "c1(NC(C(NC2CCN(C3CCN(C)CC3)CC2)=O)=O)cccc(-n2cccn2)c1",
                ][..],
            ),
            // Issue #107. P-14.4 (j): geminal amine nitrogens whose substituents
            // differ only by stereodescriptor; R, like Z, takes the lower,
            // unprimed locant (P-14.3.5).
            (
                "N^3-[(2R)-butan-2-yl]-N'^3-[(2S)-butan-2-yl]hexane-3,3-diamine",
                &[
                    "CCC(N[C@H](C)CC)(N[C@@H](C)CC)CCC",
                    "CCC(N[C@@H](C)CC)(N[C@H](C)CC)CCC",
                    "CC[C@H](C)NC(CC)(CCC)N[C@H](C)CC",
                ][..],
            ),
            (
                "N^3-[(2Z)-but-2-en-1-yl]-N'^3-[(2E)-but-2-en-1-yl]hexane-3,3-diamine",
                &[
                    "CCC(NC/C=C\\C)(NC/C=C/C)CCC",
                    "CCC(NC/C=C/C)(NC/C=C\\C)CCC",
                    "C/C=C/CNC(CC)(CCC)NC/C=C\\C",
                ][..],
            ),
            (
                "N^1-[(2R)-butan-2-yl]-N'^1-[(2S)-butan-2-yl]cyclohexane-1,1-diamine",
                &[
                    "C[C@H](CC)NC1(N[C@@H](C)CC)CCCCC1",
                    "C[C@@H](CC)NC1(N[C@H](C)CC)CCCCC1",
                ][..],
            ),
            (
                "N^1-[(2R)-butan-2-yl]-N'^1,N''^1-bis[(2S)-butan-2-yl]ethane-1,1,1-triamine",
                &[
                    "CC(N[C@H](C)CC)(N[C@@H](C)CC)N[C@@H](C)CC",
                    "CC(N[C@@H](C)CC)(N[C@@H](C)CC)N[C@H](C)CC",
                ][..],
            ),
            // A descriptor nested in the substituent name counts the same way.
            (
                "N^3-{2-[(2R)-butan-2-yl]phenyl}-N'^3-{2-[(2S)-butan-2-yl]phenyl}hexane-3,3-diamine",
                &[
                    "CCC(Nc1ccccc1[C@H](C)CC)(Nc1ccccc1[C@@H](C)CC)CCC",
                    "CCC(Nc1ccccc1[C@@H](C)CC)(Nc1ccccc1[C@H](C)CC)CCC",
                ][..],
            ),
            (
                "N^3-{[(2R)-oxolan-2-yl]methyl}-N'^3-{[(2S)-oxolan-2-yl]methyl}hexane-3,3-diamine",
                &[
                    "CCC(NC[C@H]1CCCO1)(NC[C@@H]1CCCO1)CCC",
                    "CCC(NC[C@@H]1CCCO1)(NC[C@H]1CCCO1)CCC",
                ][..],
            ),
            // ethyl, cited before methyl, stays unprimed although it is S; the R
            // prefix is still cited first (P-45.6.3).
            (
                "N'^3-[(2R)-butan-2-yl]-N^3-[(2S)-butan-2-yl]-N^3-ethyl-N'^3-methylhexane-3,3-diamine",
                &[
                    "CCC(N([C@H](C)CC)C)(N([C@@H](C)CC)CC)CCC",
                    "CCC(N([C@@H](C)CC)CC)(N([C@H](C)CC)C)CCC",
                ][..],
            ),
            // P-65.6.3.2.1: ester groups in alphanumerical order, as in ethyl
            // methyl butanedioate.
            (
                "4-(tert-butoxy)-2,2-dimethyl-4-oxobutyl methyl propanedioate",
                &[
                    "COC(=O)CC(=O)OCC(C)(C)CC(=O)OC(C)(C)C",
                    "C(OC(CC(OC)=O)=O)C(C)(C)CC(=O)OC(C)(C)C",
                ][..],
            ),
            // Issue #97. P-14.4 (c): the -ium suffix takes the
            // lower of the two equivalent locants. Issue #96, P-73.1.1.2: the
            // cation keeps the retained name of its neutral parent hydride.
            (
                "decahydro-2,7-naphthyridin-2-ium",
                &["N1CCC2C(C1)C[NH2+]CC2", "C1C[NH2+]CC2CNCCC12"][..],
            ),
            // P-14.3.5: on a symmetric assembly the chalcogenide nitrogen is
            // 1, not 1′.
            (
                "4,4'-bipyridine 1-sulfide",
                &["[S-][n+]1ccc(-c2ccncc2)cc1", "c1cc(-c2cc[n+]([S-])cc2)ccn1"][..],
            ),
            // Controls: a substituent fixes the numbering, so the oxide keeps a
            // primed locant; on the meso piperazine only a reflection relates
            // the nitrogens, so each oxide keeps its own locant.
            (
                "2-methyl-4,4'-bipyridine 1'-oxide",
                &["Cc1cc(-c2cc[n+]([O-])cc2)ccn1"][..],
            ),
            // Every nitrogen shares one Morgan rank, but N1 and N4 are not
            // equivalent.
            (
                "4-[2-(aziridin-1-yl)ethyl]-1,4,7-triazabicyclo[5.2.2]undecane 4-oxide",
                &["[N+]1([O-])(CCN2CC2)CCN3CCN(CC1)CC3"][..],
            ),
            (
                "(2R,5S)-1,2,4,5-tetramethylpiperazine 1-oxide",
                &["C[N+]1([O-])[C@H](C)CN(C)[C@@H](C)C1"][..],
            ),
            (
                "(2R,5S)-1,2,4,5-tetramethylpiperazine 4-oxide",
                &["CN1[C@H](C)C[N+](C)([O-])[C@@H](C)C1"][..],
            ),
            // Issue #128. P-14.3.5: a rotation that keeps every R/S descriptor
            // makes the nitrogens equivalent, so the oxide takes the lower locant.
            (
                "(2R,5R)-1,2,4,5-tetramethylpiperazine 1-oxide",
                &[
                    "C[C@@H]1C[N+](C)([O-])[C@H](C)CN1C",
                    "CN1C[C@@H](C)[N+](C[C@H]1C)(C)[O-]",
                ][..],
            ),
            // The three nitrogens are related only by a three-fold rotation.
            (
                "3,8,13-trimethyl-1,6,11-trioxa-3,8,13-triazacyclopentadecane 3-oxide",
                &[
                    "C[N+]1([O-])COCCN(C)COCCN(C)COCC1",
                    "C1OCN(CCOC[N+](CCOCN(C)C1)(C)[O-])C",
                ][..],
            ),
            (
                "1,4,7-trimethyl-1,4,7-triazonane 1-oxide",
                &["C[N+]1([O-])CCN(C)CCN(C)CC1"][..],
            ),
            // P-66.4.1.1, P-66.4.1.2.1.1, P-66.1.6.1.1.5: amidines rank as amides,
            // in the order of their acids, so a carboximidamide is senior to
            // guanidine (the amidine of carbonic acid).
            (
                "4-{2-[(diaminomethylidene)amino]ethyl}benzene-1-carboximidamide",
                &["NC(N)=NCCc1ccc(C(N)=N)cc1", "NC(=N)c1ccc(CCN=C(N)N)cc1"][..],
            ),
            // Issue #128. P-41 classes 11 and 12, P-43 Table 4.4 (suffixes 16 to
            // 19 and 31): amides in the order of their acids, carboxamides and
            // their chalcogen analogues, carboximidamides, then sulfonamides,
            // all before hydrazides.
            (
                "4-(hydrazinecarbonyl)benzene-1-carboximidamide",
                &["NC(=N)c1ccc(cc1)C(=O)NN", "NNC(=O)c1ccc(C(N)=N)cc1"][..],
            ),
            (
                "4-(2-acetylhydrazine-1-carbonyl)benzene-1-carboximidamide",
                &["NC(=N)c1ccc(cc1)C(=O)NNC(C)=O"][..],
            ),
            (
                "4-sulfamoylbenzene-1-carboximidamide",
                &["NC(=N)c1ccc(cc1)S(N)(=O)=O"][..],
            ),
            (
                "4-carbamimidoylbenzene-1-carbothioamide",
                &["NC(=N)c1ccc(cc1)C(N)=S"][..],
            ),
            // Issue #141. P-66.3.2.2: a chain-end hydrazide is hydrazinyl with
            // oxo, on the chain that keeps the carbonyl carbon.
            (
                "4-hydrazinyl-4-oxobutanimidamide",
                &["NC(=N)CCC(=O)NN", "NNC(=O)CCC(N)=N"][..],
            ),
            (
                "5-hydrazinyl-3-(2-hydrazinyl-2-oxoethyl)-5-oxopentanoic acid",
                &["OC(=O)CC(CC(=O)NN)CC(=O)NN", "NNC(=O)CC(CC(=O)NN)CC(O)=O"][..],
            ),
            (
                "4-hydrazinyl-3-methyl-4-oxobutanoic acid",
                &["OC(=O)CC(C)C(=O)NN", "CC(CC(O)=O)C(=O)NN"][..],
            ),
            // Issue #148. P-66.1.1.4.1.1 (1): a chain-end amide is amino with
            // oxo, on the chain that keeps the carbonyl carbon (P-45.2.1).
            (
                "4-amino-4-oxobutanoic acid",
                &["OC(=O)CCC(N)=O", "NC(=O)CCC(O)=O"][..],
            ),
            (
                "4-amino-3-methyl-4-oxobutanoic acid",
                &["OC(=O)CC(C)C(N)=O", "CC(CC(O)=O)C(N)=O"][..],
            ),
            (
                "5-amino-3-(2-amino-2-oxoethyl)-5-oxopentanoic acid",
                &["OC(=O)CC(CC(N)=O)CC(N)=O", "NC(=O)CC(CC(N)=O)CC(O)=O"][..],
            ),
            // P-66.3.2.1: the acyl prefix of hydrazinecarbothioic acid.
            (
                "4-(hydrazinecarbothioyl)benzene-1-carboximidamide",
                &["NC(=N)c1ccc(cc1)C(=S)NN"][..],
            ),
            // P-41 classes 7a and 7b, P-66.1.6.1.1.5 (N-carbamoylbenzenesulfonamide
            // in P-66.1.6.1.1.3): amides of suffix acids are senior to urea.
            (
                "4-(carbamoylamino)benzene-1-sulfonamide",
                &["NC(=O)Nc1ccc(cc1)S(N)(=O)=O"][..],
            ),
            (
                "4-(carbamoylamino)benzene-1-carbothioamide",
                &["NC(=O)Nc1ccc(cc1)C(N)=S"][..],
            ),
            (
                "4-(carbamoylamino)benzene-1-carboximidamide",
                &["NC(=O)Nc1ccc(cc1)C(N)=N"][..],
            ),
            (
                "2-(carbamothioylamino)ethane-1-sulfonamide",
                &["NC(=S)NCCS(N)(=O)=O"][..],
            ),
            // P-66.4.1.2.1.3 (N-carbamimidoylurea): urea is senior to guanidine.
            (
                "N-{2-[(diaminomethylidene)amino]ethyl}urea",
                &["NC(N)=NCCNC(N)=O"][..],
            ),
            // P-41 class 15, P-43 Table 4.4 (suffix 47): an aldehyde is senior
            // to its chalcogen analogues.
            ("4-sulfanylidenebutanal", &["O=CCCC=S", "S=CCCC=O"][..]),
            (
                "3-(methanethioyl)pentanedial",
                &["O=CCC(C=S)CC=O", "C(C=S)(CC=O)CC=O"][..],
            ),
            (
                "4-sulfanylidene-2-(2-sulfanylideneethyl)butanal",
                &["O=CC(CC=S)CC=S", "S=CCC(C=O)CC=S"][..],
            ),
            // P-66.5.4.1: nitrile oxides rank with zwitterions, ahead of nitriles.
            (
                "3-cyanopropanenitrile oxide",
                &["N#CCCC#[N+][O-]", "[O-][N+]#CCCC#N"][..],
            ),
            // P-68.3.1.2.2 ranks no hydrazone over another; the functional
            // class names take the aldehyde hydrazone whatever the atom order.
            (
                "4-(hydrazono)pentanal hydrazone",
                &["NN=CCCC(C)=NN", "C(CC(=NN)C)C=NN"][..],
            ),
            // P-44.4.1.12.2: of two butan-2-ol parents the like pair (2R,3R)
            // precedes the unlike pair (2R,3S).
            (
                "(2R,3R)-3-[(3-{[(2S,3R)-3-hydroxybutan-2-yl]amino}hexan-3-yl)amino]butan-2-ol",
                &[
                    "CCC(N[C@H](C)[C@@H](C)O)(N[C@@H](C)[C@@H](C)O)CCC",
                    "CCC(N[C@@H](C)[C@@H](C)O)(N[C@H](C)[C@@H](C)O)CCC",
                    "O[C@H](C)[C@H](C)NC(CC)(CCC)N[C@H](C)[C@@H](C)O",
                ][..],
            ),
        ] {
            let identity = crate::test_support::canonical_smiles(inputs[0]).unwrap();
            let parsed = crate::test_support::parse_strict(expected).unwrap();
            assert_eq!(
                crate::test_support::canonical_smiles(&parsed).unwrap(),
                identity,
                "{expected}"
            );
            let status = name_smiles(inputs[0], &Default::default()).unwrap().status;
            for &input in inputs {
                let prepared = smiles::prepare(input, &Default::default()).unwrap();
                let atoms = prepared.state.graph.atoms.len();
                let orders = (0..atoms).step_by(atoms.div_ceil(5)).map(|root| {
                    smiles::write::write(
                        &prepared.state,
                        smiles::write::Options {
                            canonical: false,
                            root: Some(root),
                            ..Default::default()
                        },
                    )
                    .unwrap()
                    .text
                });
                for order in std::iter::once(input.to_owned()).chain(orders) {
                    let name = name_smiles(&order, &Default::default()).unwrap();
                    assert_eq!(name.text, expected, "{input} written as {order}");
                    assert_eq!(name.status, status, "{input} written as {order}");
                }
            }
        }
    }

    #[test]
    fn symmetric_azine_ylidene_spelling_keeps_unaudited_status() {
        use nazuchem_ir::{NameReason, NameStatus};

        for (smiles, expected) in [
            ("C/C=N/N=C/C", "(1E,2E)-diethylidenehydrazine"),
            ("CC/C=N/N=C/CC", "(1E,2E)-dipropylidenehydrazine"),
            (
                "C/C(/CC)=N\\N=C(/C)\\CC",
                "(1E,2E)-di(butan-2-ylidene)hydrazine",
            ),
            ("CC(C)=NN=C(C)C", "di(propan-2-ylidene)hydrazine"),
            ("CCCCCCCCCC=NN=CCCCCCCCCC", "di(decylidene)hydrazine"),
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(
                name.status,
                NameStatus::Unverified(NameReason::UnauditedPath),
                "{smiles}"
            );
        }
        assert!(name_smiles("C/C=N/N=C\\C", &Default::default()).is_err());
    }

    #[test]
    // Provenance: fn tests::azine_spelling_preserves_alkene_control | adapted
    // Source: ReShiki | src/naming/openclatura/mod.rs::tests::azine_spelling_preserves_alkene_control
    // Source: openclatura | tests/roundtrip/test_stereo.py::test_roundtrip
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn azine_spelling_preserves_alkene_control() {
        let name = name_smiles("C/C=C/C", &Default::default()).unwrap();
        assert_eq!(name.text, "(2E)-but-2-ene");
        assert_eq!(name.status, nazuchem_ir::NameStatus::Preferred);
    }

    #[test]
    // Provenance: fn tests::alkoxy_contractions_preserve_parent_components | adapted
    // Source: ReShiki | src/naming/openclatura/mod.rs::tests::alkoxy_contractions_preserve_parent_components
    // Source: openclatura | tests/fixtures/diverse_corpus.csv::<fixture data>
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn alkoxy_contractions_preserve_parent_components() {
        for (smiles, expected) in [
            ("OC(=O)COC(C)C", "[(propan-2-yl)oxy]acetic acid"),
            ("OC(=O)COCC(C)C", "(2-methylpropoxy)acetic acid"),
            ("OC(=O)COCCCCCCl", "[(5-chloropentyl)oxy]acetic acid"),
            (
                "OC(=O)COC(C)(C)CCl",
                "[(1-chloro-2-methylpropan-2-yl)oxy]acetic acid",
            ),
            (
                "OC(=O)COc1ccc(-c2ccccc2)cc1",
                "[([1,1'-biphenyl]-4-yl)oxy]acetic acid",
            ),
            ("OC(=O)COCC1CCCCC1", "(cyclohexylmethoxy)acetic acid"),
            ("OC(=O)COCc1ccccc1", "(benzyloxy)acetic acid"),
            ("OC(=O)COC(C)CCl", "[(1-chloropropan-2-yl)oxy]acetic acid"),
            (
                "OC(=O)CSc1ccc(Cl)cc1",
                "[(4-chlorophenyl)sulfanyl]acetic acid",
            ),
            ("OC(=O)COCCO", "(2-hydroxyethoxy)acetic acid"),
            ("OC(=O)COC(C)(C)C", "tert-butoxyacetic acid"),
            ("ClCCOC(=O)c1ccccc1", "2-chloroethyl benzoate"),
            ("COc1ccccc1", "anisole"),
            ("OC(=O)COC(C)CC", "[(butan-2-yl)oxy]acetic acid"),
            ("OC(=O)C(OCCO)OCCO", "bis(2-hydroxyethoxy)acetic acid"),
            (
                "ClCOc1ccc(CCl)cc1",
                "1-(chloromethoxy)-4-(chloromethyl)benzene",
            ),
            ("ClCCOC(=O)CC(=O)O", "(2-chloroethoxycarbonyl)acetic acid"),
            ("CC(C)CO[SiH2]C", "methyl(2-methylpropoxy)silane"),
        ] {
            assert_eq!(
                name_smiles(smiles, &Default::default()).unwrap().text,
                expected,
                "{smiles}"
            );
        }
    }

    #[test]
    fn skeletal_replacement_counts_heterounits_and_the_cited_group() {
        for (smiles, required) in [
            ("COCOCOC", false),
            ("COCOCOCOC", true),
            ("CNNOOC", true),
            ("CO[SiH2]OCCOCC(=O)O", true),
            ("COCCOCCOC[SiH3]", true),
            ("COCCOCCOC[PH2]", true),
            ("COCCOCC=NCCOCC(=O)O", true),
            ("COCOCOCC(=O)NC", false),
            ("COCOCOCC(=O)NCC(=O)O", true),
            ("CC(=O)NCC(=O)NCC(=O)NCC(=O)NCC(=O)O", true),
        ] {
            let mol = super::graph_io::read_smiles(smiles, &Default::default()).unwrap();
            let atoms = mol.atoms.keys().copied().collect();
            let component =
                super::component_namer::name_component_with_citation(&mol, &atoms, None, None)
                    .unwrap();
            assert_eq!(
                super::requires_skeletal_replacement(&mol, &component.derivation),
                required,
                "{smiles}"
            );
        }
    }

    #[test]
    // Provenance: fn tests::allow_list_status_matrix | adapted
    // Source: ReShiki | src/naming/openclatura/mod.rs::tests::allow_list_status_matrix
    // Source: openclatura | tests/fixtures/diverse_corpus.golden.json::<fixture data>
    // Source: openclatura | tests/unit/test_chalcogen_functional_groups.py::test_anhydride_bridges_and_acyl_sites_are_independent
    // Source: openclatura | tests/unit/test_fused_bonds.py::test_fused_bonds_match_the_cycle_count
    // Source: openclatura | tests/roundtrip/test_replacement_parents.py::test_roundtrip
    // Source: openclatura | tests/data/parser_xml_resources/simpleGroups.json::<fixture data>
    // Source: openclatura | tests/unit/test_analysis.py::test_public_api_golden_names
    // Source: openclatura | tests/unit/test_analysis.py::test_group_13_14_central_atoms_have_a_parent_hydride
    // Source: openclatura | tests/unit/test_polyoxide_chains.py::test_substituted_oxygen_chains_use_the_oxidane_parent
    // Source: openclatura | tests/unit/test_audit_reconstruction.py::test_self_audit_confirms_modelled_names
    // Source: openclatura | tests/unit/test_chalcogen_functional_groups.py::test_ether_analogues_reuse_generic_heteroatom_subgraph_naming
    // Source: openclatura | tests/unit/test_all_chalcogen_acyl_regressions.py::test_reported_name_regressions
    // Source: openclatura | tests/unit/test_analysis.py::test_nitrogen_chain_roles_make_cyclic_hydrazines_prefixes
    // Source: openclatura | tests/unit/test_analysis.py::test_component_shortcuts_are_final_binding_audited
    // Source: openclatura | tests/unit/test_analysis.py::test_anhydride_shortcut_binds_acid_halves_and_core_separately
    // Source: openclatura | tests/unit/test_analysis.py::test_a_zwitterion_keeps_the_anion_as_the_suffix
    // Source: openclatura | tests/roundtrip/test_oxyacids.py::test_roundtrip
    // Source: openclatura | tests/unit/test_locant_elision.py::test_nonconstitutional_and_unsupported_parent_locants_are_unchanged
    // Source: openclatura | tests/roundtrip/test_polycycles.py::test_roundtrip
    // Source: openclatura | tests/roundtrip/test_role_models.py::test_roundtrip
    // Source: openclatura | tests/roundtrip/test_principal_group_breadth.py::test_roundtrip
    // Source: openclatura | tests/exploratory/test_group14_ether_ligand_parent_selection.py::test_ring_ligand_does_not_outrank_the_mononuclear_hydride_parent
    // Source: openclatura | tests/exploratory/test_group14_ether_ligand_parent_selection.py::test_a_genuine_principal_group_elsewhere_still_blocks_the_hydride_parent
    // Source: openclatura | tests/data/parser_xml_resources/arylGroups.json::<fixture data>
    // Source: openclatura | tests/data/parser_xml_resources/nonCarboxylicAcids.json::<fixture data>
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn allow_list_status_matrix() {
        use nazuchem_ir::{NameReason, NameStatus};
        for (smiles, expected_name, expected_status) in [
            ("CCO", "ethanol", NameStatus::Preferred),
            ("CC(=O)O", "acetic acid", NameStatus::Preferred),
            ("Oc1ccccc1", "phenol", NameStatus::Preferred),
            ("OC1CCCCC1", "cyclohexanol", NameStatus::Preferred),
            ("c1ccncc1", "pyridine", NameStatus::Preferred),
            ("c1ccc2ccccc2c1", "naphthalene", NameStatus::Preferred),
            (
                "c1ccc2c(c1)C=CC=CO2",
                "1-benzoxepine",
                NameStatus::Preferred,
            ),
            (
                "C1CCC2CCCCC2C1",
                "decahydronaphthalene",
                NameStatus::Preferred,
            ),
            (
                "c1ccc(cc1)-c1ccccc1",
                "1,1'-biphenyl",
                NameStatus::Preferred,
            ),
            (
                "C1CC2CCC1C2",
                "bicyclo[2.2.1]heptane",
                NameStatus::Preferred,
            ),
            (
                "O[C@@H]1C[C@H]2CC[C@@H]1C2",
                "(1R,2R,4S)-bicyclo[2.2.1]heptan-2-ol",
                NameStatus::Preferred,
            ),
            (
                "C1CCC2(CC1)CCCC2",
                "spiro[4.5]decane",
                NameStatus::Preferred,
            ),
            (
                "C1CCC2(CC1)OCCO2",
                "1,4-dioxaspiro[4.5]decane",
                NameStatus::Preferred,
            ),
            ("CCOC(C)=O", "ethyl acetate", NameStatus::Preferred),
            ("CC(=O)OC(C)=O", "acetic anhydride", NameStatus::Preferred),
            // P-65.7.6.1 example (issue #110).
            (
                "CC(=O)OC(=O)CCC(=O)OC(C)=O",
                "diacetic butanedioic dianhydride",
                NameStatus::Preferred,
            ),
            // P-65.7.6.4.2 example: two unlike linkages stay unverified.
            (
                "CC(=O)OC(=O)CCCC(=O)SC(C)=O",
                "acetic 5-(acetylsulfanyl)-5-oxopentanoic anhydride",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "c1ccc(cc1)[Si](c1ccccc1)(c1ccccc1)c1ccccc1",
                "tetraphenylsilane",
                NameStatus::Preferred,
            ),
            (
                "CP(C)(C)=O",
                "trimethyl-lambda5-phosphanone",
                NameStatus::Preferred,
            ),
            (
                "C[P+](C)(C)[O-]",
                "trimethyl-lambda5-phosphanone",
                NameStatus::Preferred,
            ),
            (
                "o1ccc2oc3ccoc3c12",
                "difuro[3,2-b:2',3'-d]furan",
                NameStatus::Unverified(NameReason::FusionMultiparent),
            ),
            (
                "O1CCC2CCCCC2C1",
                "octahydro-1H-2-benzopyran",
                NameStatus::Preferred,
            ),
            (
                "O=C(Oc1ccc(Oc2ccccc2)cc1)Oc1ccc(Oc2ccccc2)cc1",
                "bis(4-phenoxyphenyl) carbonate",
                NameStatus::Preferred,
            ),
            ("CCO[N+](=O)[O-]", "ethyl nitrate", NameStatus::Preferred),
            ("[Na+].CC(=O)[O-]", "sodium acetate", NameStatus::Preferred),
            (
                "O=C([O-])c1cc[n-]c1",
                "1H-pyrrol-1-ide-3-carboxylate",
                NameStatus::Preferred,
            ),
            (
                "C1CC2C3CCC1C2C3",
                "tricyclo[4.3.0.0^{3,9}]nonane",
                NameStatus::Unverified(NameReason::UnauditedParent),
            ),
            (
                "C1CCCCCOCCCCC1",
                "1-oxacyclododecane",
                NameStatus::Unverified(NameReason::UnauditedParent),
            ),
            (
                "COOC",
                "1,2-dimethyldioxidane",
                NameStatus::NotPreferred(NameReason::HydrideParent),
            ),
            (
                "C[Te]C",
                "dimethyltellane",
                NameStatus::NotPreferred(NameReason::HydrideParent),
            ),
            (
                "NN",
                "hydrazine",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "O[N+](=O)[O-]",
                "nitric acid",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "[SiH4]",
                "silane",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            ("CCB(O)O", "ethylboronic acid", NameStatus::Preferred),
            (
                "C12(CC3CCC1C3)C1CCC2CC1",
                "2,7'-spirobi[bicyclo[2.2.1]heptane]",
                NameStatus::Preferred,
            ),
            (
                "CC(=O)[NH-]",
                "acetamide",
                NameStatus::Unverified(NameReason::AnionParent),
            ),
            (
                "c1ccccc1C1=CC2(C3=CC=CC=C13)C=CC1=CC=CC=C12",
                "3-phenyl-1,1'-spirobi[indene]",
                NameStatus::Preferred,
            ),
            (
                "COP(=O)(OC)OCC(=O)O",
                "carboxymethyl dimethyl phosphate",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "CC1=CC2(C3=CC=CC=C13)C=CC1=CC=CC=C12",
                "3-methyl-1,1'-spirobi[indene]",
                NameStatus::Preferred,
            ),
            (
                "[O-]C(=O)c1cc(C(=O)[O-])ccc1Cc1cc(S(=O)(=O)[O-])[n-]c1",
                "4-{[2,4-bis(carboxylato)phenyl]methyl}-1H-pyrrol-1-ide-2-sulfonate",
                NameStatus::Preferred,
            ),
            (
                "C1C=Cc2cc3sccc3n12",
                "5H-thieno[2,3-b]pyrrolizine",
                NameStatus::Preferred,
            ),
            (
                "C[SH2]C",
                "dimethyl-lambda4-sulfane",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "C1=CC2(C3=CC=CC=C13)C=CC1=CC=CC=C12",
                "1,1'-spirobi[indene]",
                NameStatus::Preferred,
            ),
            (
                "CC(=O)OOC(C)=O",
                "acetic peroxyanhydride",
                NameStatus::Preferred,
            ),
            (
                "CC(=O)[Se]C(C)=O",
                "acetic selenoanhydride",
                NameStatus::Preferred,
            ),
            (
                "CN(C)C(=O)OP(=O)(O)O",
                "dimethylcarbamoyl dihydrogen phosphate",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "[Na+].[Cl-]",
                "sodium chloride",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            ("C[O-]", "methoxide", NameStatus::Preferred),
            (
                "[Na+].c1ccc2[n-]cc(C(=O)O)c2c1",
                "sodium 3-carboxy-1H-indol-1-ide",
                NameStatus::Preferred,
            ),
            (
                "[O-]C(=O)c1ccccc1-c1[n-]cc(C(=O)[O-])c1",
                "5-(2-carboxylatophenyl)-1H-pyrrol-1-ide-3-carboxylate",
                NameStatus::Preferred,
            ),
            // P-63.2.4.2: the multiplicative name is generated.
            (
                "c1ccc(Oc2ccccc2)cc1",
                "1,1'-oxydibenzene",
                NameStatus::Preferred,
            ),
            (
                "Oc1ccc(Oc2ccc(Oc3ccc(Oc4ccccc4)cc3)cc2)cc1",
                "4-[4-(4-phenoxyphenoxy)phenoxy]phenol",
                NameStatus::NotPreferred(NameReason::LinearPhane),
            ),
            (
                "C1C2CC3CC1CC2C3",
                "octahydro-2,5-methanopentalene",
                NameStatus::Unverified(NameReason::BridgedFusedParent),
            ),
            (
                "CC(=O)OOC(=O)CC",
                "acetic propanoic peroxyanhydride",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "CC(=O)OC(=O)CC",
                "acetic propanoic anhydride",
                NameStatus::Preferred,
            ),
            (
                "CP(C)Cl",
                "chlorodi(methyl)phosphane",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "COP(C)C",
                "methoxydi(methyl)phosphane",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "C[Al](C)C",
                "trimethylalumane",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "CC[Sn](CC)(CC)CC",
                "tetraethylstannane",
                NameStatus::Preferred,
            ),
            (
                "c1ccc(cc1)P(c1ccccc1)c1ccccc1",
                "triphenylphosphane",
                NameStatus::Preferred,
            ),
            (
                "CC(=S)OP(=O)(OC)OC",
                "ethanethioyl dimethyl phosphate",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            // P-45.5: 'bis' precedes 'methoxy' between the tied ethane parents.
            (
                "COCCOCCOCCOC",
                "1,2-bis(2-methoxyethoxy)ethane",
                NameStatus::Unverified(NameReason::UnauditedParent),
            ),
            (
                "CC(=O)NCCOCCOCCOCCOC",
                "N-(2-{2-[2-(2-methoxyethoxy)ethoxy]ethoxy}ethyl)acetamide",
                NameStatus::Unverified(NameReason::UnauditedParent),
            ),
            ("CCOCCOCC", "1,2-diethoxyethane", NameStatus::Preferred),
            (
                "CC(=O)NCC(=O)NCC(=O)NCC(=O)NCC(=O)O",
                "{2-[2-(2-acetamidoacetamido)acetamido]acetamido}acetic acid",
                NameStatus::Unverified(NameReason::UnauditedParent),
            ),
            (
                "C[N+](C)(C)CC(=O)[O-]",
                // P-74.1.3 (2020 correction)
                "(N,N-dimethylmethanaminiumyl)acetate",
                NameStatus::Preferred,
            ),
            (
                "C1CC2C3CCC1C2C3C4CC5C6CCC4C5C6",
                "2,7'-bi(tricyclo[4.3.0.0^{3,9}]nonane)",
                NameStatus::Unverified(NameReason::UnauditedParent),
            ),
            (
                "[O-]C(=O)c1ccccc1Cc1cc[n-]c1",
                "3-[(2-carboxylatophenyl)methyl]-1H-pyrrol-1-ide",
                NameStatus::Preferred,
            ),
            (
                "[O-]C(=O)c1ccccc1[CH2-]",
                "(2-carboxylatophenyl)methanide",
                NameStatus::Preferred,
            ),
            (
                "[Na+].[Na+].[O-]c1ccccc1Cc1ccccc1C(=O)[O-]",
                "disodium 2-[(2-oxidophenyl)methyl]benzoate",
                NameStatus::Preferred,
            ),
            // P-52.2.5.3 Example 2: the multiplicative name is generated.
            (
                "c1ccccc1Oc2ccc(Oc3ccccc3)cc2",
                "1,1'-[1,4-phenylenebis(oxy)]dibenzene",
                NameStatus::Preferred,
            ),
            (
                "[Na+]",
                "sodium(1+)",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "OCC[NH3+]",
                "2-hydroxyethan-1-aminium",
                NameStatus::Preferred,
            ),
            ("[O-]c1ccccc1", "phenoxide", NameStatus::Preferred),
        ] {
            let name = name_smiles(smiles, &Default::default())
                .unwrap_or_else(|error| panic!("{smiles}: {error}"));
            assert_eq!(name.text, expected_name, "{smiles}");
            assert_eq!(name.status, expected_status, "{smiles}: {}", name.text);
        }
    }

    #[test]
    fn polyamine_n_locants_cite_the_parent_locant() {
        // P-62.2.4.1.1, P-62.2.4.1.2: the parent locant is a superscript on N,
        // as in the PINs N^1-ethyl-N^3-methylpropane-1,3-diamine and
        // N^4,2-dimethylpentane-2,4-diamine, never N or N'.
        for (smiles, expected) in [
            ("CNCCCNCC", "N^1-ethyl-N^3-methylpropane-1,3-diamine"),
            ("CC(N)(C)CC(NC)C", "N^4,2-dimethylpentane-2,4-diamine"),
            ("CN(C)CCCNC", "N^1,N^1,N^3-trimethylpropane-1,3-diamine"),
            (
                "CCN(CC)CCCNC",
                "N^1,N^1-diethyl-N^3-methylpropane-1,3-diamine",
            ),
            (
                "CNCCCCCCCCCCCCNCC",
                "N^1-ethyl-N^12-methyldodecane-1,12-diamine",
            ),
            (
                "C1CCN(CC1)CCNC2=C(C(=NC=N2)NC3=NC4=CC=CC=C4S3)[N+](=O)[O-]",
                "N^4-(1,3-benzothiazol-2-yl)-5-nitro-N^6-[2-(piperidin-1-yl)ethyl]pyrimidine-4,6-diamine",
            ),
            // Fusion-letter and composite parent locants keep their hyphen and order.
            (
                "CNC12CCCCC1(NCC)CCCC2",
                "N^4a-ethyl-N^8a-methyldecahydronaphthalene-4a,8a-diamine",
            ),
            (
                "C1(C2CCCCCC(C3CCCCCCCCCCC3)CCCCC2)C(NC)CCCCCCCC(NC)CC1",
                "N^1^2,N^1^10-dimethyl[1^1,2^1:2^7,3^1-tercyclododecane]-1^2,1^10-diamine",
            ),
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, nazuchem_ir::NameStatus::Preferred, "{smiles}");
            let parsed = crate::test_support::parse_strict(&name.text).unwrap();
            assert_eq!(
                crate::test_support::canonical_smiles(smiles).unwrap(),
                crate::test_support::canonical_smiles(&parsed).unwrap(),
                "{expected}",
            );
        }
    }

    #[test]
    fn geminal_amine_n_locants_take_primes() {
        // P-62.2.4.1.2, P-16.9.3: geminal amine nitrogens share a superscript
        // parent locant and are told apart by primes, as in the PINs
        // N^3-ethyl-N'^3-methylhexane-3,3-diamine and
        // N^3-ethyl-N^1,N'^3-dimethylhexane-1,3,3,6-tetramine. Unprimed is
        // lower than primed (P-14.3.5), so it goes to the nitrogen with more
        // substituents, then to the one cited first (P-14.4 (f), (g)).
        for (smiles, expected) in [
            ("CCC(NCC)(NC)CCC", "N^3-ethyl-N'^3-methylhexane-3,3-diamine"),
            ("CCC(NC)(NCC)CCC", "N^3-ethyl-N'^3-methylhexane-3,3-diamine"),
            (
                "CNCCC(NCC)(NC)CCCN",
                "N^3-ethyl-N^1,N'^3-dimethylhexane-1,3,3,6-tetramine",
            ),
            ("CCC(N)(NC)CCC", "N^3-methylhexane-3,3-diamine"),
            (
                "CCC(NCC)(N(C)C)CCC",
                "N'^3-ethyl-N^3,N^3-dimethylhexane-3,3-diamine",
            ),
            (
                "CCC(N(C)CC)(NCCC)CCC",
                "N^3-ethyl-N^3-methyl-N'^3-propylhexane-3,3-diamine",
            ),
            (
                "CCNC1(NC)CCCCC1",
                "N^1-ethyl-N'^1-methylcyclohexane-1,1-diamine",
            ),
            (
                "CC(NC)(NCC)N(C)CCC",
                "N'^1-ethyl-N^1,N''^1-dimethyl-N^1-propylethane-1,1,1-triamine",
            ),
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, nazuchem_ir::NameStatus::Preferred, "{smiles}");
            let parsed = crate::test_support::parse_strict(&name.text).unwrap();
            assert_eq!(
                crate::test_support::canonical_smiles(smiles).unwrap(),
                crate::test_support::canonical_smiles(&parsed).unwrap(),
                "{expected}",
            );
        }
    }

    #[test]
    fn axial_general_names_leave_nuclide_groups_alone() {
        assert_eq!(
            super::general_axial_name("(31P)phosphanium (2M)-2-(2-hydroxyphenyl)benzen-1-olate")
                .as_deref(),
            Some("(31P)phosphanium (2Ra)-2-(2-hydroxyphenyl)benzen-1-olate")
        );
        assert_eq!(super::general_axial_name("(31P)phosphane"), None);
        assert_eq!(
            super::general_axial_name("4-[(1P,2R)-x]-(3M)-y").as_deref(),
            Some("4-[(1Sa,2R)-x]-(3Ra)-y")
        );
    }

    #[test]
    fn axial_general_alternatives_cover_ordinary_unknown_isotope_and_salt_paths() {
        use nazuchem_ir::{NameReason, NameStatus};
        for (marker, descriptor, general) in [(900_004, 'M', "Ra"), (900_008, 'P', "Sa")] {
            for (input, unknown) in [
                (
                    format!("Oc1ccc2ccccc2[c:{marker}]1-[c:{marker}]1c(O)ccc2ccccc12"),
                    false,
                ),
                (
                    format!("O[13c]1ccc2ccccc2[c:{marker}]1-[c:{marker}]1c(O)ccc2ccccc12"),
                    false,
                ),
                (
                    format!("Oc1ccc(C[C:900001](F)(Cl)Br)c[c:{marker}]1-[c:{marker}]1ccccc1O"),
                    true,
                ),
                (
                    format!("Oc1cccc[c:{marker}]1-[c:{marker}]1ccccc1[O-].[Na+]"),
                    false,
                ),
            ] {
                let name = name_smiles(&input, &Default::default())
                    .unwrap_or_else(|e| panic!("{input}: {e}"));
                assert!(name.text.contains(descriptor), "{input}: {}", name.text);
                let expected = name.text.replace(descriptor, general);
                let alternative = name
                    .alternatives
                    .iter()
                    .find(|a| a.name == expected)
                    .unwrap_or_else(|| {
                        panic!("{input}: missing {expected}: {:?}", name.alternatives)
                    });
                assert_eq!(
                    alternative.status,
                    NameStatus::NotPreferred(NameReason::AxialAlternative)
                );
                assert!(alternative.status.detail().contains("P-93.5.7.1"));
                if unknown {
                    assert_eq!(name.omitted_units, 1);
                    assert!(
                        name.alternatives
                            .iter()
                            .any(|a| a.name.contains('Ξ') && a.name.contains(general))
                    );
                }
            }
        }
    }

    #[test]
    fn complete_senior_parents_need_multiplicative_pin_status_and_round_trip() {
        use nazuchem_ir::{NameReason, NameStatus};
        // P-16.5.4 Fig. 1.3(f) and P-15.3.2.1: the first central group is
        // asymmetric; the second has a cyclic assembly between equal arms.
        for (smiles, expected) in [
            // P-45.5: 'oxopropanyl' precedes 'oxopropyl' between the tied rings.
            (
                "N#Cc1ccc(-c2ccc(OCCCCCOC(=O)C(Cc3ccc(C(=O)O)cc3)c3ccc(C(=O)O)cc3)cc2)cc1",
                "4-[3-(4-carboxyphenyl)-1-({5-[(4'-cyano[1,1'-biphenyl]-4-yl)oxy]pentyl}oxy)-1-oxopropan-2-yl]benzoic acid",
            ),
            (
                "O=C(O)CCCCCCCCCOc1ccc(-c2ccc(OCCCCCCCCCC(=O)O)cc2)cc1",
                "10-({4'-[(9-carboxynonyl)oxy][1,1'-biphenyl]-4-yl}oxy)decanoic acid",
            ),
            // P-15.3 / P-51.3.1: a more substituted core does not remove
            // identical terminal parents from multiplicative eligibility.
            (
                "c1ccccc1Oc2ccc(Oc3ccccc3)c(C)c2",
                "2-methyl-1,4-diphenoxybenzene",
            ),
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(
                name.status,
                NameStatus::NotPreferred(NameReason::Multiplicative),
                "{smiles}"
            );
            let parsed = crate::test_support::parse_strict(&name.text).unwrap();
            assert_eq!(
                crate::test_support::canonical_smiles(smiles).unwrap(),
                crate::test_support::canonical_smiles(&parsed).unwrap(),
                "{}",
                name.text,
            );
        }
    }

    #[test]
    fn multiplicative_pin_accepts_simple_cores_and_identical_concatenated_arms() {
        for smiles in [
            "c1ccccc1Oc1ccccc1",
            "O=S(c1ccccc1)c1ccccc1",
            "O=[Se](c1ccccc1)c1ccccc1",
            "O=[Te](=O)(c1ccccc1)c1ccccc1",
            "S(CC(=O)O)CC(=O)O",
            // Principal acids outrank the central phosphane (P-51.3.1).
            "P(CC(=O)O)(CC(=O)O)CC(=O)O",
            "CC(c1ccccc1)Cc1ccccc1",
            "O=C(O)c1ccc(OCCOc2ccc(C(=O)O)cc2)cc1",
        ] {
            let mol = super::graph_io::read_smiles(smiles, &Default::default()).unwrap();
            let atoms = mol.atoms.keys().copied().collect();
            let name =
                super::component_namer::name_component_with_citation(&mol, &atoms, None, None)
                    .unwrap();
            assert!(
                super::has_multiplicative_pin(&mol, &name.derivation),
                "{smiles}"
            );
        }
    }

    #[test]
    // Provenance: fn tests::senior_structural_hydrides_keep_preferred_names | adapted
    // Source: ReShiki | src/naming/openclatura/mod.rs::tests::senior_structural_hydrides_keep_preferred_names
    // Source: openclatura | tests/exploratory/test_group14_ether_ligand_parent_selection.py::test_ring_ligand_does_not_outrank_the_mononuclear_hydride_parent
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn senior_structural_hydrides_keep_preferred_names() {
        for (smiles, expected) in [
            ("P(c1ccccc1)(c1ccccc1)c1ccccc1", "triphenylphosphane"),
            (
                "[Si](c1ccccc1)(c1ccccc1)(c1ccccc1)c1ccccc1",
                "tetraphenylsilane",
            ),
            ("B(c1ccccc1)(c1ccccc1)c1ccccc1", "triphenylborane"),
            ("N(c1ccccc1)c1ccccc1", "N-phenylaniline"),
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, nazuchem_ir::NameStatus::Preferred, "{smiles}");
        }
    }

    #[test]
    // Provenance: fn tests::multiplicative_pin_requires_identical_complete_rooted_senior_units | adapted
    // Source: ReShiki | src/naming/openclatura/mod.rs::tests::multiplicative_pin_requires_identical_complete_rooted_senior_units
    // Source: openclatura | tests/unit/test_analysis.py::test_anhydride_shortcut_binds_acid_halves_and_core_separately
    // Source: openclatura | tests/unit/test_polyoxide_chains.py::test_lone_oxygen_atoms_are_untouched
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn multiplicative_pin_requires_identical_complete_rooted_senior_units() {
        for smiles in [
            // Equal parents, different attachment locants (P-51.3.1(3)).
            "O=C(O)c1ccc(Cc2cccc(C(=O)O)c2)cc1",
            // Different substituted parent units and configurations.
            "Brc1cccc(Cc2cccc(Cl)c2)c1",
            "CC[C@@H](C)c1ccc(Sc2ccc(cc2)[C@@H](C)CC)cc1",
            r"O=C(O)/C=C/COC/C=C\C(=O)O",
            // Unequal concatenated arms (P-15.3.4.1.1, P-15.3.4.1.3).
            "Oc1ccc(CNc2ccc(O)cc2)cc1",
            "O=C(O)c1ccc(CC(Cl)OCCc2ccc(C(=O)O)cc2)cc1",
            // A more senior chain contains the principal group.
            "O=C(O)CC(c1ccccc1)c1ccccc1",
            // Unfunctionalized acyclic hydrocarbons and anhydrides.
            "CCOCC",
            "CC(=O)OC(C)=O",
            // Direct ring junctions are assemblies, including their suffixes.
            "Oc1ccccc1-c1ccccc1O",
            "O=C1CCC(C2CCC(=O)CC2)CC1",
            "c1ccc(-c2ccc(-c3ccccc3)cc2)cc1",
        ] {
            let mol = super::graph_io::read_smiles(smiles, &Default::default()).unwrap();
            let atoms = mol.atoms.keys().copied().collect();
            let name =
                super::component_namer::name_component_with_citation(&mol, &atoms, None, None)
                    .unwrap();
            assert!(
                !super::has_multiplicative_pin(&mol, &name.derivation),
                "{smiles}"
            );
        }
    }

    #[test]
    fn unequal_or_junior_multiplicative_units_keep_preferred_names_and_round_trip() {
        for (smiles, expected) in [
            ("Oc1ccc(Oc2ccccc2)cc1", "4-phenoxyphenol"),
            ("Cc1ccc(Oc2ccccc2)cc1", "1-methyl-4-phenoxybenzene"),
            ("c1ccccc1Oc2ccc(Oc3ccccc3)nc2", "2,5-diphenoxypyridine"),
            ("C1CCC(CC1)c1ccccc1", "cyclohexylbenzene"),
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, nazuchem_ir::NameStatus::Preferred, "{smiles}");
            let parsed = crate::test_support::parse_strict(&name.text).unwrap();
            assert_eq!(
                crate::test_support::canonical_smiles(smiles).unwrap(),
                crate::test_support::canonical_smiles(&parsed).unwrap(),
                "{}",
                name.text
            );
        }
    }

    #[test]
    // Provenance: fn tests::anhydride_names_need_acid_halves | adapted
    // Source: ReShiki | src/naming/openclatura/mod.rs::tests::anhydride_names_need_acid_halves
    // Source: openclatura | tests/unit/test_analysis.py::test_anhydride_shortcut_binds_acid_halves_and_core_separately
    // Source: openclatura | tests/unit/test_chalcogen_citation_compatibility.py::test_mixed_site_anhydrides_keep_both_acyl_ends
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn anhydride_names_need_acid_halves() {
        // Upstream names this "2-(((hydroxy)(thioxo)methyl)sulfanyl)acetonitrile
        // 4-methylbenzoic thioanhydride", which OPSIN reads as two compounds.
        let name = name_smiles("Cc1ccc(C(=O)SC(=S)SCC#N)cc1", &Default::default());
        assert!(
            !matches!(&name, Ok(name) if name.text.contains("anhydride")),
            "{name:?}"
        );
    }

    #[test]
    // Provenance: fn tests::enclosing_marks_nest_by_what_they_enclose | adapted
    // Source: ReShiki | src/naming/openclatura/mod.rs::tests::enclosing_marks_nest_by_what_they_enclose
    // Source: openclatura | tests/unit/test_analysis.py::test_public_api_golden_names
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn enclosing_marks_nest_by_what_they_enclose() {
        for (smiles, expected) in [
            (
                "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
                "2-[4-(2-methylpropyl)phenyl]propanoic acid",
            ),
            // P-51.3.1: the multiplicative name is generated.
            (
                "OCCOCCOCCOCCO",
                "2,2'-[oxybis(ethane-2,1-diyloxy)]di(ethan-1-ol)",
            ),
            (
                "OCCOCCOCCOCCOCCO",
                "2-(2-{2-[2-(2-hydroxyethoxy)ethoxy]ethoxy}ethoxy)ethan-1-ol",
            ),
            (
                "CN(CC(=O)O)Nc1cnc(CO)cn1",
                "({[5-(hydroxymethyl)pyrazin-2-yl]amino}(methyl)amino)acetic acid",
            ),
            (
                "C[S@@](=O)CSCNN",
                "1-[({[(R)-methanesulfinyl]methyl}sulfanyl)methyl]hydrazine",
            ),
            // Stereodescriptors and compound locants count as a level; fusion
            // and von Baeyer brackets do not (P-16.5.4.1.2, P-16.5.4.1.3).
            ("C/C=C/Cc1ccc(O)cc1", "4-[(2E)-but-2-en-1-yl]phenol"),
            (
                "CN(CCO)[C@H]1CCCN1",
                "2-{methyl[(2S)-pyrrolidin-2-yl]amino}ethan-1-ol",
            ),
            (
                "OC(=O)CC1CC2CCC1C2",
                "(bicyclo[2.2.1]heptan-2-yl)acetic acid",
            ),
            (
                "CC(=O)NC12C=C1BNC2",
                "N-[3-aza-2-borabicyclo[3.1.0]hex-1(6)-en-5-yl]acetamide",
            ),
            (
                "OC(=O)Cc1ccc(cc1)-c1c[nH]c2ncccc12",
                "[4-(1H-pyrrolo[2,3-b]pyridin-3-yl)phenyl]acetic acid",
            ),
            (
                "CC(=O)N[C@@H]1CCN2C=CN[C@H]12",
                "N-[(7R,7aS)-5,6,7,7a-tetrahydro-1H-pyrrolo[1,2-a]imidazol-7-yl]acetamide",
            ),
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            if expected.ends_with("acetic acid") {
                let parsed = crate::test_support::parse_strict(&name.text).unwrap();
                assert_eq!(
                    crate::test_support::canonical_smiles(smiles).unwrap(),
                    crate::test_support::canonical_smiles(&parsed).unwrap(),
                    "{}",
                    name.text
                );
            }
        }
    }

    /// The P-16.5.4 examples, given as the namer writes them before nesting.
    #[test]
    fn enclosing_marks_follow_the_blue_book_examples() {
        for (enclosed, expected) in [
            (
                "4,4'-(3-((5-((4'-cyano[1,1'-biphenyl]-4-yl)oxy)pentyl)oxy)-3-oxopropane-1,2-diyl)dibenzoic acid",
                "4,4'-[3-({5-[(4'-cyano[1,1'-biphenyl]-4-yl)oxy]pentyl}oxy)-3-oxopropane-1,2-diyl]dibenzoic acid",
            ),
            (
                "1-(3,4-dihydroquinolin-1(2H)-yl)ethan-1-one",
                "1-(3,4-dihydroquinolin-1(2H)-yl)ethan-1-one",
            ),
            (
                "10,10'-([1,1'-biphenyl]-4,4'-diylbis(oxy))di(decanoic acid)",
                "10,10'-[[1,1'-biphenyl]-4,4'-diylbis(oxy)]di(decanoic acid)",
            ),
            (
                "4-(4-([1,1'-bi(cyclohexan)]-4-yl)phenyl)-4'-phenyl-1,1'-bi(cyclohexane)",
                "4-[4-([1,1'-bi(cyclohexan)]-4-yl)phenyl]-4'-phenyl-1,1'-bi(cyclohexane)",
            ),
            (
                "1-(tricyclo[3.3.1.1^{3,7}]decan-1-yl)ethan-1-one",
                "1-(tricyclo[3.3.1.1^{3,7}]decan-1-yl)ethan-1-one",
            ),
            (
                "2-(bicyclo[6.6.1]pentadeca-8(15)-en-1-yl)ethan-1-ol",
                "2-[bicyclo[6.6.1]pentadeca-8(15)-en-1-yl]ethan-1-ol",
            ),
            // P-16.5.4.1.5: `((2S)-` would put two parentheses together.
            (
                "(3S)-2-((2S)-2-(((2S)-1-ethoxy-1-oxo-4-phenylbutan-2-yl)amino)propanoyl)-1,2,3,4-tetrahydroisoquinoline-3-carboxylic acid",
                "(3S)-2-[(2S)-2-{[(2S)-1-ethoxy-1-oxo-4-phenylbutan-2-yl]amino}propanoyl]-1,2,3,4-tetrahydroisoquinoline-3-carboxylic acid",
            ),
        ] {
            assert_eq!(nest_enclosing_marks(enclosed), expected);
        }
    }

    #[test]
    // Provenance: fn tests::controls_keep_their_names_and_preference_status | adapted
    // Source: ReShiki | src/naming/openclatura/mod.rs::tests::controls_keep_their_names_and_preference_status
    // Source: openclatura | tests/fixtures/diverse_corpus.golden.json::<fixture data>
    // Source: openclatura | tests/data/parser_xml_resources/functionalTerms.json::<fixture data>
    // Source: openclatura | tests/unit/test_polyoxide_chains.py::test_the_other_heteroatom_chains_are_unchanged
    // Source: openclatura | tests/unit/test_audit_reconstruction.py::test_self_audit_confirms_modelled_names
    // Source: openclatura | tests/data/fusion_pin_cases.json::<fixture data>
    // Source: openclatura | tests/fixtures/diverse_corpus.csv::<fixture data>
    // Source: openclatura | tests/unit/test_analysis.py::test_sulfur_ylide_resonance_compare_accepts_lambda_fallback_graph
    // Source: openclatura | tests/unit/test_analysis.py::test_anhydride_shortcut_binds_acid_halves_and_core_separately
    // Source: openclatura | tests/unit/test_polyoxide_chains.py::test_oxygen_chains_do_not_displace_a_senior_heteroatom_chain
    // Source: openclatura | tests/roundtrip/test_alkenes_alkynes.py::test_roundtrip
    // Source: openclatura | tests/roundtrip/test_replacement_parents.py::test_roundtrip
    // Source: openclatura | tests/unit/test_chalcogen_structural_roles.py::test_acyl_pseudohalide_owns_complete_leaving_unit
    // Source: openclatura | tests/unit/test_chalcogen_functional_groups.py::test_nitrile_chalcogenides_use_zwitterion_priority_and_class_names
    // Source: openclatura | tests/unit/test_chalcogen_functional_groups.py::test_dynamic_chalcogen_rules_compose_with_existing_parent_and_modifier_paths
    // Source: openclatura | tests/data/parser_xml_resources/simpleGroups.json::<fixture data>
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn controls_keep_their_names_and_preference_status() {
        use nazuchem_ir::NameStatus;
        for (smiles, expected) in [
            ("CCO", "ethanol"),
            ("O=C(O)c1ccccc1", "benzoic acid"),
            (
                "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
                "2-[4-(2-methylpropyl)phenyl]propanoic acid",
            ),
            ("c1ccc2[nH]ccc2c1", "1H-indole"),
            // Double bonds take the lower locant once multiple bonds tie
            // (P-31.1.1.1), whichever end the input starts from.
            ("C=CCC#C", "pent-1-en-4-yne"),
            ("C#CCC=C", "pent-1-en-4-yne"),
            ("C#CCCC=C", "hex-1-en-5-yne"),
            ("C#CCC(C)(C)CC=C", "4,4-dimethylhept-1-en-6-yne"),
            ("SS", "disulfane"),
            ("SSS", "trisulfane"),
            ("CC(=O)C#N", "acetyl cyanide"),
            ("COC(=O)Nc1ccccc1", "methyl phenylcarbamate"),
            ("CCC#[N+][O-]", "propanenitrile oxide"),
            // Substitutive sulfoxide (P-63.6) and -one (P-68.3.2.3.1) names.
            ("CS(C)=O", "(methanesulfinyl)methane"),
            ("CC(C)S(C)=O", "2-(methanesulfinyl)propane"),
            ("CP(C)(C)=O", "trimethyl-lambda5-phosphanone"),
            ("C[P+](C)(C)[O-]", "trimethyl-lambda5-phosphanone"),
            ("C[As](C)(C)=O", "trimethyl-lambda5-arsanone"),
            ("CP=O", "methylphosphanone"),
            ("c1ccccc1[As]=S", "phenylarsanethione"),
            ("C[As](C)(C)=[Te]", "trimethyl-lambda5-arsanetellone"),
            ("CP(C)(C)=S", "trimethyl-lambda5-phosphanethione"),
            // Monosubstituted ethane and ethene cite no locants (P-14.3.4).
            ("CCS(C)=O", "(methanesulfinyl)ethane"),
            ("CCSC", "(methylsulfanyl)ethane"),
            ("C=CSS", "ethenedithioperoxol"),
            ("C=CS(=O)(=O)O", "ethenesulfonic acid"),
            (
                "COS(=O)(=O)c1ccc(C(=O)O)cc1",
                "4-(methoxysulfonyl)benzoic acid",
            ),
            // P-16.5.1.1: the compound prefix methoxysulfonyl is enclosed.
            ("COS(=O)(=O)NCC(=O)O", "[(methoxysulfonyl)amino]acetic acid"),
            // Sulfonyl and sulfinyl prefixes are named from the acid (P-65.3.2.2.2).
            ("CS(=O)(=O)c1ccccc1", "(methanesulfonyl)benzene"),
            ("C[S@](=O)c1ccccc1", "[(S)-methanesulfinyl]benzene"),
            ("CS(=O)(=O)N1CCCC1", "1-(methanesulfonyl)pyrrolidine"),
            (
                "O=C(O)c1ccc(cc1)S(=O)(=O)N1CCN(C)CC1",
                "4-(4-methylpiperazine-1-sulfonyl)benzoic acid",
            ),
            (
                "O=C(O)C1CCN(C1)S(=O)(=O)c1ccc(C)cc1",
                "1-(4-methylbenzene-1-sulfonyl)pyrrolidine-3-carboxylic acid",
            ),
            (
                "CC(NS(=O)(=O)c1ccc(C)cc1)C(=O)O",
                "2-(4-methylbenzene-1-sulfonamido)propanoic acid",
            ),
            (
                "CS(=O)(=O)OCCC(=O)O",
                "3-[(methanesulfonyl)oxy]propanoic acid",
            ),
            ("OC(=O)CS(=O)(=O)c1ccccc1", "(benzenesulfonyl)acetic acid"),
            ("OS(=O)(=O)ON=C1CCCCC1", "N-(sulfooxy)cyclohexanimine"),
            // Substituted anisoles take systematic PINs (P-63.2.3).
            ("COc1ccc(Cl)cc1", "1-chloro-4-methoxybenzene"),
            ("COc1ccccc1OC", "1,2-dimethoxybenzene"),
            ("O=CCl", "formyl chloride"),
            ("CNC=O", "N-methylformamide"),
            ("COC=O", "methyl formate"),
            // Chalcogen analogues of oxo and formyl (P-64.6.1, P-66.6.3),
            // multiplied with bis to avoid ambiguity (P-16.3.6).
            ("O=C(O)CC(=S)C", "3-sulfanylidenebutanoic acid"),
            ("O=C1CSC(=S)N1", "2-sulfanylidene-1,3-thiazolidin-4-one"),
            ("CC(=[Te])CC(=[Se])C", "4-tellanylidenepentane-2-selone"),
            (
                "CC(=[Se])CCc1ccc(cc1)C(=O)O",
                "4-(3-selanylidenebutyl)benzoic acid",
            ),
            ("OC(=O)CC(=S)CC(=S)C", "3,5-bis(sulfanylidene)hexanoic acid"),
            ("O=C(O)c1ccc(C=S)cc1", "4-(methanethioyl)benzoic acid"),
            ("O=C(O)c1ccc(C=[Se])cc1", "4-(methaneselenoyl)benzoic acid"),
            ("O=C(O)c1ccc(C=[Te])cc1", "4-(methanetelluroyl)benzoic acid"),
            ("O=C(O)c1ccc(C(=S)C)cc1", "4-(ethanethioyl)benzoic acid"),
            (
                "CC(=S)C(C(C)=S)C(C)=S",
                "3-(ethanethioyl)pentane-2,4-dithione",
            ),
            (
                "O=C(O)c1ccc(C(=S)c2ccccc2)cc1",
                "4-(benzenecarbothioyl)benzoic acid",
            ),
            // Amines, sulfides and hydrodisulfides are named substitutively
            // (P-62.2.2.1, P-63.2.5, P-63.4.2.1).
            ("CNC", "N-methylmethanamine"),
            ("CN(C)C", "N,N-dimethylmethanamine"),
            ("CSC", "(methylsulfanyl)methane"),
            ("SSCCOC", "2-methoxyethane-1-dithioperoxol"),
            ("SSc1ccccn1", "pyridine-2-dithioperoxol"),
            // Alkoxy prefixes other than the retained contractions (P-63.2.2.2).
            ("CC(C)Oc1ccccc1", "[(propan-2-yl)oxy]benzene"),
            ("CCC(C)Oc1ccccc1", "[(butan-2-yl)oxy]benzene"),
            ("CC(C)COc1ccccc1", "(2-methylpropoxy)benzene"),
            ("CC(C)(C)Oc1ccccc1", "tert-butoxybenzene"),
            ("CC(C)Oc1ccccc1OC(C)C", "1,2-bis[(propan-2-yl)oxy]benzene"),
            ("CC(C)COc1ccccc1OCC(C)C", "1,2-bis(2-methylpropoxy)benzene"),
            ("CC(C)OC(=O)c1ccccc1", "propan-2-yl benzoate"),
            // Retained names for general nomenclature only (P-65.1.1.2.2,
            // P-66.5.1.2.3); HCN is formonitrile (P-66.5.1.2.1).
            ("CCC#N", "propanenitrile"),
            ("CC(C)C#N", "2-methylpropanenitrile"),
            ("O=CCCC#N", "4-oxobutanenitrile"),
            ("O=C(O)CC(=O)O", "propanedioic acid"),
            ("O=C(O)CCC(=O)O", "butanedioic acid"),
            ("O=C(O)CCCC(=O)O", "pentanedioic acid"),
            ("O=C(O)CCCCC(=O)O", "hexanedioic acid"),
            ("C#N", "formonitrile"),
            // Retained PINs (P-66.6.1.2.1, P-34.1.1.4, P-66.1.1.1.2.1, P-66.5.4.1).
            ("C=O", "formaldehyde"),
            ("CC=O", "acetaldehyde"),
            ("COc1ccccc1", "anisole"),
            ("NC(=O)C(N)=O", "oxamide"),
            ("C#[N+][O-]", "formonitrile oxide"),
            ("CC#[N+][O-]", "acetonitrile oxide"),
            ("ClCC#[N+][O-]", "chloroacetonitrile oxide"),
            // Anilino (P-62.2.1.1.1) and carbamimidoylamino (P-66.4.1.2.1.3)
            // prefixes; contractions still read anilino as phenylamino.
            ("O=C(O)CNc1ccccc1", "anilinoacetic acid"),
            (
                "OC(=O)c1ccc(N(c2ccccc2)c2ccccc2)cc1",
                "4-(N-phenylanilino)benzoic acid",
            ),
            ("NC(=N)NCC(=O)O", "(carbamimidoylamino)acetic acid"),
            (
                "CN(C(=N)N)CC(=O)O",
                "[carbamimidoyl(methyl)amino]acetic acid",
            ),
            ("COC(=O)Nc1ccccc1", "methyl phenylcarbamate"),
            (
                "OC(=O)c1ccc(C(=O)Nc2ccccc2)cc1",
                "4-(phenylcarbamoyl)benzoic acid",
            ),
            // Oxazole and thiazole fusion prefixes cite their Hantzsch-Widman
            // locants (P-25.3.2.2.3).
            ("O1C=NC2=NC=CC=C21", "[1,3]oxazolo[4,5-b]pyridine"),
            ("c1cnc2scnc2c1", "[1,3]thiazolo[5,4-b]pyridine"),
            ("c1cnc2oncc2c1", "[1,2]oxazolo[5,4-b]pyridine"),
            ("c1cnc2sncc2c1", "[1,2]thiazolo[5,4-b]pyridine"),
            // Differently substituted rings stay substitutive (P-51.3.1), and
            // the retained acyl and formyl prefixes need no ring acid name.
            (
                "Brc1cccc(Cc2cccc(Cl)c2)c1",
                "1-bromo-3-[(3-chlorophenyl)methyl]benzene",
            ),
            ("O=C(O)c1ccc(cc1)C(=O)c1ccccc1", "4-benzoylbenzoic acid"),
            ("O=C(O)c1ccc(cc1)C=O", "4-formylbenzoic acid"),
            // A chalcogen or peroxy prefix encloses a ligand that names a
            // parent hydride or is complex (P-16.5.1.4, P-16.5.1.1).
            (
                "CS(=O)(=O)SCC(=O)O",
                "[(methanesulfonyl)sulfanyl]acetic acid",
            ),
            ("CC(C)OOCC(=O)O", "[(propan-2-yl)peroxy]acetic acid"),
            ("CC(=S)SCC(=O)O", "[(ethanethioyl)sulfanyl]acetic acid"),
            ("CC(=S)OOCC(=O)O", "[(ethanethioyl)peroxy]acetic acid"),
            // Prefixes are numbered and cited from their first Roman letter;
            // bracketed initial locants only break ties (P-14.5, P-14.5.4).
            (
                "OC(=O)c1cc(C)cc(c2nc3ncccc3o2)c1",
                "3-methyl-5-([1,3]oxazolo[4,5-b]pyridin-2-yl)benzoic acid",
            ),
            (
                "O1N=CC2=NC=C(C=C21)CC(CO)C=2C=C1C(=NC2)N=CO1",
                "3-([1,2]oxazolo[4,5-b]pyridin-6-yl)-2-([1,3]oxazolo[4,5-b]pyridin-6-yl)propan-1-ol",
            ),
            // Units that differ in configuration keep the substitutive PIN
            // (P-45.6.2), as does an anhydride (P-65.7.1).
            (
                "CC[C@@H](C)c1ccc(Sc2ccc(cc2)[C@@H](C)CC)cc1",
                "1-[(2R)-butan-2-yl]-4-({4-[(2S)-butan-2-yl]phenyl}sulfanyl)benzene",
            ),
            ("CC(=O)OC(C)=O", "acetic anhydride"),
            // Biphenyl is the ring assembly name itself (P-28.2.1); different
            // rings are not an assembly (P-28.7, P-61.2.2).
            ("c1ccc(cc1)-c1ccccc1", "1,1'-biphenyl"),
            ("C1CCC(CC1)c1ccccc1", "cyclohexylbenzene"),
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            let status = match smiles {
                "SS" | "SSS" => NameStatus::Unverified(nazuchem_ir::NameReason::UnauditedPath),
                _ => NameStatus::Preferred,
            };
            assert_eq!(name.status, status, "{smiles}");
            if expected.ends_with("acetic acid") || expected == "chloroacetonitrile oxide" {
                let parsed = crate::test_support::parse_strict(&name.text).unwrap();
                assert_eq!(
                    crate::test_support::canonical_smiles(smiles).unwrap(),
                    crate::test_support::canonical_smiles(&parsed).unwrap(),
                    "{}",
                    name.text
                );
            }
        }
    }

    #[test]
    // Provenance: fn tests::ring_assemblies_are_preferred_numbered_parents_and_round_trip | adapted
    // Source: ReShiki | src/naming/openclatura/mod.rs::tests::ring_assemblies_are_preferred_numbered_parents_and_round_trip
    // Source: openclatura | tests/unit/test_analysis.py::test_biphenyl_is_graph_selected_not_exact_name_rewritten
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn ring_assemblies_are_preferred_numbered_parents_and_round_trip() {
        use nazuchem_ir::NameStatus;
        for (smiles, expected) in [
            (
                "Oc1ccc2ccccc2c1-c1c(O)ccc2ccccc12",
                "[1,1'-binaphthalene]-2,2'-diol",
            ),
            ("Oc1ccccc1-c1ccccc1O", "[1,1'-biphenyl]-2,2'-diol"),
            ("Oc1ccc(cc1)-c1ccccc1", "[1,1'-biphenyl]-4-ol"),
            (
                "OC(=O)c1ccc(cc1)-c1ccccc1",
                "[1,1'-biphenyl]-4-carboxylic acid",
            ),
            ("c1ccc(cc1)-c1ccccc1", "1,1'-biphenyl"),
            ("Cc1ccc(cc1)-c1ccc(C)cc1", "4,4'-dimethyl-1,1'-biphenyl"),
            ("Cc1ccc(cc1)-c1ccccc1", "4-methyl-1,1'-biphenyl"),
            ("c1ccc(nc1)-c1ccccn1", "2,2'-bipyridine"),
            ("c1ccncc1-c1cccnc1", "3,3'-bipyridine"),
            ("c1ccncc1-c1ccncc1", "3,4'-bipyridine"),
            ("c1ccc(nc1)-c1cccnc1", "2,3'-bipyridine"),
            ("c1cnccc1-c1ccncc1", "4,4'-bipyridine"),
            ("c1ccc(s1)-c1cccs1", "2,2'-bithiophene"),
            ("C1CC1C1CC1", "1,1'-bi(cyclopropane)"),
            ("CC1(CC1)C1CC1", "1-methyl-1,1'-bi(cyclopropane)"),
            ("C1CCC(CC1)C1CCCCC1", "1,1'-bi(cyclohexane)"),
            (
                "O=C1CCC(C2CCC(=O)CC2)CC1",
                "[1,1'-bi(cyclohexane)]-4,4'-dione",
            ),
            ("c1ccc2ccccc2c1-c1ccc2ccccc2c1", "1,2'-binaphthalene"),
            ("c1ccc2ccccc2c1-c1cccc2ccccc12", "1,1'-binaphthalene"),
            (
                "c1ccc(cc1)-c1ccc(cc1)-c1ccccc1",
                "1^1,2^1:2^4,3^1-terphenyl",
            ),
            (
                "c1ccc(cc1)-c1cccc(-c2ccccc2)c1",
                "1^1,2^1:2^3,3^1-terphenyl",
            ),
            (
                "Oc1ccc(cc1)-c1ccc(cc1)-c1ccccc1",
                "[1^1,2^1:2^4,3^1-terphenyl]-1^4-ol",
            ),
            (
                "c1(-c2ccccc2)cc(-c2ccccc2)cc(-c2ccccc2)c1",
                "2^5-phenyl-1^1,2^1:2^3,3^1-terphenyl",
            ),
            // P-28.6: seven units, but only three in a row.
            (
                "c1(-c2ccccc2)c(-c2ccccc2)c(-c2ccccc2)c(-c2ccccc2)c(-c2ccccc2)c1-c1ccccc1",
                "2^3,2^4,2^5,2^6-tetraphenyl-1^1,2^1:2^2,3^1-terphenyl",
            ),
            // P-31.1.4.3.4: r before s, whichever traversal the input uses.
            (
                "C[C@H]1CC[C@@H]([C@H]2CC[C@H](C)CC2)CC1",
                "(1r,1's,4r,4's)-4,4'-dimethyl-1,1'-bi(cyclohexane)",
            ),
            (
                "C[C@H]1CC[C@H]([C@H]2CC[C@@H](C)CC2)CC1",
                "(1r,1's,4r,4's)-4,4'-dimethyl-1,1'-bi(cyclohexane)",
            ),
            ("c1c[nH]c(c1)-c1ccc[nH]1", "1H,1'H-2,2'-bipyrrole"),
            ("c1ccn(c1)-n1cccc1", "1,1'-bipyrrole"),
            (
                "OC(=O)Cc1ccc(cc1)-c1ccccc1",
                "([1,1'-biphenyl]-4-yl)acetic acid",
            ),
            (
                "C[C@H]1C[C@@H]1[C@H]1C[C@@H]1C",
                "(1S,1'S,2S,2'S)-2,2'-dimethyl-1,1'-bi(cyclopropane)",
            ),
            ("C1CCC(CC1)c1ccccc1", "cyclohexylbenzene"),
            ("c1ccc(cc1)-c1ccncc1", "4-phenylpyridine"),
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, NameStatus::Preferred, "{smiles}");
            let parsed = crate::test_support::parse_strict(&name.text).unwrap();
            assert_eq!(
                crate::test_support::canonical_smiles(smiles).unwrap(),
                crate::test_support::canonical_smiles(&parsed).unwrap(),
                "{}",
                name.text
            );
        }
    }

    #[test]
    fn linear_phane_assemblies_keep_nonpreferred_status_and_round_trip() {
        use nazuchem_ir::{NameReason, NameStatus};
        for (smiles, expected, status) in [
            // P-52.2.5.3: unsubstituted linear phanes take the phane name.
            (
                "c1ccccc1-c1ccc(-c2ccc(Oc3cccc(-c4cccc(-c5ccccn5)n4)n3)cc2)cc1",
                "4-oxa-1(2),2,3(2,6)-tripyridina-5,6(1,4),7(1)-tribenzenaheptaphane",
                NameStatus::Preferred,
            ),
            // P-51.3.2: even identical halves require the linear phane,
            // rather than a multiplicative name, for seven skeletal nodes.
            (
                "c1ccccc1-c1ccc(-c2ccc(Oc3ccc(-c4ccc(-c5ccccc5)cc4)cc3)cc2)cc1",
                "4-oxa-1,7(1),2,3,5,6(1,4)-hexabenzenaheptaphane",
                NameStatus::Preferred,
            ),
            // P-65.6.3.3.7.2.1: the carbamate separates two two-ring
            // components; its named biphenyl junction needs no demotion.
            (
                "CN(CCCCCCCCCN1CCC(CC1)OC(=O)NC2=CC=CC=C2C3=CC=CC=C3)C(=O)C4=CC(=C(C=C4)O)Cl",
                "1-[9-(3-chloro-4-hydroxy-N-methylbenzamido)nonyl]piperidin-4-yl ([1,1'-biphenyl]-2-yl)carbamate",
                NameStatus::Preferred,
            ),
            // P-52.2.5.3 applies without identical-unit ring junctions.
            (
                "Oc1ccc(Oc2ccc(Oc3ccc(Oc4ccccc4)cc3)cc2)cc1",
                "4-[4-(4-phenoxyphenoxy)phenoxy]phenol",
                NameStatus::NotPreferred(NameReason::LinearPhane),
            ),
            (
                "Oc1ccc(Cc2ccc(Cc3ccc(Cc4ccccc4)cc3)cc2)cc1",
                "4-({4-[(4-benzylphenyl)methyl]phenyl}methyl)phenol",
                NameStatus::NotPreferred(NameReason::LinearPhane),
            ),
            (
                "Oc1ccc(Oc2ccc(Sc3ccc(Oc4ccccc4)cc3)cc2)cc1",
                "4-{4-[(4-phenoxyphenyl)sulfanyl]phenoxy}phenol",
                NameStatus::NotPreferred(NameReason::LinearPhane),
            ),
            (
                "Oc1ccc(-c2ccnc(-c3ccc(-c4ccnc(-c5ccc(-c6ccnc(-c7ccccc7)c6)cc5)c4)cc3)c2)cc1",
                "4-[2-(4-{2-[4-(2-phenylpyridin-4-yl)phenyl]pyridin-4-yl}phenyl)pyridin-4-yl]phenol",
                NameStatus::NotPreferred(NameReason::LinearPhane),
            ),
            // A linear phane outranks the multiplicative name (P-51.3.2).
            (
                "c1ccccc1Oc2ccc(Oc3ccc(Oc4ccccc4)cc3)cc2",
                "2,4,6-trioxa-1,7(1),3,5(1,4)-tetrabenzenaheptaphane",
                NameStatus::Preferred,
            ),
            (
                "c1ccccc1-c1ccc(-c2ccc(-c3ccc(-c4ccc(-c5ccccc5)cc4)cc3)cc2)cc1",
                "1^1,2^1:2^4,3^1:3^4,4^1:4^4,5^1:5^4,6^1-sexiphenyl",
                NameStatus::Preferred,
            ),
            (
                "c1(-c2ccccc2)c(-c2ccccc2)c(-c2ccccc2)c(-c2ccccc2)c(-c2ccccc2)c1-c1ccccc1",
                "2^3,2^4,2^5,2^6-tetraphenyl-1^1,2^1:2^2,3^1-terphenyl",
                NameStatus::Preferred,
            ),
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, status, "{smiles}");
            let parsed = crate::test_support::parse_strict(&name.text).unwrap();
            assert_eq!(
                crate::test_support::canonical_smiles(smiles).unwrap(),
                crate::test_support::canonical_smiles(&parsed).unwrap(),
                "{}",
                name.text
            );
        }
    }

    #[test]
    fn linear_phane_status_respects_functional_class_components_and_round_trips() {
        use nazuchem_ir::{NameReason, NameStatus};
        for (smiles, expected, status) in [
            // P-65.6.3.2.1 / P-65.6.3.3.1: each ester component has two rings.
            (
                "O=C(OCc1ccc(Oc2ccccc2)cc1)c1ccc(Oc2ccccc2)cc1",
                "(4-phenoxyphenyl)methyl 4-phenoxybenzoate",
                NameStatus::Preferred,
            ),
            (
                "O=C(OCc1ccc(-c2ccccc2)cc1)c1ccc(-c2ccccc2)cc1",
                "([1,1'-biphenyl]-4-yl)methyl [1,1'-biphenyl]-4-carboxylate",
                NameStatus::Preferred,
            ),
            // P-65.6.3.3.7.2.1: carbamic acid esters use the same split.
            (
                "O=C(Oc1ccc(Oc2ccccc2)cc1)Nc1ccc(Oc2ccccc2)cc1",
                "4-phenoxyphenyl (4-phenoxyphenyl)carbamate",
                NameStatus::Preferred,
            ),
            // P-65.7.1: anhydrides cite separate acid components.
            (
                "O=C(OC(=O)c1ccc(Oc2ccccc2)cc1)c1ccc(Oc2ccccc2)cc1",
                "4-phenoxybenzoic anhydride",
                NameStatus::Preferred,
            ),
            // P-65.6.3.3.4.4: a phane within either component still counts.
            (
                "COC(=O)c1ccc(Oc2ccc(Oc3ccc(Oc4ccccc4)cc3)cc2)cc1",
                "methyl 4-[4-(4-phenoxyphenoxy)phenoxy]benzoate",
                NameStatus::NotPreferred(NameReason::LinearPhane),
            ),
            (
                "CC(=O)Oc1ccc(Oc2ccc(Oc3ccc(Oc4ccccc4)cc3)cc2)cc1",
                "4-[4-(4-phenoxyphenoxy)phenoxy]phenyl acetate",
                NameStatus::NotPreferred(NameReason::LinearPhane),
            ),
            // P-65.6.3.3.4: an ester prefix does not split the acid parent.
            (
                "OC(=O)c1ccc(Oc2ccc(C(=O)Oc3ccc(Oc4ccccc4)cc3)cc2)cc1",
                "4-{4-[(4-phenoxyphenoxy)carbonyl]phenoxy}benzoic acid",
                NameStatus::NotPreferred(NameReason::LinearPhane),
            ),
            // A same-key ester outside the selected parent remains a prefix.
            (
                "COC(=O)c1cc(C(=O)OC)cc(Oc2ccc(C(=O)Oc3ccc(Oc4ccccc4)cc3)cc2)c1",
                "dimethyl 5-{4-[(4-phenoxyphenoxy)carbonyl]phenoxy}benzene-1,3-dicarboxylate",
                NameStatus::NotPreferred(NameReason::LinearPhane),
            ),
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, status, "{smiles}");
            let parsed = crate::test_support::parse_strict(&name.text).unwrap();
            assert_eq!(
                crate::test_support::canonical_smiles(smiles).unwrap(),
                crate::test_support::canonical_smiles(&parsed).unwrap(),
                "{}",
                name.text
            );
        }
    }

    #[test]
    fn linear_phane_status_respects_oxoacid_ester_components_and_round_trips() {
        use nazuchem_ir::{NameReason, NameStatus};
        for (smiles, expected, status) in [
            // P-65.2.3 / P-65.6.3.2.1: ester organyl groups are separate components.
            (
                "O=C(Oc1ccc(Oc2ccccc2)cc1)Oc1ccc(Oc2ccccc2)cc1",
                "bis(4-phenoxyphenyl) carbonate",
                NameStatus::Preferred,
            ),
            (
                "O=P(Oc1ccc(Oc2ccccc2)cc1)(Oc1ccccc1)Oc1ccc(Oc2ccccc2)cc1",
                "bis(4-phenoxyphenyl) phenyl phosphate",
                NameStatus::Preferred,
            ),
            (
                "O=S(=O)(Oc1ccc(Oc2ccccc2)cc1)Oc1ccc(Oc2ccccc2)cc1",
                "bis(4-phenoxyphenyl) sulfate",
                NameStatus::Preferred,
            ),
            (
                "O=S(=O)(Oc1ccc(Oc2ccccc2)cc1)c1ccc(Oc2ccccc2)cc1",
                "4-phenoxyphenyl 4-phenoxybenzene-1-sulfonate",
                NameStatus::Preferred,
            ),
            (
                "CP(=O)(Oc1ccc(Oc2ccccc2)cc1)Oc1ccc(Oc2ccccc2)cc1",
                "bis(4-phenoxyphenyl) methylphosphonate",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "COP(=O)(Oc1ccc(Oc2ccccc2)cc1)c1ccc(Oc2ccccc2)cc1",
                "methyl 4-phenoxyphenyl 4-phenoxyphenylphosphonate",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "OP(=O)(Oc1ccc(Oc2ccccc2)cc1)c1ccc(Oc2ccccc2)cc1",
                "4-phenoxyphenyl hydrogen 4-phenoxyphenylphosphonate",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "[PH](=O)(Oc1ccc(Oc2ccccc2)cc1)c1ccc(Oc2ccccc2)cc1",
                "4-phenoxyphenyl (4-phenoxyphenyl)phosphinate",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            // P-52.2.5.3: a four-ring path within one component still counts.
            (
                "COC(=O)Oc1ccc(Oc2ccc(Oc3ccc(Oc4ccccc4)cc3)cc2)cc1",
                "methyl 4-[4-(4-phenoxyphenoxy)phenoxy]phenyl carbonate",
                NameStatus::NotPreferred(NameReason::LinearPhane),
            ),
            (
                "COP(=O)(c1ccc(Oc2ccccc2)cc1)c1ccc(Oc2ccccc2)cc1",
                "methyl bis(4-phenoxyphenyl)phosphinate",
                NameStatus::NotPreferred(NameReason::LinearPhane),
            ),
            (
                "OP(=O)(c1ccc(Oc2ccccc2)cc1)c1ccc(Oc2ccccc2)cc1",
                "bis(4-phenoxyphenyl)phosphinic acid",
                NameStatus::NotPreferred(NameReason::LinearPhane),
            ),
            // A phosphonate prefix under a senior nitrile keeps the full skeleton.
            (
                "N#Cc1ccc(Oc2ccc(P(=O)(OC)Oc3ccc(Oc4ccccc4)cc3)cc2)cc1",
                "4-{4-[methoxy(4-phenoxyphenoxy)phosphoryl]phenoxy}benzonitrile",
                NameStatus::NotPreferred(NameReason::LinearPhane),
            ),
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, status, "{smiles}");
            let parsed = crate::test_support::parse_strict(&name.text).unwrap();
            assert_eq!(
                crate::test_support::canonical_smiles(smiles).unwrap(),
                crate::test_support::canonical_smiles(&parsed).unwrap(),
                "{}",
                name.text
            );
        }
    }

    #[test]
    fn substituent_amine_sulfide_is_cited_as_a_prefix() {
        use nazuchem_ir::{NameReason, NameStatus};
        let name = name_smiles(
            "O=P(Oc1ccc(Oc2ccccc2)cc1)(Oc1ccccc1)Oc1ccc(Oc2ccc([N+](C)(C)[S-])cc2)cc1",
            &Default::default(),
        )
        .unwrap();
        // P-62.5 names only a parent amine; this N-sulfide is a prefix (#947).
        assert_eq!(
            name.text,
            "4-{4-[dimethyl(sulfido)ammonio]phenoxy}phenyl 4-phenoxyphenyl phenyl phosphate"
        );
        assert_eq!(
            name.status,
            NameStatus::NotPreferred(NameReason::AmineOxideZwitterion)
        );
    }

    #[test]
    fn tied_ring_assemblies_choose_the_senior_component_and_round_trip() {
        let smiles = "C1(C3CC4CCC34)C(C5CC5C6CC6C7CC7C8CC8)C2CCC12";
        let name = name_smiles(smiles, &Default::default()).unwrap();
        assert_eq!(
            name.text,
            "3-([1^1,2^1:2^2,3^1:3^2,4^1-quatercyclopropan]-1^2-yl)-2,2'-bi(bicyclo[2.2.0]hexane)"
        );
        assert_eq!(name.status, nazuchem_ir::NameStatus::Preferred);
        let parsed = crate::test_support::parse_strict(&name.text).unwrap();
        assert_eq!(
            crate::test_support::canonical_smiles(smiles).unwrap(),
            crate::test_support::canonical_smiles(&parsed).unwrap()
        );
    }

    #[test]
    fn unsupported_ring_assembly_keeps_substitutive_name_and_status() {
        let name = name_smiles("C1CCC(=C2CCCC2)C1", &Default::default()).unwrap();
        assert_eq!(name.text, "cyclopentylidenecyclopentane");
        assert_eq!(
            name.status,
            nazuchem_ir::NameStatus::NotPreferred(nazuchem_ir::NameReason::RingAssembly)
        );
    }

    #[test]
    fn composite_quaterpyridine_and_assembly_substituent_locants_round_trip() {
        for expected in [
            "1^2,2^2:2^6,3^2:3^6,4^2-quaterpyridine",
            "1^1,2^1:2^2,3^1-tercyclopropane",
            "1^1,2^1:2^4,3^1:3^4,4^1:4^4,5^1-quinquephenyl",
            "1^1,2^1:2^4,3^1:3^4,4^1:4^4,5^1:5^4,6^1-sexiphenyl",
            // Symmetric units collapse to the numberings a name can cite.
            "1^1,2^1:2^4,3^1:3^4,4^1:4^4,5^1-quinquebicyclo[2.2.2]octane",
            "6H,6'H-2,2'-bipyran",
            "1H,1'H-1,1'-biindene",
            "1H,1'H-3,3'-biindole",
            "1H,1'H-[2,2'-bipyrrol]-5-ol",
            "(1H,1'H-[2,2'-bipyrrol]-5-yl)acetic acid",
            "2,3'-bifuran",
            "2,2'-bi(bicyclo[2.2.2]octane)",
            "1,1'-bi(spiro[2.2]pentane)",
            "1,1'-bi(dispiro[2.0.2.1]heptane)",
            "([2,2'-bipyridin]-5-yl)acetic acid",
            "([1^1,2^1:2^4,3^1-terphenyl]-1^4-yl)acetic acid",
            "(4'-cyano[1,1'-biphenyl]-4-yl)acetic acid",
        ] {
            let parsed = crate::test_support::parse_strict(expected).unwrap();
            let name = name_smiles(&parsed, &Default::default()).unwrap();
            assert_eq!(name.text, expected);
            assert_eq!(name.status, nazuchem_ir::NameStatus::Preferred);
        }
    }

    #[test]
    fn symmetric_polyspiro_assembly_junctions_use_every_equivalent_numbering() {
        // P-24.2.2 fixes each dispiro parent numbering; P-28.2.1 then gives
        // both equivalent terminal rings the lowest possible junction locant.
        for smiles in [
            "C1CC12CC21CC1C1CC12CC21CC1",
            "C1C2(C1)CC12CC1C1C2(C3(CC3)C2)C1",
            "C12(CC1)C1(CC1C1C3(CC43CC4)C1)C2",
            "C1(C2C3(CC43CC4)C2)C2(C1)CC12CC1",
            "C1CC21CC21CC1C1CC21C1(CC1)C2",
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, "1,1'-bi(dispiro[2.0.2.1]heptane)", "{smiles}");
            assert_eq!(name.status, nazuchem_ir::NameStatus::Preferred, "{smiles}");
            let parsed = crate::test_support::parse_strict(&name.text).unwrap();
            assert_eq!(
                crate::test_support::canonical_smiles(smiles).unwrap(),
                crate::test_support::canonical_smiles(&parsed).unwrap(),
                "{}",
                name.text
            );
        }
    }

    #[test]
    // Provenance: fn tests::locanted_assembly_hydrides_are_enclosed_for_two_through_six_units | adapted
    // Source: ReShiki | src/naming/openclatura/mod.rs::tests::locanted_assembly_hydrides_are_enclosed_for_two_through_six_units
    // Source: openclatura | tests/data/parser_xml_resources/multipliers.json::<fixture data>
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn locanted_assembly_hydrides_are_enclosed_for_two_through_six_units() {
        // P-28.2.1 and P-16.5 separate the component's locants from those
        // identifying the junctions, independently of the unit multiplier.
        let smiles = "C1(OCCOC1)C1COCCO1";
        let name = name_smiles(smiles, &Default::default()).unwrap();
        assert_eq!(name.text, "2,2'-bi(1,4-dioxane)");
        assert_eq!(name.status, nazuchem_ir::NameStatus::Preferred);
        let parsed = crate::test_support::parse_strict(&name.text).unwrap();
        assert_eq!(
            crate::test_support::canonical_smiles(smiles).unwrap(),
            crate::test_support::canonical_smiles(&parsed).unwrap()
        );
        for expected in [
            "2,2'-bi(1,4-dioxane)",
            "1^2,2^2:2^3,3^2-ter(1,4-dioxane)",
            "1^2,2^2:2^3,3^2:3^3,4^2-quater(1,4-dioxane)",
            "1^2,2^2:2^3,3^2:3^3,4^2:4^3,5^2-quinque(1,4-dioxane)",
            "1^2,2^2:2^3,3^2:3^3,4^2:4^3,5^2:5^3,6^2-sexi(1,4-dioxane)",
        ] {
            let parsed = crate::test_support::parse_strict(expected).unwrap();
            let name = name_smiles(&parsed, &Default::default()).unwrap();
            assert_eq!(name.text, expected);
            assert_eq!(name.status, nazuchem_ir::NameStatus::Preferred);
            let reparsed = crate::test_support::parse_strict(&name.text).unwrap();
            assert_eq!(
                crate::test_support::canonical_smiles(&parsed).unwrap(),
                crate::test_support::canonical_smiles(&reparsed).unwrap()
            );
        }
        for multiplier in ["bi", "ter", "quater", "quinque", "sexi"] {
            let enclosed = format!("2-(1,1'-{multiplier}(1,4-dioxan)-2-yl)acetic acid");
            assert_eq!(nest_enclosing_marks(&enclosed), enclosed);
        }
    }

    #[test]
    fn fused_heterocycle_assemblies_keep_parent_locants_and_indicated_hydrogen() {
        for (smiles, expected) in [
            ("c1ccc2c(c1)ccc(n2)-c1ccc2ccccc2n1", "2,2'-biquinoline"),
            ("c1ccc2c(c1)oc(c2)-c1cc2ccccc2o1", "2,2'-bi(1-benzofuran)"),
            // #28, P-28.2.1: parentheses keep "bi" apart from the "di" of
            // the multiplied fusion prefix.
            (
                "c1ccc2c(c1)oc1cc(ccc12)-c1cccc2c1oc1ccccc12",
                "3,4'-bi(dibenzofuran)",
            ),
            (
                "c1ccc2c(c1)sc1cc(ccc12)-c1cccc2c1sc1ccccc12",
                "3,4'-bi(dibenzothiophene)",
            ),
            (
                "c1ccc2c(-c3c[nH]c4ccccc34)c[nH]c2c1",
                "1H,1'H-3,3'-biindole",
            ),
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, nazuchem_ir::NameStatus::Preferred, "{smiles}");
            let parsed = crate::test_support::parse_strict(&name.text).unwrap();
            assert_eq!(
                crate::test_support::canonical_smiles(smiles).unwrap(),
                crate::test_support::canonical_smiles(&parsed).unwrap(),
                "{}",
                name.text
            );
        }
        // P-16.2.4.1 also permits separating the intrinsic fusion-parent
        // locant with a hyphen. Both spellings reconstruct the same hydride;
        // generation consistently encloses locanted component names.
        let alternative = "2,2'-bi-1-benzofuran";
        let parsed = crate::test_support::parse_strict(alternative).unwrap();
        assert_eq!(
            crate::test_support::canonical_smiles("c1ccc2c(c1)oc(c2)-c1cc2ccccc2o1").unwrap(),
            crate::test_support::canonical_smiles(&parsed).unwrap()
        );
    }

    #[test]
    // Provenance: fn tests::preferred_von_baeyer_controls_do_not_report_fusion_fallback | adapted
    // Source: ReShiki | src/naming/openclatura/mod.rs::tests::preferred_von_baeyer_controls_do_not_report_fusion_fallback
    // Source: openclatura | tests/unit/test_fused_bonds.py::test_fused_bonds_match_the_cycle_count
    // Source: openclatura | tests/unit/test_von_baeyer_candidate_audit.py::test_the_preferred_candidate_is_still_offered_first
    // Source: openclatura | tests/roundtrip/test_polycycles.py::test_roundtrip
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn preferred_von_baeyer_controls_do_not_report_fusion_fallback() {
        use nazuchem_ir::NameStatus;
        for (smiles, expected) in [
            ("C1CCC2CCCCC2C1", "decahydronaphthalene"),
            ("C1CC2=CC=CC=C12", "bicyclo[4.2.0]octa-1,3,5-triene"),
            ("C1CC2CCC1CC2", "bicyclo[2.2.2]octane"),
            ("C1C2CC12", "bicyclo[1.1.0]butane"),
            ("C1CC2CCC1C2", "bicyclo[2.2.1]heptane"),
            ("C1CCC2(CC1)CCCC2", "spiro[4.5]decane"),
            // A retained fused parent bypasses the fallback assessment.
            ("c1ccc2ccccc2c1", "naphthalene"),
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(
                name.status,
                NameStatus::Preferred,
                "{smiles}: {}",
                name.text
            );
        }
    }

    #[test]
    // Provenance: fn tests::pheno_fusion_components_choose_preferred_parents_and_round_trip | original
    // Source: ReShiki | src/naming/openclatura/mod.rs::tests::pheno_fusion_components_choose_preferred_parents_and_round_trip
    fn pheno_fusion_components_choose_preferred_parents_and_round_trip() {
        use nazuchem_ir::NameStatus;
        for (smiles, expected) in [
            (
                "c1ccc2c(c1)Nc1cc3cnc4ccccc4c3cc1O2",
                "8H-quinolino[3,4-b]phenoxazine",
            ),
            (
                "C1=CC=C2C(=C1)NC1=CC3=CN=C4C=CC=CC4=C3C=C1O2",
                "8H-quinolino[3,4-b]phenoxazine",
            ),
            (
                "c1ccc2c(c1)Nc1cc3cnc4ccccc4c3cc1S2",
                "8H-quinolino[3,4-b]phenothiazine",
            ),
            (
                "c1ccc2cc3c(cc2c1)Nc1cc2ccccc2cc1O3",
                "13H-dibenzo[b,i]phenoxazine",
            ),
            ("c1ccc2c(c1)Sc1ccc3ccccc3c1S2", "benzo[a]thianthrene"),
            ("c1ccc2c(c1)Oc1ccc3ccccc3c1S2", "benzo[a]phenoxathiine"),
            ("c1ccc2c(c1)Nc1ccc3ccccc3c1S2", "7H-benzo[c]phenothiazine"),
            (
                "c1ccc2c(c1)Nc1c(ccc3ccccc13)S2",
                "12H-benzo[a]phenothiazine",
            ),
            ("c1ccc2c(c1)Nc1ccc3ccccc3c1O2", "7H-benzo[c]phenoxazine"),
            ("c1ccc2c(c1)Oc1ccc3ccccc3c1N2", "12H-benzo[a]phenoxazine"),
            ("c1ccc2c(c1)Nc1cc3ccccc3cc1S2", "12H-benzo[b]phenothiazine"),
            ("c1ccc2cc3c(cc2c1)Sc1ccccc1S3", "benzo[b]thianthrene"),
            ("c1cc2Nc3ccccc3Oc2c2ncccc12", "7H-pyrido[3,2-c]phenoxazine"),
            (
                "c1ccc2c(c1)Oc1ccc3ncccc3c1N2",
                "12H-pyrido[3,2-a]phenoxazine",
            ),
            (
                "c1ccc2nc3cc4Oc5ccccc5Nc4cc3nc2c1",
                "14H-quinoxalino[2,3-b]phenoxazine",
            ),
            (
                "c1ccc2c(c1)Oc1cc3Sc4ccccc4Sc3cc1S2",
                "[1,4]benzodithiino[2,3-b]phenoxathiine",
            ),
            (
                "c1ccc2c(c1)Nc1cc3Sc4ccccc4Nc3cc1O2",
                "7,14-dihydro[1,4]benzothiazino[2,3-b]phenoxazine",
            ),
            // P-25.2.2.2: heteranthrenes (same heteroatoms) are PINs.
            ("c1ccc2c(c1)Oc1ccccc1O2", "oxanthrene"),
            ("c1ccc2c(c1)[Se]c1ccccc1[Se]2", "selenanthrene"),
            ("c1ccc2c(c1)[Te]c1ccccc1[Te]2", "telluranthrene"),
            // P-25.2.2.3: pheno...ine names (different heteroatoms) are PINs.
            ("c1ccc2c(c1)Oc1ccccc1[Se]2", "phenoxaselenine"),
            ("C1=CC=C2C(=C1)OC1=CC=CC=C1[Se]2", "phenoxaselenine"),
            ("c1ccc2c(c1)Oc1ccccc1[Te]2", "phenoxatellurine"),
            ("c1ccc2c(c1)Nc1ccccc1[Se]2", "10H-phenoselenazine"),
            ("[Se]1c2ccccc2Nc2ccccc12", "10H-phenoselenazine"),
            ("c1ccc2c(c1)Nc1ccccc1[Te]2", "10H-phenotellurazine"),
            ("c1ccc2c(c1)Pc1ccccc1O2", "10H-phenoxaphosphinine"),
            ("c1ccc2c(c1)[AsH]c1ccccc1O2", "10H-phenoxarsinine"),
            ("c1ccc2c(c1)[SbH]c1ccccc1O2", "10H-phenoxastibinine"),
            ("c1ccc2c(c1)[AsH]c1ccccc1S2", "10H-phenothiarsinine"),
            // P-25.3.2.4: they are senior parent components in fusion names;
            // (a) senior heteroatom, (b) more rings.
            ("c1ccc2c(c1)Oc1ccc3ccccc3c1O2", "benzo[a]oxanthrene"),
            ("c1ccc2cc3c(cc2c1)Oc1ccccc1O3", "benzo[b]oxanthrene"),
            (
                "c1ccc2c(c1)Pc1ccc3ccccc3c1O2",
                "7H-benzo[c]phenoxaphosphinine",
            ),
            (
                "c1ccc2c(c1)Oc1ccc3[Se]c4ccccc4[Se]c3c1O2",
                "[1,4]benzodiselenino[2,3-a]oxanthrene",
            ),
            ("Cc1ccc2Oc3ccccc3Oc2c1", "2-methyloxanthrene"),
            ("OC(=O)c1ccc2Oc3ccccc3Oc2c1", "oxanthrene-2-carboxylic acid"),
            ("CN1c2ccccc2[Se]c2ccccc12", "10-methyl-10H-phenoselenazine"),
            (
                "C1CCC2C(C1)Oc1ccccc1O2",
                "1,2,3,4,4a,10a-hexahydrooxanthrene",
            ),
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, NameStatus::Preferred, "{smiles}");
            let parsed = crate::test_support::parse_strict(&name.text).unwrap();
            crate::test_support::verify_identity(smiles, &parsed).unwrap();
        }
    }

    #[test]
    fn pheno_fusion_derivatives_use_the_senior_parent() {
        use nazuchem_ir::{NameReason, NameStatus};
        // #24: N-substituted, anionic, quaternary and ketone derivatives keep
        // the pheno...ine parent (P-25.2.2.3, P-25.3.2.4) and its hydrogen.
        for (smiles, expected) in [
            // P-58.2.4: prefixes replace the parent's N-H (indicated or hydro).
            (
                "c1ccc2c(c1)N(C)c1cc3cnc4ccccc4c3cc1O2",
                "8-methyl-8H-quinolino[3,4-b]phenoxazine",
            ),
            (
                "c1ccc2c(c1)N(c1ccccc1)c1cc3cnc4ccccc4c3cc1O2",
                "8-phenyl-8H-quinolino[3,4-b]phenoxazine",
            ),
            (
                "c1ccc2c(c1)Oc1ccc3ccccc3c1N2C",
                "12-methyl-12H-benzo[a]phenoxazine",
            ),
            (
                "CCN1c2ccccc2Sc2cc3ccccc3cc21",
                "12-ethyl-12H-benzo[b]phenothiazine",
            ),
            (
                "c1ccc2c(c1)N(C)c1cc3Sc4ccccc4Nc3cc1O2",
                "14-methyl-7,14-dihydro[1,4]benzothiazino[2,3-b]phenoxazine",
            ),
            (
                "c1ccc2c(c1)Nc1cc3Sc4ccccc4N(C)c3cc1O2",
                "7-methyl-7,14-dihydro[1,4]benzothiazino[2,3-b]phenoxazine",
            ),
            (
                "c1ccc2c(c1)N(C)c1cc3Sc4ccccc4N(C)c3cc1O2",
                "7,14-dimethyl-7,14-dihydro[1,4]benzothiazino[2,3-b]phenoxazine",
            ),
            // P-72.2.2.1 / P-73.1.1.2: '-ide' and '-ium' modify the parent
            // hydride, which keeps its indicated hydrogen (P-58.2.1.1).
            (
                "c1ccc2c(c1)Oc1ccc3ccccc3c1[N-]2",
                "12H-benzo[a]phenoxazin-12-ide",
            ),
            (
                "c1ccc2c(c1)Oc1ccc3ccccc3c1[N+]2(C)C",
                "12,12-dimethyl-12H-benzo[a]phenoxazin-12-ium",
            ),
            (
                "c1ccc2c(c1)Sc1ccc3ccccc3c1[N-]2",
                "12H-benzo[a]phenothiazin-12-ide",
            ),
            (
                "c1ccc2c(c1)Sc1ccc3ccccc3c1[N+]2(C)C",
                "12,12-dimethyl-12H-benzo[a]phenothiazin-12-ium",
            ),
            (
                "c1ccc2c(c1)[N-]c1cc3cnc4ccccc4c3cc1O2",
                "8H-quinolino[3,4-b]phenoxazin-8-ide",
            ),
            (
                "c1ccc2c(c1)[N+](C)(C)c1cc3cnc4ccccc4c3cc1O2",
                "8,8-dimethyl-8H-quinolino[3,4-b]phenoxazin-8-ium",
            ),
            // P-58.2.3.1.1: the ketone takes the indicated hydrogen
            // ('7H-1-benzopyran-7-one').
            (
                "CCN(CC)c1ccc2nc3c(cc(=O)c4ccccc34)oc2c1",
                "9-(diethylamino)-5H-benzo[a]phenoxazin-5-one",
            ),
            (
                "O=C1C=CC2=Nc3ccc4ccccc4c3OC2=C1",
                "10H-benzo[c]phenoxazin-10-one",
            ),
            // P-58.2.2: added hydrogen; the N-H stays as indicated or hydro
            // hydrogen of the parent hydride (P-58.2.4).
            (
                "O=C1CC=C2Oc3cc4Nc5ccccc5Sc4cc3NC2=C1",
                "7,14-dihydro[1,4]benzothiazino[2,3-b]phenoxazin-2(3H)-one",
            ),
            (
                "O=C1CC=C2Oc3cc4N(C)c5ccccc5Sc4cc3NC2=C1",
                "7-methyl-7,14-dihydro[1,4]benzothiazino[2,3-b]phenoxazin-2(3H)-one",
            ),
            (
                "O=C1CC=C2Oc3cc4Nc5ccccc5Sc4cc3N(C)C2=C1",
                "14-methyl-7,14-dihydro[1,4]benzothiazino[2,3-b]phenoxazin-2(3H)-one",
            ),
            // P-58.2.3.1.1: one indicated hydrogen, one suffix, as in
            // '1,2-dihydro-3H-indol-3-one' (P-64.3.1).
            (
                "O=C1CC=C2Nc3cc4ccccc4cc3OC2=C1",
                "2,12-dihydro-3H-benzo[b]phenoxazin-3-one",
            ),
            (
                "O=C1CC=C2N(C)c3cc4ccccc4cc3OC2=C1",
                "12-methyl-2,12-dihydro-3H-benzo[b]phenoxazin-3-one",
            ),
            (
                "O=C1C=Cc2ccc3Sc4ccccc4Sc3c2C1",
                "benzo[a]thianthren-2(1H)-one",
            ),
            // P-58.2.3.1.3: fewer indicated hydrogens than suffixes; the
            // indicated hydrogen takes the lowest mancude position, as in
            // '1H-cyclopenta[a]naphthalene-1,2(3H)-dione'.
            (
                "O=C1C=CC(=O)c2c1ccc1Oc3ccccc3Nc21",
                "4H-benzo[a]phenoxazine-1,4(12H)-dione",
            ),
            // #86: the unfused pheno...ine is the retained parent itself
            // (P-25.2.2.3), its indicated hydrogen goes to the ketone
            // (P-58.2.3.1.1, '1,2-dihydro-3H-indol-3-one' in P-64.3.1) or,
            // for a dione, to the lowest position (P-58.2.3.1.3).
            (
                "O=C1CC=C2Nc3ccccc3OC2=C1",
                "2,10-dihydro-3H-phenoxazin-3-one",
            ),
            (
                "O=C1CC=C2N(C)c3ccccc3OC2=C1",
                "10-methyl-2,10-dihydro-3H-phenoxazin-3-one",
            ),
            (
                "O=C1CC=C2Nc3ccc(Cl)cc3OC2=C1",
                "7-chloro-2,10-dihydro-3H-phenoxazin-3-one",
            ),
            (
                "O=C1CC=C2Nc3ccccc3SC2=C1",
                "2,10-dihydro-3H-phenothiazin-3-one",
            ),
            (
                "O=C1C=CC(=O)c2Nc3ccccc3Oc12",
                "1H-phenoxazine-1,4(10H)-dione",
            ),
            (
                "O=C1C(=O)C2=C(CC1)Nc1ccccc1O2",
                "1H-phenoxazine-3,4(2H,10H)-dione",
            ),
            (
                "O=C1C=CC2=C(C1=O)Oc1ccccc1N2",
                "3H-phenoxazine-3,4(10H)-dione",
            ),
            // P-58.2.3.1.4: a ketone that cannot take the indicated hydrogen
            // leaves it at the lowest permissible site (P-58.2.3.1.3).
            ("O=C1C=CC2=C(Nc3ccccc3O2)C1", "1H-phenoxazin-2(10H)-one"),
            (
                "O=C1CCC=C2Nc3ccccc3OC21",
                "2,10-dihydro-3H-phenoxazin-4(4aH)-one",
            ),
            ("O=C1C=CC2=Nc3ccccc3OC2=C1", "3H-phenoxazin-3-one"),
            // Controls.
            ("O=C1CNc2ccccc12", "1,2-dihydro-3H-indol-3-one"),
            ("O=C1CC=CC2=C1c1ccccc1N2", "3,9-dihydro-4H-carbazol-4-one"),
            ("O=c1[nH]cnc2nc[nH]c12", "1,7-dihydro-6H-purin-6-one"),
            (
                "CN1c2ccccc2Nc2ccncc21",
                "10-methyl-5,10-dihydropyrido[3,4-b]quinoxaline",
            ),
            (
                "O=S1c2ccccc2S(=O)c2ccccc21",
                "5lambda^4,10lambda^4-thianthrene-5,10-dione",
            ),
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, NameStatus::Preferred, "{smiles}");
            let parsed = crate::test_support::parse_strict(&name.text).unwrap();
            crate::test_support::verify_identity(smiles, &parsed).unwrap();
        }
        // Not composed yet, so these keep main's names labelled not preferred:
        // an ion whose ring also has a ketone cites it as 'oxo' (P-73.1.2.1),
        // an 'oxo' prefix needs its own indicated hydrogen, and an ion at a
        // hydro N takes added hydrogen (P-72.2.2.1.1).
        let fallback = NameStatus::NotPreferred(NameReason::VonBaeyerFusionFallback);
        let unregistered = NameStatus::NotPreferred(NameReason::UnregisteredRetainedComponent);
        for (smiles, status) in [
            ("O=C1CC=C2Oc3cc4Nc5ccccc5Sc4cc3[N+](C)(C)C2=C1", &fallback),
            ("O=C1CC=C2Oc3cc4[N+](C)(C)c5ccccc5Sc4cc3NC2=C1", &fallback),
            ("OC(=O)c1ccc2nc3c(cc(=O)c4ccccc34)oc2c1", &unregistered),
            ("C[N+]1(C)c2ccccc2Sc2cc3Oc4ccccc4N(C)c3cc21", &fallback),
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(&name.status, status, "{smiles}: {}", name.text);
        }
    }

    #[test]
    // Provenance: fn tests::fusion_retained_component_controls_remain_preferred | adapted
    // Source: ReShiki | src/naming/openclatura/mod.rs::tests::fusion_retained_component_controls_remain_preferred
    // Source: openclatura | tests/roundtrip/test_retained_fused_production.py::test_expected_retained_rings_are_registered
    // Source: openclatura | tests/roundtrip/test_fused_parents.py::test_roundtrip
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn fusion_retained_component_controls_remain_preferred() {
        use nazuchem_ir::NameStatus;
        for (smiles, expected) in [
            (
                "c1ccc2c(c1)Cc1cc3ncccc3cc1O2",
                "11H-[1]benzopyrano[2,3-g]quinoline",
            ),
            ("S1c2ccccc2Cc2ccccc21", "9H-thioxanthene"),
            ("c1ccc2c(c1)Cc1ccccc1O2", "9H-xanthene"),
            ("c1ccc2c(c1)Nc1ccccc1O2", "10H-phenoxazine"),
            ("c1ccc2c(c1)Nc1ccccc1S2", "10H-phenothiazine"),
            (
                "C1=CC(=CC=2SC3=CC=CC=C3CC12)C(=O)O",
                "9H-thioxanthene-3-carboxylic acid",
            ),
            (
                "N1=CC=CC2=CC=C3C=C4C(=NC3=C12)C=CC=C4",
                "benzo[b][1,10]phenanthroline",
            ),
            ("c1c3ncncc3cc4nc2ccccc2cc14", "pyrimido[4,5-b]acridine"),
            ("c1ccc2cc3nccnc3cc2c1", "benzo[g]quinoxaline"),
            ("O1C=CC2=NC=3C(C=C21)=COC3", "difuro[3,2-b:3',4'-e]pyridine"),
            ("c1ccc2c(c1)oc1ccccc12", "dibenzofuran"),
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, NameStatus::Preferred, "{smiles}");
        }
    }

    #[test]
    // Provenance: fn tests::fusion_intrinsic_carbon_hydrogen_relocation_is_preferred_and_round_trips | adapted
    // Source: ReShiki | src/naming/openclatura/mod.rs::tests::fusion_intrinsic_carbon_hydrogen_relocation_is_preferred_and_round_trips
    // Source: openclatura | tests/unit/test_retained_bridge_h_relocation.py::test_retained_wrapper_consumes_matched_h_metadata
    // Source: openclatura | tests/roundtrip/test_retained_fused_production.py::test_expected_retained_rings_are_registered
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn fusion_intrinsic_carbon_hydrogen_relocation_is_preferred_and_round_trips() {
        use nazuchem_ir::NameStatus;
        // Issue #969: P-25.3.2.4 selects the senior retained component;
        // P-25.3.1.3 and P-25.3.3 assign its descriptor and completed locants.
        // In the two pyrido cases, parent side b is component 2 -> 3,
        // giving [2,3-b] and [3,2-b], respectively.
        for (smiles, expected) in [
            ("C1C=Cc2cc3sccc3n12", "5H-thieno[2,3-b]pyrrolizine"),
            ("C1C=Cc2cc3ccsc3n12", "7H-thieno[3,2-b]pyrrolizine"),
            ("C1=CCc2c3ccsc3cn12", "4H-thieno[3,2-a]pyrrolizine"),
            (
                "CC1C=Cc2cc3sccc3n12",
                "5-methyl-5H-thieno[2,3-b]pyrrolizine",
            ),
            ("C1C=Cc2cc3occc3n12", "5H-furo[2,3-b]pyrrolizine"),
            ("C1C=Cc2cc3ncccc3n12", "6H-pyrido[2,3-b]pyrrolizine"),
            ("C1C=Cc2cc3cccnc3n12", "8H-pyrido[3,2-b]pyrrolizine"),
            ("c1ccc2c(c1)Cc1cccn12", "9H-pyrrolo[1,2-a]indole"),
            ("C1C=CN2C=CC=C12", "1H-pyrrolizine"),
            ("C1=CC2=CC=CN2C1", "3H-pyrrolizine"),
            ("C1CC2=CC=CN2C1", "2,3-dihydro-1H-pyrrolizine"),
            ("CC1=CC2=CC=CN2C1", "2-methyl-3H-pyrrolizine"),
            ("C1C=Cc2ccc3sccc3c2N1", "1,2-dihydrothieno[2,3-h]quinoline"),
            ("c1cc2cc3occc3n2cc1", "furo[2,3-b]indolizine"),
            ("c1cn2cnc3[nH]cnc3c2n1", "3H-imidazo[2,1-i]purine"),
            (
                "S1CC=C2C1=NC=1N2C=2C(=C3C1N=C1C=CC=CC1=C3)N=CC2",
                "2H-benzo[b]pyrrolo[2,3-f]thieno[2',3':4,5]imidazo[1,2-h][1,7]naphthyridine",
            ),
            ("c1ccc2c(c1)Cc1c(ccc3ccccc13)O2", "12H-benzo[a]xanthene"),
            ("C1C=Cn2c1cc1sccc12", "7H-thieno[2,3-b]pyrrolizine"),
            // Quinolizine also has a neutral trivalent donor N5. Its
            // component C4 (input atom 13) releases H to C6 (atom 11).
            // Parent 2 -> 3 and thiophene 2 -> 3 both map to input 5 -> 4,
            // so P-25.3.1.3 gives [2,3-b]. Completed input-order locants
            // are 1,2,3,3a,10a,10,9a,9,8,7,6,5,4 (P-25.3.3).
            ("S1C=CC=2C1=CC1=CC=CCN1C2", "6H-thieno[2,3-b]quinolizine"),
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, NameStatus::Preferred, "{smiles}");
            let parsed = crate::test_support::parse_strict(&name.text).unwrap();
            assert_eq!(
                crate::test_support::canonical_smiles(smiles).unwrap(),
                crate::test_support::canonical_smiles(&parsed).unwrap(),
                "{}",
                name.text,
            );
        }
    }

    #[test]
    fn fusion_registered_carbazole_keeps_existing_layout_status() {
        use nazuchem_ir::{NameReason, NameStatus};
        // Carbazole is registered, but this first-order fusion still uses an
        // OPSIN-entry layout; its existing numbering status remains unverified.
        let name = name_smiles("c1ccc2c(c1)ccc1c2[nH]c2ccccc12", &Default::default()).unwrap();
        assert_eq!(name.text, "11H-benzo[a]carbazole");
        assert_eq!(
            name.status,
            NameStatus::Unverified(NameReason::OpsinFusionLayout)
        );
    }

    #[test]
    fn compound_prefixes_inside_heteroatom_prefixes_are_enclosed() {
        // P-16.5.1.1, P-29.1.2, P-35.4.1: a compound prefix is enclosed within a larger
        // prefix, as in (chloromethyl)amino and (methoxyamino)phosphanyl; a simple one is not.
        for (smiles, expected) in [
            ("ClCSCC(=O)O", "[(chloromethyl)sulfanyl]acetic acid"),
            ("ClCNCC(=O)O", "[(chloromethyl)amino]acetic acid"),
            ("COCSCC(=O)O", "[(methoxymethyl)sulfanyl]acetic acid"),
            ("C1CC1CNCC(=O)O", "[(cyclopropylmethyl)amino]acetic acid"),
            ("ClC[Se]CC(=O)O", "[(chloromethyl)selanyl]acetic acid"),
            ("CONSCC(=O)O", "[(methoxyamino)sulfanyl]acetic acid"),
            (
                "ClCN=C1CCCC1C(=O)O",
                "2-[(chloromethyl)imino]cyclopentane-1-carboxylic acid",
            ),
            ("CSCC(=O)O", "(methylsulfanyl)acetic acid"),
            ("CNCC(=O)O", "(methylamino)acetic acid"),
            // Formerly declined: (cyclohexylperoxyamino) lost its carbamate contraction.
            ("COC(=O)NOOC1CCCCC1", "methyl (cyclohexylperoxy)carbamate"),
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, nazuchem_ir::NameStatus::Preferred, "{smiles}");
            let parsed = crate::test_support::parse_strict(&name.text).unwrap();
            assert_eq!(
                crate::test_support::canonical_smiles(smiles).unwrap(),
                crate::test_support::canonical_smiles(&parsed).unwrap(),
                "{expected}"
            );
        }
    }

    #[test]
    fn first_locantless_simple_prefix_is_not_enclosed() {
        // P-16.5.1.3.1, P-16.5.1.3.2: the first cited simple prefix is not enclosed
        // (a compound one is) and later ones are, as in cyclopropyl(phenyl)methanol and
        // hydroxydi(phenyl)acetic acid (P-65.1.1.2.4); a first prefix with a
        // locant keeps its parentheses. The tellanyl prefixes follow
        // hydroxy(oxo)-lambda^5-azanylidene (P-61.5.3) and
        // methyldi(oxo)-lambda^5-phosphane (P-64.1.2.2).
        for (smiles, expected) in [
            ("OC(C(=O)O)c1ccccc1", "hydroxy(phenyl)acetic acid"),
            (
                "OC(C(=O)O)(c1ccccc1)c1ccccc1",
                "hydroxydi(phenyl)acetic acid",
            ),
            (
                "OC(C(=O)O)c1ccccc1Cl",
                "(2-chlorophenyl)(hydroxy)acetic acid",
            ),
            (
                "OC(c1ccccc1)c1ccc(cc1)C(=O)O",
                "4-[hydroxy(phenyl)methyl]benzoic acid",
            ),
            ("OC(C1CC1)c1ccccc1", "cyclopropyl(phenyl)methanol"),
            (
                "OC(=O)CO[Te](=O)O",
                "{[hydroxy(oxo)-lambda^4-tellanyl]oxy}acetic acid",
            ),
            (
                "OC(=O)CO[Te](=O)(=O)O",
                "{[hydroxydi(oxo)-lambda^6-tellanyl]oxy}acetic acid",
            ),
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, nazuchem_ir::NameStatus::Preferred, "{smiles}");
            let parsed = crate::test_support::parse_strict(&name.text).unwrap();
            assert_eq!(
                crate::test_support::canonical_smiles(smiles).unwrap(),
                crate::test_support::canonical_smiles(&parsed).unwrap(),
                "{expected}"
            );
        }
        // A first compound prefix keeps its parentheses, as the compound heterone
        // prefixes oxophosphanyl, oxoarsanyl and oxostibanyl (P-61.6) do.
        for (smiles, expected) in [
            ("OC(=O)C(P=O)c1ccccc1", "(oxophosphanyl)(phenyl)acetic acid"),
            ("OC(=O)C([As]=O)c1ccccc1", "(oxoarsanyl)(phenyl)acetic acid"),
            (
                "OC(=O)C([Sb]=O)c1ccccc1",
                "(oxostibanyl)(phenyl)acetic acid",
            ),
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
        }
    }

    #[test]
    fn halogen_acyl_amino_and_triazene_prefixes_take_iupac_enclosing_marks() {
        // P-16.5.1.3.1, P-16.5.1.3.2: later simple prefixes are enclosed, halogens
        // included: bromo(chloro)acetic acid, (R)-bromo(chloro)(fluoro)methane
        // (P-91.3), bromo(chloro)(phenyl)borane (P-67.1.2.5.2). P-16.5.1.1: the compound
        // prefix methoxycarbonyl is enclosed inside amino. P-25.4.2.1.4 and
        // P-68.3.1.4.1: H2N-N=N- is triaz-1-en-1-yl, HN=N-NH- triaz-2-en-1-yl and
        // H2N-NH-NH- triazan-1-yl.
        for (smiles, expected) in [
            ("FC(Cl)Br", "bromo(chloro)(fluoro)methane"),
            ("F[C@H](Cl)Br", "(R)-bromo(chloro)(fluoro)methane"),
            ("OC(=O)C(Cl)Br", "bromo(chloro)acetic acid"),
            ("N#CC(F)C1CC1", "cyclopropyl(fluoro)acetonitrile"),
            ("ClC(Cl)(Br)C(=O)O", "bromodi(chloro)acetic acid"),
            ("ClC(Br)c1ccccc1", "[bromo(chloro)methyl]benzene"),
            ("ClC(Cl)C(=O)O", "dichloroacetic acid"),
            ("OC(=O)CNC(=O)OC", "[(methoxycarbonyl)amino]acetic acid"),
            (
                "OC(=O)C(NC(=O)OC)c1ccccc1",
                "[(methoxycarbonyl)amino](phenyl)acetic acid",
            ),
            (
                "CC(NC(=O)OC(C)(C)C)C(=O)O",
                "2-[(tert-butoxycarbonyl)amino]propanoic acid",
            ),
            ("OC(=O)CNS(=O)(=O)OC", "[(methoxysulfonyl)amino]acetic acid"),
            (
                "OC(=O)CN(C)C(=O)OC",
                "[(methoxycarbonyl)(methyl)amino]acetic acid",
            ),
            (
                "OC(=O)CN(CC)C(=O)OC",
                "[ethyl(methoxycarbonyl)amino]acetic acid",
            ),
            (
                "OC(=O)C(Cl)(Br)N=NN",
                "bromo(chloro)(triaz-1-en-1-yl)acetic acid",
            ),
            (
                "OC(=O)C(N=NN)c1ccccc1",
                "phenyl(triaz-1-en-1-yl)acetic acid",
            ),
            ("NN=Nc1ccc(cc1)C(=O)O", "4-(triaz-1-en-1-yl)benzoic acid"),
            ("OC(=O)CNN=N", "(triaz-2-en-1-yl)acetic acid"),
            ("OC(=O)CNNN", "(triazan-1-yl)acetic acid"),
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, nazuchem_ir::NameStatus::Preferred, "{smiles}");
            let parsed = crate::test_support::parse_strict(&name.text).unwrap();
            assert_eq!(
                crate::test_support::canonical_smiles(smiles).unwrap(),
                crate::test_support::canonical_smiles(&parsed).unwrap(),
                "{expected}"
            );
        }
    }

    #[test]
    // Provenance: fn tests::known_departures_keep_their_names_and_report_the_rule_family | adapted
    // Source: ReShiki | src/naming/openclatura/mod.rs::tests::known_departures_keep_their_names_and_report_the_rule_family
    // Source: openclatura | tests/unit/test_chalcogen_functional_groups.py::test_peroxide_linkages_preserve_both_chalcogen_sites
    // Source: openclatura | tests/unit/test_rules_tables.py::test_anhydride_half_citations_are_structured_rule_data
    // Source: openclatura | tests/unit/test_audit_reconstruction.py::test_self_audit_confirms_modelled_names
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn known_departures_keep_their_names_and_report_the_rule_family() {
        use nazuchem_ir::{NameReason as R, NameStatus};
        for (smiles, expected, reason) in [
            (
                "C1=CN=C2C=CC=C3C2=C1C=C1C2=C(C3)C=CC=C2C=C1",
                "17-azapentacyclo[10.7.1.1^{3,6}.0^{16,20}.0^{10,21}]henicosa-1(19),2,4,6(21),7,9,12(20),13,15,17-decaene",
                R::VonBaeyerFusionFallback,
            ),
            (
                "OC(=O)CS(=O)(=O)CS(=O)(=O)O",
                "(sulfomethylsulfonyl)acetic acid",
                R::SulfonylPrefix,
            ),
            (
                "CNC(=N)NCC(=O)O",
                "(3-methylguanidino)acetic acid",
                R::AminoPrefix,
            ),
            (
                "N/C(=N/[H])NCC(=O)O",
                "[(Z)-guanidino]acetic acid",
                R::AminoPrefix,
            ),
            (
                "CC[S@](C)=O",
                "(S)-ethyl methyl sulfoxide",
                R::FunctionalClass,
            ),
            (
                "CP(C)(F)=O",
                "fluorodi(methyl)phosphane oxide",
                R::FunctionalClass,
            ),
            ("CN(C)Cl", "chlorodi(methyl)azane", R::HydrideParent),
            ("CC=NN", "acetaldehyde hydrazone", R::FunctionalClass),
            (
                "ClC#[N+][O-]",
                "chloromethanenitrile oxide",
                R::UnretainedParent,
            ),
            (
                "CCC(C)=NN=CC",
                "butan-2-one ethan-1-ylidenehydrazone",
                R::FunctionalClass,
            ),
            ("CO[Se]C", "dimethyl selenoperoxide", R::FunctionalClass),
            ("CSSC", "1,2-dimethyldisulfane", R::HydrideParent),
            ("CCOSC", "ethoxy(methyl)sulfane", R::HydrideParent),
            ("CCCOSS", "1-propoxydisulfane", R::HydrideParent),
            (
                "CNC(=O)C(N)=O",
                "N-methylethanediamide",
                R::UnretainedParent,
            ),
            (
                "CN(C=O)C(=O)O",
                "(N-methylformamido)formic acid",
                R::UnretainedParent,
            ),
            ("O=C(O)Cl", "hydroxyformyl chloride", R::UnretainedParent),
            ("O=C(OC)Cl", "methoxyformyl chloride", R::UnretainedParent),
            (
                "CNC(=O)Cl",
                "1-chloro-N-methylformamide",
                R::UnretainedParent,
            ),
            (
                "O=CO[PH](=O)O",
                "(hydroxyphosphinoyloxy)methanal",
                R::UnretainedParent,
            ),
            (
                "O=C(O)OP(=O)(O)c1ccccn1",
                "{[hydroxy(pyridin-2-yl)phosphoryl]oxy}formic acid",
                R::UnretainedParent,
            ),
            (
                "COC(=O)N(C=O)c1ccsc1",
                "methyl [N-(thiophen-3-yl)formamido]formate",
                R::UnretainedParent,
            ),
            (
                "CCNC(=O)CC(=O)NCC",
                "N,N'-diethylpropanediamide",
                R::DiamideLocants,
            ),
            (
                "CNC(=O)CCC(N)=O",
                "N-methylbutanediamide",
                R::DiamideLocants,
            ),
            (
                "c1cc2cccc3conc3c-2c1",
                "azuleno[4,5-c]isoxazole",
                R::FusionComponent,
            ),
            // P-15.3.2.4.1, P-15.3.1.2.2: substituted units and concatenated
            // linking groups are not multiplied yet.
            // P-15.3.1.3: the unit numbering that puts the linking group
            // at locant 1 is not built, so the substitutive name stays.
            (
                "Brc1ccc(Oc2ccc(Br)c(Br)c2)cc1Br",
                "1,2-dibromo-4-(3,4-dibromophenoxy)benzene",
                R::Multiplicative,
            ),
            (
                "O=C(O)c1ccc(cc1)C(=O)C1CCCCC1",
                "4-(cyclohexylcarbonyl)benzoic acid",
                R::AcylPrefix,
            ),
            (
                "O=C(O)c1ccc(cc1)C(=S)C1CCCCC1",
                "4-(cyclohexylcarbonothioyl)benzoic acid",
                R::AcylPrefix,
            ),
            (
                "O=C(O)c1ccc(cc1)P(C)(C)=O",
                "4-(dimethyloxophosphanyl)benzoic acid",
                R::AcylPrefix,
            ),
            ("C[PH](C)=O", "(methyloxophosphanyl)methane", R::AcylPrefix),
            // P-67.1.4.1.1.5: the preferred prefix follows
            // hydroxy(methylphosphonoyl).
            (
                "O[As](C)(=O)CC(=O)O",
                "[hydroxy(methyl)arsoryl]acetic acid",
                R::AcylPrefix,
            ),
            (
                "CO[As](C)(=O)CC(=O)O",
                "[methoxy(methyl)arsoryl]acetic acid",
                R::AcylPrefix,
            ),
            (
                "O=C(O)c1ccc([As](=O)(C)C)cc1",
                "4-(dimethyloxoarsanyl)benzoic acid",
                R::AcylPrefix,
            ),
            (
                "CC[C@@H](C)c1ccc(Sc2ccc(cc2)[C@H](C)CC)cc1",
                "1-[(2R)-butan-2-yl]-4-({4-[(2R)-butan-2-yl]phenyl}sulfanyl)benzene",
                R::Multiplicative,
            ),
            // P-51.4.1.1: four heterounits; the skeletal replacement name
            // is not built.
            (
                "OC(=O)COCCOCCOCCOCC(=O)O",
                "(2-{2-[2-(carboxymethoxy)ethoxy]ethoxy}ethoxy)acetic acid",
                R::Multiplicative,
            ),
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, NameStatus::NotPreferred(reason), "{smiles}");
            if expected.ends_with("acetic acid") || expected.ends_with("formic acid") {
                let parsed = crate::test_support::parse_strict(&name.text).unwrap();
                assert_eq!(
                    crate::test_support::canonical_smiles(smiles).unwrap(),
                    crate::test_support::canonical_smiles(&parsed).unwrap(),
                    "{}",
                    name.text
                );
            }
        }
    }

    #[test]
    // Provenance: fn tests::unchecked_paths_keep_their_names_and_report_the_uncertainty | adapted
    // Source: ReShiki | src/naming/openclatura/mod.rs::tests::unchecked_paths_keep_their_names_and_report_the_uncertainty
    // Source: openclatura | tests/integration/test_systematic_fusion.py::test_retained_complete_system_matrix_precedes_systematic_fusion
    // Source: openclatura | tests/unit/test_fusion_wrappers.py::test_bridged_fusion_wrapper_uses_retained_parent_locants
    // Source: openclatura | tests/unit/test_fusion_wrappers.py::test_bridge_wrapper_composes_with_shared_parent_derivative_operations
    // Source: openclatura | tests/unit/test_analysis.py::test_public_api_golden_names
    // Source: openclatura | tests/roundtrip/test_p13_fusion.py::test_roundtrip
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn unchecked_paths_keep_their_names_and_report_the_uncertainty() {
        use nazuchem_ir::{NameReason as R, NameStatus};
        for (smiles, expected, reason) in [
            ("NCC=O", "2-aminoethanal", R::SubstitutedRetainedParent),
            // P-52.2.8: with one principal group each, the ring is the parent.
            (
                "OC(=O)Cc1ccc(cc1)C(=O)O",
                "4-(carboxymethyl)benzoic acid",
                R::PrincipalGroupTie,
            ),
            (
                "NC(=O)CNC(=O)c1ccccc1",
                "N-(2-amino-2-oxoethyl)benzamide",
                R::PrincipalGroupTie,
            ),
            (
                "C1=Cc2cc3ccc(cc4nc(cc5ccc(cc1n2)[nH]5)C=C4)[nH]3",
                "porphyrin",
                R::UnauditedParent,
            ),
            (
                "C1=C2CCC(=N2)C=C2CCC(N2)C2CCC(=N2)C=C2CCC1=N2",
                "corrin",
                R::UnauditedParent,
            ),
            (
                "C12C(C=C(C3=CC=CC=C13)O2)=O",
                "1,4-epoxynaphthalen-2(1H)-one",
                R::AtomOrder,
            ),
            (
                "C12=CC=C(C3=CC=CC=C13)O2",
                "1,4-epoxynaphthalene",
                R::BridgedFusedParent,
            ),
            // P-44.4.1.12.2: three centres, and RS against SR, need CIP
            // rules 4 and 5.
            (
                "CCC(N[C@H](C)[C@@H](Cl)[C@@H](C)O)(N[C@@H](C)[C@@H](Cl)[C@@H](C)O)CCC",
                "(2R,3R,4R)-3-chloro-4-[(3-{[(2S,3R,4R)-3-chloro-4-hydroxypentan-2-yl]amino}hexan-3-yl)amino]pentan-2-ol",
                R::ParentConfiguration,
            ),
            (
                "CCC(N[C@H](C)[C@H](C)O)(N[C@@H](C)[C@@H](C)O)CCC",
                "(2S,3R)-3-[(3-{[(2S,3R)-3-hydroxybutan-2-yl]amino}hexan-3-yl)amino]butan-2-ol",
                R::ParentConfiguration,
            ),
            (
                "c1oc2cc3nc4ccccc4nc3cc2c1",
                "furo[2,3-b]phenazine",
                R::OpsinFusionLayout,
            ),
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, NameStatus::Unverified(reason), "{smiles}");
            if expected.ends_with("benzoic acid") {
                let parsed = crate::test_support::parse_strict(&name.text).unwrap();
                assert_eq!(
                    crate::test_support::canonical_smiles(smiles).unwrap(),
                    crate::test_support::canonical_smiles(&parsed).unwrap(),
                    "{}",
                    name.text
                );
            }
        }
    }

    #[test]
    // Provenance: fn tests::retained_and_systematic_fused_parents_follow_parent_precedence | adapted
    // Source: ReShiki | src/naming/openclatura/mod.rs::tests::retained_and_systematic_fused_parents_follow_parent_precedence
    // Source: openclatura | tests/roundtrip/test_fused_parents.py::test_fused_parent_suite_pin_claims_are_explicit_and_evidence_backed
    // Source: openclatura | tests/roundtrip/test_fused_hydrocarbon_retained.py::test_hydrogenated_parent_uses_data_backed_preferred_name
    // Source: openclatura | tests/integration/test_systematic_fusion.py::test_retained_complete_system_matrix_precedes_systematic_fusion
    // Source: openclatura | tests/integration/test_systematic_fusion.py::test_multiple_indicated_hydrogens_are_graph_derived_and_roundtrip
    // Source: openclatura | tests/integration/test_systematic_fusion.py::test_positive_heteroatom_charge_is_audited_as_a_locanted_parent_operation
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn retained_and_systematic_fused_parents_follow_parent_precedence() {
        for (smiles, expected) in [
            ("c1ccc2ccccc2c1", "naphthalene"),
            ("c1ccc2cc3ccccc3cc2c1", "anthracene"),
            ("c1ccc2[nH]ccc2c1", "1H-indole"),
            ("c1ccc2c(c1)[nH]c1ccccc12", "9H-carbazole"),
            ("C1CCC2=CC=CC=C2C1", "1,2,3,4-tetrahydronaphthalene"),
            ("c1ccc2c(c1)CCN2", "2,3-dihydro-1H-indole"),
            ("c1oc2ncccc2c1", "furo[2,3-b]pyridine"),
            ("C1OCC2COCC12", "tetrahydro-1H,3H-furo[3,4-c]furan"),
        ] {
            assert_eq!(
                name_smiles(smiles, &Default::default())
                    .map(|name| name.text)
                    .as_deref(),
                Ok(expected),
                "{smiles}"
            );
        }
        assert_eq!(
            name_smiles("C1CCC2(C1)Oc1ccccc1O2", &Default::default())
                .map(|name| name.text)
                .as_deref(),
            Ok("spiro[[1,3]benzodioxole-2,1'-cyclopentane]"),
        );
    }

    #[test]
    fn saturated_fused_parents_cite_only_the_indicated_hydrogen_the_mancude_parent_needs() {
        // Issue #87. P-14.7.1: indicated hydrogen marks the positions the
        // mancude parent leaves without a double bond; P-31.2.2: hydro
        // prefixes saturate double bonds, so they add hydrogen in pairs.
        for (smiles, expected) in [
            // Nine ring atoms: one indicated hydrogen, not 1H,6H plus hexahydro.
            ("N1CCCC2C1CNC2", "octahydro-1H-pyrrolo[3,4-b]pyridine"),
            (
                "N1CCC(C)C2C1CNC2",
                "4-methyloctahydro-1H-pyrrolo[3,4-b]pyridine",
            ),
            // Eleven ring atoms.
            (
                "N1CCC2C1C1CCNC1C2",
                "decahydro-1H-cyclopenta[1,2-b:4,5-b']dipyrrole",
            ),
            // Even ring sizes need no indicated hydrogen.
            ("C1NCC2C1CC1CNCC1C2", "dodecahydropyrrolo[3,4-f]isoindole"),
            // Controls that were already correct.
            ("N1CCCC2C1NCC2", "octahydro-1H-pyrrolo[2,3-b]pyridine"),
            ("N1CCC2C1CNC2", "octahydropyrrolo[3,4-b]pyrrole"),
            ("N1CCCC2C1CNCC2", "decahydro-1,7-naphthyridine"),
            // Issue #108. P-14.7.1, P-58.2.1.2: a fusion N takes no double
            // bond, so it can split the mancude parent into parts that each
            // need indicated hydrogen, also in an even-membered system
            // ('1H,5H-pyrido[3,2,1-ij]quinoline' in P-58.2.3.1.4).
            ("C1CN2CCCN2C1", "tetrahydro-1H,5H-pyrazolo[1,2-a]pyrazole"),
            (
                "C1CC2CCCC3CCCN(C1)C23",
                "decahydro-1H,5H-pyrido[3,2,1-ij]quinoline",
            ),
            (
                "C1CC2=C3C(=CC=C2)CCCN3C1",
                "2,3,6,7-tetrahydro-1H,5H-pyrido[3,2,1-ij]quinoline",
            ),
            // The fusion N leaves a nine-atom chain; position 1 pairs up, so
            // the lowest indicated hydrogen is 2H (P-58.2.1.2).
            ("C1CCN2CCNCC2C1", "octahydro-2H-pyrido[1,2-a]pyrazine"),
            ("C1CN2CCNCC2CN1", "octahydro-2H-pyrazino[1,2-a]pyrazine"),
            // Controls with a fusion N that were already correct.
            ("C1CCN2CCCN2C1", "hexahydro-1H-pyrazolo[1,2-a]pyridazine"),
            ("C1CCN2CCCCN2C1", "octahydropyridazino[1,2-a]pyridazine"),
            ("C1CCN2CCCCC2C1", "octahydro-2H-quinolizine"),
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, nazuchem_ir::NameStatus::Preferred, "{smiles}");
            let parsed = crate::test_support::parse_strict(&name.text).unwrap();
            assert_eq!(
                crate::test_support::canonical_smiles(smiles).unwrap(),
                crate::test_support::canonical_smiles(&parsed).unwrap(),
                "{expected}"
            );
        }
    }

    #[test]
    fn polyhelicenes_take_their_retained_names() {
        for (smiles, expected) in [
            // P-25.1.2.6: six or more helically ortho-fused benzene rings.
            ("c1ccc2c(c1)ccc1ccc3ccc4ccc5ccccc5c4c3c12", "hexahelicene"),
            (
                "c1ccc2c(c1)ccc1ccc3ccc4ccc5ccc6ccccc6c5c4c3c12",
                "heptahelicene",
            ),
            (
                "c1ccc2c(c1)ccc1ccc3ccc4ccc5ccc6ccc7ccc8ccc9ccccc9c8c7c6c5c4c3c12",
                "decahelicene",
            ),
            // P-25.3.3.1.1: numbering starts at the inner CH of a terminal ring.
            (
                "Clc1cccc2ccc3ccc4ccc5ccc6ccccc6c5c4c3c12",
                "1-chlorohexahelicene",
            ),
            (
                "c1ccc2c(c1)ccc1ccc3ccc4ccc5ccc(Cl)cc5c4c3c12",
                "2-chlorohexahelicene",
            ),
            (
                "c1ccc2c(c1)ccc1ccc3ccc4cc(Cl)c5ccccc5c4c3c12",
                "5-chlorohexahelicene",
            ),
            (
                "c1ccc2c(c1)ccc1ccc3c(Cl)cc4ccc5ccccc5c4c3c12",
                "8-chlorohexahelicene",
            ),
            // Hydro prefixes are kept at every ring count.
            (
                "C1Cc2ccc3ccc4ccc5ccc6ccccc6c5c4c3c2C=C1",
                "3,4-dihydrohexahelicene",
            ),
            (
                "C1Cc2ccc3ccc4ccc5ccc6ccc7ccc8ccc9ccc%10ccccc%10c9c8c7c6c5c4c3c2C=C1",
                "3,4-dihydrodecahelicene",
            ),
            // Five rings are not a helicene.
            (
                "c1ccc2c(c1)ccc1ccc3ccc4ccccc4c3c12",
                "dibenzo[c,g]phenanthrene",
            ),
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, nazuchem_ir::NameStatus::Preferred, "{smiles}");
            let parsed = crate::test_support::parse_strict(&name.text).unwrap();
            assert_eq!(
                crate::test_support::canonical_smiles(smiles).unwrap(),
                crate::test_support::canonical_smiles(&parsed).unwrap(),
                "{expected}"
            );
        }
    }

    #[test]
    fn substituted_pyrrolic_nitrogen_keeps_systematic_fusion_indicated_hydrogen() {
        // #16, #29: a PIN always cites the parent hydride's indicated hydrogen
        // (P-14.7.1, P-58.2.1.1); substituent prefixes are added after it
        // (P-58.2.4), so a substituted N-H keeps it, as in 1-methyl-1H-indole.
        for (smiles, expected) in [
            ("Cn1ccc2sccc21", "4-methyl-4H-thieno[3,2-b]pyrrole"),
            ("Cn1cnc2ccncc21", "3-methyl-3H-imidazo[4,5-c]pyridine"),
            (
                "Cn1nnc2ccncc21",
                "3-methyl-3H-[1,2,3]triazolo[4,5-c]pyridine",
            ),
            ("Cn1ccc2ccc3ccccc3c21", "1-methyl-1H-benzo[g]indole"),
            ("Cn1ccc2c1ccc1ccccc12", "3-methyl-3H-benzo[e]indole"),
            ("Cn1c2ccncc2c2ccccc21", "5-methyl-5H-pyrido[4,3-b]indole"),
            ("Cn1c2cccnc2c2ccccc21", "5-methyl-5H-pyrido[3,2-b]indole"),
            ("CN1N=CC=2C1=NC=CC2", "1-methyl-1H-pyrazolo[3,4-b]pyridine"),
            ("CN1C=CC=2C=NC=CC21", "1-methyl-1H-pyrrolo[3,2-c]pyridine"),
            (
                "O=C1C=Cc2ccn(C)c21",
                "1-methylcyclopenta[b]pyrrol-6(1H)-one",
            ),
            // Controls the issues report as already correct.
            ("Cn1ccc2ncccc21", "1-methyl-1H-pyrrolo[3,2-b]pyridine"),
            ("Cn1c2cnccc2c2ccccc21", "9-methyl-9H-pyrido[3,4-b]indole"),
            (
                "C[NH+]1C=CC=2C=NC=CC21",
                "1-methyl-1H-pyrrolo[3,2-c]pyridin-1-ium",
            ),
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, nazuchem_ir::NameStatus::Preferred, "{smiles}");
            let parsed = crate::test_support::parse_strict(&name.text).unwrap();
            assert_eq!(
                crate::test_support::canonical_smiles(smiles).unwrap(),
                crate::test_support::canonical_smiles(&parsed).unwrap(),
                "{expected}"
            );
        }
    }

    #[test]
    fn fused_parents_place_indicated_hydrogen_before_hydro_prefixes_and_on_ketones() {
        use nazuchem_ir::NameStatus;
        for (smiles, expected) in [
            // #95. P-31.2.2, P-58.2.1.2: indicated hydrogen takes lower locants
            // than hydro prefixes ('2,3-dihydro-1H-indene', not
            // '1,3-dihydro-2H-indene'; '3,4-dihydro-2H-1-benzopyran').
            ("C1NCc2ncccc21", "6,7-dihydro-5H-pyrrolo[3,4-b]pyridine"),
            (
                "C1CNCc2c1[nH]c1ccccc12",
                "2,3,4,5-tetrahydro-1H-pyrido[4,3-b]indole",
            ),
            ("C1CC2OCC=CC2O1", "3,3a,5,7a-tetrahydro-2H-furo[3,2-b]pyran"),
            // The same for retained parents whose name states the hydrogen.
            (
                "C1CNCc2c1c1ccccc1[nH]2",
                "2,3,4,9-tetrahydro-1H-pyrido[3,4-b]indole",
            ),
            (
                "C1CCc2c(C1)Oc1ccccc1N2",
                "2,3,4,10-tetrahydro-1H-phenoxazine",
            ),
            ("C1CSc2ccccc2N1", "3,4-dihydro-2H-1,4-benzothiazine"),
            ("C1CNc2ncncc21", "6,7-dihydro-5H-pyrrolo[2,3-d]pyrimidine"),
            // P-58.2.4: prefixes are added after the hydrogen of the N-H
            // parent, whether that N-H is indicated or hydro hydrogen (#52).
            (
                "CN1Cc2ncccc2C1",
                "6-methyl-6,7-dihydro-5H-pyrrolo[3,4-b]pyridine",
            ),
            (
                "CN1CCc2c(C1)c1ccccc1n2C",
                "2,5-dimethyl-2,3,4,5-tetrahydro-1H-pyrido[4,3-b]indole",
            ),
            (
                "CN1CCc2sccc21",
                "4-methyl-5,6-dihydro-4H-thieno[3,2-b]pyrrole",
            ),
            (
                "CN1Cc2sccc2C1",
                "5-methyl-5,6-dihydro-4H-thieno[2,3-c]pyrrole",
            ),
            (
                "CN1CCC2=C1C=NN2",
                "4-methyl-1,4,5,6-tetrahydropyrrolo[3,2-c]pyrazole",
            ),
            (
                "N1CCCC2C1CN(C)C2",
                "6-methyloctahydro-1H-pyrrolo[3,4-b]pyridine",
            ),
            (
                "CN1CCCC2C1CN(C)C2",
                "1,6-dimethyloctahydro-1H-pyrrolo[3,4-b]pyridine",
            ),
            ("CN1CC2CCNC2C1", "5-methyloctahydropyrrolo[3,4-b]pyrrole"),
            // P-72.2.2.1, P-73.1.1.2: '-ide' and '-ium' modify the parent
            // hydride, which keeps its hydrogen ('5H-...-5-ide' in P-74.1.1).
            ("c1ccc2c(c1)[n-]c1ccncc12", "5H-pyrido[4,3-b]indol-5-ide"),
            ("[n-]1ccc2sccc21", "4H-thieno[3,2-b]pyrrol-4-ide"),
            (
                "C[NH+]1Cc2ncccc2C1",
                "6-methyl-6,7-dihydro-5H-pyrrolo[3,4-b]pyridin-6-ium",
            ),
            (
                "C1C[NH2+]Cc2c1c1ccccc1[nH]2",
                "2,3,4,9-tetrahydro-1H-pyrido[3,4-b]indol-2-ium",
            ),
            (
                "C[NH+]1CCc2c(C1)[nH]c1ccccc21",
                "2-methyl-2,3,4,9-tetrahydro-1H-pyrido[3,4-b]indol-2-ium",
            ),
            (
                "C1C[NH2+]Cc2c1[nH]c1ccccc12",
                "2,3,4,5-tetrahydro-1H-pyrido[4,3-b]indol-2-ium",
            ),
            // P-58.2.3.1.1: with as many indicated hydrogens as ketones the
            // ketone takes the hydrogen and the other saturated positions are
            // hydro positions ('7H-1-benzopyran-7-one'; P-64.3.1
            // '1,2-dihydro-3H-indol-3-one, not 1H-indol-3(2H)-one').
            ("O=c1ccoc2ncccc12", "4H-pyrano[2,3-b]pyridin-4-one"),
            ("O=c1ccc2cccnc2o1", "2H-pyrano[2,3-b]pyridin-2-one"),
            ("O=c1ccoc2ccsc12", "7H-thieno[3,2-b]pyran-7-one"),
            (
                "O=C1NCc2ncccc21",
                "6,7-dihydro-5H-pyrrolo[3,4-b]pyridin-5-one",
            ),
            (
                "O=C1CNc2[nH]c3ccccc3c2C1",
                "1,2,4,9-tetrahydro-3H-pyrido[2,3-b]indol-3-one",
            ),
            (
                "O=c1[nH]ccc2[nH]ccc12",
                "1,5-dihydro-4H-pyrrolo[3,2-c]pyridin-4-one",
            ),
            (
                "O=C1CN=C(c2ccccc2)c2cc(Cl)ccc2N1",
                "7-chloro-5-phenyl-1,3-dihydro-2H-1,4-benzodiazepin-2-one",
            ),
            (
                "O=c1[nH]ccc2c1[nH]c1ccccc12",
                "2,9-dihydro-1H-pyrido[3,4-b]indol-1-one",
            ),
            (
                "O=C1NCCc2c1[nH]c1ccccc12",
                "2,3,4,9-tetrahydro-1H-pyrido[3,4-b]indol-1-one",
            ),
            ("O=C1CCCc2c[nH]cc21", "2,5,6,7-tetrahydro-4H-isoindol-4-one"),
            (
                "O=c1[nH]cnc2[nH]ccc12",
                "3,7-dihydro-4H-pyrrolo[2,3-d]pyrimidin-4-one",
            ),
            (
                "O=c1[nH]cnc2c1cnn2C",
                "1-methyl-1,5-dihydro-4H-pyrazolo[3,4-d]pyrimidin-4-one",
            ),
            (
                "O=C1CC=C2Nc3ccccc3OC2=C1",
                "2,10-dihydro-3H-phenoxazin-3-one",
            ),
            // P-58.2.3.1.3: with fewer indicated hydrogens than ketones it
            // takes the lowest position and the rest use added hydrogen
            // ('1H-cyclopenta[a]naphthalene-1,2(3H)-dione').
            (
                "O=C1NC(=O)c2ncccc12",
                "5H-pyrrolo[3,4-b]pyridine-5,7(6H)-dione",
            ),
            (
                "O=C1C=CC(=O)c2Nc3ccccc3Oc12",
                "1H-phenoxazine-1,4(10H)-dione",
            ),
            // P-58.2.3.1.4: a ketone position that cannot take the hydrogen
            // leaves it at the lowest permissible one.
            ("O=C1C=CC2=C(Nc3ccccc3O2)C1", "1H-phenoxazin-2(10H)-one"),
            ("O=C1COc2ccccc2N1", "2H-1,4-benzoxazin-3(4H)-one"),
            // Controls: already lowest, a hydro N, no indicated hydrogen in
            // the parent, and parents whose name does not state the hydrogen.
            (
                "C1CCc2c(C1)n(C)c1ccncc21",
                "5-methyl-6,7,8,9-tetrahydro-5H-pyrido[4,3-b]indole",
            ),
            (
                "CN1CCCc2sccc21",
                "4-methyl-4,5,6,7-tetrahydrothieno[3,2-b]pyridine",
            ),
            (
                "O=c1c2ccccc2nc2n1CCC2",
                "2,3-dihydropyrrolo[2,1-b]quinazolin-9(1H)-one",
            ),
            (
                "C1CCc2c(C1)[nH]c1ccccc12",
                "2,3,4,9-tetrahydro-1H-carbazole",
            ),
            ("O=C1CC=CC2=C1c1ccccc1N2", "3,9-dihydro-4H-carbazol-4-one"),
            ("O=c1[nH]cnc2nc[nH]c12", "1,7-dihydro-6H-purin-6-one"),
            (
                "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
                "1,3,7-trimethyl-3,7-dihydro-1H-purine-2,6-dione",
            ),
            ("O=C1NC(=O)c2ccccc12", "1H-isoindole-1,3(2H)-dione"),
            // Issue #95 leftovers. P-25.7.1.1, P-14.7.1: a fusion N takes no
            // double bond, so the parent's indicated hydrogen is counted
            // without it.
            // P-58.2.3.1.4, the Blue Book example: neither ketone position
            // can take the hydrogen, and the two ketones remove one double
            // bond, so no added hydrogen is cited (P-58.2.2.3).
            (
                "C1C=CN2CC(=O)C=C3C(=O)C=CC1=C23",
                "1H,5H-pyrido[3,2,1-ij]quinoline-6,8-dione",
            ),
            (
                "O=C1C=C2C(=O)C=CC3=C2N(C1)C=CC3",
                "1H,5H-pyrido[3,2,1-ij]quinoline-6,8-dione",
            ),
            // P-58.2.3.1.4 with P-58.2.3.1.3 (1)-(3): 5,6 cannot both take
            // the hydrogen; one stays at the lowest position (3), the other
            // goes to ketone 5, and ketone 6 takes added hydrogen.
            (
                "O=C1C(=O)CN2N=CCN12",
                "3H,5H-pyrazolo[1,2-a][1,2,3]triazole-5,6(7H)-dione",
            ),
            // P-58.2.3.1.1: as many indicated hydrogens as ketones.
            (
                "O=C1C=CN2C=CC(=O)N12",
                "1H,7H-pyrazolo[1,2-a]pyrazole-1,7-dione",
            ),
            ("O=C1C=CN2C=CC=CC2=N1", "2H-pyrido[1,2-a]pyrimidin-2-one"),
            ("O=C1C=CN=C2C=CC=CN12", "4H-pyrido[1,2-a]pyrimidin-4-one"),
            (
                "O=C1C=CN2C=CSC2=N1",
                "7H-[1,3]thiazolo[3,2-a]pyrimidin-7-one",
            ),
            // P-58.2.3.1.2: the hydrogen the ketone does not use takes the
            // lowest locant.
            (
                "C1CN2CCCN2C1=O",
                "tetrahydro-1H,5H-pyrazolo[1,2-a]pyrazol-1-one",
            ),
            // P-58.2.3.1.4: position 3 cannot take the hydrogen (S and the
            // fusion N isolate C2-C3), so it stays at 5 and the ketone takes
            // added hydrogen.
            (
                "O=C1CSC2=NC=CCN12",
                "5H-[1,3]thiazolo[3,2-a]pyrimidin-3(2H)-one",
            ),
            ("O=C1CCCN2C=CC=CC12", "3,4-dihydro-2H-quinolizin-1(9aH)-one"),
            ("O=C1C=CC=C2C=CC=CN12", "4H-quinolizin-4-one"),
            // P-58.2.4: groups cited as prefixes replace hydrogen of hydro
            // positions; only a pair of suffixes removes the double bond
            // without them (P-58.2.2.3).
            (
                "COC(=O)C1=CN=C2SC(=C)C(=O)N2C1",
                "methyl 2-methylidene-3-oxo-2,3-dihydro-5H-[1,3]thiazolo[3,2-a]pyrimidine-6-carboxylate",
            ),
            (
                "C=C1SC2=NC=CCN2C1=C",
                "2,3-dimethylidene-2,3-dihydro-5H-[1,3]thiazolo[3,2-a]pyrimidine",
            ),
            (
                "O=C1SC2=NC=CCN2C1=O",
                "5H-[1,3]thiazolo[3,2-a]pyrimidine-2,3-dione",
            ),
            // P-58.2.3.1.1: 'locants for hydro prefixes are those of the
            // saturated positions', also when the double bonds between them
            // move ('1,2-dihydro-3H-indol-3-one' in P-64.3.1).
            ("O=C1CC=CC2=CNC=C12", "2,5-dihydro-4H-isoindol-4-one"),
            ("O=C1CC=Cc2c[nH]cc12", "2,5-dihydro-4H-isoindol-4-one"),
            // P-58.2.2.2: without indicated hydrogen in the parent, the added
            // hydrogen takes the lowest position the double bonds allow.
            (
                "O=C1CCc2c1[nH]c1ccccc12",
                "1,4-dihydrocyclopenta[b]indol-3(2H)-one",
            ),
            (
                "O=C1CCc2c1c1ccccc1[nH]2",
                "3,4-dihydrocyclopenta[b]indol-1(2H)-one",
            ),
            // P-58.2.3.1.4: 2H-perimidine has one double bond fewer than
            // 1H-perimidine, so the ketone cannot take the hydrogen.
            ("O=C1Nc2cccc3cccc(N1)c23", "1H-perimidin-2(3H)-one"),
            // P-25.7.1.1 ('pyrrolo[3,2-b]pyrrole'): mancude
            // pyrrolo[2,3-b]indole has six double bonds (N1=C2, C3=C3a,
            // N8=C8a and the benzo ring) and no indicated hydrogen; both N-H
            // are hydro positions (P-31.2.2).
            (
                "c1ccc2c(c1)[nH]c1[nH]ccc12",
                "1,8-dihydropyrrolo[2,3-b]indole",
            ),
            (
                "CN1c2ccccc2C2=C1NC(=O)C2",
                "8-methyl-3,8-dihydropyrrolo[2,3-b]indol-2(1H)-one",
            ),
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, NameStatus::Preferred, "{smiles}");
            let parsed = crate::test_support::parse_strict(&name.text).unwrap();
            crate::test_support::verify_identity(smiles, &parsed).unwrap();
        }
        // A bridged fused parent replays the plan, so the relocated hydrogen
        // must carry its own bond proof; otherwise the plan is left as it was.
        for (smiles, expected) in [
            (
                "C1NCc2nc3cc(C3)c21",
                "6,7-dihydro-2,4-methano-5H-pyrrolo[3,4-b]pyridine",
            ),
            (
                "C1NCc2nc3cc(CC3)c21",
                "6,7-dihydro-2,4-ethano-5H-pyrrolo[3,4-b]pyridine",
            ),
            (
                "C1CNCc2c1[nH]c1cc3cc(C3)c12",
                "1,2,3,4-tetrahydro-7,9-methano-5H-pyrido[4,3-b]indole",
            ),
            (
                "CN1CCC23C4C1CC5=C2C(=CC=C5)OC3CC=C4",
                "3-methyl-2,3,4,4a,7,7a-hexahydro-4,12-methano-1H-[1]benzofuro[3,2-e]isoquinoline",
            ),
            // A bridge on the ring N is not a substituent of the N-H parent.
            (
                "C1CN3Cc2c1[nH]c1cc(CC3)ccc12",
                "1,2,3,4-tetrahydro-2,7-ethano-5H-pyrido[4,3-b]indole",
            ),
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            let parsed = crate::test_support::parse_strict(&name.text).unwrap();
            crate::test_support::verify_identity(smiles, &parsed).unwrap();
        }
    }

    #[test]
    fn angular_oxazole_thiazole_acetic_acid_round_trips() {
        let smiles = "O=C(O)Cc1nc2ccc3ocnc3c2s1";
        let name = name_smiles(smiles, &Default::default()).unwrap().text;
        assert_eq!(
            name,
            "([1,3]thiazolo[5,4-e][1,3]benzoxazol-7-yl)acetic acid"
        );
        let parsed = crate::test_support::parse_strict(&name).unwrap();
        assert_eq!(
            crate::test_support::canonical_smiles(smiles).unwrap(),
            crate::test_support::canonical_smiles(&parsed).unwrap()
        );
    }

    #[test]
    // Provenance: fn tests::spiro_parents_and_substituents_follow_oracle_names | adapted
    // Source: ReShiki | src/naming/openclatura/mod.rs::tests::spiro_parents_and_substituents_follow_oracle_names
    // Source: openclatura | tests/unit/test_fused_bonds.py::test_fused_bonds_match_the_cycle_count
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn spiro_parents_and_substituents_follow_oracle_names() {
        for (smiles, expected) in [
            ("C1CCC2(CC1)CCCC2", "spiro[4.5]decane"),
            ("C1CCC2(CC1)OCCO2", "1,4-dioxaspiro[4.5]decane"),
            ("O=C1CCCC12CCCCC2", "spiro[4.5]decan-1-one"),
            ("C1=CC2(C1)CC1(CCC1)C2", "dispiro[3.1.3.1]dec-1-ene"),
            ("C1=CC12NC21CC1", "7-azadispiro[2.0.2.1]hept-1-ene"),
            (
                "CN/C=C\\C1CC2(C1)CCC2",
                "(1Z)-N-methyl-2-(spiro[3.3]heptan-2-yl)ethen-1-amine",
            ),
            (
                "CC1=CC=C(C=C1)N2C(=O)C[C@@]3(C2=O)C=C(NN3)C(=O)OCC4=CC=CC=C4",
                "benzyl (5S)-7-(4-methylphenyl)-6,8-dioxo-1,2,7-triazaspiro[4.4]non-3-ene-3-carboxylate",
            ),
        ] {
            assert_eq!(
                name_smiles(smiles, &Default::default())
                    .map(|name| name.text)
                    .as_deref(),
                Ok(expected),
                "{smiles}"
            );
        }
    }

    #[test]
    fn unsupported_spiro_topologies_decline() {
        for smiles in ["C1CC12CC23CC34CC4", "S123(CC1)(CC2)CC3"] {
            assert!(
                name_smiles(smiles, &Default::default()).is_err(),
                "{smiles}"
            );
        }
    }

    #[test]
    // Provenance: fn tests::declines_unrendered_group_and_attachment_features | adapted
    // Source: ReShiki | src/naming/openclatura/mod.rs::tests::declines_unrendered_group_and_attachment_features
    // Source: openclatura | tests/unit/test_chalcogen_functional_groups.py::test_polonium_is_parseable_but_has_an_explicit_unvalidated_boundary
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn declines_unrendered_group_and_attachment_features() {
        for smiles in [
            "CC(=[Po])[OH]",
            "C/N=N/c1cccnc1",
            "C/N=N\\c1cccnc1",
            "CC(=O)NOC(=O)Cc1nc(C(=O)O)cs1",
            "O=C(O)[C@H](O)[C@@H](O)C(=O)OO",
            "CN(C(=O)[O-])C1(B(O)O)CCC1",
            "O[N@]1CC12CC2",
            "CN(C)[N@]1CC1(C)CO",
            "C[P@]1CCCC1C",
            "C[P@]1CCCC12CCCC2",
            "C[As@]1CCCC1C",
            "C[S@+]1CCCC1C",
        ] {
            assert!(
                name_smiles(smiles, &Default::default()).is_err(),
                "{smiles}"
            );
        }
    }

    #[test]
    fn tetrahedral_group_four_ring_hydrides_have_no_lone_pair() {
        for element in ["Si", "Ge", "Sn"] {
            let smiles = format!("C[{element}@H]1CCCC1F");
            assert!(
                name_smiles(&smiles, &Default::default()).is_ok(),
                "{smiles}"
            );
        }
        for smiles in ["C[P@]1CCCC1F", "C[As@]1CCCC1F", "C[S@+]1CCCC1F"] {
            assert_eq!(
                name_smiles(smiles, &Default::default()).unwrap_err(),
                "openclatura: N and ring lone-pair stereocentres are not validated",
                "{smiles}",
            );
        }
    }

    #[test]
    #[ignore = "development parity against the external openclatura oracle TSV"]
    fn parity_file() {
        let path = std::env::var("OPENCLATURA_PARITY").expect("OPENCLATURA_PARITY TSV path");
        let out =
            std::env::var("OPENCLATURA_PARITY_OUT").expect("OPENCLATURA_PARITY_OUT diff path");
        let input = std::fs::read_to_string(path).expect("read parity input");
        let mut output =
            std::io::BufWriter::new(std::fs::File::create(out).expect("create parity diff"));
        writeln!(output, "class\tsmiles\texpected\tactual\tstatus\treason").unwrap();
        let (mut exact, mut different, mut declined) = (0, 0, 0);
        for (index, line) in input.lines().filter(|line| !line.is_empty()).enumerate() {
            let (smiles, expected) = line.split_once('\t').expect("SMILES and expected name");
            if std::env::var_os("OPENCLATURA_PARITY_VERBOSE").is_some() {
                eprintln!("parity row {index}: {smiles}");
            }
            match name_smiles(smiles, &Default::default()) {
                Ok(actual) if actual.text == expected => exact += 1,
                Ok(actual) => {
                    different += 1;
                    let wire = serde_json::to_value(actual.status).unwrap();
                    let status = wire["status"].as_str().unwrap();
                    let reason = wire["reason"].as_str().unwrap_or("");
                    writeln!(
                        output,
                        "different\t{smiles}\t{expected}\t{}\t{status}\t{reason}",
                        actual.text
                    )
                    .unwrap();
                }
                Err(reason) => {
                    declined += 1;
                    writeln!(output, "declined\t{smiles}\t{expected}\t{reason}").unwrap();
                }
            }
            output.flush().unwrap();
        }
        output.flush().unwrap();
        println!(
            "openclatura parity: exact={exact} different={different} declined={declined} total={}",
            exact + different + declined
        );
    }
}
