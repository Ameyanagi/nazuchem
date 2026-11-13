// SPDX-License-Identifier: MIT AND (MIT OR Apache-2.0)
// Origin: mixed
// Source: openclatura | src/openclatura/fusion/planner.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: ReShiki | src/naming/openclatura/fusion/planner.rs | 5f559fd482e237307a589dfcb7ea333cee0a6e13
// Source: openclatura | tests/roundtrip/test_bridges.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/integration/test_systematic_fusion.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: Nazuchem | crates/nazuchem-s2n/src/fusion/planner.rs | issue #95
// Copyright (c) 2026 lamalab-org
// Copyright (c) 2026 Ameyanagi and ReShiki contributors
// Changes: Openclatura audited fusion planner with separate local benzo-parent planning, spiro entry and localized-pi refusal proof.
// Local contributions: MIT OR Apache-2.0; see NOTICE and PROVENANCE-nazuchem.toml.

//! Audited systematic fusion planning in upstream candidate and rule order.
use super::super::{
    locants::retained_locant_sort_key as locant_key,
    molecule::{AtomSet, Molecule},
    polycycle_topology::{bicyclo_proof, edges_within_atoms},
    retained_fused_templates::is_ium_nitrogen,
};
use super::{
    audit::audit_fusion_plan,
    charges::{fusion_charge_operations, fusion_component_charge_parent},
    citation_numbering::numbered_parent_graphs_agree,
    citation_pi::citation_pi_dead_end,
    citation_pi_donor::charged_donor_citation_pi_dead_end,
    config::fusion_nomenclature_config,
    descriptor::{DescriptorOptions, iter_fusion_name_asts, render_fusion_name},
    faces::{BoundedFaceModel, cached_bounded_face_model, typed_face_model},
    indicated_hydrogen::{
        component_carbon_h_relocation_scope, component_h_locants, component_parent_atoms,
        component_parent_graph, intrinsic_carbon_candidate_atoms, intrinsic_carbon_parent_model,
        intrinsic_parent_lone_pair_sites, parent_hydrogen_candidates,
    },
    layout::{
        can_use_component_entry_layout, component_entry_layouts, preferred_intrinsic_layouts,
    },
    mancude::{
        ParentDerivativeOptions, ParentDerivativeState, all_single,
        cited_saturated_nitrogen_hydrogen_sites, has_complete_saturated_hydrogenation,
        indicated_hydrogen_parent_bond_model, parent_degree, parent_derivative_state,
        protonated_hydrogenation_is_complete, saturated_nitrogen_hydrogen_sites,
        single_site_parent_model, spiro_carbon_sites,
    },
    model::*,
    numbering::{
        CompletedNumberingSelection, completed_system_numbering_selection,
        observed_parent_matches_bond_model, parent_bond_model,
    },
    registry::{FusionComponentRegistry, fusion_component_registry},
    rules::{fusion_ring_size_gate, fusion_supported},
    valence::fusion_lambda_descriptors,
};
use std::collections::BTreeMap;

type PlanResult = Result<Option<FusionParentPlan>, String>;
type Specs = BTreeMap<usize, FusionComponentSpec>;

/// Only a confirmed plan enters naming. An audited abstention preserves the
/// legacy route, while an unported branch cannot silently choose a different name.
// Provenance: fn plan_fusion_parent | adapted
// Source: openclatura | src/openclatura/fusion/planner.py::plan_fusion_parent
pub fn plan_fusion_parent(mol: &Molecule, atoms: &AtomSet) -> PlanResult {
    if let Some(cached) = mol.fusion_plan_cache.borrow().get(atoms) {
        return cached.clone();
    }
    let result = plan_uncached(mol, atoms, false, true);
    mol.fusion_plan_cache
        .borrow_mut()
        .insert(atoms.clone(), result.clone());
    result
}

/// Component-local indicated H cannot discard a valid lower spiro junction.
/// This context is kept outside the ordinary isolated-parent plan cache.
// Provenance: fn plan_fusion_parent_for_spiro | original
// Source: ReShiki | src/naming/openclatura/fusion/planner.rs::plan_fusion_parent_for_spiro
pub(in super::super) fn plan_fusion_parent_for_spiro(
    mol: &Molecule,
    atoms: &AtomSet,
) -> PlanResult {
    plan_uncached(mol, atoms, true, true)
}

// Provenance: fn rejected | adapted
fn rejected(error: String) -> PlanResult {
    if error.starts_with("openclatura:") && error.contains("not ported") {
        Err(error)
    } else {
        Ok(None)
    }
}

// Provenance: fn plan_uncached | adapted
// Source: openclatura | src/openclatura/fusion/planner.py::_plan_uncached
fn plan_uncached(
    mol: &Molecule,
    atoms: &AtomSet,
    spiro_component: bool,
    neutral_parent: bool,
) -> PlanResult {
    let config = fusion_nomenclature_config();
    let support = &config.rules.support;
    if atoms.len() < 6 || atoms.iter().any(|a| !mol.atoms.contains_key(a)) {
        return Ok(None);
    }
    // P-58.2.4, P-72.2.2.1, P-73.1.1.2: substituents, '-ide' and '-ium' are
    // added to the neutral parent hydride, so an N-substituted, anionic or
    // quaternary ring N keeps the indicated or hydro hydrogen of its N-H
    // parent. Plan that parent; the plan proved on the input is kept when it
    // cites the same hydrogen. A bridge on the N is cut with the substituents,
    // so a bridged parent uses the N-H parent only for a pheno...ine
    // component (P-25.2.2.3), which the input structure does not keep.
    let n_h_sites: Vec<_> = atoms
        .iter()
        .copied()
        .filter(|&id| {
            let atom = &mol.atoms[&id];
            let exocyclic = mol
                .get_neighbors(id)
                .iter()
                .filter(|n| !atoms.contains(n))
                .count();
            atom.symbol == "N"
                && atom.total_h_count == 0
                && all_single(mol, id)
                && mol.degree(id) == exocyclic + 2
                && exocyclic as i32 == 1 + atom.charge
        })
        .collect();
    // An ion whose ring also carries a suffix group cites it as a prefix
    // (P-73.1.2.1, '2,2-dimethyl-1,3-dioxo-...-2-ium'); that is not composed here.
    let ring_suffix = atoms.iter().any(|&a| {
        mol.get_neighbors(a)
            .iter()
            .any(|n| !atoms.contains(n) && mol.bond_order(a, *n) > 1)
    });
    if neutral_parent
        && !n_h_sites.is_empty()
        && !(ring_suffix && n_h_sites.iter().any(|id| mol.atoms[id].charge != 0))
    {
        let mut kept: AtomSet = mol.atoms.keys().copied().collect();
        for &id in &n_h_sites {
            kept.retain(|a| atoms.contains(a) || !mol.get_neighbors(id).contains(a));
        }
        let mut parent = mol.subgraph(&kept);
        for &id in &n_h_sites {
            if let Some(atom) = parent.atoms.get_mut(&id) {
                atom.charge = 0;
                atom.total_h_count = 1;
                atom.explicit_h_count = 1;
            }
        }
        // An ion at a hydro N is named with added hydrogen instead
        // ('pyridin-1(2H)-ide', P-72.2.2.1.1), which is not composed here: an
        // ionic N must carry the parent's only indicated hydrogen.
        let charged: Vec<_> = n_h_sites
            .iter()
            .filter(|id| mol.atoms[id].charge != 0)
            .collect();
        if let Ok(Some(mut plan)) = plan_uncached(&parent, atoms, spiro_component, true)
            && (!bridged(mol, atoms)
                || plan.ast.component_occurrences.iter().any(|m| {
                    fusion_component_registry()
                        .spec_for_match(m)
                        .is_ok_and(|s| s.parent_hydride_only)
                }))
            && (charged.is_empty()
                || plan.indicated_hydrogens.len() == 1
                    && charged.iter().all(|id| {
                        plan.numbering
                            .abstract_atom_to_locant
                            .iter()
                            .any(|(a, l)| a == *id && plan.indicated_hydrogens.contains(l))
                    }))
        {
            if let Ok(Some(own)) = plan_uncached(mol, atoms, spiro_component, false)
                && own.rendered_base_name == plan.rendered_base_name
                && (own.derivative_state.hydro_operations.iter())
                    .map(|operation| &operation.locants)
                    .eq((plan.derivative_state.hydro_operations.iter())
                        .map(|operation| &operation.locants))
            {
                return Ok(Some(own));
            }
            plan.charge_operations =
                fusion_charge_operations(mol, &plan.abstract_parent_graph, &plan.numbering)?;
            for variant in &mut plan.numbering_variants {
                variant.charge_operations = fusion_charge_operations(
                    mol,
                    &variant.abstract_parent_graph,
                    &variant.numbering,
                )?;
            }
            return Ok(Some(plan));
        }
    }
    // P-73.1.1.2 / P-14.4(b): a cation formed by adding a hydron to a
    // saturated ring N is named from its neutral parent hydride plus -ium, so
    // plan its indicated and hydro hydrogen on the unsubstituted neutral
    // parent (ring N substituents become H). A quaternary ring N is the same
    // cation with its hydrogen substituted; one with four ring bonds has no
    // neutral parent. Ring atoms with exocyclic multiple bonds (suffix
    // composition) stay on the charged path below.
    let protonated: Vec<_> = atoms
        .iter()
        .copied()
        .filter(|&id| is_ium_nitrogen(mol, id) && parent_degree(mol, atoms, id) < 4)
        .collect();
    if !protonated.is_empty()
        && !spiro_component
        && atoms.iter().all(|&a| {
            mol.get_neighbors(a)
                .iter()
                .all(|n| atoms.contains(n) || mol.bond_order(a, *n) == 1)
        })
    {
        let mut kept: AtomSet = mol.atoms.keys().copied().collect();
        for &id in &protonated {
            kept.retain(|a| atoms.contains(a) || !mol.get_neighbors(id).contains(a));
        }
        let mut neutral = mol.subgraph(&kept);
        for &id in &protonated {
            let h = 3 - neutral.degree(id);
            if let Some(atom) = neutral.atoms.get_mut(&id) {
                atom.charge = 0;
                atom.total_h_count = h;
                atom.explicit_h_count = h;
            }
        }
        if let Ok(Some(mut plan)) = plan_uncached(&neutral, atoms, false, true) {
            let locants: BTreeMap<_, _> = plan
                .numbering
                .abstract_atom_to_locant
                .iter()
                .cloned()
                .collect();
            let lowest_saturated = atoms
                .iter()
                .filter(|&&a| all_single(mol, a))
                .map(|a| locant_key(&locants[a]))
                .min();
            // P-31.2.2: indicated hydrogen takes lower locants than hydro
            // prefixes; a saturated parent's plan already cites the lowest
            // positions its mancude parent permits (P-58.2.1.2), which a
            // fusion N can move off the lowest locant. P-25.2.1: a retained
            // neutral parent (quinolizine, naphthyridine) outranks the
            // planner's systematic fusion name.
            let lowest_indicated = atoms.iter().all(|&a| all_single(&neutral, a))
                || plan
                    .indicated_hydrogens
                    .iter()
                    .all(|l| lowest_saturated.as_ref() == Some(&locant_key(l)));
            let retained_parent = !neutral.use_benzoheterocycle_names()
                || super::super::rules::retained::get_retained_ring(
                    &neutral,
                    &atoms.iter().copied().collect::<Vec<_>>(),
                )
                .is_none_or(|retained| {
                    plan.rendered_base_name.ends_with(
                        retained
                            .name
                            .split_once("H-")
                            .map_or(retained.name.as_str(), |(_, name)| name),
                    )
                });
            if lowest_indicated
                && retained_parent
                && let Ok(operations) =
                    fusion_charge_operations(mol, &plan.abstract_parent_graph, &plan.numbering)
            {
                plan.charge_operations = operations;
                for variant in &mut plan.numbering_variants {
                    variant.charge_operations = fusion_charge_operations(
                        mol,
                        &variant.abstract_parent_graph,
                        &variant.numbering,
                    )?;
                }
                return Ok(Some(plan));
            }
        }
    }
    if (!support.charged_parents && atoms.iter().any(|a| mol.atoms[a].charge != 0))
        || atoms
            .iter()
            .any(|a| !fusion_supported(&mol.atoms[a].symbol))
        || (!support.nonstandard_valence && !standard_valence_parent(mol, atoms))
    {
        return Ok(None);
    }
    let bounded = match cached_bounded_face_model(mol, atoms) {
        Ok(Some(bounded)) => bounded,
        Ok(None) => return Ok(None),
        Err(error) => return rejected(error),
    };
    if bounded.cycle_rank < 2 {
        return Ok(None);
    }
    if !bounded.interior_atoms().is_empty()
        && (bicyclo_proof(atoms, &edges_within_atoms(mol, atoms)).is_some()
            || !support.interior_atoms)
    {
        return Ok(None);
    }
    if !fusion_ring_size_gate(
        &bounded
            .faces
            .iter()
            .map(|f| f.atoms.len())
            .collect::<Vec<_>>(),
    ) {
        return Ok(None);
    }
    let registry = fusion_component_registry();
    let matching_parent = match fusion_component_charge_parent(mol, atoms) {
        Ok(parent) => parent,
        Err(error) => return rejected(error),
    };
    let face_model = typed_face_model(mol, &bounded);
    let matches = match registry.match_faces(mol, &matching_parent, &face_model.faces, None) {
        Ok(matches) => matches,
        Err(error) => return rejected(error),
    };
    if mol.use_benzoheterocycle_names()
        && face_model.faces.len() == 2
        && let Some(plan) = plan_benzoheterocycle(mol, atoms, &matches, registry, &face_model)?
    {
        return Ok(Some(plan));
    }
    let mut intrinsic_layouts = None;
    if !can_use_component_entry_layout(&face_model) {
        match preferred_intrinsic_layouts(&face_model, config.search.layout_states) {
            Ok(layouts) if !layouts.is_empty() => intrinsic_layouts = Some(layouts),
            Ok(_) => return Ok(None),
            Err(error) => return rejected(error),
        }
    }
    let options = DescriptorOptions::default();
    let candidates = match iter_fusion_name_asts(mol, &matches, registry, &options) {
        Ok(candidates) => candidates,
        Err(error) => return rejected(error),
    };
    let mut numbering_cache = BTreeMap::new();
    for candidate in candidates.take(config.search.maximum_name_candidates) {
        let ast = match candidate {
            Ok(ast) => ast,
            Err(error) => return rejected(error),
        };
        match plan_numbered_candidate(
            mol,
            atoms,
            &ast,
            registry,
            &bounded,
            &face_model,
            &mut intrinsic_layouts,
            &mut numbering_cache,
            spiro_component,
        ) {
            Ok(Some(plan)) => return Ok(Some(plan)),
            Ok(None) => (),
            Err(error) => {
                rejected(error)?;
            }
        }
    }
    Ok(None)
}

#[allow(clippy::too_many_arguments)]
// Provenance: fn plan_numbered_candidate | adapted
// Source: openclatura | src/openclatura/fusion/planner.py::_plan_numbered_candidate
pub(super) fn plan_numbered_candidate(
    mol: &Molecule,
    atoms: &AtomSet,
    ast: &FusionNameAst,
    registry: &FusionComponentRegistry,
    bounded: &BoundedFaceModel,
    face_model: &FaceModel,
    intrinsic_layouts: &mut Option<Vec<FusedLayout>>,
    numbering_cache: &mut BTreeMap<bool, CompletedNumberingSelection>,
    spiro_component: bool,
) -> PlanResult {
    let specs = component_specs(ast, registry)?;
    let directed = component_entry_layouts(face_model, ast, &specs)?;
    let component_directed = directed.is_some();
    let layouts = if let Some(layouts) = &directed {
        layouts.as_slice()
    } else {
        if intrinsic_layouts.is_none() {
            *intrinsic_layouts = Some(preferred_intrinsic_layouts(
                face_model,
                fusion_nomenclature_config().search.layout_states,
            )?);
        }
        intrinsic_layouts.as_deref().unwrap_or_default()
    };
    if layouts.is_empty() {
        return Ok(None);
    }
    let mut changed_component_parent = false;
    for spec in specs.values() {
        changed_component_parent |= component_parent_atoms(spec)? != spec.atoms();
    }
    let carbon_h_candidates = intrinsic_carbon_candidate_atoms(ast, &specs, mol)?;
    let suffix_carbon_h: AtomSet = carbon_h_candidates
        .iter()
        .copied()
        .filter(|&a| {
            mol.atoms[&a].symbol == "C"
                && mol.atoms[&a].charge == 0
                && mol.get_neighbors(a).iter().any(|&n| {
                    !atoms.contains(&n)
                        && mol.atoms[&n].symbol == "O"
                        && mol.atoms[&n].charge == 0
                        && mol.get_bond(a, n).is_some_and(|b| b.order == 2)
                })
        })
        .collect();
    let prove_carbon_h = !carbon_h_candidates.is_empty()
        || (changed_component_parent
            && specs
                .values()
                .all(|s| s.atoms().iter().all(|a| a.symbol == "C" && a.charge == 0)));
    let defer_indicated_hydrogen = prove_carbon_h || spiro_component;
    let selection = if !component_directed {
        if let Some(cached) = numbering_cache.get(&defer_indicated_hydrogen) {
            cached.clone()
        } else {
            let selection = completed_system_numbering_selection(
                mol,
                bounded,
                Some(face_model),
                layouts,
                defer_indicated_hydrogen,
            )?;
            numbering_cache.insert(defer_indicated_hydrogen, selection.clone());
            selection
        }
    } else {
        completed_system_numbering_selection(
            mol,
            bounded,
            Some(face_model),
            layouts,
            defer_indicated_hydrogen,
        )?
    };
    if selection.accepted.is_empty() {
        return Ok(None);
    }
    if selection.accepted.len() > 1
        && !numbered_parent_graphs_agree(
            mol,
            &selection
                .accepted
                .iter()
                .map(|n| n.string_map())
                .collect::<Vec<_>>(),
        )
    {
        return Ok(None);
    }
    // Keep the intrinsic FR-5 proof when every preferred numbering is also
    // witnessed by the parser's component-entry perimeter.
    let independent = if ast.joins.iter().any(|join| join.order > 1)
        && layouts
            .iter()
            .all(|layout| layout.kind == LayoutKind::OpsinEntry)
    {
        preferred_intrinsic_layouts(
            face_model,
            fusion_nomenclature_config().search.layout_states,
        )
        .ok()
        .and_then(|layouts| {
            let intrinsic = completed_system_numbering_selection(
                mol,
                bounded,
                Some(face_model),
                &layouts,
                defer_indicated_hydrogen,
            )
            .ok()?;
            let directed_maps: Vec<_> = selection.accepted.iter().map(|n| n.string_map()).collect();
            (!intrinsic.accepted.is_empty()
                && intrinsic
                    .accepted
                    .iter()
                    .all(|n| directed_maps.contains(&n.string_map())))
            .then_some((layouts, intrinsic))
        })
    } else {
        None
    };
    let (layouts, selection) = independent
        .as_ref()
        .map_or((layouts, &selection), |(layouts, selection)| {
            (layouts.as_slice(), selection)
        });
    let graph = component_parent_graph(
        ast,
        &specs,
        component_carbon_h_relocation_scope(ast, &specs),
    )?;
    let intrinsic_sites = intrinsic_parent_lone_pair_sites(mol, &graph);
    let initial_model = if intrinsic_sites.is_empty() {
        parent_bond_model(&graph)?
    } else {
        indicated_hydrogen_parent_bond_model(&graph, &intrinsic_sites)?
    };
    let mut alternatives = Vec::new();
    for candidate in &selection.accepted {
        let Some(index) = candidate.layout_index else {
            return Ok(None);
        };
        let layout = &layouts[index];
        let proof = FusionNumberingProof {
            selected_face_model: face_model.clone(),
            selected_layout: layout.clone(),
            orientation_score: ComparisonKey::Tuple(vec![
                layout.orientation_score.clone().into(),
                candidate.score.clone(),
            ]),
            abstract_atom_to_locant: candidate.atom_to_locant.clone(),
            input_locant_maps: vec![candidate.atom_to_locant.clone()],
            rejected_numberings: selection.rejected.clone(),
        };
        match complete_fusion_plan(
            mol,
            atoms,
            ast,
            registry,
            &specs,
            proof,
            &graph,
            &initial_model,
        ) {
            Ok(Some(plan))
                if citation_pi_dead_end(&plan).is_none()
                    && charged_donor_citation_pi_dead_end(mol, &plan).is_none() =>
            {
                alternatives.push(plan)
            }
            Err(error) => {
                rejected(error)?;
            }
            _ => (),
        }
    }
    if alternatives.is_empty() {
        return Ok(None);
    }
    let mut suffix_composition_proved = false;
    if !suffix_carbon_h.is_empty() {
        for plan in &alternatives {
            let Some(locants) = plan.numbering.string_input_locant_maps().into_iter().next() else {
                continue;
            };
            let carbon_h = |cited: &[String]| {
                locants
                    .iter()
                    .any(|(a, l)| mol.atoms[a].symbol == "C" && cited.contains(l))
            };
            if carbon_h(&plan.indicated_hydrogens)
                && let Some((indicated, _, _)) =
                    super::super::additive::systematic_suffix_hydrogen_state(
                        mol,
                        plan,
                        &suffix_carbon_h,
                        &locants,
                    )?
                && carbon_h(&indicated)
            {
                suffix_composition_proved = true;
                break;
            }
        }
    }
    // A ketone suffix can consume intrinsic carbon H at a different site.
    // Keep tied skeletal numberings for suffix composition in the assembly
    // pipeline before comparing indicated-H and suffix locants (P-14.4).
    if prove_carbon_h && !suffix_composition_proved && !spiro_component {
        let best = alternatives
            .iter()
            .map(|p| {
                p.indicated_hydrogens
                    .iter()
                    .map(|s| locant_key(s))
                    .collect::<Vec<_>>()
            })
            .min()
            .unwrap_or_default();
        alternatives.retain(|p| {
            p.indicated_hydrogens
                .iter()
                .map(|s| locant_key(s))
                .collect::<Vec<_>>()
                == best
        });
    }
    let mut plan = alternatives[0].clone();
    plan.numbering.input_locant_maps = alternatives
        .iter()
        .map(|p| p.numbering.input_locant_maps[0].clone())
        .collect();
    plan.numbering_variants = alternatives;
    Ok(Some(plan))
}

// Provenance: fn plan_benzoheterocycle | adapted
// Source: openclatura | src/openclatura/fusion/planner.py::_complete_fusion_plan
fn plan_benzoheterocycle(
    mol: &Molecule,
    atoms: &AtomSet,
    matches: &[FusionComponentMatch],
    registry: &FusionComponentRegistry,
    face_model: &FaceModel,
) -> PlanResult {
    let mut variants = Vec::new();
    for benzo in matches.iter().filter(|m| m.spec_key == "benzene") {
        for hetero in matches.iter().filter(|m| m.spec_key != "benzene") {
            for matched in registry.benzoheterocycle_parent_matches(benzo, hetero)? {
                if variants
                    .iter()
                    .any(|p: &FusionParentPlan| p.ast.component_occurrences == [matched.clone()])
                {
                    continue;
                }
                // P-25.2.2.4 fixes numbering on the heterocyclic perimeter;
                // the generated component supplies both admissible directions.
                let locants: Vec<_> = matched
                    .local_to_input_atom
                    .iter()
                    .map(|(l, a)| (*a, l.clone()))
                    .collect();
                let ast = FusionNameAst {
                    plan_kind: "benzoheterocycle".into(),
                    parent_occurrences: vec![matched.occurrence_id],
                    citation_tree: Some(FusionCitationNode {
                        occurrence_id: matched.occurrence_id,
                        children: vec![],
                    }),
                    component_occurrences: vec![matched],
                    ..Default::default()
                };
                let specs = component_specs(&ast, registry)?;
                let graph = component_parent_graph(
                    &ast,
                    &specs,
                    component_carbon_h_relocation_scope(&ast, &specs),
                )?;
                let intrinsic_sites = intrinsic_parent_lone_pair_sites(mol, &graph);
                let model = if intrinsic_sites.is_empty() {
                    parent_bond_model(&graph)?
                } else {
                    indicated_hydrogen_parent_bond_model(&graph, &intrinsic_sites)?
                };
                let numbering = FusionNumberingProof {
                    selected_face_model: face_model.clone(),
                    abstract_atom_to_locant: locants.clone(),
                    input_locant_maps: vec![locants],
                    ..Default::default()
                };
                if let Some(plan) = complete_fusion_plan(
                    mol, atoms, &ast, registry, &specs, numbering, &graph, &model,
                )? {
                    variants.push(plan);
                }
            }
        }
    }
    let Some(mut plan) = variants.first().cloned() else {
        return Ok(None);
    };
    plan.numbering.input_locant_maps = variants
        .iter()
        .map(|p| p.numbering.input_locant_maps[0].clone())
        .collect();
    plan.numbering_variants = variants;
    Ok(Some(plan))
}

// Provenance: fn component_specs | adapted
// Source: openclatura | src/openclatura/fusion/planner.py::inline in _abstract_graph
fn component_specs(
    ast: &FusionNameAst,
    registry: &FusionComponentRegistry,
) -> Result<Specs, String> {
    ast.component_occurrences
        .iter()
        .map(|m| Ok((m.occurrence_id, registry.spec_for_match(m)?)))
        .collect()
}

#[allow(clippy::too_many_arguments)]
// Provenance: fn complete_fusion_plan | adapted
// Source: openclatura | src/openclatura/fusion/planner.py::_complete_fusion_plan
fn complete_fusion_plan(
    mol: &Molecule,
    atoms: &AtomSet,
    ast: &FusionNameAst,
    registry: &FusionComponentRegistry,
    specs: &Specs,
    numbering: FusionNumberingProof,
    graph: &FusionGraph,
    initial_model: &ParentBondModel,
) -> PlanResult {
    let charge_operations = fusion_charge_operations(mol, graph, &numbering)?;
    let intrinsic_n_h = intrinsic_parent_lone_pair_sites(mol, graph);
    let input_locants: BTreeMap<_, _> = numbering.input_locant_maps[0].iter().cloned().collect();
    let initial_citations = cited_indicated_hydrogens(mol, ast, specs, &numbering, initial_model)?;
    let cited_n_h = input_locants
        .iter()
        .filter(|(a, l)| {
            initial_citations.contains(l) && mol.atoms[a].symbol == "N" && mol.atoms[a].charge == 0
        })
        .map(|(&a, _)| a)
        .collect();
    let carbon_candidates = intrinsic_carbon_candidate_atoms(ast, specs, mol)?;
    let (mut bond_model, mut intrinsic_carbon_h) = intrinsic_carbon_parent_model(
        mol,
        graph,
        initial_model,
        &input_locants,
        &carbon_candidates,
        &intrinsic_n_h,
        &cited_n_h,
    )?;
    if ast.plan_kind == "benzoheterocycle" && !cited_n_h.is_subset(&intrinsic_n_h) {
        // P-31.1.4.2.4 / P-58.2: indicated hydrogen precedes hydro locants,
        // including parent N-H replaced by a substituent. Either may instead
        // be supplied by hydrogenation of C=N.
        let (candidate, carbon_h) = intrinsic_carbon_parent_model(
            mol,
            graph,
            initial_model,
            &input_locants,
            &carbon_candidates,
            &intrinsic_n_h,
            &AtomSet::new(),
        )?;
        let hydrogen_key = |sites: AtomSet| {
            let mut locants: Vec<_> = sites
                .iter()
                .map(|a| locant_key(&input_locants[a]))
                .collect();
            locants.sort();
            locants
        };
        if !carbon_h.is_empty()
            && hydrogen_key(carbon_h.union(&intrinsic_n_h).copied().collect())
                < hydrogen_key(intrinsic_carbon_h.union(&cited_n_h).copied().collect())
        {
            bond_model = candidate;
            intrinsic_carbon_h = carbon_h;
        }
    }
    let mut indicated_h = cited_indicated_hydrogens(mol, ast, specs, &numbering, &bond_model)?;
    indicated_h.extend(intrinsic_carbon_h.iter().map(|a| input_locants[a].clone()));
    indicated_h.sort_by_key(|l| locant_key(l));
    indicated_h.dedup();
    let mut indicated_atoms: AtomSet = input_locants
        .iter()
        .filter(|(_, l)| indicated_h.contains(l))
        .map(|(&a, _)| a)
        .collect();
    let Some(mut derivative_state) = parent_derivative_state(
        mol,
        atoms,
        &bond_model,
        &input_locants,
        &ParentDerivativeOptions {
            indicated_hydrogen_atom_ids: Some(&indicated_atoms),
            ..Default::default()
        },
    )?
    else {
        return Ok(None);
    };
    let hydrogenated = derivative_state.bond_delta.hydrogenated_atom_ids();
    if !indicated_atoms.is_empty()
        && atoms.len().is_multiple_of(2)
        && specs
            .values()
            .all(|s| s.template.default_indicated_h.is_empty())
        && (atoms
            .difference(&indicated_atoms)
            .filter(|a| !hydrogenated.contains(a))
            .count()
            > 1
            || hydrogenated.iter().any(|a| mol.atoms[a].symbol == "N"))
        && !saturated_nitrogen_hydrogen_sites(
            mol,
            atoms,
            &atoms
                .iter()
                .filter(|a| mol.atoms[a].symbol == "N")
                .copied()
                .collect(),
        )
        .is_empty()
        && let Some(complete_state) = parent_derivative_state(
            mol,
            atoms,
            &bond_model,
            &input_locants,
            &ParentDerivativeOptions::default(),
        )?
        && has_complete_saturated_hydrogenation(mol, atoms, &complete_state)
    {
        derivative_state = complete_state;
        indicated_h.clear();
    }
    // P-14.7.1, P-58.2.1.2: a neutral ring atom whose ring bonds use up its
    // standard valence (a fusion N) takes no double bond, so such atoms can
    // split the mancude parent into parts that each leave one position
    // unpaired ('tetrahydro-1H,5H-pyrazolo[1,2-a]pyrazole'). A saturated
    // system cites exactly the indicated hydrogen that parent needs, at the
    // lowest nonfusion locants; hydro prefixes pair up the rest (P-31.2.2).
    let fixed: AtomSet = atoms
        .iter()
        .copied()
        .filter(|&a| {
            mol.atoms[&a].charge == 0
                && mol.atoms[&a]
                    .element()
                    .is_some_and(|e| e.standard_valence <= parent_degree(mol, atoms, a))
        })
        .collect();
    if atoms
        .iter()
        .all(|&a| mol.atoms[&a].charge == 0 && !mol.atoms[&a].is_aromatic && all_single(mol, a))
        && spiro_carbon_sites(mol, atoms).is_empty()
    {
        let mancude = single_site_parent_model(initial_model, &fixed)?;
        if atoms.len() - fixed.len() - 2 * mancude.maximum_non_cumulative_double_bonds
            != indicated_h.len()
        {
            let key = |a: &usize| locant_key(&input_locants[a]);
            let Some(sites) = mancude
                .allowed_kekule_assignments
                .iter()
                .map(|assignment| {
                    let paired = assignment.paired_atoms();
                    let mut sites: Vec<_> = atoms
                        .iter()
                        .copied()
                        .filter(|a| !fixed.contains(a) && !paired.contains(a))
                        .collect();
                    sites.sort_by_key(key);
                    sites
                })
                .min_by_key(|sites| {
                    (
                        sites
                            .iter()
                            .filter(|&&a| parent_degree(mol, atoms, a) > 2)
                            .count(),
                        sites.iter().map(key).collect::<Vec<_>>(),
                    )
                })
            else {
                return Ok(None);
            };
            indicated_h = sites.iter().map(|a| input_locants[a].clone()).collect();
            indicated_atoms = sites.into_iter().collect();
            bond_model = single_site_parent_model(&mancude, &indicated_atoms)?;
            let Some(state) = parent_derivative_state(
                mol,
                atoms,
                &bond_model,
                &input_locants,
                &ParentDerivativeOptions {
                    indicated_hydrogen_atom_ids: Some(&indicated_atoms),
                    ..Default::default()
                },
            )?
            else {
                return Ok(None);
            };
            derivative_state = state;
        }
    }
    // Upstream divergence: removing a reduced edge incident to indicated N-H
    // can lose the carbon endpoint's hydrogen. Require the cited H/hydro sites
    // to preserve the ring pi count and the pi degree of each nonaromatic carbon.
    if !protonated_hydrogenation_is_complete(
        mol,
        atoms,
        &bond_model,
        &derivative_state,
        &indicated_atoms,
    )? {
        return Err("openclatura: protonated fused hydrogenation state not ported".into());
    }
    if fusion_nomenclature_config()
        .rules
        .support
        .maximum_indicated_hydrogens
        .is_some_and(|n| indicated_h.len() > n)
    {
        return Ok(None);
    }
    let rendered_core = render_fusion_name(ast, registry)?;
    let lambda_descriptors = fusion_lambda_descriptors(mol, atoms, &input_locants)?;
    let mut rendered = String::new();
    if !indicated_h.is_empty() {
        rendered.push_str(
            &indicated_h
                .iter()
                .map(|l| format!("{l}H"))
                .collect::<Vec<_>>()
                .join(","),
        );
        rendered.push('-');
    }
    let benzo_name = ast
        .component_occurrences
        .first()
        .filter(|m| m.standalone_benzo_parent)
        .and_then(|m| specs[&m.occurrence_id].benzo_parent_name.as_ref());
    if let Some(name) = benzo_name {
        rendered.push_str(&name.render(&lambda_descriptors));
    } else if !lambda_descriptors.is_empty() {
        rendered.push_str(
            &lambda_descriptors
                .iter()
                .map(|d| d.text())
                .collect::<Vec<_>>()
                .join(","),
        );
        rendered.push('-');
    }
    if benzo_name.is_none() {
        rendered.push_str(&rendered_core);
    }
    let mut plan = FusionParentPlan {
        ast: ast.clone(),
        rendered_base_name: rendered,
        abstract_parent_graph: graph.clone(),
        numbering,
        bond_model,
        indicated_hydrogens: indicated_h,
        derivative_state,
        charge_operations,
        lambda_descriptors,
        ..Default::default()
    };
    plan.audit = audit_fusion_plan(mol, atoms, &plan, registry, &rendered_core)?;
    if !plan.audit.confirmed() {
        return Ok(None);
    }
    // Upstream divergence: aromatic equivalence can hide ring triples and
    // cumulative double bonds during parent-hydride matching. Reject a confirmed
    // fusion citation that would erase them; a planner abstention must still
    // permit the legacy route to cite bridged-parent unsaturation explicitly.
    if mol.bonds.values().any(|bond| {
        atoms.contains(&bond.u) && atoms.contains(&bond.v) && !matches!(bond.order, 1 | 2)
    }) || atoms.iter().any(|&atom| {
        mol.get_neighbors(atom)
            .iter()
            .filter(|&&neighbor| {
                atoms.contains(&neighbor)
                    && mol
                        .get_bond(atom, neighbor)
                        .is_some_and(|bond| bond.order == 2)
            })
            .count()
            > 1
    }) {
        return Err(
            "openclatura: fused-parent cumulative and higher bond orders are not ported".into(),
        );
    }
    if unresolved_localized_pi_cycle(
        mol,
        &plan.bond_model,
        &plan.derivative_state,
        &indicated_atoms,
    ) {
        return Err("openclatura: localized nonaromatic fused pi state not ported".into());
    }
    if !super::layout::parser_numbering_is_supported(&plan.numbering.selected_face_model) {
        return Err("openclatura: nonhexagonal fusion parser numbering not ported".into());
    }
    // P-73.1.1.2: a protonated N-H⁺ cites the hydrogenation of its neutral
    // parent hydride plus -ium. Indicated H at such an N lets the mancude
    // delta drop the hydro citation of its saturated neighbour, naming an
    // impossible parent (3,4-dihydro-2H-isoquinolin-2-ium) or another tautomer.
    // Every sp3 ring carbon must then be cited; correct spelling is not ported.
    let indicated: AtomSet = plan
        .numbering
        .abstract_atom_to_locant
        .iter()
        .filter(|(_, l)| plan.indicated_hydrogens.contains(l))
        .map(|(id, _)| *id)
        .collect();
    let delta = &plan.derivative_state.bond_delta;
    let hydro = delta.hydrogenated_atom_ids();
    let protonated_indicated = indicated.iter().any(|&id| {
        let atom = &mol.atoms[&id];
        atom.symbol == "N"
            && atom.charge == 1
            && atom.total_h_count > 0
            && super::mancude::all_single(mol, id)
    });
    if protonated_indicated && !hydro.is_empty() {
        let cited: AtomSet = indicated
            .iter()
            .chain(&hydro)
            .copied()
            .chain(
                plan.derivative_state
                    .hydro_operations
                    .iter()
                    .chain(&delta.added_hydrogen_operations)
                    .chain(&delta.intrinsic_hydro_operations)
                    .flat_map(|op| op.atom_ids.iter().copied()),
            )
            .collect();
        if atoms.iter().any(|&id| {
            mol.atoms[&id].is_carbon()
                && super::mancude::all_single(mol, id)
                && !cited.contains(&id)
        }) {
            return Err(
                "openclatura: protonated ring N-H hydrogen citation from the neutral parent hydride not ported"
                    .into(),
            );
        }
    }
    // P-31.2.2, P-58.2.1.2: indicated hydrogen takes lower locants than hydro
    // prefixes and prefers nonfusion atoms ('2,3-dihydro-1H-indene', not
    // '1,3-dihydro-2H-indene'; '3a,5-dihydro-4H-indene', not
    // '4,5-dihydro-3aH-indene'). Among the saturated positions it is cited
    // at the lowest ones that still leave a mancude parent.
    let state = &plan.derivative_state;
    if !indicated.is_empty()
        && state.hydro_operations.len() == 1
        && state.external_pi_operations().next().is_none()
        && state.unsaturation_operations.is_empty()
        && delta.added_hydrogen_operations.is_empty()
        && delta.intrinsic_hydro_operations.is_empty()
        && atoms.iter().all(|a| mol.atoms[a].charge == 0)
        && spiro_carbon_sites(mol, atoms).is_empty()
    {
        let key = |a: &usize| {
            (
                parent_degree(mol, atoms, *a) > 2,
                locant_key(&input_locants[a]),
            )
        };
        let mut sites: Vec<_> = indicated
            .iter()
            .chain(&state.hydro_operations[0].atom_ids)
            .copied()
            .collect();
        sites.sort_by_key(key);
        let mut cited = indicated.clone();
        for site in &sites {
            let Some(&old) = cited
                .iter()
                .filter(|a| key(a) > key(site))
                .max_by_key(|a| key(a))
                .filter(|_| !cited.contains(site))
            else {
                continue;
            };
            let mut candidate = cited.clone();
            candidate.remove(&old);
            candidate.insert(*site);
            let single = candidate.union(&fixed).copied().collect();
            if indicated_hydrogen_parent_bond_model(graph, &single).is_ok_and(|model| {
                model.maximum_non_cumulative_double_bonds
                    == plan.bond_model.maximum_non_cumulative_double_bonds
            }) {
                cited = candidate;
            }
        }
        if cited != indicated {
            sites.retain(|a| !cited.contains(a));
            sites.sort_by_key(|a| locant_key(&input_locants[a]));
            let mut locants: Vec<_> = cited.iter().map(|a| input_locants[a].clone()).collect();
            locants.sort_by_key(|l| locant_key(l));
            let mut relocated = plan.clone();
            relocated.rendered_base_name =
                super::super::retained_name_policy::render_retained_hydrogen_state(
                    &plan.rendered_base_name,
                    &plan.indicated_hydrogens,
                    Some(&locants),
                );
            relocated.indicated_hydrogens = locants;
            let model = indicated_hydrogen_parent_bond_model(graph, &cited)?;
            // Prove the relocated hydrogen like the original plan. The bond
            // model cannot place hydro prefixes that move double bonds
            // through an aromatic ring (N-H of '2,3,4,5-tetrahydro-1H-
            // pyrido[4,3-b]indole'); those keep the audited plan and only
            // exchange the cited positions. A bridged parent recomputes this
            // state, so it keeps the plan as it is.
            let proved = parent_derivative_state(
                mol,
                atoms,
                &model,
                &input_locants,
                &ParentDerivativeOptions {
                    indicated_hydrogen_atom_ids: Some(&cited),
                    ..Default::default()
                },
            )?
            .filter(|state| {
                state.unsaturation_operations.is_empty()
                    && state.bond_delta.added_hydrogen_operations.is_empty()
                    && state.bond_delta.intrinsic_hydro_operations.is_empty()
                    && (state.hydro_operations.iter())
                        .map(|operation| &operation.atom_ids)
                        .eq([&sites])
            })
            .map(|state| FusionParentPlan {
                bond_model: model,
                derivative_state: state,
                ..relocated.clone()
            })
            .filter(|proved| {
                audit_fusion_plan(mol, atoms, proved, registry, &rendered_core)
                    .is_ok_and(|audit| audit.confirmed())
            });
            if let Some(proved) = proved {
                plan = proved;
            } else if !bridged(mol, atoms) {
                let operation = &mut relocated.derivative_state.hydro_operations[0];
                operation.locants = sites.iter().map(|a| input_locants[a].clone()).collect();
                operation.atom_ids = sites;
                plan = relocated;
            }
        }
    }
    Ok(Some(plan))
}

/// Whether atoms outside the ring system join two of its atoms.
// Provenance: fn bridged | original
// Source: Nazuchem | crates/nazuchem-s2n/src/fusion/planner.rs::bridged
fn bridged(mol: &Molecule, atoms: &AtomSet) -> bool {
    let mut seen = AtomSet::new();
    for &start in atoms {
        for &outside in mol.get_neighbors(start) {
            let mut joined = AtomSet::new();
            let mut pending = vec![outside];
            while let Some(atom) = pending.pop() {
                if atoms.contains(&atom) {
                    joined.insert(atom);
                } else if seen.insert(atom) {
                    pending.extend(mol.get_neighbors(atom));
                }
            }
            if joined.len() > 1 && !atoms.contains(&outside) {
                return true;
            }
        }
    }
    false
}

// Provenance: fn cited_indicated_hydrogens | adapted
// Source: openclatura | src/openclatura/fusion/planner.py::_cited_indicated_hydrogens
fn cited_indicated_hydrogens(
    mol: &Molecule,
    ast: &FusionNameAst,
    specs: &Specs,
    numbering: &FusionNumberingProof,
    bond_model: &ParentBondModel,
) -> Result<Vec<String>, String> {
    let locants: BTreeMap<_, _> = numbering.abstract_atom_to_locant.iter().cloned().collect();
    let candidates = parent_hydrogen_candidates(mol, ast, &locants, bond_model);
    let mut unpaired = AtomSet::new();
    for assignment in &bond_model.allowed_kekule_assignments {
        let paired = assignment.paired_atoms();
        unpaired.extend(candidates.difference(&paired));
    }
    let mut component_h_roles = AtomSet::new();
    for occurrence in &ast.component_occurrences {
        let map = occurrence.input_atom_by_locant();
        component_h_roles.extend(
            component_h_locants(&specs[&occurrence.occurrence_id], "N")?
                .iter()
                .map(|l| map[l]),
        );
    }
    let mut cited: AtomSet = candidates
        .iter()
        .filter(|a| {
            mol.atoms[a].symbol != "C"
                && (mol.atoms[a].is_aromatic
                    || mol.atoms[a].charge != 0
                    || unpaired.contains(a)
                    || ast.plan_kind != "benzoheterocycle" && component_h_roles.contains(a))
        })
        .copied()
        .collect();
    cited.extend(cited_saturated_nitrogen_hydrogen_sites(
        mol,
        &locants.keys().copied().collect(),
        &candidates
            .iter()
            .filter(|a| mol.atoms[a].symbol == "N")
            .copied()
            .collect(),
        bond_model,
    ));
    let has_attached_carbocycle = ast
        .component_occurrences
        .iter()
        .filter(|m| !ast.parent_occurrences.contains(&m.occurrence_id))
        .any(|m| {
            let spec = &specs[&m.occurrence_id];
            spec.rings.len() == 1
                && spec.usable_as_attached
                && !spec.usable_as_parent
                && spec.atoms().iter().all(|a| a.symbol == "C")
        });
    let has_polycyclic_parent = ast
        .parent_occurrences
        .iter()
        .any(|id| specs[id].rings.len() > 1);
    if has_attached_carbocycle
        && has_polycyclic_parent
        && observed_parent_matches_bond_model(mol, bond_model)
    {
        cited.extend(candidates.iter().filter(|a| mol.atoms[a].symbol == "C"));
    }
    let mut result = cited.iter().map(|a| locants[a].clone()).collect::<Vec<_>>();
    result.sort_by_key(|l| locant_key(l));
    Ok(result)
}

// Provenance: fn standard_valence_parent | adapted
// Source: openclatura | src/openclatura/fusion/planner.py::_standard_valence_parent
fn standard_valence_parent(mol: &Molecule, atoms: &AtomSet) -> bool {
    atoms.iter().all(|&id| {
        let atom = &mol.atoms[&id];
        let h = if atom.total_h_count != 0 {
            atom.total_h_count
        } else {
            atom.explicit_h_count
        };
        atom.element().is_some_and(|e| {
            h + mol
                .get_neighbors(id)
                .iter()
                .filter_map(|&n| mol.get_bond(id, n))
                .map(|b| b.order)
                .sum::<usize>()
                <= e.standard_valence
        })
    })
}

// Upstream divergence: maximum matching proves that the observed pi state
// exists, but not that an unqualified parent citation determines it. An
// alternating 4n cycle with a nonaromatic bond can retain its skeleton and
// hydrogens while reconstructing a different localized double-bond arrangement.
// Aromatic endpoints do not make an explicit inter-ring bond aromatic.
// Provenance: fn unresolved_localized_pi_cycle | original
// Source: ReShiki | src/naming/openclatura/fusion/planner.rs::unresolved_localized_pi_cycle
fn unresolved_localized_pi_cycle(
    mol: &Molecule,
    model: &ParentBondModel,
    state: &ParentDerivativeState,
    indicated: &AtomSet,
) -> bool {
    let model = state.bond_delta.composition_model.as_ref().unwrap_or(model);
    let chosen = &state.bond_delta.assignment;
    if model.allowed_kekule_assignments.len() < 2
        || !model.allowed_kekule_assignments.contains(chosen)
    {
        return false;
    }
    let mut protected = indicated.clone();
    protected.extend(state.bond_delta.hydrogenated_atom_ids());
    for operation in state
        .hydro_operations
        .iter()
        .chain(&state.bond_delta.added_hydrogen_operations)
        .chain(&state.bond_delta.intrinsic_hydro_operations)
    {
        protected.extend(operation.atom_ids.iter().copied());
    }
    protected.extend(
        state
            .unsaturation_operations
            .iter()
            .flat_map(|op| op.atom_ids.iter().copied()),
    );
    protected.extend(state.external_pi_operations().map(|op| op.parent_atom_id));
    if let Some(redistribution) = &state.pi_redistribution {
        protected.extend(&redistribution.hydrogenated_atom_ids);
        for id in redistribution
            .removed_bond_ids
            .union(&redistribution.added_bond_ids)
        {
            if let Some(bond) = mol.bonds.get(id) {
                protected.extend([bond.u, bond.v]);
            }
        }
    }
    let chosen: BTreeMap<_, _> = chosen.orders.iter().copied().collect();
    for alternative in &model.allowed_kekule_assignments {
        let alternative_orders: BTreeMap<_, _> = alternative.orders.iter().copied().collect();
        let mut adjacency: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for &(edge, order) in &alternative.orders {
            if matches!((chosen.get(&edge), order), (Some(1), 2) | (Some(2), 1)) {
                adjacency.entry(edge.0).or_default().push(edge.1);
                adjacency.entry(edge.1).or_default().push(edge.0);
            }
        }
        let mut seen = AtomSet::new();
        for &start in adjacency.keys() {
            if seen.contains(&start) {
                continue;
            }
            let mut component = AtomSet::new();
            let mut pending = vec![start];
            while let Some(atom) = pending.pop() {
                if component.insert(atom) {
                    pending.extend(adjacency[&atom].iter().copied());
                }
            }
            seen.extend(&component);
            if component.len() < 4
                || !component.len().is_multiple_of(4)
                || !component.is_disjoint(&protected)
                || component.iter().any(|id| {
                    let atom = &mol.atoms[id];
                    adjacency[id].len() != 2
                        || atom.charge != 0
                        || !matches!(atom.symbol.as_str(), "C" | "N")
                })
            {
                continue;
            }
            if !component.iter().all(|&atom| {
                adjacency[&atom]
                    .iter()
                    .filter(|&&other| chosen[&normalize_edge(atom, other)] == 2)
                    .count()
                    == 1
            }) {
                continue;
            }
            let localized: Vec<_> = component
                .iter()
                .flat_map(|&atom| adjacency[&atom].iter().map(move |&other| (atom, other)))
                .filter(|&(atom, other)| {
                    atom < other
                        && mol
                            .get_bond(atom, other)
                            .is_some_and(|bond| !bond.is_aromatic)
                })
                .collect();
            if !localized.is_empty()
                && [&chosen, &alternative_orders]
                    .into_iter()
                    .any(|assignment| {
                        localized.iter().all(|&(atom, other)| {
                            mol.get_bond(atom, other)
                                .is_some_and(|bond| bond.order == assignment[&(atom, other)])
                        })
                    })
            {
                return true;
            }
        }
    }
    false
}

#[cfg(test)]
// Provenance: mod tests | original
// Source: ReShiki | src/naming/openclatura/fusion/planner.rs::tests
mod tests {
    use super::super::super::graph_io::read_smiles;
    use super::*;

    #[test]
    // Provenance: fn tests::higher_order_components_preserve_intrinsic_numbering_and_identity | adapted
    fn higher_order_components_preserve_intrinsic_numbering_and_identity() {
        for (smiles, expected) in [
            (
                "C1=NN=C2C1=C1C(=NC=3C=CC=CC13)C=CO2",
                "pyrazolo[4',3':6,7]oxepino[4,5-b]indole",
            ),
            (
                "O1C=CC=C2C1=C1C(=C3N=NC=C3O1)C=CS2",
                "pyrano[2'',3'':6',7']thiepino[4',5':4,5]furo[3,2-c]pyrazole",
            ),
            (
                "S1C=CC=C2C1=C1C(=C3N=NC=C3O1)C=CS2",
                "thiopyrano[2'',3'':6',7']thiepino[4',5':4,5]furo[3,2-c]pyrazole",
            ),
            (
                "[Se]1C=CC=C2C1=C1C(=C3N=NC=C3O1)C=CS2",
                "selenopyrano[2'',3'':6',7']thiepino[4',5':4,5]furo[3,2-c]pyrazole",
            ),
            (
                "[Te]1C=CC=C2C1=C1C(=C3N=NC=C3O1)C=CS2",
                "telluropyrano[2'',3'':6',7']thiepino[4',5':4,5]furo[3,2-c]pyrazole",
            ),
        ] {
            let mol = read_smiles(smiles, &Default::default()).unwrap();
            let atoms = mol.atoms.keys().copied().collect();
            let plan = plan_fusion_parent(&mol, &atoms).unwrap().unwrap();
            assert!(plan.audit.confirmed(), "{:?}", plan.audit);
            assert_eq!(plan.numbering.selected_layout.kind, LayoutKind::Intrinsic);
            let name = crate::name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected);
            assert_eq!(name.status, nazuchem_ir::NameStatus::Preferred);
            let parsed = nazuchem_n2s::Parser::new()
                .unwrap()
                .parse(&name.text, &nazuchem_n2s::ParseOptions::strict());
            assert_eq!(
                parsed.status,
                nazuchem_n2s::Status::Success,
                "{}",
                parsed.message
            );
            assert!(parsed.warnings.is_empty(), "{:?}", parsed.warnings);
            let round_trip = parsed.structure.unwrap().semantic_cxsmiles().unwrap();
            crate::test_support::verify_identity(smiles, &round_trip).unwrap();
        }
    }

    #[test]
    fn protonated_fused_citations_preserve_carbon_hydrogenation() {
        for smiles in [
            "C1=CC2CCCCC2[NH2+]1",
            "C1=CCC2CC[NH2+]CCC2CC1",
            // P-25.2.1: the neutral parent is retained quinolizine, so the
            // planner's pyrido[1,2-a]pyridine name is not used for the cation.
            "C1CC[NH+]2CCCCC2C1",
        ] {
            let mol = read_smiles(smiles, &Default::default()).unwrap();
            let atoms = mol.atoms.keys().copied().collect();
            let error = plan_fusion_parent(&mol, &atoms).unwrap_err();
            assert!(error.contains("protonated"), "{smiles}: {error}");
        }
        // P-73.1.1.2: the cation cites its neutral parent's hydro prefixes,
        // 1,2,3,4,4a,5,6,8a-octahydroisoquinolin-2-ium (not 2H-isoquinoline).
        let mol = read_smiles("C1=CC2C[NH2+]CCC2CC1", &Default::default()).unwrap();
        let atoms = mol.atoms.keys().copied().collect();
        let plan = plan_fusion_parent(&mol, &atoms).unwrap().unwrap();
        assert_eq!(plan.rendered_base_name, "isoquinoline");
        assert!(plan.indicated_hydrogens.is_empty());
        assert_eq!(plan.charge_operations.len(), 1);
        // #95: the neutral parent now cites its indicated hydrogen at the
        // lowest saturated locant (P-31.2.2), so these cations are planned too.
        for (smiles, expected) in [
            (
                "C1=CCC2CC[NH2+]CC2CC1",
                "2,3,4,4a,5,8,9,9a-octahydro-1H-cyclohepta[c]pyridin-2-ium",
            ),
            (
                "C1=CC2CCCC2C[NH2+]1",
                "2,4a,5,6,7,7a-hexahydro-1H-cyclopenta[c]pyridin-2-ium",
            ),
            (
                "C1=CC2CCCCCC2C[NH2+]1",
                "2,4a,5,6,7,8,9,9a-octahydro-1H-cyclohepta[c]pyridin-2-ium",
            ),
        ] {
            let name = crate::name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected);
            let parsed = crate::test_support::parse_strict(&name.text).unwrap();
            crate::test_support::verify_identity(smiles, &parsed).unwrap();
        }
        // #96, P-73.1.1.2: a quaternary ring N is the same cation with its
        // hydrogen substituted, so it is planned on the neutral parent too.
        for (smiles, expected) in [
            (
                "C[N+]1(C)CCC2CCCC2C1",
                "2,2-dimethyloctahydro-1H-cyclopenta[c]pyridin-2-ium",
            ),
            (
                "C[N+]1(C)CCc2sccc2C1",
                "5,5-dimethyl-4,5,6,7-tetrahydrothieno[3,2-c]pyridin-5-ium",
            ),
            // #108, P-58.2.1.2: the fusion N takes no double bond, so the
            // lowest indicated hydrogen of pyrido[1,2-a]pyrazine is 2H.
            (
                "C[N+]1(C)CCN2CCCCC2C1",
                "2,2-dimethyloctahydro-2H-pyrido[1,2-a]pyrazin-2-ium",
            ),
            (
                "C1CCN2CC[NH2+]CC2C1",
                "octahydro-2H-pyrido[1,2-a]pyrazin-2-ium",
            ),
            (
                "C1CCN2CC[NH+](C)CC2C1",
                "2-methyloctahydro-2H-pyrido[1,2-a]pyrazin-2-ium",
            ),
        ] {
            let name = crate::name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected);
            assert_eq!(name.status, nazuchem_ir::NameStatus::Preferred);
            let parsed = crate::test_support::parse_strict(&name.text).unwrap();
            crate::test_support::verify_identity(smiles, &parsed).unwrap();
        }
        // #96 controls, unchanged names: an N⁺ with four ring bonds has no
        // neutral parent hydride, a spiro N⁺ stays with spiro nomenclature
        // (P-24.5.2), and a bridgehead hydron on a von Baeyer parent is
        // spelled 'azonia' (P-73.4), which the parser reads as written.
        for (smiles, expected) in [
            (
                "[N+]123CCC1CCC2CCC3",
                "1-azoniatricyclo[5.3.0.0^{1,4}]decane",
            ),
            (
                "C1CC[N+]2(CC1)CCc1ccccc1C2",
                "spiro[3-azoniabicyclo[4.4.0]deca-1(10),6,8-triene-3,1'-(piperidin-1-ium)]",
            ),
            ("C1C[NH+]2CCC1CC2", "1-azoniabicyclo[2.2.2]octane"),
        ] {
            let name = crate::name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected);
            assert_ne!(
                name.status,
                nazuchem_ir::NameStatus::Unverified(nazuchem_ir::NameReason::FusionNitrogenHydron)
            );
        }
        // #96: this cation's hydro prefixes do not yet match its neutral
        // parent hydride (1-oxodecahydroisoquinolin-2-ium, P-73.1.2.1), so
        // the name keeps its locants and stays unverified.
        let name = crate::name_smiles("O=C1[NH2+]CCC2CCCCC12", &Default::default()).unwrap();
        assert_eq!(
            name.status,
            nazuchem_ir::NameStatus::Unverified(nazuchem_ir::NameReason::HydroPrefixLocants)
        );
        // #96, P-73.1.1.2: a hydron at a fusion N is cited as -ium, but the
        // parser reads that spelling as the λ5-ylium cation (P-73.3.1), so
        // the name stays unverified.
        for (smiles, expected) in [
            ("C1CC2CCC[NH+]2C1", "hexahydro-1H-pyrrolizin-4-ium"),
            ("C1=CC=C2C=CC=C[NH+]2C1", "4H-quinolizin-5-ium"),
        ] {
            let name = crate::name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected);
            assert_eq!(
                name.status,
                nazuchem_ir::NameStatus::Unverified(nazuchem_ir::NameReason::FusionNitrogenHydron)
            );
        }
    }

    #[test]
    // Provenance: fn tests::localized_nonaromatic_pi_alternatives_require_a_citation_witness | adapted
    // Source: ReShiki | src/naming/openclatura/fusion/planner.rs::tests::localized_nonaromatic_pi_alternatives_require_a_citation_witness
    // Source: openclatura | tests/roundtrip/test_bridges.py::test_roundtrip
    // Source: openclatura | tests/integration/test_systematic_fusion.py::test_multiple_indicated_hydrogens_are_graph_derived_and_roundtrip
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn localized_nonaromatic_pi_alternatives_require_a_citation_witness() {
        for text in [
            "C1=CC2=CC3=C(C=CN4CN=CCC34)C2=N1",
            "ClC1=CC=CC2=CC=CC=CC2=C1",
            "c1cn2c3cc4ccnc-4c=3ccn2cn1",
        ] {
            let mol = read_smiles(text, &Default::default()).unwrap();
            let atoms = mol
                .atoms
                .keys()
                .copied()
                .filter(|id| mol.atoms[id].symbol != "Cl")
                .collect();
            assert_eq!(
                plan_fusion_parent(&mol, &atoms).unwrap_err(),
                "openclatura: localized nonaromatic fused pi state not ported"
            );
        }
        for text in ["c1oc2ncccc2c1", "C1OCC2COCC12"] {
            let mol = read_smiles(text, &Default::default()).unwrap();
            let atoms = mol.atoms.keys().copied().collect();
            assert!(
                plan_fusion_parent(&mol, &atoms).unwrap().is_some(),
                "{text}"
            );
        }
    }
}
