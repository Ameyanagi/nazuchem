// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: ReShiki | src/naming/openclatura/heterospiro.rs | 327d32fed63a5b3df15b7b8dc5ab84a7c36e0865
// Copyright (c) 2026 Ameyanagi and ReShiki contributors

//! P-24.5: hydrogen and suffix composition of two different spiro components.
use super::{
    assembly_parts::{AssemblyParts, HydroOperation, SubstituentItem},
    assembly_spiro::{spiro_component, spiro_component_sort_key},
    molecule::{AtomSet, Molecule},
    nomenclature::RULES,
    spirobi::{Component, completed_graph, component, prime_locant, project},
};
use std::collections::BTreeMap;

fn nonstandard_bonding_number(mol: &Molecule, id: usize) -> Option<usize> {
    let atom = &mol.atoms[&id];
    let bonding = atom.total_h_count
        + mol
            .get_neighbors(id)
            .iter()
            .map(|&neighbor| mol.bond_order(id, neighbor))
            .sum::<usize>();
    (atom.charge == 0
        && atom
            .element()
            .is_some_and(|element| bonding > element.standard_valence))
    .then_some(bonding)
}

pub(super) fn has_nonstandard_bonding(mol: &Molecule, atoms: &AtomSet) -> bool {
    atoms
        .iter()
        .any(|&id| nonstandard_bonding_number(mol, id).is_some())
}

fn identity(mol: &Molecule, parts: &AssemblyParts) -> Result<Option<(String, bool)>, String> {
    match component(parts)? {
        Some(Component::Fused(name)) => {
            let indicated: Vec<_> = super::additive::name_indicated_hydrogen_locants(Some(&name))
                .into_iter()
                .collect();
            Ok(Some((
                super::retained_name_policy::render_retained_hydrogen_state(
                    &name,
                    &indicated,
                    Some(&[]),
                ),
                true,
            )))
        }
        Some(Component::VonBaeyer(_)) => Ok(None),
        None if parts.is_simple_ring() => {
            // The side projection can use skeletal replacements internally;
            // recover the retained/Hantzsch–Widman citation from its ring graph.
            if !parts.a_prefixes.is_empty() {
                let mut numbered: Vec<_> = parts.parent_atom_ids_by_locant.iter().collect();
                numbered.sort_by_key(|(l, _)| super::assembly_utils::parse_locant(l));
                let path: Vec<_> = numbered.into_iter().map(|(_, &id)| id).collect();
                return Ok(super::rules::retained::get_retained_ring(mol, &path)
                    .map(|parent| (parent.name, false)));
            }
            let mut parent = parts.clone();
            parent.principal_group = None;
            let (stem, terminal) = super::assembly_parent::parent_stem_and_terminal(&parent)?;
            let tail = super::assembly_parent::format_parent_tail(&parent, &stem, &terminal)?;
            Ok(Some((
                format!("{}{}{}", tail.stem, tail.unsaturation, tail.terminal),
                false,
            )))
        }
        None => Ok(None),
    }
}

pub(super) fn compose(
    mol: &Molecule,
    parts: &AssemblyParts,
) -> Result<Option<(String, String)>, String> {
    let sides: Vec<_> = parts
        .substituents
        .iter()
        .filter_map(|s| s.spiro.as_deref())
        .collect();
    let [side] = sides.as_slice() else {
        return Ok(None);
    };
    let Some(other) = side.side_parts.as_deref() else {
        return Ok(None);
    };
    if side.continuation.is_some() || other.substituents.iter().any(|s| s.spiro.is_some()) {
        return Ok(None);
    }
    let (Some(first), Some(second)) = (identity(mol, parts)?, identity(mol, other)?) else {
        return Ok(None);
    };
    if first.0 == second.0 || !(first.1 || second.1) {
        return Ok(None);
    }
    let shared: Vec<_> = parts
        .parent_atom_ids
        .intersection(&other.parent_atom_ids)
        .copied()
        .collect();
    let [junction] = shared.as_slice() else {
        return Ok(None);
    };
    let junction = *junction;
    if super::spiro_subgraph::heterogeneous_spiro_junction(mol) != Some(junction) {
        return Ok(None);
    }
    // A spiro prefix has one free valence on the component that was rooted.
    if other.is_substituent
        || parts.is_substituent && (parts.attachment_order != 1 || parts.is_acyl_substituent)
    {
        return Ok(None);
    }
    // A prefix cites every characteristic group as a prefix (P-58.2.4).
    let principal_key = parts
        .principal_group
        .as_ref()
        .or(other.principal_group.as_ref())
        .filter(|_| !parts.is_substituent)
        .map(|group| group.key.as_str());
    if principal_key.is_some_and(|key| key != "ketone") {
        return Ok(None);
    }
    for c in [parts, other] {
        if !c.principal_suffix_modifiers.is_empty()
            // P-41: only the senior class is the suffix; a junior group on
            // the other component is cited as a prefix.
            || c.principal_group.as_ref().is_some_and(|g| {
                Some(g.key.as_str()) != principal_key
                    && RULES.functional_groups.prefix_for(&g.key).is_none()
            })
            || !c.stereo_features.is_empty()
            || !c.endo_exo_descriptors.is_empty()
            || !c.relative_stereo_prefixes.is_empty()
            || !c.front_modifiers.is_empty()
        {
            return Ok(None);
        }
    }
    let groups: Vec<_> = super::perception::perceive_groups(mol)
        .into_iter()
        .filter(|group| Some(group.key.as_str()) == principal_key)
        .collect();
    let mut components = [(parts, first), (other, second)];
    components.sort_by_key(|(_, (name, _))| spiro_component_sort_key(name));
    let mut projected = AssemblyParts::new(parts.parent_length + other.parent_length - 1);
    projected.is_ring = true;
    projected.omit_redundant_locants = false;
    projected.spiro_prefix_composition = true;
    let mut names = Vec::new();
    let mut junctions = Vec::new();
    let mut lambdas = Vec::new();
    // P-58.2.3.1: indicated and added hydrogen are allocated over the whole
    // spiro system, so the fused components share one completed graph.
    let mut fused = BTreeMap::new();
    for (index, (c, (name, is_fused))) in components.iter().enumerate() {
        let map: BTreeMap<_, _> = c
            .parent_atom_ids_by_locant
            .iter()
            .map(|(l, &id)| {
                (
                    id,
                    if index == 0 {
                        l.clone()
                    } else {
                        prime_locant(l)
                    },
                )
            })
            .collect();
        if map.keys().copied().collect::<AtomSet>() != c.parent_atom_ids {
            return Ok(None);
        }
        let Some(mut local) = project(mol, c, &map, false) else {
            return Ok(None);
        };
        // The monocycle citation already contains its skeletal heteroatoms.
        local.a_prefixes.clear();
        if let Some(group) = local
            .principal_group
            .take_if(|g| Some(g.key.as_str()) != principal_key)
        {
            let mut item = SubstituentItem::new(
                RULES
                    .functional_groups
                    .prefix_for(&group.key)
                    .unwrap_or_default(),
                group.locants,
            );
            item.atom_ids = group.atom_ids;
            item.bond_ids = group.bond_ids;
            local.substituents.push(item);
        }
        super::spirobi::promote_suffixes(mol, &mut local, &groups, &[&map]);
        junctions.push(map[&junction].clone());
        names.push(spiro_component(name));
        lambdas.extend(map.iter().filter_map(|(&id, locant)| {
            Some((locant.clone(), nonstandard_bonding_number(mol, id)?))
        }));
        if *is_fused {
            fused.extend(map);
        }
        projected.substituents.append(&mut local.substituents);
        if let Some(group) = local.principal_group {
            if let Some(principal) = &mut projected.principal_group {
                if principal.key != group.key {
                    return Ok(None);
                }
                principal.locants.extend(group.locants);
                principal.atom_ids.extend(group.atom_ids);
            } else {
                projected.principal_group = Some(group);
            }
        }
    }
    // The numbering proof ranks a fused ketone as a suffix site; as an oxo
    // prefix of a substituent it keeps the unverified fallback.
    if parts.is_substituent
        && super::perception::perceive_groups(mol)
            .iter()
            .any(|g| g.key == "ketone" && fused.contains_key(&g.attachment_carbon))
    {
        return Ok(None);
    }
    let suffix: AtomSet = fused
        .iter()
        .filter(|(_, locant)| {
            projected
                .principal_group
                .as_ref()
                .is_some_and(|g| g.locants.contains(locant))
        })
        .map(|(&id, _)| id)
        .collect();
    // P-58.2.3.1: indicated hydrogen also accommodates the free valence.
    let free_valences: AtomSet = parts
        .parent_atom_ids_by_locant
        .get(&parts.attachment_locant)
        .filter(|id| parts.is_substituent && fused.contains_key(id))
        .into_iter()
        .chain(&suffix)
        .copied()
        .collect();
    let atoms: AtomSet = fused.keys().copied().collect();
    let graph = completed_graph(mol, &atoms, &fused, junction);
    let Some((indicated, _, state)) = super::additive::graph_suffix_hydrogen_state(
        mol,
        &graph,
        &suffix,
        &free_valences,
        &fused,
        Some(junction),
    )?
    else {
        return Ok(None);
    };
    if !state.bond_delta.compatible || !state.unsaturation_operations.is_empty() {
        return Ok(None);
    }
    projected.hydro_operations.push(HydroOperation {
        key: "indicated_hydrogen".into(),
        atom_ids: indicated
            .iter()
            .filter_map(|l| fused.iter().find(|(_, v)| *v == l).map(|(&id, _)| id))
            .collect(),
        locants: indicated,
        operation_kind: "indicated_hydrogen".into(),
        ..Default::default()
    });
    projected.hydro_operations.extend(state.hydro_operations);
    projected
        .hydro_operations
        .extend(state.bond_delta.added_hydrogen_operations);
    let mut suffix = String::new();
    if parts.is_substituent {
        let locant = if std::ptr::eq(components[0].0, parts) {
            parts.attachment_locant.clone()
        } else {
            prime_locant(&parts.attachment_locant)
        };
        if names[1].ends_with('e') {
            names[1].pop();
        }
        suffix = format!("-{locant}-yl");
    } else if projected.principal_group.is_some() {
        let (terminal, rendered, _) =
            super::assembly_parent::format_principal_suffix(&projected, "e")?;
        if terminal.is_empty() && names[1].ends_with('e') {
            names[1].pop();
        }
        suffix = rendered;
    }
    // P-24.8.4.1: lambda descriptors precede `spiro`, after indicated hydrogen.
    lambdas.sort_by_key(|(locant, _)| super::assembly_utils::parse_locant(locant));
    let mut core: String = lambdas
        .iter()
        .map(|(locant, bonding)| format!("{locant}lambda^{bonding}"))
        .collect::<Vec<_>>()
        .join(",");
    if !core.is_empty() {
        core.push('-');
    }
    core.push_str(&format!(
        "spiro[{}-{},{}-{}]",
        names[0], junctions[0], junctions[1], names[1]
    ));
    super::assembler::assemble_completed_parent(&mut projected, &core, &suffix).map(Some)
}

/// Component-local H and mancude ketones require a completed-parent proof.
pub(super) fn needs_hydrogen_proof(mol: &Molecule, parts: &AssemblyParts) -> bool {
    let component_needs_proof = |c: &AssemblyParts| {
        has_nonstandard_bonding(mol, &c.parent_atom_ids)
            || (c.is_simple_ring() && !c.unsaturations.is_empty())
            || !c.hydro_operations.is_empty()
            || !c.indicated_hydrogens.is_empty()
            || c.retained_name.as_deref().is_some_and(|n| {
                n.contains("hydro")
                    || !super::additive::name_indicated_hydrogen_locants(Some(n)).is_empty()
            })
            || ((c.fusion_plan.is_some()
                || c.retained_parent_metadata
                    .as_ref()
                    .is_some_and(|m| m.mancude_double_bonds > 0))
                && super::perception::perceive_groups(mol).iter().any(|group| {
                    group.key == "ketone" && c.parent_atom_ids.contains(&group.attachment_carbon)
                }))
    };
    parts
        .substituents
        .iter()
        .filter_map(|s| s.spiro.as_deref())
        .any(|s| {
            if s.continuation.is_none()
                && let (Ok(Some(first)), Some(other)) = (component(parts), s.side_parts.as_deref())
                && component(other)
                    .ok()
                    .flatten()
                    .is_some_and(|second| first == second)
            {
                // Identical-component fallbacks retain their existing assessment.
                return false;
            }
            component_needs_proof(parts)
                || s.side_parts
                    .as_deref()
                    .is_some_and(|c| component_needs_proof(c) || needs_hydrogen_proof(mol, c))
        })
}
