// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: ReShiki | crates/model/src/molecule_names.rs | 2df788dc456793c687b934afd6d4b317e3c73ef8
// Changes: copy pure identity helpers onto kernel State and forward an explicit Budget; add odd-cumulene E/Z markers.
// Intentionally mirrored on different State types until ReShiki's chemistry and the kernel converge; kept in step by crates/model/src/molecule_names/identity_parity.rs.

//! Canonical molecular identity, independent of document selection and captions.
use crate::{self as chemistry, budget::Budget};
mod axis_markers;
pub use axis_markers::{AXIS_M, AXIS_P, apply_axis_markers, encode_axis_markers};

/// Chemical semantics that cannot be discarded when associating a name.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IdentityRejection {
    UnsupportedAtomsOrBonds,
    NonclassicalStereo,
    /// Preserve the SMILES writer's own diagnostic at both call sites.
    Write(String),
}

/// Reserved naming map for an E odd cumulene; both terminal atoms carry it.
pub const CUMULENE_E: i32 = 900_016;
/// Reserved naming map for a Z odd cumulene; both terminal atoms carry it.
pub const CUMULENE_Z: i32 = 900_032;

/// Record each SMILES-specified odd cumulene as its CIP descriptor
/// (P-93.4.2.3) in the maps of its two terminal atoms.
pub fn mark_cumulenes(prepared: &mut chemistry::smiles::Prepared) -> Result<(), IdentityRejection> {
    let mut markers = Vec::new();
    for ez in &prepared.cumulene_ez {
        let e = chemistry::stereo::cumulene_is_e(&prepared.state, ez)
            .ok()
            .flatten()
            .ok_or(IdentityRejection::NonclassicalStereo)?;
        markers.push((ez.ends, if e { CUMULENE_E } else { CUMULENE_Z }));
    }
    for (end, marker) in markers
        .into_iter()
        .flat_map(|(ends, marker)| ends.map(|end| (end, marker)))
    {
        let atom = prepared
            .state
            .metadata
            .atoms
            .get_mut(end)
            .ok_or(IdentityRejection::NonclassicalStereo)?;
        atom.map_number = marker;
        atom.map_present = true;
    }
    Ok(())
}

/// Marked odd cumulenes and whether each is E. A marker is valid only on both
/// terminal atoms of one stereogenic odd cumulene; any other use fails.
pub fn cumulene_markers(
    state: &chemistry::stereo::perception::State,
) -> Result<Vec<(chemistry::stereo::Cumulene, bool)>, IdentityRejection> {
    let marker = |atom: usize| {
        state
            .metadata
            .atoms
            .get(atom)
            .map(|meta| meta.map_number)
            .filter(|map| matches!(*map, CUMULENE_E | CUMULENE_Z))
    };
    let marked = (0..state.metadata.atoms.len())
        .filter(|&atom| marker(atom).is_some())
        .count();
    if marked == 0 {
        return Ok(Vec::new());
    }
    let mut found = Vec::new();
    for chain in chemistry::stereo::stereogenic_odd_cumulenes(state)
        .map_err(|_| IdentityRejection::NonclassicalStereo)?
    {
        let [(a, _), (b, _)] = chain.ends;
        if let (Some(x), Some(y)) = (marker(a), marker(b))
            && x == y
        {
            found.push((chain, x == CUMULENE_E));
        }
    }
    if marked != 2 * found.len() {
        return Err(IdentityRejection::NonclassicalStereo);
    }
    Ok(found)
}

/// Shared canonical isomeric identity, with private maps for M/P biaryl axes.
/// Ordinary atom maps are excluded. The naming adapter and caption identity
/// retain their own document bounds and error messages.
pub fn canonical_identity(
    state: &chemistry::stereo::perception::State,
    budget: &Budget,
) -> Result<String, IdentityRejection> {
    identity(state, false, budget)
}

/// [`canonical_identity`] of a state read from SMILES: its odd-cumulene E/Z
/// markers are molecular identity and stay in the result.
pub fn canonical_smiles_identity(
    state: &chemistry::stereo::perception::State,
    budget: &Budget,
) -> Result<String, IdentityRejection> {
    identity(state, true, budget)
}

fn identity(
    state: &chemistry::stereo::perception::State,
    cumulenes: bool,
    budget: &Budget,
) -> Result<String, IdentityRejection> {
    if state
        .graph
        .atoms
        .iter()
        .any(|atom| atom.atomic_number == 0 || atom.radical_electrons > 2)
        || state
            .graph
            .bonds
            .iter()
            .any(|bond| !(1..=4).contains(&bond.order))
        || !state.metadata.groups.is_empty()
    {
        return Err(IdentityRejection::UnsupportedAtomsOrBonds);
    }
    let mut encoded = encode_axis_markers(state, budget)?;
    if cumulenes {
        for (chain, e) in cumulene_markers(state)? {
            for (end, _) in chain.ends {
                let atom = encoded
                    .metadata
                    .atoms
                    .get_mut(end)
                    .ok_or(IdentityRejection::NonclassicalStereo)?;
                atom.map_number = if e { CUMULENE_E } else { CUMULENE_Z };
                atom.map_present = true;
            }
        }
    }
    let state = unspecified_unknown_stereo(&encoded);
    let smiles = chemistry::smiles::write::write(
        &state,
        chemistry::smiles::write::Options {
            ignore_maps: false,
            ..Default::default()
        },
    )
    .map_err(|error| IdentityRejection::Write(error.to_string()))?;
    if state
        .graph
        .atoms
        .iter()
        .any(|atom| atom.radical_electrons != 0)
    {
        // Plain SMILES encodes radicals through a valence deficit, rather than
        // an explicit electron field. Admit a radical only when that encoding
        // reconstructs every atom's exact electron count, including atoms that
        // have no radicals. The reader folds ordinary H vertices, so compare
        // against the same detached projection and its aggregate H counts.
        // Keep the caller's graph and atom identities intact; surviving source
        // indices map the writer's traversal onto the projected reader order.
        let projected = chemistry::hydrogens::remove(&chemistry::hydrogens::Input {
            graph: state.graph.clone(),
            metadata: state.metadata.clone(),
            directions: state.directions.clone(),
            unknown_atoms: state
                .properties
                .atoms
                .iter()
                .map(|atom| atom.unknown)
                .collect(),
            annotations: chemistry::hydrogens::Annotations {
                // Only neutral, closed-shell H may become a count. Retaining
                // other H vertices makes a lossy reader reconstruction fail.
                protected_atoms: state
                    .graph
                    .atoms
                    .iter()
                    .enumerate()
                    .filter_map(|(index, atom)| {
                        (atom.atomic_number == 1
                            && (atom.charge != 0 || atom.radical_electrons != 0))
                            .then_some(index)
                    })
                    .collect(),
                ..Default::default()
            },
        })
        .map_err(|_| IdentityRejection::UnsupportedAtomsOrBonds)?;
        let projected_valences = projected
            .graph
            .valences()
            .map_err(|_| IdentityRejection::UnsupportedAtomsOrBonds)?;
        let mut projected_indices = vec![None; state.graph.atoms.len()];
        for (projected_index, &source) in projected.kept_atoms.iter().enumerate() {
            *projected_indices
                .get_mut(source)
                .ok_or(IdentityRejection::UnsupportedAtomsOrBonds)? = Some(projected_index);
        }
        let projected_order = smiles
            .atom_order
            .iter()
            .filter_map(|&source| projected_indices.get(source).copied().flatten())
            .collect::<Vec<_>>();
        let reconstructed = chemistry::smiles::prepare(&smiles.text, budget)
            .map_err(|_| IdentityRejection::UnsupportedAtomsOrBonds)?;
        if projected_order.len() != reconstructed.state.graph.atoms.len()
            || projected_order
                .iter()
                .zip(&reconstructed.state.graph.atoms)
                .enumerate()
                .any(|(target_index, (&source, target))| {
                    projected.graph.atoms.get(source).is_none_or(|atom| {
                        let source_hydrogens = projected_valences.get(source).map(|valence| {
                            u32::from(atom.explicit_hydrogens) + valence.implicit_hydrogens
                        });
                        let target_hydrogens =
                            reconstructed
                                .state
                                .valences
                                .get(target_index)
                                .map(|valence| {
                                    u32::from(target.explicit_hydrogens)
                                        + valence.implicit_hydrogens
                                });
                        atom.atomic_number != target.atomic_number
                            || atom.charge != target.charge
                            || atom.isotope != target.isotope
                            || atom.radical_electrons != target.radical_electrons
                            || source_hydrogens != target_hydrogens
                    })
                })
        {
            return Err(IdentityRejection::UnsupportedAtomsOrBonds);
        }
    }
    Ok(smiles.text)
}

/// The identity named by an unknown-configuration drawing. Clear only the
/// explicitly unknown units; specified neighbouring stereocentres remain part
/// of the identity checked by naming and generated captions.
fn unspecified_unknown_stereo(
    state: &chemistry::stereo::perception::State,
) -> chemistry::stereo::perception::State {
    use chemistry::kekulize::Direction;
    let unknown_atoms = unknown_atom_units(state);
    let mut state = state.clone();
    for ((property, metadata), unknown) in state
        .properties
        .atoms
        .iter_mut()
        .zip(&mut state.metadata.atoms)
        .zip(unknown_atoms)
    {
        if unknown {
            property.unknown = false;
            property.cip_code = None;
            metadata.chiral_tag = 0;
            metadata.chiral_permutation = None;
            metadata.ring_stereo = false;
        }
    }
    for (idx, metadata) in state.metadata.bonds.iter_mut().enumerate() {
        if metadata.unknown_stereo || metadata.stereo == 1 {
            metadata.unknown_stereo = false;
            metadata.stereo = 0;
            metadata.stereo_atoms.clear();
            if let Some(code) = state.properties.bond_codes.get_mut(idx) {
                *code = None;
            }
            if let Some(direction) = state.directions.get_mut(idx) {
                *direction = Direction::None;
            }
        }
    }
    state
}

/// Unknown tetrahedral anchors can be explicit properties or the first
/// endpoint of a native wavy single bond. Double-bond unknowns are separate.
fn unknown_atom_units(state: &chemistry::stereo::perception::State) -> Vec<bool> {
    let mut unknown: Vec<_> = state
        .properties
        .atoms
        .iter()
        .map(|atom| atom.unknown)
        .collect();
    for (idx, bond) in state.graph.bonds.iter().enumerate() {
        if bond.order == 1
            && (state.directions.get(idx) == Some(&chemistry::kekulize::Direction::Unknown)
                || state
                    .metadata
                    .bonds
                    .get(idx)
                    .is_some_and(|meta| meta.unknown_stereo))
            && let Some(anchor) = unknown.get_mut(bond.a)
        {
            *anchor = true;
        }
    }
    unknown
}
