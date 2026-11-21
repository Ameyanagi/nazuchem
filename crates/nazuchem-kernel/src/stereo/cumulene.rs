// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: ReShiki | crates/chemistry/src/stereo/cumulene.rs | be63a799aa03a83dc6ec74bb294f4d42b89298d3

//! Odd cumulene chains with distinguishable ligands at both ends.
// Copied from crates/chemistry/src/stereo/cumulene.rs at be63a799aa03a83dc6ec74bb294f4d42b89298d3; changed: add the specified E/Z record and its CIP descriptor
use super::{at, cip::label, double_bond_directions, perception, perception::State};
use crate::ranking;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Cumulene {
    /// Each terminal atom and its inner chain neighbor, in graph atom indices.
    pub ends: [(usize, usize); 2],
    /// Graph bond indices in order from the first end to the second.
    pub bonds: Vec<usize>,
}

/// Cis/trans geometry that SMILES directions specify across an odd cumulene.
/// `refs[i]` is a ligand of the terminal atom `ends[i]`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CumuleneEz {
    pub ends: [usize; 2],
    pub refs: [usize; 2],
    pub cis: bool,
}

/// Whether a specified odd cumulene is E (P-93.4.2.3): the CIP-senior ligands
/// of its two terminal atoms are trans. `None` when CIP cannot rank an end.
pub fn cumulene_is_e(state: &State, ez: &CumuleneEz) -> Result<Option<bool>, String> {
    let double = |atom: usize, skip: usize| {
        state
            .graph
            .bonds
            .iter()
            .enumerate()
            .find_map(|(id, bond)| match (bond.a, bond.b) {
                (x, y) | (y, x) if bond.order == 2 && x == atom && y != skip => Some((id, y)),
                _ => None,
            })
            .ok_or("Invalid cumulene end")
    };
    // Label each terminal double bond as cis from the reference to the next
    // chain atom: it is Z exactly when the reference is the senior ligand.
    let mut probe = state.clone();
    let mut bonds = [0; 2];
    for ((&end, &reference), slot) in ez.ends.iter().zip(&ez.refs).zip(&mut bonds) {
        let (bond, inner) = double(end, usize::MAX)?;
        let (_, next) = double(inner, end)?;
        let meta = probe
            .metadata
            .bonds
            .get_mut(bond)
            .ok_or("Invalid cumulene bond")?;
        meta.stereo = 4;
        meta.stereo_atoms = if at(&state.graph.bonds, bond)?.a == end {
            vec![reference, next]
        } else {
            vec![next, reference]
        };
        *slot = bond;
    }
    let codes = label::assign(&probe, &label::Options::default())
        .map_err(|error| error.to_string())?
        .state
        .properties
        .bond_codes;
    let mut senior = [false; 2];
    for (bond, senior) in bonds.into_iter().zip(&mut senior) {
        *senior = match at(&codes, bond)?.as_deref() {
            Some("Z") => true,
            Some("E") => false,
            _ => return Ok(None),
        };
    }
    Ok(Some(ez.cis != (senior[0] == senior[1])))
}

#[derive(PartialEq, Eq)]
enum Ligand {
    Atom(u32),
    Hydrogen,
    LonePair,
}

fn without_maps(state: &State) -> Result<State, String> {
    let mut unmapped = state.clone();
    for atom in &mut unmapped.metadata.atoms {
        atom.map_number = 0;
        atom.map_present = false;
    }
    // Directions reconstructed for adjacent specified bonds can flank an
    // unspecified double bond. Temporarily mark it unknown to prevent the
    // second perception from inventing a configuration there.
    for (bond, meta) in unmapped
        .graph
        .bonds
        .iter()
        .zip(&mut unmapped.metadata.bonds)
    {
        if bond.order == 2 && meta.stereo == 0 {
            meta.stereo = 1;
        }
    }
    let geometry = double_bond_directions(
        &unmapped.graph,
        &unmapped.metadata,
        &unmapped.directions,
        None,
        &unmapped.rings.atoms,
    )?;
    unmapped.metadata = geometry.metadata;
    unmapped.directions = geometry.directions;
    let mut unmapped = perception::perceive(
        &unmapped,
        perception::Options {
            clean: true,
            force: true,
            flag_possible: false,
        },
    )?;
    for ((bond, original), meta) in state
        .graph
        .bonds
        .iter()
        .zip(&state.metadata.bonds)
        .zip(&mut unmapped.metadata.bonds)
    {
        if bond.order == 2 && original.stereo == 0 {
            *meta = original.clone();
        }
    }
    Ok(unmapped)
}

/// Maximal chains of an odd number (at least three) of double bonds, whatever
/// their terminal ligands.
pub fn odd_cumulenes(state: &State) -> Result<Vec<Cumulene>, String> {
    Ok(chains(state)?.0)
}

type Neighbors = Vec<Vec<(usize, usize)>>;

fn chains(state: &State) -> Result<(Vec<Cumulene>, Neighbors), String> {
    let mut neighbors = vec![Vec::new(); state.graph.atoms.len()];
    let mut doubles = neighbors.clone();
    for (id, bond) in state.graph.bonds.iter().enumerate() {
        for (end, other) in [(bond.a, bond.b), (bond.b, bond.a)] {
            neighbors
                .get_mut(end)
                .ok_or("Invalid cumulene atom index")?
                .push((other, id));
            if bond.order == 2 {
                doubles
                    .get_mut(end)
                    .ok_or("Invalid cumulene atom index")?
                    .push((other, id));
            }
        }
    }
    let mut chains = Vec::new();
    for (start, adjacent) in doubles.iter().enumerate() {
        let &[(first, bond)] = adjacent.as_slice() else {
            continue;
        };
        let (mut previous, mut end) = (start, first);
        let mut bonds = vec![bond];
        while let &[a, b] = at(&doubles, end)?.as_slice() {
            let (next, bond) = if a.0 == previous { b } else { a };
            bonds.push(bond);
            (previous, end) = (end, next);
        }
        if start < end && bonds.len() >= 3 && bonds.len() % 2 == 1 && at(&doubles, end)?.len() == 1
        {
            chains.push(Cumulene {
                ends: [(start, first), (end, previous)],
                bonds,
            });
        }
    }
    Ok((chains, neighbors))
}

/// Find maximal chains of an odd number (at least three) of double bonds whose
/// ends each have two different ligands. The state must be prepared, with
/// perceived ordinary stereo and caches belonging to the same graph. Isotopes
/// and stereo distinguish ligands; atom maps do not. No configuration is assigned.
pub fn stereogenic_odd_cumulenes(state: &State) -> Result<Vec<Cumulene>, String> {
    let (chains, neighbors) = chains(state)?;
    if chains.is_empty() {
        return Ok(chains);
    }
    // Legacy perception itself uses maps in its priorities. Merely excluding
    // maps from ranking would retain map-dependent E/Z or tetrahedral labels.
    let unmapped = state
        .metadata
        .atoms
        .iter()
        .any(|a| a.map_present || a.map_number != 0)
        .then(|| without_maps(state))
        .transpose()?;
    let state = unmapped.as_ref().unwrap_or(state);
    let ranks = ligand_ranks(state, false)?;
    let mut result = Vec::new();
    for chain in chains {
        if distinct_ends(state, &neighbors, &ranks, chain.ends)? {
            result.push(chain);
        }
    }
    Ok(result)
}

fn ligand_ranks(state: &State, cumulene_labels: bool) -> Result<Vec<u32>, String> {
    let ranks = ranking::rank(
        &state.graph,
        &state.rings.atoms,
        &state.metadata,
        ranking::Options {
            break_ties: false,
            include_chirality: true,
            include_isotopes: true,
            include_maps: cumulene_labels,
            // Slash geometry refers to the depicted member of an enhanced group.
            include_stereo_groups: false,
            ..ranking::Options::default()
        },
    )?;
    Ok(ranks)
}

fn distinct_ends(
    state: &State,
    neighbors: &Neighbors,
    ranks: &[u32],
    ends: [(usize, usize); 2],
) -> Result<bool, String> {
    let distinct = |(end, inner)| -> Result<bool, String> {
        let atom = at(&state.graph.atoms, end)?;
        let hydrogens =
            u32::from(atom.explicit_hydrogens) + at(&state.valences, end)?.implicit_hydrogens;
        let mut ligands = Vec::new();
        for &(other, _) in at(neighbors, end)? {
            if other != inner {
                let atom = at(&state.graph.atoms, other)?;
                // An explicit 1H is the same ligand as H (P-92.3).
                ligands.push(if atom.atomic_number == 1 && atom.isotope <= 1 {
                    Ligand::Hydrogen
                } else {
                    Ligand::Atom(*at(ranks, other)?)
                });
            }
        }
        ligands.extend((0..hydrogens).map(|_| Ligand::Hydrogen));
        while ligands.len() < 2 {
            ligands.push(Ligand::LonePair);
        }
        Ok(matches!(ligands.as_slice(), [a, b] if a != b))
    };
    Ok(distinct(ends[0])? && distinct(ends[1])?)
}

/// Unconfigured ordinary double bonds that could acquire a slash on both ends.
pub(crate) fn unspecified_directed_alkenes(
    state: &State,
    marked: &[Vec<usize>],
) -> Result<Vec<[usize; 2]>, String> {
    let mut candidates = vec![false; state.graph.bonds.len()];
    for ((bond, meta), candidate) in state
        .graph
        .bonds
        .iter()
        .zip(&state.metadata.bonds)
        .zip(&mut candidates)
    {
        *candidate = bond.order == 2
            && meta.stereo <= 1
            && !at(marked, bond.a)?.is_empty()
            && !at(marked, bond.b)?.is_empty();
    }
    if !candidates.contains(&true) {
        return Ok(Vec::new());
    }
    let unmapped = state
        .metadata
        .atoms
        .iter()
        .any(|a| a.map_present || a.map_number != 0)
        .then(|| without_maps(state))
        .transpose()?;
    let state = unmapped.as_ref().unwrap_or(state);
    let mut result = Vec::new();
    for id in perception::stereogenic_units(state)?.1 {
        if *at(&candidates, id)? {
            let bond = at(&state.graph.bonds, id)?;
            result.push([bond.a, bond.b]);
        }
    }
    Ok(result)
}

/// Reject contradictory directions around an ordinary stereogenic double bond
/// and report directed units whose ligands require another cumulene's descriptor.
pub(crate) fn check_stereo_alongside_cumulenes(
    state: &State,
    configurations: &[CumuleneEz],
    directed: &[bool],
    conflicts: &[bool],
) -> Result<bool, String> {
    let (chains, neighbors) = chains(state)?;
    let both_directed = |ends: [(usize, usize); 2]| {
        ends.iter()
            .all(|&(end, _)| directed.get(end) == Some(&true))
    };
    let mut unresolved: Vec<_> = chains
        .into_iter()
        .map(|chain| chain.ends)
        .filter(|&ends| {
            both_directed(ends)
                && !configurations
                    .iter()
                    .any(|ez| ez.ends == ends.map(|(end, _)| end))
        })
        .collect();
    let mut contradictory = Vec::new();
    let ordinary_end = |end: usize| -> Result<bool, String> {
        let double_count = at(&neighbors, end)?
            .iter()
            .filter(|&&(_, id)| {
                state
                    .graph
                    .bonds
                    .get(id)
                    .is_some_and(|bond| bond.order == 2)
            })
            .count();
        Ok(double_count == 1)
    };
    for (bond, meta) in state.graph.bonds.iter().zip(&state.metadata.bonds) {
        let ends = [(bond.a, bond.b), (bond.b, bond.a)];
        if bond.order != 2
            || !both_directed(ends)
            || !ordinary_end(bond.a)?
            || !ordinary_end(bond.b)?
        {
            continue;
        }
        if ends
            .iter()
            .any(|&(end, _)| conflicts.get(end) == Some(&true))
        {
            contradictory.push(ends);
        }
        if meta.stereo == 0 {
            unresolved.push(ends);
        }
    }
    if contradictory.is_empty() && unresolved.is_empty() {
        return Ok(false);
    }
    let baseline = without_maps(state)?;
    let ranks = ligand_ranks(&baseline, false)?;
    for ends in contradictory {
        if distinct_ends(&baseline, &neighbors, &ranks, ends)? {
            return Err("contradictory SMILES directions around an ordinary double bond".into());
        }
    }
    if unresolved.is_empty() {
        return Ok(false);
    }
    let mut labelled = baseline.clone();
    for ez in configurations {
        let label = match cumulene_is_e(&baseline, ez)? {
            Some(true) => 1,
            Some(false) => 2,
            // A known geometry with unknown CIP still merits a conservative
            // check for dependency.
            None => 3,
        };
        for &end in &ez.ends {
            let atom = labelled
                .metadata
                .atoms
                .get_mut(end)
                .ok_or("Invalid cumulene atom index")?;
            atom.map_number = label;
            atom.map_present = true;
        }
    }
    let with_descriptors = ligand_ranks(&labelled, true)?;
    for ends in unresolved {
        if !distinct_ends(&baseline, &neighbors, &ranks, ends)?
            && distinct_ends(&baseline, &neighbors, &with_descriptors, ends)?
        {
            return Ok(true);
        }
    }
    Ok(false)
}

#[cfg(test)]
mod tests;
