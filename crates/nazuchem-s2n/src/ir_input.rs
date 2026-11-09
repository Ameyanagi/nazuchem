// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | crates/nazuchem-s2n/src/ir_input.rs | Phase B step 8b
// Copyright (c) 2026 Ameyanagi and Nazuchem contributors

//! Direct IR admission and chemistry preparation. Atom and bond positions are
//! identities throughout; this adapter never serializes or reparses a graph.

use nazuchem_ir::{
    AxialWinding, BondOrder, Molecule, Ref, Relation, Stereo, StereoGroupKind, Winding,
};
use nazuchem_kernel::{
    budget::Budget,
    graph::{Atom, Bond, Graph},
    kekulize::Direction,
    ranking::{Metadata, StereoGroup},
    sanitize,
    stereo::{self, cip::configuration::parity4, perception},
};

const TETRAHEDRAL_LOSS: &str = "Native import could not retain every specified tetrahedral center; no name or graph was assigned";
const CUMULENE_STEREO: &str =
    "Naming does not support cumulene E/Z configuration yet; no less specific name was substituted";
const NONCLASSICAL_STEREO: &str = "Naming supports tetrahedral, ordinary alkene and biaryl axis stereochemistry; other stereo classes cannot be discarded";

/// Build the kernel state directly from the public, localized graph.
///
/// IR hydrogen counts are authoritative. The kernel's inferred representation
/// is recovered where it gives exactly that count; otherwise the count remains
/// frozen. The naming graph reads hydrogen *origin* from the original IR.
pub(crate) fn prepare(molecule: &Molecule, budget: &Budget) -> Result<perception::State, String> {
    if molecule.atoms.len() > 512 {
        return Err(
            "Naming supports at most 512 input atom tokens; no oversized graph was imported".into(),
        );
    }
    molecule.validate().map_err(|error| error.to_string())?;
    if !molecule.position_variation.is_empty()
        || !molecule.polymers.is_empty()
        || molecule
            .atoms
            .iter()
            .any(|atom| atom.label.is_some() || atom.attachment_point.is_some())
    {
        return Err(
            "Naming requires ordinary molecular SMILES without reaction or CX extensions".into(),
        );
    }
    let graph = Graph {
        atoms: molecule
            .atoms
            .iter()
            .map(|atom| Atom {
                atomic_number: atom.element,
                isotope: atom.isotope.unwrap_or_default(),
                charge: atom.charge,
                explicit_hydrogens: atom.hydrogens,
                no_implicit: true,
                aromatic: false,
                radical_electrons: atom.radical_electrons,
            })
            .collect(),
        bonds: molecule
            .bonds
            .iter()
            .map(|bond| {
                Ok(Bond {
                    a: bond.atoms[0] as usize,
                    b: bond.atoms[1] as usize,
                    order: match bond.order {
                        BondOrder::Single => 1,
                        BondOrder::Double => 2,
                        BondOrder::Triple => 3,
                        BondOrder::Quadruple => 6,
                        _ => return Err("Unsupported IR bond order".into()),
                    },
                    aromatic: false,
                })
            })
            .collect::<Result<_, String>>()?,
    };
    let adjacent = neighbors(&graph);
    let mut metadata = Metadata::unspecified(&graph);
    for (atom, meta) in molecule.atoms.iter().zip(&mut metadata.atoms) {
        if let Some(map) = atom.atom_map {
            // The native SMILES integer reader rejects values at or above this
            // bound before its i32 cast. Preserve that naming admission limit;
            // the public IR itself continues to retain the complete u32 range.
            if map >= 2_147_483_640 {
                return Err("SMILES exceeds parser resource limits".into());
            }
            meta.map_number = map as i32;
            meta.map_present = true;
        }
    }
    for entry in &molecule.stereo {
        match entry {
            Stereo::Tetrahedral {
                center,
                neighbors,
                winding,
            } => {
                let center = *center as usize;
                if *winding != Winding::Unknown {
                    metadata.atoms[center].chiral_tag =
                        tetrahedral_tag(center, neighbors, *winding, &adjacent)?;
                }
            }
            Stereo::DoubleBond { bond, relation, .. } => {
                let bond = *bond as usize;
                let axis = &graph.bonds[bond];
                if *relation != Relation::Unknown
                    && graph.bonds.iter().enumerate().any(|(id, other)| {
                        id != bond
                            && other.order == 2
                            && [axis.a, axis.b]
                                .iter()
                                .any(|end| other.a == *end || other.b == *end)
                    })
                {
                    return Err(CUMULENE_STEREO.into());
                }
            }
            Stereo::Atropisomer {
                bond,
                refs,
                winding,
            } => {
                let bond = *bond as usize;
                metadata.bonds[bond].stereo = axial_tag(&graph, &adjacent, bond, refs, *winding)?;
                metadata.bonds[bond].stereo_atoms = axial_carriers(&graph, &adjacent, bond)?;
            }
            _ => return Err(NONCLASSICAL_STEREO.into()),
        }
    }
    let sanitized = sanitize::sanitize(
        &graph,
        &metadata,
        &vec![Direction::None; graph.bonds.len()],
        budget,
    )
    .map_err(|error| error.to_string())?;
    let mut state = perception::State {
        properties: perception::Properties::unspecified(&sanitized.graph),
        graph: sanitized.graph,
        metadata: sanitized.metadata,
        directions: sanitized.directions,
        valences: sanitized.valences,
        conjugated: sanitized.conjugated,
        hybridizations: sanitized.hybridizations,
        rings: perception::RingCache {
            kind: perception::RingKind::Symmetric,
            atoms: sanitized.rings,
        },
    };
    if state.graph.atoms.len() != molecule.atoms.len()
        || state.graph.bonds.len() != molecule.bonds.len()
        || state
            .graph
            .bonds
            .iter()
            .zip(&molecule.bonds)
            .any(|(bond, source)| [bond.a as u32, bond.b as u32] != source.atoms)
    {
        return Err("Native import could not retain IR atom and bond identities".into());
    }
    // Sanitization may normalize resonance but must never silently reinterpret
    // a complete IR hydrogen or electron count as an unspecified SMILES atom.
    for (index, ((source, atom), valence)) in molecule
        .atoms
        .iter()
        .zip(&state.graph.atoms)
        .zip(&state.valences)
        .enumerate()
    {
        if u32::from(atom.explicit_hydrogens) + valence.implicit_hydrogens
            != u32::from(source.hydrogens)
        {
            return Err(format!(
                "Native import could not retain IR hydrogen count at atom {index}"
            ));
        }
        if atom.radical_electrons != source.radical_electrons {
            return Err(format!(
                "Native import could not retain IR radical electron count at atom {index}"
            ));
        }
    }
    // Thaw one atom at a time: each inference is local to that atom and must
    // reproduce the stored count, including aromatic heteroatom hydrogens.
    for (index, atom) in molecule.atoms.iter().enumerate() {
        if atom.hydrogens_stated {
            continue;
        }
        let previous = state.graph.atoms[index].clone();
        state.graph.atoms[index].explicit_hydrogens = 0;
        state.graph.atoms[index].no_implicit = false;
        match state.graph.valences() {
            Ok(valences) if valences[index].implicit_hydrogens == u32::from(atom.hydrogens) => {
                state.valences = valences
            }
            _ => state.graph.atoms[index] = previous,
        }
    }
    state.valences = state.graph.valences()?;
    if molecule.stereo.iter().any(|entry| {
        matches!(
            entry,
            Stereo::Atropisomer {
                winding: AxialWinding::Unknown,
                ..
            }
        )
    }) {
        // Unknown winding does not relax the shared naming topology contract.
        // At this point every axis still carries a detached trial winding, so
        // the ordinary biaryl encoder can validate it before Unknown becomes
        // native stereo=1. Its result is discarded: no configuration or atom
        // map is asserted in the returned state.
        nazuchem_kernel::identity::encode_axis_markers(&state, budget)
            .map_err(|_| NONCLASSICAL_STEREO.to_owned())?;
    }
    for entry in &molecule.stereo {
        match entry {
            Stereo::Tetrahedral {
                center,
                winding: Winding::Unknown,
                ..
            } => state.properties.atoms[*center as usize].unknown = true,
            Stereo::Atropisomer { bond, winding, .. } => {
                let bond = *bond as usize;
                if !matches!(state.metadata.bonds[bond].stereo, 6 | 7) {
                    return Err(NONCLASSICAL_STEREO.into());
                }
                if *winding == AxialWinding::Unknown {
                    state.metadata.bonds[bond].stereo = 1;
                    state.metadata.bonds[bond].unknown_stereo = true;
                }
            }
            _ => {}
        }
    }
    // The ordinary clean pass validates tetrahedral winding. Double-bond IR
    // entries are installed afterward: legacy cleanup otherwise discards their
    // authoritative cis/trans relation merely because no slash carriers exist.
    state = perception::perceive(
        &state,
        perception::Options {
            clean: true,
            force: true,
            flag_possible: true,
        },
    )?;
    let (_, stereogenic_bonds) = perception::stereogenic_units(&state)?;
    for entry in &molecule.stereo {
        match entry {
            Stereo::Tetrahedral {
                center, winding, ..
            } => {
                let center = *center as usize;
                if *winding != Winding::Unknown && state.metadata.atoms[center].chiral_tag == 0 {
                    return Err(TETRAHEDRAL_LOSS.into());
                }
            }
            Stereo::DoubleBond {
                bond,
                refs,
                relation,
            } => {
                let bond = *bond as usize;
                // Explicit Unknown retains the same structural annotation as
                // the existing private-marker import, even for equivalent
                // ligands. It never asserts a configuration.
                if *relation != Relation::Unknown && !stereogenic_bonds.contains(&bond) {
                    return Err(NONCLASSICAL_STEREO.into());
                }
                state.metadata.bonds[bond].stereo = match relation {
                    Relation::Cis => 4,
                    Relation::Trans => 5,
                    Relation::Unknown => 1,
                    _ => return Err(NONCLASSICAL_STEREO.into()),
                };
                state.metadata.bonds[bond].stereo_atoms =
                    refs.iter().map(|&atom| atom as usize).collect();
                state.metadata.bonds[bond].unknown_stereo = *relation == Relation::Unknown;
            }
            Stereo::Atropisomer { bond, .. } => {
                if !stereo::axis_ends_distinct(&state.graph, *bond as usize, budget)? {
                    return Err(NONCLASSICAL_STEREO.into());
                }
            }
            _ => return Err(NONCLASSICAL_STEREO.into()),
        }
    }
    // Enhanced groups are mapped without flattening their relative/racemic
    // meaning. The shared naming identity gate decides their admission.
    for source in &molecule.stereo_groups {
        let (kind, id) = match source.kind {
            StereoGroupKind::And(id) => (2, id),
            StereoGroupKind::Or(id) => (1, id),
            _ => return Err(NONCLASSICAL_STEREO.into()),
        };
        let mut group = StereoGroup {
            kind,
            read_id: id,
            write_id: id,
            ..StereoGroup::default()
        };
        for &member in &source.stereo {
            match molecule.stereo[member as usize] {
                Stereo::Tetrahedral { center, .. } => group.atoms.push(center as usize),
                Stereo::DoubleBond { bond, .. } | Stereo::Atropisomer { bond, .. } => {
                    group.bonds.push(bond as usize)
                }
                _ => return Err(NONCLASSICAL_STEREO.into()),
            }
        }
        state.metadata.groups.push(group);
    }
    state.metadata.validate(&state.graph)?;
    Ok(state)
}

fn neighbors(graph: &Graph) -> Vec<Vec<usize>> {
    let mut adjacent = vec![Vec::new(); graph.atoms.len()];
    for bond in &graph.bonds {
        adjacent[bond.a].push(bond.b);
        adjacent[bond.b].push(bond.a);
    }
    adjacent
}

fn tetrahedral_tag(
    center: usize,
    refs: &[Ref; 4],
    winding: Winding,
    adjacent: &[Vec<usize>],
) -> Result<u8, String> {
    let mut carriers = adjacent[center]
        .iter()
        .copied()
        .map(Some)
        .collect::<Vec<_>>();
    if carriers.len() < 4 {
        carriers.push(Some(center));
    }
    if carriers.len() < 4 {
        carriers.push(None);
    }
    let refs = refs
        .map(|reference| match reference {
            Ref::Atom(atom) => Ok(Some(atom as usize)),
            Ref::ImplicitHydrogen => Ok(Some(center)),
            Ref::LonePair if adjacent[center].len() == 3 => Ok(Some(center)),
            Ref::LonePair => Ok(None),
            _ => Err(NONCLASSICAL_STEREO.to_owned()),
        })
        .into_iter()
        .collect::<Result<Vec<_>, _>>()?;
    let parity = parity4(&refs, &carriers).map_err(|error| error.to_string())?;
    if parity == 0 {
        return Err(TETRAHEDRAL_LOSS.into());
    }
    match winding {
        Winding::Clockwise => Ok(if parity == 2 { 1 } else { 2 }),
        Winding::Anticlockwise => Ok(if parity == 2 { 2 } else { 1 }),
        _ => Err(NONCLASSICAL_STEREO.into()),
    }
}

fn axial_carriers(
    graph: &Graph,
    adjacent: &[Vec<usize>],
    bond: usize,
) -> Result<Vec<usize>, String> {
    let axis = &graph.bonds[bond];
    [(axis.a, axis.b), (axis.b, axis.a)]
        .into_iter()
        .map(|(end, other)| {
            let refs = adjacent[end]
                .iter()
                .copied()
                .filter(|&atom| atom != other)
                .collect::<Vec<_>>();
            if refs.len() != 2 {
                return Err(NONCLASSICAL_STEREO.into());
            }
            refs.into_iter()
                .min()
                .ok_or_else(|| NONCLASSICAL_STEREO.into())
        })
        .collect()
}

fn axial_tag(
    graph: &Graph,
    adjacent: &[Vec<usize>],
    bond: usize,
    refs: &[Ref; 4],
    winding: AxialWinding,
) -> Result<u8, String> {
    let carriers = axial_carriers(graph, adjacent, bond)?;
    if refs
        .iter()
        .any(|reference| !matches!(reference, Ref::Atom(_)))
    {
        return Err(NONCLASSICAL_STEREO.into());
    }
    let mut clockwise = match winding {
        AxialWinding::Clockwise | AxialWinding::Unknown => true,
        AxialWinding::Anticlockwise => false,
        _ => return Err(NONCLASSICAL_STEREO.into()),
    };
    for side in 0..2 {
        if let Ref::Atom(reference) = refs[side * 2] {
            clockwise ^= reference as usize != carriers[side];
        }
    }
    Ok(if clockwise { 6 } else { 7 })
}

/// Test-only bridge for comparing existing SMILES corpora with the direct IR
/// route. Prepared atom/bond order is retained; no writer determines identity.
#[cfg(test)]
pub(crate) fn from_smiles(text: &str, budget: &Budget) -> Result<Molecule, String> {
    use nazuchem_ir::{Atom as IrAtom, Bond as IrBond, StereoGroup as IrGroup};
    use nazuchem_kernel::{identity, kekulize, ranking, smiles};
    let prepared = smiles::read(text, budget).map_err(|error| error.to_string())?;
    if prepared.dropped_cumulene_ez {
        return Err(CUMULENE_STEREO.into());
    }
    let mut state = identity::apply_axis_markers(&prepared.prepared.state, budget)
        .map_err(|error| format!("{error:?}"))?;
    let ranks = ranking::rank_cached(
        &state.graph,
        &state.rings.atoms,
        &state.metadata,
        ranking::Options {
            fragment: true,
            ..Default::default()
        },
        &state.valences,
    )?;
    let localized = kekulize::assign(
        &state.graph,
        &state.rings.atoms,
        &state.directions,
        kekulize::Options {
            ranks: Some(&ranks),
            ..Default::default()
        },
    )?;
    state.graph = localized.graph;
    state.directions = localized.directions;
    state.valences = state.graph.valences()?;
    let adjacent = neighbors(&state.graph);
    let mut molecule = Molecule {
        atoms: state
            .graph
            .atoms
            .iter()
            .zip(&state.valences)
            .zip(&state.metadata.atoms)
            .map(|((atom, valence), meta)| {
                Ok(IrAtom {
                    element: atom.atomic_number,
                    isotope: (atom.isotope != 0).then_some(atom.isotope),
                    charge: atom.charge,
                    hydrogens: u8::try_from(
                        u32::from(atom.explicit_hydrogens) + valence.implicit_hydrogens,
                    )
                    .map_err(|_| "IR hydrogen count overflow")?,
                    hydrogens_stated: atom.no_implicit,
                    radical_electrons: atom.radical_electrons,
                    atom_map: (meta.map_present || meta.map_number != 0)
                        .then_some(meta.map_number as u32),
                    ..IrAtom::default()
                })
            })
            .collect::<Result<_, String>>()?,
        bonds: state
            .graph
            .bonds
            .iter()
            .map(|bond| {
                Ok(IrBond {
                    atoms: [bond.a as u32, bond.b as u32],
                    order: match bond.order {
                        1 => BondOrder::Single,
                        2 => BondOrder::Double,
                        3 => BondOrder::Triple,
                        6 => BondOrder::Quadruple,
                        _ => return Err("Unlocalized test bond".to_owned()),
                    },
                })
            })
            .collect::<Result<_, String>>()?,
        ..Molecule::default()
    };
    let mut atom_entries = vec![None; molecule.atoms.len()];
    let mut bond_entries = vec![None; molecule.bonds.len()];
    let unknown_atoms = state
        .metadata
        .atoms
        .iter()
        .map(|meta| matches!(meta.map_number, 900_001 | 900_003))
        .collect::<Vec<_>>();
    let unknown_bond_ends = state
        .metadata
        .atoms
        .iter()
        .map(|meta| matches!(meta.map_number, 900_002 | 900_003))
        .collect::<Vec<_>>();
    for (index, meta) in state.metadata.atoms.iter().enumerate() {
        if matches!(meta.map_number, 900_001..=900_003) {
            molecule.atoms[index].atom_map = None;
        }
        let unknown = unknown_atoms[index] || state.properties.atoms[index].unknown;
        if !matches!(meta.chiral_tag, 1 | 2) && !unknown {
            if meta.chiral_tag != 0 {
                return Err(NONCLASSICAL_STEREO.into());
            }
            continue;
        }
        let mut refs = adjacent[index]
            .iter()
            .map(|&atom| Ref::Atom(atom as u32))
            .collect::<Vec<_>>();
        if refs.len() < 4 {
            refs.push(if molecule.atoms[index].hydrogens > 0 {
                Ref::ImplicitHydrogen
            } else {
                Ref::LonePair
            });
        }
        if refs.len() < 4 {
            refs.push(Ref::LonePair);
        }
        let neighbors = refs.try_into().map_err(|_| TETRAHEDRAL_LOSS.to_owned())?;
        atom_entries[index] = Some(molecule.stereo.len() as u32);
        molecule.stereo.push(Stereo::Tetrahedral {
            center: index as u32,
            neighbors,
            winding: if unknown {
                Winding::Unknown
            } else if meta.chiral_tag == 1 {
                Winding::Clockwise
            } else {
                Winding::Anticlockwise
            },
        });
    }
    for (index, (bond, meta)) in state
        .graph
        .bonds
        .iter()
        .zip(&state.metadata.bonds)
        .enumerate()
    {
        let entry = if bond.order == 2
            && (meta.stereo != 0 || unknown_bond_ends[bond.a] && unknown_bond_ends[bond.b])
        {
            let refs = if meta.stereo_atoms.len() == 2 {
                [meta.stereo_atoms[0] as u32, meta.stereo_atoms[1] as u32]
            } else {
                let reference = |end: usize, other: usize| {
                    adjacent[end]
                        .iter()
                        .copied()
                        .find(|&atom| atom != other)
                        .map(|atom| atom as u32)
                        .ok_or_else(|| NONCLASSICAL_STEREO.to_owned())
                };
                [reference(bond.a, bond.b)?, reference(bond.b, bond.a)?]
            };
            Some(Stereo::DoubleBond {
                bond: index as u32,
                refs,
                relation: if meta.stereo == 1
                    || unknown_bond_ends[bond.a] && unknown_bond_ends[bond.b]
                {
                    Relation::Unknown
                } else if matches!(meta.stereo, 2 | 4) {
                    Relation::Cis
                } else {
                    Relation::Trans
                },
            })
        } else if matches!(meta.stereo, 6 | 7)
            || bond.order == 1 && meta.stereo == 1 && meta.unknown_stereo
        {
            let mut refs = Vec::new();
            for (end, other) in [(bond.a, bond.b), (bond.b, bond.a)] {
                let mut pair = adjacent[end]
                    .iter()
                    .copied()
                    .filter(|&atom| atom != other)
                    .collect::<Vec<_>>();
                if pair.len() != 2 {
                    return Err(NONCLASSICAL_STEREO.into());
                }
                pair.sort_unstable();
                refs.extend(pair.into_iter().map(|atom| Ref::Atom(atom as u32)));
            }
            Some(Stereo::Atropisomer {
                bond: index as u32,
                refs: refs
                    .try_into()
                    .map_err(|_| NONCLASSICAL_STEREO.to_owned())?,
                winding: if meta.stereo == 1 {
                    AxialWinding::Unknown
                } else if meta.stereo == 6 {
                    AxialWinding::Clockwise
                } else {
                    AxialWinding::Anticlockwise
                },
            })
        } else {
            None
        };
        if let Some(entry) = entry {
            bond_entries[index] = Some(molecule.stereo.len() as u32);
            molecule.stereo.push(entry);
        }
    }
    for group in &state.metadata.groups {
        let kind = match group.kind {
            0 => continue,
            1 => StereoGroupKind::Or(group.read_id),
            2 => StereoGroupKind::And(group.read_id),
            _ => return Err(NONCLASSICAL_STEREO.into()),
        };
        let entries = group
            .atoms
            .iter()
            .map(|&atom| atom_entries[atom])
            .chain(group.bonds.iter().map(|&bond| bond_entries[bond]))
            .map(|entry| entry.ok_or_else(|| NONCLASSICAL_STEREO.to_owned()))
            .collect::<Result<_, _>>()?;
        molecule.stereo_groups.push(IrGroup {
            kind,
            stereo: entries,
        });
    }
    molecule.validate().map_err(|error| error.to_string())?;
    Ok(molecule)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nazuchem_ir::{NameReason, NameStatus, StereoGroup as IrGroup};
    use nazuchem_kernel::stereo::cip::label;

    fn input(text: &str) -> Molecule {
        from_smiles(text, &Budget::default()).unwrap()
    }

    fn identity(molecule: &Molecule) -> String {
        crate::identity::canonical_state(
            &prepare(molecule, &Budget::default()).unwrap(),
            &Budget::default(),
        )
        .unwrap()
    }

    fn complete_name(molecule: &Molecule) -> serde_json::Value {
        serde_json::to_value(
            crate::name_molecule(
                molecule,
                nazuchem_ir::IsotopeStyle::default(),
                &Budget::default(),
            )
            .unwrap(),
        )
        .unwrap()
    }

    #[test]
    fn direct_preparation_retains_identity_order_hydrogen_origin_and_maps() {
        for text in [
            "[13CH3:7][C@@H:0](F)Cl",
            "[2H]O[CH2:2147483639][CH2]",
            "c1cc[nH]c1",
        ] {
            let molecule = input(text);
            let state = prepare(&molecule, &Budget::default()).unwrap();
            assert_eq!(state.graph.atoms.len(), molecule.atoms.len(), "{text}");
            assert_eq!(state.graph.bonds.len(), molecule.bonds.len(), "{text}");
            for (index, source) in molecule.atoms.iter().enumerate() {
                let atom = &state.graph.atoms[index];
                assert_eq!(atom.atomic_number, source.element, "{text}: atom {index}");
                assert_eq!(atom.isotope, source.isotope.unwrap_or_default());
                assert_eq!(atom.charge, source.charge);
                assert_eq!(atom.radical_electrons, source.radical_electrons);
                assert_eq!(
                    u32::from(atom.explicit_hydrogens) + state.valences[index].implicit_hydrogens,
                    u32::from(source.hydrogens)
                );
                assert_eq!(
                    state.metadata.atoms[index].map_present,
                    source.atom_map.is_some()
                );
                assert_eq!(
                    state.metadata.atoms[index].map_number as u32,
                    source.atom_map.unwrap_or_default()
                );
            }
            for (bond, source) in state.graph.bonds.iter().zip(&molecule.bonds) {
                assert_eq!([bond.a as u32, bond.b as u32], source.atoms);
            }
            assert_eq!(
                identity(&molecule),
                crate::canonical_smiles(text, &Budget::default()).unwrap()
            );
        }
        for (stated, explicit) in [(false, 0), (true, 4)] {
            let mut methane = input("C");
            methane.atoms[0].hydrogens_stated = stated;
            let state = prepare(&methane, &Budget::default()).unwrap();
            assert_eq!(state.graph.atoms[0].explicit_hydrogens, explicit);
            assert_eq!(state.graph.atoms[0].no_implicit, stated);
        }
        // Pyrrolic H is inferred by the name parser, but cannot be inferred
        // from native aromatic N alone. Its count remains frozen and exact.
        let mut pyrrole = input("c1cc[nH]c1");
        pyrrole.atoms[3].hydrogens_stated = false;
        let state = prepare(&pyrrole, &Budget::default()).unwrap();
        assert_eq!(state.graph.atoms[3].explicit_hydrogens, 1);
        assert!(state.graph.atoms[3].no_implicit);
    }

    #[test]
    fn ordinary_reserved_maps_do_not_assert_ir_stereochemistry() {
        for map in [
            0,
            900_001,
            900_002,
            900_003,
            900_004,
            900_008,
            2_147_483_639,
        ] {
            let mut molecule = input("FC(Cl)(Br)I");
            molecule.atoms[1].atom_map = Some(map);
            let state = prepare(&molecule, &Budget::default()).unwrap();
            assert_eq!(state.metadata.atoms[1].map_number as u32, map);
            assert_eq!(state.metadata.atoms[1].chiral_tag, 0);
            assert!(!state.properties.atoms[1].unknown);
            assert!(state.metadata.bonds.iter().all(|meta| meta.stereo == 0));
            assert_eq!(identity(&molecule), identity(&input("FC(Cl)(Br)I")));
        }
        // Values outside native import's range are refused, never truncated.
        for map in [2_147_483_640, i32::MAX as u32, u32::MAX] {
            let text = format!("[CH4:{map}]");
            assert!(nazuchem_kernel::smiles::parse(&text).is_err());
            let mut molecule = input("C");
            molecule.atoms[0].atom_map = Some(map);
            assert_eq!(
                prepare(&molecule, &Budget::default()).unwrap_err(),
                "SMILES exceeds parser resource limits"
            );
        }
    }

    #[test]
    fn typed_stereo_conventions_match_full_cip_and_reference_permutations() {
        let parser = nazuchem_n2s::Parser::new().unwrap();
        for (name, expected) in [
            ("(2R)-butan-2-ol", "R"),
            ("(2S)-butan-2-ol", "S"),
            ("(R)-bromochlorofluoroiodomethane", "R"),
            ("(S)-bromochlorofluoroiodomethane", "S"),
            ("(E)-but-2-ene", "E"),
            ("(Z)-but-2-ene", "Z"),
        ] {
            let parsed = parser.parse(name, &nazuchem_n2s::ParseOptions::strict());
            let source = parsed
                .structure
                .unwrap_or_else(|| panic!("{name}: {}", parsed.message));
            let expected_name = serde_json::to_value(
                crate::name_smiles(&source.semantic_cxsmiles().unwrap(), &Budget::default())
                    .unwrap(),
            )
            .unwrap();
            let molecule = source.to_ir().unwrap();
            let state = label::assign(
                &prepare(&molecule, &Budget::default()).unwrap(),
                &Default::default(),
            )
            .unwrap()
            .state;
            let labels = state
                .properties
                .atoms
                .iter()
                .filter_map(|atom| atom.cip_code.as_deref())
                .chain(
                    state
                        .properties
                        .bond_codes
                        .iter()
                        .filter_map(|code| code.as_deref()),
                )
                .collect::<Vec<_>>();
            assert_eq!(labels, [expected], "{name}");
            let mut permuted = molecule.clone();
            for entry in &mut permuted.stereo {
                match entry {
                    Stereo::Tetrahedral {
                        neighbors, winding, ..
                    } => {
                        neighbors.swap(0, 1);
                        *winding = if *winding == Winding::Clockwise {
                            Winding::Anticlockwise
                        } else {
                            Winding::Clockwise
                        };
                    }
                    Stereo::DoubleBond { bond, refs, .. } => {
                        permuted.bonds[*bond as usize].atoms.swap(0, 1);
                        refs.swap(0, 1);
                    }
                    _ => panic!("ordinary named stereochemistry"),
                }
            }
            assert_eq!(identity(&molecule), identity(&permuted), "{name}");
            assert_eq!(complete_name(&molecule), expected_name, "{name}: original");
            assert_eq!(complete_name(&permuted), expected_name, "{name}: permuted");
        }
        for (marker, descriptor) in [(900_004, "M"), (900_008, "P")] {
            let text = format!("Oc1ccc2ccccc2[c:{marker}]1-[c:{marker}]1c(O)ccc2ccccc12");
            let molecule = input(&text);
            let expected_name =
                serde_json::to_value(crate::name_smiles(&text, &Budget::default()).unwrap())
                    .unwrap();
            let labeled = label::assign(
                &prepare(&molecule, &Budget::default()).unwrap(),
                &Default::default(),
            )
            .unwrap()
            .state;
            assert_eq!(
                labeled
                    .properties
                    .bond_codes
                    .iter()
                    .filter_map(|code| code.as_deref())
                    .collect::<Vec<_>>(),
                [descriptor]
            );
            let mut permuted = molecule.clone();
            for entry in &mut permuted.stereo {
                if let Stereo::Atropisomer { refs, winding, .. } = entry {
                    refs.swap(0, 1);
                    *winding = if *winding == AxialWinding::Clockwise {
                        AxialWinding::Anticlockwise
                    } else {
                        AxialWinding::Clockwise
                    };
                }
            }
            assert_eq!(identity(&molecule), identity(&permuted));
            let mut reversed = molecule.clone();
            for entry in &mut reversed.stereo {
                if let Stereo::Atropisomer { bond, refs, .. } = entry {
                    reversed.bonds[*bond as usize].atoms.swap(0, 1);
                    refs.rotate_left(2);
                }
            }
            assert_eq!(identity(&molecule), identity(&reversed));
            assert_eq!(
                complete_name(&molecule),
                expected_name,
                "{descriptor}: original"
            );
            assert_eq!(
                complete_name(&permuted),
                expected_name,
                "{descriptor}: permuted"
            );
            assert_eq!(
                complete_name(&reversed),
                expected_name,
                "{descriptor}: reversed"
            );
        }
    }

    #[test]
    fn named_ir_retains_hydrogen_vertex_identities_while_naming_recovers_the_same_parent() {
        let parser = nazuchem_n2s::Parser::new().unwrap();
        for name in ["(R)-butan-2-ol", "D-(+)-glucose"] {
            let parsed = parser.parse(name, &nazuchem_n2s::ParseOptions::strict());
            let source = parsed
                .structure
                .unwrap_or_else(|| panic!("{name}: {}", parsed.message));
            let molecule = source.to_ir().unwrap();
            let vertices = molecule
                .atoms
                .iter()
                .enumerate()
                .filter_map(|(index, atom)| {
                    (atom.element == 1 && atom.isotope.is_none()).then_some(index)
                })
                .collect::<Vec<_>>();
            assert!(!vertices.is_empty(), "{name}: expected protected stereo H");
            let state = prepare(&molecule, &Budget::default()).unwrap();
            assert_eq!(state.graph.atoms.len(), molecule.atoms.len());
            assert_eq!(state.graph.bonds.len(), molecule.bonds.len());
            for &hydrogen in &vertices {
                assert_eq!(state.graph.atoms[hydrogen].atomic_number, 1);
            }
            let text = source.semantic_cxsmiles().unwrap();
            let expected =
                serde_json::to_value(crate::name_smiles(&text, &Budget::default()).unwrap())
                    .unwrap();
            assert_eq!(complete_name(&molecule), expected, "{name}: original");
            for target in [0, molecule.atoms.len() / 2] {
                // Rotate the protected H into the front or middle so folding
                // leaves sparse naming IDs, rather than only removing a tail.
                let count = molecule.atoms.len();
                let shift = (target + count - vertices[0]) % count;
                let remap_atom = |atom: u32| ((atom as usize + shift) % count) as u32;
                let remap_ref = |reference: &mut Ref| {
                    if let Ref::Atom(atom) = reference {
                        *atom = remap_atom(*atom);
                    }
                };
                let mut rotated = molecule.clone();
                rotated.atoms.rotate_right(shift);
                for bond in &mut rotated.bonds {
                    bond.atoms = bond.atoms.map(remap_atom);
                }
                for entry in &mut rotated.stereo {
                    match entry {
                        Stereo::Tetrahedral {
                            center, neighbors, ..
                        } => {
                            *center = remap_atom(*center);
                            neighbors.iter_mut().for_each(remap_ref);
                        }
                        Stereo::DoubleBond { refs, .. } => *refs = refs.map(remap_atom),
                        Stereo::Atropisomer { refs, .. } => refs.iter_mut().for_each(remap_ref),
                        _ => panic!("test stereo variant"),
                    }
                }
                rotated
                    .parent
                    .iter_mut()
                    .for_each(|atom| *atom = remap_atom(*atom));
                rotated.validate().unwrap();
                assert_eq!(rotated.atoms[target].element, 1);
                let prepared = prepare(&rotated, &Budget::default()).unwrap();
                for &hydrogen in &vertices {
                    assert_eq!(
                        prepared.graph.atoms[remap_atom(hydrogen as u32) as usize].atomic_number,
                        1
                    );
                }
                let graph = crate::graph_io::read_molecule(&rotated, &Budget::default()).unwrap();
                let heavy_ids = rotated
                    .atoms
                    .iter()
                    .enumerate()
                    .filter_map(|(index, atom)| (atom.element != 1).then_some(index))
                    .collect::<Vec<_>>();
                assert_eq!(graph.atoms.keys().copied().collect::<Vec<_>>(), heavy_ids);
                assert_eq!(
                    complete_name(&rotated),
                    expected,
                    "{name}: hydrogen at {target}"
                );
            }
        }
    }

    #[test]
    fn authoritative_double_bonds_survive_without_direction_carriers() {
        let molecule = input("F/C(Cl)=C(/Br)I");
        let expected = identity(&molecule);
        let mut permuted = molecule.clone();
        let Stereo::DoubleBond {
            bond,
            refs,
            relation,
        } = &mut permuted.stereo[0]
        else {
            panic!("double bond")
        };
        let end = permuted.bonds[*bond as usize].atoms[0];
        refs[0] = permuted
            .bonds
            .iter()
            .flat_map(|bond| {
                [
                    (bond.atoms[0], bond.atoms[1]),
                    (bond.atoms[1], bond.atoms[0]),
                ]
            })
            .find_map(|(a, b)| {
                (a == end && b != refs[0] && b != permuted.bonds[*bond as usize].atoms[1])
                    .then_some(b)
            })
            .unwrap();
        *relation = if *relation == Relation::Cis {
            Relation::Trans
        } else {
            Relation::Cis
        };
        let state = prepare(&permuted, &Budget::default()).unwrap();
        assert!(
            state
                .directions
                .iter()
                .all(|direction| *direction == Direction::None)
        );
        assert_eq!(identity(&permuted), expected);
        let mut reversed = molecule.clone();
        let Stereo::DoubleBond { bond, refs, .. } = &mut reversed.stereo[0] else {
            panic!("double bond")
        };
        reversed.bonds[*bond as usize].atoms.swap(0, 1);
        refs.swap(0, 1);
        assert_eq!(identity(&reversed), expected);
    }

    #[test]
    fn partial_conjugated_configuration_needs_no_cx_projection() {
        let parser = nazuchem_n2s::Parser::new().unwrap();
        for name in [
            "(2E,6E)-octa-2,4,6-triene",
            "(2E,8E)-deca-2,4,5,6,8-pentaene",
        ] {
            let parsed = parser.parse(name, &nazuchem_n2s::ParseOptions::strict());
            let molecule = parsed
                .structure
                .unwrap_or_else(|| panic!("{name}: {}", parsed.message))
                .to_ir()
                .unwrap();
            assert!(matches!(
                molecule.to_cxsmiles(),
                Err(nazuchem_ir::CxsmilesError::StereoNotExpressible)
            ));
            let state = prepare(&molecule, &Budget::default()).unwrap();
            let specified = molecule
                .stereo
                .iter()
                .filter_map(|entry| match entry {
                    Stereo::DoubleBond { bond, relation, .. } if *relation != Relation::Unknown => {
                        Some(*bond as usize)
                    }
                    _ => None,
                })
                .collect::<Vec<_>>();
            assert_eq!(specified.len(), 2);
            for (index, bond) in state
                .graph
                .bonds
                .iter()
                .enumerate()
                .filter(|(_, bond)| bond.order == 2)
            {
                assert_eq!(
                    state.metadata.bonds[index].stereo != 0,
                    specified.contains(&index),
                    "{name}: bond {}–{}",
                    bond.a,
                    bond.b
                );
            }
        }
    }

    #[test]
    fn unknown_units_and_enhanced_groups_retain_structural_meaning() {
        for text in ["F[C@](Cl)(Br)I", "F/C=C/F"] {
            let mut molecule = input(text);
            match &mut molecule.stereo[0] {
                Stereo::Tetrahedral { winding, .. } => *winding = Winding::Unknown,
                Stereo::DoubleBond { relation, .. } => *relation = Relation::Unknown,
                _ => panic!("ordinary stereo"),
            }
            let state = prepare(&molecule, &Budget::default()).unwrap();
            match molecule.stereo[0] {
                Stereo::Tetrahedral { center, .. } => {
                    assert!(state.properties.atoms[center as usize].unknown);
                    assert_eq!(state.metadata.atoms[center as usize].chiral_tag, 0);
                }
                Stereo::DoubleBond { bond, .. } => {
                    assert_eq!(state.metadata.bonds[bond as usize].stereo, 1);
                    assert!(state.metadata.bonds[bond as usize].unknown_stereo);
                }
                _ => unreachable!(),
            }
            let named = crate::name_molecule(
                &molecule,
                nazuchem_ir::IsotopeStyle::default(),
                &Budget::default(),
            )
            .unwrap();
            assert_eq!(
                named.status,
                NameStatus::NotPreferred(NameReason::UnknownConfiguration)
            );
            assert_eq!(named.omitted_units, 1);
        }
        let mut axis = input("Oc1ccc2ccccc2[c:900004]1-[c:900004]1c(O)ccc2ccccc12");
        let Stereo::Atropisomer { bond, winding, .. } = &mut axis.stereo[0] else {
            panic!("axis")
        };
        *winding = AxialWinding::Unknown;
        let bond = *bond as usize;
        let state = prepare(&axis, &Budget::default()).unwrap();
        assert_eq!(state.metadata.bonds[bond].stereo, 1);
        assert!(state.metadata.bonds[bond].unknown_stereo);
        let named = crate::name_molecule(
            &axis,
            nazuchem_ir::IsotopeStyle::default(),
            &Budget::default(),
        )
        .unwrap();
        assert_eq!(
            named.status,
            NameStatus::NotPreferred(NameReason::UnknownConfiguration)
        );
        assert_eq!(named.omitted_units, 1);
        for text in [
            "F[C@](Cl)(Br)I",
            "F/C=C/F",
            "Oc1ccc2ccccc2[c:900004]1-[c:900004]1c(O)ccc2ccccc12",
        ] {
            for kind in [StereoGroupKind::And(19), StereoGroupKind::Or(23)] {
                let mut molecule = input(text);
                molecule.stereo_groups.push(IrGroup {
                    kind,
                    stereo: vec![0],
                });
                let state = prepare(&molecule, &Budget::default()).unwrap();
                let group = &state.metadata.groups[0];
                assert_eq!(
                    (group.kind, group.read_id, group.write_id),
                    match kind {
                        StereoGroupKind::And(id) => (2, id, id),
                        StereoGroupKind::Or(id) => (1, id, id),
                        _ => unreachable!(),
                    }
                );
                match molecule.stereo[0] {
                    Stereo::Tetrahedral { center, .. } => {
                        assert_eq!(group.atoms, [center as usize])
                    }
                    Stereo::DoubleBond { bond, .. } | Stereo::Atropisomer { bond, .. } => {
                        assert_eq!(group.bonds, [bond as usize])
                    }
                    _ => unreachable!(),
                }
                assert!(
                    crate::name_molecule(
                        &molecule,
                        nazuchem_ir::IsotopeStyle::default(),
                        &Budget::default()
                    )
                    .is_err()
                );
            }
        }
    }

    #[test]
    fn direct_admission_does_not_clean_unsupported_or_inconsistent_input() {
        let mut invalid = input("FCCl");
        invalid.atoms[1].hydrogens = 1;
        invalid.stereo.push(Stereo::Tetrahedral {
            center: 1,
            neighbors: [
                Ref::Atom(0),
                Ref::Atom(2),
                Ref::ImplicitHydrogen,
                Ref::LonePair,
            ],
            winding: Winding::Clockwise,
        });
        assert!(prepare(&invalid, &Budget::default()).is_err());
        let mut cumulene = input("CC=C=CC");
        cumulene.stereo.push(Stereo::DoubleBond {
            bond: 1,
            refs: [0, 3],
            relation: Relation::Trans,
        });
        assert_eq!(
            prepare(&cumulene, &Budget::default()).unwrap_err(),
            CUMULENE_STEREO
        );
        let mut invalid_radical = input("[CH3]");
        invalid_radical.atoms[0].radical_electrons = 0;
        assert!(
            prepare(&invalid_radical, &Budget::default())
                .unwrap_err()
                .contains("IR radical electron count")
        );
        let mut oversized = Molecule::default();
        oversized.atoms.resize(
            513,
            nazuchem_ir::Atom {
                element: 6,
                hydrogens: 4,
                ..Default::default()
            },
        );
        assert!(
            prepare(&oversized, &Budget::default())
                .unwrap_err()
                .contains("512 input atom tokens")
        );
        let budget = Budget::new(0, usize::MAX, None, usize::MAX, |_| Some(Box::new(())));
        assert!(
            prepare(&input("C"), &budget)
                .unwrap_err()
                .contains("estimated memory budget exceeded")
        );
    }

    #[test]
    fn unknown_axes_require_the_same_biaryl_topology_as_known_axes() {
        let mut molecule = input("FC(=O)C(=O)F");
        let axis = molecule
            .bonds
            .iter()
            .position(|bond| bond.atoms == [1, 3])
            .unwrap() as u32;
        molecule.stereo.push(Stereo::Atropisomer {
            bond: axis,
            refs: [Ref::Atom(0), Ref::Atom(2), Ref::Atom(4), Ref::Atom(5)],
            winding: AxialWinding::Unknown,
        });
        assert_eq!(
            prepare(&molecule, &Budget::default()).unwrap_err(),
            NONCLASSICAL_STEREO
        );
        // Maps that happen to match the legacy private transport have no role
        // in the typed IR eligibility check.
        for atom in [1, 3] {
            molecule.atoms[atom].atom_map = Some(900_004);
        }
        assert_eq!(
            prepare(&molecule, &Budget::default()).unwrap_err(),
            NONCLASSICAL_STEREO
        );
    }

    #[test]
    fn unknown_symmetric_ligands_keep_existing_marker_admission_and_name_results() {
        for text in [
            "C[CH:900001](C)O",
            "C[CH:900002]=[C:900002](C)C",
            "C[CH:900002]=[C:900002]=CC",
        ] {
            let molecule = input(text);
            let legacy = crate::name_smiles(text, &Budget::default()).unwrap();
            let direct = crate::name_molecule(
                &molecule,
                nazuchem_ir::IsotopeStyle::default(),
                &Budget::default(),
            )
            .unwrap();
            assert_eq!(
                serde_json::to_value(direct).unwrap(),
                serde_json::to_value(legacy).unwrap(),
                "{text}"
            );
        }
    }
}
