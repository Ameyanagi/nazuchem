// SPDX-License-Identifier: MIT
// Origin: original
// Source: Nazuchem | crates/nazuchem-n2s/tests/support/ir_stereo.rs | Phase B step 8a

//! Independent identity-based audits, before any kernel stereo cleanup.

use std::collections::{BTreeMap, BTreeSet};

use nazuchem_ir::{AxialWinding, Molecule, Ref, Relation, Stereo, Winding};
use nazuchem_kernel::{smiles, stereo::bond_stereo_from_directions};
use nazuchem_n2s::{
    Structure,
    graph::{AtomId, BondStereoValue, StereoGroupType, StereoReference},
};

#[derive(Debug, Default)]
pub struct Counts {
    pub tetrahedral: usize,
    pub double_bond: usize,
    pub atropisomer: usize,
}

/// Parity of an explicit reference permutation, independent of the writer.
fn odd_permutation<T: Eq>(from: &[T], to: &[T]) -> Result<bool, String> {
    if from.len() != to.len() {
        return Err("reference counts differ".into());
    }
    let permutation = from
        .iter()
        .map(|reference| {
            to.iter()
                .position(|other| other == reference)
                .ok_or("reference sets differ".to_owned())
        })
        .collect::<Result<Vec<_>, _>>()?;
    let mut odd = false;
    for (i, &a) in permutation.iter().enumerate() {
        for &b in &permutation[i + 1..] {
            odd ^= a > b;
        }
    }
    Ok(odd)
}

/// Compare every source assertion against IR under the exact source identity map.
/// Also check complete locants and principal-parent membership before H folding.
pub fn source(
    structure: &Structure,
    molecule: &Molecule,
    order: &[AtomId],
) -> Result<Counts, String> {
    let graph = &structure.graph;
    let map: BTreeMap<_, _> = order
        .iter()
        .enumerate()
        .map(|(i, atom)| (atom.0, i as u32))
        .collect();
    let index = |atom: AtomId| {
        map.get(&atom.0)
            .copied()
            .ok_or_else(|| format!("source identity {} was folded", atom.0))
    };
    for &atom in &graph.fragment(structure.fragment).atoms {
        if !graph.atom(atom).locants.is_empty() && !map.contains_key(&atom.0) {
            return Err(format!("locanted source atom {} was folded", atom.0));
        }
    }
    for (&source, atom) in order.iter().zip(&molecule.atoms) {
        if graph.atom(source).locants != atom.locants {
            return Err(format!("complete locants changed at {}", source.0));
        }
    }
    let expected_parent = structure
        .principal_parent_atoms
        .iter()
        .copied()
        .filter(|id| graph.fragment(structure.fragment).atoms.contains(id))
        .map(index)
        .collect::<Result<BTreeSet<_>, _>>()?;
    if expected_parent != molecule.parent.iter().copied().collect() {
        return Err("complete principal-parent membership changed".into());
    }
    let mut counts = Counts::default();
    for &id in &graph.fragment(structure.fragment).atoms {
        let Some(parity) = &graph.atom(id).parity else {
            continue;
        };
        let center = index(id)?;
        let expected = parity
            .atom_refs
            .iter()
            .map(|reference| match reference {
                Some(StereoReference::Atom(other)) if *other == id => Ok(Ref::LonePair),
                Some(StereoReference::Atom(other)) => index(*other).map(Ref::Atom),
                Some(StereoReference::ImplicitHydrogen) => Ok(Ref::ImplicitHydrogen),
                _ => Err("unresolved source parity reference".into()),
            })
            .collect::<Result<Vec<_>, String>>()?;
        let Some(Stereo::Tetrahedral {
            neighbors, winding, ..
        }) = molecule
            .stereo
            .iter()
            .find(|stereo| matches!(stereo, Stereo::Tetrahedral {center:c,..} if *c == center))
        else {
            return Err("source tetrahedral assertion omitted".into());
        };
        let odd = odd_permutation(&expected, neighbors)?;
        if parity.stereo_group.kind == StereoGroupType::Unknown {
            if *winding != Winding::Unknown {
                return Err("unknown source tetrahedral assertion became known".into());
            }
        } else {
            // OPSIN parity +1 describes clockwise rotation viewed from ref 0.
            let clockwise = (parity.parity == 1) ^ odd;
            if *winding
                != if clockwise {
                    Winding::Clockwise
                } else {
                    Winding::Anticlockwise
                }
            {
                return Err(format!("source tetrahedral winding changed at {center}"));
            }
        }
        counts.tetrahedral += 1;
    }
    for &id in &graph.fragment(structure.fragment).bonds {
        let bond = graph.bond(id);
        if bond.stereo.is_none() && bond.axial_stereo.is_none() {
            continue;
        }
        let ends = [index(bond.from)?, index(bond.to)?];
        let ir_bond = molecule
            .bonds
            .iter()
            .position(|b| b.atoms == ends || b.atoms == [ends[1], ends[0]])
            .ok_or("stereo bond omitted")?;
        let reversed = molecule.bonds[ir_bond].atoms != ends;
        if let Some(stereo) = &bond.stereo {
            let mut expected = [index(stereo.atom_refs[0])?, index(stereo.atom_refs[3])?];
            if stereo.atom_refs[1] != bond.from {
                expected.swap(0, 1);
            }
            if reversed {
                expected.swap(0, 1);
            }
            let Some(Stereo::DoubleBond { refs, relation, .. }) = molecule
                .stereo
                .iter()
                .find(|s| matches!(s, Stereo::DoubleBond {bond,..} if *bond as usize == ir_bond))
            else {
                return Err("source double-bond assertion omitted".into());
            };
            let cis = (stereo.value == BondStereoValue::Cis)
                ^ (expected[0] != refs[0])
                ^ (expected[1] != refs[1]);
            if *relation != if cis { Relation::Cis } else { Relation::Trans } {
                return Err("source double-bond configuration changed".into());
            }
            counts.double_bond += 1;
        }
        if let Some(axis) = &bond.axial_stereo {
            let mut expected = [
                Ref::Atom(index(axis.refs[0])?),
                Ref::Atom(index(axis.refs[1])?),
            ];
            if reversed {
                expected.swap(0, 1);
            }
            let Some(Stereo::Atropisomer { refs, winding, .. }) = molecule
                .stereo
                .iter()
                .find(|s| matches!(s, Stereo::Atropisomer {bond,..} if *bond as usize == ir_bond))
            else {
                return Err("source axis assertion omitted".into());
            };
            let clockwise = axis.clockwise ^ (expected[0] != refs[0]) ^ (expected[1] != refs[2]);
            if *winding
                != if clockwise {
                    AxialWinding::P
                } else {
                    AxialWinding::M
                }
            {
                return Err("source axis configuration changed".into());
            }
            counts.atropisomer += 1;
        }
    }
    if counts.tetrahedral + counts.double_bond + counts.atropisomer != molecule.stereo.len() {
        return Err("IR introduced a source-absent stereo assertion".into());
    }
    Ok(counts)
}

/// Reparse the emitted standard SMILES with the kernel, without perception that
/// could erase unsupported source assertions. Every atom and bond is checked,
/// including those with no IR stereo entry (not just the specified subset).
pub fn written(molecule: &Molecule, text: &str, order: &[u32]) -> Result<(), String> {
    let parsed = smiles::parse(super::cx::bare(text)).map_err(|e| e.to_string())?;
    let metadata =
        bond_stereo_from_directions(&parsed.graph, &parsed.metadata, &parsed.directions)?;
    // Parse the actual CX extension too, without preparation/cleanup. In
    // particular, w/ctu must describe only unknown units present in the IR.
    let mut unknown_atoms = BTreeSet::new();
    let mut unknown_bonds = BTreeSet::new();
    if let Some((_, extension)) = text.split_once(" |") {
        use nazuchem_kernel::cx::{self, Event};
        let topology = cx::Topology {
            atoms: parsed.graph.atoms.len(),
            bonds: parsed
                .graph
                .bonds
                .iter()
                .zip(&parsed.bond_indices)
                .map(|(bond, &index)| cx::ParseBond {
                    a: bond.a,
                    b: bond.b,
                    index: Some(index),
                })
                .collect(),
        };
        let events =
            cx::read(&format!("|{extension}"), &topology).map_err(|error| error.to_string())?;
        for event in events.events {
            match event {
                Event::Wedge { bond, begin, kind } => {
                    if kind != 2 {
                        return Err("unexpected CX wedge assertion".into());
                    }
                    let center = order[begin];
                    let ends = [
                        order[parsed.graph.bonds[bond].a],
                        order[parsed.graph.bonds[bond].b],
                    ];
                    let expected = molecule.stereo.iter().any(|unit| match unit {
                        Stereo::Tetrahedral {
                            center: c,
                            winding: Winding::Unknown,
                            ..
                        } => *c == center,
                        Stereo::DoubleBond {
                            bond,
                            relation: Relation::Unknown,
                            ..
                        } => {
                            let actual = molecule.bonds[*bond as usize].atoms;
                            actual == ends || actual == [ends[1], ends[0]]
                        }
                        _ => false,
                    });
                    if !expected {
                        return Err("CX introduced unknown stereo outside IR".into());
                    }
                    unknown_atoms.insert(center);
                }
                Event::BondStereo { bond, stereo } => {
                    let ends = [
                        order[parsed.graph.bonds[bond].a],
                        order[parsed.graph.bonds[bond].b],
                    ];
                    let expected = molecule.stereo.iter().any(|unit| matches!(unit, Stereo::DoubleBond {bond, relation:Relation::Unknown, ..} if { let actual = molecule.bonds[*bond as usize].atoms; actual == ends || actual == [ends[1], ends[0]] }));
                    if stereo != 1 || !expected {
                        return Err("CX introduced a bond configuration outside IR".into());
                    }
                    unknown_bonds.insert([ends[0].min(ends[1]), ends[0].max(ends[1])]);
                }
                _ => {}
            }
        }
    }
    for unit in &molecule.stereo {
        match unit {
            Stereo::Tetrahedral {
                center,
                winding: Winding::Unknown,
                ..
            } if !unknown_atoms.contains(center) => {
                return Err("writer dropped unknown tetrahedral state".into());
            }
            Stereo::DoubleBond {
                bond,
                relation: Relation::Unknown,
                ..
            } => {
                let [a, b] = molecule.bonds[*bond as usize].atoms;
                if !unknown_bonds.contains(&[a.min(b), a.max(b)]) {
                    return Err("writer dropped unknown double-bond state".into());
                }
            }
            _ => {}
        }
    }
    for (atom, meta) in metadata.atoms.iter().enumerate() {
        let center = order[atom];
        let expected = molecule.stereo.iter().find_map(|s| match s {
            Stereo::Tetrahedral {
                center: c,
                neighbors,
                winding,
            } if *c == center && *winding != Winding::Unknown => Some((neighbors, winding)),
            _ => None,
        });
        let Some((expected, winding)) = expected else {
            if meta.chiral_tag != 0 {
                return Err(format!("writer introduced tetrahedral stereo at {center}"));
            }
            continue;
        };
        if !matches!(meta.chiral_tag, 1 | 2) {
            return Err(format!("writer dropped tetrahedral stereo at {center}"));
        }
        let mut neighbors = parsed
            .graph
            .bonds
            .iter()
            .filter_map(|b| {
                if b.a == atom {
                    Some(Ref::Atom(order[b.b]))
                } else if b.b == atom {
                    Some(Ref::Atom(order[b.a]))
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();
        if neighbors.len() == 3 {
            neighbors.push(if expected.contains(&Ref::ImplicitHydrogen) {
                Ref::ImplicitHydrogen
            } else {
                Ref::LonePair
            });
        }
        let clockwise = (*winding == Winding::Clockwise) ^ odd_permutation(expected, &neighbors)?;
        if (meta.chiral_tag == 1) != clockwise {
            return Err(format!(
                "writer changed tetrahedral configuration at {center}"
            ));
        }
    }
    for (id, bond) in parsed.graph.bonds.iter().enumerate() {
        let ends = [order[bond.a], order[bond.b]];
        let ir_bond = molecule
            .bonds
            .iter()
            .position(|b| b.atoms == ends || b.atoms == [ends[1], ends[0]])
            .ok_or("writer introduced bond")?;
        let expected = molecule.stereo.iter().find_map(|s| match s {
            Stereo::DoubleBond {
                bond,
                refs,
                relation,
            } if *bond as usize == ir_bond && *relation != Relation::Unknown => {
                Some((refs, relation))
            }
            _ => None,
        });
        let meta = &metadata.bonds[id];
        let Some((refs, relation)) = expected else {
            if meta.stereo != 0 {
                return Err(format!("writer introduced bond stereo at {ir_bond}"));
            }
            continue;
        };
        if !matches!(meta.stereo, 4 | 5) || meta.stereo_atoms.len() != 2 {
            return Err(format!("writer dropped double-bond stereo at {ir_bond}"));
        }
        let mut actual = [order[meta.stereo_atoms[0]], order[meta.stereo_atoms[1]]];
        if molecule.bonds[ir_bond].atoms != ends {
            actual.swap(0, 1);
        }
        let cis = (meta.stereo == 4) ^ (actual[0] != refs[0]) ^ (actual[1] != refs[1]);
        if cis != (*relation == Relation::Cis) {
            return Err(format!(
                "writer changed double-bond configuration at {ir_bond}"
            ));
        }
    }
    if molecule.stereo.iter().any(
        |s| matches!(s, Stereo::Atropisomer { winding, .. } if *winding != AxialWinding::Unknown),
    ) {
        return Err("known axis cannot be omitted by standard projection".into());
    }
    Ok(())
}
