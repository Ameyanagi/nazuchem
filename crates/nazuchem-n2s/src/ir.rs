// SPDX-License-Identifier: MIT
// Origin: original
// Source: Nazuchem | crates/nazuchem-n2s/src/ir.rs | Phase B step 8a

//! Projection of the finalized OPSIN graph into the versioned data contract.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::fmt;

use nazuchem_ir::{
    Atom, AxialWinding, Bond, BondOrder, Connectivity, Molecule, PositionVariation, Ref, Relation,
    RepeatUnit, Stereo, StereoGroup, StereoGroupKind, Winding,
};

use crate::graph::{AtomId, BondStereoValue, Element, StereoGroupType, StereoReference};
use crate::{SerializationError, Structure};

/// A graph cannot be projected without discarding information.
#[non_exhaustive]
#[derive(Debug)]
pub enum IrError {
    /// The existing semantic exporter rejected the construction graph.
    Serialization(SerializationError),
    /// Unfinished construction state or a value outside the IR's numeric range.
    InvalidGraph(String),
    /// The projected graph failed the IR's structural validation.
    Validation(nazuchem_ir::ValidationError),
}

impl fmt::Display for IrError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Serialization(error) => error.fmt(f),
            Self::InvalidGraph(reason) => f.write_str(reason),
            Self::Validation(error) => error.fmt(f),
        }
    }
}

impl std::error::Error for IrError {}

fn invalid(reason: impl Into<String>) -> IrError {
    IrError::InvalidGraph(reason.into())
}

impl Structure {
    /// Convert a finalized graph without chemistry perception or a SMILES parse.
    ///
    /// Atom identity follows the semantic CXSMILES writer, followed by hydrogen
    /// vertices whose metadata or stereochemical role requires an identity.
    /// Only ordinary hydrogen vertices without metadata become counts.
    pub fn to_ir(&self) -> Result<Molecule, IrError> {
        self.to_ir_with_atom_mapping().map(|(molecule, _)| molecule)
    }

    /// Project the graph and return the source atom at each IR atom identity.
    ///
    /// The mapping is independent of any subsequent CXSMILES traversal. It also
    /// allows corpus gates to audit source metadata and stereo before hydrogen
    /// folding, without reconstructing source identities from output strings.
    #[doc(hidden)]
    pub fn to_ir_with_atom_mapping(&self) -> Result<(Molecule, Vec<AtomId>), IrError> {
        let graph = &self.graph;
        let fragment = graph
            .fragments
            .get(self.fragment.0)
            .filter(|fragment| fragment.active)
            .ok_or_else(|| invalid("Missing or inactive structure fragment"))?;
        let mut order = crate::smiles::semantic_atom_order(graph, self.fragment)
            .map_err(IrError::Serialization)?;
        let mut retained: HashSet<_> = order.iter().copied().collect();
        let mut protected: HashSet<_> = self.principal_parent_atoms.iter().copied().collect();
        for &id in &fragment.atoms {
            let atom = graph.atom(id);
            if atom.parity.is_some() {
                protected.insert(id);
            }
            if let Some(parity) = &atom.parity {
                protected.extend(
                    parity
                        .atom_refs
                        .iter()
                        .filter_map(|reference| match reference {
                            Some(StereoReference::Atom(id)) => Some(*id),
                            _ => None,
                        }),
                );
            }
            if let Some(endpoints) = &atom.properties.position_variation_bond {
                protected.insert(id);
                protected.extend(endpoints.iter().copied());
            }
        }
        for &id in &fragment.bonds {
            let bond = graph.bond(id);
            if let Some(stereo) = &bond.stereo {
                protected.extend(stereo.atom_refs);
            }
            if let Some(stereo) = &bond.axial_stereo {
                protected.extend([bond.from, bond.to]);
                protected.extend(stereo.refs);
                protected.extend(graph.neighbours(bond.from));
                protected.extend(graph.neighbours(bond.to));
            }
        }
        // Stated unknown configurations keep the same identities as defined ones.
        for &id in &self.unknown_stereo_atoms {
            protected.insert(id);
            protected.extend(graph.neighbours(id));
        }
        for &id in &self.unknown_stereo_bonds {
            let bond = graph.bond(id);
            protected.extend(graph.neighbours(bond.from));
            protected.extend(graph.neighbours(bond.to));
        }
        if let Some(points) = &fragment.polymer_attachment_points {
            protected.extend(points.iter().copied());
        }
        // Keep the legacy heavy-atom traversal, but never let its display-only
        // H folding discard source identity or structural/name metadata.
        for &id in &fragment.atoms {
            if retained.contains(&id) {
                continue;
            }
            let atom = graph.atom(id);
            let foldable = atom.element == Element::H
                && atom.isotope.is_none()
                && atom.charge == 0
                && atom.locants.is_empty()
                && atom.properties.homology_group.is_none()
                && atom.properties.atom_class.is_none()
                && atom.properties.smiles_hydrogen_count.is_none()
                && atom.explicit_hydrogens == 0
                && atom.bonds.len() == 1
                && graph.bond(atom.bonds[0]).order == 1
                && !protected.contains(&id);
            if !foldable {
                retained.insert(id);
                order.push(id);
            }
        }
        let indices: HashMap<_, _> = order
            .iter()
            .enumerate()
            .map(|(index, &atom)| {
                u32::try_from(index)
                    .map(|index| (atom, index))
                    .map_err(|_| invalid("Too many atoms for IR v1"))
            })
            .collect::<Result<_, _>>()?;
        let index = |atom: AtomId| {
            indices
                .get(&atom)
                .copied()
                .ok_or_else(|| invalid("Reference to an unwritten atom"))
        };
        let polymer_points: HashSet<_> = fragment
            .polymer_attachment_points
            .as_ref()
            .into_iter()
            .flatten()
            .copied()
            .collect();
        let polymer = !polymer_points.is_empty();
        let mut used_classes = HashSet::new();
        let mut next_class = 1;
        let mut atoms = Vec::with_capacity(order.len());
        for &id in &order {
            let atom = graph.atom(id);
            let folded = graph
                .neighbours(id)
                .into_iter()
                .filter(|other| !indices.contains_key(other))
                .count();
            let hydrogens = u8::try_from(folded as u64 + u64::from(atom.explicit_hydrogens))
                .map_err(|_| invalid("Hydrogen count exceeds IR v1 range"))?;
            let attachment_point = if atom.element == Element::R
                && !polymer_points.contains(&id)
                && atom.properties.homology_group.is_none()
            {
                let class = if let Some(class) = atom.properties.atom_class {
                    used_classes.insert(class);
                    class
                } else {
                    while used_classes.contains(&next_class) {
                        next_class += 1;
                    }
                    let class = next_class;
                    next_class += 1;
                    class
                };
                Some(class)
            } else {
                None
            };
            atoms.push(Atom {
                element: atom.element.atomic_number(),
                isotope: atom
                    .isotope
                    .map(u16::try_from)
                    .transpose()
                    .map_err(|_| invalid("Isotope exceeds IR v1 range"))?,
                charge: i8::try_from(atom.charge)
                    .map_err(|_| invalid("Charge exceeds IR v1 range"))?,
                hydrogens,
                hydrogens_stated: atom.properties.smiles_hydrogen_count.is_some(),
                // Out valency tracks construction attachment sites. Elementary
                // atoms can also have unpaired electrons without any out atom.
                radical_electrons: 0,
                atom_map: attachment_point
                    .is_none()
                    .then_some(atom.properties.atom_class)
                    .flatten(),
                attachment_point,
                // Homology/R-group names use the literal standard CX `_p`
                // label suffix. Attachment and polymer roles stay typed.
                label: atom
                    .properties
                    .homology_group
                    .as_ref()
                    .map(|group| format!("{group}_p")),
                locants: atom.locants.clone(),
            });
        }
        let mut bonds = Vec::new();
        let mut bond_indices = HashMap::new();
        for &id in &fragment.bonds {
            let bond = graph.bond(id);
            let (Some(&from), Some(&to)) = (indices.get(&bond.from), indices.get(&bond.to)) else {
                continue; // The writer folded an ordinary hydrogen vertex.
            };
            let order = match bond.order {
                1 => BondOrder::Single,
                2 => BondOrder::Double,
                3 => BondOrder::Triple,
                4 => BondOrder::Quadruple,
                _ => return Err(invalid("Non-Kekule bond order")),
            };
            bond_indices.insert(id, bonds.len() as u32);
            bonds.push(Bond {
                atoms: [from, to],
                order,
            });
        }
        assign_radical_electrons(&mut atoms, &bonds)?;
        let mut stereo = Vec::new();
        let mut groups = BTreeMap::<_, Vec<u32>>::new();
        for &id in &order {
            let atom = graph.atom(id);
            let Some(parity) = &atom.parity else {
                continue;
            };
            let mut neighbors = [Ref::LonePair; 4];
            for (target, reference) in neighbors.iter_mut().zip(parity.atom_refs) {
                *target = match reference {
                    Some(StereoReference::Atom(other)) if other == id => Ref::LonePair,
                    Some(StereoReference::Atom(other)) => match indices.get(&other) {
                        Some(&index) => Ref::Atom(index),
                        None if graph.atom(other).element == Element::H => Ref::ImplicitHydrogen,
                        None => return Err(invalid("Unwritten tetrahedral reference")),
                    },
                    Some(StereoReference::ImplicitHydrogen) => Ref::ImplicitHydrogen,
                    _ => return Err(invalid("Unresolved tetrahedral reference")),
                };
            }
            let winding = match parity.parity {
                _ if parity.stereo_group.kind == StereoGroupType::Unknown => Winding::Unknown,
                1 => Winding::Clockwise,
                -1 => Winding::Anticlockwise,
                _ => return Err(invalid("Tetrahedral winding must be +1 or -1")),
            };
            let stereo_index = stereo.len() as u32;
            stereo.push(Stereo::Tetrahedral {
                center: index(id)?,
                neighbors,
                winding,
            });
            if !matches!(
                parity.stereo_group.kind,
                StereoGroupType::Unknown | StereoGroupType::Absolute
            ) {
                let kind = match parity.stereo_group.kind {
                    StereoGroupType::Racemic => StereoGroupKind::And(parity.stereo_group.number),
                    StereoGroupType::Relative => StereoGroupKind::Or(parity.stereo_group.number),
                    StereoGroupType::Unknown | StereoGroupType::Absolute => unreachable!(),
                };
                groups.entry(kind).or_default().push(stereo_index);
            }
        }
        for &id in &self.unknown_stereo_atoms {
            let mut neighbors = [Ref::LonePair; 4];
            let atoms = graph.neighbours(id);
            if !(3..=4).contains(&atoms.len()) {
                return Err(invalid("Unknown stereocentre is not tetrahedral"));
            }
            for (target, atom) in neighbors.iter_mut().zip(atoms) {
                *target = Ref::Atom(index(atom)?);
            }
            stereo.push(Stereo::Tetrahedral {
                center: index(id)?,
                neighbors,
                winding: Winding::Unknown,
            });
        }
        for &id in &fragment.bonds {
            let bond = graph.bond(id);
            let Some(configuration) = &bond.stereo else {
                continue;
            };
            let refs = configuration.atom_refs;
            if bond.other_atom(refs[1]) != Some(refs[2]) {
                return Err(invalid("IR v1 cannot represent odd-cumulene E/Z"));
            }
            let refs = if refs[1] == bond.from && refs[2] == bond.to {
                [index(refs[0])?, index(refs[3])?]
            } else if refs[2] == bond.from && refs[1] == bond.to {
                [index(refs[3])?, index(refs[0])?]
            } else {
                return Err(invalid("Double-bond stereo endpoints are inconsistent"));
            };
            stereo.push(Stereo::DoubleBond {
                bond: *bond_indices
                    .get(&id)
                    .ok_or_else(|| invalid("Folded stereo bond"))?,
                refs,
                relation: match configuration.value {
                    BondStereoValue::Cis => Relation::Cis,
                    BondStereoValue::Trans => Relation::Trans,
                },
            });
        }
        for &id in &self.unknown_stereo_bonds {
            let bond = graph.bond(id);
            let substituent = |end: AtomId, partner: AtomId| {
                graph
                    .neighbours(end)
                    .into_iter()
                    .find(|&other| other != partner)
                    .ok_or_else(|| invalid("Unknown double bond end has no substituent"))
                    .and_then(index)
            };
            stereo.push(Stereo::DoubleBond {
                bond: *bond_indices
                    .get(&id)
                    .ok_or_else(|| invalid("Folded stereo bond"))?,
                refs: [
                    substituent(bond.from, bond.to)?,
                    substituent(bond.to, bond.from)?,
                ],
                relation: Relation::Unknown,
            });
        }
        for &id in &fragment.bonds {
            let bond = graph.bond(id);
            let Some(configuration) = &bond.axial_stereo else {
                continue;
            };
            let mut refs = [Ref::ImplicitHydrogen; 4];
            for (end, (atom, partner, reference)) in [
                (bond.from, bond.to, configuration.refs[0]),
                (bond.to, bond.from, configuration.refs[1]),
            ]
            .into_iter()
            .enumerate()
            {
                let neighbors: Vec<_> = graph
                    .neighbours(atom)
                    .into_iter()
                    .filter(|&other| other != partner)
                    .collect();
                if !neighbors.contains(&reference) || !(1..=2).contains(&neighbors.len()) {
                    return Err(invalid(
                        "Axial stereo references do not match axis neighbours",
                    ));
                }
                refs[end * 2] = Ref::Atom(index(reference)?);
                if let Some(other) = neighbors.into_iter().find(|&other| other != reference) {
                    refs[end * 2 + 1] = Ref::Atom(index(other)?);
                } else if graph.atom(atom).explicit_hydrogens == 0 {
                    return Err(invalid("Axial endpoint has no second ligand"));
                }
            }
            stereo.push(Stereo::Atropisomer {
                bond: *bond_indices
                    .get(&id)
                    .ok_or_else(|| invalid("Folded axial bond"))?,
                refs,
                // Source refs[0] and refs[1] become IR refs[0] and refs[2],
                // with the same bond direction. Preserve their relative turn;
                // the requested CIP descriptor does not label this winding.
                winding: if configuration.clockwise {
                    AxialWinding::Clockwise
                } else {
                    AxialWinding::Anticlockwise
                },
            });
        }
        let stereo_groups = groups
            .into_iter()
            .map(|(kind, stereo)| StereoGroup { kind, stereo })
            .collect();
        let mut position_variation = Vec::new();
        for &id in &order {
            if let Some(endpoints) = &graph.atom(id).properties.position_variation_bond {
                let candidates: Vec<_> = graph
                    .atom(id)
                    .bonds
                    .iter()
                    .filter_map(|bond| bond_indices.get(bond).copied())
                    .collect();
                if candidates.len() != 1 {
                    return Err(invalid("Position variation requires one attachment bond"));
                }
                let bond_index = candidates[0];
                let attachment = index(id)?;
                let bond = &mut bonds[bond_index as usize];
                if bond.atoms[0] != attachment {
                    if bond.atoms[1] != attachment {
                        return Err(invalid(
                            "Position-variation bond is not incident to its attachment",
                        ));
                    }
                    // IR v1 identifies the source of m: through the first bond
                    // endpoint. Keep double-bond substituents in the same
                    // corresponding-end order when orienting this bond.
                    bond.atoms.swap(0, 1);
                    for entry in &mut stereo {
                        match entry {
                            Stereo::DoubleBond { bond, refs, .. } if *bond == bond_index => {
                                refs.swap(0, 1)
                            }
                            Stereo::Atropisomer { bond, refs, .. } if *bond == bond_index => {
                                refs.rotate_left(2)
                            }
                            _ => {}
                        }
                    }
                }
                position_variation.push(PositionVariation {
                    bond: bond_index,
                    endpoints: endpoints
                        .iter()
                        .map(|&atom| index(atom))
                        .collect::<Result<_, _>>()?,
                });
            }
        }
        let polymers = if polymer {
            let atoms: Vec<_> = order
                .iter()
                .enumerate()
                .filter_map(|(i, id)| (!polymer_points.contains(id)).then_some(i as u32))
                .collect();
            let included: HashSet<_> = atoms.iter().copied().collect();
            let crossing_bonds = bonds
                .iter()
                .enumerate()
                .filter_map(|(i, bond)| {
                    (included.contains(&bond.atoms[0]) != included.contains(&bond.atoms[1]))
                        .then_some(i as u32)
                })
                .collect();
            vec![RepeatUnit {
                atoms,
                crossing_bonds,
                attachment_points: fragment
                    .polymer_attachment_points
                    .as_ref()
                    .unwrap()
                    .iter()
                    .map(|&atom| index(atom))
                    .collect::<Result<_, _>>()?,
                connectivity: Connectivity::HeadToTail,
                label: None,
            }]
        } else {
            Vec::new()
        };
        let parent_set: HashSet<_> = self.principal_parent_atoms.iter().copied().collect();
        let parent = order
            .iter()
            .enumerate()
            .filter_map(|(i, atom)| parent_set.contains(atom).then_some(i as u32))
            .collect();
        let molecule = Molecule {
            atoms,
            bonds,
            stereo,
            stereo_groups,
            position_variation,
            polymers,
            parent,
        };
        molecule.validate().map_err(IrError::Validation)?;
        Ok((molecule, order))
    }
}

/// Assign only radicals on the already finalized, localized graph. Every H count
/// is fixed, including counts originating from unbracketed resource atoms. This
/// avoids SMILES parsing, implicit-H inference, aromaticity, or stereo perception.
fn assign_radical_electrons(atoms: &mut [Atom], bonds: &[Bond]) -> Result<(), IrError> {
    use nazuchem_kernel::graph::{Atom as KernelAtom, Bond as KernelBond, Graph as KernelGraph};

    let graph = KernelGraph {
        atoms: atoms
            .iter()
            .map(|atom| KernelAtom {
                atomic_number: atom.element,
                isotope: atom.isotope.unwrap_or(0),
                charge: atom.charge,
                explicit_hydrogens: atom.hydrogens,
                no_implicit: true,
                ..KernelAtom::default()
            })
            .collect(),
        bonds: bonds
            .iter()
            .map(|bond| {
                Ok(KernelBond {
                    a: bond.atoms[0] as usize,
                    b: bond.atoms[1] as usize,
                    order: match bond.order {
                        BondOrder::Single => 1,
                        BondOrder::Double => 2,
                        BondOrder::Triple => 3,
                        BondOrder::Quadruple => 6,
                        _ => return Err(invalid("Unsupported IR bond order")),
                    },
                    aromatic: false,
                })
            })
            .collect::<Result<_, IrError>>()?,
    };
    let radicals = graph.assign_radicals().map_err(|reason| {
        invalid(format!(
            "Cannot assign finalized radical electrons: {reason}"
        ))
    })?;
    for (atom, radical_electrons) in atoms.iter_mut().zip(radicals) {
        atom.radical_electrons = radical_electrons;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Graph;

    #[test]
    fn isotope_vertices_all_locants_and_parent_projection_survive() {
        let mut graph = Graph::default();
        let fragment = graph.add_fragment("test");
        let carbon = graph.add_atom(fragment, Element::C);
        let deuterium = graph.add_atom(fragment, Element::H);
        graph.atom_mut(deuterium).isotope = Some(2);
        graph.atom_mut(carbon).locants = vec!["1".into(), "alpha".into()];
        graph.atom_mut(deuterium).locants = vec!["D".into()];
        graph.add_bond(carbon, deuterium, 1).unwrap();
        let structure = Structure {
            graph,
            fragment,
            principal_parent_atoms: vec![carbon],
            unknown_stereo_atoms: Vec::new(),
            unknown_stereo_bonds: Vec::new(),
        };
        let molecule = structure.to_ir().unwrap();
        assert_eq!(molecule.atoms.len(), 2);
        assert_eq!(molecule.atoms[0].locants, ["1", "alpha"]);
        assert_eq!(molecule.atoms[1].element, 1);
        assert_eq!(molecule.atoms[1].isotope, Some(2));
        assert_eq!(molecule.atoms[1].locants, ["D"]);
        assert_eq!(molecule.parent, [0]);
    }

    #[test]
    fn multi_attachment_and_wildcard_roles_are_structural_data() {
        let mut graph = Graph::default();
        let fragment = graph.add_fragment("test");
        let attachment = graph.add_atom(fragment, Element::R);
        let first = graph.add_atom(fragment, Element::C);
        let second = graph.add_atom(fragment, Element::C);
        graph.add_bond(attachment, first, 1).unwrap();
        graph.add_bond(first, second, 1).unwrap();
        graph
            .atom_mut(attachment)
            .properties
            .position_variation_bond = Some(vec![first, second]);
        let structure = Structure {
            graph,
            fragment,
            principal_parent_atoms: Vec::new(),
            unknown_stereo_atoms: Vec::new(),
            unknown_stereo_bonds: Vec::new(),
        };
        let molecule = structure.to_ir().unwrap();
        assert_eq!(molecule.atoms[0].label, None);
        assert_eq!(molecule.atoms[0].attachment_point, Some(1));
        assert_eq!(molecule.atoms[0].radical_electrons, 0);
        assert_eq!(
            molecule.position_variation,
            [PositionVariation {
                bond: 0,
                endpoints: vec![1, 2]
            }]
        );
    }

    #[test]
    fn folding_preserves_every_source_locant_parent_and_hydrogen_role() {
        for role in [
            "ordinary", "locants", "parent", "isotope", "charge", "label", "class", "stated",
        ] {
            let mut graph = Graph::default();
            let fragment = graph.add_fragment("test");
            let carbon = graph.add_atom(fragment, Element::C);
            let hydrogen = graph.add_atom(fragment, Element::H);
            graph.add_bond(carbon, hydrogen, 1).unwrap();
            graph.atom_mut(carbon).locants = vec!["1".into(), "alpha".into()];
            let mut principal_parent_atoms = vec![carbon];
            match role {
                "locants" => graph.atom_mut(hydrogen).locants = vec!["H".into(), "H1".into()],
                "parent" => principal_parent_atoms.push(hydrogen),
                "isotope" => graph.atom_mut(hydrogen).isotope = Some(2),
                "charge" => graph.atom_mut(hydrogen).charge = 1,
                "label" => {
                    graph.atom_mut(hydrogen).properties.homology_group = Some("hydrogen".into())
                }
                "class" => graph.atom_mut(hydrogen).properties.atom_class = Some(17),
                "stated" => graph.atom_mut(hydrogen).properties.smiles_hydrogen_count = Some(0),
                _ => {}
            }
            let structure = Structure {
                graph,
                fragment,
                principal_parent_atoms,
                unknown_stereo_atoms: Vec::new(),
                unknown_stereo_bonds: Vec::new(),
            };
            let (molecule, source_atoms) = structure.to_ir_with_atom_mapping().unwrap();
            assert_eq!(
                source_atoms.contains(&hydrogen),
                role != "ordinary",
                "{role}"
            );
            for (index, &source) in source_atoms.iter().enumerate() {
                let source = structure.graph.atom(source);
                assert_eq!(molecule.atoms[index].locants, source.locants, "{role}");
                assert_eq!(
                    molecule.parent.contains(&(index as u32)),
                    structure.principal_parent_atoms.contains(&source.id),
                    "{role}"
                );
                assert_eq!(
                    molecule.atoms[index].label,
                    source
                        .properties
                        .homology_group
                        .as_ref()
                        .map(|group| format!("{group}_p")),
                    "{role}"
                );
                assert_eq!(
                    molecule.atoms[index].atom_map, source.properties.atom_class,
                    "{role}"
                );
                assert_eq!(
                    molecule.atoms[index].hydrogens_stated,
                    source.properties.smiles_hydrogen_count.is_some(),
                    "{role}"
                );
            }
            assert_eq!(
                molecule.atoms[0].hydrogens,
                u8::from(role == "ordinary"),
                "{role}"
            );
        }
    }

    #[test]
    fn explicit_hydrido_locants_survive_before_the_legacy_writer_folds_them() {
        let parser = crate::Parser::new().unwrap();
        let result = parser.parse(
            "dihydrido(naphthalen-2-yl)rhenium",
            &crate::ParseOptions::strict(),
        );
        let structure = result.structure.unwrap();
        let hydrides: Vec<_> = structure
            .graph
            .fragment(structure.fragment)
            .atoms
            .iter()
            .copied()
            .filter(|&id| {
                structure.graph.atom(id).element == Element::H
                    && !structure.graph.atom(id).locants.is_empty()
            })
            .collect();
        assert_eq!(hydrides.len(), 2);
        let (molecule, source_atoms) = structure.to_ir_with_atom_mapping().unwrap();
        for hydrogen in hydrides {
            let index = source_atoms.iter().position(|&id| id == hydrogen).unwrap();
            assert_eq!(molecule.atoms[index].element, 1);
            assert_eq!(
                molecule.atoms[index].locants,
                structure.graph.atom(hydrogen).locants
            );
        }
        let parents: HashSet<_> = molecule
            .parent
            .iter()
            .map(|&index| source_atoms[index as usize])
            .collect();
        assert_eq!(
            parents,
            structure.principal_parent_atoms.iter().copied().collect()
        );
        for (index, &source) in source_atoms.iter().enumerate() {
            assert_eq!(
                molecule.atoms[index].locants,
                structure.graph.atom(source).locants
            );
        }
    }

    #[test]
    fn mapped_atoms_attachment_points_and_polymer_stars_have_typed_roles() {
        let mut graph = Graph::default();
        let fragment = graph.add_fragment("test");
        let first = graph.add_atom(fragment, Element::R);
        let carbon = graph.add_atom(fragment, Element::C);
        let last = graph.add_atom(fragment, Element::R);
        graph.add_bond(first, carbon, 1).unwrap();
        graph.add_bond(carbon, last, 1).unwrap();
        graph.atom_mut(first).properties.atom_class = Some(27);
        graph.atom_mut(carbon).properties.atom_class = Some(42);
        let mut structure = Structure {
            graph,
            fragment,
            principal_parent_atoms: vec![carbon],
            unknown_stereo_atoms: Vec::new(),
            unknown_stereo_bonds: Vec::new(),
        };
        let (molecule, sources) = structure.to_ir_with_atom_mapping().unwrap();
        let first_index = sources.iter().position(|&id| id == first).unwrap();
        let carbon_index = sources.iter().position(|&id| id == carbon).unwrap();
        assert_eq!(molecule.atoms[first_index].attachment_point, Some(27));
        assert_eq!(molecule.atoms[first_index].atom_map, None);
        assert_eq!(molecule.atoms[carbon_index].atom_map, Some(42));
        assert_eq!(molecule.atoms[carbon_index].attachment_point, None);
        assert!(molecule.atoms.iter().all(|atom| atom.label.is_none()));

        structure
            .graph
            .fragment_mut(fragment)
            .polymer_attachment_points = Some(vec![first, last]);
        let (molecule, sources) = structure.to_ir_with_atom_mapping().unwrap();
        assert!(
            molecule
                .atoms
                .iter()
                .all(|atom| atom.label.is_none() && atom.attachment_point.is_none())
        );
        assert_eq!(molecule.polymers.len(), 1);
        assert_eq!(
            molecule.polymers[0]
                .attachment_points
                .iter()
                .map(|&index| sources[index as usize])
                .collect::<Vec<_>>(),
            [first, last]
        );
        assert_eq!(molecule.polymers[0].atoms.len(), 1);
        assert_eq!(molecule.polymers[0].crossing_bonds.len(), 2);

        // A different wildcard in the same fragment remains inside the repeat
        // unit and keeps its ordinary attachment role.
        let pendant = structure.graph.add_atom(fragment, Element::R);
        structure.graph.add_bond(carbon, pendant, 1).unwrap();
        structure.graph.atom_mut(pendant).properties.atom_class = Some(91);
        let (molecule, sources) = structure.to_ir_with_atom_mapping().unwrap();
        let pendant = sources.iter().position(|&id| id == pendant).unwrap() as u32;
        assert_eq!(molecule.atoms[pendant as usize].attachment_point, Some(91));
        assert_eq!(molecule.atoms[pendant as usize].label, None);
        assert!(molecule.polymers[0].atoms.contains(&pendant));
        assert!(!molecule.polymers[0].attachment_points.contains(&pendant));
        assert_eq!(molecule.polymers[0].crossing_bonds.len(), 2);
    }

    #[test]
    fn hydrogen_references_of_double_bonds_and_axes_remain_vertices() {
        use crate::graph::{AxialStereo, BondStereo};
        for axial in [false, true] {
            let mut graph = Graph::default();
            let fragment = graph.add_fragment("test");
            let left = graph.add_atom(fragment, Element::C);
            let right = graph.add_atom(fragment, Element::C);
            let hydrogen = graph.add_atom(fragment, Element::H);
            let fluorine = graph.add_atom(fragment, Element::F);
            let chlorine = graph.add_atom(fragment, Element::Cl);
            let bromine = graph.add_atom(fragment, Element::Br);
            let bond = graph
                .add_bond(left, right, if axial { 1 } else { 2 })
                .unwrap();
            graph.add_bond(left, hydrogen, 1).unwrap();
            graph.add_bond(left, fluorine, 1).unwrap();
            graph.add_bond(right, chlorine, 1).unwrap();
            graph.add_bond(right, bromine, 1).unwrap();
            if axial {
                // H is the second ligand, rather than one of OPSIN's two
                // selected references. Its axis role still needs an identity.
                graph.bond_mut(bond).axial_stereo = Some(AxialStereo {
                    refs: [fluorine, chlorine],
                    clockwise: true,
                    descriptor: "P".into(),
                });
            } else {
                graph.bond_mut(bond).stereo = Some(BondStereo {
                    atom_refs: [hydrogen, left, right, chlorine],
                    value: BondStereoValue::Trans,
                });
            }
            let structure = Structure {
                graph,
                fragment,
                principal_parent_atoms: Vec::new(),
                unknown_stereo_atoms: Vec::new(),
                unknown_stereo_bonds: Vec::new(),
            };
            let (molecule, sources) = structure.to_ir_with_atom_mapping().unwrap();
            let hydrogen = sources.iter().position(|&id| id == hydrogen).unwrap() as u32;
            assert_eq!(molecule.atoms[hydrogen as usize].element, 1);
            assert!(molecule.stereo.iter().any(|stereo| match stereo {
                Stereo::Atropisomer { refs, .. } => refs.contains(&Ref::Atom(hydrogen)),
                Stereo::DoubleBond { refs, .. } => refs.contains(&hydrogen),
                _ => false,
            }));
        }
    }

    /// Read the IR's documented structural convention directly into the kernel
    /// metadata. No adapter helper, source graph parity, or SMILES perception is
    /// involved in this independent CIP check.
    fn ir_cip_labels(molecule: &Molecule) -> nazuchem_kernel::stereo::perception::State {
        use nazuchem_kernel::{
            electronic::Hybridization,
            graph::{Atom as KernelAtom, Bond as KernelBond, Graph as KernelGraph},
            kekulize::Direction,
            ranking::Metadata,
            stereo::{
                cip::{configuration::parity4, label},
                perception::{Properties, RingCache, State},
            },
        };
        let graph = KernelGraph {
            atoms: molecule
                .atoms
                .iter()
                .map(|atom| KernelAtom {
                    atomic_number: atom.element,
                    isotope: atom.isotope.unwrap_or_default(),
                    charge: atom.charge,
                    explicit_hydrogens: atom.hydrogens,
                    no_implicit: true,
                    ..KernelAtom::default()
                })
                .collect(),
            bonds: molecule
                .bonds
                .iter()
                .map(|bond| KernelBond {
                    a: bond.atoms[0] as usize,
                    b: bond.atoms[1] as usize,
                    order: match bond.order {
                        BondOrder::Single => 1,
                        BondOrder::Double => 2,
                        BondOrder::Triple => 3,
                        _ => panic!("test order"),
                    },
                    aromatic: false,
                })
                .collect(),
        };
        let mut metadata = Metadata::unspecified(&graph);
        for stereo in &molecule.stereo {
            match stereo {
                Stereo::Tetrahedral {
                    center,
                    neighbors,
                    winding,
                } => {
                    let center = *center as usize;
                    let mut carriers: Vec<_> = graph
                        .bonds
                        .iter()
                        .filter_map(|bond| {
                            if bond.a == center {
                                Some(bond.b)
                            } else if bond.b == center {
                                Some(bond.a)
                            } else {
                                None
                            }
                        })
                        .collect();
                    if carriers.len() == 3 {
                        carriers.push(center);
                    }
                    let refs = neighbors.map(|reference| match reference {
                        Ref::Atom(index) => index as usize,
                        Ref::ImplicitHydrogen | Ref::LonePair => center,
                        _ => panic!("test reference"),
                    });
                    let permutation = parity4(&refs, &carriers).unwrap();
                    assert_ne!(permutation, 0);
                    // Kernel 1 is clockwise viewed from carrier[0] toward the
                    // centre; an odd reference permutation reverses winding.
                    metadata.atoms[center].chiral_tag = match winding {
                        Winding::Clockwise => {
                            if permutation == 2 {
                                1
                            } else {
                                2
                            }
                        }
                        Winding::Anticlockwise => {
                            if permutation == 2 {
                                2
                            } else {
                                1
                            }
                        }
                        Winding::Unknown => 0,
                        _ => panic!("test winding"),
                    };
                }
                Stereo::Atropisomer {
                    bond,
                    refs,
                    winding,
                } => {
                    let bond_index = *bond as usize;
                    let bond = &graph.bonds[bond_index];
                    let mut clockwise = match winding {
                        AxialWinding::Clockwise => true,
                        AxialWinding::Anticlockwise => false,
                        _ => panic!("test axial winding"),
                    };
                    let mut references = Vec::new();
                    for (end, partner, reference) in
                        [(bond.a, bond.b, refs[0]), (bond.b, bond.a, refs[2])]
                    {
                        let Ref::Atom(reference) = reference else {
                            panic!("test explicit reference")
                        };
                        let carrier = graph
                            .bonds
                            .iter()
                            .filter_map(|bond| {
                                if bond.a == end {
                                    Some(bond.b)
                                } else if bond.b == end {
                                    Some(bond.a)
                                } else {
                                    None
                                }
                            })
                            .filter(|&other| other != partner)
                            .min()
                            .unwrap();
                        clockwise ^= carrier != reference as usize;
                        references.push(reference as usize);
                    }
                    metadata.bonds[bond_index].stereo = if clockwise { 6 } else { 7 };
                    metadata.bonds[bond_index].stereo_atoms = references;
                }
                _ => {}
            }
        }
        let state = State {
            valences: graph.valences().unwrap(),
            directions: vec![Direction::None; graph.bonds.len()],
            conjugated: vec![false; graph.bonds.len()],
            hybridizations: vec![Hybridization::Unspecified; graph.atoms.len()],
            properties: Properties::unspecified(&graph),
            rings: RingCache::default(),
            graph,
            metadata,
        };
        label::assign(&state, &Default::default()).unwrap().state
    }

    #[test]
    fn tetrahedral_first_reference_view_matches_named_r_and_s_cip() {
        let parser = crate::Parser::new().unwrap();
        for (name, expected) in [
            ("(2R)-butan-2-ol", "R"),
            ("(2S)-butan-2-ol", "S"),
            ("(R)-bromochlorofluoroiodomethane", "R"),
            ("(S)-bromochlorofluoroiodomethane", "S"),
        ] {
            let result = parser.parse(name, &crate::ParseOptions::strict());
            let structure = result
                .structure
                .unwrap_or_else(|| panic!("{name}: {}", result.message));
            let (molecule, source_atoms) = structure.to_ir_with_atom_mapping().unwrap();
            let labeled = ir_cip_labels(&molecule);
            let labels: Vec<_> = labeled
                .properties
                .atoms
                .iter()
                .filter_map(|atom| atom.cip_code.as_deref())
                .collect();
            assert_eq!(labels, [expected], "{name}");
            for &id in &source_atoms {
                if let Some(parity) = &structure.graph.atom(id).parity {
                    for reference in parity.atom_refs {
                        if let Some(StereoReference::Atom(reference)) = reference {
                            assert!(source_atoms.contains(&reference), "{name}: stereo H folded");
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn both_binaphthol_hands_and_ra_sa_keep_independent_axial_cip() {
        let parser = crate::Parser::new().unwrap();
        for (descriptor, expected) in [("M", "M"), ("P", "P"), ("Ra", "M"), ("Sa", "P")] {
            let name = format!("(1{descriptor})-[1,1'-binaphthalene]-2,2'-diol");
            let structure = parser
                .parse(&name, &crate::ParseOptions::strict())
                .structure
                .unwrap();
            let molecule = structure.to_ir().unwrap();
            let labeled = ir_cip_labels(&molecule);
            let codes: Vec<_> = labeled
                .properties
                .bond_codes
                .iter()
                .filter_map(|code| code.as_deref())
                .collect();
            assert_eq!(codes, [expected], "{name}");
            assert_eq!(
                Molecule::from_json(&molecule.to_json().unwrap()).unwrap(),
                molecule
            );
        }
    }

    #[test]
    fn radical_and_attachment_parse_options_remain_distinct() {
        let parser = crate::Parser::new().unwrap();
        for wildcard in [false, true] {
            let result = parser.parse(
                "methyl",
                &crate::ParseOptions {
                    allow_radicals: true,
                    output_radicals_as_wildcard_atoms: wildcard,
                    ..crate::ParseOptions::strict()
                },
            );
            let molecule = result.structure.unwrap().to_ir().unwrap();
            assert_eq!(
                molecule
                    .atoms
                    .iter()
                    .map(|atom| atom.radical_electrons)
                    .sum::<u8>(),
                u8::from(!wildcard)
            );
            assert_eq!(
                molecule
                    .atoms
                    .iter()
                    .filter(|atom| atom.element == 0)
                    .count(),
                usize::from(wildcard)
            );
        }
    }

    #[test]
    fn position_variation_keeps_a_reverse_source_with_any_endpoint_elements() {
        for elements in [
            [Element::R, Element::R],
            [Element::C, Element::C],
            [Element::R, Element::C],
            [Element::C, Element::R],
        ] {
            let mut graph = Graph::default();
            let fragment = graph.add_fragment("test");
            let first = graph.add_atom(fragment, elements[0]);
            let source = graph.add_atom(fragment, elements[1]);
            let alternative = graph.add_atom(fragment, Element::C);
            graph.add_bond(first, source, 1).unwrap();
            graph.atom_mut(source).properties.position_variation_bond =
                Some(vec![first, alternative]);
            let structure = Structure {
                graph,
                fragment,
                principal_parent_atoms: Vec::new(),
                unknown_stereo_atoms: Vec::new(),
                unknown_stereo_bonds: Vec::new(),
            };
            let molecule = structure.to_ir().unwrap();
            let source_index = crate::smiles::semantic_atom_order(&structure.graph, fragment)
                .unwrap()
                .iter()
                .position(|&atom| atom == source)
                .unwrap() as u32;
            let variation = &molecule.position_variation[0];
            assert_eq!(
                molecule.bonds[variation.bond as usize].atoms[0],
                source_index
            );
            assert_eq!(
                molecule.semantic_cxsmiles_for_tests().unwrap(),
                structure.semantic_cxsmiles().unwrap(),
                "{elements:?}"
            );
            let restored = Molecule::from_json(&molecule.to_json().unwrap()).unwrap();
            assert_eq!(restored, molecule);
            assert_eq!(
                restored.semantic_cxsmiles_for_tests().unwrap(),
                structure.semantic_cxsmiles().unwrap()
            );
        }
    }

    #[test]
    fn elementary_radicals_do_not_require_construction_attachment_sites() {
        let parser = crate::Parser::new().unwrap();
        for (name, electrons) in [
            ("oxygen", 2),
            ("nitrogen", 3),
            ("carbon", 4),
            ("sulfur", 2),
            ("sodium", 1),
            ("chloride", 0),
            ("iron", 0),
        ] {
            let result = parser.parse(name, &crate::ParseOptions::strict());
            let structure = result.structure.unwrap();
            assert!(
                structure
                    .graph
                    .fragment(structure.fragment)
                    .out_atoms
                    .is_empty(),
                "{name}"
            );
            let molecule = structure.to_ir().unwrap();
            assert_eq!(molecule.atoms[0].radical_electrons, electrons, "{name}");
            assert_eq!(
                molecule.semantic_cxsmiles_for_tests().unwrap(),
                structure.semantic_cxsmiles().unwrap(),
                "{name}"
            );
        }
    }
}
