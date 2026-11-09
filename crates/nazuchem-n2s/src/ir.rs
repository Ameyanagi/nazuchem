// SPDX-License-Identifier: MIT
// Origin: original
// Source: Nazuchem | crates/nazuchem-n2s/src/ir.rs | Phase B step 8a

//! Projection of the finalized OPSIN graph into the versioned data contract.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::fmt;

use nazuchem_ir::{
    Atom, Bond, BondOrder, Connectivity, Molecule, PositionVariation, Ref, Relation, RepeatUnit,
    Stereo, StereoGroup, StereoGroupKind, Winding,
};

use crate::graph::{AtomId, BondStereoValue, Element, StereoGroupType, StereoReference};
use crate::{SerializationError, Structure};

/// A graph cannot be projected without discarding information.
#[non_exhaustive]
#[derive(Debug)]
pub enum IrError {
    /// The existing semantic exporter rejected the construction graph.
    Serialization(SerializationError),
    /// IR v1 reserves axial configurations for a later stereo variant.
    UnsupportedAxialStereo,
    /// Unfinished construction state or a value outside the IR's numeric range.
    InvalidGraph(String),
    /// The projected graph failed the IR's structural validation.
    Validation(nazuchem_ir::ValidationError),
}

impl fmt::Display for IrError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Serialization(error) => error.fmt(f),
            Self::UnsupportedAxialStereo => {
                f.write_str("IR v1 does not represent axial stereochemistry")
            }
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
    /// Atom identity follows the semantic CXSMILES writer; ordinary foldable H
    /// vertices become counts, while isotopic H and stereo-required H remain
    /// vertices. All locants on written atoms retain their original order.
    /// Axial configurations return [`IrError::UnsupportedAxialStereo`].
    pub fn to_ir(&self) -> Result<Molecule, IrError> {
        let graph = &self.graph;
        let fragment = graph
            .fragments
            .get(self.fragment.0)
            .filter(|fragment| fragment.active)
            .ok_or_else(|| invalid("Missing or inactive structure fragment"))?;
        if fragment
            .bonds
            .iter()
            .any(|id| graph.bond(*id).axial_stereo.is_some())
        {
            return Err(IrError::UnsupportedAxialStereo);
        }
        let order = crate::smiles::semantic_atom_order(graph, self.fragment)
            .map_err(IrError::Serialization)?;
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
        let polymer = fragment
            .polymer_attachment_points
            .as_ref()
            .is_some_and(|points| !points.is_empty());
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
            let label = if let Some(group) = &atom.properties.homology_group {
                // The legacy writer escapes every non-alphanumeric UTF-16
                // unit before adding this homology suffix (including '_').
                Some(format!("{group}_p"))
            } else if atom.element == Element::R {
                Some(if polymer {
                    "star_e".into()
                } else {
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
                    format!("_AP{class}")
                })
            } else {
                if atom.properties.atom_class.is_some() {
                    return Err(invalid(
                        "IR v1 has no atom-map field for a non-wildcard atom",
                    ));
                }
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
                label,
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
                1 => Winding::Clockwise,
                -1 => Winding::Anticlockwise,
                _ => return Err(invalid("Tetrahedral winding must be +1 or -1")),
            };
            let stereo_index = stereo.len() as u32;
            stereo.push(Stereo::Tetrahedral {
                center: index(id)?,
                neighbors,
                winding,
                known: parity.stereo_group.kind != StereoGroupType::Unknown,
            });
            if parity.stereo_group.kind != StereoGroupType::Unknown {
                let kind = match parity.stereo_group.kind {
                    StereoGroupType::Absolute => StereoGroupKind::Abs,
                    StereoGroupType::Racemic => StereoGroupKind::And(parity.stereo_group.number),
                    StereoGroupType::Relative => StereoGroupKind::Or(parity.stereo_group.number),
                    StereoGroupType::Unknown => unreachable!(),
                };
                groups.entry(kind).or_default().push(stereo_index);
            }
        }
        for &id in &fragment.bonds {
            let bond = graph.bond(id);
            let Some(configuration) = &bond.stereo else {
                continue;
            };
            let refs = configuration.atom_refs;
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
                known: true,
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
                        if let Stereo::DoubleBond { bond, refs, .. } = entry
                            && *bond == bond_index
                        {
                            refs.swap(0, 1);
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
                .filter_map(|(i, &id)| (graph.atom(id).element != Element::R).then_some(i as u32))
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
        Ok(molecule)
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
    fn multi_attachment_and_wildcard_labels_are_structural_data() {
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
        };
        let molecule = structure.to_ir().unwrap();
        assert_eq!(molecule.atoms[0].label.as_deref(), Some("_AP1"));
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
                molecule.semantic_cxsmiles().unwrap(),
                structure.semantic_cxsmiles().unwrap(),
                "{elements:?}"
            );
            let restored = Molecule::from_json(&molecule.to_json().unwrap()).unwrap();
            assert_eq!(restored, molecule);
            assert_eq!(
                restored.semantic_cxsmiles().unwrap(),
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
                molecule.semantic_cxsmiles().unwrap(),
                structure.semantic_cxsmiles().unwrap(),
                "{name}"
            );
        }
    }
}
