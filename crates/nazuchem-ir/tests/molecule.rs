// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | crates/nazuchem-ir/tests/molecule.rs | Phase B step 8a

use nazuchem_ir::{
    Atom, Bond, BondOrder, Connectivity, JsonError, Molecule, PositionVariation, Ref, Relation,
    RepeatUnit, Stereo, StereoGroup, StereoGroupKind, ValidationError, Winding,
};
use serde_json::json;

fn carbon(hydrogens: u8) -> Atom {
    Atom {
        element: 6,
        hydrogens,
        ..Atom::default()
    }
}

fn bond(a: u32, b: u32, order: BondOrder) -> Bond {
    Bond {
        atoms: [a, b],
        order,
    }
}

fn tetrahedral() -> Molecule {
    let mut molecule = Molecule::new(
        vec![
            carbon(1),
            carbon(3),
            Atom {
                element: 9,
                ..Atom::default()
            },
            Atom {
                element: 17,
                ..Atom::default()
            },
        ],
        vec![
            bond(0, 1, BondOrder::Single),
            bond(0, 2, BondOrder::Single),
            bond(0, 3, BondOrder::Single),
        ],
    )
    .unwrap();
    molecule.stereo.push(Stereo::Tetrahedral {
        center: 0,
        neighbors: [
            Ref::Atom(1),
            Ref::Atom(2),
            Ref::Atom(3),
            Ref::ImplicitHydrogen,
        ],
        winding: Winding::Anticlockwise,
        known: false,
    });
    molecule
}

fn alkene() -> Molecule {
    let mut molecule = Molecule::new(
        vec![carbon(1), carbon(1), carbon(3), carbon(3)],
        vec![
            bond(0, 1, BondOrder::Double),
            bond(0, 2, BondOrder::Single),
            bond(1, 3, BondOrder::Single),
        ],
    )
    .unwrap();
    molecule.stereo.push(Stereo::DoubleBond {
        bond: 0,
        refs: [2, 3],
        relation: Relation::Trans,
        known: true,
    });
    molecule
}

#[test]
fn all_annotations_round_trip_without_reordering_or_folding_hydrogen() {
    let mut molecule = tetrahedral();
    molecule.atoms[0].locants = vec!["2".into(), "2′".into()];
    molecule.atoms[1].hydrogens_stated = true;
    molecule.atoms.push(Atom {
        element: 1,
        isotope: Some(2),
        ..Atom::default()
    });
    molecule.atoms.push(Atom {
        label: Some("R;λ".into()),
        ..Atom::default()
    });
    molecule.atoms.push(Atom {
        element: 8,
        charge: -1,
        radical_electrons: 1,
        ..Atom::default()
    });
    molecule.bonds.push(bond(1, 4, BondOrder::Single));
    molecule.bonds.push(bond(1, 5, BondOrder::Single));
    molecule.parent = vec![1, 0];
    molecule.stereo_groups.push(StereoGroup {
        kind: StereoGroupKind::And(7),
        stereo: vec![0],
    });
    molecule.position_variation.push(PositionVariation {
        bond: 4,
        endpoints: vec![1, 2, 3],
    });
    molecule.polymers.push(RepeatUnit {
        atoms: vec![1, 4],
        crossing_bonds: vec![0, 4],
        connectivity: Connectivity::HeadToHead,
        label: Some("n".into()),
    });

    let encoded = molecule.to_json().unwrap();
    assert_eq!(
        serde_json::from_str::<serde_json::Value>(&encoded).unwrap()["nazuchem_ir"],
        1
    );
    assert_eq!(Molecule::from_json(&encoded).unwrap(), molecule);
    assert_eq!(
        serde_json::from_str::<Molecule>(&encoded).unwrap(),
        molecule
    );
    assert_eq!(
        Molecule::from_json(&encoded).unwrap().atoms[4].isotope,
        Some(2)
    );
}

#[test]
fn lone_pairs_and_double_bond_reference_orientation_are_validated() {
    let mut molecule = tetrahedral();
    let Stereo::Tetrahedral { neighbors, .. } = &mut molecule.stereo[0] else {
        unreachable!()
    };
    neighbors[3] = Ref::LonePair;
    molecule.atoms[0].hydrogens = 0;
    molecule.validate().unwrap();

    let mut molecule = alkene();
    molecule.validate().unwrap();
    let Stereo::DoubleBond { refs, .. } = &mut molecule.stereo[0] else {
        unreachable!()
    };
    refs.swap(0, 1);
    assert!(matches!(
        molecule.validate(),
        Err(ValidationError::InvalidStereo { .. })
    ));
}

#[test]
fn graph_rejects_invalid_elements_isotopes_and_bond_identities() {
    assert!(matches!(
        Molecule::new(
            vec![Atom {
                element: 119,
                ..Atom::default()
            }],
            vec![]
        ),
        Err(ValidationError::InvalidElement { .. })
    ));
    assert!(matches!(
        Molecule::new(
            vec![Atom {
                isotope: Some(0),
                ..carbon(0)
            }],
            vec![]
        ),
        Err(ValidationError::InvalidIsotope { .. })
    ));
    assert!(matches!(
        Molecule::new(vec![carbon(0)], vec![bond(0, 1, BondOrder::Single)]),
        Err(ValidationError::IndexOutOfBounds { .. })
    ));
    assert!(matches!(
        Molecule::new(vec![carbon(0)], vec![bond(0, 0, BondOrder::Single)]),
        Err(ValidationError::SelfBond { .. })
    ));
    assert!(matches!(
        Molecule::new(
            vec![carbon(0), carbon(0)],
            vec![bond(0, 1, BondOrder::Single), bond(1, 0, BondOrder::Double)]
        ),
        Err(ValidationError::DuplicateBond { .. })
    ));
    // Valence is the kernel's job: the IR preserves explicitly stated data.
    Molecule::new(vec![carbon(255)], vec![]).unwrap();
}

#[test]
fn stereo_rejects_duplicate_nonadjacent_and_unavailable_references() {
    let mut duplicate = tetrahedral();
    let Stereo::Tetrahedral { neighbors, .. } = &mut duplicate.stereo[0] else {
        unreachable!()
    };
    neighbors[1] = neighbors[0];
    assert!(matches!(
        duplicate.validate(),
        Err(ValidationError::InvalidStereo { .. })
    ));

    let mut missing_hydrogen = tetrahedral();
    missing_hydrogen.atoms[0].hydrogens = 0;
    assert!(matches!(
        missing_hydrogen.validate(),
        Err(ValidationError::InvalidStereo { .. })
    ));

    let mut omitted_neighbor = tetrahedral();
    omitted_neighbor.atoms.push(carbon(0));
    omitted_neighbor.bonds.push(bond(0, 4, BondOrder::Single));
    assert!(matches!(
        omitted_neighbor.validate(),
        Err(ValidationError::InvalidStereo { .. })
    ));

    let mut duplicate_entry = alkene();
    duplicate_entry
        .stereo
        .push(duplicate_entry.stereo[0].clone());
    assert!(matches!(
        duplicate_entry.validate(),
        Err(ValidationError::InvalidStereo { .. })
    ));

    let mut wrong_order = alkene();
    wrong_order.bonds[0].order = BondOrder::Single;
    assert!(matches!(
        wrong_order.validate(),
        Err(ValidationError::InvalidStereo { .. })
    ));
}

#[test]
fn group_parent_position_variation_and_polymer_references_are_checked() {
    let mut invalid_group = tetrahedral();
    invalid_group.stereo_groups.push(StereoGroup {
        kind: StereoGroupKind::Abs,
        stereo: vec![1],
    });
    assert!(matches!(
        invalid_group.validate(),
        Err(ValidationError::IndexOutOfBounds { .. })
    ));
    invalid_group.stereo_groups[0].stereo = vec![0];
    invalid_group.stereo_groups.push(StereoGroup {
        kind: StereoGroupKind::Or(1),
        stereo: vec![0],
    });
    assert!(matches!(
        invalid_group.validate(),
        Err(ValidationError::InvalidStereoGroup { .. })
    ));

    let mut invalid_parent = Molecule::new(vec![carbon(0)], vec![]).unwrap();
    invalid_parent.parent = vec![0, 0];
    assert!(matches!(
        invalid_parent.validate(),
        Err(ValidationError::DuplicateIndex { .. })
    ));
    invalid_parent.parent = vec![1];
    assert!(matches!(
        invalid_parent.validate(),
        Err(ValidationError::IndexOutOfBounds { .. })
    ));

    let mut invalid_variation = tetrahedral();
    invalid_variation
        .position_variation
        .push(PositionVariation {
            bond: 0,
            endpoints: vec![1, 1],
        });
    assert!(matches!(
        invalid_variation.validate(),
        Err(ValidationError::DuplicateIndex { .. })
    ));
    invalid_variation.position_variation[0].endpoints.clear();
    assert!(matches!(
        invalid_variation.validate(),
        Err(ValidationError::InvalidPositionVariation { .. })
    ));

    let mut invalid_polymer = tetrahedral();
    invalid_polymer.polymers.push(RepeatUnit {
        atoms: vec![0, 1],
        crossing_bonds: vec![0],
        connectivity: Connectivity::Unknown,
        label: None,
    });
    assert!(matches!(
        invalid_polymer.validate(),
        Err(ValidationError::InvalidRepeatUnit { .. })
    ));
    invalid_polymer.polymers[0].crossing_bonds = vec![1, 2];
    invalid_polymer.validate().unwrap();
}

#[test]
fn json_rejects_version_unknown_fields_variants_and_invalid_graphs() {
    let molecule = tetrahedral();
    let original = serde_json::to_value(&molecule).unwrap();
    let mut version = original.clone();
    version["nazuchem_ir"] = json!(2);
    assert!(matches!(
        Molecule::from_json(&version.to_string()),
        Err(JsonError::UnsupportedVersion { version: 2 })
    ));
    assert!(serde_json::from_value::<Molecule>(version).is_err());

    let mut unknown_stereo_field = original.clone();
    unknown_stereo_field["stereo"][0]["future_field"] = json!(true);
    assert!(matches!(
        Molecule::from_json(&unknown_stereo_field.to_string()),
        Err(JsonError::Json(_))
    ));
    let mut unknown_field = original.clone();
    unknown_field["future_field"] = json!(true);
    assert!(Molecule::from_json(&unknown_field.to_string()).is_err());
    let mut unknown_atom_field = original.clone();
    unknown_atom_field["atoms"][0]["aromatic"] = json!(true);
    assert!(Molecule::from_json(&unknown_atom_field.to_string()).is_err());
    let mut unknown_variant = original.clone();
    unknown_variant["stereo"][0]["kind"] = json!("atropisomer");
    assert!(matches!(
        Molecule::from_json(&unknown_variant.to_string()),
        Err(JsonError::Json(_))
    ));
    let mut invalid_graph = original;
    invalid_graph["bonds"][0]["atoms"][1] = json!(999);
    assert!(matches!(
        Molecule::from_json(&invalid_graph.to_string()),
        Err(JsonError::Validation(
            ValidationError::IndexOutOfBounds { .. }
        ))
    ));
    assert!(serde_json::from_value::<Molecule>(invalid_graph).is_err());
}

#[test]
fn direct_serde_serialization_also_rejects_mutated_invalid_graphs() {
    let mut molecule = tetrahedral();
    molecule.parent.push(999);
    assert!(matches!(molecule.to_json(), Err(JsonError::Validation(_))));
    assert!(serde_json::to_value(molecule).is_err());
}

#[test]
fn empty_graph_is_valid_but_version_and_collections_are_required_on_input() {
    let molecule = Molecule::default();
    assert_eq!(
        Molecule::from_json(&molecule.to_json().unwrap()).unwrap(),
        molecule
    );
    assert!(Molecule::from_json("{}").is_err());
    assert!(Molecule::from_json("{\"nazuchem_ir\":1}").is_err());
    let mut missing_collection = serde_json::to_value(molecule).unwrap();
    missing_collection.as_object_mut().unwrap().remove("parent");
    assert!(Molecule::from_json(&missing_collection.to_string()).is_err());
}

#[test]
fn all_enum_types_reject_unknown_wire_variants_and_tag_fields() {
    assert!(serde_json::from_value::<BondOrder>(json!("aromatic")).is_err());
    assert!(serde_json::from_value::<Winding>(json!("unknown")).is_err());
    assert!(serde_json::from_value::<Relation>(json!("unspecified")).is_err());
    assert!(serde_json::from_value::<Connectivity>(json!("tail_to_tail")).is_err());
    assert!(serde_json::from_value::<Ref>(json!({"kind":"future"})).is_err());
    assert!(serde_json::from_value::<StereoGroupKind>(json!({"kind":"future"})).is_err());
    assert!(serde_json::from_value::<Ref>(json!({"kind":"atom","atom":0,"future":true})).is_err());
    assert!(serde_json::from_value::<Ref>(json!({"kind":"lone_pair","future":true})).is_err());
    assert!(
        serde_json::from_value::<StereoGroupKind>(json!({"kind":"and","id":1,"future":true}))
            .is_err()
    );
    assert!(
        serde_json::from_value::<StereoGroupKind>(json!({"kind":"abs","future":true})).is_err()
    );
}
