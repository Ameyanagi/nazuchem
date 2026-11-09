// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | crates/nazuchem-ir/tests/molecule.rs | Phase B step 8a

use nazuchem_ir::{
    Atom, AxialWinding, Bond, BondOrder, Connectivity, JsonError, Molecule, PositionVariation, Ref,
    Relation, RepeatUnit, Stereo, StereoGroup, StereoGroupKind, ValidationError, Winding,
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
        winding: Winding::Unknown,
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
    });
    molecule
}

#[test]
fn all_annotations_round_trip_without_reordering_or_folding_hydrogen() {
    let mut molecule = tetrahedral();
    let Stereo::Tetrahedral { winding, .. } = &mut molecule.stereo[0] else {
        unreachable!()
    };
    *winding = Winding::Clockwise;
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
        attachment_points: vec![5],
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
    let Stereo::Tetrahedral { winding, .. } = &mut invalid_group.stereo[0] else {
        unreachable!()
    };
    *winding = Winding::Clockwise;
    invalid_group.stereo_groups.push(StereoGroup {
        kind: StereoGroupKind::And(1),
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
        attachment_points: vec![],
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
        Err(JsonError::UnknownField { .. })
    ));
    let mut unknown_field = original.clone();
    unknown_field["future_field"] = json!(true);
    assert!(Molecule::from_json(&unknown_field.to_string()).is_err());
    let mut unknown_atom_field = original.clone();
    unknown_atom_field["atoms"][0]["aromatic"] = json!(true);
    assert!(Molecule::from_json(&unknown_atom_field.to_string()).is_err());
    let mut unknown_variant = original.clone();
    unknown_variant["stereo"][0]["kind"] = json!("extended_tetrahedral");
    assert!(matches!(
        Molecule::from_json(&unknown_variant.to_string()),
        Err(JsonError::UnknownVariant { .. })
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
fn optional_annotations_default_and_empty_fields_are_omitted() {
    let molecule = Molecule::default();
    assert_eq!(
        Molecule::from_json(&molecule.to_json().unwrap()).unwrap(),
        molecule
    );
    assert!(Molecule::from_json("{}").is_err());
    assert!(Molecule::from_json("{\"nazuchem_ir\":1}").is_err());
    let minimal = json!({"nazuchem_ir":1,"atoms":[],"bonds":[]});
    assert_eq!(serde_json::to_value(&molecule).unwrap(), minimal);
    assert_eq!(Molecule::from_json(&minimal.to_string()).unwrap(), molecule);
    assert_eq!(
        serde_json::from_value::<Molecule>(minimal).unwrap(),
        molecule
    );
    let atom = serde_json::to_value(carbon(4)).unwrap();
    for field in [
        "isotope",
        "charge",
        "hydrogens_stated",
        "radical_electrons",
        "label",
        "locants",
        "atom_map",
        "attachment_point",
    ] {
        assert!(atom.get(field).is_none(), "{field}");
    }
    assert_eq!(serde_json::from_value::<Atom>(atom).unwrap(), carbon(4));
}

#[test]
fn atom_defaults_are_optional_on_read_and_omitted_on_write() {
    let minimal = json!({"element":6,"hydrogens":4});
    assert_eq!(serde_json::to_value(carbon(4)).unwrap(), minimal);
    assert_eq!(
        serde_json::from_value::<Atom>(minimal.clone()).unwrap(),
        carbon(4)
    );
    let document = json!({"nazuchem_ir":1,"atoms":[minimal],"bonds":[]});
    let molecule = Molecule::new(vec![carbon(4)], vec![]).unwrap();
    assert_eq!(
        Molecule::from_json(&document.to_string()).unwrap(),
        molecule
    );
    assert_eq!(
        serde_json::from_value::<Molecule>(document).unwrap(),
        molecule
    );

    let mut stated = carbon(4);
    stated.charge = -1;
    stated.radical_electrons = 2;
    stated.hydrogens_stated = true;
    let value = serde_json::to_value(&stated).unwrap();
    assert_eq!(
        value,
        json!({"element":6,"hydrogens":4,"charge":-1,"radical_electrons":2,"hydrogens_stated":true})
    );
    assert_eq!(serde_json::from_value::<Atom>(value).unwrap(), stated);

    let explicit_defaults = json!({"element":6,"hydrogens":4,"charge":0,"radical_electrons":0,"hydrogens_stated":false});
    let read: Atom = serde_json::from_value(explicit_defaults).unwrap();
    assert_eq!(
        serde_json::to_value(read).unwrap(),
        json!({"element":6,"hydrogens":4})
    );
    for missing_required in [json!({"element":6}), json!({"hydrogens":4})] {
        assert!(serde_json::from_value::<Atom>(missing_required).is_err());
    }
}

#[test]
fn all_enum_types_reject_unknown_wire_variants_and_tag_fields() {
    assert!(serde_json::from_value::<BondOrder>(json!("aromatic")).is_err());
    assert!(serde_json::from_value::<BondOrder>(json!({"single":null})).is_err());
    assert!(serde_json::from_value::<Winding>(json!({"clockwise":null})).is_err());
    assert!(serde_json::from_value::<Relation>(json!({"cis":null})).is_err());
    assert!(serde_json::from_value::<AxialWinding>(json!({"clockwise":null})).is_err());
    assert!(serde_json::from_value::<AxialWinding>(json!("m")).is_err());
    assert!(serde_json::from_value::<AxialWinding>(json!("p")).is_err());
    assert!(serde_json::from_value::<Connectivity>(json!({"unknown":null})).is_err());
    assert!(serde_json::from_value::<Winding>(json!("future")).is_err());
    assert!(serde_json::from_value::<Relation>(json!("unspecified")).is_err());
    assert!(serde_json::from_value::<Connectivity>(json!("tail_to_tail")).is_err());
    assert!(serde_json::from_value::<Ref>(json!({"kind":"future"})).is_err());
    assert!(serde_json::from_value::<Ref>(json!({"lone_pair":null})).is_err());
    assert!(serde_json::from_value::<Ref>(json!({"implicit_hydrogen":null})).is_err());
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

fn axis() -> Molecule {
    let mut molecule = Molecule::new(
        vec![carbon(0); 6],
        vec![
            bond(0, 1, BondOrder::Single),
            bond(0, 2, BondOrder::Single),
            bond(0, 3, BondOrder::Single),
            bond(1, 4, BondOrder::Single),
            bond(1, 5, BondOrder::Single),
        ],
    )
    .unwrap();
    molecule.stereo.push(Stereo::Atropisomer {
        bond: 0,
        refs: [Ref::Atom(2), Ref::Atom(3), Ref::Atom(4), Ref::Atom(5)],
        winding: AxialWinding::Clockwise,
    });
    molecule
}

#[test]
fn untagged_refs_flattened_groups_and_typed_atom_annotations_have_exact_json() {
    let mut molecule = tetrahedral();
    let Stereo::Tetrahedral { winding, .. } = &mut molecule.stereo[0] else {
        unreachable!()
    };
    *winding = Winding::Clockwise;
    molecule.atoms[0].atom_map = Some(12);
    molecule.atoms.push(Atom {
        attachment_point: Some(2),
        ..Atom::default()
    });
    molecule.stereo_groups.push(StereoGroup {
        kind: StereoGroupKind::And(7),
        stereo: vec![0],
    });
    let value = serde_json::to_value(&molecule).unwrap();
    assert_eq!(
        value["stereo"][0]["neighbors"],
        json!([1, 2, 3, "implicit_hydrogen"])
    );
    assert_eq!(value["stereo"][0]["winding"], "clockwise");
    assert!(value["stereo"][0].get("known").is_none());
    assert_eq!(
        value["stereo_groups"],
        json!([{"kind":"and","id":7,"stereo":[0]}])
    );
    assert_eq!(value["atoms"][0]["atom_map"], 12);
    assert_eq!(value["atoms"][4]["attachment_point"], 2);
    assert_eq!(Molecule::from_json(&value.to_string()).unwrap(), molecule);
    assert_eq!(
        serde_json::from_value::<Ref>(json!("lone_pair")).unwrap(),
        Ref::LonePair
    );
    assert_eq!(
        serde_json::to_value(Ref::Atom(u32::MAX)).unwrap(),
        json!(u32::MAX)
    );
}

#[test]
fn axial_unknown_and_known_stereo_validate_and_round_trip() {
    for winding in [
        AxialWinding::Anticlockwise,
        AxialWinding::Clockwise,
        AxialWinding::Unknown,
    ] {
        let mut molecule = axis();
        let Stereo::Atropisomer { winding: value, .. } = &mut molecule.stereo[0] else {
            unreachable!()
        };
        *value = winding;
        assert_eq!(
            Molecule::from_json(&molecule.to_json().unwrap()).unwrap(),
            molecule
        );
    }
    let mut wrong_end = axis();
    let Stereo::Atropisomer { refs, .. } = &mut wrong_end.stereo[0] else {
        unreachable!()
    };
    refs.swap(1, 2);
    assert!(matches!(
        wrong_end.validate(),
        Err(ValidationError::InvalidStereo { .. })
    ));
    let mut duplicate = axis();
    duplicate.stereo.push(duplicate.stereo[0].clone());
    assert!(matches!(
        duplicate.validate(),
        Err(ValidationError::InvalidStereo { .. })
    ));
    let mut wrong_order = axis();
    wrong_order.bonds[0].order = BondOrder::Double;
    assert!(matches!(
        wrong_order.validate(),
        Err(ValidationError::InvalidStereo { .. })
    ));

    let mut virtual_refs = axis();
    virtual_refs.bonds.remove(4);
    virtual_refs.bonds.remove(2);
    virtual_refs.atoms[0].hydrogens = 1;
    virtual_refs.atoms[1].hydrogens = 1;
    let Stereo::Atropisomer { refs, .. } = &mut virtual_refs.stereo[0] else {
        unreachable!()
    };
    *refs = [
        Ref::Atom(2),
        Ref::ImplicitHydrogen,
        Ref::Atom(4),
        Ref::ImplicitHydrogen,
    ];
    virtual_refs.validate().unwrap();
    virtual_refs.atoms[1].hydrogens = 0;
    assert!(matches!(
        virtual_refs.validate(),
        Err(ValidationError::InvalidStereo { .. })
    ));

    let mut alkene = alkene();
    let Stereo::DoubleBond { relation, .. } = &mut alkene.stereo[0] else {
        unreachable!()
    };
    *relation = Relation::Unknown;
    assert_eq!(
        Molecule::from_json(&alkene.to_json().unwrap()).unwrap(),
        alkene
    );
}

#[test]
fn groups_require_known_configuration_for_every_stereo_kind() {
    for mut molecule in [tetrahedral(), alkene(), axis()] {
        match &mut molecule.stereo[0] {
            Stereo::Tetrahedral { winding, .. } => *winding = Winding::Clockwise,
            Stereo::DoubleBond { relation, .. } => *relation = Relation::Cis,
            Stereo::Atropisomer { winding, .. } => *winding = AxialWinding::Clockwise,
            _ => unreachable!(),
        }
        for kind in [StereoGroupKind::And(1), StereoGroupKind::Or(1)] {
            let mut grouped = molecule.clone();
            grouped.stereo_groups.push(StereoGroup {
                kind,
                stereo: vec![0],
            });
            grouped.validate().unwrap();
            let mut document = serde_json::to_value(&grouped).unwrap();
            let field = match &mut grouped.stereo[0] {
                Stereo::Tetrahedral { winding, .. } => {
                    *winding = Winding::Unknown;
                    "winding"
                }
                Stereo::DoubleBond { relation, .. } => {
                    *relation = Relation::Unknown;
                    "relation"
                }
                Stereo::Atropisomer { winding, .. } => {
                    *winding = AxialWinding::Unknown;
                    "winding"
                }
                _ => unreachable!(),
            };
            assert!(matches!(
                grouped.validate(),
                Err(ValidationError::InvalidStereoGroup { group: 0, .. })
            ));
            assert!(matches!(
                grouped.to_json(),
                Err(JsonError::Validation(ValidationError::InvalidStereoGroup {
                    group: 0,
                    ..
                }))
            ));
            assert!(serde_json::to_value(&grouped).is_err());
            document["stereo"][0][field] = json!("unknown");
            assert!(matches!(
                Molecule::from_json(&document.to_string()),
                Err(JsonError::Validation(ValidationError::InvalidStereoGroup {
                    group: 0,
                    ..
                }))
            ));
            assert!(serde_json::from_value::<Molecule>(document).is_err());
            grouped.stereo_groups.clear();
            grouped.validate().unwrap();
        }
    }
}

#[test]
fn repeat_units_require_the_complete_boundary_and_external_wildcard_end_roles() {
    let mut molecule = Molecule::new(
        vec![Atom::default(), carbon(2), carbon(2), Atom::default()],
        vec![
            bond(0, 1, BondOrder::Single),
            bond(1, 2, BondOrder::Single),
            bond(2, 3, BondOrder::Single),
        ],
    )
    .unwrap();
    molecule.polymers.push(RepeatUnit {
        atoms: vec![1, 2],
        crossing_bonds: vec![0],
        attachment_points: vec![0, 3],
        connectivity: Connectivity::HeadToTail,
        label: None,
    });
    assert!(matches!(
        molecule.validate(),
        Err(ValidationError::InvalidRepeatUnit { .. })
    ));
    molecule.polymers[0].crossing_bonds.push(2);
    molecule.validate().unwrap();
    molecule.polymers[0].attachment_points = vec![1];
    assert!(matches!(
        molecule.validate(),
        Err(ValidationError::InvalidRepeatUnit { .. })
    ));
    molecule.polymers[0].attachment_points = vec![0, 0];
    assert!(matches!(
        molecule.validate(),
        Err(ValidationError::DuplicateIndex { .. })
    ));
}

#[test]
fn unsupported_version_is_checked_before_new_fields_and_shapes() {
    for value in [
        json!({"nazuchem_ir":2,"future":true}),
        json!({"nazuchem_ir":2,"atoms":[["new format"]],"bonds":true}),
        json!({"nazuchem_ir":2,"stereo":[{"kind":"future"}]}),
        json!({"nazuchem_ir":0}),
    ] {
        let version = value["nazuchem_ir"].as_u64().unwrap();
        assert!(matches!(Molecule::from_json(&value.to_string()),
            Err(JsonError::UnsupportedVersion { version: actual }) if actual == version));
        let error = serde_json::from_value::<Molecule>(value)
            .unwrap_err()
            .to_string();
        assert!(error.contains("unsupported nazuchem IR version"), "{error}");
    }
}

#[test]
fn unknown_fields_and_variants_are_typed_and_per_variant_fields_are_exact() {
    let original = serde_json::to_value(tetrahedral()).unwrap();
    for (pointer, value, field) in [
        ("/stereo/0/known", json!(false), "known"),
        ("/stereo/0/relation", json!("cis"), "relation"),
        ("/atoms/0/aromatic", json!(true), "aromatic"),
    ] {
        let mut document = original.clone();
        let (object, _) = pointer.rsplit_once('/').unwrap();
        document.pointer_mut(object).unwrap()[field] = value;
        assert!(
            matches!(Molecule::from_json(&document.to_string()), Err(JsonError::UnknownField { field: actual, .. }) if actual == field)
        );
        assert!(serde_json::from_value::<Molecule>(document).is_err());
    }
    for (pointer, value) in [
        ("/stereo/0/winding", "future"),
        ("/stereo/0/neighbors/0", "virtual"),
        ("/bonds/0/order", "dative"),
        ("/stereo/0/kind", "future"),
    ] {
        let mut document = original.clone();
        *document.pointer_mut(pointer).unwrap() = json!(value);
        assert!(
            matches!(Molecule::from_json(&document.to_string()), Err(JsonError::UnknownVariant { variant, .. }) if variant == value)
        );
        assert!(serde_json::from_value::<Molecule>(document).is_err());
    }
}

// Use the same independent, fail-closed schema interpreter as the corpus gate.
#[path = "../../nazuchem-n2s/tests/support/ir_schema.rs"]
mod schema;

#[test]
fn schema_and_both_deserializers_agree_on_representative_wire_shapes() {
    let validator =
        schema::Validator::new(serde_json::from_str(nazuchem_ir::JSON_SCHEMA).unwrap()).unwrap();
    let tetra = serde_json::to_value(tetrahedral()).unwrap();
    let axis = serde_json::to_value(axis()).unwrap();
    let alkene_document = serde_json::to_value(alkene()).unwrap();
    let grouped = {
        let mut m = tetrahedral();
        let Stereo::Tetrahedral { winding, .. } = &mut m.stereo[0] else {
            unreachable!()
        };
        *winding = Winding::Anticlockwise;
        m.stereo_groups.push(StereoGroup {
            kind: StereoGroupKind::Or(2),
            stereo: vec![0],
        });
        serde_json::to_value(m).unwrap()
    };
    let polymer = {
        let mut m = Molecule::new(vec![carbon(0)], vec![]).unwrap();
        m.polymers.push(RepeatUnit {
            atoms: vec![0],
            crossing_bonds: vec![],
            attachment_points: vec![],
            connectivity: Connectivity::Unknown,
            label: None,
        });
        serde_json::to_value(m).unwrap()
    };
    let pv = {
        let mut m = alkene();
        m.position_variation.push(PositionVariation {
            bond: 1,
            endpoints: vec![1, 2],
        });
        serde_json::to_value(m).unwrap()
    };
    let mut integral_float = tetra.clone();
    integral_float["atoms"][0]["hydrogens"] = json!(1.0);
    integral_float["atoms"][0]["charge"] = json!(-0.0);
    integral_float["stereo"][0]["neighbors"][0] = json!(1.0);
    let mut documents = vec![
        (integral_float, true),
        (json!({"nazuchem_ir":1,"atoms":[],"bonds":[]}), true),
        (tetra.clone(), true),
        (axis.clone(), true),
        (alkene_document.clone(), true),
        (grouped.clone(), true),
        (polymer.clone(), true),
        (pv.clone(), true),
        (json!([1, [], [], [], [], [], [], []]), false),
    ];
    for (base, pointer, invalid) in [
        (
            &tetra,
            "/atoms/0",
            json!([6, null, 0, 1, false, 0, null, []]),
        ),
        (&tetra, "/bonds/0", json!([[0, 1], "single"])),
        (
            &tetra,
            "/stereo/0",
            json!(["tetrahedral", 0, [1, 2, 3, "implicit_hydrogen"], "unknown"]),
        ),
        (&grouped, "/stereo_groups/0", json!(["or", 2, [0]])),
        (&polymer, "/polymers/0", json!([[0], [], "unknown", null])),
        (&pv, "/position_variation/0", json!([1, [1, 2]])),
        (
            &tetra,
            "/stereo/0/neighbors/3",
            json!({"kind":"implicit_hydrogen","atom":null}),
        ),
        (
            &tetra,
            "/stereo/0/neighbors/3",
            json!({"kind":"lone_pair","atom":null}),
        ),
        (&tetra, "/stereo/0/neighbors/3", json!({"lone_pair":null})),
        (
            &tetra,
            "/stereo/0/neighbors/3",
            json!({"implicit_hydrogen":null}),
        ),
        (&tetra, "/stereo/0/neighbors/0", json!(-1)),
        (&tetra, "/stereo/0/neighbors/0", json!(4294967296_u64)),
        (&tetra, "/stereo/0/neighbors/0", json!(1.5)),
        (&tetra, "/stereo/0/neighbors/0", json!("future")),
        (
            &grouped,
            "/stereo_groups/0",
            json!({"kind":"abs","id":null,"stereo":[0]}),
        ),
        (
            &grouped,
            "/stereo_groups/0",
            json!({"kind":{"kind":"or","id":2},"stereo":[0]}),
        ),
        (
            &tetra,
            "/stereo/0",
            json!({"kind":"tetrahedral","center":0,"neighbors":[1,2,3,"implicit_hydrogen"],"winding":"unknown","known":false}),
        ),
        (
            &alkene_document,
            "/stereo/0",
            json!({"kind":"double_bond","bond":0,"refs":[2,3],"relation":"trans","neighbors":null}),
        ),
        (
            &axis,
            "/stereo/0",
            json!({"kind":"atropisomer","bond":0,"refs":[2,3,4,5],"winding":"clockwise","relation":null}),
        ),
        (&tetra, "/bonds/0/order", json!({"single":null})),
        (&tetra, "/stereo/0/winding", json!({"unknown":null})),
        (&alkene_document, "/stereo/0/relation", json!({"cis":null})),
        (&axis, "/stereo/0/winding", json!({"clockwise":null})),
        (
            &polymer,
            "/polymers/0/connectivity",
            json!({"unknown":null}),
        ),
        (&tetra, "/bonds/0/order", json!("dative")),
        (&tetra, "/bonds/0/order", json!("zero")),
        (&tetra, "/nazuchem_ir", json!(2)),
    ] {
        let mut document = base.clone();
        *document.pointer_mut(pointer).unwrap() = invalid;
        documents.push((document, false));
    }
    let mut unknown_field = tetra.clone();
    unknown_field["coordinates"] = json!([]);
    documents.push((unknown_field, false));
    for (document, expected) in documents {
        let schema = validator.validate(&document).is_ok();
        let api = Molecule::from_json(&document.to_string()).is_ok();
        let serde = serde_json::from_value::<Molecule>(document.clone()).is_ok();
        assert_eq!(schema, expected, "schema: {document}");
        assert_eq!(api, schema, "from_json: {document}");
        assert_eq!(serde, schema, "serde: {document}");
    }
}

#[test]
fn standalone_struct_deserializers_also_require_objects_and_exact_variant_fields() {
    assert!(serde_json::from_value::<Atom>(json!([6, null, 0, 4, false, 0, null, []])).is_err());
    assert!(serde_json::from_value::<Bond>(json!([[0, 1], "single"])).is_err());
    assert!(serde_json::from_value::<RepeatUnit>(json!([[0], [], "unknown", null])).is_err());
    assert!(serde_json::from_value::<PositionVariation>(json!([0, [1]])).is_err());
    assert!(serde_json::from_value::<StereoGroup>(json!(["and", 7, [0]])).is_err());
    assert!(serde_json::from_value::<StereoGroupKind>(json!(["and", 7])).is_err());
    assert!(
        serde_json::from_value::<StereoGroupKind>(json!({"kind":"and","id":7,"stereo":null}))
            .is_err()
    );
    assert!(
        serde_json::from_value::<Stereo>(
            json!({"kind":"double_bond","bond":0,"refs":[2,3],"relation":"trans","winding":null})
        )
        .is_err()
    );
}

#[test]
fn integer_spellings_accept_only_exact_values_within_the_identity_space() {
    assert_eq!(
        serde_json::from_value::<Ref>(json!(1.0)).unwrap(),
        Ref::Atom(1)
    );
    assert_eq!(
        serde_json::from_value::<Ref>(json!(-0.0)).unwrap(),
        Ref::Atom(0)
    );
    assert_eq!(
        serde_json::from_value::<Ref>(json!(4294967295.0)).unwrap(),
        Ref::Atom(u32::MAX)
    );
    for value in [
        json!(-1.0),
        json!(1.5),
        json!(4294967296.0),
        json!(18446744073709551616.0),
        json!(1e100),
    ] {
        assert!(serde_json::from_value::<Ref>(value).is_err());
    }
    assert_eq!(
        Molecule::from_json(r#"{"nazuchem_ir":1e0,"atoms":[],"bonds":[]}"#).unwrap(),
        Molecule::default()
    );
    assert!(matches!(
        Molecule::from_json(r#"{"nazuchem_ir":2.0,"future":true}"#),
        Err(JsonError::UnsupportedVersion { version: 2 })
    ));
    assert!(
        Molecule::from_json(r#"{"nazuchem_ir":18446744073709551616.0,"atoms":[],"bonds":[]}"#)
            .is_err()
    );
}
