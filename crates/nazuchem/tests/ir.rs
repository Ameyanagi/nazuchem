// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | crates/nazuchem/tests/ir.rs | Phase B step 8a
// Source: Nazuchem | crates/nazuchem/tests/ir.rs | Phase B step 8b

use nazuchem::{AxialWinding, CxsmilesError, IR_VERSION, JSON_SCHEMA, Molecule};

#[cfg(feature = "n2s")]
use nazuchem::{n2s, name_to_molecule, name_to_smiles};

#[cfg(any(feature = "n2s", feature = "s2n"))]
use nazuchem::ErrorKind;

#[cfg(feature = "s2n")]
use nazuchem::{
    Atom, Bond, BondOrder, Budget, IsotopeStyle, molecule_to_name, s2n, smiles_to_name,
};

#[test]
fn ir_contract_is_available_without_engine_features() {
    let schema: serde_json::Value = serde_json::from_str(JSON_SCHEMA).unwrap();
    assert_eq!(schema["properties"]["nazuchem_ir"]["const"], IR_VERSION);
    assert_eq!(
        serde_json::to_value(AxialWinding::Clockwise).unwrap(),
        "clockwise"
    );
    let output: Result<String, CxsmilesError> = Molecule::default().to_cxsmiles();
    assert_eq!(output.unwrap(), "");
}

#[cfg(feature = "n2s")]
#[test]
fn facade_returns_ethanol_ir_with_parent_locants_and_json_roundtrip() {
    let molecule = name_to_molecule("ethanol").unwrap();
    molecule.validate().unwrap();
    assert_eq!(molecule.atoms.len(), 3);
    assert_eq!(molecule.bonds.len(), 2);
    assert_eq!(
        molecule.displayed_locants(),
        [Some("1".into()), Some("2".into()), None]
    );
    assert_eq!(
        molecule.to_cxsmiles().unwrap().split(" |").next().unwrap(),
        name_to_smiles("ethanol").unwrap()
    );
    let restored = Molecule::from_json(&molecule.to_json().unwrap()).unwrap();
    assert_eq!(restored, molecule);
    assert_eq!(
        restored.to_cxsmiles().unwrap(),
        molecule.to_cxsmiles().unwrap()
    );
}

#[cfg(feature = "n2s")]
#[test]
fn strict_ir_facade_preserves_warning_and_parse_error_classification() {
    let parser = n2s::Parser::new().unwrap();
    for (name, status, kind) in [
        ("(+)-lactic acid", n2s::Status::Warning, ErrorKind::Warning),
        (
            "not a chemical name",
            n2s::Status::Failure,
            ErrorKind::Parse,
        ),
    ] {
        let parsed = parser.parse(name, &n2s::ParseOptions::strict());
        assert_eq!(parsed.status, status);
        let error = name_to_molecule(name).unwrap_err();
        assert_eq!(error.kind(), kind);
        assert_eq!(error.to_string(), parsed.message);
    }
}

#[cfg(feature = "n2s")]
#[test]
fn axial_ir_is_available_through_the_facade() {
    let name = "(1M)-[1,1'-binaphthalene]-2,2'-diol";
    let parsed = n2s::Parser::new()
        .unwrap()
        .parse(name, &n2s::ParseOptions::strict());
    assert_eq!(parsed.status, n2s::Status::Success, "{}", parsed.message);
    let molecule = name_to_molecule(name).unwrap();
    assert_eq!(molecule, parsed.structure.unwrap().to_ir().unwrap());
    assert!(
        molecule
            .stereo
            .iter()
            .any(|stereo| matches!(stereo, nazuchem::Stereo::Atropisomer { .. }))
    );
    assert_eq!(
        Molecule::from_json(&molecule.to_json().unwrap()).unwrap(),
        molecule
    );
}

#[cfg(feature = "s2n")]
fn ethanol_ir() -> Molecule {
    Molecule {
        atoms: [(6, 3), (6, 2), (8, 1)]
            .into_iter()
            .map(|(element, hydrogens)| Atom {
                element,
                hydrogens,
                ..Atom::default()
            })
            .collect(),
        bonds: vec![
            Bond {
                atoms: [0, 1],
                order: BondOrder::Single,
            },
            Bond {
                atoms: [1, 2],
                order: BondOrder::Single,
            },
        ],
        ..Molecule::default()
    }
}

#[cfg(feature = "s2n")]
#[test]
fn s2n_only_facade_names_ir_and_preserves_input_annotations() {
    let mut molecule = ethanol_ir();
    // Ordinary maps have no stereo meaning, including values reserved by the
    // legacy string transport. Locants do not constrain generated numbering.
    molecule.atoms[0].atom_map = Some(900_001);
    molecule.atoms[1].locants = vec!["1".into(), "alpha".into()];
    molecule.parent = vec![0, 1];
    molecule.validate().unwrap();
    let original = molecule.clone();

    let name = molecule_to_name(&molecule).unwrap();
    assert_eq!(name.text, "ethanol");
    assert_eq!(
        serde_json::to_value(name).unwrap(),
        serde_json::to_value(smiles_to_name("CCO").unwrap()).unwrap()
    );
    assert_eq!(molecule, original);
}

#[cfg(feature = "s2n")]
#[test]
fn ir_engine_reexport_accepts_style_and_facade_uses_iupac() {
    let mut molecule = ethanol_ir();
    molecule.atoms[0].isotope = Some(13);
    molecule.atoms[0].hydrogens_stated = true;
    let budget = Budget::default();
    let facade = molecule_to_name(&molecule).unwrap();
    let iupac = s2n::name_molecule(&molecule, IsotopeStyle::Iupac, &budget).unwrap();
    let boughton = s2n::name_molecule(&molecule, IsotopeStyle::Boughton, &budget).unwrap();
    assert_eq!(
        serde_json::to_value(facade).unwrap(),
        serde_json::to_value(iupac).unwrap()
    );
    assert_eq!(
        serde_json::to_value(boughton).unwrap(),
        serde_json::to_value(
            s2n::name_smiles_with_style("[13CH3]CO", IsotopeStyle::Boughton, &budget).unwrap()
        )
        .unwrap()
    );
}

#[cfg(feature = "s2n")]
#[test]
fn ir_facade_preserves_admission_errors_as_naming_errors() {
    let mut invalid_bond = ethanol_ir();
    invalid_bond.bonds[0].atoms[1] = 99;
    let mut unsupported_label = ethanol_ir();
    unsupported_label.atoms[0].label = Some("R".into());
    for molecule in [Molecule::default(), invalid_bond, unsupported_label] {
        let expected =
            s2n::name_molecule(&molecule, IsotopeStyle::Iupac, &Budget::default()).unwrap_err();
        let error = molecule_to_name(&molecule).unwrap_err();
        assert_eq!(error.kind(), ErrorKind::Naming);
        assert_eq!(error.to_string(), expected);
    }
}

#[cfg(feature = "s2n")]
#[test]
fn ir_facade_rejects_nonbiaryl_axes_even_when_configuration_is_unknown() {
    // Oxalyl fluoride has sp2 axis ends with distinct O/F substituents, but
    // this acyclic bond is outside the naming engine's biaryl stereo scope.
    let mut molecule = Molecule {
        atoms: [6, 6, 8, 9, 8, 9]
            .into_iter()
            .map(|element| Atom {
                element,
                ..Atom::default()
            })
            .collect(),
        bonds: [
            ([0, 1], BondOrder::Single),
            ([0, 2], BondOrder::Double),
            ([0, 3], BondOrder::Single),
            ([1, 4], BondOrder::Double),
            ([1, 5], BondOrder::Single),
        ]
        .into_iter()
        .map(|(atoms, order)| Bond { atoms, order })
        .collect(),
        ..Molecule::default()
    };
    for winding in [AxialWinding::Clockwise, AxialWinding::Unknown] {
        molecule.stereo = vec![nazuchem::Stereo::Atropisomer {
            bond: 0,
            refs: [
                nazuchem::Ref::Atom(2),
                nazuchem::Ref::Atom(3),
                nazuchem::Ref::Atom(4),
                nazuchem::Ref::Atom(5),
            ],
            winding,
        }];
        molecule.validate().unwrap();
        let error = molecule_to_name(&molecule).unwrap_err();
        assert_eq!(error.kind(), ErrorKind::Naming);
        assert!(error.to_string().contains("biaryl axis stereochemistry"));
    }
}

#[cfg(feature = "s2n")]
#[test]
fn ir_unknown_units_preserve_legacy_admission_without_a_stereogenicity_requirement() {
    let tetra = Molecule {
        atoms: [(6, 3), (6, 1), (6, 3), (8, 1)]
            .into_iter()
            .map(|(element, hydrogens)| Atom {
                element,
                hydrogens,
                ..Atom::default()
            })
            .collect(),
        bonds: [[0, 1], [1, 2], [1, 3]]
            .into_iter()
            .map(|atoms| Bond {
                atoms,
                order: BondOrder::Single,
            })
            .collect(),
        stereo: vec![nazuchem::Stereo::Tetrahedral {
            center: 1,
            neighbors: [
                nazuchem::Ref::Atom(0),
                nazuchem::Ref::Atom(2),
                nazuchem::Ref::Atom(3),
                nazuchem::Ref::ImplicitHydrogen,
            ],
            winding: nazuchem::Winding::Unknown,
        }],
        ..Molecule::default()
    };
    let double = Molecule {
        atoms: [3, 1, 0, 3, 3]
            .into_iter()
            .map(|hydrogens| Atom {
                element: 6,
                hydrogens,
                ..Atom::default()
            })
            .collect(),
        bonds: [
            ([0, 1], BondOrder::Single),
            ([1, 2], BondOrder::Double),
            ([2, 3], BondOrder::Single),
            ([2, 4], BondOrder::Single),
        ]
        .into_iter()
        .map(|(atoms, order)| Bond { atoms, order })
        .collect(),
        stereo: vec![nazuchem::Stereo::DoubleBond {
            bond: 1,
            refs: [0, 3],
            relation: nazuchem::Relation::Unknown,
        }],
        ..Molecule::default()
    };
    let mut cumulene = double.clone();
    cumulene.atoms[3].hydrogens = 1;
    cumulene.bonds[2].order = BondOrder::Double;
    cumulene.bonds[3].atoms = [3, 4];
    for (molecule, smiles) in [
        (tetra, "C[CH:900001](C)O"),
        (double, "C[CH:900002]=[C:900002](C)C"),
        (cumulene, "C[CH:900002]=[C:900002]=CC"),
    ] {
        molecule.validate().unwrap();
        assert_eq!(
            serde_json::to_value(molecule_to_name(&molecule).unwrap()).unwrap(),
            serde_json::to_value(smiles_to_name(smiles).unwrap()).unwrap(),
            "{smiles}"
        );
    }
}

#[cfg(all(feature = "n2s", feature = "s2n"))]
#[test]
fn ir_unknown_salts_preserve_legacy_rejections_without_dropping_components() {
    let mut salt = name_to_molecule("sodium (2R)-2-fluoropropanoate").unwrap();
    let winding = salt
        .stereo
        .iter_mut()
        .find_map(|stereo| match stereo {
            nazuchem::Stereo::Tetrahedral { winding, .. } => Some(winding),
            _ => None,
        })
        .unwrap();
    *winding = nazuchem::Winding::Unknown;
    for with_neutral_component in [false, true] {
        let mut molecule = salt.clone();
        let mut smiles = "C[CH:900001](F)C(=O)[O-].[Na+]".to_owned();
        if with_neutral_component {
            molecule.atoms.push(Atom {
                element: 6,
                hydrogens: 4,
                ..Atom::default()
            });
            smiles.push_str(".C");
        }
        molecule.validate().unwrap();
        let original = molecule.clone();
        let expected = smiles_to_name(&smiles).unwrap_err();
        let error = molecule_to_name(&molecule).unwrap_err();
        assert_eq!(error.kind(), ErrorKind::Naming);
        assert_eq!(error, expected, "{smiles}");
        if with_neutral_component {
            // Admission of the unspecified mixture fails before the later
            // connectedness check for rendering unknown stereo descriptors.
            assert!(
                error
                    .to_string()
                    .contains("neutral components and mixtures")
            );
        } else {
            assert_eq!(
                error.to_string(),
                "Openclatura requires a single connected molecule"
            );
        }
        assert_eq!(molecule, original, "{smiles}");
    }
}

#[cfg(all(feature = "n2s", feature = "s2n"))]
#[test]
fn name_ir_name_facade_composes_for_localized_aromatic_and_chiral_inputs() {
    for name in ["benzene", "(2R)-2-chlorobutane"] {
        let molecule = name_to_molecule(name).unwrap();
        let restored = Molecule::from_json(&molecule.to_json().unwrap()).unwrap();
        assert_eq!(
            serde_json::to_value(molecule_to_name(&restored).unwrap_or_else(|error| {
                panic!("IR naming of {name} failed: {error}; IR {restored:?}")
            }))
            .unwrap(),
            serde_json::to_value(
                smiles_to_name(&name_to_smiles(name).unwrap())
                    .unwrap_or_else(|error| { panic!("SMILES naming of {name} failed: {error}") })
            )
            .unwrap(),
            "{name}"
        );
    }
}

#[cfg(feature = "n2s")]
#[test]
fn large_fusion_letters_and_single_large_components_do_not_panic_or_overrun_the_limit() {
    // An implicit fusion edge beyond 'z' used to unwrap an invalid char.
    assert!(name_to_smiles("benzo[55202]annulene").is_err());
    // A count of one copies nothing, so a large component must not hit the
    // copy-expansion limit (36,001 atoms alone are accepted).
    let alone = name_to_molecule("tetrakis(nonalian-1-yl)methane").unwrap();
    let paired =
        name_to_molecule("tetrakis(nonalian-1-yl)methane compound with methane (1:1)").unwrap();
    assert_eq!(paired.atoms.len(), alone.atoms.len() + 1); // methane: one C, hydrogens folded
}

#[cfg(feature = "n2s")]
#[test]
fn extreme_name_numbers_return_parse_errors_on_every_pointer_width() {
    // Native tests exercise 64-bit arithmetic; the Wasm binding suite runs the
    // same boundaries on wasm32. Every call must return, without unwinding.
    for number in [
        "2147483647",
        "4294967295",
        "18446744073709551615",
        "999999999999999999999999999999999999999",
    ] {
        for name in [
            format!("spiro[{number}.5]pentane"),
            format!("spiro[{number}.{number}.5]pentane"),
            format!("bicyclo[{number}.2.1]heptane"),
            format!("[{number}]annulene"),
            format!("{number}methane"),
            format!("{number}-methylpentane"),
            format!("[1,{number}'-biphenyl]"),
            format!("methanol compound with ethanol (1:{number})"),
            format!("lambda^{number}-phosphane"),
        ] {
            let smiles = name_to_smiles(&name).unwrap_err();
            let molecule = name_to_molecule(&name).unwrap_err();
            assert_eq!(smiles.kind(), ErrorKind::Parse, "{name}: {smiles}");
            assert_eq!(molecule.kind(), ErrorKind::Parse, "{name}: {molecule}");
        }
    }
    assert_eq!(name_to_smiles("ethanol").unwrap(), "C(C)O");
    assert_eq!(name_to_molecule("ethanol").unwrap().atoms.len(), 3);
    assert_eq!(name_to_smiles("lambda^5-phosphane").unwrap(), "[PH5]");
    for name in [
        "lambda^2147483647-phosphane",
        "methanol compound with ethanol (1:2147483647)",
    ] {
        assert!(
            name_to_smiles(name)
                .unwrap_err()
                .to_string()
                .contains("expansion limit")
        );
        assert!(
            name_to_molecule(name)
                .unwrap_err()
                .to_string()
                .contains("expansion limit")
        );
    }
}
