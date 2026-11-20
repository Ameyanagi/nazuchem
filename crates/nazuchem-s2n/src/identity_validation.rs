// SPDX-License-Identifier: MIT AND (MIT OR Apache-2.0)
// Origin: mixed
// Source: Nazuchem | crates/nazuchem-s2n/src/identity_validation.rs | amide nitrogen charge follow-up
// Source: Nazuchem | crates/nazuchem-s2n/src/identity_validation.rs | issue #186
// Source: ReShiki | src/naming/openclatura/identity_validation.rs | 327d32fed63a5b3df15b7b8dc5ab84a7c36e0865
// Copyright (c) 2026 lamalab-org
// Source: openclatura | tests/data/parser_xml_resources/heteroAtoms.json | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/fixtures/diverse_corpus.csv | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/unit/test_analysis.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/unit/test_chalcogen_functional_groups.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/data/parser_xml_resources/wordRules.json | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/data/parser_xml_resources/groupStemsAllowingInlineSuffixes.json | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Copyright (c) 2026 Ameyanagi and ReShiki contributors

//! Input states that the plain-name renderers cannot preserve.
use super::molecule::Molecule;
use super::perception::PerceivedGroup;
use std::collections::BTreeSet;

// Provenance: fn validate_principal_charge_roles | original
// Source: ReShiki | src/naming/openclatura/identity_validation.rs::validate_principal_charge_roles
pub fn validate_principal_charge_roles(
    mol: &Molecule,
    groups: &[PerceivedGroup],
    key: Option<&str>,
) -> Result<(), String> {
    let selected = || groups.iter().filter(|g| Some(g.key.as_str()) == key);
    let has_derivative = |derivative: &str| {
        selected().any(|g| {
            g.descriptor
                .as_ref()
                .is_some_and(|d| d.derivative == derivative)
        })
    };
    // Upstream divergence: a shared ester/anion suffix cannot assign separate charge states to mixed principal groups.
    if has_derivative("anion") && has_derivative("ester") {
        return Err(
            "openclatura: mixed ester and anion principal citations are not supported".into(),
        );
    }
    // P-72.2.2.2.4: a neutral suffix such as amide would drop the charge of
    // its nitrogen (acetamide for acetylazanide).
    if selected().any(|g| {
        g.is_principal_candidate
            && !super::perception::is_anionic_group_key(&g.key)
            && g.atoms_involved
                .iter()
                .any(|i| mol.atoms[i].symbol == "N" && mol.atoms[i].charge < 0)
            && g.atom_ids()
                .iter()
                .map(|i| mol.atoms[i].charge)
                .sum::<i32>()
                < 0
    }) {
        return Err("openclatura: anionic nitrogen in a neutral suffix is not supported".into());
    }
    Ok(())
}

// Provenance: fn validate_input | original
// Source: ReShiki | src/naming/openclatura/identity_validation.rs::validate_input
pub fn validate_input(mol: &Molecule) -> Result<(), String> {
    for (&idx, atom) in &mol.atoms {
        // Upstream divergence: elemental ion shortcuts require their exact charge and hydrogen state.
        if mol.atoms.len() == 1
            && atom.charge != 0
            && !matches!(
                (atom.symbol.as_str(), atom.charge, atom.total_h_count),
                ("Li" | "Na" | "K", 1, 0)
                    | ("Mg" | "Ca", 2, 0)
                    | ("F" | "Cl" | "Br" | "I", -1, 0)
                    | ("N", 1, 4)
                    | ("N", -1, 2)
                    | ("O" | "S", -1, 1)
                    | ("O" | "S" | "Se", -2, 0)
                    | ("O", 1, 3)
                    | ("P", 1, 4)
                    | ("H", -1 | 1, 0)
                    | ("C" | "Si", 1, 3)
            )
        {
            return Err(format!(
                "openclatura: monatomic {} charge/hydrogen state is not supported",
                atom.symbol
            ));
        }
        // Upstream divergence: unit-charge names cannot represent unsupported charge magnitudes.
        if atom.charge.unsigned_abs() > 1
            && !(mol.atoms.len() == 1
                && atom.total_h_count == 0
                && matches!(
                    (atom.symbol.as_str(), atom.charge),
                    ("Mg" | "Ca", 2) | ("O" | "S" | "Se", -2)
                ))
        {
            return Err(format!(
                "openclatura: charge {} on {} atom {idx} is not supported",
                atom.charge, atom.symbol
            ));
        }
        if atom.charge != 0
            || atom.total_h_count == 0
            || !matches!(
                atom.symbol.as_str(),
                "P" | "As" | "Sb" | "Bi" | "S" | "Se" | "Te"
            )
            || mol.degree(idx) != 1
        {
            continue;
        }
        let neighbor = mol.get_neighbors(idx)[0];
        let bond = mol
            .get_bond(idx, neighbor)
            .ok_or("openclatura: missing terminal hydride bond")?;
        let expected = atom
            .element()
            .ok_or("openclatura: unsupported terminal hydride element")?
            .standard_valence;
        // Upstream divergence: ordinary terminal hydride/group names require their exact bonding number and H count.
        if (bond.order == 1 || mol.atoms[&neighbor].is_carbon())
            && bond.order + atom.total_h_count != expected
        {
            return Err(format!(
                "openclatura: terminal {} atom {idx} with {} hydrogens and bonding number {} is not supported",
                atom.symbol,
                atom.total_h_count,
                bond.order + atom.total_h_count
            ));
        }
    }
    Ok(())
}

/// Upstream divergence: CIP implementations disagree on the duplicate atom of
/// a mancude ring bond (RDKit's labels and the OPSIN reading of the same name
/// differ). Decline a stereocentre whose ligand ranking, compared sphere by
/// sphere, is decided or tied where that duplicate's element is ambiguous.
// Provenance: fn validate_mancude_duplicate_ranking | original
// Source: ReShiki | src/naming/openclatura/identity_validation.rs::validate_mancude_duplicate_ranking
pub fn validate_mancude_duplicate_ranking(mol: &Molecule) -> Result<(), String> {
    const SPHERES: usize = 8;
    let z = |idx: usize| mol.atoms[&idx].element().map_or(0, |e| e.atomic_number);
    // Per sphere, descending atomic numbers with each ambiguous duplicate at
    // its lowest and at its highest element.
    let spheres = |center: usize, ligand: usize| -> Vec<[Vec<usize>; 2]> {
        let mut seen = BTreeSet::from([center, ligand]);
        let mut frontier = vec![(ligand, center)];
        let mut result = Vec::new();
        for _ in 0..SPHERES {
            let mut sphere: [Vec<usize>; 2] = Default::default();
            let mut next = Vec::new();
            for &(atom_idx, parent) in &frontier {
                let atom = &mol.atoms[&atom_idx];
                let mut entries = vec![1; atom.total_h_count];
                let mut aromatic_double = false;
                for &n in mol.get_neighbors(atom_idx) {
                    let order = mol.get_bond(atom_idx, n).map_or(1, |b| b.order);
                    let aromatic = atom.is_aromatic && mol.atoms[&n].is_aromatic;
                    aromatic_double |= aromatic && order == 2;
                    let duplicates = if aromatic { 0 } else { order - 1 };
                    entries.extend(std::iter::repeat_n(z(n), duplicates));
                    if n == parent {
                        continue;
                    }
                    entries.push(z(n));
                    if seen.insert(n) {
                        next.push((n, atom_idx));
                    }
                }
                let variants: Vec<_> = mol
                    .get_neighbors(atom_idx)
                    .iter()
                    .filter(|&&n| aromatic_double && mol.atoms[&n].is_aromatic)
                    .map(|&n| z(n))
                    .collect();
                for (side, duplicate) in sphere
                    .iter_mut()
                    .zip([variants.iter().min(), variants.iter().max()])
                {
                    side.extend(&entries);
                    side.extend(duplicate);
                }
            }
            for side in &mut sphere {
                side.sort_unstable_by(|a, b| b.cmp(a));
            }
            result.push(sphere);
            frontier = next;
        }
        result
    };
    for (&center, atom) in &mol.atoms {
        if atom.stereo.is_none() {
            continue;
        }
        let ligands: Vec<_> = mol
            .get_neighbors(center)
            .iter()
            .map(|&l| (z(l), spheres(center, l)))
            .collect();
        for (i, (zi, left)) in ligands.iter().enumerate() {
            for (zj, right) in &ligands[i + 1..] {
                if zi != zj {
                    continue;
                }
                for (l, r) in left.iter().zip(right) {
                    let orders: BTreeSet<_> = l
                        .iter()
                        .flat_map(|l| r.iter().map(move |r| l.cmp(r)))
                        .collect();
                    if orders.len() > 1 {
                        return Err(format!(
                            "openclatura: CIP rank at atom {center} depends on a mancude duplicate atom"
                        ));
                    }
                    if !orders.contains(&std::cmp::Ordering::Equal) {
                        break;
                    }
                }
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::name_smiles;

    #[test]
    // Provenance: fn tests::charged_characteristic_groups_preserve_identity | original
    // Source: Nazuchem | crates/nazuchem-s2n/src/identity_validation.rs::tests::charged_characteristic_groups_preserve_identity
    fn charged_characteristic_groups_preserve_identity() {
        let budget = Default::default();
        for (smiles, expected) in [
            (
                "CC[N+](C)=C=S",
                "N-ethyl-N-methyl-1-sulfanylidenemethaniminium",
            ),
            ("CC[N+](C)=C=O", "N-ethyl-N-methyl-1-oxomethaniminium"),
            (
                "C1CCC[N+]1=C=S",
                "1-(sulfanylidenemethylidene)pyrrolidin-1-ium",
            ),
            (
                "[O-]C(=O)CC[N+](C)=C=S",
                "3-[methyl(sulfanylidenemethylidene)ammonio]propanoate",
            ),
            (
                "C[NH+](C)[N+](C)=C=S",
                "N-methyl-N-[methyl(sulfanylidenemethylidene)ammonio]methanaminium",
            ),
            (
                "CC[N+](C)=C=S.[Cl-]",
                "N-ethyl-N-methyl-1-sulfanylidenemethaniminium chloride",
            ),
            ("CCN=C=[OH+]", "N-ethyl-1-oxidaniumylidenemethanimine"),
            ("CCN=C=S", "isothiocyanatoethane"),
            ("CNNC(C)=O", "N'-methylacetohydrazide"),
            ("[NH-]C(=O)C", "acetylazanide"),
        ] {
            let ir = crate::ir_input::from_smiles(smiles, &budget).unwrap();
            let result = [
                name_smiles(smiles, &budget),
                crate::name_molecule(&ir, nazuchem_ir::IsotopeStyle::Iupac, &budget),
            ];
            for name in result {
                assert_eq!(name.unwrap().text, expected, "{smiles}");
            }
            let parsed = crate::test_support::parse_strict(expected).unwrap();
            crate::test_support::verify_identity(smiles, &parsed).unwrap();
        }
        for smiles in ["CN[N-]C(C)=O", "CCN=C=[SH+]", "C[NH+](C)N=C=[SH+]"] {
            let ir = crate::ir_input::from_smiles(smiles, &budget).unwrap();
            assert!(name_smiles(smiles, &budget).is_err(), "{smiles}");
            assert!(
                crate::name_molecule(&ir, nazuchem_ir::IsotopeStyle::Iupac, &budget).is_err(),
                "{smiles}"
            );
        }
    }

    #[test]
    // Provenance: fn tests::amidium_names_preserve_charge_in_both_public_naming_paths | original
    // Source: Nazuchem | crates/nazuchem-s2n/src/identity_validation.rs::tests::amidium_names_preserve_charge_in_both_public_naming_paths
    fn amidium_names_preserve_charge_in_both_public_naming_paths() {
        use crate::test_support::{canonical_smiles, parse_strict};

        let budget = Default::default();
        for (smiles, expected) in [
            // P-73.1.2.1 and Table 7.4: use amidium for positive nitrogen.
            ("C[NH+](C)C=S", "N,N-dimethylmethanethioamidium"),
            (
                "C[NH+](C)C=S.[Cl-]",
                "N,N-dimethylmethanethioamidium chloride",
            ),
            ("C[NH+](C)C=O", "N,N-dimethylformamidium"),
            ("C[NH+](C)C=O.[Cl-]", "N,N-dimethylformamidium chloride"),
            (
                "O=C(O)CC=[N+](C)C=S",
                "N-(3-hydroxy-3-oxopropylidene)-N-methylmethanethioamidium",
            ),
            (
                "C[NH+](C)C(=S)c1ccccc1",
                "N,N-dimethylbenzenecarbothioamidium",
            ),
            ("C[N+](C)(C)C=S", "N,N,N-trimethylmethanethioamidium"),
            ("C[NH+](C)C=[Se]", "N,N-dimethylmethaneselenoamidium"),
            ("C[NH+](C)C=[Te]", "N,N-dimethylmethanetelluroamidium"),
            ("CC(=[OH+])N", "1-oxidaniumylideneethan-1-amine"),
            (
                "CC(=[OH+])N.[Cl-]",
                "1-oxidaniumylideneethan-1-amine chloride",
            ),
            (
                "C[NH+](C)C(=S)CC(=O)[O-]",
                "3-(N-methylmethanaminiumyl)-3-sulfanylidenepropanoate",
            ),
            (
                "[NH3+]C(=S)CC(=O)[O-]",
                "3-azaniumyl-3-sulfanylidenepropanoate",
            ),
            ("[NH3+]C(=O)CC(=O)[O-]", "3-azaniumyl-3-oxopropanoate"),
            ("[NH3+]C(=S)C(=O)[O-]", "azaniumyl(sulfanylidene)acetate"),
            ("[NH3+]C(=O)C(=O)[O-]", "azaniumyl(oxo)acetate"),
            ("CN(C)C=S", "N,N-dimethylmethanethioamide"),
            ("CC(=S)NN=C(C)C", "N'-(propan-2-ylidene)ethanethiohydrazide"),
            ("C[N+]([O-])(C)C=S", "N,N-dimethyl-N-oxidomethanethioamide"),
        ] {
            let ir = crate::ir_input::from_smiles(smiles, &budget).unwrap();
            for result in [
                name_smiles(smiles, &budget),
                crate::name_molecule(&ir, nazuchem_ir::IsotopeStyle::Iupac, &budget),
            ] {
                let name = result.unwrap_or_else(|error| panic!("{smiles}: {error}"));
                assert_eq!(name.text, expected, "{smiles}");
                let reconstructed = parse_strict(&name.text).unwrap();
                assert_eq!(
                    canonical_smiles(&reconstructed).unwrap(),
                    canonical_smiles(smiles).unwrap(),
                    "{smiles} -> {} -> {reconstructed}",
                    name.text,
                );
            }
        }
    }

    #[test]
    // Provenance: fn tests::charged_amidium_variants_preserve_protons_and_charge_pairs | original
    // Source: Nazuchem | crates/nazuchem-s2n/src/identity_validation.rs::tests::charged_amidium_variants_preserve_protons_and_charge_pairs
    fn charged_amidium_variants_preserve_protons_and_charge_pairs() {
        use crate::test_support::{canonical_smiles, parse_strict};
        use nazuchem_ir::{NameReason, NameStatus};

        let budget = Default::default();
        let fallback = NameStatus::Unverified(NameReason::CationParent);
        for (smiles, expected, status) in [
            // The retained uronium names parse as C=N+ tautomers.
            ("NC(=O)[NH3+]", "aminoformamidium", fallback),
            ("[NH3+]C(=O)N", "aminoformamidium", fallback),
            ("NC(=S)[NH3+]", "aminomethanethioamidium", fallback),
            ("[NH3+]C(=S)N", "aminomethanethioamidium", fallback),
            (
                "[NH3+]C(=O)C(=O)[NH3+]",
                "ethanebis(amidium)",
                NameStatus::Preferred,
            ),
            (
                "[NH3+]C(=S)CC(=S)[NH3+]",
                "propanebis(thioamidium)",
                NameStatus::Preferred,
            ),
            ("[NH2+]([O-])C(=O)C", "N-oxidoacetamidium", fallback),
            ("[NH2+]([S-])C(=S)C", "N-sulfidoethanethioamidium", fallback),
            (
                "[NH2+]([Se-])C(=[Se])C",
                "N-selenidoethaneselenoamidium",
                fallback,
            ),
            (
                "[NH2+]([Te-])C(=[Te])C",
                "N-telluridoethanetelluroamidium",
                fallback,
            ),
            (
                "C[NH+]([O-])C(=O)C",
                "N-methyl-N-oxidoacetamidium",
                fallback,
            ),
        ] {
            let ir = crate::ir_input::from_smiles(smiles, &budget).unwrap();
            let reconstructed = parse_strict(expected).unwrap();
            assert_eq!(
                canonical_smiles(&reconstructed).unwrap(),
                canonical_smiles(smiles).unwrap(),
                "{smiles} -> {expected} -> {reconstructed}",
            );
            for result in [
                name_smiles(smiles, &budget),
                crate::name_molecule(&ir, nazuchem_ir::IsotopeStyle::Iupac, &budget),
            ] {
                let name = result.unwrap_or_else(|error| panic!("{smiles}: {error}"));
                assert_eq!(name.text, expected, "{smiles}");
                assert_eq!(name.status, status, "{smiles}");
            }
        }
    }

    #[test]
    // Provenance: fn tests::amide_nitrogen_chalcogenides_keep_the_neutral_amide_suffix | original
    // Source: Nazuchem | crates/nazuchem-s2n/src/identity_validation.rs::tests::amide_nitrogen_chalcogenides_keep_the_neutral_amide_suffix
    fn amide_nitrogen_chalcogenides_keep_the_neutral_amide_suffix() {
        use crate::test_support::{canonical_smiles, parse_strict};
        use nazuchem_ir::{NameReason, NameStatus};

        let budget = Default::default();
        for (element, prefix) in [
            ("O", "oxido"),
            ("S", "sulfido"),
            ("Se", "selenido"),
            ("Te", "tellurido"),
        ] {
            let expected = format!("N,N,N',N'-tetramethyl-N-{prefix}ethanediamide");
            let reconstructed = parse_strict(&expected).unwrap();
            let canonical = canonical_smiles(&reconstructed).unwrap();
            for smiles in [
                format!("C[N+]([{element}-])(C)C(=O)C(=O)N(C)C"),
                format!("CN(C)C(=O)C(=O)[N+](C)(C)[{element}-]"),
            ] {
                let ir = crate::ir_input::from_smiles(&smiles, &budget).unwrap();
                assert_eq!(canonical_smiles(&smiles).unwrap(), canonical, "{smiles}");
                for result in [
                    name_smiles(&smiles, &budget),
                    crate::name_molecule(&ir, nazuchem_ir::IsotopeStyle::Iupac, &budget),
                ] {
                    let name = result.unwrap_or_else(|error| panic!("{smiles}: {error}"));
                    assert_eq!(name.text, expected, "{smiles}");
                    assert_eq!(
                        name.status,
                        NameStatus::NotPreferred(NameReason::UnretainedParent),
                        "{smiles}",
                    );
                }
            }
        }
    }

    #[test]
    // Provenance: fn tests::multiplied_amidium_n_locants_are_independent_of_atom_order | original
    // Source: Nazuchem | crates/nazuchem-s2n/src/identity_validation.rs::tests::multiplied_amidium_n_locants_are_independent_of_atom_order
    fn multiplied_amidium_n_locants_are_independent_of_atom_order() {
        use crate::test_support::{canonical_smiles, parse_strict};

        let budget = Default::default();
        for (smiles_pair, expected) in [
            (
                ["C[NH2+]C(=O)C(=O)[NH3+]", "[NH3+]C(=O)C(=O)[NH2+]C"],
                "N-methylethanebis(amidium)",
            ),
            (
                ["C[NH2+]C(=S)C(=S)[NH3+]", "[NH3+]C(=S)C(=S)[NH2+]C"],
                "N-methylethanebis(thioamidium)",
            ),
            (
                ["C[NH2+]C(=S)C(=S)[NH2+]CC", "CC[NH2+]C(=S)C(=S)[NH2+]C"],
                "N-ethyl-N'-methylethanebis(thioamidium)",
            ),
            (
                ["C[NH2+]C(=O)C(=S)[NH2+]CC", "CC[NH2+]C(=S)C(=O)[NH2+]C"],
                "2-(ethanaminiumyl)-N-methyl-2-sulfanylideneacetamidium",
            ),
        ] {
            let reconstructed = parse_strict(expected).unwrap();
            let canonical = canonical_smiles(&reconstructed).unwrap();
            for smiles in smiles_pair {
                let ir = crate::ir_input::from_smiles(smiles, &budget).unwrap();
                assert_eq!(canonical_smiles(smiles).unwrap(), canonical, "{smiles}");
                for result in [
                    name_smiles(smiles, &budget),
                    crate::name_molecule(&ir, nazuchem_ir::IsotopeStyle::Iupac, &budget),
                ] {
                    let name = result.unwrap_or_else(|error| panic!("{smiles}: {error}"));
                    assert_eq!(name.text, expected, "{smiles}");
                    assert_eq!(name.status, nazuchem_ir::NameStatus::Preferred, "{smiles}");
                }
            }
        }
    }

    #[test]
    // Provenance: fn tests::shared_nitrogen_amidium_parent_is_independent_of_atom_order | original
    // Source: Nazuchem | crates/nazuchem-s2n/src/identity_validation.rs::tests::shared_nitrogen_amidium_parent_is_independent_of_atom_order
    fn shared_nitrogen_amidium_parent_is_independent_of_atom_order() {
        use crate::test_support::{canonical_smiles, parse_strict};

        let budget = Default::default();
        for (smiles_pair, expected) in [
            (
                ["CC(=O)[NH+](C)C(=O)CC", "CCC(=O)[NH+](C)C(=O)C"],
                "N-acetyl-N-methylpropanamidium",
            ),
            (
                ["CC(=S)[NH+](C)C(=S)CC", "CCC(=S)[NH+](C)C(=S)C"],
                "N-(ethanethioyl)-N-methylpropanethioamidium",
            ),
            (
                ["CC(=O)[NH+](C)C(=S)CC", "CCC(=S)[NH+](C)C(=O)C"],
                "N-methyl-N-(propanethioyl)acetamidium",
            ),
        ] {
            let reconstructed = parse_strict(expected).unwrap();
            let canonical = canonical_smiles(&reconstructed).unwrap();
            for smiles in smiles_pair {
                let ir = crate::ir_input::from_smiles(smiles, &budget).unwrap();
                assert_eq!(canonical_smiles(smiles).unwrap(), canonical, "{smiles}");
                for result in [
                    name_smiles(smiles, &budget),
                    crate::name_molecule(&ir, nazuchem_ir::IsotopeStyle::Iupac, &budget),
                ] {
                    let name = result.unwrap_or_else(|error| panic!("{smiles}: {error}"));
                    assert_eq!(name.text, expected, "{smiles}");
                    assert_eq!(name.status, nazuchem_ir::NameStatus::Preferred, "{smiles}");
                }
            }
        }
    }

    #[test]
    // Provenance: fn tests::amidium_substituent_prefixes_without_a_charge_citation_are_declined | original
    // Source: Nazuchem | crates/nazuchem-s2n/src/identity_validation.rs::tests::amidium_substituent_prefixes_without_a_charge_citation_are_declined
    fn amidium_substituent_prefixes_without_a_charge_citation_are_declined() {
        let budget = Default::default();
        let smiles = "C[NH+](C)C(=S)c1ccc(C(=O)[O-])cc1";
        let ir = crate::ir_input::from_smiles(smiles, &budget).unwrap();
        for result in [
            name_smiles(smiles, &budget),
            crate::name_molecule(&ir, nazuchem_ir::IsotopeStyle::Iupac, &budget),
        ] {
            assert_eq!(
                result.unwrap_err(),
                "openclatura: amidium substituent prefix is not supported",
            );
        }
    }

    #[test]
    // Provenance: fn tests::public_entry_declines_kekule_dependent_cip_ranks | original
    // Source: ReShiki | src/naming/openclatura/identity_validation.rs::tests::public_entry_declines_kekule_dependent_cip_ranks
    fn public_entry_declines_kekule_dependent_cip_ranks() {
        for smiles in [
            "NC(N)[C@@H](C)c1ncncc1",
            "O=[N+]([O-])[C@@H]1NCC[C@@H]1c1ncccc1",
            "O=[N+]([O-])[C@@H]1NCC[C@@H]1c1ncncc1Br",
            "C[C@H](Cc1ncccc1)CC=NC1CCCCC1",
            "C[C@H](Cc1ncccc1)CC=NC1CCC2(CC1)CCCC2",
        ] {
            assert!(
                name_smiles(smiles, &Default::default()).is_err(),
                "{smiles}"
            );
        }
        for (smiles, expected) in [
            (
                "O=[N+]([O-])[C@@H]1NCC[C@@H]1c1ccncc1",
                "4-[(2S,3R)-2-nitropyrrolidin-3-yl]pyridine",
            ),
            (
                "O=[N+]([O-])C[C@@H](C)c1ncccc1",
                "2-[(2R)-1-nitropropan-2-yl]pyridine",
            ),
        ] {
            assert_eq!(
                name_smiles(smiles, &Default::default())
                    .map(|name| name.text)
                    .as_deref(),
                Ok(expected),
                "{smiles}"
            );
        }
    }

    #[test]
    // Provenance: fn tests::public_entry_declines_unrepresented_charge_states | adapted
    // Source: ReShiki | src/naming/openclatura/identity_validation.rs::tests::public_entry_declines_unrepresented_charge_states
    // Source: openclatura | tests/data/parser_xml_resources/heteroAtoms.json::<fixture data>
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn public_entry_declines_unrepresented_charge_states() {
        for smiles in [
            "C[NH+2]C",
            "C[Si-](C)(C)C",
            "[Mg+]",
            "[NaH+]",
            "[ClH-]",
            "COC(=O)C1CC(C(=O)[O-])NC(C(N)=S)C1",
            "CCC[PH-]",
            "CCC[As-](C)",
            "CCC[SH2-]",
            "CCC[TeH2-]",
            "CCC[ClH+]",
            "CCC[IH-]",
        ] {
            assert!(
                name_smiles(smiles, &Default::default()).is_err(),
                "{smiles}"
            );
        }
    }

    #[test]
    // Provenance: fn tests::public_entry_declines_unrepresented_parent_elements_and_hydrogens | original
    // Source: ReShiki | src/naming/openclatura/identity_validation.rs::tests::public_entry_declines_unrepresented_parent_elements_and_hydrogens
    fn public_entry_declines_unrepresented_parent_elements_and_hydrogens() {
        for smiles in [
            "[Tl]1CC1",
            "[InH]1CCCC1",
            "C[SH3]",
            "C[SH5]",
            "C[SeH3]",
            "C[PH4]",
            "CO[SH3]",
            "CC(=O)[SeH3]",
            "[SH2]1CCCC1",
        ] {
            assert!(
                name_smiles(smiles, &Default::default()).is_err(),
                "{smiles}"
            );
        }
    }

    #[test]
    // Provenance: fn tests::public_entry_preserves_supported_charge_and_hydride_names | adapted
    // Source: ReShiki | src/naming/openclatura/identity_validation.rs::tests::public_entry_preserves_supported_charge_and_hydride_names
    // Source: openclatura | tests/fixtures/diverse_corpus.csv::<fixture data>
    // Source: openclatura | tests/unit/test_analysis.py::test_group_13_14_central_atoms_have_a_parent_hydride
    // Source: openclatura | tests/unit/test_chalcogen_functional_groups.py::test_simple_chalcogen_analogues_are_named_from_graph_roles
    // Source: openclatura | tests/data/parser_xml_resources/wordRules.json::<fixture data>
    // Source: openclatura | tests/data/parser_xml_resources/groupStemsAllowingInlineSuffixes.json::<fixture data>
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn public_entry_preserves_supported_charge_and_hydride_names() {
        for (smiles, expected) in [
            ("C", "methane"),
            ("CC", "ethane"),
            ("C[NH2+]C", "N-methylmethanaminium"),
            ("C[N+](C)(C)C", "N,N,N-trimethylmethanaminium"),
            ("CC(=O)[O-]", "acetate"),
            ("C[S-]", "methanethiolate"),
            ("C[SeH]", "methaneselenol"),
            ("CCC[SeH2+]", "propylselanium"),
            ("CS", "methanethiol"),
            ("C[PH2]", "phosphanylmethane"),
            ("C[SH2]C", "dimethyl-lambda4-sulfane"),
            ("[CH2-][S+](C)C", "(dimethylsulfaniumyl)methanide"),
            ("C[Si](C)(C)C", "tetramethylsilane"),
            ("[NH4+]", "ammonium"),
            ("[O-2]", "oxide"),
            ("[Ca+2]", "calcium(2+)"),
            ("[N-]1[NH+]=CC=C1", "1H-pyrazol-2-ium-1-ide"),
        ] {
            assert_eq!(
                name_smiles(smiles, &Default::default())
                    .map(|name| name.text)
                    .as_deref(),
                Ok(expected),
                "{smiles}"
            );
        }
    }
}
