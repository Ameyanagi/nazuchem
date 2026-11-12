// SPDX-License-Identifier: MIT AND (MIT OR Apache-2.0)
// Origin: mixed
// Source: ReShiki | src/naming/openclatura/radical_naming.rs | 5f559fd482e237307a589dfcb7ea333cee0a6e13
// Source: Nazuchem | crates/nazuchem-s2n/src/radical_naming.rs | issue #27
// Copyright (c) 2026 lamalab-org
// Source: openclatura | tests/exploratory/test_parenthesis_functions_direct.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/unit/test_analysis.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/data/parser_xml_resources/simpleSubstituents.json | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/fuzz/test_ligand_parentheses.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/unit/test_chalcogen_functional_groups.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/exploratory/test_charge_representation_gaps.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Copyright (c) 2026 Ameyanagi and ReShiki contributors

//! Neutral radicals whose subtractive P-71 names have exact parser identities.
//!
//! <https://iupac.qmul.ac.uk/BlueBook/P7.html>, P-71.2.1–P-71.3.4.
//! A virtual excluded bond reuses the ordinary substituent renderer without
//! discarding the free-valence order. More elaborate radical characteristic
//! groups and radical ions are deliberately rejected before ordinary naming.
use super::{
    Name,
    formatting::strip_outer_parentheses,
    molecule::{Atom, AtomSet, Bond, Molecule},
};

// Provenance: fn unsupported | original
// Source: ReShiki | src/naming/openclatura/radical_naming.rs::unsupported
fn unsupported() -> String {
    "openclatura: this neutral radical's parent or characteristic group is not supported".into()
}

// Provenance: fn name_radical | original
// Source: ReShiki | src/naming/openclatura/radical_naming.rs::name_radical
pub(super) fn name_radical(mol: &Molecule) -> Result<Option<Name>, String> {
    let centers: Vec<_> = mol
        .atoms
        .values()
        .filter(|atom| atom.radical_electrons != 0)
        .collect();
    if centers.is_empty() {
        return Ok(ylium(mol));
    }
    if mol.atoms.values().any(|atom| atom.charge != 0) {
        return Err("openclatura: radical ions are not supported".into());
    }
    if centers
        .iter()
        .map(|atom| atom.radical_electrons)
        .sum::<usize>()
        > 2
        || mol
            .atoms
            .values()
            .any(|atom| atom.has_stereo() || atom.isotope.is_some())
        || mol
            .bonds
            .values()
            .any(|bond| bond.stereo.is_some() || !(1..=3).contains(&bond.order))
        || super::graph_io::get_connected_components(mol).len() != 1
    {
        return Err(unsupported());
    }
    let text = if centers.len() == 2 {
        carbon_diyl(mol)?
    } else {
        let center = centers[0];
        match center.symbol.as_str() {
            "C" if carbon_parent(mol, center.idx, center.radical_electrons == 1) => {
                substituent_radical(mol, center.idx, center.radical_electrons)?
            }
            // P-71.2.1.2: any position of a parent hydride, so an arene carbon
            // keeps the ring as parent and cites everything else as prefixes.
            "C" if center.radical_electrons == 1 && center.is_aromatic => {
                substituent_radical(mol, center.idx, 1)?
            }
            // P-71.2.1.1 silyl and P-71.2.1.2 boranyl, with substituent prefixes.
            "Si" | "B"
                if center.radical_electrons == 1
                    && mononuclear_center(mol, center.idx)
                    && center.element().is_some_and(|e| {
                        mol.degree(center.idx) + center.total_h_count + 1 == e.standard_valence
                    }) =>
            {
                let name = substituent_radical(mol, center.idx, 1)?;
                // P-71.2.1.2: boranyl, not the boryl prefix.
                match name.strip_suffix("boryl") {
                    Some(stem) => format!("{stem}boranyl"),
                    None => name,
                }
            }
            // P-71.2.1.2 with the lambda-convention (P-14.1): a free valence on a
            // centre two above its standard bonding number, as in
            // tert-butoxytri(phenyl)-lambda5-phosphanyl.
            "P" | "As" | "S" | "Se"
                if center.radical_electrons == 1
                    && mononuclear_center(mol, center.idx)
                    && center.element().is_some_and(|e| {
                        mol.degree(center.idx) + center.total_h_count + 1 == e.standard_valence + 2
                    }) =>
            {
                let name = substituent_radical(mol, center.idx, 1)?;
                // The unsubstituted hydride radical is rendered without its lambda.
                if !name.contains("lambda") {
                    return Err(unsupported());
                }
                name
            }
            "N" => amine_radical(mol, center.idx, center.radical_electrons)?,
            "O" | "S" | "Se" if center.radical_electrons == 1 => {
                chalcogen_radical(mol, center.idx)?
            }
            _ => return Err(unsupported()),
        }
    };
    Ok(Some(Name {
        text: super::nest_enclosing_marks(&text),
        status: mol.name_status(),
        alternatives: vec![],
        omitted_stereo: vec![],
        omitted_units: 0,
        verification: None,
    }))
}

/// A singly bonded centre whose molecule has no second parent-hydride
/// heteroatom: disilanyl, disiloxanyl and the like need their own parent.
// Provenance: fn mononuclear_center | original
// Source: Nazuchem | crates/nazuchem-s2n/src/radical_naming.rs::mononuclear_center
fn mononuclear_center(mol: &Molecule, center: usize) -> bool {
    mol.get_neighbors(center)
        .iter()
        .all(|&n| mol.bond_order(center, n) == 1)
        && mol.atoms.values().all(|atom| {
            atom.idx == center
                || matches!(
                    atom.symbol.as_str(),
                    "C" | "N" | "O" | "S" | "F" | "Cl" | "Br" | "I"
                )
        })
}

/// Only hydrocarbons and a carbonyl at the free carbon are audited here.
// Provenance: fn carbon_parent | original
// Source: ReShiki | src/naming/openclatura/radical_naming.rs::carbon_parent
fn carbon_parent(mol: &Molecule, root: usize, allow_acyl: bool) -> bool {
    let carbon_degree = |id| {
        mol.get_neighbors(id)
            .iter()
            .filter(|n| mol.atoms[n].is_carbon())
            .count()
    };
    let unbranched = mol
        .atoms
        .values()
        .filter(|atom| atom.is_carbon())
        .all(|atom| carbon_degree(atom.idx) <= 2);
    // The only branched parent audited here is the saturated star centered on
    // the free valence (in particular 2-methylpropan-2-yl). Other branching and
    // ring side chains need a separate P-46/P-71 parent-selection audit.
    let centered_star = mol.atoms.values().all(|atom| atom.is_carbon())
        && mol.bonds.values().all(|bond| bond.order == 1)
        && mol
            .atoms
            .keys()
            .all(|&id| id == root || mol.degree(id) == 1)
        && mol.degree(root) + 1 == mol.atoms.len();
    (unbranched || centered_star)
        && mol.bonds.len() <= mol.atoms.len()
        && mol.atoms.values().all(|atom| {
            atom.is_carbon()
                || allow_acyl
                    && atom.symbol == "O"
                    && mol.degree(atom.idx) == 1
                    && mol.bond_order(root, atom.idx) == 2
                    && atom.total_h_count == 0
        })
}

/// An excluded virtual parent makes the free valence indistinguishable from
/// the attachment bond already handled by alkyl/alkylidene and acyl naming.
// Provenance: fn substituent_radical | original
// Source: ReShiki | src/naming/openclatura/radical_naming.rs::substituent_radical
fn substituent_radical(mol: &Molecule, root: usize, order: usize) -> Result<String, String> {
    let mut attached = mol.clone();
    for atom in attached.atoms.values_mut() {
        atom.radical_electrons = 0;
    }
    // The cation centre of an ylium name is the same free valence.
    if let Some(atom) = attached.atoms.get_mut(&root) {
        atom.charge = 0;
    }
    let parent = attached.atoms.keys().next_back().copied().unwrap_or(0) + 1;
    let bond = attached.bonds.keys().next_back().copied().unwrap_or(0) + 1;
    let mut atom = Atom::new(parent, "C");
    atom.total_h_count = 4 - order;
    attached.add_atom(atom)?;
    attached.add_bond(Bond::new(bond, root, parent, order))?;
    let name = super::namer::name_subgraph_rendered(
        &attached,
        root,
        &AtomSet::from([parent]),
        Some(parent),
        false,
    )?;
    // P-82.2.1: the name cites isotope labels as the prefix does; the isotope
    // route declines whatever is left uncited.
    let mut uncited = mol.isotope_sites();
    uncited.retain(|site| !name.isotope_atoms.contains(&site.anchor));
    mol.set_isotope_sites(uncited);
    let name = name.text;
    // P-71.2.1.2: 2-methylpropan-2-yl (PIN), not the tert-butyl prefix; the
    // "yl" suffix goes on a parent hydride, which the benzyl prefix is not.
    Ok(match strip_outer_parentheses(&name).as_str() {
        "tert-butyl" => "2-methylpropan-2-yl".into(),
        "benzyl" => "phenylmethyl".into(),
        name => name.into(),
    })
}

/// P-73.2.2.1: a carbon or silicon cation formed by loss of hydride is named
/// like the radical at that position, with "ylium" for "yl" (propan-2-ylium,
/// triphenylsilylium). A cationic centre is senior to every other class
/// present here (P-41), so the rest of the structure is cited as prefixes.
/// Centres with multiple bonds (furan-2-ylium, P-73.2.3.1 acylium) decline.
// Provenance: fn ylium | original
// Source: Nazuchem | crates/nazuchem-s2n/src/radical_naming.rs::ylium
fn ylium(mol: &Molecule) -> Option<Name> {
    let mut charged = mol.atoms.values().filter(|atom| atom.charge != 0);
    let (Some(center), None) = (charged.next(), charged.next()) else {
        return None;
    };
    if center.charge != 1
        || !matches!(center.symbol.as_str(), "C" | "Si")
        || center.is_aromatic
        || mol.degree(center.idx) + center.total_h_count != 3
        || !(mononuclear_center(mol, center.idx)
            || center.is_carbon()
                && mol
                    .get_neighbors(center.idx)
                    .iter()
                    .all(|&n| mol.bond_order(center.idx, n) == 1))
        || mol.atoms.values().any(|atom| atom.has_stereo())
        || mol.bonds.values().any(|bond| bond.stereo.is_some())
        || super::graph_io::get_connected_components(mol).len() != 1
    {
        return None;
    }
    let name = substituent_radical(mol, center.idx, 1).ok()?;
    name.ends_with("yl").then(|| Name {
        text: super::nest_enclosing_marks(&format!("{name}ium")),
        status: mol.name_status(),
        alternatives: vec![],
        omitted_stereo: vec![],
        omitted_units: 0,
        verification: None,
    })
}

// Provenance: fn amine_radical | original
// Source: ReShiki | src/naming/openclatura/radical_naming.rs::amine_radical
fn amine_radical(mol: &Molecule, root: usize, electrons: usize) -> Result<String, String> {
    let atom = &mol.atoms[&root];
    if mol.atoms.len() == 1 {
        return match (atom.total_h_count, electrons) {
            (2, 1) => Ok("azanyl".into()),
            (1, 2) => Ok("azanylidene".into()),
            _ => Err(unsupported()),
        };
    }
    // The pinned parser cannot interpret substituted aminylidene suffixes.
    // Decline these rather than return an unreconstructable preferred form.
    if electrons != 1
        || mol.bonds.len() >= mol.atoms.len()
        || mol
            .atoms
            .values()
            .any(|atom| atom.idx != root && !atom.is_carbon())
        || mol
            .get_neighbors(root)
            .iter()
            .any(|&n| mol.bond_order(root, n) != 1)
        || atom.total_h_count + mol.degree(root) + electrons != 3
    {
        return Err(unsupported());
    }
    let mut capped = mol.clone();
    let center = capped.atoms.get_mut(&root).ok_or_else(unsupported)?;
    center.radical_electrons = 0;
    center.total_h_count += electrons;
    center.explicit_h_count += electrons;
    let component = capped.atoms.keys().copied().collect::<AtomSet>();
    let name = super::component_namer::name_component(&capped, &component)?;
    let stem = name.strip_suffix("amine").ok_or_else(unsupported)?;
    Ok(format!(
        "{stem}{}",
        if electrons == 1 {
            "aminyl"
        } else {
            "aminylidene"
        }
    ))
}

// Provenance: fn chalcogen_radical | original
// Source: ReShiki | src/naming/openclatura/radical_naming.rs::chalcogen_radical
fn chalcogen_radical(mol: &Molecule, root: usize) -> Result<String, String> {
    let atom = &mol.atoms[&root];
    if mol.atoms.len() == 1 {
        return match (atom.symbol.as_str(), atom.total_h_count) {
            ("O", 1) => Ok("hydroxyl".into()),
            ("S", 1) => Ok("sulfanyl".into()),
            ("Se", 1) => Ok("selanyl".into()),
            _ => Err(unsupported()),
        };
    }
    if atom.symbol == "O"
        && mol.atoms.len() == 2
        && atom.total_h_count == 0
        && mol.atoms.values().any(|other| {
            other.idx != root
                && other.symbol == "O"
                && other.total_h_count == 1
                && mol.bond_order(root, other.idx) == 1
        })
    {
        return Ok("hydroperoxyl".into());
    }
    if atom.symbol == "O" && atom.total_h_count == 0 && mol.degree(root) == 1 {
        let nitrogen = mol.get_neighbors(root)[0];
        if mol.atoms[&nitrogen].symbol == "N"
            && mol.bond_order(root, nitrogen) == 1
            && mol
                .get_neighbors(nitrogen)
                .iter()
                .all(|&n| mol.bond_order(nitrogen, n) == 1)
            && mol.degree(nitrogen) + mol.atoms[&nitrogen].total_h_count == 3
            && mol
                .atoms
                .values()
                .all(|other| other.idx == root || other.idx == nitrogen || other.is_carbon())
            && mol.bonds.len() < mol.atoms.len()
        {
            let prefix = super::namer::name_subgraph_with_upstream(
                mol,
                nitrogen,
                &AtomSet::from([root]),
                Some(root),
            )?;
            let prefix = strip_outer_parentheses(&prefix);
            let prefix = super::formatting::phenylamino_spelling(&prefix);
            let stem = prefix.strip_suffix("amino").ok_or_else(unsupported)?;
            return Ok(format!("{stem}aminoxyl"));
        }
    }
    if atom.total_h_count != 0
        || mol.degree(root) != 1
        || mol.bonds.len() > mol.atoms.len()
        || mol
            .get_neighbors(root)
            .iter()
            .any(|&n| mol.bond_order(root, n) != 1)
        || mol
            .atoms
            .values()
            .any(|other| other.idx != root && !other.is_carbon())
    {
        return Err(unsupported());
    }
    let name = substituent_radical(mol, root, 1)?;
    if atom.symbol == "O" {
        // P-71.3.4: methoxy -> methoxyl, including retained ethoxy/phenoxy.
        if !name.ends_with("oxy") {
            return Err(unsupported());
        }
        Ok(format!("{name}l"))
    } else {
        Ok(name)
    }
}

/// The two independent free valences must both receive explicit locants.
/// This path supports unbranched saturated chains, saturated monocycles and
/// benzene; a branch or further characteristic group needs a richer assembly.
// Provenance: fn carbon_diyl | original
// Source: ReShiki | src/naming/openclatura/radical_naming.rs::carbon_diyl
fn carbon_diyl(mol: &Molecule) -> Result<String, String> {
    if mol
        .atoms
        .values()
        .any(|atom| !atom.is_carbon() || atom.radical_electrons > 1)
        || mol.atoms.len() < 2
        || mol.atoms.keys().any(|&i| mol.degree(i) > 2)
    {
        return Err(unsupported());
    }
    let cyclic = mol.bonds.len() == mol.atoms.len();
    let aromatic = mol.atoms.values().all(|atom| atom.is_aromatic);
    if !(mol.bonds.values().all(|bond| bond.order == 1)
        || cyclic && aromatic && mol.atoms.len() == 6)
    {
        return Err(unsupported());
    }
    let mut best = None::<Vec<usize>>;
    for &start in mol.atoms.keys().filter(|&&i| cyclic || mol.degree(i) == 1) {
        for &second in mol.get_neighbors(start) {
            let mut path = vec![start, second];
            while path.len() < mol.atoms.len() {
                let current = *path.last().ok_or_else(unsupported)?;
                let previous = path[path.len() - 2];
                let Some(&next) = mol
                    .get_neighbors(current)
                    .iter()
                    .find(|&&n| n != previous && n != start)
                else {
                    break;
                };
                path.push(next);
            }
            if path.len() != mol.atoms.len() {
                continue;
            }
            let locants: Vec<_> = path
                .iter()
                .enumerate()
                .filter_map(|(index, id)| {
                    (mol.atoms[id].radical_electrons == 1).then_some(index + 1)
                })
                .collect();
            if best.as_ref().is_none_or(|current| locants < *current) {
                best = Some(locants);
            }
        }
    }
    let locants = best.ok_or_else(unsupported)?;
    let parent = if aromatic {
        "benzene".into()
    } else {
        format!(
            "{}{}ane",
            if cyclic { "cyclo" } else { "" },
            super::rules::stems::stem_for(mol.atoms.len())
        )
    };
    // P-71.2.3: benzene-1,4-diyl and ethane-1,2-diyl are the PINs.
    Ok(format!("{parent}-{},{}-diyl", locants[0], locants[1]))
}

#[cfg(test)]
mod tests {
    use super::*;
    use nazuchem_ir::{NameReason, NameStatus};

    // Provenance: fn tests::roundtrip | original
    // Source: ReShiki | src/naming/openclatura/radical_naming.rs::tests::roundtrip
    fn roundtrip(smiles: &str, expected: &str) -> Result<(), String> {
        let mol = super::super::graph_io::read_smiles(smiles, &Default::default())?;
        let name = name_radical(&mol)?.ok_or("Expected radical name")?;
        assert_eq!(name.text, expected, "{smiles}");
        assert_eq!(
            name.status,
            NameStatus::Unverified(NameReason::UnauditedPath),
            "{smiles}"
        );
        let parsed = nazuchem_n2s::Parser::new()
            .map_err(|error| error.to_string())?
            .parse(
                &name.text,
                &nazuchem_n2s::ParseOptions {
                    allow_radicals: true,
                    ..Default::default()
                },
            );
        let structure = parsed
            .structure
            .ok_or_else(|| format!("{}: {}", name.text, parsed.message))?;
        let reconstructed = structure
            .semantic_cxsmiles()
            .map_err(|error| error.to_string())?;
        crate::test_support::verify_identity(smiles, &reconstructed)
    }

    #[test]
    // Provenance: fn tests::carbon_monoradicals_and_single_center_diradicals_roundtrip | adapted
    // Source: ReShiki | src/naming/openclatura/radical_naming.rs::tests::carbon_monoradicals_and_single_center_diradicals_roundtrip
    // Source: openclatura | tests/exploratory/test_parenthesis_functions_direct.py::test_is_substituted_alkyl_ligand
    // Source: openclatura | tests/unit/test_analysis.py::test_rooted_alkyl_ligands_are_derived_as_systematic_substituents
    // Source: openclatura | tests/unit/test_analysis.py::test_is_composite_prefix
    // Source: openclatura | tests/unit/test_analysis.py::test_acetyloxy_substituent_tokens_bind_the_acyl_and_the_ester_oxygen
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn carbon_monoradicals_and_single_center_diradicals_roundtrip() -> Result<(), String> {
        for (smiles, name) in [
            ("[CH3]", "methyl"),
            ("C[CH2]", "ethyl"),
            ("C[CH]C", "propan-2-yl"),
            ("CCC[CH]CCC", "heptan-4-yl"),
            ("C[CH]CCCC", "hexan-2-yl"),
            ("C[C](C)C", "2-methylpropan-2-yl"),
            ("[CH2]C=C", "prop-2-en-1-yl"),
            ("[CH]1CCCCC1", "cyclohexyl"),
            ("[c]1ccccc1", "phenyl"),
            ("[CH2]", "methylidene"),
            ("C[CH]", "ethylidene"),
            ("C[C]C", "propan-2-ylidene"),
            ("[C]1CCCCC1", "cyclohexylidene"),
            ("C[C]=O", "acetyl"),
        ] {
            roundtrip(smiles, name)?;
        }
        Ok(())
    }

    #[test]
    // Provenance: fn tests::separate_radical_centers_roundtrip_with_lowest_locants | original
    // Source: ReShiki | src/naming/openclatura/radical_naming.rs::tests::separate_radical_centers_roundtrip_with_lowest_locants
    fn separate_radical_centers_roundtrip_with_lowest_locants() -> Result<(), String> {
        for (smiles, name) in [
            ("[CH2][CH2]", "ethane-1,2-diyl"),
            ("[CH2]CC[CH2]", "butane-1,4-diyl"),
            ("C[CH]C[CH2]", "butane-1,3-diyl"),
            ("[CH]1C[CH]CCC1", "cyclohexane-1,3-diyl"),
            ("[c]1cc[c]cc1", "benzene-1,4-diyl"),
        ] {
            roundtrip(smiles, name)?;
        }
        Ok(())
    }

    #[test]
    // Provenance: fn tests::heteroatom_radicals_use_p71_forms_and_roundtrip | adapted
    // Source: ReShiki | src/naming/openclatura/radical_naming.rs::tests::heteroatom_radicals_use_p71_forms_and_roundtrip
    // Source: openclatura | tests/data/parser_xml_resources/simpleSubstituents.json::<fixture data>
    // Source: openclatura | tests/unit/test_analysis.py::test_is_composite_prefix
    // Source: openclatura | tests/fuzz/test_ligand_parentheses.py::branched_functional_center_pairs
    // Source: openclatura | tests/unit/test_chalcogen_functional_groups.py::test_neutral_chalcogen_radical_is_not_hydrogen_bearing
    // Source: openclatura | tests/unit/test_chalcogen_functional_groups.py::test_all_peroxy_acyl_site_combinations_have_typed_graph_descriptors
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn heteroatom_radicals_use_p71_forms_and_roundtrip() -> Result<(), String> {
        for (smiles, name) in [
            ("[NH2]", "azanyl"),
            ("[NH]", "azanylidene"),
            ("C[NH]", "methanaminyl"),
            ("C[N]C", "N-methylmethanaminyl"),
            ("[OH]", "hydroxyl"),
            ("O[O]", "hydroperoxyl"),
            ("C[O]", "methoxyl"),
            ("N[O]", "aminoxyl"),
            ("CN(C)[O]", "dimethylaminoxyl"),
            ("CC[O]", "ethoxyl"),
            ("CCC[O]", "propoxyl"),
            ("CCCC[O]", "butoxyl"),
            ("c1ccccc1[O]", "phenoxyl"),
            ("[SH]", "sulfanyl"),
            ("C[S]", "methylsulfanyl"),
            ("[SeH]", "selanyl"),
            ("C[Se]", "methylselanyl"),
        ] {
            roundtrip(smiles, name)?;
        }
        Ok(())
    }

    #[test]
    // Provenance: fn tests::silyl_boranyl_and_substituted_aryl_radicals_roundtrip | original
    // Source: Nazuchem | crates/nazuchem-s2n/src/radical_naming.rs::tests::silyl_boranyl_and_substituted_aryl_radicals_roundtrip
    fn silyl_boranyl_and_substituted_aryl_radicals_roundtrip() -> Result<(), String> {
        // Issue #27. P-71.2.1.1: silyl; P-71.2.1.2: boranyl (not boryl) and
        // "yl" at any position of a parent hydride, other groups as prefixes.
        for (smiles, name) in [
            ("[SiH3]", "silyl"),
            ("C[Si](C)C", "trimethylsilyl"),
            (
                "CCCCO[Si](OCCCC)OCC(CC)CC=C",
                "dibutoxy[(2-ethylpent-4-en-1-yl)oxy]silyl",
            ),
            ("[BH2]", "boranyl"),
            ("C[B]C", "dimethylboranyl"),
            ("[B]1C2C1C=CC2", "6-borabicyclo[3.1.0]hex-2-en-6-yl"),
            ("C1CC1COC2=CC=[C]C=C2", "4-(cyclopropylmethoxy)phenyl"),
            ("Oc1cc[c]cc1", "4-hydroxyphenyl"),
            ("c1cc[c]nc1", "pyridin-2-yl"),
        ] {
            roundtrip(smiles, name)?;
        }
        // Disilanyl, disiloxanyl and silylidene parents are not composed.
        for smiles in ["[SiH3][SiH2]", "C[Si](C)O[Si](C)(C)C", "C[Si]C"] {
            assert!(
                super::super::name_smiles(smiles, &Default::default()).is_err(),
                "{smiles}"
            );
        }
        Ok(())
    }

    #[test]
    // Provenance: fn tests::lambda_convention_radicals_roundtrip | original
    // Source: Nazuchem | crates/nazuchem-s2n/src/radical_naming.rs::tests::lambda_convention_radicals_roundtrip
    fn lambda_convention_radicals_roundtrip() -> Result<(), String> {
        // Issue #27. P-71.2.1.2 with the lambda-convention of P-14.1; the Blue
        // Book example is tert-butoxytri(phenyl)-lambda5-phosphanyl.
        for (smiles, name) in [
            (
                "CC(C)(C)O[P](c1ccccc1)(c1ccccc1)c1ccccc1",
                "(tert-butoxy)tri(phenyl)-lambda^5-phosphanyl",
            ),
            (
                "CC(C)(C)O[PH](OC(C)(C)C)OC(C)(C)C",
                "tris(tert-butoxy)-lambda^5-phosphanyl",
            ),
            ("C[SeH]CCO", "(2-hydroxyethyl)(methyl)-lambda^4-selanyl"),
            ("C[S](C)C", "trimethyl-lambda^4-sulfanyl"),
            ("C[As](C)(C)C", "tetramethyl-lambda^5-arsanyl"),
        ] {
            roundtrip(smiles, name)?;
        }
        // Unsubstituted hydrides, standard bonding numbers and centres more
        // than two above the standard bonding number are not composed.
        for smiles in ["[PH4]", "[SH3]", "C[PH]", "C[S](C)(C)(C)C", "C[PH3]"] {
            assert!(
                super::super::name_smiles(smiles, &Default::default()).is_err(),
                "{smiles}"
            );
        }
        Ok(())
    }

    #[test]
    // Provenance: fn tests::hydride_loss_cations_take_ylium_and_roundtrip | original
    // Source: Nazuchem | crates/nazuchem-s2n/src/radical_naming.rs::tests::hydride_loss_cations_take_ylium_and_roundtrip
    fn hydride_loss_cations_take_ylium_and_roundtrip() -> Result<(), String> {
        // Issue #27. P-73.2.2.1.1: methylium, propylium, cyclobutylium,
        // triphenylsilylium; P-73.2.2.1.2: "ylium" at any position
        // (cyclopenta-2,4-dien-1-ylium), other groups as prefixes.
        for (smiles, name) in [
            ("[CH3+]", "methylium"),
            ("[SiH3+]", "silylium"),
            ("CC[CH2+]", "propylium"),
            ("C[CH+]C", "propan-2-ylium"),
            ("C[C+](C)C", "2-methylpropan-2-ylium"),
            ("[CH+]1CCC1", "cyclobutylium"),
            ("[CH+]1C=CC=C1", "cyclopenta-2,4-dien-1-ylium"),
            ("C1CC2CCC1[CH+]2", "bicyclo[2.2.1]heptan-7-ylium"),
            (
                "CC(C)C[C+](C)CC(C)C(C)C",
                "2,4,6,7-tetramethyloctan-4-ylium",
            ),
            ("[CH2+]c1ccccc1", "phenylmethylium"),
            ("OC(=O)C[CH+]C", "1-carboxypropan-2-ylium"),
            ("c1ccccc1[Si+](c1ccccc1)c1ccccc1", "triphenylsilylium"),
        ] {
            roundtrip(smiles, name)?;
        }
        assert_eq!(
            super::super::name_smiles("C[C+](C)C.[Cl-]", &Default::default())?.text,
            "2-methylpropan-2-ylium chloride"
        );
        // Multiply bonded centres, a second charge and polysilanes still decline.
        for smiles in [
            "C[C+]=O",
            "[c+]1ccccc1",
            "C[N+](C)(C)C[CH+]C",
            "C[Si+](C)[Si](C)(C)C",
        ] {
            assert!(
                super::super::name_smiles(smiles, &Default::default()).is_err(),
                "{smiles}"
            );
        }
        Ok(())
    }

    #[test]
    // Provenance: fn tests::unsupported_radicals_never_reach_closed_shell_naming | adapted
    // Source: ReShiki | src/naming/openclatura/radical_naming.rs::tests::unsupported_radicals_never_reach_closed_shell_naming
    // Source: openclatura | tests/unit/test_analysis.py::test_public_api_golden_names
    // Source: openclatura | tests/exploratory/test_charge_representation_gaps.py::test_bare_hydrogen_less_atoms_are_not_named_as_their_hydride
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn unsupported_radicals_never_reach_closed_shell_naming() -> Result<(), String> {
        for smiles in [
            "[CH2+]",
            "[O]",
            "C[N]",
            "CC(=O)[NH]",
            "ON[O]",
            "[CH2]O",
            "C[CH]C[CH]C(C)C",
            "[CH2]C[CH][CH2]",
            "[PH2]",
            "[CH2]CC(C)C",
            "C[CH]CC(C)C",
            "[CH]1C2CCCCC2CC1",
        ] {
            assert!(
                super::super::name_smiles(smiles, &Default::default()).is_err(),
                "{smiles}"
            );
        }
        let closed = super::super::graph_io::read_smiles("CCC", &Default::default())?;
        assert!(name_radical(&closed)?.is_none());
        Ok(())
    }
}
