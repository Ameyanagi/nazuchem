// SPDX-License-Identifier: MIT AND (MIT OR Apache-2.0)
// Origin: mixed
// Source: openclatura | src/openclatura/special_cases.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | src/openclatura/nitrogen_roles.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | src/openclatura/oxoacid_roles.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | src/openclatura/component_namer.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | src/openclatura/chalcogen_vocabulary.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: ReShiki | src/naming/openclatura/special_cases.rs | 5f559fd482e237307a589dfcb7ea333cee0a6e13
// Source: Nazuchem | crates/nazuchem-s2n/src/special_cases.rs | issue #17
// Source: Nazuchem | crates/nazuchem-s2n/src/special_cases.rs | issue #158
// Source: openclatura | tests/exploratory/test_group14_ether_ligand_parent_selection.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/unit/test_analysis.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/roundtrip/test_substituents.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/roundtrip/test_replacement_parents.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/roundtrip/test_oxyacids.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/roundtrip/test_role_models.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/roundtrip/test_principal_group_breadth.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/unit/test_polyoxide_chains.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/unit/test_chalcogen_functional_groups.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/data/parser_xml_resources/functionalTerms.json | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/unit/test_all_chalcogen_acyl_regressions.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Copyright (c) 2026 lamalab-org
// Copyright (c) 2026 Ameyanagi and ReShiki contributors
// Changes: Large openclatura shortcut/ligand/functional-class adaptation, consolidated across role modules, with local shortcut safety/certification extensions and tests that include adapted upstream cases.
// Local contributions: MIT OR Apache-2.0; see NOTICE and PROVENANCE-nazuchem.toml.

//! Whole-component shortcuts from Openclatura's `special_cases.py`.
//!
//! Matching is graph based; bindings and audit-only metadata are omitted from
//! this plain-name port.

use std::collections::{BTreeMap, BTreeSet, VecDeque};

use super::assembly_prefixes::substituent_sort_key;
use super::chains::{canonical_ranks, get_cyclic_atoms};
use super::formatting::{
    descriptors_omitted, format_center_ligands_by, format_multiplier, is_complex_prefix,
    needs_complex_multiplier, needs_simple_multiplier_parentheses, strip_outer_parentheses,
};
use super::graph_kernel::{Edge, biconnected_edge_components, normalize_edge};
use super::molecule::{Atom, AtomSet, Bond, Molecule, component_atoms_until_blocked};
use super::namer::{
    name_subgraph_rendered, name_subgraph_with_upstream, name_subgraph_without_attachment_stereo,
};
use super::nomenclature::RULES;
use super::perception::{PerceivedGroup, perceive_groups};
use super::rules::{elements, multipliers, retained, stems};
use super::subgraph_tools::subgraph_component;
use nazuchem_ir::{NameReason, NameStatus};

const SHORTCUT_STEREO_ERROR: &str =
    "openclatura: shortcut cannot preserve configured stereochemistry";

// Provenance: fn set | adapted
// Source: openclatura | src/openclatura/special_cases.py::peroxide_linkage_result
// Source: openclatura | src/openclatura/special_cases.py::organophosphinic_acid_result
fn set<const N: usize>(atoms: [usize; N]) -> AtomSet {
    atoms.into_iter().collect()
}
// Provenance: fn all_atoms | adapted
// Source: openclatura | src/openclatura/special_cases.py::simple_azine_parent_name
// Source: openclatura | src/openclatura/special_cases.py::phosphane_borane_zwitterion_result
fn all_atoms(mol: &Molecule) -> AtomSet {
    mol.atoms.keys().copied().collect()
}
// Provenance: fn outside | adapted
// Source: openclatura | src/openclatura/special_cases.py::peroxide_linkage_result
// Source: openclatura | src/openclatura/special_cases.py::simple_azine_parent_name
// Source: openclatura | src/openclatura/special_cases.py::phosphane_borane_zwitterion_result
fn outside(mol: &Molecule, atoms: &AtomSet) -> AtomSet {
    all_atoms(mol).difference(atoms).copied().collect()
}
// Provenance: fn require_stereo_free_core | adapted
// Source: openclatura | src/openclatura/special_cases.py::phosphane_borane_zwitterion_result
// Source: openclatura | src/openclatura/special_cases.py::sulfonium_ylide_name
// Source: openclatura | src/openclatura/special_cases.py::oxoacid_ester_result
fn require_stereo_free_core(mol: &Molecule, atoms: &AtomSet) -> Result<(), String> {
    // Upstream divergence: manually spelled shortcut cores must not discard configured atom or bond stereo.
    if atoms.iter().any(|&idx| mol.atoms[&idx].has_stereo())
        || mol
            .bonds
            .values()
            .any(|bond| atoms.contains(&bond.u) && atoms.contains(&bond.v) && bond.stereo.is_some())
    {
        return Err(SHORTCUT_STEREO_ERROR.into());
    }
    Ok(())
}
// Provenance: fn center_ligands | adapted
// Source: openclatura | src/openclatura/special_cases.py::_grouped_ligand_prefix
fn center_ligands(names: &[String]) -> String {
    format_center_ligands_by(names, substituent_sort_key)
}
// Provenance: fn branch_name | adapted
// Source: openclatura | src/openclatura/special_cases.py::phosphane_borane_zwitterion_result
// Source: openclatura | src/openclatura/special_cases.py::sulfonium_ylide_name
fn branch_name(
    mol: &Molecule,
    component: &AtomSet,
    root: usize,
    blocked: &AtomSet,
    upstream: usize,
) -> Result<String, String> {
    let mut exclude = outside(mol, component);
    exclude.extend(blocked);
    name_subgraph_with_upstream(mol, root, &exclude, Some(upstream))
}

// Provenance: fn single_atom_component_name | adapted
// Source: openclatura | src/openclatura/special_cases.py::single_atom_component_name
pub(super) fn single_atom_component_name(
    mol: &Molecule,
    component_atoms: &AtomSet,
) -> Option<String> {
    if component_atoms.len() != 1 {
        return None;
    }
    let atom = &mol.atoms[component_atoms.first()?];
    // Upstream divergence: a retained single-atom hydride does not cite atom stereo.
    if mol.atoms[&atom.idx].has_stereo() {
        return None;
    }
    if atom.charge > 0 && RULES.ions.single_atom_cations.contains(&atom.symbol) {
        return elements::get(&atom.symbol).map(|e| e.name.clone());
    }
    if atom.charge < 0
        && let Some(name) = RULES.ions.single_atom_anions.get(&atom.symbol)
    {
        return Some(name.clone());
    }
    let ion = match (atom.symbol.as_str(), atom.charge, atom.total_h_count) {
        ("N", 1, 4) => "ammonium",
        ("N", -1, 2) => "azanide",
        ("O", -1, 1) => "hydroxide",
        ("O", -2, 0) => "oxide",
        ("O", 1, 3) => "oxonium",
        ("S", -1, 1) => "sulfanide",
        ("S", -2, 0) => "sulfide",
        ("P", 1, 4) => "phosphanium",
        ("H", 1, 0) => "hydron",
        ("H", -1, 0) => "hydride",
        ("Se", -2, 0) => "selenide",
        _ => "",
    };
    if !ion.is_empty() {
        return Some(ion.into());
    }
    if atom.charge != 0 || atom.total_h_count != elements::get(&atom.symbol)?.standard_valence {
        return None;
    }
    let retained = match atom.symbol.as_str() {
        "N" => "ammonia",
        "O" => "water",
        "F" => "hydrogen fluoride",
        "Cl" => "hydrogen chloride",
        "Br" => "hydrogen bromide",
        "I" => "hydrogen iodide",
        _ => "",
    };
    if !retained.is_empty() {
        return Some(retained.into());
    }
    RULES
        .components
        .mononuclear_parent_hydrides
        .get(&atom.symbol)
        .cloned()
}

#[derive(PartialEq, Eq)]
struct SmallSignature {
    atoms: Vec<(String, i32)>,
    bonds: Vec<(String, String, usize)>,
}
type SmallParentSpec<'a> = (
    &'a str,
    &'a [(&'a str, i32)],
    &'a [(&'a str, &'a str, usize)],
);
// Provenance: fn small_signature | adapted
// Source: openclatura | src/openclatura/special_cases.py::_small_inorganic_signature
fn small_signature(mol: &Molecule, atoms: &AtomSet) -> SmallSignature {
    let mut atom_sig: Vec<_> = atoms
        .iter()
        .map(|i| (mol.atoms[i].symbol.clone(), mol.atoms[i].charge))
        .collect();
    atom_sig.sort();
    let mut bonds: Vec<_> = mol
        .bonds
        .values()
        .filter(|b| atoms.contains(&b.u) && atoms.contains(&b.v))
        .map(|b| {
            let (a, z) = (&mol.atoms[&b.u].symbol, &mol.atoms[&b.v].symbol);
            (a.min(z).clone(), a.max(z).clone(), b.order)
        })
        .collect();
    bonds.sort();
    SmallSignature {
        atoms: atom_sig,
        bonds,
    }
}

// Provenance: fn small_inorganic_parent_name | adapted
// Source: openclatura | src/openclatura/special_cases.py::small_inorganic_parent_result
fn small_inorganic_parent_name(mol: &Molecule, atoms: &AtomSet) -> Option<String> {
    if !(2..=9).contains(&atoms.len())
        || atoms
            .iter()
            .any(|i| mol.atoms[i].is_carbon() && mol.atoms[i].total_h_count != 0)
    {
        return None;
    }
    // The upstream retained table is a multiset of charged atoms and bonds.
    let specs: &[SmallParentSpec<'_>] = &[
        ("hydroxylamine", &[("N", 0), ("O", 0)], &[("N", "O", 1)]),
        ("hydrogen peroxide", &[("O", 0), ("O", 0)], &[("O", "O", 1)]),
        (
            "phosphonic acid",
            &[("O", 0), ("O", 0), ("O", 0), ("P", 0)],
            &[("O", "P", 1), ("O", "P", 1), ("O", "P", 2)],
        ),
        (
            "phosphinic acid",
            &[("O", 0), ("O", 0), ("P", 0)],
            &[("O", "P", 1), ("O", "P", 2)],
        ),
        (
            "sulfuric diamide",
            &[("N", 0), ("N", 0), ("O", 0), ("O", 0), ("S", 0)],
            &[("N", "S", 1), ("N", "S", 1), ("O", "S", 2), ("O", "S", 2)],
        ),
        (
            "diphosphoric acid",
            &[
                ("O", 0),
                ("O", 0),
                ("O", 0),
                ("O", 0),
                ("O", 0),
                ("O", 0),
                ("O", 0),
                ("P", 0),
                ("P", 0),
            ],
            &[
                ("O", "P", 1),
                ("O", "P", 1),
                ("O", "P", 1),
                ("O", "P", 1),
                ("O", "P", 1),
                ("O", "P", 1),
                ("O", "P", 2),
                ("O", "P", 2),
            ],
        ),
        (
            "carbon dioxide",
            &[("C", 0), ("O", 0), ("O", 0)],
            &[("C", "O", 2), ("C", "O", 2)],
        ),
        (
            "carbon disulfide",
            &[("C", 0), ("S", 0), ("S", 0)],
            &[("C", "S", 2), ("C", "S", 2)],
        ),
        (
            "carbonyl sulfide",
            &[("C", 0), ("O", 0), ("S", 0)],
            &[("C", "O", 2), ("C", "S", 2)],
        ),
        ("carbon monoxide", &[("C", -1), ("O", 1)], &[("C", "O", 3)]),
        (
            "carbonyl dichloride",
            &[("C", 0), ("O", 0), ("Cl", 0), ("Cl", 0)],
            &[("C", "O", 2), ("C", "Cl", 1), ("C", "Cl", 1)],
        ),
        (
            "carbonyl difluoride",
            &[("C", 0), ("O", 0), ("F", 0), ("F", 0)],
            &[("C", "O", 2), ("C", "F", 1), ("C", "F", 1)],
        ),
        (
            "carbonyl dibromide",
            &[("C", 0), ("O", 0), ("Br", 0), ("Br", 0)],
            &[("C", "O", 2), ("C", "Br", 1), ("C", "Br", 1)],
        ),
        (
            "thiocyanic acid",
            &[("C", 0), ("N", 0), ("S", 0)],
            &[("C", "N", 3), ("C", "S", 1)],
        ),
        (
            "isocyanic acid",
            &[("C", 0), ("N", 0), ("O", 0)],
            &[("C", "N", 2), ("C", "O", 2)],
        ),
        (
            "isothiocyanic acid",
            &[("C", 0), ("N", 0), ("S", 0)],
            &[("C", "N", 2), ("C", "S", 2)],
        ),
        (
            "cyanate",
            &[("C", 0), ("N", 0), ("O", -1)],
            &[("C", "N", 3), ("C", "O", 1)],
        ),
        (
            "cyanate",
            &[("C", 0), ("N", -1), ("O", 0)],
            &[("C", "N", 2), ("C", "O", 2)],
        ),
        (
            "thiocyanate",
            &[("C", 0), ("N", 0), ("S", -1)],
            &[("C", "N", 3), ("C", "S", 1)],
        ),
        (
            "thiocyanate",
            &[("C", 0), ("N", -1), ("S", 0)],
            &[("C", "N", 2), ("C", "S", 2)],
        ),
        ("cyanide", &[("C", -1), ("N", 0)], &[("C", "N", 3)]),
        (
            "azide",
            &[("N", -1), ("N", 1), ("N", -1)],
            &[("N", "N", 2), ("N", "N", 2)],
        ),
        (
            "nitrate",
            &[("N", 1), ("O", 0), ("O", -1), ("O", -1)],
            &[("N", "O", 2), ("N", "O", 1), ("N", "O", 1)],
        ),
        (
            "nitrite",
            &[("N", 0), ("O", 0), ("O", -1)],
            &[("N", "O", 2), ("N", "O", 1)],
        ),
        ("peroxide", &[("O", -1), ("O", -1)], &[("O", "O", 1)]),
    ];
    let signature = small_signature(mol, atoms);
    for &(name, atom_spec, bond_spec) in specs {
        let mut a: Vec<_> = atom_spec.iter().map(|&(s, c)| (s.to_owned(), c)).collect();
        a.sort();
        let mut b: Vec<_> = bond_spec
            .iter()
            .map(|&(s, t, o)| (s.min(t).to_owned(), s.max(t).to_owned(), o))
            .collect();
        b.sort();
        if signature == (SmallSignature { atoms: a, bonds: b }) {
            return Some(name.to_owned());
        }
    }
    None
}

// Provenance: fn carbon_ligand_atoms | adapted
// Source: openclatura | src/openclatura/special_cases.py::_carbon_ligand_atoms
fn carbon_ligand_atoms(mol: &Molecule, atoms: &AtomSet, root: usize, central: usize) -> AtomSet {
    let mut found = AtomSet::new();
    let mut queue = VecDeque::from([root]);
    while let Some(idx) = queue.pop_front() {
        if found.contains(&idx) {
            continue;
        }
        if !atoms.contains(&idx) || !mol.atoms[&idx].is_carbon() {
            return AtomSet::new();
        }
        found.insert(idx);
        queue.extend(
            mol.get_neighbors(idx)
                .iter()
                .copied()
                .filter(|n| *n != central && atoms.contains(n)),
        );
    }
    found
}

// Provenance: fn carbon_path | adapted
// Source: openclatura | src/openclatura/special_cases.py::_carbon_path_between
fn carbon_path(
    mol: &Molecule,
    atoms: &AtomSet,
    start: usize,
    end: usize,
) -> Result<Vec<usize>, String> {
    let mut queue = VecDeque::from([(start, vec![start])]);
    let mut visited = AtomSet::new();
    while let Some((idx, path)) = queue.pop_front() {
        // Upstream divergence: shortcut path searches share the request's exhaustive-work budget.
        mol.spend_naming_work(path.len())?;
        if idx == end {
            return Ok(path);
        }
        if !visited.insert(idx) {
            continue;
        }
        let mut neighbors: Vec<_> = mol
            .get_neighbors(idx)
            .iter()
            .copied()
            .filter(|n| atoms.contains(n) && !path.contains(n))
            .collect();
        neighbors.sort();
        for n in neighbors {
            let mut p = path.clone();
            p.push(n);
            queue.push_back((n, p));
        }
    }
    Ok(Vec::new())
}

// Provenance: fn alkyl_ligand_name | adapted
// Source: openclatura | src/openclatura/special_cases.py::_alkyl_ligand_name
fn alkyl_ligand_name(
    mol: &Molecule,
    atoms: &AtomSet,
    root: usize,
    central: usize,
) -> Result<String, String> {
    let carbons = carbon_ligand_atoms(mol, atoms, root, central);
    if carbons.is_empty() {
        return Ok(String::new());
    }
    let mut degree_sum = 0;
    for &idx in &carbons {
        for &n in mol.get_neighbors(idx) {
            if !atoms.contains(&n) || n == central {
                continue;
            }
            if mol.bond_order(idx, n) != 1 || !mol.atoms[&n].is_carbon() {
                return Ok(String::new());
            }
            degree_sum += 1;
        }
    }
    if degree_sum / 2 != carbons.len() - 1 {
        return Ok(String::new());
    }
    // Upstream divergence: the direct alkyl spelling cannot represent ligand stereo or formal charge.
    require_stereo_free_core(mol, &carbons)?;
    if carbons.iter().any(|idx| mol.atoms[idx].charge != 0) {
        return Err("openclatura: shortcut cannot preserve charged alkyl ligand".into());
    }
    systematic_alkyl_name(mol, &carbons, root)
}

// Provenance: fn systematic_alkyl_name | adapted
// Source: openclatura | src/openclatura/special_cases.py::_systematic_alkyl_ligand_name
fn systematic_alkyl_name(mol: &Molecule, atoms: &AtomSet, root: usize) -> Result<String, String> {
    let _depth = mol.enter_naming()?;
    if atoms.len() == 1 {
        return Ok("methyl".into());
    }
    let atom_list: Vec<_> = atoms.iter().copied().collect();
    let mut choices = Vec::new();
    for (i, &left) in atom_list.iter().enumerate() {
        for &right in &atom_list[i..] {
            let path = carbon_path(mol, atoms, left, right)?;
            if !path.contains(&root) {
                continue;
            }
            for p in [path.clone(), path.into_iter().rev().collect()] {
                mol.spend_naming_work(p.len())?;
                let Some(root_offset) = p.iter().position(|i| *i == root) else {
                    continue;
                };
                let root_loc = root_offset + 1;
                let locs: Vec<_> = p
                    .iter()
                    .enumerate()
                    .filter(|(_, a)| {
                        mol.get_neighbors(**a)
                            .iter()
                            .any(|n| atoms.contains(n) && !p.contains(n))
                    })
                    .map(|(i, _)| i + 1)
                    .collect();
                choices.push((std::cmp::Reverse(p.len()), root_loc, locs, p));
            }
        }
    }
    let Some((_, root_loc, _, path)) = choices.into_iter().min() else {
        return Ok(String::new());
    };
    let path_set: AtomSet = path.iter().copied().collect();
    let mut locants: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (loc, &idx) in path.iter().enumerate() {
        for &n in mol.get_neighbors(idx) {
            if !atoms.contains(&n) || path_set.contains(&n) {
                continue;
            }
            let branch = carbon_ligand_atoms(mol, atoms, n, idx);
            let name = systematic_alkyl_name(mol, &branch, n)?;
            if name.is_empty() {
                return Ok(String::new());
            }
            locants.entry(name).or_default().push(loc + 1);
        }
    }
    let mut prefix = String::new();
    for (name, locs) in locants {
        prefix.push_str(&format!(
            "{}-{}",
            locs.iter()
                .map(usize::to_string)
                .collect::<Vec<_>>()
                .join(","),
            format_multiplier(&name, locs.len(), false)
        ));
    }
    let stem = stems::stem_for(path.len());
    let parent = if root_loc == 1 {
        format!("{stem}yl")
    } else {
        format!("{stem}an-{root_loc}-yl")
    };
    Ok(format!("{prefix}{parent}"))
}

// Provenance: fn ligand_name | adapted
// Source: openclatura | src/openclatura/special_cases.py::_phosphorus_ligand_name
fn ligand_name(
    mol: &Molecule,
    atoms: &AtomSet,
    root: usize,
    central: usize,
    blocked: &AtomSet,
) -> Result<String, String> {
    let alkyl = alkyl_ligand_name(mol, atoms, root, central)?;
    if !alkyl.is_empty() {
        return Ok(alkyl);
    }
    Ok(strip_outer_parentheses(&branch_name(
        mol, atoms, root, blocked, central,
    )?))
}

#[derive(Clone, Debug)]
struct AzineSide {
    carbon: usize,
    nitrogen: usize,
    atoms: AtomSet,
}

/// The graph roles from `nitrogen_roles.azine_roles`: the two carbon sides
/// must partition the component after the C=N-N=C nitrogens are blocked.
// Provenance: fn azine_roles | adapted
// Source: openclatura | src/openclatura/nitrogen_roles.py::azine_roles
fn azine_roles(mol: &Molecule, atoms: &AtomSet) -> Vec<[AzineSide; 2]> {
    let mut roles = Vec::new();
    for &n1 in atoms {
        if mol.atoms[&n1].symbol != "N" {
            continue;
        }
        for &n2 in mol.get_neighbors(n1) {
            if n2 <= n1
                || !atoms.contains(&n2)
                || mol.atoms[&n2].symbol != "N"
                || mol.bond_order(n1, n2) != 1
            {
                continue;
            }
            let double_carbon = |nitrogen, blocked| {
                let candidates: Vec<_> = mol
                    .get_neighbors(nitrogen)
                    .iter()
                    .copied()
                    .filter(|n| {
                        *n != blocked
                            && mol.atoms[n].is_carbon()
                            && mol.bond_order(nitrogen, *n) == 2
                    })
                    .collect();
                (candidates.len() == 1).then(|| candidates[0])
            };
            let (Some(c1), Some(c2)) = (double_carbon(n1, n2), double_carbon(n2, n1)) else {
                continue;
            };
            let core = set([n1, n2]);
            let side1 = component_atoms_until_blocked(mol, atoms, c1, &core);
            let side2 = component_atoms_until_blocked(mol, atoms, c2, &core);
            if side1.is_empty() || side2.is_empty() || !side1.is_disjoint(&side2) {
                continue;
            }
            if side1
                .iter()
                .chain(&side2)
                .chain(&core)
                .copied()
                .collect::<AtomSet>()
                != *atoms
            {
                continue;
            }
            roles.push([
                AzineSide {
                    carbon: c1,
                    nitrogen: n1,
                    atoms: side1,
                },
                AzineSide {
                    carbon: c2,
                    nitrogen: n2,
                    atoms: side2,
                },
            ]);
        }
    }
    roles
}

// Provenance: fn is_amidino_carbon | adapted
// Source: openclatura | src/openclatura/special_cases.py::_is_amidino_carbon
fn is_amidino_carbon(mol: &Molecule, carbon: usize, nitrogen: usize) -> bool {
    let others: Vec<_> = mol
        .get_neighbors(carbon)
        .iter()
        .copied()
        .filter(|n| *n != nitrogen)
        .collect();
    others.len() == 2
        && others
            .iter()
            .all(|n| mol.atoms[n].symbol == "N" && mol.bond_order(carbon, *n) == 1)
}

// Provenance: fn hydrazone_name | adapted
// Source: openclatura | src/openclatura/special_cases.py::inline in simple_azine_parent_name
fn hydrazone_name(
    mol: &Molecule,
    carbon: usize,
    nitrogen: usize,
    parent: &str,
    ylidene: &str,
) -> String {
    let name = format!("{parent} {ylidene}hydrazone");
    let Some(bond) = mol.get_bond(carbon, nitrogen) else {
        return name;
    };
    let Some(stereo) = bond
        .stereo
        .as_deref()
        .filter(|stereo| matches!(*stereo, "E" | "Z" | "Ξ"))
    else {
        return name;
    };
    // Functional-class words obey the same positional-equivalence proof.
    let locant = parent
        .strip_suffix("-one")
        .and_then(|s| s.rsplit_once("an-"))
        .map(|(_, locant)| locant)
        .filter(|locant| !locant.is_empty() && locant.bytes().all(|b| b.is_ascii_digit()))
        .unwrap_or(
            if name.contains(|c: char| c.is_ascii_digit())
                || !super::locant_elision::bond_stereo_locant_is_unambiguous(mol, bond.idx)
            {
                "1"
            } else {
                ""
            },
        );
    if descriptors_omitted() {
        return name;
    }
    format!("({locant}{stereo})-{name}")
}

/// Port of the ordered `simple_azine_parent_name` shortcut. In particular,
/// substituted amidino sides use their recursively named ylidenes before an
/// anhydride elsewhere in a ligand can choose the whole-component citation.
// Provenance: fn simple_azine_parent_name | adapted
// Source: openclatura | src/openclatura/special_cases.py::simple_azine_parent_name
fn simple_azine_parent_name(mol: &Molecule, atoms: &AtomSet) -> Result<Option<String>, String> {
    for sides in azine_roles(mol, atoms) {
        let status = mol.name_status();
        require_stereo_free_core(mol, &set([sides[0].nitrogen, sides[1].nitrogen]))?;
        let parents = [
            carbonyl_equivalent_side_name(mol, &sides[0], false)?,
            carbonyl_equivalent_side_name(mol, &sides[1], false)?,
        ];
        // Only hydrazone names cite a parent, and those are marked as a whole.
        mol.set_name_status(status);
        let ylidenes = [
            carbonyl_equivalent_side_name(mol, &sides[0], true)?,
            carbonyl_equivalent_side_name(mol, &sides[1], true)?,
        ];
        if parents.iter().chain(&ylidenes).any(String::is_empty) {
            for (parent_side, parent) in parents.iter().enumerate() {
                if parent.is_empty() {
                    continue;
                }
                let other = &sides[1 - parent_side];
                let ylidene = strip_outer_parentheses(&branch_name(
                    mol,
                    atoms,
                    other.carbon,
                    &set([other.nitrogen]),
                    other.nitrogen,
                )?);
                if !ylidene.is_empty() && ylidene.ends_with("ylidene") {
                    let side = &sides[parent_side];
                    let name = hydrazone_name(mol, side.carbon, side.nitrogen, parent, &ylidene);
                    mol.mark_name(NameStatus::NotPreferred(NameReason::FunctionalClass));
                    return Ok(Some(name));
                }
            }
            let amidino = sides
                .each_ref()
                .map(|side| is_amidino_carbon(mol, side.carbon, side.nitrogen));
            let aldehyde = sides.each_ref().map(|side| {
                mol.get_neighbors(side.carbon)
                    .iter()
                    .filter(|n| **n != side.nitrogen && mol.atoms[n].symbol != "H")
                    .count()
                    <= 1
            });
            let substituted_amidino = sides.iter().zip(amidino).any(|(side, amidino)| {
                amidino
                    && mol.get_neighbors(side.carbon).iter().any(|ligand| {
                        mol.atoms[ligand].symbol == "N"
                            && *ligand != side.nitrogen
                            && mol.get_neighbors(*ligand).iter().any(|other| {
                                *other != side.carbon && mol.atoms[other].symbol != "H"
                            })
                    })
            });
            if amidino.into_iter().filter(|yes| *yes).count() == 1
                && (substituted_amidino
                    || !aldehyde
                        .into_iter()
                        .zip(amidino)
                        .any(|(aldehyde, amidino)| aldehyde && !amidino))
            {
                let mut names = Vec::new();
                for side in &sides {
                    names.push(strip_outer_parentheses(&branch_name(
                        mol,
                        atoms,
                        side.carbon,
                        &set([side.nitrogen]),
                        side.nitrogen,
                    )?));
                }
                if names
                    .iter()
                    .all(|name| !name.is_empty() && name.ends_with("ylidene"))
                {
                    names.sort_by_key(|name| substituent_sort_key(name));
                    return Ok(Some(format!("1-({})-2-({})hydrazine", names[0], names[1])));
                }
            }
            // Side names of a skipped role are not cited.
            mol.set_name_status(status);
            continue;
        }
        let ranks = canonical_ranks(mol, None);
        let left_bond = mol.get_bond(sides[0].carbon, sides[0].nitrogen);
        let right_bond = mol.get_bond(sides[1].carbon, sides[1].nitrogen);
        if ranks[&sides[0].carbon] == ranks[&sides[1].carbon]
            && ranks[&sides[0].nitrogen] == ranks[&sides[1].nitrogen]
            && let (Some(left), Some(right)) = (left_bond, right_bond)
            && left.stereo == right.stereo
        {
            let mut exclude = outside(mol, atoms);
            exclude.insert(sides[0].nitrogen);
            let stereo = left.stereo.as_deref().filter(|s| matches!(*s, "E" | "Z"));
            // The hydrazine attachment locants cite the outward C=N stereo.
            let ylidene = if stereo.is_some() {
                name_subgraph_without_attachment_stereo(
                    mol,
                    sides[0].carbon,
                    &exclude,
                    sides[0].nitrogen,
                )?
            } else {
                name_subgraph_with_upstream(
                    mol,
                    sides[0].carbon,
                    &exclude,
                    Some(sides[0].nitrogen),
                )?
            };
            let ylidene = strip_outer_parentheses(&ylidene);
            mol.set_name_status(status);
            let complex = needs_complex_multiplier(&ylidene);
            let multiplier = if complex { "bis" } else { "di" };
            let prefix = if complex
                || is_complex_prefix(&ylidene)
                || needs_simple_multiplier_parentheses(&ylidene)
            {
                format!("{multiplier}({ylidene})")
            } else {
                format!("{multiplier}{ylidene}")
            };
            return Ok(Some(
                if let Some(stereo) = stereo.filter(|_| !descriptors_omitted()) {
                    format!("(1{stereo},2{stereo})-{prefix}hydrazine")
                } else {
                    format!("{prefix}hydrazine")
                },
            ));
        }
        let parent_side =
            usize::from((sides[1].atoms.len(), &parents[1]) > (sides[0].atoms.len(), &parents[0]));
        let side = &sides[parent_side];
        // Upstream divergence: an asymmetric hydrazone citation does not render the other C=N descriptor.
        if mol
            .get_bond(
                sides[1 - parent_side].carbon,
                sides[1 - parent_side].nitrogen,
            )
            .is_some_and(|bond| bond.stereo.is_some())
        {
            return Err(SHORTCUT_STEREO_ERROR.into());
        }
        let name = hydrazone_name(
            mol,
            side.carbon,
            side.nitrogen,
            &parents[parent_side],
            &ylidenes[1 - parent_side],
        );
        mol.mark_name(NameStatus::NotPreferred(NameReason::FunctionalClass));
        return Ok(Some(name));
    }
    Ok(None)
}

// Provenance: fn ordered_simple_ring | adapted
// Source: openclatura | src/openclatura/special_cases.py::_ordered_simple_ring
fn ordered_simple_ring(mol: &Molecule, atoms: &AtomSet) -> Vec<usize> {
    if atoms.len() < 3
        || atoms.iter().any(|i| {
            mol.get_neighbors(*i)
                .iter()
                .filter(|n| atoms.contains(n))
                .count()
                != 2
        })
    {
        return Vec::new();
    }
    let Some(&start) = atoms.first() else {
        return Vec::new();
    };
    let mut path = vec![start];
    let mut previous = None;
    loop {
        let current = *path.last().unwrap_or(&start);
        let next = mol
            .get_neighbors(current)
            .iter()
            .copied()
            .filter(|n| atoms.contains(n) && Some(*n) != previous)
            .min();
        let Some(next) = next else {
            return Vec::new();
        };
        if next == start {
            return if path.len() == atoms.len() {
                path
            } else {
                Vec::new()
            };
        }
        if path.contains(&next) {
            return Vec::new();
        }
        previous = Some(current);
        path.push(next);
    }
}

// Provenance: fn simple_cycle_unsaturation_locants | adapted
// Source: openclatura | src/openclatura/special_cases.py::_simple_cycle_unsaturation_locants
fn simple_cycle_unsaturation_locants(mol: &Molecule, path: &[usize]) -> Vec<usize> {
    path.iter()
        .enumerate()
        .filter_map(|(i, atom)| {
            (mol.bond_order(*atom, path[(i + 1) % path.len()]) == 2).then_some(i + 1)
        })
        .collect()
}

// Provenance: fn simple_replacement_cycle_name | adapted
// Source: openclatura | src/openclatura/special_cases.py::_simple_replacement_cycle_name
fn simple_replacement_cycle_name(mol: &Molecule, path: &[usize]) -> String {
    let mut hetero = BTreeMap::<(usize, String), Vec<usize>>::new();
    for (offset, atom) in path.iter().enumerate() {
        if mol.atoms[atom].is_carbon() {
            continue;
        }
        let Some(element) = elements::get(&mol.atoms[atom].symbol) else {
            return String::new();
        };
        let Some(prefix) = element.hw_stem.as_ref().filter(|s| !s.is_empty()) else {
            return String::new();
        };
        hetero
            .entry((element.hw_priority.unwrap_or(999), prefix.clone()))
            .or_default()
            .push(offset + 1);
    }
    let mut prefix = String::new();
    for ((_, name), locants) in hetero {
        let multiplier = if locants.len() == 1 {
            String::new()
        } else {
            multipliers::basic(locants.len())
        };
        prefix.push_str(&format!(
            "{}-{multiplier}{name}",
            locants
                .iter()
                .map(usize::to_string)
                .collect::<Vec<_>>()
                .join(",")
        ));
    }
    let stem = stems::stem_for(path.len());
    let locants = simple_cycle_unsaturation_locants(mol, path);
    if locants.is_empty() {
        return format!("{prefix}cyclo{stem}an");
    }
    let multiplier = if locants.len() == 1 {
        String::new()
    } else {
        multipliers::basic(locants.len())
    };
    format!(
        "{prefix}cyclo{stem}a-{}-{multiplier}en",
        locants
            .iter()
            .map(usize::to_string)
            .collect::<Vec<_>>()
            .join(",")
    )
}

// Provenance: fn carbonyl_side_locant_sort_key | adapted
// Source: openclatura | src/openclatura/special_cases.py::_retained_locant_sort_key
fn carbonyl_side_locant_sort_key(locant: &str) -> (usize, &str) {
    let digits: String = locant.chars().take_while(char::is_ascii_digit).collect();
    (digits.parse().unwrap_or(999), locant)
}

// Provenance: fn carbaldehyde_side_ring | adapted
// Source: openclatura | src/openclatura/special_cases.py::_carbaldehyde_side_ring
fn carbaldehyde_side_ring(mol: &Molecule, side: &AzineSide) -> Option<(AtomSet, usize)> {
    if !side.atoms.contains(&side.carbon) || !mol.atoms[&side.carbon].is_carbon() {
        return None;
    }
    let ring: AtomSet = side
        .atoms
        .difference(&set([side.carbon]))
        .copied()
        .collect();
    if ring.len() < 3 {
        return None;
    }
    let neighbors: Vec<_> = mol
        .get_neighbors(side.carbon)
        .iter()
        .copied()
        .filter(|n| ring.contains(n) && mol.bond_order(side.carbon, *n) == 1)
        .collect();
    (neighbors.len() == 1).then(|| (ring, neighbors[0]))
}

// Provenance: fn carbonyl_equivalent_side_name | adapted
// Source: openclatura | src/openclatura/special_cases.py::_carbonyl_equivalent_side_name
fn carbonyl_equivalent_side_name(
    mol: &Molecule,
    side: &AzineSide,
    as_ylidene: bool,
) -> Result<String, String> {
    let name = carbonyl_equivalent_side_name_unchecked(mol, side, as_ylidene)?;
    if !name.is_empty() {
        require_stereo_free_core(mol, &side.atoms)?;
        // Upstream divergence: carbonyl-equivalent side citations do not render formal charges.
        if side.atoms.iter().any(|idx| mol.atoms[idx].charge != 0) {
            return Err("openclatura: shortcut cannot preserve charged carbonyl side".into());
        }
    }
    Ok(name)
}

// Provenance: fn carbonyl_equivalent_side_name_unchecked | adapted
// Source: openclatura | src/openclatura/special_cases.py::_carbonyl_equivalent_side_name
fn carbonyl_equivalent_side_name_unchecked(
    mol: &Molecule,
    side: &AzineSide,
    as_ylidene: bool,
) -> Result<String, String> {
    if let Some((ring, attachment)) = carbaldehyde_side_ring(mol, side) {
        let path = ordered_simple_ring(mol, &ring);
        if !path.is_empty() {
            let homocyclic = path
                .iter()
                .all(|i| mol.atoms[i].symbol == mol.atoms[&path[0]].symbol);
            if let Some(parent) = retained::get_retained_ring(mol, &path) {
                let mut locant = parent
                    .locant_maps
                    .iter()
                    .filter_map(|map| map.get(&attachment))
                    .min_by_key(|s| carbonyl_side_locant_sort_key(s))
                    .cloned()
                    .unwrap_or_else(|| {
                        if homocyclic {
                            "1".into()
                        } else {
                            String::new()
                        }
                    });
                if locant.is_empty() {
                    let status = mol.name_status();
                    let name = branch_name(mol, &ring, attachment, &AtomSet::new(), side.carbon);
                    // This recursive citation only supplies a locant for the retained parent.
                    mol.set_name_status(status);
                    let name = strip_outer_parentheses(&name?);
                    if let Some((_, candidate)) =
                        name.strip_suffix("-yl").and_then(|s| s.rsplit_once('-'))
                        && candidate.as_bytes().first().is_some_and(u8::is_ascii_digit)
                        && {
                            let tail = candidate.trim_start_matches(|c: char| c.is_ascii_digit());
                            tail.is_empty()
                                || (tail.len() == 1 && tail.as_bytes()[0].is_ascii_lowercase())
                        }
                    {
                        locant = candidate.to_owned();
                    }
                }
                if !locant.is_empty() {
                    return Ok(if parent.name == "benzene" {
                        if as_ylidene {
                            "benzylidene"
                        } else {
                            "benzaldehyde"
                        }
                        .into()
                    } else if as_ylidene {
                        let stem = RULES
                            .retained
                            .substituent_stems
                            .get(&parent.name)
                            .map(|(stem, _)| stem.as_str())
                            .filter(|stem| !stem.is_empty())
                            .unwrap_or_else(|| parent.name.trim_end_matches('e'));
                        format!("{stem}-{locant}-ylmethylidene")
                    } else {
                        format!("{}-{locant}-carbaldehyde", parent.name)
                    });
                }
            }
            // A heterocycle uses its retained/Hantzsch-Widman citation rather
            // than numbering replacement prefixes from the attachment.
            if homocyclic {
                let offset = path.iter().position(|i| *i == attachment).unwrap_or(0);
                let forward: Vec<_> = path[offset..]
                    .iter()
                    .chain(&path[..offset])
                    .copied()
                    .collect();
                let reverse: Vec<_> = std::iter::once(forward[0])
                    .chain(forward[1..].iter().copied().rev())
                    .collect();
                let oriented = [forward, reverse]
                    .into_iter()
                    .min_by_key(|path| (simple_cycle_unsaturation_locants(mol, path), path.clone()))
                    .unwrap_or_default();
                let name = simple_replacement_cycle_name(mol, &oriented);
                if !name.is_empty() {
                    return Ok(if as_ylidene {
                        format!("{name}-1-ylmethylidene")
                    } else {
                        format!("{name}-1-carbaldehyde")
                    });
                }
            }
        }
    }
    if as_ylidene {
        let path = ordered_simple_ring(mol, &side.atoms);
        if !path.is_empty()
            && path.iter().all(|i| {
                mol.atoms[i].is_carbon()
                    || elements::get(&mol.atoms[i].symbol).is_some_and(|e| e.hw_stem.is_some())
            })
        {
            type RingChoice = ((Vec<usize>, usize, Vec<usize>), String);
            let mut choices: Vec<RingChoice> = Vec::new();
            let priority = path
                .iter()
                .filter(|i| !mol.atoms[i].is_carbon())
                .filter_map(|i| elements::get(&mol.atoms[i].symbol)?.hw_priority)
                .min();
            for offset in 0..path.len() {
                let forward: Vec<_> = path[offset..]
                    .iter()
                    .chain(&path[..offset])
                    .copied()
                    .collect();
                let reverse: Vec<_> = std::iter::once(forward[0])
                    .chain(forward[1..].iter().copied().rev())
                    .collect();
                for oriented in [forward, reverse] {
                    if let Some(priority) = priority
                        && (mol.atoms[&oriented[0]].is_carbon()
                            || elements::get(&mol.atoms[&oriented[0]].symbol)
                                .and_then(|e| e.hw_priority)
                                != Some(priority))
                    {
                        continue;
                    }
                    let locant = oriented.iter().position(|i| *i == side.carbon).unwrap_or(0) + 1;
                    let hetero = oriented
                        .iter()
                        .enumerate()
                        .filter_map(|(i, atom)| (!mol.atoms[atom].is_carbon()).then_some(i + 1))
                        .collect();
                    let name = simple_replacement_cycle_name(mol, &oriented);
                    if !name.is_empty() {
                        choices.push((
                            (hetero, locant, oriented),
                            format!("{name}-{locant}-ylidene"),
                        ));
                    }
                }
            }
            if let Some((_, name)) = choices.into_iter().min() {
                return Ok(name);
            }
        }
    }
    simple_carbonyl_side_name(mol, side, as_ylidene)
}

// Provenance: fn simple_carbonyl_side_name | adapted
// Source: openclatura | src/openclatura/special_cases.py::_simple_carbonyl_side_name
fn simple_carbonyl_side_name(
    mol: &Molecule,
    side: &AzineSide,
    as_ylidene: bool,
) -> Result<String, String> {
    if side.atoms.iter().any(|i| !mol.atoms[i].is_carbon()) {
        return Ok(String::new());
    }
    let edges: Vec<_> = side
        .atoms
        .iter()
        .flat_map(|a| {
            mol.get_neighbors(*a)
                .iter()
                .filter(move |b| side.atoms.contains(b) && a < b)
                .map(move |b| (*a, *b))
        })
        .collect();
    if edges.len() != side.atoms.len() - 1
        || edges
            .iter()
            .any(|(a, b)| !matches!(mol.bond_order(*a, *b), 1..=3))
    {
        return Ok(String::new());
    }
    let ends: Vec<_> = side
        .atoms
        .iter()
        .copied()
        .filter(|i| {
            mol.get_neighbors(*i)
                .iter()
                .filter(|n| side.atoms.contains(n))
                .count()
                <= 1
        })
        .collect();
    let mut path = Vec::new();
    for &start in &ends {
        for &end in &ends {
            let candidate = carbon_path(mol, &side.atoms, start, end)?;
            if candidate.contains(&side.carbon) && candidate.len() > path.len() {
                path = candidate;
            }
        }
    }
    if path.len() != side.atoms.len() {
        return Ok(String::new());
    }
    // Upstream takes the first endpoint from a Python set, then retains that
    // path orientation for unsaturation while minimizing the carbonyl locant.
    // Only name sides whose spelling is independent of that incidental order.
    let bond_orders: Vec<_> = path
        .windows(2)
        .map(|pair| mol.bond_order(pair[0], pair[1]))
        .collect();
    if !bond_orders.iter().eq(bond_orders.iter().rev()) {
        return Err("openclatura: azine unsaturation depends on upstream set ordering".into());
    }
    let locant = path.iter().position(|i| *i == side.carbon).unwrap_or(0) + 1;
    let locant = locant.min(path.len() - locant + 1);
    let stem = stems::stem_for(path.len());
    let mut parts = Vec::new();
    for (bond_order, ending) in [(2, "ene"), (3, "yne")] {
        let locants: Vec<_> = path
            .windows(2)
            .enumerate()
            .filter_map(|(i, pair)| {
                (mol.bond_order(pair[0], pair[1]) == bond_order).then_some(i + 1)
            })
            .collect();
        if !locants.is_empty() {
            let multiplier = if locants.len() == 1 {
                String::new()
            } else {
                multipliers::basic(locants.len())
            };
            parts.push(format!(
                "{}-{multiplier}{ending}",
                locants
                    .iter()
                    .map(usize::to_string)
                    .collect::<Vec<_>>()
                    .join(",")
            ));
        }
    }
    let unsaturation = if parts.is_empty() {
        "an".into()
    } else {
        format!("-{}", parts.join("-"))
    };
    Ok(if as_ylidene {
        if path.len() == 1 {
            "methylidene".into()
        } else {
            format!("{stem}{unsaturation}-{locant}-ylidene")
        }
    } else if path.len() == 1 {
        "formaldehyde".into()
    } else if locant == 1 {
        format!("{stem}{unsaturation}al")
    } else {
        format!("{stem}{unsaturation}-{locant}-one")
    })
}

// Provenance: fn sulfoxide_parent_name | adapted
// Source: openclatura | src/openclatura/special_cases.py::sulfoxide_parent_result
fn sulfoxide_parent_name(mol: &Molecule, atoms: &AtomSet) -> Result<Option<String>, String> {
    let sulfurs: Vec<_> = atoms
        .iter()
        .copied()
        .filter(|i| mol.atoms[i].symbol == "S")
        .collect();
    if sulfurs.len() != 1 {
        return Ok(None);
    }
    let central = sulfurs[0];
    let mut oxygens = Vec::new();
    let mut roots = Vec::new();
    for &n in mol.get_neighbors(central) {
        if !atoms.contains(&n) {
            return Ok(None);
        }
        if mol.atoms[&n].symbol == "O"
            && (mol.bond_order(central, n) == 2 || mol.atoms[&n].charge == -1)
        {
            oxygens.push(n);
        } else if mol.atoms[&n].is_carbon() && mol.bond_order(central, n) == 1 {
            roots.push(n);
        } else {
            return Ok(None);
        }
    }
    if oxygens.len() != 1 || roots.len() != 2 {
        return Ok(None);
    }
    let left = carbon_ligand_atoms(mol, atoms, roots[0], central);
    let right = carbon_ligand_atoms(mol, atoms, roots[1], central);
    if left.is_empty() || right.is_empty() || !left.is_disjoint(&right) {
        return Ok(None);
    }
    let ligand_results: Vec<_> = roots
        .iter()
        .map(|r| alkyl_ligand_name(mol, atoms, *r, central))
        .collect();
    if ligand_results
        .iter()
        .any(|result| result.as_ref().is_ok_and(String::is_empty))
    {
        return Ok(None);
    }
    let mut ligands: Vec<_> = ligand_results.into_iter().collect::<Result<_, _>>()?;
    // Upstream divergence: neutral sulfoxide class names require a neutral or charge-separated S=O group.
    let oxygen = oxygens[0];
    if !matches!(
        (
            mol.atoms[&central].charge,
            mol.atoms[&oxygen].charge,
            mol.bond_order(central, oxygen)
        ),
        (0, 0, 2) | (1, -1, 1)
    ) {
        return Err("openclatura: shortcut cannot preserve sulfoxide formal charges".into());
    }
    require_stereo_free_core(mol, &oxygens.into_iter().collect())?;
    // Upstream divergence: an unassigned sulfoxide raw tag is not a rendered absolute descriptor.
    if mol.atoms[&central].raw_stereo.is_some()
        && matches!(
            mol.atoms[&central].stereo.as_deref(),
            None | Some("r" | "s")
        )
    {
        return Err(SHORTCUT_STEREO_ERROR.into());
    }
    // Achiral sulfoxides take substitutive names (P-63.6); a chiral one keeps the
    // class name, as its PIN cites the descriptor before the whole name (P-93.3.4.1).
    let Some(stereo) = &mol.atoms[&central].stereo else {
        return Ok(None);
    };
    ligands.sort();
    mol.mark_name(NameStatus::NotPreferred(NameReason::FunctionalClass));
    if descriptors_omitted() {
        return Ok(Some(format!("{} sulfoxide", ligands.join(" "))));
    }
    Ok(Some(format!("({stereo})-{} sulfoxide", ligands.join(" "))))
}

// Provenance: fn has_principal_group | adapted
// Source: openclatura | src/openclatura/special_cases.py::_has_principal_group
fn has_principal_group(mol: &Molecule, atoms: &AtomSet) -> bool {
    perceive_groups(mol)
        .iter()
        .any(|g| g.is_principal_candidate && atoms.contains(&g.attachment_carbon))
}
// Provenance: fn competing_principal_group | adapted
// Source: openclatura | src/openclatura/special_cases.py::_has_principal_group
fn competing_principal_group(mol: &Molecule, atoms: &AtomSet, central: usize) -> bool {
    perceive_groups(mol).iter().any(|g| {
        g.is_principal_candidate
            && atoms.contains(&g.attachment_carbon)
            && !g.atoms_involved.contains(&central)
            && g.attachment_carbon != central
    })
}

// Provenance: fn terminal_ligand_name | adapted
// Source: openclatura | src/openclatura/special_cases.py::_terminal_ligand_name
fn terminal_ligand_name(mol: &Molecule, idx: usize) -> String {
    // Upstream divergence: neutral terminal prefixes do not represent formal charge or atom stereo.
    if mol.atoms[&idx].charge != 0 || mol.atoms[&idx].has_stereo() {
        return String::new();
    }
    match mol.atoms[&idx].symbol.as_str() {
        "F" => "fluoro".into(),
        "Cl" => "chloro".into(),
        "Br" => "bromo".into(),
        "I" => "iodo".into(),
        "C" if mol.degree(idx) == 1 => "methyl".into(),
        _ => String::new(),
    }
}

// Provenance: fn central_ligand_name | adapted
// Source: openclatura | src/openclatura/special_cases.py::_central_ligand_name
fn central_ligand_name(
    mol: &Molecule,
    atoms: &AtomSet,
    central: usize,
    neighbor: usize,
) -> Result<String, String> {
    if mol.bond_order(central, neighbor) != 1 {
        return Ok(String::new());
    }
    let name = terminal_ligand_name(mol, neighbor);
    if !name.is_empty() {
        return Ok(name);
    }
    if mol.atoms[&neighbor].symbol != "O" || mol.atoms[&neighbor].charge != 0 {
        return Ok(String::new());
    }
    let roots: Vec<_> = mol
        .get_neighbors(neighbor)
        .iter()
        .copied()
        .filter(|n| atoms.contains(n) && *n != central)
        .collect();
    if roots.len() != 1 {
        return Ok(String::new());
    }
    let alkyl = alkyl_ligand_name(mol, atoms, roots[0], neighbor)?;
    if alkyl.is_empty() {
        Ok(String::new())
    } else {
        let mut exclude = outside(mol, atoms);
        exclude.insert(neighbor);
        Ok(
            name_subgraph_rendered(mol, roots[0], &exclude, Some(neighbor), false)?
                .into_oxy_prefix(),
        )
    }
}

// Provenance: fn simple_central_parent_hydride_name | adapted
// Source: openclatura | src/openclatura/special_cases.py::simple_central_parent_hydride_result
fn simple_central_parent_hydride_name(
    mol: &Molecule,
    atoms: &AtomSet,
) -> Result<Option<StructuralParentName>, String> {
    let mut centers: Vec<_> = atoms
        .iter()
        .copied()
        .filter(|i| {
            // P-72.3: iodane takes a hydride as well (diphenyliodanuide).
            (RULES
                .components
                .mononuclear_parent_hydrides
                .contains_key(&mol.atoms[i].symbol)
                || mol.atoms[i].symbol == "I" && mol.atoms[i].charge == -1)
                && (mol.atoms[i].symbol != "O" || mol.atoms[i].charge == 1)
                && (mol.degree(*i) >= 2
                    || mol.degree(*i) == 1
                        && mol.atoms[i].charge == 1
                        && matches!(mol.atoms[i].symbol.as_str(), "O" | "S" | "Se" | "Te"))
        })
        .collect();
    // P-41, P-73.1.2: a cation outranks every neutral class, so a cationic
    // chalcogen centre is the parent over neutral heteroatoms such as an amine.
    // An anion is senior to cations and neutral classes alike.
    if centers.len() > 1 {
        centers.retain(|i| {
            mol.atoms[i].charge == -1
                || mol.atoms[i].charge == 1
                    && matches!(mol.atoms[i].symbol.as_str(), "O" | "S" | "Se" | "Te")
        });
    }
    if centers.len() != 1 {
        return Ok(None);
    }
    let c = centers[0];
    let a = &mol.atoms[&c];
    if mol.atoms[&c].has_stereo() {
        return Ok(None);
    }
    let senior = matches!(
        a.symbol.as_str(),
        "P" | "As" | "Sb" | "Bi" | "Si" | "Ge" | "Sn" | "Pb" | "B" | "Al" | "Ga" | "In" | "Tl"
    );
    let is_chalcogen = |n: usize| matches!(mol.atoms[&n].symbol.as_str(), "O" | "S" | "Se" | "Te");
    let organyl = |n: usize| {
        mol.get_neighbors(c)
            .iter()
            .all(|&o| o == n || mol.atoms[&o].is_carbon())
    };
    // P+-O- and its S, Se and Te analogues stand for P=O.
    let oxide = a.charge == 1
        && matches!(a.symbol.as_str(), "P" | "As" | "Sb")
        && mol
            .get_neighbors(c)
            .iter()
            .any(|&n| is_chalcogen(n) && mol.atoms[&n].charge == -1 && mol.degree(n) == 1);
    // P-73.1.1.1/.2: protonated chalcogen hydrides retain the systematic
    // oxidanium, sulfanium, selanium or tellanium parent, including remaining
    // hydrogens. P-73.1.2.1: ylidene and ylidyne ligands substitute the same
    // parent, as in ethyl(propan-2-ylidene)oxidanium.
    // Other ionic or radical centers require the ordinary parent selection.
    // The isotope adapter needs atom locants, which this shortcut lacks.
    let chalcogen_onium = !mol.is_isotope_parent()
        && a.charge == 1
        && matches!(a.symbol.as_str(), "O" | "S" | "Se" | "Te")
        && mol
            .get_neighbors(c)
            .iter()
            .map(|n| mol.bond_order(c, *n))
            .sum::<usize>()
            + a.total_h_count.max(a.explicit_h_count)
            == 3
        // A heteroatom on an ylidene or ylidyne carbon can need a retained
        // parent such as uronium (P-73.1.2.2) or meet prefixes that assume a
        // single attachment, so those keep the ordinary route.
        && mol.get_neighbors(c).iter().all(|&n| {
            mol.atoms[&n].is_carbon()
                && (mol.bond_order(c, n) == 1
                    || mol
                        .get_neighbors(n)
                        .iter()
                        .all(|&m| m == c || mol.atoms[&m].is_carbon()))
        })
        // Ligand names cannot carry E/Z of a double bond at the centre.
        && !mol
            .bonds
            .values()
            .any(|b| (b.u == c || b.v == c) && b.stereo.is_some())
        && mol
            .atoms
            .values()
            .all(|atom| atom.radical_electrons == 0 && (atom.idx == c || atom.charge == 0))
        && !get_cyclic_atoms(mol, &outside(mol, atoms)).contains(&c);
    let onium = chalcogen_onium
        || a.charge == 1
            && senior
            && !oxide
            && a.total_h_count == 0
            && mol
                .get_neighbors(c)
                .iter()
                .all(|n| mol.atoms[n].is_carbon());
    // P-72.3: a hydride added to the neutral hydride gives the suffix "uide"
    // (tetramethylboranuide).
    let uide = a.charge == -1
        && (senior || a.symbol == "I")
        && mol
            .get_neighbors(c)
            .iter()
            .all(|&n| atoms.contains(&n) && mol.bond_order(c, n) == 1 && mol.atoms[&n].charge == 0)
        && elements::get(&a.symbol)
            .is_some_and(|e| e.standard_valence + 1 == mol.degree(c) + a.total_h_count);
    if a.charge != 0 && !oxide && !onium && !uide {
        return Ok(None);
    }
    // P-41: the anion outranks any characteristic group on its ligands.
    let hydrocarbyl = chalcogen_onium || uide || senior && !has_principal_group(mol, atoms);
    let mut names = Vec::new();
    let mut chalcogen = None;
    for &n in mol.get_neighbors(c) {
        if !atoms.contains(&n) {
            continue;
        }
        if hydrocarbyl
            && chalcogen.is_none()
            && matches!(a.symbol.as_str(), "P" | "As" | "Sb")
            && a.total_h_count == 0
            && is_chalcogen(n)
            && mol.degree(n) == 1
            && (mol.bond_order(c, n) == 2
                || (mol.bond_order(c, n) == 1 && mol.atoms[&n].charge == -1))
            && (mol.atoms[&n].symbol == "O" || organyl(n))
        {
            chalcogen = Some(n);
            continue;
        }
        let mut name = central_ligand_name(mol, atoms, c, n)?;
        if name.is_empty()
            && hydrocarbyl
            && mol.atoms[&n].is_carbon()
            && (chalcogen_onium || mol.bond_order(c, n) == 1)
        {
            let ligand = component_atoms_until_blocked(mol, atoms, n, &set([c]));
            if !ligand.is_empty()
                && !ligand.iter().any(|i| mol.atoms[i].charge != 0)
                && !mol
                    .get_neighbors(c)
                    .iter()
                    .any(|other| *other != n && ligand.contains(other))
            {
                name = if chalcogen_onium {
                    let mut exclude = outside(mol, atoms);
                    exclude.insert(c);
                    name_subgraph_rendered(mol, n, &exclude, Some(c), true)?.text
                } else {
                    branch_name(mol, atoms, n, &set([c]), c)?
                };
            }
        }
        if name.is_empty() {
            return Ok(None);
        }
        names.push(name);
    }
    if names.is_empty() {
        return Ok(None);
    }
    let prefix = center_ligands(&names);
    // Halides, pseudohalides and esters are named from the acids (P-68.3.2.3.1),
    // so they keep the oxide class name.
    let class_oxide = chalcogen.is_some_and(|n| !organyl(n));
    let valence: usize = mol
        .get_neighbors(c)
        .iter()
        .map(|n| mol.bond_order(c, *n))
        .sum::<usize>()
        + a.total_h_count.max(a.explicit_h_count)
        + usize::from(oxide)
        - if class_oxide { 2 } else { 0 };
    let lambda =
        if onium || uide || elements::get(&a.symbol).is_some_and(|e| e.standard_valence == valence)
        {
            String::new()
        } else {
            format!("lambda{valence}-")
        };
    let mut parent = format!(
        "{lambda}{}",
        RULES
            .components
            .mononuclear_parent_hydrides
            .get(&a.symbol)
            .map_or("iodane", String::as_str)
    );
    if class_oxide {
        parent.push_str(" oxide");
        mol.mark_name(NameStatus::NotPreferred(NameReason::FunctionalClass));
    } else if let Some(n) = chalcogen {
        // The suffix is preferred to the oxide or sulfide class name (P-68.3.2.3.1).
        match mol.atoms[&n].symbol.as_str() {
            "O" => {
                parent.pop();
                parent.push_str("one");
            }
            "S" => parent.push_str("thione"),
            "Se" => parent.push_str("selone"),
            _ => parent.push_str("tellone"),
        }
    } else if onium {
        if parent.ends_with('e') {
            parent.pop();
        }
        parent.push_str("ium");
    } else if uide {
        parent.pop();
        parent.push_str("uide");
    }
    if a.symbol == "N" || a.symbol == "S" && !onium && lambda.is_empty() {
        // Methyl-only amines and sulfides take substitutive names (P-62.2.2.1, P-63.2.5).
        if mol
            .get_neighbors(c)
            .iter()
            .all(|n| mol.atoms[n].is_carbon())
        {
            return Ok(None);
        }
        mol.mark_name(NameStatus::NotPreferred(NameReason::HydrideParent));
    }
    // P-63.2.5: (methyltellanyl)methane, not dimethyltellane; λ4/λ6 hydrides
    // and onium cations keep the hydride name (P-68.4.3.1, P-73.1.1).
    if matches!(a.symbol.as_str(), "S" | "Se" | "Te")
        && !onium
        && lambda.is_empty()
        && mol
            .get_neighbors(c)
            .iter()
            .any(|n| mol.atoms[n].is_carbon())
    {
        mol.mark_name(NameStatus::NotPreferred(NameReason::HydrideParent));
    }
    // P-44.1.2 and P-68.3.2.3.1: senior hydrides with organyl ligands,
    // including a terminal chalcogen suffix on P, As or Sb. The charge-separated
    // P+-O- representation (and its analogues) has the same ligand role as P=O.
    // An anion cited over its ligands' characteristic groups is not audited.
    let certified = chalcogen_onium
        || !(uide && has_principal_group(mol, atoms))
            && matches!(
                a.symbol.as_str(),
                "B" | "Si" | "Ge" | "Sn" | "Pb" | "P" | "As" | "Sb" | "Bi"
            )
            && mol.get_neighbors(c).iter().all(|&n| {
                mol.bond_order(c, n) == 1 && mol.atoms[&n].is_carbon()
                    || (mol.bond_order(c, n) == 2
                        || oxide && mol.bond_order(c, n) == 1 && mol.atoms[&n].charge == -1)
                        && matches!(a.symbol.as_str(), "P" | "As" | "Sb")
                        && is_chalcogen(n)
                        && mol.degree(n) == 1
            });
    let name = if !prefix.is_empty() && !lambda.is_empty() {
        format!("{prefix}-{parent}")
    } else {
        format!("{prefix}{parent}")
    };
    Ok(Some(StructuralParentName {
        certified,
        ..StructuralParentName::ordinary(name)
    }))
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum OxoLigandRole {
    Oxo,
    Oxido,
    Hydroxy,
    Alkoxy,
    Peroxy,
}

#[derive(Clone, Debug)]
struct OxoLigand {
    oxygen: usize,
    role: OxoLigandRole,
    attachment: Option<usize>,
}
#[derive(Clone, Debug)]
struct OxoRole {
    central: usize,
    ligands: Vec<OxoLigand>,
}
impl OxoRole {
    // Provenance: fn OxoRole::count | adapted
    // Source: openclatura | src/openclatura/oxoacid_roles.py::CentralOxoRole.count
    fn count(&self, role: OxoLigandRole) -> usize {
        self.ligands.iter().filter(|l| l.role == role).count()
    }
    // Provenance: fn OxoRole::core | adapted
    // Source: openclatura | src/openclatura/oxoacid_roles.py::CentralOxoRole fields central/ligands.oxygen
    fn core(&self) -> AtomSet {
        std::iter::once(self.central)
            .chain(self.ligands.iter().map(|l| l.oxygen))
            .collect()
    }
}
// Provenance: fn oxo_roles | adapted
// Source: openclatura | src/openclatura/oxoacid_roles.py::central_oxo_roles
fn oxo_roles(mol: &Molecule, atoms: &AtomSet) -> Vec<OxoRole> {
    let mut roles = Vec::new();
    for &central in atoms {
        if mol.get_neighbors(central).is_empty()
            || mol
                .get_neighbors(central)
                .iter()
                .any(|n| !atoms.contains(n) || mol.atoms[n].symbol != "O")
        {
            continue;
        }
        let mut ligands = Vec::new();
        let mut valid = true;
        for &oxygen in mol.get_neighbors(central) {
            let normalized = matches!(mol.atoms[&central].symbol.as_str(), "Cl" | "Br" | "I")
                && mol.bond_order(central, oxygen) == 1
                && mol.atoms[&oxygen].charge < 0
                && mol.degree(oxygen) == 1;
            let (role, attachment) = if mol.bond_order(central, oxygen) == 2 || normalized {
                (OxoLigandRole::Oxo, None)
            } else {
                let neighbors: Vec<_> = mol
                    .get_neighbors(oxygen)
                    .iter()
                    .copied()
                    .filter(|n| *n != central)
                    .collect();
                if neighbors.is_empty() {
                    (
                        if mol.atoms[&oxygen].charge < 0 {
                            OxoLigandRole::Oxido
                        } else {
                            OxoLigandRole::Hydroxy
                        },
                        None,
                    )
                } else if neighbors.len() == 1 && atoms.contains(&neighbors[0]) {
                    match mol.atoms[&neighbors[0]].symbol.as_str() {
                        "C" => (OxoLigandRole::Alkoxy, Some(neighbors[0])),
                        "O" => (OxoLigandRole::Peroxy, Some(neighbors[0])),
                        _ => {
                            valid = false;
                            break;
                        }
                    }
                } else {
                    valid = false;
                    break;
                }
            };
            ligands.push(OxoLigand {
                oxygen,
                role,
                attachment,
            });
        }
        if valid {
            roles.push(OxoRole { central, ligands });
        }
    }
    roles
}
// Provenance: fn matching_oxo_spec | adapted
// Source: openclatura | src/openclatura/special_cases.py::_matching_oxoacid_spec
fn matching_oxo_spec<'a>(
    mol: &Molecule,
    role: &OxoRole,
) -> Option<&'a super::nomenclature::OxoacidSpec>
where
    'static: 'a,
{
    let atom = &mol.atoms[&role.central];
    if atom.total_h_count != 0 {
        return None;
    }
    let single = role
        .ligands
        .iter()
        .filter(|l| l.role != OxoLigandRole::Oxo)
        .count();
    let double = role.count(OxoLigandRole::Oxo);
    RULES
        .components
        .replacement_parent_oxoacid_specs
        .iter()
        .find(|spec| {
            spec.central == atom.symbol
                && spec.single_o == single
                && spec.double_o == double
                && spec.charge.is_none_or(|c| c == atom.charge)
        })
}
// Provenance: fn oxoacid_parent_name | adapted
// Source: openclatura | src/openclatura/special_cases.py::oxoacid_parent_result
fn oxoacid_parent_name(mol: &Molecule, atoms: &AtomSet) -> Option<String> {
    let roles = oxo_roles(mol, atoms);
    if roles.len() != 1 {
        return None;
    }
    let r = &roles[0];
    if r.count(OxoLigandRole::Alkoxy) > 0
        || r.count(OxoLigandRole::Peroxy) > 0
        || r.core() != *atoms
    {
        return None;
    }
    let spec = matching_oxo_spec(mol, r)?;
    if r.count(OxoLigandRole::Oxido) == 0 || spec.ester_suffix.is_empty() {
        return Some(spec.name.clone());
    }
    // P-65.6.2.3.2, P-72.2.2.2.1.2: an anion that keeps acid hydrogen, such
    // as HCO3-, is `hydrogen carbonate`, and H2PO4- `dihydrogen phosphate`.
    // A zwitterion such as O[N+](=O)[O-] has no net charge and stays the acid.
    Some(match r.count(OxoLigandRole::Hydroxy) {
        0 => spec.ester_suffix.clone(),
        _ if atoms.iter().map(|a| mol.atoms[a].charge).sum::<i32>() >= 0 => spec.name.clone(),
        1 => format!("hydrogen {}", spec.ester_suffix),
        n => format!("{}hydrogen {}", multipliers::basic(n), spec.ester_suffix),
    })
}
// Provenance: fn modifier_phrase | adapted
// Source: openclatura | src/openclatura/special_cases.py::_ester_modifier_phrase
fn modifier_phrase(names: &[String], safe: bool) -> String {
    let mut counts: BTreeMap<String, usize> = BTreeMap::new();
    for name in names {
        *counts.entry(name.clone()).or_default() += 1;
    }
    let mut counted: Vec<_> = counts.into_iter().collect();
    counted.sort_by_key(|(name, _)| substituent_sort_key(name));
    counted
        .into_iter()
        .map(|(name, count)| {
            if count == 1 {
                name
            } else {
                format_multiplier(&name, count, safe)
            }
        })
        .collect::<Vec<_>>()
        .join(" ")
}
// Provenance: fn ester_modifier_name | adapted
// Source: openclatura | src/openclatura/special_cases.py::_ester_modifier_name
fn ester_modifier_name(
    mol: &Molecule,
    atoms: &AtomSet,
    root: usize,
    oxygen: usize,
    acid: &AtomSet,
) -> Result<String, String> {
    let name = branch_name(mol, atoms, root, acid, oxygen)?;
    if !name.is_empty() {
        return Ok(strip_outer_parentheses(&name));
    }
    alkyl_ligand_name(mol, atoms, root, oxygen)
}
// Provenance: fn oxoacid_ester_name | adapted
// Source: openclatura | src/openclatura/special_cases.py::oxoacid_ester_result
fn oxoacid_ester_name(
    mol: &Molecule,
    atoms: &AtomSet,
) -> Result<Option<StructuralParentName>, String> {
    let status = mol.name_status();
    let cyclic_bonds: BTreeSet<_> = biconnected_edge_components(mol, atoms)
        .into_iter()
        .filter(|block| block.len() > 1)
        .flatten()
        .collect();
    let groups = perceive_groups(mol);
    // P-67.1.3.2 / P-67.1.3.3: acid-derived groups form anhydrides, not esters;
    // a principal group inside the ligand takes the suffix instead (P-67.1.4.1.3).
    let alcohol_derived = |root: usize, ligand: &AtomSet| {
        mol.atoms[&root].is_carbon()
            && mol.get_neighbors(root).iter().all(|&n| {
                mol.bond_order(root, n) == 1
                    || mol.atoms[&n].is_carbon()
                    || cyclic_bonds.contains(&normalize_edge(root, n))
            })
            && !groups
                .iter()
                .any(|g| g.is_principal_candidate && g.atom_ids().is_subset(ligand))
    };
    let mut matches = Vec::new();
    for role in oxo_roles(mol, atoms) {
        mol.set_name_status(status);
        if mol.atoms[&role.central].has_stereo() {
            continue;
        }
        let core = role.core();
        let central = &mol.atoms[&role.central];
        if role.count(OxoLigandRole::Peroxy) > 0 {
            // Verified charge-normalized halogen peroxy ester templates.
            if role.count(OxoLigandRole::Peroxy) != 1
                || role.count(OxoLigandRole::Oxo) != 2
                || role.ligands.len() != 3
                || central.charge != 2
            {
                continue;
            }
            let suffix = match central.symbol.as_str() {
                "Cl" => "peroxychlorate",
                "Br" => "peroxybromate",
                "I" => "peroxyiodate",
                _ => continue,
            };
            let Some(peroxy) = role
                .ligands
                .iter()
                .find(|l| l.role == OxoLigandRole::Peroxy)
            else {
                continue;
            };
            let Some(distal) = peroxy.attachment else {
                continue;
            };
            let roots: Vec<_> = mol
                .get_neighbors(distal)
                .iter()
                .copied()
                .filter(|n| atoms.contains(n) && *n != role.central && *n != peroxy.oxygen)
                .collect();
            if roots.len() != 1 {
                continue;
            }
            let mut acid = core.clone();
            acid.insert(distal);
            let modifier_atoms = component_atoms_until_blocked(mol, atoms, roots[0], &acid);
            if modifier_atoms.is_empty()
                || modifier_atoms.union(&acid).copied().collect::<AtomSet>() != *atoms
            {
                continue;
            }
            let modifier = ester_modifier_name(mol, atoms, roots[0], distal, &acid)?;
            if !modifier.is_empty() {
                require_stereo_free_core(mol, &acid)?;
                matches.push((
                    StructuralParentName {
                        certified: alcohol_derived(roots[0], &modifier_atoms),
                        component_splits: BTreeSet::from([normalize_edge(distal, roots[0])]),
                        ..StructuralParentName::ordinary(format!("{modifier} {suffix}"))
                    },
                    mol.name_status(),
                ));
            }
            continue;
        }
        let Some(spec) = matching_oxo_spec(mol, &role) else {
            continue;
        };
        let base = spec.ester_suffix.as_str();
        if base.is_empty() {
            continue;
        }
        let esters: Vec<_> = role
            .ligands
            .iter()
            .filter(|l| l.role == OxoLigandRole::Alkoxy)
            .collect();
        if esters.is_empty() {
            continue;
        }
        let mut covered = AtomSet::new();
        let mut names = Vec::new();
        let mut component_splits = BTreeSet::new();
        let mut certified = true;
        let mut valid = true;
        for ligand in esters {
            let Some(root) = ligand.attachment else {
                valid = false;
                break;
            };
            let ligand_atoms = component_atoms_until_blocked(mol, atoms, root, &core);
            if ligand_atoms.is_empty() || !ligand_atoms.is_disjoint(&covered) {
                valid = false;
                break;
            }
            let name = ester_modifier_name(mol, atoms, root, ligand.oxygen, &core)?;
            if name.is_empty() {
                valid = false;
                break;
            }
            names.push(name);
            certified &= alcohol_derived(root, &ligand_atoms);
            component_splits.insert(normalize_edge(ligand.oxygen, root));
            covered.extend(ligand_atoms);
        }
        if !valid || covered.union(&core).copied().collect::<AtomSet>() != *atoms {
            continue;
        }
        require_stereo_free_core(mol, &core)?;
        let hydroxy = role.count(OxoLigandRole::Hydroxy);
        let suffix = if hydroxy == 0 || spec.charge.is_some() {
            base.to_owned()
        } else {
            let hydrogen = if hydroxy == 1 {
                "hydrogen".into()
            } else {
                format!("{}hydrogen", multipliers::basic(hydroxy))
            };
            format!("{hydrogen} {base}")
        };
        matches.push((
            StructuralParentName {
                certified,
                component_splits,
                ..StructuralParentName::ordinary(format!(
                    "{} {suffix}",
                    modifier_phrase(&names, false)
                ))
            },
            mol.name_status(),
        ));
    }
    mol.set_name_status(status);
    Ok(if matches.len() == 1 {
        matches.pop().map(|(name, status)| {
            mol.set_name_status(status);
            name
        })
    } else {
        None
    })
}

// Provenance: fn organophosphorus_acid_name | adapted
// Source: openclatura | src/openclatura/special_cases.py::organophosphinic_acid_result
// Source: openclatura | src/openclatura/special_cases.py::organophosphonic_acid_result
/// P-67.1.1.2: organyl groups substitute the hydrogen on the central atom of
/// the phosphorus, arsenic and antimony acids of P-67.1.1.1. Their anions and
/// esters change `ic acid` to `ate` and `ous acid` to `ite` (P-67.1.3.1,
/// P-67.1.3.2, P-72.2.2.2.1).
fn organophosphorus_acid_name(
    mol: &Molecule,
    atoms: &AtomSet,
) -> Result<Option<StructuralParentName>, String> {
    let centers: Vec<_> = atoms
        .iter()
        .copied()
        .filter(|i| matches!(mol.atoms[i].symbol.as_str(), "P" | "As" | "Sb"))
        .collect();
    if centers.len() != 1 {
        return Ok(None);
    }
    let c = centers[0];
    let a = &mol.atoms[&c];
    if get_cyclic_atoms(mol, &AtomSet::new()).contains(&c)
        || mol.atoms[&c].has_stereo()
        || competing_principal_group(mol, atoms, c)
    {
        return Ok(None);
    }
    let (mut double, mut hydroxy, mut oxido, mut esters, mut roots) =
        (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
    for &n in mol.get_neighbors(c) {
        if !atoms.contains(&n) {
            return Ok(None);
        }
        let atom = &mol.atoms[&n];
        let bond = mol.bond_order(c, n);
        if atom.symbol == "O" && bond == 2 {
            double.push(n);
        } else if atom.symbol == "O" && bond == 1 && atom.charge == -1 && mol.degree(n) == 1 {
            oxido.push(n);
        } else if atom.symbol == "O" && bond == 1 && atom.charge == 0 && mol.degree(n) == 1 {
            hydroxy.push(n);
        } else if atom.symbol == "O" && bond == 1 && atom.charge == 0 && mol.degree(n) == 2 {
            esters.push(n);
        } else if atom.is_carbon() && bond == 1 {
            roots.push(n);
        } else {
            return Ok(None);
        }
    }
    // The charge-separated P(+)-O(-) form of a phosphonic acid or ester.
    if a.charge == 1 && a.symbol == "P" && double.is_empty() && oxido.len() == 1 {
        double.append(&mut oxido);
    } else if a.charge != 0 {
        return Ok(None);
    }
    let acidic = hydroxy.len() + oxido.len() + esters.len();
    // One acidic oxygen is the `in` acid and two the `on` acid; the trivalent
    // acids and esters without an anionic oxygen keep their substitutive names.
    if !(1..=2).contains(&acidic)
        || roots.is_empty()
        || roots.len() + a.total_h_count.max(a.explicit_h_count) + acidic != 3
        || double.len() > 1
        || a.charge == 1 && acidic != 2
        || double.is_empty() && oxido.is_empty()
    {
        return Ok(None);
    }
    let stem = format!(
        "{}{}",
        match a.symbol.as_str() {
            "P" => "phosph",
            "As" => "ars",
            _ => "stib",
        },
        if acidic == 1 { "in" } else { "on" }
    );
    let oxo = !double.is_empty();
    let acid: AtomSet = std::iter::once(c)
        .chain(double)
        .chain(hydroxy.iter().copied())
        .chain(oxido.iter().copied())
        .chain(esters.iter().copied())
        .collect();
    require_stereo_free_core(mol, &acid)?;
    let mut ligands = Vec::new();
    for root in roots {
        let name = ligand_name(mol, atoms, root, c, &acid)?;
        if name.is_empty() {
            return Ok(None);
        }
        ligands.push(name);
    }
    // The `in` acids group their carbon ligands as a centre prefix. The
    // single `on` acid ligand is cited directly by the upstream route.
    let ligand = if acidic == 1 {
        center_ligands(&ligands)
    } else {
        ligands.remove(0)
    };
    let mut modifiers = Vec::new();
    let mut component_splits = BTreeSet::new();
    for oxygen in esters {
        let Some(root) = mol.get_neighbors(oxygen).iter().copied().find(|n| *n != c) else {
            return Ok(None);
        };
        let name = ester_modifier_name(mol, atoms, root, oxygen, &acid)?;
        if name.is_empty() {
            return Ok(None);
        }
        modifiers.push(name);
        component_splits.insert(normalize_edge(oxygen, root));
    }
    if modifiers.is_empty() && oxido.is_empty() {
        return Ok(Some(StructuralParentName::ordinary(format!(
            "{ligand}{stem}ic acid"
        ))));
    }
    // P-67.1.3.2, P-72.2.2.2.1.2: ester groups, then `hydrogen` for each
    // remaining hydroxy group, then the anion name.
    let mut words = Vec::new();
    if !modifiers.is_empty() {
        words.push(modifier_phrase(&modifiers, acidic == 2));
    }
    if !hydroxy.is_empty() {
        words.push("hydrogen".into());
    }
    words.push(format!("{ligand}{stem}{}", if oxo { "ate" } else { "ite" }));
    Ok(Some(StructuralParentName {
        component_splits,
        ..StructuralParentName::ordinary(words.join(" "))
    }))
}

// Provenance: fn organoboronic_acid_name | adapted
// Source: openclatura | src/openclatura/special_cases.py::organoboronic_acid_result
fn organoboronic_acid_name(
    mol: &Molecule,
    atoms: &AtomSet,
) -> Result<Option<StructuralParentName>, String> {
    let centers: Vec<_> = atoms
        .iter()
        .copied()
        .filter(|i| mol.atoms[i].symbol == "B")
        .collect();
    if centers.len() != 1 {
        return Ok(None);
    }
    let c = centers[0];
    let a = &mol.atoms[&c];
    if a.charge != 0
        || a.total_h_count != 0
        || get_cyclic_atoms(mol, &AtomSet::new()).contains(&c)
        || competing_principal_group(mol, atoms, c)
    {
        return Ok(None);
    }
    let (mut hydroxy, mut roots) = (Vec::new(), Vec::new());
    for &n in mol.get_neighbors(c) {
        if !atoms.contains(&n) || mol.bond_order(c, n) != 1 {
            return Ok(None);
        }
        let a = &mol.atoms[&n];
        if a.symbol == "O" && a.charge == 0 && mol.degree(n) == 1 {
            hydroxy.push(n);
        } else if a.is_carbon() {
            roots.push(n);
        } else {
            return Ok(None);
        }
    }
    if !matches!((hydroxy.len(), roots.len()), (2, 1) | (1, 2)) {
        return Ok(None);
    }
    let certified = hydroxy.iter().all(|i| mol.atoms[i].total_h_count == 1);
    let acid: AtomSet = std::iter::once(c).chain(hydroxy).collect();
    require_stereo_free_core(mol, &acid)?;
    let mut names = Vec::new();
    for root in roots {
        let name = ligand_name(mol, atoms, root, c, &acid)?;
        if name.is_empty() {
            return Ok(None);
        }
        names.push(name);
    }
    Ok(Some(StructuralParentName {
        // P-68.1.4.1: the non-carbon ligands are hydroxy groups.
        certified,
        ..StructuralParentName::ordinary(format!(
            "{}{} acid",
            center_ligands(&names),
            if names.len() == 1 {
                "boronic"
            } else {
                "borinic"
            }
        ))
    }))
}

// Provenance: fn connected_branches | adapted
// Source: openclatura | src/openclatura/special_cases.py::_connected_branches
fn connected_branches(mol: &Molecule, atoms: &AtomSet) -> Vec<AtomSet> {
    let mut branches = Vec::new();
    let mut seen = AtomSet::new();
    for &start in atoms {
        if seen.contains(&start) {
            continue;
        }
        let mut branch = set([start]);
        seen.insert(start);
        let mut stack = vec![start];
        while let Some(idx) = stack.pop() {
            for &n in mol.get_neighbors(idx) {
                if atoms.contains(&n) && seen.insert(n) {
                    branch.insert(n);
                    stack.push(n);
                }
            }
        }
        branches.push(branch);
    }
    branches
}
// Provenance: fn ordered_backbone | adapted
// Source: openclatura | src/openclatura/special_cases.py::_ordered_backbone_chain
fn ordered_backbone(mol: &Molecule, atoms: &AtomSet) -> Option<Vec<usize>> {
    let degrees: BTreeMap<_, _> = atoms
        .iter()
        .map(|i| {
            (
                *i,
                mol.get_neighbors(*i)
                    .iter()
                    .filter(|n| atoms.contains(n))
                    .count(),
            )
        })
        .collect();
    if degrees.values().any(|d| *d > 2) {
        return None;
    }
    let ends: Vec<_> = degrees
        .iter()
        .filter(|(_, d)| **d <= 1)
        .map(|(i, _)| *i)
        .collect();
    if ends.len() != 2 {
        return None;
    }
    let mut chain = vec![*ends.iter().min()?];
    let mut previous = None;
    while chain.len() < atoms.len() {
        let current = *chain.last()?;
        let neighbors: Vec<_> = mol
            .get_neighbors(current)
            .iter()
            .copied()
            .filter(|n| atoms.contains(n) && Some(*n) != previous)
            .collect();
        if neighbors.len() != 1 || chain.contains(&neighbors[0]) {
            return None;
        }
        previous = Some(current);
        chain.push(neighbors[0]);
    }
    Some(chain)
}
// Provenance: fn chain_oxide_oxygens | adapted
// Source: openclatura | src/openclatura/special_cases.py::_chain_oxide_oxygens
fn chain_oxide_oxygens(
    mol: &Molecule,
    chain: &[usize],
    atoms: &AtomSet,
) -> Option<BTreeMap<usize, usize>> {
    let mut oxides = BTreeMap::new();
    for &i in chain {
        if mol.atoms[&i].charge == 0 {
            continue;
        }
        if mol.atoms[&i].charge != 1 {
            return None;
        }
        let oxygens: Vec<_> = mol
            .get_neighbors(i)
            .iter()
            .copied()
            .filter(|n| {
                atoms.contains(n)
                    && mol.atoms[n].symbol == "O"
                    && mol.atoms[n].charge == -1
                    && mol.degree(*n) == 1
                    && mol.bond_order(i, *n) == 1
            })
            .collect();
        if oxygens.len() != 1 {
            return None;
        }
        oxides.insert(i, oxygens[0]);
    }
    (!oxides.is_empty()).then_some(oxides)
}
// Provenance: fn homonuclear_backbone | adapted
// Source: openclatura | src/openclatura/special_cases.py::_homonuclear_backbone
fn homonuclear_backbone(mol: &Molecule, atoms: &AtomSet) -> Option<(String, Vec<usize>)> {
    let mut candidates = Vec::new();
    let symbols: BTreeSet<_> = atoms.iter().map(|i| mol.atoms[i].symbol.clone()).collect();
    for symbol in symbols {
        if !RULES
            .components
            .mononuclear_parent_hydrides
            .contains_key(&symbol)
            || matches!(symbol.as_str(), "F" | "Cl" | "Br" | "I")
        {
            continue;
        }
        let backbone: AtomSet = atoms
            .iter()
            .copied()
            .filter(|i| mol.atoms[i].symbol == symbol)
            .collect();
        if backbone.len() < 2 {
            continue;
        }
        if let Some(chain) = ordered_backbone(mol, &backbone) {
            let rank = ["N", "P", "Si", "B", "O", "S", "C"]
                .iter()
                .position(|s| *s == symbol)
                .unwrap_or(8);
            candidates.push((rank, symbol, chain));
        }
    }
    if let Some((_, symbol, chain)) = candidates.into_iter().min() {
        return Some((symbol, chain));
    }
    // A lone formal oxide centre may lie within a connected same-element chain
    // even when other disconnected atoms of that element occur in its ligands.
    let centers: Vec<_> = atoms
        .iter()
        .copied()
        .filter(|i| {
            RULES
                .components
                .mononuclear_parent_hydrides
                .contains_key(&mol.atoms[i].symbol)
                && !matches!(mol.atoms[i].symbol.as_str(), "O" | "F" | "Cl" | "Br" | "I")
                && mol.atoms[i].charge == 1
                && chain_oxide_oxygens(mol, &[*i], atoms).is_some()
        })
        .collect();
    if centers.len() != 1 {
        return None;
    }
    let c = centers[0];
    let symbol = mol.atoms[&c].symbol.clone();
    let mut connected = set([c]);
    let mut stack = vec![c];
    while let Some(i) = stack.pop() {
        for &n in mol.get_neighbors(i) {
            if atoms.contains(&n) && mol.atoms[&n].symbol == symbol && connected.insert(n) {
                stack.push(n);
            }
        }
    }
    if connected.len() < 2 {
        return None;
    }
    Some((symbol, ordered_backbone(mol, &connected)?))
}
// Provenance: fn same_element_parent_name | adapted
// Source: openclatura | src/openclatura/special_cases.py::_same_element_parent_name
fn same_element_parent_name(symbol: &str, length: usize, orders: &[usize]) -> Option<String> {
    let parent = RULES.components.mononuclear_parent_hydrides.get(symbol)?;
    if length < 2 || !parent.ends_with("ane") {
        return None;
    }
    if !orders.iter().any(|o| *o > 1)
        && let Some(retained) = RULES
            .components
            .retained_homonuclear_chain_names
            .get(&format!("{symbol}{length}"))
    {
        return Some(retained.clone());
    }
    let base = format!("{}{parent}", multipliers::basic(length));
    let unsaturated: Vec<_> = orders.iter().enumerate().filter(|(_, o)| **o > 1).collect();
    if unsaturated.is_empty() {
        return Some(base);
    }
    if unsaturated.len() != 1 {
        return None;
    }
    let stem = &base[..base.len() - 3];
    let ending = match unsaturated[0].1 {
        2 => "ene",
        3 => "yne",
        _ => return None,
    };
    Some(if length == 2 {
        format!("{stem}{ending}")
    } else {
        format!("{stem}-{}-{ending}", unsaturated[0].0 + 1)
    })
}
// Provenance: fn locanted_ligand_prefix | adapted
// Source: openclatura | src/openclatura/special_cases.py::_locanted_ligand_prefix
fn locanted_ligand_prefix(items: &[(usize, String)]) -> String {
    let mut counts: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (loc, name) in items {
        counts.entry(name.clone()).or_default().push(*loc);
    }
    let mut groups: Vec<_> = counts.into_iter().collect();
    groups.sort_by_key(|(name, _)| substituent_sort_key(name));
    groups
        .into_iter()
        .map(|(name, mut locs)| {
            locs.sort();
            format!(
                "{}-{}",
                locs.iter()
                    .map(usize::to_string)
                    .collect::<Vec<_>>()
                    .join(","),
                format_multiplier(&name, locs.len(), false)
            )
        })
        .collect::<Vec<_>>()
        .join("-")
}
// Provenance: fn homonuclear_chain_parent_name | adapted
// Source: openclatura | src/openclatura/special_cases.py::homonuclear_chain_parent_result
fn homonuclear_chain_parent_name(
    mol: &Molecule,
    atoms: &AtomSet,
) -> Result<Option<String>, String> {
    let Some((symbol, chain)) = homonuclear_backbone(mol, atoms) else {
        return Ok(None);
    };
    let chain_set: AtomSet = chain.iter().copied().collect();
    let oxide = chain_oxide_oxygens(mol, &chain, atoms);
    if chain.iter().any(|i| mol.atoms[i].charge != 0) && oxide.is_none() {
        return Ok(None);
    }
    let oxide = oxide.unwrap_or_default();
    let oxide_atoms: AtomSet = oxide.values().copied().collect();
    let ligand_atoms: AtomSet = atoms.difference(&oxide_atoms).copied().collect();
    if perceive_groups(mol).iter().any(|g| {
        g.is_principal_candidate
            && g.atoms_involved
                .intersection(&ligand_atoms)
                .any(|i| !chain_set.contains(i))
    }) {
        return Ok(None);
    }
    // P-72.2.2.2.4: an anionic oxygen on a chalcogen chain is named from the
    // anionic parent hydride, as in tert-butyltrioxidanide and
    // (ethylsulfanyl)oxidanide (P-72.2.2.2.2). `oxido` on the chain would be
    // read as a chalcogen oxide.
    if matches!(symbol.as_str(), "S" | "Se" | "Te")
        && oxide.is_empty()
        && let Some(&o) = ligand_atoms.iter().find(|&&i| {
            mol.atoms[&i].symbol == "O" && mol.atoms[&i].charge == -1 && mol.degree(i) == 1
        })
        && let Some(&center) = mol.get_neighbors(o).first()
        && chain_set.contains(&center)
    {
        // The prefix carries no other charge and no backbone stereo.
        if atoms.iter().any(|&i| i != o && mol.atoms[&i].charge != 0) {
            return Err("openclatura: second charge beside a chalcogen-chain oxidanide".into());
        }
        if chain_set.iter().any(|i| mol.atoms[i].has_stereo()) {
            return Err(SHORTCUT_STEREO_ERROR.into());
        }
        let prefix = branch_name(mol, atoms, center, &set([o]), o)?;
        return Ok(Some(format!(
            "({})oxidanide",
            strip_outer_parentheses(&prefix)
        )));
    }
    let mut ligands = Vec::new();
    let mut covered = chain_set.clone();
    for branch in connected_branches(mol, &ligand_atoms.difference(&chain_set).copied().collect()) {
        let mut links = Vec::new();
        for &i in &branch {
            for &n in mol.get_neighbors(i) {
                if chain_set.contains(&n) {
                    links.push((n, i));
                }
            }
        }
        links.sort();
        if links.len() != 1 {
            return Ok(None);
        }
        let (attachment, root) = links[0];
        let bond = mol.bond_order(attachment, root);
        if !matches!(bond, 1 | 2) || (bond == 2 && chain.len() <= 2) {
            return Ok(None);
        }
        let mut name = if bond == 2 {
            String::new()
        } else {
            terminal_ligand_name(mol, root)
        };
        if name.is_empty() {
            name = branch_name(mol, &ligand_atoms, root, &chain_set, attachment)?;
        }
        if name.is_empty() || bond == 2 && !strip_outer_parentheses(&name).ends_with("ylidene") {
            return Ok(None);
        }
        ligands.push((attachment, name));
        covered.extend(branch);
    }
    if covered != ligand_atoms {
        return Ok(None);
    }
    let forward: Vec<_> = chain
        .windows(2)
        .map(|p| mol.bond_order(p[0], p[1]))
        .collect();
    struct Choice {
        key: (usize, Vec<usize>, Vec<usize>, Vec<usize>),
        parent: String,
        placed: Vec<(usize, String)>,
        oriented: Vec<usize>,
    }
    let mut best: Option<Choice> = None;
    for (oriented, orders) in [
        (chain.clone(), forward.clone()),
        (
            chain.iter().copied().rev().collect(),
            forward.iter().copied().rev().collect(),
        ),
    ] {
        let Some(parent) = same_element_parent_name(&symbol, oriented.len(), &orders) else {
            continue;
        };
        let locants: BTreeMap<_, _> = oriented
            .iter()
            .enumerate()
            .map(|(offset, atom)| (*atom, offset + 1))
            .collect();
        let mut placed: Vec<_> = ligands
            .iter()
            .map(|(i, name)| (locants[i], name.clone()))
            .collect();
        placed.sort_by_key(|(loc, name)| (substituent_sort_key(name), *loc));
        let mut oxide_locs: Vec<_> = oxide.keys().map(|i| locants[i]).collect();
        oxide_locs.sort();
        let mut sorted_locs: Vec<_> = placed.iter().map(|(l, _)| *l).collect();
        sorted_locs.sort();
        let key = (
            orders.iter().position(|o| *o > 1).unwrap_or(orders.len()),
            oxide_locs,
            sorted_locs,
            placed.iter().map(|(l, _)| *l).collect(),
        );
        if best.as_ref().is_none_or(|best| key < best.key) {
            best = Some(Choice {
                key,
                parent,
                placed,
                oriented,
            });
        }
    }
    let Some(Choice {
        parent,
        placed,
        oriented,
        ..
    }) = best
    else {
        return Ok(None);
    };
    // Upstream divergence: homonuclear chain assembly renders bond descriptors but no backbone atom stereo.
    if chain_set
        .iter()
        .chain(&oxide_atoms)
        .any(|&idx| mol.atoms[&idx].has_stereo())
    {
        return Err(SHORTCUT_STEREO_ERROR.into());
    }
    let prefix = locanted_ligand_prefix(&placed);
    let mut name = format!("{prefix}{parent}");
    if !oxide.is_empty() {
        let locants: BTreeMap<_, _> = oriented
            .iter()
            .enumerate()
            .map(|(offset, atom)| (*atom, offset + 1))
            .collect();
        let mut locs: Vec<_> = oxide.keys().map(|i| locants[i]).collect();
        locs.sort();
        let mult = if locs.len() > 1 {
            multipliers::basic(locs.len())
        } else {
            String::new()
        };
        name.push_str(&format!(
            " {}-{mult}oxide",
            locs.iter()
                .map(usize::to_string)
                .collect::<Vec<_>>()
                .join(",")
        ));
    }
    let stereo: Vec<_> = oriented
        .windows(2)
        .enumerate()
        .filter_map(|(i, p)| {
            let b = mol.get_bond(p[0], p[1])?;
            let s = b.stereo.as_deref()?;
            (b.order == 2 && matches!(s, "E" | "Z")).then(|| format!("{}{s}", i + 1))
        })
        .collect();
    if !stereo.is_empty() && !descriptors_omitted() {
        name = format!("({})-{name}", stereo.join(","));
    }
    if symbol == "S" && chain.len() <= 2 && !placed.is_empty() {
        // R-SSH with a carbon R takes the dithioperoxol suffix (P-63.4.2.1), unless a
        // phosphorus acid, which the generic path cites as a prefix, outranks it.
        if ligands.len() == 1
            && oxide.is_empty()
            && mol
                .get_neighbors(ligands[0].0)
                .iter()
                .any(|n| !chain_set.contains(n) && mol.atoms[n].is_carbon())
            && !atoms.iter().any(|i| mol.atoms[i].symbol == "P")
        {
            return Ok(None);
        }
        mol.mark_name(NameStatus::NotPreferred(NameReason::HydrideParent));
    }
    if symbol == "O"
        && chain.len() == 2
        && chain.iter().any(|&i| {
            mol.get_neighbors(i)
                .iter()
                .any(|n| !chain_set.contains(n) && mol.atoms[n].is_carbon())
        })
    {
        // P-63.3.1 prefers a peroxy prefix for organyl-substituted dioxidanes.
        mol.mark_name(NameStatus::NotPreferred(NameReason::HydrideParent));
    }
    Ok(Some(name))
}

/// Hydrazine as the parent hydride of an acid, ester, acid halide, amide or
/// hydrazide carbon on one of its nitrogens: N-N-C(=X)-Y is named with a
/// `carbo...` suffix on hydrazine, numbered 1 and 2, as in
/// 2-methylhydrazine-1-carboxylic acid, N,1-dimethylhydrazine-1-carboxamide
/// and hydrazinecarbohydrazide (P-65.1.2.2.2, P-68.3.1.2.1, P-68.3.1.2.4,
/// P-68.3.1.2.6). Hydrazine is the parent only while it carries at least as
/// many principal characteristic groups as any other parent (P-44.1.1) and
/// no ring nitrogen or second hydrazine competes (P-44.1.2, P-44.1.2.2);
/// anything this cannot decide is left to the general namer.
// Provenance: fn hydrazine_carbonyl_parent_name | original
// Source: Nazuchem | crates/nazuchem-s2n/src/special_cases.rs::hydrazine_carbonyl_parent_name
fn hydrazine_carbonyl_parent_name(
    mol: &Molecule,
    atoms: &AtomSet,
) -> Result<Option<String>, String> {
    // A suffix carbon as seen from one hydrazine nitrogen. `code` is the
    // class seniority used by the functional-group rules (20 acid, 22 peroxy
    // acid, 40 ester, 50-53 acid halide, 60 amide, 68 hydrazide); `plain`
    // marks the unreplaced suffix, the only one compared with other groups.
    struct Suffix {
        carbon: usize,
        nitrogen: usize,
        code: i64,
        plain: bool,
        text: String,
        // Italic element locant of a chalcogen acid or ester (P-65.1.5.1).
        site: &'static str,
        y: usize,
    }
    if mol.is_isotope_parent() {
        return Ok(None);
    }
    let cyclic = get_cyclic_atoms(mol, &AtomSet::new());
    let symbol = |i: usize| mol.atoms[&i].symbol.as_str();
    let heavy = |i: usize| -> Vec<usize> {
        mol.get_neighbors(i)
            .iter()
            .copied()
            .filter(|&x| symbol(x) != "H")
            .collect()
    };
    let chalcogen = |i: usize| matches!(symbol(i), "O" | "S" | "Se" | "Te");
    let open = |i: usize| mol.atoms[&i].charge == 0 && !cyclic.contains(&i);
    // The anion of an acid on hydrazine: N-N-C(=X)-Y(-), named with its
    // suffix (P-72.2.2.2.1).
    let anion = |y: usize| {
        mol.atoms[&y].charge == -1
            && chalcogen(y)
            && matches!(heavy(y)[..], [c] if symbol(c) == "C"
                && heavy(c).iter().any(|&x| chalcogen(x) && mol.bond_order(c, x) == 2)
                && heavy(c).iter().any(|&n| symbol(n) == "N"
                    && heavy(n).iter().any(|&m| symbol(m) == "N")))
    };
    for &i in atoms {
        let atom = &mol.atoms[&i];
        // Other ions, other elements and sulfur oxoacids have their own
        // classes.
        if !matches!(
            symbol(i),
            "C" | "H" | "N" | "O" | "S" | "Se" | "Te" | "F" | "Cl" | "Br" | "I"
        ) || atom.charge != 0
            && !anion(i)
            && !mol
                .get_neighbors(i)
                .iter()
                .any(|&x| mol.atoms[&x].charge == -atom.charge)
            || chalcogen(i)
                && heavy(i).len() > 2
                && heavy(i)
                    .iter()
                    .any(|&x| mol.bond_order(i, x) == 1 && !matches!(symbol(x), "C" | "N"))
        {
            return Ok(None);
        }
    }
    // The other nitrogen of an acyclic N-N single bond whose nitrogens carry
    // only carbon.
    let partner = |n: usize| -> Option<usize> {
        let carbon_only = |a: usize, b: usize| {
            symbol(a) == "N" && open(a) && heavy(a).iter().all(|&x| x == b || symbol(x) == "C")
        };
        let ns: Vec<_> = heavy(n).into_iter().filter(|&x| symbol(x) == "N").collect();
        match ns[..] {
            [m] if mol.bond_order(n, m) == 1 && carbon_only(n, m) && carbon_only(m, n) => Some(m),
            _ => None,
        }
    };
    let suffix_at = |carbon: usize, nitrogen: usize| -> Option<Suffix> {
        let others: Vec<_> = heavy(carbon)
            .into_iter()
            .filter(|&x| x != nitrogen)
            .collect();
        let [p, q] = others[..] else { return None };
        let (x, y) = if mol.bond_order(carbon, p) == 2 {
            (p, q)
        } else {
            (q, p)
        };
        if symbol(carbon) != "C"
            || !open(carbon)
            || mol.bond_order(carbon, nitrogen) != 1
            || mol.bond_order(carbon, x) != 2
            || mol.bond_order(carbon, y) != 1
            || mol.atoms[&x].charge != 0
            || heavy(x).len() != 1
            || !(open(y) || anion(y))
        {
            return None;
        }
        let infix = match symbol(x) {
            "O" => "",
            "S" => "thio",
            "Se" => "seleno",
            "Te" => "telluro",
            _ => return None,
        };
        let plain = infix.is_empty();
        let beyond: Vec<_> = heavy(y).into_iter().filter(|&z| z != carbon).collect();
        let (code, plain, text, site) = match symbol(y) {
            "F" | "Cl" | "Br" | "I" if beyond.is_empty() => {
                let (code, word) = match symbol(y) {
                    "F" => (50, "fluoride"),
                    "Cl" => (51, "chloride"),
                    "Br" => (52, "bromide"),
                    _ => (53, "iodide"),
                };
                let acyl = if plain { "carbonyl" } else { "carbo" };
                let yl = if plain { "" } else { "yl" };
                (code, plain, format!("{acyl}{infix}{yl} {word}"), "")
            }
            "N" if beyond.iter().all(|&z| mol.bond_order(y, z) == 1) => {
                match beyond
                    .iter()
                    .filter(|&&z| symbol(z) == "N")
                    .collect::<Vec<_>>()[..]
                {
                    [] if plain => (60, true, "carboxamide".into(), ""),
                    [] => (60, false, format!("carbo{infix}amide"), ""),
                    [_] if partner(y).is_some() => {
                        (68, plain, format!("carbo{infix}hydrazide"), "")
                    }
                    _ => return None,
                }
            }
            "O" | "S" | "Se" | "Te" => {
                let ending = match beyond[..] {
                    [] if mol.atoms[&y].total_h_count == 1 => "ic acid",
                    [] if anion(y) => "ate",
                    // An acyl carbon here would make an anhydride.
                    [r] if symbol(r) == "C"
                        && heavy(r)
                            .iter()
                            .all(|&z| symbol(z) == "C" || mol.bond_order(r, z) == 1) =>
                    {
                        "ate"
                    }
                    // P-65.1.4.1: -CO-OOH, carboperoxoic acid.
                    [o] if plain
                        && symbol(y) == "O"
                        && symbol(o) == "O"
                        && mol.atoms[&o].charge == 0
                        && heavy(o).len() == 1
                        && mol.atoms[&o].total_h_count == 1 =>
                    {
                        return Some(Suffix {
                            carbon,
                            nitrogen,
                            code: 22,
                            plain: false,
                            text: "carboperoxoic acid".into(),
                            site: "",
                            y,
                        });
                    }
                    _ => return None,
                };
                let code = if beyond.is_empty() { 20 } else { 40 };
                let replaced = match symbol(y) {
                    "O" => "",
                    "S" => "thio",
                    "Se" => "seleno",
                    _ => "telluro",
                };
                // P-65.1.5.1: one replaced chalcogen needs the O/S/Se/Te
                // locant of the acid or ester site; mixed pairs are left out.
                match (infix, replaced) {
                    ("", "") if ending == "ate" => (code, true, "carboxylate".into(), ""),
                    ("", "") => (code, true, "carboxylic acid".into(), ""),
                    (a, b) if a == b => (code, false, format!("carbodi{a}{ending}"), ""),
                    // The charge of a mixed anion is not placed.
                    _ if anion(y) => return None,
                    (a, "") | ("", a) => {
                        let site = match symbol(y) {
                            "O" => "O-",
                            "S" => "S-",
                            "Se" => "Se-",
                            _ => "Te-",
                        };
                        (code, false, format!("carbo{a}{ending}"), site)
                    }
                    _ => return None,
                }
            }
            _ => return None,
        };
        Some(Suffix {
            carbon,
            nitrogen,
            code,
            plain,
            text,
            site,
            y,
        })
    };
    // A carbon substituent the prefix namer can cite: alkyl, aryl or acyl.
    let citable = |root: usize, from: usize| {
        let multiple: Vec<_> = heavy(root)
            .into_iter()
            .filter(|&x| x != from && symbol(x) != "C" && mol.bond_order(root, x) > 1)
            .collect();
        symbol(root) == "C"
            && match multiple[..] {
                [] => true,
                [x] => {
                    mol.bond_order(root, x) == 2
                        && chalcogen(x)
                        && heavy(root)
                            .iter()
                            .all(|&z| z == from || z == x || symbol(z) == "C")
                }
                _ => false,
            }
    };
    let mut suffixes = Vec::new();
    for &n in atoms {
        if partner(n).is_some() {
            suffixes.extend(heavy(n).into_iter().filter_map(|c| suffix_at(c, n)));
        }
    }
    // Every acid-type carbon on an N-N nitrogen must be one of those suffixes.
    for &c in atoms {
        if symbol(c) == "C"
            && !suffixes.iter().any(|s: &Suffix| s.carbon == c)
            && open(c)
            && heavy(c).iter().any(|&n| {
                symbol(n) == "N" && heavy(n).iter().any(|&m| symbol(m) == "N") && !citable(c, n)
            })
        {
            return Ok(None);
        }
    }
    let Some(code) = suffixes.iter().map(|s| s.code).min() else {
        return Ok(None);
    };
    // P-41 classes: acids, anhydrides, esters, acid halides, amides,
    // hydrazides, then everything junior.
    let class = |code: i64| match code {
        ..20 => 0,
        20..30 => 1,
        30..40 => 2,
        40..50 => 3,
        50..60 => 4,
        60..68 => 5,
        68 => 6,
        _ => 7,
    };
    let principal: Vec<_> = suffixes.iter().filter(|s| s.code == code).collect();
    let first = principal[0];
    if suffixes.iter().any(|s| {
        class(s.code) == class(code) && (s.code, &s.text, s.site) != (code, &first.text, first.site)
    }) {
        return Ok(None);
    }
    // Groups elsewhere: a senior class wins outright and the same suffix
    // competes by count (P-44.1.1). Within a class the unreplaced carbon
    // suffix leads (P-41, classes 7a and 9 to 12); other sub-class orders are
    // not decided here.
    let mut outside = 0;
    for group in perceive_groups(mol) {
        let own = group.atom_ids();
        if !atoms.contains(&group.attachment_carbon)
            || suffixes.iter().any(|s| own.contains(&s.carbon))
        {
            continue;
        }
        let Some(seniority) = group
            .resolved_rule
            .as_ref()
            .and_then(|r| r.seniority.clone())
        else {
            continue;
        };
        if class(seniority[0]) < class(code) {
            return Ok(None);
        }
        if class(seniority[0]) == class(code) {
            if first.plain && seniority == [code] {
                outside += 1;
            } else if !first.plain || seniority[0] <= code {
                return Ok(None);
            }
        }
    }
    let mut units: Vec<[usize; 2]> = principal
        .iter()
        .filter_map(|s| Some([s.nitrogen, partner(s.nitrogen)?]))
        .map(|[a, b]| [a.min(b), a.max(b)])
        .collect();
    units.sort();
    units.dedup();
    let on_unit = |unit: &[usize; 2]| {
        principal
            .iter()
            .filter(|s| unit.contains(&s.nitrogen))
            .count()
    };
    let count = units.iter().map(on_unit).max().unwrap_or(0);
    units.retain(|unit| on_unit(unit) == count);
    // Two hydrazines tie only as the two ends of one carbohydrazide.
    let shared = units.len() == 2
        && count == 1
        && principal.len() == 2
        && principal[0].carbon == principal[1].carbon;
    if units.len() > 1 && !shared
        || outside > count
        || outside == count
            && atoms
                .iter()
                .any(|&i| symbol(i) == "N" && cyclic.contains(&i))
        || count > 2
        || count == 2 && !(first.plain || first.code == 60)
    {
        return Ok(None);
    }
    let rank = |locant: &str| ["N", "N'", "1", "2"].iter().position(|l| *l == locant);
    type Key = (usize, Vec<usize>, Vec<Option<usize>>, Vec<Option<usize>>);
    let mut best: Option<(Key, String)> = None;
    for unit in &units {
        let own: Vec<_> = principal
            .iter()
            .filter(|s| unit.contains(&s.nitrogen))
            .collect();
        if own.len() == 2 && own[0].nitrogen == own[1].nitrogen {
            return Ok(None);
        }
        for [one, two] in [[unit[0], unit[1]], [unit[1], unit[0]]] {
            let mut core: AtomSet = unit.iter().copied().collect();
            // Substituent prefixes with their locants, and ester words.
            let mut cited: Vec<(&str, String)> = Vec::new();
            let mut esters = Vec::new();
            let mut sites: Vec<(&str, usize, usize)> = vec![("1", one, two), ("2", two, one)];
            for s in &own {
                core.extend([s.carbon, s.y]);
                core.extend(heavy(s.carbon));
                match s.code {
                    40 => {
                        let root = heavy(s.y).into_iter().find(|&r| r != s.carbon);
                        let root = root.ok_or("openclatura: ester without alkyl group")?;
                        esters.push(ester_modifier_name(
                            mol,
                            atoms,
                            root,
                            s.y,
                            &set([s.carbon, s.y]),
                        )?);
                    }
                    60 => sites.push(("N", s.y, s.carbon)),
                    68 => {
                        let far = partner(s.y).ok_or("openclatura: hydrazide without N'")?;
                        core.insert(far);
                        sites.extend([("N", s.y, s.carbon), ("N'", far, s.y)]);
                    }
                    _ => {}
                }
            }
            let mut bonds = 0;
            for &(locant, n, back) in &sites {
                for root in heavy(n) {
                    if root == back || core.contains(&root) {
                        continue;
                    }
                    let order = mol.bond_order(n, root);
                    let numeric = rank(locant) > Some(1);
                    // No Blue Book example places an E/Z descriptor here.
                    let ylidene = order == 2
                        && heavy(root).iter().all(|&z| z == n || symbol(z) == "C")
                        && mol.get_bond(n, root).is_some_and(|b| b.stereo.is_none());
                    // A junior acid-type group on hydrazine is a prefix.
                    let junior = suffix_at(root, n).filter(|_| numeric && order == 1);
                    if !(ylidene || junior.is_some() || order == 1 && citable(root, n)) {
                        return Ok(None);
                    }
                    let name = match &junior {
                        // P-65.2.1.5: Cl-CO- is carbonochloridoyl.
                        Some(s) if matches!(s.code, 50..=53) && s.plain => format!(
                            "carbono{}oyl",
                            s.text["carbonyl ".len()..].trim_end_matches('e')
                        ),
                        // P-66.3.2.1: H2N-NH-CO- is hydrazinecarbonyl.
                        Some(s)
                            if s.code == 68
                                && s.plain
                                && heavy(s.y).len() == 2
                                && partner(s.y).is_some_and(|far| heavy(far).len() == 1) =>
                        {
                            "(hydrazinecarbonyl)".into()
                        }
                        Some(s) if s.code != 60 && !(s.code == 40 && s.plain) => {
                            return Ok(None);
                        }
                        _ => branch_name(mol, atoms, root, &set([n]), n)?,
                    };
                    if name.is_empty()
                        || ylidene != strip_outer_parentheses(&name).ends_with("ylidene")
                    {
                        return Ok(None);
                    }
                    // P-29.6.2.1: unsubstituted benzylidene is preferred.
                    let name = if name == "(phenylmethylidene)" {
                        "benzylidene".into()
                    } else {
                        name
                    };
                    bonds += order;
                    cited.push((locant, name));
                }
            }
            if require_stereo_free_core(mol, &core).is_err() || esters.iter().any(String::is_empty)
            {
                return Ok(None);
            }
            let mut names: Vec<_> = cited.iter().map(|(_, name)| name.clone()).collect();
            names.sort_by_key(|name| substituent_sort_key(name));
            names.dedup();
            cited
                .sort_by_key(|(locant, name)| (names.iter().position(|n| n == name), rank(locant)));
            let numbered = cited.iter().any(|(locant, _)| rank(locant) > Some(1));
            // Not spelled here: N^1/N^2 locants and unlike ester words on two
            // suffixes, and two suffixes with every position substituted
            // alike (P-14.3.4.5).
            if own.len() == 2
                && (cited.iter().any(|(locant, _)| rank(locant) < Some(2))
                    || esters.iter().any(|e| *e != esters[0])
                    || names.len() == 1 && bonds == 2)
            {
                return Ok(None);
            }
            // P-14.3.4.5: no locants when every position carries the same
            // substituent.
            let full = own.len() == 1
                && names.len() == 1
                && bonds == 2 * sites.len() - 1 - usize::from(code == 68);
            let prefixes = names
                .iter()
                .map(|name| {
                    let locants: Vec<_> = cited
                        .iter()
                        .filter(|(_, n)| n == name)
                        .map(|(locant, _)| *locant)
                        .collect();
                    let text = format_multiplier(name, locants.len(), false);
                    if full {
                        text
                    } else {
                        format!("{}-{text}", locants.join(","))
                    }
                })
                .collect::<Vec<_>>()
                .join("-");
            let mut locants: Vec<_> = own
                .iter()
                .map(|s| if s.nitrogen == one { "1" } else { "2" })
                .collect();
            locants.sort();
            // P-14.3.4.2 (b): monosubstituted hydrazine needs no locant, as
            // in N-phenylhydrazinecarboxamide.
            let position = if own.len() == 1 && (full || !numbered) {
                String::new()
            } else {
                format!("-{}-", locants.join(","))
            };
            let suffix = if own.len() == 1 {
                first.text.clone()
            } else {
                // hydrazine-1,2-dicarbonyl dichloride
                first
                    .text
                    .split(' ')
                    .map(|word| {
                        if word == "acid" {
                            word.into()
                        } else {
                            format!("di{word}")
                        }
                    })
                    .collect::<Vec<_>>()
                    .join(" ")
            };
            let ester = if esters.is_empty() {
                String::new()
            } else if first.site.is_empty() {
                format!("{} ", modifier_phrase(&esters, false))
            } else {
                format!("{}{} ", first.site, format_multiplier(&esters[0], 1, false))
            };
            let suffix = suffix.replace(" acid", &format!(" {}acid", first.site));
            let key: Key = (
                usize::MAX - cited.iter().filter(|(l, _)| rank(l) > Some(1)).count(),
                locants.iter().map(|l| usize::from(*l == "2")).collect(),
                {
                    let mut all: Vec<_> = cited.iter().map(|(l, _)| rank(l)).collect();
                    all.sort();
                    all
                },
                cited.iter().map(|(l, _)| rank(l)).collect(),
            );
            let name = format!("{ester}{prefixes}hydrazine{position}{suffix}");
            if best.as_ref().is_none_or(|(k, _)| key < *k) {
                best = Some((key, name));
            }
        }
    }
    Ok(best.map(|(_, name)| name))
}

/// Cyanohydrazide, NC-NH-NH2, the hydrazide of cyanic acid (P-66.3.1.2.1,
/// P-66.3.5.1). Its carbon substituents take the locants N, the nitrogen on
/// the cyano carbon, and N' (P-66.3.3.1). Left to the general namer: a
/// second cyanohydrazide, and a group that is senior to a hydrazide or ties
/// with it (P-41).
// Provenance: fn cyanohydrazide_name | original
// Source: Nazuchem | crates/nazuchem-s2n/src/special_cases.rs::cyanohydrazide_name
fn cyanohydrazide_name(mol: &Molecule, atoms: &AtomSet) -> Result<Option<String>, String> {
    if mol.is_isotope_parent() || atoms.iter().any(|i| mol.atoms[i].charge != 0) {
        return Ok(None);
    }
    let symbol = |i: usize| mol.atoms[&i].symbol.as_str();
    let cyclic = get_cyclic_atoms(mol, &AtomSet::new());
    let heavy = |i: usize| -> Vec<usize> {
        mol.get_neighbors(i)
            .iter()
            .copied()
            .filter(|&x| symbol(x) != "H")
            .collect()
    };
    // N#C-N-N with an acyclic two-nitrogen chain.
    let mut cores = Vec::new();
    for &c in atoms {
        let [p, q] = heavy(c)[..] else { continue };
        for (nitrile, n) in [(p, q), (q, p)] {
            let chain: Vec<_> = heavy(n).into_iter().filter(|&x| symbol(x) == "N").collect();
            if let [far] = chain[..]
                && symbol(c) == "C"
                && symbol(nitrile) == "N"
                && mol.bond_order(c, nitrile) == 3
                && symbol(n) == "N"
                && mol.bond_order(c, n) == 1
                && mol.bond_order(n, far) == 1
                && !cyclic.contains(&n)
                && !cyclic.contains(&far)
                && heavy(far).iter().all(|&x| x == n || symbol(x) != "N")
            {
                cores.push([c, nitrile, n, far]);
            }
        }
    }
    let [[c, nitrile, n, far]] = cores[..] else {
        return Ok(None);
    };
    for group in perceive_groups(mol) {
        if atoms.contains(&group.attachment_carbon)
            && !group.atom_ids().contains(&c)
            && group
                .resolved_rule
                .as_ref()
                .and_then(|r| r.seniority.as_ref())
                .is_some_and(|s| s[0] <= 68)
        {
            return Ok(None);
        }
    }
    let core = set([c, nitrile, n, far]);
    // Not every senior class is perceived as a group: carbonic, phosphorus
    // and nitric acids and their esters are not. They all have a carbon
    // with two heteroatoms or a bond between two heteroatoms.
    let hetero =
        |i: usize| !core.contains(&i) && !matches!(symbol(i), "C" | "H" | "F" | "Cl" | "Br" | "I");
    if atoms.iter().any(|&i| {
        !core.contains(&i)
            && heavy(i).into_iter().filter(|&x| hetero(x)).count() > usize::from(!hetero(i))
    }) {
        return Ok(None);
    }
    let mut cited: Vec<(String, &str)> = Vec::new();
    for (locant, nitrogen) in [("N", n), ("N'", far)] {
        for root in heavy(nitrogen) {
            if core.contains(&root) {
                continue;
            }
            // No Blue Book example places an E/Z descriptor here.
            let order = mol.bond_order(nitrogen, root);
            if symbol(root) != "C"
                || order > 2
                || mol
                    .get_bond(nitrogen, root)
                    .is_some_and(|b| b.stereo.is_some())
            {
                return Ok(None);
            }
            let name = branch_name(mol, atoms, root, &set([nitrogen]), nitrogen)?;
            // A name that does not carry the double bond, such as
            // methanethioyl for =C=S, would lose it.
            if name.is_empty()
                || (order == 2) != strip_outer_parentheses(&name).ends_with("ylidene")
            {
                return Ok(None);
            }
            // P-29.6.2.1: unsubstituted benzylidene is preferred.
            let name = if name == "(phenylmethylidene)" {
                "benzylidene".into()
            } else {
                name
            };
            cited.push((name, locant));
        }
    }
    if require_stereo_free_core(mol, &core).is_err() {
        return Ok(None);
    }
    let mut names: Vec<_> = cited.iter().map(|(name, _)| name.clone()).collect();
    names.sort_by_cached_key(|name| (substituent_sort_key(name), name.clone()));
    names.dedup();
    // P-14.3.4.5: no locants when all three positions carry the same
    // substituent, as in tetrafluorourea.
    let full = cited.len() == 3 && names.len() == 1;
    let prefixes = names
        .iter()
        .map(|name| {
            let locants: Vec<_> = cited
                .iter()
                .filter(|(n, _)| n == name)
                .map(|(_, locant)| *locant)
                .collect();
            let text = format_multiplier(name, locants.len(), false);
            if full {
                text
            } else {
                format!("{}-{text}", locants.join(","))
            }
        })
        .collect::<Vec<_>>()
        .join("-");
    Ok(Some(format!("{prefixes}cyanohydrazide")))
}

// Provenance: fn nitrogen_ligands | adapted
// Source: openclatura | src/openclatura/special_cases.py::_urea_n_ligand_names
fn nitrogen_ligands(
    mol: &Molecule,
    atoms: &AtomSet,
    nitrogen: usize,
    core: &AtomSet,
    carbon_only: bool,
) -> Result<Option<Vec<(AtomSet, String)>>, String> {
    let mut ligands = Vec::new();
    for &root in mol.get_neighbors(nitrogen) {
        if core.contains(&root) || !atoms.contains(&root) || mol.atoms[&root].symbol == "H" {
            continue;
        }
        if carbon_only && !mol.atoms[&root].is_carbon() {
            return Ok(None);
        }
        let branch = component_atoms_until_blocked(mol, atoms, root, core);
        if branch.is_empty() {
            return Ok(None);
        }
        if carbon_only
            && (branch.iter().any(|i| !mol.atoms[i].is_carbon())
                || branch.iter().any(|i| {
                    mol.get_neighbors(*i)
                        .iter()
                        .any(|n| branch.contains(n) && mol.bond_order(*i, *n) != 1)
                }))
        {
            return Ok(None);
        }
        let mut name = strip_outer_parentheses(&branch_name(mol, atoms, root, core, nitrogen)?);
        if name.is_empty() {
            name = alkyl_ligand_name(mol, atoms, root, nitrogen)?;
        }
        if name.is_empty() {
            return Ok(None);
        }
        ligands.push((branch, name));
    }
    ligands.sort_by(|a, b| a.1.cmp(&b.1));
    Ok(Some(ligands))
}
// Provenance: fn hydroxyurea_parent_name | adapted
// Source: openclatura | src/openclatura/special_cases.py::hydroxyurea_parent_result
fn hydroxyurea_parent_name(mol: &Molecule, atoms: &AtomSet) -> Result<Option<String>, String> {
    let centers: Vec<_> = atoms
        .iter()
        .copied()
        .filter(|i| {
            mol.atoms[i].is_carbon()
                && mol
                    .get_neighbors(*i)
                    .iter()
                    .filter(|n| {
                        atoms.contains(n)
                            && mol.atoms[n].symbol == "O"
                            && mol.bond_order(*i, **n) == 2
                    })
                    .count()
                    == 1
                && mol
                    .get_neighbors(*i)
                    .iter()
                    .filter(|n| {
                        atoms.contains(n)
                            && mol.atoms[n].symbol == "N"
                            && mol.bond_order(*i, **n) == 1
                    })
                    .count()
                    == 2
        })
        .collect();
    if centers.len() != 1 {
        return Ok(None);
    }
    let c = centers[0];
    let Some(oxygen) = mol
        .get_neighbors(c)
        .iter()
        .copied()
        .find(|n| atoms.contains(n) && mol.atoms[n].symbol == "O")
    else {
        return Ok(None);
    };
    let nitrogens: Vec<_> = mol
        .get_neighbors(c)
        .iter()
        .copied()
        .filter(|n| atoms.contains(n) && mol.atoms[n].symbol == "N" && mol.bond_order(c, *n) == 1)
        .collect();
    if nitrogens
        .iter()
        .any(|n| mol.get_bond(c, *n).is_some_and(|b| b.in_small_ring))
    {
        return Ok(None);
    }
    let mut hydroxy = Vec::new();
    for &n in &nitrogens {
        let oxygens: Vec<_> = mol
            .get_neighbors(n)
            .iter()
            .copied()
            .filter(|i| {
                atoms.contains(i)
                    && *i != c
                    && mol.atoms[i].symbol == "O"
                    && mol.bond_order(n, *i) == 1
            })
            .collect();
        if !oxygens.is_empty() {
            if oxygens.len() != 1
                || mol
                    .get_neighbors(oxygens[0])
                    .iter()
                    .any(|i| *i != n && atoms.contains(i))
            {
                return Ok(None);
            }
            hydroxy.push((n, oxygens[0]));
        }
    }
    if hydroxy.len() != 1 {
        return Ok(None);
    }
    let (n_prime, oh) = hydroxy[0];
    let Some(&n) = nitrogens.iter().find(|i| **i != n_prime) else {
        return Ok(None);
    };
    let core = set([c, oxygen, n, n_prime, oh]);
    require_stereo_free_core(mol, &core)?;
    let Some(first) = nitrogen_ligands(mol, atoms, n, &core, true)? else {
        return Ok(None);
    };
    let Some(second) = nitrogen_ligands(mol, atoms, n_prime, &core, true)? else {
        return Ok(None);
    };
    let mut represented = core;
    for (branch, _) in first.iter().chain(&second) {
        represented.extend(branch);
    }
    if represented != *atoms {
        return Ok(None);
    }
    let mut prefixes: Vec<_> = first
        .iter()
        .map(|(_, name)| format!("N-{}", format_multiplier(name, 1, false)))
        .collect();
    prefixes.extend(
        second
            .iter()
            .map(|(_, name)| format!("N'-{}", format_multiplier(name, 1, false))),
    );
    prefixes.push("N'-hydroxy".into());
    Ok(Some(format!("{}urea", prefixes.join("-"))))
}
// Provenance: fn sulfamic_acid_name | adapted
// Source: openclatura | src/openclatura/special_cases.py::sulfamic_acid_result
fn sulfamic_acid_name(mol: &Molecule, atoms: &AtomSet) -> Result<Option<String>, String> {
    let cyclic = get_cyclic_atoms(mol, &AtomSet::new());
    let mut centers = Vec::new();
    for &c in atoms {
        if mol.atoms[&c].symbol != "S" || mol.atoms[&c].charge != 0 {
            continue;
        }
        let neighbors: Vec<_> = mol
            .get_neighbors(c)
            .iter()
            .copied()
            .filter(|n| atoms.contains(n))
            .collect();
        let oxygens: Vec<_> = neighbors
            .iter()
            .copied()
            .filter(|n| mol.atoms[n].symbol == "O")
            .collect();
        let nitrogens: Vec<_> = neighbors
            .iter()
            .copied()
            .filter(|n| mol.atoms[n].symbol == "N")
            .collect();
        if neighbors.len() != 4
            || oxygens.len() != 3
            || nitrogens.len() != 1
            || cyclic.contains(&nitrogens[0])
        {
            continue;
        }
        let double: Vec<_> = oxygens
            .iter()
            .copied()
            .filter(|o| mol.bond_order(c, *o) == 2)
            .collect();
        let single: Vec<_> = oxygens
            .iter()
            .copied()
            .filter(|o| !double.contains(o))
            .collect();
        if double.len() != 2
            || single.len() != 1
            || mol.atoms[&single[0]].charge != 0
            || mol.degree(single[0]) != 1
        {
            continue;
        }
        centers.push((c, nitrogens[0], oxygens));
    }
    if centers.len() != 1 {
        return Ok(None);
    }
    let Some((c, n, oxygens)) = centers.pop() else {
        return Ok(None);
    };
    if mol.atoms[&n].charge != 0 {
        return Ok(None);
    }
    let core: AtomSet = [c, n].into_iter().chain(oxygens).collect();
    require_stereo_free_core(mol, &core)?;
    let Some(ligands) = nitrogen_ligands(mol, atoms, n, &core, false)? else {
        return Ok(None);
    };
    let mut represented = core;
    for (branch, _) in &ligands {
        represented.extend(branch);
    }
    if represented != *atoms {
        return Ok(None);
    }
    let names: Vec<_> = ligands.into_iter().map(|(_, n)| n).collect();
    Ok(Some(format!("{}sulfamic acid", center_ligands(&names))))
}
/// P-67.1.2.4.1.1: HO-SO2-CN is sulfurocyanidic acid, sulfuric acid with the
/// infix cyanido (P-67.1.2.3.2). Sulfurous and phosphoric acid take infixes
/// too (P-67.1.2.1), so (HO)2P(O)-CN is phosphorocyanidic acid, like
/// phosphorodicyanatidic acid (P-67.1.2.4.1.3). Anions end in ate or ite and
/// keep `hydrogen` for an acid hydrogen (P-67.1.3.1); esters cite the alkyl
/// groups before that, as in sodium methyl hydrogen phosphate (P-67.1.3.2).
/// A substitutive name on methane would have the nitrile carbon as its
/// whole parent.
// Provenance: fn cyanidic_acid_name | original
// Source: Nazuchem | crates/nazuchem-s2n/src/special_cases.rs::cyanidic_acid_name
fn cyanidic_acid_name(
    mol: &Molecule,
    atoms: &AtomSet,
) -> Result<Option<StructuralParentName>, String> {
    let mut found = None;
    for &centre in atoms {
        let atom = &mol.atoms[&centre];
        if atom.charge != 0 || atom.total_h_count != 0 || atom.has_stereo() {
            continue;
        }
        let around = mol.get_neighbors(centre);
        let cyano: Vec<_> = around
            .iter()
            .copied()
            .filter(|&c| {
                mol.atoms[&c].is_carbon()
                    && mol.degree(c) == 2
                    && mol.get_neighbors(c).iter().any(|&n| {
                        mol.atoms[&n].symbol == "N"
                            && mol.bond_order(c, n) == 3
                            && mol.atoms[&n].charge == 0
                    })
            })
            .collect();
        let oxygen = |order| {
            around
                .iter()
                .copied()
                .filter(|&o| {
                    mol.atoms[&o].symbol == "O"
                        && mol.bond_order(centre, o) == order
                        && (order == 1 || mol.degree(o) == 1 && mol.atoms[&o].charge == 0)
                })
                .collect::<Vec<_>>()
        };
        let (oxo, single) = (oxygen(2), oxygen(1));
        let names = match (atom.symbol.as_str(), oxo.len(), single.len()) {
            ("S", 2, 1) => ("sulfurocyanidic acid", "sulfurocyanidate"),
            ("S", 1, 1) => ("sulfurocyanidous acid", "sulfurocyanidite"),
            ("P", 1, 2) => ("phosphorocyanidic acid", "phosphorocyanidate"),
            _ => continue,
        };
        if cyano.len() == 1 && around.len() == oxo.len() + single.len() + 1 {
            if found.is_some() {
                return Ok(None);
            }
            found = Some((centre, cyano[0], oxo, single, names));
        }
    }
    let Some((centre, c, oxo, single, (acid, anion))) = found else {
        return Ok(None);
    };
    let mut core: AtomSet = [centre, c].into_iter().chain(oxo).collect();
    core.extend(&single);
    core.extend(mol.get_neighbors(c));
    let mut name = StructuralParentName {
        certified: true,
        ..StructuralParentName::ordinary(String::new())
    };
    let mut covered = core.clone();
    let (mut alkyls, mut hydroxy) = (Vec::new(), 0);
    for o in single {
        let roots: Vec<_> = mol
            .get_neighbors(o)
            .iter()
            .copied()
            .filter(|&r| r != centre)
            .collect();
        match (&roots[..], mol.atoms[&o].charge) {
            ([], -1) => {}
            ([], 0) => hydroxy += 1,
            (&[root], 0) if mol.atoms[&root].is_carbon() => {
                let ligand = component_atoms_until_blocked(mol, atoms, root, &core);
                // An acyl group makes an anhydride (P-67.1.3.3) and a principal
                // group or a charge in the ligand may be senior to the ester
                // (P-41); these are left to the other routes.
                if !ligand.is_disjoint(&covered)
                    || ligand.iter().any(|a| mol.atoms[a].charge != 0)
                    || mol
                        .get_neighbors(root)
                        .iter()
                        .any(|&n| mol.bond_order(root, n) != 1 && !mol.atoms[&n].is_carbon())
                    || perceive_groups(mol)
                        .iter()
                        .any(|g| g.is_principal_candidate && g.atom_ids().is_subset(&ligand))
                {
                    return Ok(None);
                }
                let alkyl = ester_modifier_name(mol, atoms, root, o, &core)?;
                if alkyl.is_empty() {
                    return Ok(None);
                }
                name.component_splits.insert(normalize_edge(o, root));
                alkyls.push(alkyl);
                covered.extend(ligand);
            }
            _ => return Ok(None),
        }
    }
    if covered != *atoms {
        return Ok(None);
    }
    let anionic = atoms.iter().any(|a| mol.atoms[a].charge != 0);
    name.name = if alkyls.is_empty() && !anionic {
        acid.into()
    } else {
        let hydrogen = match hydroxy {
            0 => None,
            1 => Some("hydrogen".into()),
            n => Some(format!("{}hydrogen", multipliers::basic(n))),
        };
        let esters = Some(modifier_phrase(&alkyls, false)).filter(|e| !e.is_empty());
        esters
            .into_iter()
            .chain(hydrogen)
            .chain([anion.to_owned()])
            .collect::<Vec<_>>()
            .join(" ")
    };
    require_stereo_free_core(mol, &core)?;
    Ok(Some(name))
}
// Provenance: fn azinic_acid_name | adapted
// Source: openclatura | src/openclatura/special_cases.py::azinic_acid_result
fn azinic_acid_name(mol: &Molecule, atoms: &AtomSet) -> Result<Option<String>, String> {
    let mut centers = Vec::new();
    for &n in atoms {
        if mol.atoms[&n].symbol != "N" || mol.atoms[&n].charge != 1 {
            continue;
        }
        let neighbors: Vec<_> = mol
            .get_neighbors(n)
            .iter()
            .copied()
            .filter(|i| atoms.contains(i))
            .collect();
        let oxygens: Vec<_> = neighbors
            .iter()
            .copied()
            .filter(|i| mol.atoms[i].symbol == "O")
            .collect();
        let others: Vec<_> = neighbors
            .iter()
            .copied()
            .filter(|i| !oxygens.contains(i))
            .collect();
        if oxygens.len() != 2 || others.len() != 1 {
            continue;
        }
        let oxido = oxygens
            .iter()
            .filter(|o| mol.atoms[o].charge == -1 && mol.degree(**o) == 1)
            .count();
        let hydroxy = oxygens
            .iter()
            .filter(|o| {
                mol.atoms[o].charge == 0 && mol.degree(**o) == 1 && mol.atoms[o].total_h_count > 0
            })
            .count();
        if oxido == 1 && hydroxy == 1 {
            centers.push((n, others[0], oxygens));
        }
    }
    if centers.len() != 1 {
        return Ok(None);
    }
    let Some((n, root, oxygens)) = centers.pop() else {
        return Ok(None);
    };
    let bond = mol.bond_order(n, root);
    if !matches!(bond, 1 | 2) {
        return Ok(None);
    }
    let core: AtomSet = std::iter::once(n).chain(oxygens).collect();
    require_stereo_free_core(mol, &core)?;
    let name = strip_outer_parentheses(&branch_name(mol, atoms, root, &core, n)?);
    let branch = component_atoms_until_blocked(mol, atoms, root, &core);
    if name.is_empty()
        || branch.is_empty()
        || branch.union(&core).copied().collect::<AtomSet>() != *atoms
    {
        return Ok(None);
    }
    Ok(Some(if bond == 2 {
        format!("{}azinic acid", format_multiplier(&name, 1, true))
    } else {
        format!("N-{}azinic acid", format_multiplier(&name, 1, false))
    }))
}

// Provenance: fn phosphane_borane_zwitterion_name | adapted
// Source: openclatura | src/openclatura/special_cases.py::phosphane_borane_zwitterion_result
fn phosphane_borane_zwitterion_name(
    mol: &Molecule,
    atoms: &AtomSet,
) -> Result<Option<String>, String> {
    let mut pairs = Vec::new();
    for &p in atoms {
        if mol.atoms[&p].symbol != "P" || mol.atoms[&p].charge <= 0 {
            continue;
        }
        for &b in mol.get_neighbors(p) {
            if atoms.contains(&b)
                && mol.atoms[&b].symbol == "B"
                && mol.atoms[&b].charge < 0
                && mol.atoms[&b].total_h_count == 3
                && mol
                    .get_neighbors(p)
                    .iter()
                    .filter(|n| atoms.contains(n) && **n != b)
                    .count()
                    == 3
            {
                pairs.push((p, b));
            }
        }
    }
    let Some(&(p, b)) = pairs.first() else {
        return Ok(None);
    };
    if mol.atoms[&p].has_stereo() {
        return Ok(None);
    }
    let roots: Vec<_> = mol
        .get_neighbors(p)
        .iter()
        .copied()
        .filter(|n| atoms.contains(n) && *n != b && mol.atoms[n].is_carbon())
        .collect();
    if roots.len() != 3 {
        return Ok(None);
    }
    let core = set([p, b]);
    require_stereo_free_core(mol, &core)?;
    let mut seen = AtomSet::new();
    let mut names = Vec::new();
    for root in roots {
        let branch = component_atoms_until_blocked(mol, atoms, root, &core);
        if branch.is_empty() || !branch.is_disjoint(&seen) {
            return Ok(None);
        }
        seen.extend(branch);
        let mut name = strip_outer_parentheses(&branch_name(mol, atoms, root, &core, p)?);
        if name.is_empty() {
            name = alkyl_ligand_name(mol, atoms, root, p)?;
        }
        if name.is_empty() {
            return Ok(None);
        }
        names.push(name);
    }
    if seen.union(&core).copied().collect::<AtomSet>() != *atoms {
        return Ok(None);
    }
    Ok(Some(format!(
        "({}phosphaniumyl)boranuide",
        center_ligands(&names)
    )))
}

// Provenance: fn sulfonium_ylide_name | adapted
// Source: openclatura | src/openclatura/special_cases.py::sulfonium_ylide_name
fn sulfonium_ylide_name(mol: &Molecule, atoms: &AtomSet) -> Result<Option<String>, String> {
    let sulfurs: Vec<_> = atoms
        .iter()
        .copied()
        .filter(|i| mol.atoms[i].symbol == "S" && mol.atoms[i].charge > 0)
        .collect();
    let carbons: Vec<_> = atoms
        .iter()
        .copied()
        .filter(|i| mol.atoms[i].is_carbon() && mol.atoms[i].charge < 0)
        .collect();
    if sulfurs.len() != 1 || carbons.len() != 1 {
        return Ok(None);
    }
    let (s, c) = (sulfurs[0], carbons[0]);
    if mol.atoms[&s].has_stereo()
        || mol.atoms[&s].total_h_count != 0
        || mol.atoms[&s].explicit_h_count != 0
        || mol.bond_order(s, c) != 1
    {
        return Ok(None);
    }
    let roots: Vec<_> = mol
        .get_neighbors(s)
        .iter()
        .copied()
        .filter(|n| *n != c && atoms.contains(n) && mol.atoms[n].is_carbon())
        .collect();
    if roots.len() != 2 {
        return Ok(None);
    }
    require_stereo_free_core(mol, &set([s, c]))?;
    let mut sulfur_atoms = AtomSet::new();
    let mut names = Vec::new();
    for root in roots {
        let branch = carbon_ligand_atoms(mol, atoms, root, s);
        if branch.is_empty() || !branch.is_disjoint(&sulfur_atoms) {
            return Ok(None);
        }
        sulfur_atoms.extend(branch);
        let name = alkyl_ligand_name(mol, atoms, root, s)?;
        if name.is_empty() {
            return Ok(None);
        }
        names.push(name);
    }
    let subs: Vec<_> = mol
        .get_neighbors(c)
        .iter()
        .copied()
        .filter(|n| *n != s && atoms.contains(n) && mol.atoms[n].symbol != "H")
        .collect();
    if subs.len() > 1 {
        return Ok(None);
    }
    let (mut sub_atoms, mut sub_name) = (AtomSet::new(), String::new());
    if let Some(&root) = subs.first() {
        let mut blocked = set([s, c]);
        blocked.extend(&sulfur_atoms);
        sub_atoms = component_atoms_until_blocked(mol, atoms, root, &blocked);
        if sub_atoms.is_empty() || !sub_atoms.is_disjoint(&sulfur_atoms) {
            return Ok(None);
        }
        sub_name = strip_outer_parentheses(&branch_name(mol, atoms, root, &blocked, c)?);
        if sub_name.is_empty() {
            sub_name = alkyl_ligand_name(mol, atoms, root, c)?;
        }
        if sub_name.is_empty() {
            return Ok(None);
        }
    }
    let mut represented = set([s, c]);
    represented.extend(sulfur_atoms);
    represented.extend(&sub_atoms);
    if represented != *atoms {
        return Ok(None);
    }
    let prefix = format!("{}sulfaniumyl", center_ligands(&names));
    let parent = if subs.is_empty() {
        "methanide".to_owned()
    } else {
        let degrees: Vec<_> = sub_atoms
            .iter()
            .map(|i| {
                mol.get_neighbors(*i)
                    .iter()
                    .filter(|n| sub_atoms.contains(n))
                    .count()
            })
            .collect();
        let linear = sub_atoms.iter().all(|i| mol.atoms[i].is_carbon())
            && degrees.iter().sum::<usize>() / 2 == sub_atoms.len() - 1
            && degrees.iter().all(|d| *d <= 2)
            && sub_atoms.iter().all(|i| {
                mol.get_neighbors(*i)
                    .iter()
                    .all(|n| !sub_atoms.contains(n) || mol.bond_order(*i, *n) == 1)
            });
        if linear {
            return Ok(Some(format!(
                "1-({prefix}){}an-1-ide",
                stems::stem_for(sub_atoms.len() + 1)
            )));
        }
        if sub_name.ends_with("yl") {
            format!("{sub_name}methanide")
        } else {
            format!("({sub_name})methanide")
        }
    };
    // This shortcut fixes the parent at the ylide carbon. Another anionic
    // centre in the substituent may belong to a senior parent (P-72.7(a)).
    if sub_atoms.iter().map(|i| mol.atoms[i].charge).sum::<i32>() < 0 {
        mol.mark_name(NameStatus::Unverified(NameReason::AnionParent));
    }
    Ok(Some(format!("({prefix}){parent}")))
}

// Provenance: fn nitrogen_chalcogenide_name | adapted
// Source: openclatura | src/openclatura/special_cases.py::nitrogen_chalcogenide_result
fn nitrogen_chalcogenide_name(
    mol: &Molecule,
    atoms: &AtomSet,
) -> Result<Option<StructuralParentName>, String> {
    use super::nomenclature::NAMER_RULES;
    let mut matches = Vec::new();
    for &n in atoms {
        let a = &mol.atoms[&n];
        if a.symbol != "N" || a.charge <= 0 || a.is_aromatic {
            continue;
        }
        for &ch in mol.get_neighbors(n) {
            let e = &mol.atoms[&ch];
            if !atoms.contains(&ch)
                || e.charge >= 0
                || mol.degree(ch) != 1
                || mol.bond_order(n, ch) != 1
                || NAMER_RULES
                    .chalcogen_elements
                    .get(&e.symbol)
                    .and_then(|c| c.nitrogen_chalcogenide_citation.as_deref())
                    != Some("functional_class")
            {
                continue;
            }
            let others: Vec<_> = mol
                .get_neighbors(n)
                .iter()
                .copied()
                .filter(|i| *i != ch)
                .collect();
            if others.is_empty() || others.iter().any(|i| !mol.atoms[i].is_carbon()) {
                continue;
            }
            let orders: Vec<_> = others.iter().map(|i| mol.bond_order(n, *i)).collect();
            let imine = orders.iter().filter(|o| **o == 2).count() == 1
                && orders.iter().all(|o| matches!(o, 1 | 2));
            let amine = orders.iter().all(|o| *o == 1);
            if !imine && !amine {
                continue;
            }
            matches.push((n, ch, imine));
        }
    }
    if matches.len() != 1 {
        return Ok(None);
    }
    let (n, ch, imine) = matches[0];
    require_stereo_free_core(mol, &set([n, ch]))?;
    let reduced_atoms: AtomSet = atoms.iter().copied().filter(|i| *i != ch).collect();
    let mut reduced = mol.subgraph(&reduced_atoms);
    reduced.set_atom_charge(n, 0)?;
    let parent =
        super::component_namer::name_component_with_citation(&reduced, &reduced_atoms, None, None)?;
    // P-62.5: the class name only follows a parent amine or imine; a
    // substituent N-chalcogenide is cited as a prefix (P-73.6).
    let names_the_amine = match &parent.derivation {
        super::component_namer::Derivation::Plan {
            parent_atoms,
            principal_group,
            ..
        } => match principal_group {
            None => parent_atoms.contains(&n),
            Some((key, group)) => matches!(key.as_str(), "amine" | "imine") && group.contains(&n),
        },
        super::component_namer::Derivation::Shortcut { .. } => false,
    };
    if parent.text.is_empty() || !names_the_amine {
        return Ok(None);
    }
    let Some(class) = NAMER_RULES
        .chalcogenide_class_names
        .get(&mol.atoms[&ch].symbol)
    else {
        return Ok(None);
    };
    mol.mark_name(NameStatus::Unverified(NameReason::AmineChalcogenide));
    Ok(Some(StructuralParentName {
        component_splits: parent.derivation.component_splits().clone(),
        ..StructuralParentName::ordinary(if imine {
            format!("{} N-{class}", parent.text)
        } else {
            format!("{} {class}", parent.text)
        })
    }))
}

// Provenance: fn peroxide_linkage_name | adapted
// Source: openclatura | src/openclatura/special_cases.py::peroxide_linkage_result
fn peroxide_linkage_name(mol: &Molecule, atoms: &AtomSet) -> Result<Option<String>, String> {
    use super::nomenclature::NAMER_RULES;
    let is_chalcogen = |i: usize| {
        NAMER_RULES
            .chalcogen_elements
            .contains_key(&mol.atoms[&i].symbol)
    };
    let mut matches = Vec::new();
    for bond in mol.bonds.values() {
        let (l, r) = (bond.u, bond.v);
        if !atoms.contains(&l)
            || !atoms.contains(&r)
            || bond.order != 1
            || !is_chalcogen(l)
            || !is_chalcogen(r)
            || [l, r]
                .iter()
                .any(|i| mol.atoms[i].charge != 0 || mol.atoms[i].radical_electrons != 0)
        {
            continue;
        }
        let left: Vec<_> = mol
            .get_neighbors(l)
            .iter()
            .copied()
            .filter(|n| *n != r)
            .collect();
        let right: Vec<_> = mol
            .get_neighbors(r)
            .iter()
            .copied()
            .filter(|n| *n != l)
            .collect();
        if left.len() != 1
            || right.len() != 1
            || !mol.atoms[&left[0]].is_carbon()
            || !mol.atoms[&right[0]].is_carbon()
        {
            continue;
        }
        let blocked = set([l, r]);
        let branch = subgraph_component(mol, left[0], &blocked);
        if branch.contains(&right[0]) {
            continue;
        }
        if [(left[0], l), (right[0], r)]
            .iter()
            .any(|(attachment, linker)| {
                mol.get_neighbors(*attachment).iter().any(|n| {
                    n != linker && is_chalcogen(*n) && mol.bond_order(*attachment, *n) == 2
                })
            })
        {
            continue;
        }
        matches.push((l, r, left[0], right[0]));
    }
    if matches.len() != 1 {
        return Ok(None);
    }
    let (l, r, la, ra) = matches[0];
    let (le, re) = (&mol.atoms[&l].symbol, &mol.atoms[&r].symbol);
    if le == re {
        return Ok(None);
    }
    // The O/S pair deliberately reuses the ordinary substituent pipeline.
    if [le.as_str(), re.as_str()]
        .iter()
        .all(|s| matches!(*s, "O" | "S"))
    {
        return Ok(None);
    }
    let mut unique = vec![le.as_str(), re.as_str()];
    unique.sort();
    unique.dedup();
    let Some(class) = NAMER_RULES
        .chalcogen_peroxide_class_names
        .iter()
        .find(|row| {
            let mut e: Vec<_> = row.elements.iter().map(String::as_str).collect();
            e.sort();
            e == unique
        })
        .map(|row| &row.name)
    else {
        return Ok(None);
    };
    let blocked = set([l, r]);
    require_stereo_free_core(mol, &blocked)?;
    let left = strip_outer_parentheses(&branch_name(mol, atoms, la, &blocked, l)?);
    let right = strip_outer_parentheses(&branch_name(mol, atoms, ra, &blocked, r)?);
    if left.is_empty() || right.is_empty() {
        return Ok(None);
    }
    let phrase = if left == right {
        format_multiplier(&left, 2, false)
    } else {
        let mut located = [format!("{le}-{left}"), format!("{re}-{right}")];
        located.sort_by_key(|s| substituent_sort_key(s));
        located.join(" ")
    };
    mol.mark_name(NameStatus::NotPreferred(NameReason::FunctionalClass));
    Ok(Some(format!("{phrase} {class}")))
}

// Provenance: fn anhydride_half_atoms | adapted
// Source: openclatura | src/openclatura/special_cases.py::_anhydride_half_atoms
fn anhydride_half_atoms(mol: &Molecule, start: usize, bridge: &AtomSet) -> AtomSet {
    subgraph_component(mol, start, bridge)
}
// Provenance: fn anhydride_half_name | adapted
// Source: openclatura | src/openclatura/special_cases.py::anhydride_half_name
/// Names the acid of one anhydride half, with a hydroxy group on each of
/// `starts`. Also tells whether every one of them is cited as a suffix, and
/// whether the acid is unsubstituted (all atoms in its parent or suffixes).
fn anhydride_half_name(
    mol: &Molecule,
    starts: &[usize],
    bridge: &AtomSet,
) -> Result<(String, bool, bool), String> {
    let mut atoms = anhydride_half_atoms(mol, starts[0], bridge);
    let mut fragment = mol.subgraph(&atoms);
    let first_oh = mol.atoms.keys().max().copied().unwrap_or(0) + 100;
    let mut hydroxy = Vec::new();
    for (oh, &start) in (first_oh..).zip(starts) {
        let mut oxygen = Atom::new(oh, "O");
        oxygen.total_h_count = 1;
        fragment.add_atom(oxygen)?;
        let bond = fragment.bonds.keys().max().copied().unwrap_or(0) + 1;
        fragment.add_bond(Bond::new(bond, start, oh, 1))?;
        atoms.insert(oh);
        hydroxy.push(oh);
    }
    let name = super::component_namer::name_component_with_citation(
        &fragment,
        &atoms,
        Some("anhydride_half"),
        None,
    )?;
    let cited = hydroxy
        .iter()
        .all(|oh| name.principal_group_atoms.contains(oh));
    let unsubstituted = matches!(&name.derivation, super::component_namer::Derivation::Plan { parent_atoms, .. }
        if atoms.iter().all(|a| parent_atoms.contains(a) || name.principal_group_atoms.contains(a)));
    // An acid named as a whole, such as hydrazinecarboxylic acid, keeps its
    // class word, which an anhydride name omits (P-65.7.2). Only for a
    // carbon acid: the anhydrides of sulfurocyanidic and sulfurocyanidous
    // acid (issue #158) are not generated, so those halves keep the word
    // and the callers decline them.
    let text = if mol.atoms[&starts[0]].is_carbon() {
        name.text.strip_suffix(" acid").unwrap_or(&name.text)
    } else {
        &name.text
    };
    Ok((text.to_owned(), cited, unsubstituted))
}
// Provenance: fn anhydride_class_name | adapted
// Source: openclatura | src/openclatura/chalcogen_vocabulary.py::anhydride_class_name
fn anhydride_class_name(mol: &Molecule, path: &[usize]) -> Option<String> {
    use super::nomenclature::{NAMER_RULES, chalcogen_replacement_infix};
    let elements: Vec<_> = path[1..path.len() - 1]
        .iter()
        .map(|i| mol.atoms[i].symbol.as_str())
        .collect();
    let (template, replacement) = match elements.as_slice() {
        ["O"] => ("anhydride", String::new()),
        [element] => (
            "single_replacement_anhydride",
            chalcogen_replacement_infix(element)?.to_owned(),
        ),
        ["O", "O"] => ("peroxyanhydride", String::new()),
        [a, b] => {
            let replacements: Vec<_> = [a, b].into_iter().copied().filter(|e| *e != "O").collect();
            if replacements.len() == 1 {
                (
                    "single_replacement_peroxyanhydride",
                    chalcogen_replacement_infix(replacements[0])?.to_owned(),
                )
            } else if a == b {
                (
                    "double_replacement_peroxyanhydride",
                    chalcogen_replacement_infix(a)?.to_owned(),
                )
            } else {
                let mut infixes = Vec::new();
                for e in replacements {
                    infixes.push(chalcogen_replacement_infix(e)?);
                }
                infixes.sort();
                ("mixed_replacement_peroxyanhydride", infixes.concat())
            }
        }
        _ => return None,
    };
    Some(
        NAMER_RULES
            .chalcogen_templates
            .get(template)?
            .replace("{replacement}", &replacement),
    )
}
type NameLocant = (bool, u32, String, usize);
pub(super) type AlphanumericalKey = (String, String, Vec<NameLocant>, Vec<NameLocant>, Vec<u8>);

/// The alphanumerical order of a name (P-14.5): its nonitalic Roman letters,
/// its italic letters, its locants outside enclosing marks, all its locants,
/// then its enclosing marks and punctuation.
// Provenance: fn alphanumerical_order_key | original
// Source: ReShiki | src/naming/openclatura/special_cases.rs::try_name_anhydride_component
pub(super) fn alphanumerical_order_key(name: &str) -> AlphanumericalKey {
    const ITALIC_WORDS: [&str; 10] = [
        "tert", "sec", "cis", "trans", "endo", "exo", "syn", "anti", "as", "asym",
    ];
    let chars: Vec<char> = name.chars().collect();
    let (mut roman, mut italic) = (String::new(), String::new());
    let (mut locants, mut marks) = (Vec::new(), Vec::new());
    let (mut depth, mut i) = (0usize, 0);
    while i < chars.len() {
        let c = chars[i];
        let run = |from: usize, f: fn(&char) -> bool| {
            from + chars[from..].iter().take_while(|c| f(c)).count()
        };
        if c == '(' {
            // A descriptor such as (S), (2E) or (1r,4r) has a
            // letter, no hyphen and no word of two letters; it
            // is dropped with its hyphen.
            let end = chars[i..]
                .iter()
                .position(|&c| c == ')')
                .map_or(chars.len(), |n| i + n);
            let inner = &chars[i + 1..end];
            if inner.iter().any(char::is_ascii_alphabetic)
                && !inner.contains(&'-')
                && !inner
                    .windows(2)
                    .any(|w| w.iter().all(char::is_ascii_lowercase))
            {
                i = end + 1 + usize::from(chars.get(end + 1) == Some(&'-'));
                continue;
            }
        }
        let start = i;
        if c.is_ascii_digit() {
            let digits = run(i, char::is_ascii_digit);
            let number: u32 = chars[i..digits]
                .iter()
                .collect::<String>()
                .parse()
                .unwrap_or(u32::MAX);
            // Letters and primes written onto a locant: 4a, 3'a,
            // 1H, 3aH, 1lambda.
            i = run(digits, |c| c.is_ascii_alphabetic() || *c == '\'');
            let tail = &chars[digits..i];
            let letters: Vec<char> = tail
                .iter()
                .filter(|c| c.is_ascii_alphabetic())
                .copied()
                .collect();
            let letter = match letters[..] {
                [l] | [l, 'A'..='Z', ..] if l.is_ascii_lowercase() => l.into(),
                _ => String::new(),
            };
            let primes = tail.iter().filter(|&&c| c == '\'').count();
            locants.push((depth, (true, number, letter, primes)));
        } else if c.is_ascii_alphabetic() {
            let end = run(i, char::is_ascii_alphabetic);
            i = run(end, |c| *c == '\'');
            let word: String = chars[start..end].iter().collect();
            if word.contains(|c: char| c.is_ascii_uppercase()) {
                // An element locant such as N or N'.
                locants.push((depth, (false, 0, word, i - end)));
            } else {
                // Fusion letters run in alphabetical sequence: c, cd.
                let fusion = word.as_bytes().windows(2).all(|w| w[1] == w[0] + 1);
                if fusion || ITALIC_WORDS.contains(&word.as_str()) {
                    italic.push_str(&word);
                } else if word != "lambda" {
                    // The Greek lambda of lambda^4 takes no part.
                    roman.push_str(&word);
                }
            }
        } else {
            match c {
                '(' | '[' | '{' => depth += 1,
                ')' | ']' | '}' => depth = depth.saturating_sub(1),
                _ => {}
            }
            i += 1;
            // The bonding number of lambda^4 is no locant.
            if c == '^' {
                i = run(i, char::is_ascii_digit);
            }
        }
        marks.extend(chars[start..i].iter().map(|c| match c {
            '{' | '}' => 1,
            '[' | ']' => 2,
            '(' | ')' => 3,
            '.' => 4,
            ',' => 5,
            ';' => 6,
            ':' => 7,
            '-' => 8,
            _ => 0,
        }));
    }
    let primary: Vec<_> = locants
        .iter()
        .filter(|l| l.0 == 0)
        .map(|l| l.1.clone())
        .collect();
    let all: Vec<_> = locants.into_iter().map(|l| l.1).collect();
    (roman, italic, primary, all, marks)
}
// Provenance: fn try_name_anhydride_component | adapted
// Source: openclatura | src/openclatura/special_cases.py::try_name_anhydride_component_result
pub(super) fn try_name_anhydride_component(
    mol: &Molecule,
    groups: &[PerceivedGroup],
    principal_key: Option<&str>,
) -> Result<Option<StructuralParentName>, String> {
    // The caller supplies component-scoped groups. Count all anhydride types,
    // including those with a different principal key (P-65.7.6).
    let anhydride_count = groups
        .iter()
        .filter(|group| {
            group
                .resolved_rule
                .as_ref()
                .is_some_and(|rule| rule.has_capability("anhydride"))
        })
        .count();
    // P-65.7.6.1: a dibasic acid whose two acid groups each form an anhydride
    // with a monobasic acid gives a dianhydride, as in "acetic
    // methanedisulfonic propanoic dianhydride".
    let links: Vec<_> = groups
        .iter()
        .filter(|g| {
            g.resolved_rule
                .as_ref()
                .is_some_and(|r| r.has_capability("anhydride"))
        })
        .collect();
    if let [g1, g2] = links[..]
        && [g1, g2]
            .iter()
            .any(|g| Some(g.key.as_str()) == principal_key)
        && let (Some(d1), Some(d2)) = (&g1.descriptor, &g2.descriptor)
        && let (Some(p1), Some(p2)) = (d1.linker_paths.first(), d2.linker_paths.first())
    {
        let result = mol.candidate(|| {
            let paths = [p1, p2];
            // Single-atom bridges of one element (P-65.7.6.4.1); P-65.7.6.4.2
            // names unlike bridges substitutively.
            if paths.iter().any(|p| p.len() != 3)
                || mol.atoms[&paths[0][1]].symbol != mol.atoms[&paths[1][1]].symbol
            {
                return Ok(None);
            }
            let bridges: AtomSet = paths.iter().map(|p| p[1]).collect();
            let component = |atom| subgraph_component(mol, atom, &bridges);
            // The dibasic acid holds one center of each link; each end holds
            // a single acid group.
            let Some((m1, m2, e1, e2)) = (0..2)
                .flat_map(|i| (0..2).map(move |j| (i, j)))
                .find(|&(i, j)| component(d1.centers[i]).contains(&d2.centers[j]))
                .map(|(i, j)| {
                    (
                        d1.centers[i],
                        d2.centers[j],
                        d1.centers[1 - i],
                        d2.centers[1 - j],
                    )
                })
            else {
                return Ok(None);
            };
            let (end1, end2, central) = (component(e1), component(e2), component(m1));
            if end1.contains(&m1) || end1.contains(&e2) || end2.contains(&m2) {
                return Ok(None);
            }
            let (left, _, plain_end) = anhydride_half_name(mol, &[e1], &bridges)?;
            let (right, ..) = anhydride_half_name(mol, &[e2], &bridges)?;
            let (dibasic, cited, plain) = anhydride_half_name(mol, &[m1, m2], &bridges)?;
            // A multiplied thioic suffix is cited as bis(thioic) (P-65.7.6.4.3),
            // which the acid namer does not yet write.
            let thioxo = [m1, m2].iter().any(|&m| {
                mol.get_neighbors(m)
                    .iter()
                    .any(|n| mol.bond_order(m, *n) == 2 && mol.atoms[n].symbol != "O")
            });
            if !cited || thioxo || !left.ends_with("ic") || !right.ends_with("ic") {
                return Ok(None);
            }
            // Different ends on unlike acid groups of the dibasic acid need
            // locants (P-65.7.6.2, P-65.7.6.3), which are not built here. The
            // ranks ignore configuration, so any stereo counts as unlike.
            let ranks = canonical_ranks(mol, Some(&central));
            if left != right
                && (ranks[&m1] != ranks[&m2] || require_stereo_free_core(mol, &central).is_err())
            {
                return Ok(None);
            }
            let Some(class) = anhydride_class_name(mol, paths[0]) else {
                return Ok(None);
            };
            require_stereo_free_core(mol, &bridges)?;
            let class = if class == "anhydride" {
                "dianhydride".to_owned()
            } else {
                format!("bis({class})")
            };
            // The acids are cited in their order of occurrence, beginning
            // with the end acid lower in alphanumerical order (P-14.5). Two
            // equal ends take "bis" when substituted (P-16.3.5), else "di",
            // enclosed when the acid begins with a numerical term or cites
            // locants or brackets, as in di(decanoic) (P-16.3.4 (c)-(f)).
            let name = if left == right && !plain_end {
                format!("bis({left}) {dibasic} {class}")
            } else if left == right
                && left.bytes().all(|b| b.is_ascii_lowercase())
                && !super::formatting::starts_with_multiplier(&left)
            {
                format!("di{left} {dibasic} {class}")
            } else if left == right {
                format!("di({left}) {dibasic} {class}")
            } else {
                let mut ends = [left, right];
                ends.sort_by_cached_key(|end| alphanumerical_order_key(end));
                format!("{} {dibasic} {} {class}", ends[0], ends[1])
            };
            // P-65.7.6.1 cites "2-methylbutanedioic" where P-14.3.4.3 omits
            // the locant of a substituted dibasic acid; such names stay
            // unverified.
            // The di/bis choice for a carbamic or hydrazinecarboxylic end
            // (P-16.3.3, P-16.3.5) is not made here.
            let nitrogen_end = [e1, e2].iter().any(|&e| {
                mol.get_neighbors(e)
                    .iter()
                    .any(|n| mol.atoms[n].symbol == "N")
            });
            Ok(Some(StructuralParentName {
                certified: class == "dianhydride"
                    && plain
                    && !nitrogen_end
                    && bridges.is_disjoint(&get_cyclic_atoms(mol, &AtomSet::new())),
                component_splits: paths
                    .iter()
                    .flat_map(|p| p.windows(2).map(|w| normalize_edge(w[0], w[1])))
                    .collect(),
                ..StructuralParentName::ordinary(name)
            }))
        })?;
        if result.is_some() {
            return Ok(result);
        }
    }
    for group in groups {
        let result = mol.candidate(|| {
            if Some(group.key.as_str()) != principal_key
                || !group
                    .resolved_rule
                    .as_ref()
                    .is_some_and(|r| r.has_capability("anhydride"))
            {
                return Ok(None);
            }
            let Some(descriptor) = group.descriptor.as_ref() else {
                return Ok(None);
            };
            if descriptor.centers.len() != 2 || descriptor.linker_paths.is_empty() {
                return Ok(None);
            }
            let path = &descriptor.linker_paths[0];
            if !(3..=4).contains(&path.len()) {
                return Ok(None);
            }
            let bridge: AtomSet = path[1..path.len() - 1].iter().copied().collect();
            if bridge.is_empty() {
                return Ok(None);
            }
            let (left, ..) = anhydride_half_name(mol, &descriptor.centers[..1], &bridge)?;
            let (right, ..) = anhydride_half_name(mol, &descriptor.centers[1..], &bridge)?;
            // Anhydride names cite acid names (P-65.7.1). A half named for another
            // principal group, such as a nitrile, would read as a separate compound.
            if !left.ends_with("ic") || !right.ends_with("ic") {
                return Ok(None);
            }
            let Some(class) = anhydride_class_name(mol, path) else {
                return Ok(None);
            };
            require_stereo_free_core(mol, &bridge)?;
            // P-65.7.1–P-65.7.4: O bridges, or chalcogen/peroxy bridges between
            // graph-equivalent acid halves.
            let ranks = canonical_ranks(mol, None);
            let certified = anhydride_count == 1
                && bridge.is_disjoint(&get_cyclic_atoms(mol, &AtomSet::new()))
                && (path.len() == 3 && mol.atoms[&path[1]].symbol == "O"
                    || ranks[&descriptor.centers[0]] == ranks[&descriptor.centers[1]]);
            let name = if left == right {
                format!("{left} {class}")
            } else {
                let mut halves = [(left.clone(), left), (right.clone(), right)];
                if path.len() == 4 && mol.atoms[&path[1]].symbol != mol.atoms[&path[2]].symbol {
                    halves[0].1 = format!("{}-{}", mol.atoms[&path[1]].symbol, halves[0].1);
                    halves[1].1 = format!("{}-{}", mol.atoms[&path[2]].symbol, halves[1].1);
                }
                // P-65.7.2 cites the acids in alphabetical order, the
                // alphanumerical order of P-14.5: nonitalic Roman letters first,
                // then italic letters (P-14.5.3), then the locants, lowest at
                // the first point of difference (P-14.5.4): first those outside
                // enclosing marks, then all, in order of citation, with element
                // locants (N) before numerals. Locant letters (4a, 1H, 1lambda)
                // and stereodescriptors such as (S) or (1r,4r) give no Roman
                // letters; fusion letters ([3,2-c], [cd]) and words such as tert
                // are italic. Remaining ties compare enclosing marks and
                // punctuation, stereodescriptors excluded (P-14.6).
                halves.sort_by_cached_key(|half| alphanumerical_order_key(&half.0));
                format!("{} {} {class}", halves[0].1, halves[1].1)
            };
            Ok(Some(StructuralParentName {
                certified,
                component_splits: path
                    .windows(2)
                    .map(|w| normalize_edge(w[0], w[1]))
                    .collect(),
                ..StructuralParentName::ordinary(name)
            }))
        })?;
        if result.is_some() {
            return Ok(result);
        }
    }
    Ok(None)
}

#[derive(Clone, Debug)]
pub(super) struct StructuralParentName {
    pub(super) name: String,
    pub(super) certified: bool,
    /// Boundaries between the components of this functional-class name.
    pub(super) component_splits: BTreeSet<Edge>,
    /// Only the homonuclear-chain route supplies the chain audit plan used by
    /// the upstream component-level suffix preference exception.
    pub(super) is_chain: bool,
}

impl StructuralParentName {
    // Provenance: fn StructuralParentName::ordinary | adapted
    // Source: openclatura | src/openclatura/special_cases.py::_component_name_result
    // Source: openclatura | src/openclatura/special_cases.py::SpecialComponentName
    pub(super) fn ordinary(name: String) -> Self {
        Self {
            name,
            certified: false,
            component_splits: BTreeSet::new(),
            is_chain: false,
        }
    }
}

// Provenance: fn blocking_suffix_defers_to_chain_parent | adapted
// Source: openclatura | src/openclatura/component_namer.py::_blocking_suffix_defers_to_chain_parent
pub(super) fn blocking_suffix_defers_to_chain_parent(
    mol: &Molecule,
    atoms: &AtomSet,
    groups: &[PerceivedGroup],
    principal_key: Option<&str>,
    result: Option<&StructuralParentName>,
) -> bool {
    if !result.is_some_and(|result| result.is_chain) {
        return false;
    }
    let characteristic: AtomSet = groups
        .iter()
        .filter(|group| Some(group.key.as_str()) == principal_key)
        .filter(|group| {
            group
                .resolved_rule
                .as_ref()
                .is_some_and(|rule| rule.has_capability("blocks_structural_parent"))
        })
        .flat_map(|group| {
            group
                .descriptor
                .as_ref()
                .map_or_else(|| group.atom_ids(), |descriptor| descriptor.atom_ids())
        })
        .collect();
    atoms.iter().any(|idx| {
        !characteristic.contains(idx)
            && !mol.atoms[idx].is_carbon()
            && mol.atoms[idx].symbol != "H"
            && mol.degree(*idx) > 1
    })
}

#[cfg(test)]
// Provenance: fn structural_replacement_parent_name | adapted
// Source: openclatura | src/openclatura/special_cases.py::structural_replacement_parent_name
fn structural_replacement_parent_name(
    mol: &Molecule,
    atoms: &AtomSet,
) -> Result<Option<String>, String> {
    Ok(structural_replacement_parent_result(mol, atoms)?.map(|result| result.name))
}

// Provenance: fn structural_replacement_parent_result | adapted
// Source: openclatura | src/openclatura/special_cases.py::structural_replacement_parent_result
pub(super) fn structural_replacement_parent_result(
    mol: &Molecule,
    atoms: &AtomSet,
) -> Result<Option<StructuralParentName>, String> {
    // Keep the Python shortcut precedence before ordinary chain and ring parents.
    if let Some(name) = small_inorganic_parent_name(mol, atoms) {
        require_stereo_free_core(mol, atoms)?;
        return Ok(Some(StructuralParentName::ordinary(name)));
    }
    if let Some(name) = mol.candidate(|| hydrazine_carbonyl_parent_name(mol, atoms))? {
        return Ok(Some(StructuralParentName {
            certified: true,
            ..StructuralParentName::ordinary(name)
        }));
    }
    if let Some(name) = mol.candidate(|| cyanohydrazide_name(mol, atoms))? {
        return Ok(Some(StructuralParentName {
            certified: true,
            ..StructuralParentName::ordinary(name)
        }));
    }
    if let Some(name) = mol.candidate(|| simple_azine_parent_name(mol, atoms))? {
        return Ok(Some(StructuralParentName::ordinary(name)));
    }
    if let Some(name) = mol.candidate(|| phosphane_borane_zwitterion_name(mol, atoms))? {
        return Ok(Some(StructuralParentName::ordinary(name)));
    }
    if let Some(name) = mol.candidate(|| sulfonium_ylide_name(mol, atoms))? {
        return Ok(Some(StructuralParentName::ordinary(name)));
    }
    if let Some(name) = mol.candidate(|| nitrogen_chalcogenide_name(mol, atoms))? {
        return Ok(Some(name));
    }
    if let Some(name) = mol.candidate(|| peroxide_linkage_name(mol, atoms))? {
        return Ok(Some(StructuralParentName::ordinary(name)));
    }
    if let Some(name) = mol.candidate(|| hydroxyurea_parent_name(mol, atoms))? {
        return Ok(Some(StructuralParentName::ordinary(name)));
    }
    if let Some(name) = mol.candidate(|| sulfamic_acid_name(mol, atoms))? {
        return Ok(Some(StructuralParentName::ordinary(name)));
    }
    if let Some(name) = mol.candidate(|| cyanidic_acid_name(mol, atoms))? {
        return Ok(Some(name));
    }
    if let Some(name) = mol.candidate(|| azinic_acid_name(mol, atoms))? {
        return Ok(Some(StructuralParentName::ordinary(name)));
    }
    if let Some(name) = mol.candidate(|| oxoacid_ester_name(mol, atoms))? {
        return Ok(Some(name));
    }
    if let Some(name) = oxoacid_parent_name(mol, atoms) {
        require_stereo_free_core(mol, atoms)?;
        return Ok(Some(StructuralParentName::ordinary(name)));
    }
    if let Some(name) = mol.candidate(|| organophosphorus_acid_name(mol, atoms))? {
        return Ok(Some(name));
    }
    if let Some(name) = mol.candidate(|| organoboronic_acid_name(mol, atoms))? {
        return Ok(Some(name));
    }
    if let Some(name) = mol.candidate(|| sulfoxide_parent_name(mol, atoms))? {
        return Ok(Some(StructuralParentName::ordinary(name)));
    }
    if let Some(name) = mol.candidate(|| homonuclear_chain_parent_name(mol, atoms))? {
        return Ok(Some(StructuralParentName {
            is_chain: true,
            ..StructuralParentName::ordinary(name)
        }));
    }
    mol.candidate(|| simple_central_parent_hydride_name(mol, atoms))
}

#[cfg(test)]
// Provenance: mod tests | original
// Source: ReShiki | src/naming/openclatura/special_cases.rs::tests
mod tests {
    use super::super::graph_io::read_smiles;
    use super::*;

    fn special(smiles: &str) -> Option<String> {
        let mol = read_smiles(smiles, &Default::default()).unwrap();
        structural_replacement_parent_name(&mol, &all_atoms(&mol)).unwrap()
    }

    #[test]
    fn hydrazine_carries_carbon_acid_suffixes() {
        // Issue #17. An acid, ester, acid halide, amide or hydrazide carbon
        // on hydrazine is a suffix on hydrazine, numbered 1 and 2
        // (P-65.1.2.2.2, P-68.3.1.2.1).
        for (smiles, expected) in [
            ("NNC(=O)O", "hydrazinecarboxylic acid"),
            ("CNNC(=O)O", "2-methylhydrazine-1-carboxylic acid"),
            ("OC(=O)N(C)N", "1-methylhydrazine-1-carboxylic acid"),
            ("OC(=O)NN(C)C", "2,2-dimethylhydrazine-1-carboxylic acid"),
            ("OC(=O)NNC(=O)O", "hydrazine-1,2-dicarboxylic acid"),
            (
                "CN(C(=O)O)NC(=O)O",
                "1-methylhydrazine-1,2-dicarboxylic acid",
            ),
            // An acyl hydrazide is junior to the acid (P-41).
            ("CC(=O)NNC(=O)O", "2-acetylhydrazine-1-carboxylic acid"),
            // P-14.3.4.5: every position substituted alike.
            ("CN(C)N(C)C(=O)O", "trimethylhydrazinecarboxylic acid"),
            ("CN(C)N(C)C(=O)N(C)C", "pentamethylhydrazinecarboxamide"),
            // P-65.1.4.1: peroxy acid.
            ("NNC(=O)OO", "hydrazinecarboperoxoic acid"),
            // P-65.1.5.1: chalcogen acids, with the O/S locant of the site.
            ("NNC(=S)O", "hydrazinecarbothioic O-acid"),
            ("NNC(=O)S", "hydrazinecarbothioic S-acid"),
            ("NNC(=S)S", "hydrazinecarbodithioic acid"),
            ("NNC(=O)[SeH]", "hydrazinecarboselenoic Se-acid"),
            // P-65.6.3.2.1: esters.
            ("NNC(=O)OC", "methyl hydrazinecarboxylate"),
            ("CNNC(=O)OC", "methyl 2-methylhydrazine-1-carboxylate"),
            ("NNC(=O)OC(C)(C)C", "tert-butyl hydrazinecarboxylate"),
            ("NNC(=S)OC", "O-methyl hydrazinecarbothioate"),
            ("NNC(=O)SC", "S-methyl hydrazinecarbothioate"),
            ("NNC(=S)SC", "methyl hydrazinecarbodithioate"),
            // P-65.5.1: acid halides.
            ("NNC(=O)Cl", "hydrazinecarbonyl chloride"),
            ("NNC(=S)Cl", "hydrazinecarbothioyl chloride"),
            ("ClC(=O)NNC(=O)Cl", "hydrazine-1,2-dicarbonyl dichloride"),
            // P-68.3.1.2.4: semicarbazides, with the Blue Book examples
            // N-phenylhydrazinecarboxamide and
            // N,1-dimethylhydrazine-1-carboxamide; P-66.1.4.1.1 for the
            // chalcogen analogues.
            ("NNC(N)=O", "hydrazinecarboxamide"),
            ("NNC(=O)Nc1ccccc1", "N-phenylhydrazinecarboxamide"),
            ("CN(N)C(=O)NC", "N,1-dimethylhydrazine-1-carboxamide"),
            ("NC(=O)N(C)N", "1-methylhydrazine-1-carboxamide"),
            ("CNNC(=S)N", "2-methylhydrazine-1-carbothioamide"),
            ("NNC(N)=[Se]", "hydrazinecarboselenoamide"),
            ("NC(=S)NNC(N)=S", "hydrazine-1,2-dicarbothioamide"),
            (
                "NC(=O)NN(C(N)=O)C([N+](=O)[O-])[N+](=O)[O-]",
                "1-(dinitromethyl)hydrazine-1,2-dicarboxamide",
            ),
            // P-68.3.1.2.5: semicarbazones; benzylidene is the preferred
            // prefix (P-29.6.2.1).
            (
                "CCCC(CC)=NNC(=O)N(c1ccccc1)c1ccccc1",
                "2-(hexan-3-ylidene)-N,N-diphenylhydrazine-1-carboxamide",
            ),
            (
                "CC(C)=NNC(N)=S",
                "2-(propan-2-ylidene)hydrazine-1-carbothioamide",
            ),
            ("NC(=O)NN=Cc1ccccc1", "2-benzylidenehydrazine-1-carboxamide"),
            // P-68.3.1.2.6, P-66.3.4: carbohydrazides; the hydrazine with
            // more substituents is the parent.
            ("NNC(=O)NN", "hydrazinecarbohydrazide"),
            ("NNC(=S)NN", "hydrazinecarbothiohydrazide"),
            ("NNC(=O)NNC", "2-methylhydrazine-1-carbohydrazide"),
            ("CNNC(=O)N(C)N", "N,2-dimethylhydrazine-1-carbohydrazide"),
            ("CC(=O)NNC(=O)NN", "2-acetylhydrazine-1-carbohydrazide"),
            // A junior group on the same hydrazine is a prefix: carbamoyl
            // (P-68.3.1.2.4), carbonochloridoyl (P-65.2.1.5),
            // hydrazinecarbonyl (P-66.3.2.1).
            ("NC(=O)NNC(=O)O", "2-carbamoylhydrazine-1-carboxylic acid"),
            (
                "COC(=O)NNC(=O)O",
                "2-(methoxycarbonyl)hydrazine-1-carboxylic acid",
            ),
            (
                "ClC(=O)NNC(=O)O",
                "2-carbonochloridoylhydrazine-1-carboxylic acid",
            ),
            (
                "NNC(=O)NNC(=O)O",
                "2-(hydrazinecarbonyl)hydrazine-1-carboxylic acid",
            ),
            // P-44.1.1, P-44.1.2: with equal numbers of the principal group,
            // nitrogen makes hydrazine senior to a carbon chain or ring, as
            // in 2-(2-carboxyethyl)disilane-1-carboxylic acid (P-44.1.2.1).
            (
                "COC(=O)CNNC(=O)O",
                "2-(2-methoxy-2-oxoethyl)hydrazine-1-carboxylic acid",
            ),
            (
                "O=C(Cl)CNNC(=O)O",
                "2-(2-chloro-2-oxoethyl)hydrazine-1-carboxylic acid",
            ),
            (
                "NC(=O)CCNNC(=O)O",
                "2-(3-amino-3-oxopropyl)hydrazine-1-carboxylic acid",
            ),
            (
                "OC(=O)c1ccccc1NNC(=O)O",
                "2-(2-carboxyphenyl)hydrazine-1-carboxylic acid",
            ),
            (
                "NC(=O)CCNNC(N)=O",
                "2-(3-amino-3-oxopropyl)hydrazine-1-carboxamide",
            ),
            (
                "O=C(O)NNc1ccccn1",
                "2-(pyridin-2-yl)hydrazine-1-carboxylic acid",
            ),
            // P-41, class 11: amides rank as their acids, carboxylic before
            // sulfonic.
            (
                "NS(=O)(=O)CCNNC(N)=O",
                "2-(2-sulfamoylethyl)hydrazine-1-carboxamide",
            ),
        ] {
            let name = super::super::name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, NameStatus::Preferred, "{smiles}");
            let parsed = crate::test_support::parse_strict(&name.text).unwrap();
            assert_eq!(
                crate::test_support::canonical_smiles(smiles).unwrap(),
                crate::test_support::canonical_smiles(&parsed).unwrap(),
                "{expected}",
            );
        }
        // P-72.2.2.2.1.1: the anion replaces 'ic acid' by 'ate'. The charge
        // of a mixed chalcogen anion is not placed, and an acid beside an
        // anion is left to the general namer.
        for (smiles, expected) in [
            (
                "CNNC(=O)[O-].[Na+]",
                "sodium 2-methylhydrazine-1-carboxylate",
            ),
            (
                "NN(C)C(=O)[O-].[Na+]",
                "sodium 1-methylhydrazine-1-carboxylate",
            ),
            (
                "CNNC(=S)[S-].[Na+]",
                "sodium 2-methylhydrazine-1-carbodithioate",
            ),
            (
                "[O-]C(=O)NNC(=O)[O-].[Na+].[Na+]",
                "disodium hydrazine-1,2-dicarboxylate",
            ),
            (
                "OC(=O)CNNC(=O)[O-].[Na+]",
                "sodium 2-(carboxymethyl)hydrazine-1-carboxylate",
            ),
        ] {
            let name = super::super::name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, NameStatus::Preferred, "{smiles}");
            let parsed = crate::test_support::parse_strict(&name.text).unwrap();
            assert_eq!(
                crate::test_support::canonical_smiles(smiles).unwrap(),
                crate::test_support::canonical_smiles(&parsed).unwrap(),
                "{expected}",
            );
        }
        for smiles in ["CNNC(=O)[S-]", "[O-]C(=O)NNC(=O)O", "[O-]c1ccccc1NNC(=O)O"] {
            let mol = read_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(
                hydrazine_carbonyl_parent_name(&mol, &all_atoms(&mol)).unwrap(),
                None,
                "{smiles}"
            );
        }
        // A suffix needs a terminal =X and a terminal halogen.
        for smiles in ["NNC(=[S](C)C)O", "NNC(=S(=O)=O)O", "NNC(=O)I(C)C"] {
            let mol = read_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(
                hydrazine_carbonyl_parent_name(&mol, &all_atoms(&mol)).unwrap(),
                None,
                "{smiles}"
            );
        }
        // Left to the general namer: a parent with more of the principal
        // group (P-44.1.1), a nitrogen ring with as many (P-44.1.2.2), two
        // hydrazines with one each, a senior class elsewhere (P-41), and
        // acylhydrazides (P-66.3.1.1).
        for smiles in [
            "OC(=O)CC(C(=O)O)NNC(=O)O",
            "OC(=O)c1ccncc1NNC(=O)O",
            "OC(=O)NNCCNNC(=O)O",
            "OC(=O)CCNNC(N)=O",
            "NNC=O",
            "CC(=O)NNC(C)=O",
        ] {
            let mol = read_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(
                hydrazine_carbonyl_parent_name(&mol, &all_atoms(&mol)).unwrap(),
                None,
                "{smiles}"
            );
        }
    }

    #[test]
    fn cyanohydrazide_is_the_hydrazide_of_cyanic_acid() {
        // Issue #113. P-66.3.1.2.1, P-66.3.5.1: NC-NH-NH2 is cyanohydrazide,
        // substituted like an amide with the locants N and N' (P-66.3.3.1).
        for (smiles, expected) in [
            ("NNC#N", "cyanohydrazide"),
            ("CNNC#N", "N'-methylcyanohydrazide"),
            ("CN(N)C#N", "N-methylcyanohydrazide"),
            ("CNN(C)C#N", "N,N'-dimethylcyanohydrazide"),
            ("CN(C)NC#N", "N',N'-dimethylcyanohydrazide"),
            // P-14.3.4.5: no locants when every position is substituted
            // alike.
            ("CN(C)N(C)C#N", "trimethylcyanohydrazide"),
            ("N#CN(C)Nc1ccccc1", "N-methyl-N'-phenylcyanohydrazide"),
            ("N#CNN=C(C)C", "N'-(propan-2-ylidene)cyanohydrazide"),
            ("N#CNN=Cc1ccccc1", "N'-benzylidenecyanohydrazide"),
            // P-41: a hydrazide is senior to an alcohol and a nitrile.
            ("N#CNNCCO", "N'-(2-hydroxyethyl)cyanohydrazide"),
            ("N#CNNCC#N", "N'-(cyanomethyl)cyanohydrazide"),
        ] {
            let name = super::super::name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, NameStatus::Preferred, "{smiles}");
            let parsed = crate::test_support::parse_strict(&name.text).unwrap();
            assert_eq!(
                crate::test_support::canonical_smiles(smiles).unwrap(),
                crate::test_support::canonical_smiles(&parsed).unwrap(),
                "{expected}",
            );
        }
        // Left to the general namer and not labelled preferred: a senior or
        // equal class (P-41), also a carbonic, nitric, phosphorus or sulfur
        // acid or ester, a second cyano group on the hydrazine, a heteroatom
        // substituent, and a configured hydrazone.
        for smiles in [
            "N#CNNCOC(=O)O",
            "N#CNNCOC(=O)OC",
            "N#CNNCON(=O)=O",
            "N#CNNCP(=O)(O)O",
            "N#CNNCP(=O)(C)O",
            "N#CNNCS(=O)(=O)O",
            "N#CNNCC(=O)O",
            "N#CNNC(N)=O",
            "N#CNNC(C)=O",
            "N#CNNC#N",
            "N#CNNCl",
            "N#CN/N=C/c1ccccc1",
            // methanethioyl would not carry the N=C bond.
            "N#CNN=C=S",
        ] {
            let mol = read_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(
                cyanohydrazide_name(&mol, &all_atoms(&mol)).unwrap(),
                None,
                "{smiles}"
            );
            let name = super::super::name_smiles(smiles, &Default::default()).unwrap();
            assert_ne!(
                name.status,
                NameStatus::Preferred,
                "{smiles}: {}",
                name.text
            );
            let parsed = crate::test_support::parse_strict(&name.text).unwrap();
            assert_eq!(
                crate::test_support::canonical_smiles(smiles).unwrap(),
                crate::test_support::canonical_smiles(&parsed).unwrap(),
                "{}",
                name.text
            );
        }
        // A ring nitrogen takes the suffix carbonitrile (P-66.5.1.1.3).
        assert_eq!(
            super::super::name_smiles("N#CN1CCCN1", &Default::default())
                .unwrap()
                .text,
            "pyrazolidine-1-carbonitrile"
        );
    }

    #[test]
    fn nitrogen_acids_form_mixed_anhydrides() {
        // Issue #113. P-65.7.2: the acids of a mixed anhydride are cited in
        // alphabetical order; hydrazinecarboxylic acid (P-68.3.1.2.1) and
        // carbamic acid are monobasic acids like any other.
        for (smiles, expected) in [
            ("NNC(=O)OC(C)=O", "acetic hydrazinecarboxylic anhydride"),
            (
                "CNNC(=O)OC(C)=O",
                "acetic 2-methylhydrazine-1-carboxylic anhydride",
            ),
            (
                "CN(N)C(=O)OC(=O)CC",
                "1-methylhydrazine-1-carboxylic propanoic anhydride",
            ),
            (
                "NNC(=O)OC(=O)c1ccccc1",
                "benzoic hydrazinecarboxylic anhydride",
            ),
            ("NNC(=O)OC(=O)NN", "hydrazinecarboxylic anhydride"),
            ("NNC(=O)OC(=O)N", "carbamic hydrazinecarboxylic anhydride"),
            ("NC(=O)OC(C)=O", "acetic carbamic anhydride"),
            ("CNC(=O)OC(C)=O", "acetic methylcarbamic anhydride"),
            ("c1ccccc1NC(=O)OC(C)=O", "acetic phenylcarbamic anhydride"),
        ] {
            let name = super::super::name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, NameStatus::Preferred, "{smiles}");
            let parsed = crate::test_support::parse_strict(&name.text).unwrap();
            assert_eq!(
                crate::test_support::canonical_smiles(smiles).unwrap(),
                crate::test_support::canonical_smiles(&parsed).unwrap(),
                "{expected}",
            );
        }
        // Not named as anhydrides, and not labelled preferred: a thioacyl
        // half and a chalcogen or peroxy bridge. Nor is a dianhydride with
        // carbamic ends, whose di/bis choice is not made (P-16.3.3, P-16.3.5).
        for smiles in [
            "NNC(=S)OC(C)=O",
            "NNC(=O)SC(C)=O",
            "NNC(=O)OOC(C)=O",
            "NC(=O)OC(=O)CCC(=O)OC(N)=O",
        ] {
            let name = super::super::name_smiles(smiles, &Default::default()).unwrap();
            assert_ne!(
                name.status,
                NameStatus::Preferred,
                "{smiles}: {}",
                name.text
            );
        }
    }

    #[test]
    fn discarded_shortcut_candidates_restore_the_previous_status() {
        let mol = read_smiles("CC[S@](C)=O", &Default::default()).unwrap();
        let atoms = all_atoms(&mol);
        let previous = NameStatus::Unverified(NameReason::UnauditedParent);
        mol.mark_name(previous);
        let declined: Result<Option<String>, String> = mol.candidate(|| {
            assert_eq!(
                sulfoxide_parent_name(&mol, &atoms)?.as_deref(),
                Some("(S)-ethyl methyl sulfoxide")
            );
            assert_eq!(
                mol.name_status(),
                NameStatus::NotPreferred(NameReason::FunctionalClass)
            );
            Ok(None)
        });
        assert_eq!(declined.unwrap(), None);
        assert_eq!(mol.name_status(), previous);

        let failed: Result<Option<String>, String> = mol.candidate(|| {
            sulfoxide_parent_name(&mol, &atoms)?;
            Err("discarded shortcut".into())
        });
        assert_eq!(failed.unwrap_err(), "discarded shortcut");
        assert_eq!(mol.name_status(), previous);
    }

    #[test]
    fn functional_class_splits_follow_only_selected_component_names() {
        use super::super::{
            component_namer::name_component_with_citation, phane::has_linear_skeleton,
        };
        let mol = read_smiles(
            "O=C(Oc1ccc(Oc2ccccc2)cc1)Oc1ccc(Oc2ccccc2)cc1",
            &Default::default(),
        )
        .unwrap();
        let atoms = all_atoms(&mol);
        let expected = "bis(4-phenoxyphenyl) carbonate";
        let trial = structural_replacement_parent_result(&mol, &atoms)
            .unwrap()
            .unwrap();
        assert_eq!(trial.name, expected);
        assert!(has_linear_skeleton(&mol, &BTreeSet::new()));
        assert!(!has_linear_skeleton(&mol, &trial.component_splits));

        for (smiles, expected, split_count) in [
            (
                "O=C(OCc1ccc(Oc2ccccc2)cc1)c1ccc(Oc2ccccc2)cc1",
                "(4-phenoxyphenyl)methyl 4-phenoxybenzoate",
                1,
            ),
            (
                "O=C(Oc1ccc(Oc2ccccc2)cc1)Nc1ccc(Oc2ccccc2)cc1",
                "4-phenoxyphenyl (4-phenoxyphenyl)carbamate",
                1,
            ),
            (
                "O=C(OC(=O)c1ccc(Oc2ccccc2)cc1)c1ccc(Oc2ccccc2)cc1",
                "4-phenoxybenzoic anhydride",
                2,
            ),
            (
                "O=C(Oc1ccc(Oc2ccccc2)cc1)Oc1ccc(Oc2ccccc2)cc1",
                "bis(4-phenoxyphenyl) carbonate",
                2,
            ),
            (
                "CP(=O)(Oc1ccc(Oc2ccccc2)cc1)Oc1ccc(Oc2ccccc2)cc1",
                "bis(4-phenoxyphenyl) methylphosphonate",
                2,
            ),
        ] {
            let mol = read_smiles(smiles, &Default::default()).unwrap();
            let atoms = all_atoms(&mol);
            for fail in [false, true] {
                let discarded: Result<Option<()>, String> = mol.candidate(|| {
                    let name = name_component_with_citation(&mol, &atoms, None, None)?;
                    assert_eq!(name.text, expected);
                    assert_eq!(name.derivation.component_splits().len(), split_count);
                    assert!(!has_linear_skeleton(
                        &mol,
                        name.derivation.component_splits()
                    ));
                    if fail {
                        Err("discarded functional-class name".into())
                    } else {
                        Ok(None)
                    }
                });
                assert_eq!(discarded.is_err(), fail);
                assert!(has_linear_skeleton(&mol, &BTreeSet::new()));
            }

            let fragment = mol.subgraph(&atoms);
            let fragment_name =
                name_component_with_citation(&fragment, &atoms, None, None).unwrap();
            assert_eq!(fragment_name.text, expected);
            let name = name_component_with_citation(&mol, &atoms, None, None).unwrap();
            assert_eq!(name.text, expected);
            assert_eq!(
                name.derivation.component_splits(),
                fragment_name.derivation.component_splits()
            );
            assert!(!has_linear_skeleton(
                &mol,
                name.derivation.component_splits()
            ));
            assert!(has_linear_skeleton(&mol.clone(), &BTreeSet::new()));
        }
    }

    #[test]
    // Provenance: fn tests::single_atom_hydrides_require_complete_hydrogen_count | adapted
    // Source: ReShiki | src/naming/openclatura/special_cases.rs::tests::single_atom_hydrides_require_complete_hydrogen_count
    // Source: openclatura | tests/exploratory/test_group14_ether_ligand_parent_selection.py::test_a_genuine_principal_group_elsewhere_still_blocks_the_hydride_parent
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn single_atom_hydrides_require_complete_hydrogen_count() {
        let mut mol = Molecule::new();
        let atom = Atom::new(0, "Si");
        mol.add_atom(atom).unwrap();
        assert_eq!(single_atom_component_name(&mol, &set([0])), None);
        mol.atoms.get_mut(&0).unwrap().total_h_count = 4;
        assert_eq!(
            single_atom_component_name(&mol, &set([0])).as_deref(),
            Some("silane")
        );
    }

    #[test]
    // Provenance: fn tests::graph_shortcuts_preserve_ligand_and_oxygen_roles | adapted
    // Source: ReShiki | src/naming/openclatura/special_cases.rs::tests::graph_shortcuts_preserve_ligand_and_oxygen_roles
    // Source: openclatura | tests/unit/test_analysis.py::test_oxoacid_ester_suffixes_are_role_template_backed
    // Source: openclatura | tests/roundtrip/test_substituents.py::test_roundtrip
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn graph_shortcuts_preserve_ligand_and_oxygen_roles() {
        assert_eq!(special("CS(C)=O"), None);
        assert_eq!(
            special("CC[S@@](C)=O").as_deref(),
            Some("(R)-ethyl methyl sulfoxide")
        );
        assert_eq!(
            special("CP(=O)(O)O").as_deref(),
            Some("methylphosphonic acid")
        );
        assert_eq!(
            special("CO[P](=O)(O)O").as_deref(),
            Some("methyl dihydrogen phosphate")
        );
        assert_eq!(special("CB(O)O").as_deref(), Some("methylboronic acid"));
    }

    #[test]
    fn chalcogen_onium_hydrides_use_preferred_systematic_parents() {
        for (smiles, expected) in [
            ("C[S+](C)C", "trimethylsulfanium"),
            ("C[SH+]C", "dimethylsulfanium"),
            ("C[SH2+]", "methylsulfanium"),
            ("C[S+](C)c1ccccc1", "dimethyl(phenyl)sulfanium"),
            ("CC[S+](CC)CC", "triethylsulfanium"),
            ("C[S+](C)CC", "ethyldi(methyl)sulfanium"),
            ("CC[S+](C)CCC", "ethyl(methyl)(propyl)sulfanium"),
            ("C[S+](C)CCO", "(2-hydroxyethyl)di(methyl)sulfanium"),
            ("C[S+](C)CC(=O)O", "(carboxymethyl)di(methyl)sulfanium"),
            ("C[Se+](C)C", "trimethylselanium"),
            ("CCC[SeH2+]", "propylselanium"),
            ("C[Te+](C)C", "trimethyltellanium"),
            ("C[O+](C)C", "trimethyloxidanium"),
            ("C[OH+]C", "dimethyloxidanium"),
            ("C[OH2+]", "methyloxidanium"),
            ("CC[O+](CC)CC", "triethyloxidanium"),
            ("C[O+](C)CC", "ethyldi(methyl)oxidanium"),
            ("C[O+](C)c1ccccc1", "dimethyl(phenyl)oxidanium"),
            ("C[O+](C)CC(=O)O", "(carboxymethyl)di(methyl)oxidanium"),
            ("C[O+](C)C(C)=O", "acetyldi(methyl)oxidanium"),
            // Issue #18. P-73.1.1.1, P-73.1.2.1: ylidene and ylidyne ligands
            // (ethylideneoxidanium, ethyl(propan-2-ylidene)oxidanium,
            // ethylidyneoxidanium are the Blue Book PINs).
            ("C=[O+]C", "methyl(methylidene)oxidanium"),
            ("CC=[O+]C", "ethylidene(methyl)oxidanium"),
            ("CC=[OH+]", "ethylideneoxidanium"),
            ("CC(C)=[O+]CC", "ethyl(propan-2-ylidene)oxidanium"),
            ("CC#[O+]", "ethylidyneoxidanium"),
            ("C1CCCCC1=[O+]C", "cyclohexylidene(methyl)oxidanium"),
            ("C=[S+]C", "methyl(methylidene)sulfanium"),
            ("CC=[Se+]C", "ethylidene(methyl)selanium"),
            ("NCC=[OH+]", "(2-aminoethylidene)oxidanium"),
            // Issue #18. P-41, P-73.1.2: the cation outranks a neutral amine or
            // other heteroatom centre, which becomes a prefix.
            (
                "C[O+](C)CCN(C)C",
                "[2-(dimethylamino)ethyl]di(methyl)oxidanium",
            ),
            ("C=[O+]CCN", "(2-aminoethyl)(methylidene)oxidanium"),
            ("C[O+](C)C1CCNCC1", "dimethyl(piperidin-4-yl)oxidanium"),
            (
                "C[S+](C)CCN(C)C",
                "[2-(dimethylamino)ethyl]di(methyl)sulfanium",
            ),
        ] {
            let mol = read_smiles(smiles, &Default::default()).unwrap();
            let result = simple_central_parent_hydride_name(&mol, &all_atoms(&mol))
                .unwrap()
                .unwrap();
            // Enclosing marks are normalized to brackets after this shortcut.
            assert_eq!(
                result.name,
                expected.replace('[', "(").replace(']', ")"),
                "{smiles}"
            );
            assert!(result.certified, "{smiles}");
            let name = super::super::name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, NameStatus::Preferred, "{smiles}");
            let parsed = nazuchem_n2s::Parser::new()
                .unwrap()
                .parse(&name.text, &nazuchem_n2s::ParseOptions::strict());
            let reconstructed = parsed.structure.unwrap().semantic_cxsmiles().unwrap();
            crate::test_support::verify_identity(smiles, &reconstructed)
                .unwrap_or_else(|error| panic!("{smiles} -> {expected}: {error}"));
        }
    }

    #[test]
    fn chalcogen_onium_shortcut_declines_rings_and_other_ionic_centers() {
        for smiles in [
            "C[S+]1CCCC1",
            "C[O+]1CCCC1",
            // The ligand names cannot carry E/Z of the C=O+ bond.
            "C/C=[O+]/C",
            // P-73.1.2.2: protonated (thio)urea is named from uronium.
            "[OH+]=C(N)N",
            "[SH+]=C(N)N",
            // Methanethioyl would cite a single bond to the centre.
            "S=C=[OH+]",
            "C[O+](C)CC[S+](C)C",
            "C[S+](C)[CH2-]",
            "C[S+](C)C[O-]",
            "C[S+](C)CCCCCC[N+](C)(C)C",
        ] {
            let mol = read_smiles(smiles, &Default::default()).unwrap();
            assert!(
                simple_central_parent_hydride_name(&mol, &all_atoms(&mol))
                    .unwrap()
                    .is_none(),
                "{smiles}"
            );
        }
        for radical_center in [0, 1] {
            let mut mol = read_smiles("C[S+](C)C", &Default::default()).unwrap();
            mol.atoms
                .get_mut(&radical_center)
                .unwrap()
                .radical_electrons = 1;
            assert!(
                simple_central_parent_hydride_name(&mol, &all_atoms(&mol))
                    .unwrap()
                    .is_none()
            );
        }
        for (smiles, expected) in [
            (
                "C[S+](C)CCCCCC[N+](C)(C)C",
                "6-(dimethylsulfaniumyl)-N,N,N-trimethylhexan-1-aminium",
            ),
            ("[CH2-][S+](C)C", "(dimethylsulfaniumyl)methanide"),
            ("C[N+](C)(C)C", "N,N,N-trimethylmethanaminium"),
        ] {
            assert_eq!(
                super::super::name_smiles(smiles, &Default::default())
                    .unwrap()
                    .text,
                expected,
                "{smiles}"
            );
        }
    }

    #[test]
    fn mononuclear_hydrides_enclose_later_repeated_simple_ligands() {
        for (smiles, expected) in [
            ("CC[Si](C)(C)C", "ethyltri(methyl)silane"),
            ("CC[P+](C)(C)C", "ethyltri(methyl)phosphanium"),
            ("C[Si](C)(C)c1ccccc1", "trimethyl(phenyl)silane"),
        ] {
            assert_eq!(special(smiles).as_deref(), Some(expected), "{smiles}");
        }
    }

    #[test]
    // Provenance: fn tests::shortcut_certificates_require_audited_ligands | adapted
    // Source: ReShiki | src/naming/openclatura/special_cases.rs::tests::shortcut_certificates_require_audited_ligands
    // Source: openclatura | tests/unit/test_analysis.py::test_group_13_14_central_atoms_have_a_parent_hydride
    // Source: openclatura | tests/roundtrip/test_replacement_parents.py::test_roundtrip
    // Source: openclatura | tests/unit/test_analysis.py::test_nitrogen_chain_roles_make_cyclic_hydrazines_prefixes
    // Source: openclatura | tests/roundtrip/test_oxyacids.py::test_roundtrip
    // Source: openclatura | tests/roundtrip/test_role_models.py::test_roundtrip
    // Source: openclatura | tests/roundtrip/test_principal_group_breadth.py::test_roundtrip
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn shortcut_certificates_require_audited_ligands() {
        for (smiles, expected, certified) in [
            ("CP(C)(C)=O", "trimethyl-lambda5-phosphanone", true),
            ("C[P+](C)(C)[O-]", "trimethyl-lambda5-phosphanone", true),
            ("C[As+](C)(C)[O-]", "trimethyl-lambda5-arsanone", true),
            ("C[Sb+](C)(C)[O-]", "trimethyl-lambda5-stibanone", true),
            ("CP(C)(C)=S", "trimethyl-lambda5-phosphanethione", true),
            ("CP(C)(C)=[Se]", "trimethyl-lambda5-phosphaneselone", true),
            ("CP(C)(C)=[Te]", "trimethyl-lambda5-phosphanetellone", true),
            ("CC[Sn](CC)(CC)CC", "tetraethylstannane", true),
            ("CP(C)Cl", "chlorodi(methyl)phosphane", false),
            ("COP(C)C", "methoxydi(methyl)phosphane", false),
            ("C[Al](C)C", "trimethylalumane", false),
            ("CCB(O)O", "ethylboronic acid", true),
            ("CCO[N+](=O)[O-]", "ethyl nitrate", true),
            (
                "CN(C)C(=O)OP(=O)(O)O",
                "dimethylcarbamoyl dihydrogen phosphate",
                false,
            ),
            (
                "CC(=S)OP(=O)(OC)OC",
                "ethanethioyl dimethyl phosphate",
                false,
            ),
            (
                "CC(=N)OP(=O)(OC)OC",
                "1-iminoethyl dimethyl phosphate",
                false,
            ),
            ("NN", "hydrazine", false),
        ] {
            let mol = read_smiles(smiles, &Default::default()).unwrap();
            let result = structural_replacement_parent_result(&mol, &all_atoms(&mol))
                .unwrap()
                .unwrap();
            assert_eq!(result.name, expected, "{smiles}");
            assert_eq!(result.certified, certified, "{smiles}");
        }
    }

    #[test]
    fn alkyl_shortcuts_decline_unrendered_ligand_stereo() {
        for smiles in [
            "C[C@H](CC)P(=O)(O)O",
            "C[C@@H](CC)P(=O)(O)O",
            "C[C@H](CC)B(O)O",
            "C[C@@H](CC)B(O)O",
            "C[C@H](CC)S(C)=O",
            "CCC[C@H](C)[S@@](C)=O",
        ] {
            assert_eq!(
                super::super::name_smiles(smiles, &Default::default()).unwrap_err(),
                "openclatura: shortcut cannot preserve configured stereochemistry",
                "{smiles}"
            );
        }
        let left = super::super::name_smiles("C[S@](=O)CC", &Default::default())
            .unwrap()
            .text;
        let right = super::super::name_smiles("C[S@@](=O)CC", &Default::default())
            .unwrap()
            .text;
        assert_ne!(left, right);
        assert!(left.starts_with("(R)-") || left.starts_with("(S)-"));
        assert!(right.starts_with("(R)-") || right.starts_with("(S)-"));
    }

    #[test]
    fn manual_azine_shortcuts_decline_unrendered_atom_and_bond_stereo() {
        for smiles in [
            "C=NN=C[C@H]1CCOC1",
            "C/C=C/C=NN=C",
            "CCC/C=N/N=C/C",
            "C/C=N/N=CCC",
            "C/C=N/N=C\\C",
        ] {
            assert!(
                super::super::name_smiles(smiles, &Default::default()).is_err(),
                "{smiles}"
            );
        }
    }

    #[test]
    fn charged_alkyl_ligands_are_not_rendered_as_neutral() {
        for smiles in ["C[CH-]P(=O)(O)O", "C[CH-]B(O)O", "C[CH-]S(C)=O"] {
            assert!(
                super::super::name_smiles(smiles, &Default::default()).is_err(),
                "{smiles}"
            );
        }
    }

    #[test]
    fn public_naming_bounds_shortcut_alkyl_path_search() {
        use std::time::{Duration, Instant};

        let mut tree = "C".to_owned();
        for _ in 0..7 {
            tree = format!("C({tree}){tree}");
        }
        let started = Instant::now();
        assert_eq!(
            super::super::name_smiles(&format!("P(=O)(O)(O){tree}"), &Default::default())
                .map(|name| name.text),
            Err("openclatura: naming complexity limit exceeded".into())
        );
        assert!(started.elapsed() < Duration::from_secs(10));
        assert_eq!(
            super::super::name_smiles("CCC(CC)P(=O)(O)O", &Default::default())
                .map(|name| name.text)
                .as_deref(),
            Ok("pentan-3-ylphosphonic acid")
        );
    }

    #[test]
    fn phosphorus_functional_parents_keep_their_own_ligand_grouping() {
        for (smiles, expected) in [
            ("CCC(CC)P(=O)(OC)OC", "dimethyl pentan-3-ylphosphonate"),
            ("CCC(CC)P(=O)(O)O", "pentan-3-ylphosphonic acid"),
            ("CCC(CC)[PH](=O)O", "(pentan-3-yl)phosphinic acid"),
        ] {
            assert_eq!(special(smiles).as_deref(), Some(expected), "{smiles}");
        }
    }

    #[test]
    fn azine_shortcut_uses_symmetric_carbonyl_and_ring_side_roles() {
        for (smiles, expected) in [
            ("C=N-N=C", "dimethylidenehydrazine"),
            ("CC(C)=NN=C(C)C", "di(propan-2-ylidene)hydrazine"),
            ("CCC(C)=NN=CC", "butan-2-one ethan-1-ylidenehydrazone"),
            ("c1ccccc1C=NN=C", "benzaldehyde methylidenehydrazone"),
            (
                "c1ccncc1C=NN=C",
                "pyridine-3-carbaldehyde methylidenehydrazone",
            ),
            ("C1CCCCC1=NN=C", "formaldehyde cyclohexylidenehydrazone"),
        ] {
            let mol = read_smiles(smiles, &Default::default()).unwrap();
            let actual =
                super::super::component_namer::name_component(&mol, &all_atoms(&mol)).unwrap();
            assert_eq!(
                super::super::nest_enclosing_marks(&actual),
                expected,
                "{smiles}"
            );
        }
    }

    #[test]
    fn azine_unsaturation_requires_order_independent_locants() {
        for smiles in [
            "CCCC=NN=CCC=C",
            "CCCC=NN=CCC#C",
            "C=CCC=NN=CCCC",
            "C#CCC=NN=CCCC",
        ] {
            assert_eq!(
                super::super::name_smiles(smiles, &Default::default()).unwrap_err(),
                "openclatura: azine unsaturation depends on upstream set ordering",
                "{smiles}"
            );
        }
        for smiles in ["CCCC=NN=CC=CC", "CC=CC=NN=CCCC"] {
            assert_eq!(
                super::super::name_smiles(smiles, &Default::default())
                    .map(|name| name.text)
                    .as_deref(),
                Ok("butanal but-2-ene-1-ylidenehydrazone"),
                "{smiles}"
            );
        }
    }

    #[test]
    fn substituted_amidino_azine_precedes_a_ligand_anhydride() {
        let smiles = "CC(=O)OC(=O)N(/C(N)=N/N=C=O)C(C)C";
        let mol = read_smiles(smiles, &Default::default()).unwrap();
        let atoms = all_atoms(&mol);
        let expected = "1-((E)-(((acetyloxy)carbonyl)(propan-2-yl)amino)(amino)methylidene)-2-(oxomethylidene)hydrazine";
        assert_eq!(
            simple_azine_parent_name(&mol, &atoms).unwrap().as_deref(),
            Some(expected)
        );
        assert_eq!(
            super::super::component_namer::name_component(&mol, &atoms).unwrap(),
            expected
        );
    }

    #[test]
    // Provenance: fn tests::homonuclear_parent_keeps_chain_and_substituent_locants | adapted
    // Source: ReShiki | src/naming/openclatura/special_cases.rs::tests::homonuclear_parent_keeps_chain_and_substituent_locants
    // Source: openclatura | tests/unit/test_polyoxide_chains.py::test_the_other_heteroatom_chains_are_unchanged
    // Source: openclatura | tests/unit/test_chalcogen_functional_groups.py::test_peroxide_linkages_preserve_both_chalcogen_sites
    // Source: openclatura | tests/data/parser_xml_resources/functionalTerms.json::<fixture data>
    // Source: openclatura | tests/unit/test_analysis.py::test_nitrogen_chain_roles_make_cyclic_hydrazines_prefixes
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn homonuclear_parent_keeps_chain_and_substituent_locants() {
        assert_eq!(special("CSSC").as_deref(), Some("1,2-dimethyldisulfane"));
        for (smiles, expected) in [("SS", "disulfane"), ("SSS", "trisulfane")] {
            let mol = read_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(
                homonuclear_chain_parent_name(&mol, &all_atoms(&mol))
                    .unwrap()
                    .as_deref(),
                Some(expected)
            );
            assert_eq!(
                mol.name_status(),
                NameStatus::Unverified(NameReason::UnauditedPath),
                "{smiles}"
            );
        }
        assert_eq!(special("NN").as_deref(), Some("hydrazine"));
        assert_eq!(
            special("Cl[SiH2][SiH2]F").as_deref(),
            Some("1-chloro-2-fluorodisilane")
        );
    }

    #[test]
    // Provenance: fn tests::anhydride_halves_use_acid_citation_suffixes | adapted
    // Source: ReShiki | src/naming/openclatura/special_cases.rs::tests::anhydride_halves_use_acid_citation_suffixes
    // Source: openclatura | tests/unit/test_all_chalcogen_acyl_regressions.py::test_reported_name_regressions
    // Source: openclatura | tests/unit/test_chalcogen_functional_groups.py::test_anhydride_bridges_and_acyl_sites_are_independent
    // Source: openclatura | tests/unit/test_analysis.py::test_anhydride_shortcut_binds_acid_halves_and_core_separately
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn anhydride_halves_use_acid_citation_suffixes() {
        for (smiles, expected) in [
            ("CC(=O)OC(C)=O", "acetic anhydride"),
            ("CCC(=S)OC(=O)", "formic propanethioic anhydride"),
            ("CCC(=S)OSC(=S)CC", "propanethioic thioperoxyanhydride"),
            ("CC(=S)OC(C)=S", "ethanethioic anhydride"),
        ] {
            let mol = read_smiles(smiles, &Default::default()).unwrap();
            let actual =
                super::super::component_namer::name_component(&mol, &all_atoms(&mol)).unwrap();
            assert_eq!(actual, expected, "{smiles}");
        }
    }

    #[test]
    fn mixed_anhydride_acids_follow_alphanumerical_order() {
        // Issues #9 and #19. P-65.7.2: the acids are cited in alphanumerical
        // order (P-14.5).
        let preferred = NameStatus::Preferred;
        for (smiles, expected, status) in [
            // P-14.5: Roman letters first; locants take no part.
            (
                "CC(=O)OC(=O)C(Cl)C",
                "acetic 2-chloropropanoic anhydride",
                preferred,
            ),
            (
                "CCCC(=O)OC(=O)C(C)Cl",
                "butanoic 2-chloropropanoic anhydride",
                preferred,
            ),
            (
                "ClCC(=O)OC(=O)C",
                "acetic chloroacetic anhydride",
                preferred,
            ),
            // P-14.5: the letter of a compound locant takes no part.
            (
                "O=C(OC(C)=O)C1CC2(Br)CCCCC2CC1",
                "acetic 8a-bromodecahydronaphthalene-2-carboxylic anhydride",
                preferred,
            ),
            // P-14.5.4: 4 is lower than 4a.
            (
                "BrC1C(CC(C2CCCCC12)C(=O)OC(=O)C1C(CCC2(CCCCC12)Br)Cl)Cl",
                "4-bromo-3-chlorodecahydronaphthalene-1-carboxylic \
                 4a-bromo-2-chlorodecahydronaphthalene-1-carboxylic anhydride",
                preferred,
            ),
            // P-14.5: enclosing marks take no part ('phenoxy' before 'phenyls').
            (
                "O=C(OC(=O)c1ccc(Sc2ccccc2)cc1)c1ccc(Oc2ccccc2)cc1",
                "4-phenoxybenzoic 4-(phenylsulfanyl)benzoic anhydride",
                preferred,
            ),
            (
                "CC(=O)OC(=O)c1ccc(Oc2ccc(Oc3ccc(Oc4ccccc4)cc3)cc2)cc1",
                "acetic 4-[4-(4-phenoxyphenoxy)phenoxy]benzoic anhydride",
                NameStatus::NotPreferred(NameReason::LinearPhane),
            ),
            // P-14.5.4: equal letters; lowest locants at the first point of
            // difference, in order of citation.
            (
                "O=C(OC(=O)CCCCCCCCC(Cl)C)C(Cl)CCCCCCCCC",
                "2-chloroundecanoic 10-chloroundecanoic anhydride",
                preferred,
            ),
            (
                "O=C(OC(=O)C(Cl)C(Br)C)C(Br)C(Cl)C",
                "2-bromo-3-chlorobutanoic 3-bromo-2-chlorobutanoic anhydride",
                preferred,
            ),
            // P-14.5: locants of the primary substituents are compared before
            // those inside enclosing marks ({2,4} before {2,5}).
            (
                "O=C(OC(=O)c1c(C(C)CCC)ccc(S)c1)c1c(C(CC)CC)cc(S)cc1",
                "2-(pentan-3-yl)-4-sulfanylbenzoic 2-(pentan-2-yl)-5-sulfanylbenzoic anhydride",
                preferred,
            ),
            // P-14.5: Greek letters take no part; lambda^4 gives no Roman
            // letters ('phosphanyl' before 'sulfanyl').
            (
                "C[SH](C)c1ccc(C(=O)OC(=O)c2ccc(P(C)C)cc2)cc1",
                "4-(dimethylphosphanyl)benzoic 4-(dimethyl-lambda^4-sulfanyl)benzoic anhydride",
                preferred,
            ),
            // P-14.5: an italic element locant is lower than a numeral.
            (
                "CN(CCC(=O)OC(=O)CCNC(=O)c1ccc(C)cc1C)C(=O)c1ccc(C)cc1",
                "3-(N,4-dimethylbenzamido)propanoic 3-(2,4-dimethylbenzamido)propanoic anhydride",
                preferred,
            ),
            // P-14.5.3: s of s-indacene is italic, set off by a hyphen
            // (P-16.2.4.1(d)); 'chloroindacene' before 'chloroindene'.
            (
                "ClC=1C=C(C2=CC3=CC=CC3=CC12)C(=O)OC(=O)C=1CC2=CC=CC=C2C1Cl",
                "3-chloro-s-indacene-1-carboxylic 3-chloro-1H-indene-2-carboxylic anhydride",
                preferred,
            ),
            // P-14.5.3: fusion letters are italic and count only after the
            // Roman letters ('carbot' before 'carbox', 'cyclopentan' before
            // 'cyclopentap').
            (
                "BrC1=CC=C2C(=N1)C=C(O2)C(=O)OC(=S)C2=C(C=1C=NC=CC1O2)Br",
                "3-bromofuro[3,2-c]pyridine-2-carbothioic \
                 5-bromofuro[3,2-b]pyridine-2-carboxylic anhydride",
                preferred,
            ),
            (
                "C1(=CC2=C3C(=CC4=CC=CC5=CC=C1C3=C45)C=C2)C(=O)OC(=O)C2CCCC2",
                "cyclopentanecarboxylic cyclopenta[cd]pyrene-1-carboxylic anhydride",
                preferred,
            ),
            // P-14.5: the indicated hydrogen 1H is a locant ('iodo' before
            // 'pyrrole').
            (
                "ICC(=O)OC(=O)c1ccc[nH]1",
                "iodoacetic 1H-pyrrole-2-carboxylic anhydride",
                preferred,
            ),
            // P-14.5.1, P-14.5.3: tert is italic ('butyl' before 'chloro').
            (
                "CC(C)(C)c1ccc(cc1)C(=O)OC(=O)c1ccc(Cl)cc1",
                "4-(tert-butyl)benzoic 4-chlorobenzoic anhydride",
                preferred,
            ),
            // P-14.5: stereodescriptors take no part ('bromob' before 'bromoc').
            (
                "O=C(OC(=O)c1c(Br)cc([C@H](Br)Cl)cc1)c1cc(Br)c([C@H](F)Cl)cc1",
                "2-bromo-4-[(S)-bromo(chloro)methyl]benzoic \
                 3-bromo-4-[(R)-chloro(fluoro)methyl]benzoic anhydride",
                preferred,
            ),
            // P-14.6: equal alphanumerical characters; a letter is preferred
            // to an opening parenthesis, whichever half the input gives first.
            (
                "N[SiH](C)c1ccc(cc1)C(=O)OC(=O)c1ccc(cc1)[SiH2]CN",
                "4-[amino(methyl)silyl]benzoic 4-[(aminomethyl)silyl]benzoic anhydride",
                preferred,
            ),
            (
                "NC[SiH2]c1ccc(cc1)C(=O)OC(=O)c1ccc(cc1)[SiH](C)N",
                "4-[amino(methyl)silyl]benzoic 4-[(aminomethyl)silyl]benzoic anhydride",
                preferred,
            ),
            // P-14.6 with the stereodescriptor (R) and its hyphen excluded.
            (
                "N[P@](C)c1ccc(cc1)C(=O)OC(=O)c1ccc(cc1)[PH]CN",
                "4-[(R)-amino(methyl)phosphanyl]benzoic 4-[(aminomethyl)phosphanyl]benzoic anhydride",
                preferred,
            ),
        ] {
            let name = super::super::name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, status, "{smiles}");
            let parsed = crate::test_support::parse_strict(&name.text).unwrap();
            assert_eq!(
                crate::test_support::canonical_smiles(smiles).unwrap(),
                crate::test_support::canonical_smiles(&parsed).unwrap(),
                "{expected}"
            );
        }
    }

    #[test]
    fn sulfonic_and_sulfinic_anhydrides_are_named_as_anhydrides() {
        // Issue #93. P-65.7.1 and P-65.7.2 name anhydrides of sulfonic and
        // sulfinic acids like those of carboxylic acids, not as acyl esters.
        for (smiles, expected) in [
            // P-65.7.2 example, acids in alphanumerical order (P-14.5).
            (
                "O=C(OS(=O)(=O)c1ccc([N+](=O)[O-])cc1)CCl",
                "chloroacetic 4-nitrobenzene-1-sulfonic anhydride",
            ),
            ("CC(=O)OS(=O)(=O)C", "acetic methanesulfonic anhydride"),
            ("CS(=O)(=O)OC(C)=O", "acetic methanesulfonic anhydride"),
            // P-65.6.3.1.1: an ester of an alcohol stays an ester.
            ("COS(=O)(=O)C", "methyl methanesulfonate"),
            ("CS(=O)OC", "methyl methanesulfinate"),
            (
                "O=S(=O)(OC(=O)c1ccccc1)c1ccc(C)cc1",
                "benzoic 4-methylbenzene-1-sulfonic anhydride",
            ),
            // P-65.7.2 example.
            (
                "O=S(c1ccccc1)OS(=O)(=O)CC",
                "benzenesulfinic ethanesulfonic anhydride",
            ),
            // P-65.7.1 examples.
            (
                "O=S(=O)(c1ccccc1)OS(=O)(=O)c1ccccc1",
                "benzenesulfonic anhydride",
            ),
            ("ClCCS(=O)OS(=O)CCCl", "2-chloroethane-1-sulfinic anhydride"),
            // P-41: the acid is senior to the anhydride, which becomes the
            // enclosed compound prefix (acetyloxy)sulfonyl (P-16.5.1.1, P-65.3.2.3).
            (
                "OC(=O)CCS(=O)(=O)OC(C)=O",
                "3-[(acetyloxy)sulfonyl]propanoic acid",
            ),
            (
                "OC(=O)CCS(=O)(=O)OC(=O)CCl",
                "3-{[(2-chloroacetyl)oxy]sulfonyl}propanoic acid",
            ),
        ] {
            let name = super::super::name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, NameStatus::Preferred, "{smiles}");
            let parsed = crate::test_support::parse_strict(&name.text).unwrap();
            assert_eq!(
                crate::test_support::canonical_smiles(smiles).unwrap(),
                crate::test_support::canonical_smiles(&parsed).unwrap(),
                "{expected}"
            );
        }
    }

    #[test]
    // Provenance: fn tests::anhydride_certificate_counts_linkages_within_the_component | adapted
    // Source: ReShiki | src/naming/openclatura/special_cases.rs::tests::anhydride_certificate_counts_linkages_within_the_component
    // Source: openclatura | tests/unit/test_analysis.py::test_anhydride_shortcut_binds_acid_halves_and_core_separately
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn anhydride_certificate_counts_linkages_within_the_component() {
        use super::super::principal_groups::{component_groups, component_principal_key};

        let mut mol = Molecule::new();
        let mut cases = Vec::new();
        for (smiles, expected, certified) in [
            ("CC(=O)OC(C)=O", "acetic anhydride", true),
            // P-65.7.6.1 example (issue #110).
            (
                "CC(=O)OC(=O)CCC(=O)OC(C)=O",
                "diacetic butanedioic dianhydride",
                true,
            ),
            // P-65.7.6.4.2 example: unlike linkages are not a dianhydride.
            (
                "CC(=O)OC(=O)CCCC(=O)SC(C)=O",
                "acetic 5-(acetylsulfanyl)-5-oxopentanoic anhydride",
                false,
            ),
        ] {
            let fragment = read_smiles(smiles, &Default::default()).unwrap();
            let atom_offset = mol.atoms.keys().max().map_or(0, |idx| idx + 1);
            let bond_offset = mol.bonds.keys().max().map_or(0, |idx| idx + 1);
            let mut atoms = AtomSet::new();
            for mut atom in fragment.atoms.into_values() {
                atom.idx += atom_offset;
                atoms.insert(atom.idx);
                mol.add_atom(atom).unwrap();
            }
            for mut bond in fragment.bonds.into_values() {
                bond.idx += bond_offset;
                bond.u += atom_offset;
                bond.v += atom_offset;
                mol.add_bond(bond).unwrap();
            }
            cases.push((atoms, expected, certified));
        }
        assert_eq!(
            super::super::graph_io::get_connected_components(&mol).len(),
            3
        );
        for (atoms, expected, certified) in cases {
            let groups = component_groups(&mol, &atoms);
            let key = component_principal_key(&mol, &groups, false);
            let result = try_name_anhydride_component(&mol, &groups, key.as_deref())
                .unwrap()
                .unwrap();
            assert_eq!(result.name, expected);
            assert_eq!(result.certified, certified, "{expected}");
        }
    }

    #[test]
    // Provenance: fn tests::peroxol_suffix_defers_only_for_a_senior_heteroatom_chain | adapted
    // Source: ReShiki | src/naming/openclatura/special_cases.rs::tests::peroxol_suffix_defers_only_for_a_senior_heteroatom_chain
    // Source: openclatura | tests/unit/test_chalcogen_functional_groups.py::test_peroxol_table_is_order_sensitive
    // Source: openclatura | tests/unit/test_polyoxide_chains.py::test_oxygen_chains_do_not_displace_a_senior_heteroatom_chain
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn peroxol_suffix_defers_only_for_a_senior_heteroatom_chain() {
        for (smiles, expected, chain_wins) in [
            ("COO", "methaneperoxol", false),
            ("CS[SH]", "methanedithioperoxol", false),
            ("C[Se][SeH]", "methanediselenoperoxol", false),
            ("C[Te][TeH]", "methaneditelluroperoxol", false),
            ("SSCCOC", "2-methoxyethane-1-dithioperoxol", false),
            (
                "O=P(O)(O)c1ccccc1SS",
                "1-(2-phosphonophenyl)disulfane",
                true,
            ),
        ] {
            let mol = read_smiles(smiles, &Default::default()).unwrap();
            let atoms = all_atoms(&mol);
            let groups = super::super::principal_groups::component_groups(&mol, &atoms);
            let key = super::super::principal_groups::component_principal_key(&mol, &groups, false);
            let result = structural_replacement_parent_result(&mol, &atoms).unwrap();
            assert_eq!(
                blocking_suffix_defers_to_chain_parent(
                    &mol,
                    &atoms,
                    &groups,
                    key.as_deref(),
                    result.as_ref()
                ),
                chain_wins,
                "{smiles}"
            );
            let actual = super::super::component_namer::name_component(&mol, &atoms).unwrap();
            assert_eq!(actual, expected, "{smiles}");
        }
    }
}
