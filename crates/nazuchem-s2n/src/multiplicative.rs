// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | crates/nazuchem-s2n/src/multiplicative.rs | issue #4
// Copyright (c) 2026 Ameyanagi and Nazuchem contributors

//! Multiplicative names for identical units on a simple or concatenated
//! linking group (P-15.3.1, P-51.3.1).
use super::{
    chains::{self, RingSystem, get_cyclic_atoms},
    component_group_rules, component_namer,
    graph_kernel::{EdgeSet, connected_components, normalize_edge},
    molecule::{AtomSet, Molecule},
    parent_selection::{self, ParentSelection},
    perception, principal_groups,
    rules::{multipliers, stems},
    subgraph_tools::subgraph_component,
};
use nazuchem_ir::{NameReason, NameStatus};
use std::collections::BTreeMap;

/// Chains and ring systems that parents are chosen from.
type Skeletons = (Vec<Vec<usize>>, Vec<RingSystem>);

/// The equally senior parents that may be multiplied, each with the atoms of
/// its suffix groups, and the skeletons they were chosen from.
pub(super) fn parents(mol: &Molecule) -> Option<(Vec<(ParentSelection, AtomSet)>, Skeletons)> {
    let cyclic = get_cyclic_atoms(mol, &AtomSet::new());
    let mut groups = perception::perceive_groups(mol);
    let key = principal_groups::component_principal_key(mol, &groups, false);
    let mut excluded = AtomSet::new();
    component_group_rules::retarget_external_carbonyl_groups(
        mol,
        &mut groups,
        key.as_deref(),
        &mut excluded,
        &cyclic,
    );
    let (principal, _) =
        principal_groups::partition_principal_and_prefix_groups(&groups, key.as_deref());
    component_group_rules::exclude_nonparent_group_atoms(mol, &groups, &mut excluded, &cyclic);
    let paths = chains::find_all_carbon_paths(mol, &excluded).ok()?;
    let rings = chains::find_ring_systems(mol, &excluded).ok()?;
    let parents = parent_selection::senior_parent_candidates(
        mol, &paths, &rings, &principal, &groups, &excluded,
    )
    .into_iter()
    // P-15.3.1.1 / P-15.3.4.2 exclude acyclic hydrocarbons without suffixes.
    .filter(|p| p.is_ring || principal.iter().any(|a| p.primary_path().contains(a)))
    .map(|p| {
        let mut complete = p.atom_set();
        complete.extend(component_group_rules::principal_involved_atoms(
            mol,
            &groups,
            key.as_deref(),
            p.primary_path(),
        ));
        (p, complete)
    })
    .collect();
    Some((parents, (paths, rings)))
}

/// Two to four identical units on a linking group built from single atoms,
/// unbranched alkane chains and benzene rings, where each unit is a ring or
/// functionalized chain bearing nothing but suffix groups. Everything else
/// keeps its substitutive name.
pub(super) fn name(mol: &Molecule, component: &AtomSet) -> Option<String> {
    if mol
        .atoms
        .values()
        .any(|a| a.charge != 0 || a.has_stereo() || a.radical_electrons != 0)
        || mol.bonds.values().any(|b| b.stereo.is_some())
    {
        return None;
    }
    // A reason recorded while the unit is planned and rendered would hide
    // behind the pending multiplicative status, so look for it from the default.
    let status = mol.name_status();
    mol.set_name_status(NameStatus::Unverified(NameReason::UnauditedPath));
    let name = parts(mol, component);
    let clean = mol.name_status() == NameStatus::Unverified(NameReason::UnauditedPath);
    mol.set_name_status(status);
    name.filter(|_| clean)
}

/// A senior parent with its suffix groups and substituents, and its one bond
/// to the linking group.
struct Unit<'a> {
    parent: &'a ParentSelection,
    suffixed: &'a AtomSet,
    atoms: AtomSet,
    root: usize,
    link: usize,
}

/// One component of a linking group: a single atom with its oxo groups, an
/// unbranched alkane chain or a benzene ring.
struct Segment {
    atoms: AtomSet,
    /// Bonds to the neighbouring components and units, inner atom first.
    ports: Vec<(usize, usize)>,
}

impl Segment {
    /// The prefix for the given attachment atoms. In an arm they are the atom
    /// towards the central group, then the atom towards the unit.
    fn name(&self, mol: &Molecule, ports: &[usize], central: bool) -> Option<String> {
        let inside = |a: usize| {
            mol.get_neighbors(a)
                .iter()
                .filter(|n| self.atoms.contains(n))
                .count()
        };
        let first = *self.atoms.first()?;
        if let Some(prefix) = self.atom_prefix(mol, ports) {
            return Some(prefix.into());
        }
        if !self.atoms.iter().all(|a| mol.atoms[a].is_carbon()) || self.atoms.len() < 2 {
            return None;
        }
        let ring = mol.atoms[&first].is_aromatic;
        if ring && (self.atoms.len() != 6 || self.atoms.iter().any(|&a| inside(a) != 2)) {
            return None;
        }
        // Walk the ring from any atom, or the chain from an end.
        let start = if ring {
            first
        } else {
            *self.atoms.iter().find(|&&a| inside(a) == 1)?
        };
        let mut path = vec![start];
        while let Some(&next) = mol
            .get_neighbors(path[path.len() - 1])
            .iter()
            .find(|n| self.atoms.contains(n) && !path.contains(n))
        {
            path.push(next);
        }
        if path.len() != self.atoms.len() {
            return None;
        }
        let size = path.len();
        let mut numberings = vec![path.clone()];
        if ring {
            numberings = (0..size)
                .map(|i| path[i..].iter().chain(&path[..i]).copied().collect())
                .collect();
        }
        numberings.extend(
            numberings
                .clone()
                .into_iter()
                .map(|n: Vec<usize>| n.into_iter().rev().collect()),
        );
        let locants = numberings
            .iter()
            .filter_map(|numbering| {
                let mut locants = ports
                    .iter()
                    .map(|p| Some(numbering.iter().position(|a| a == p)? + 1))
                    .collect::<Option<Vec<_>>>()?;
                if central {
                    // P-29.3.2.2: lowest locants to the free valences.
                    locants.sort();
                    Some((locants.clone(), locants))
                } else {
                    // P-15.3.1.2.2.4: the lowest locant is nearest to the
                    // unit and is cited last, as in ethane-2,1-diyl.
                    Some((locants.iter().rev().copied().collect(), locants))
                }
            })
            .min()?
            .1;
        let cited = locants
            .iter()
            .map(usize::to_string)
            .collect::<Vec<_>>()
            .join(",");
        Some(if !ring {
            // P-29.3.2.2: propane-1,2-diyl, ethane-1,1,2-triyl.
            format!(
                "{}ane-{cited}-{}yl",
                stems::stem_for(size),
                multipliers::basic(ports.len())
            )
        } else if ports.len() == 2 {
            // P-29.6.1: the retained prefix phenylene.
            format!("{cited}-phenylene")
        } else {
            format!("benzene-{cited}-{}yl", multipliers::basic(ports.len()))
        })
    }

    /// The prefix of a single linking atom. P-15.3.1.2.1.1 (methylene, oxy,
    /// sulfanediyl, nitrilo), P-15.3.1.2.1.2 (azanediyl), P-15.3.3.2.1
    /// (methanetriyl), P-51.3.1 (phosphanetriyl), P-63.2.5 (selanediyl,
    /// tellanediyl), P-63.6 and P-35.2.3 (sulfinyl to telluronyl).
    fn atom_prefix(&self, mol: &Molecule, ports: &[usize]) -> Option<&'static str> {
        let &centre = ports.first()?;
        let atom = &mol.atoms[&centre];
        if ports.iter().any(|&p| p != centre) {
            return None;
        }
        Some(
            match (
                atom.symbol.as_str(),
                self.atoms.len() - 1,
                ports.len(),
                atom.total_h_count,
            ) {
                ("O", 0, 2, 0) => "oxy",
                ("S", 0, 2, 0) => "sulfanediyl",
                ("S", 1, 2, 0) => "sulfinyl",
                ("S", 2, 2, 0) => "sulfonyl",
                ("Se", 0, 2, 0) => "selanediyl",
                ("Se", 1, 2, 0) => "seleninyl",
                ("Se", 2, 2, 0) => "selenonyl",
                ("Te", 0, 2, 0) => "tellanediyl",
                ("Te", 1, 2, 0) => "tellurinyl",
                ("Te", 2, 2, 0) => "telluronyl",
                ("C", 0, 2, 2) => "methylene",
                ("C", 0, 3, 1) => "methanetriyl",
                ("N", 0, 2, 1) => "azanediyl",
                ("N", 0, 3, 0) => "nitrilo",
                ("P", 0, 3, 0) => "phosphanetriyl",
                _ => return None,
            },
        )
    }
}

/// The multiplicative name, or nothing where it is not built.
fn parts(mol: &Molecule, component: &AtomSet) -> Option<String> {
    let heavy = |n: &usize| mol.atoms[n].symbol != "H";
    let plain = |u: usize, v: usize| {
        mol.get_bond(u, v)
            .is_some_and(|b| b.order == 1 && !b.is_aromatic)
    };
    // A unit is a senior parent with its suffix groups and one bond towards
    // the other parents; P-51.3.2.1 multiplies all of them. Its branches
    // that hold no senior parent are its substituents (P-15.3.2.4.1).
    let (parents, _) = parents(mol)?;
    let mut units: Vec<Unit> = Vec::new();
    for (parent, suffixed) in &parents {
        let mut atoms = suffixed.clone();
        let mut links = Vec::new();
        for &u in suffixed {
            for &n in mol.get_neighbors(u) {
                if suffixed.contains(&n) || !heavy(&n) {
                    continue;
                }
                let branch = subgraph_component(mol, n, suffixed);
                if parents
                    .iter()
                    .any(|(_, other)| other != suffixed && !other.is_disjoint(&branch))
                {
                    links.push((u, n));
                } else {
                    atoms.extend(branch);
                }
            }
        }
        if let [(root, link)] = links[..]
            && units.iter().all(|u| u.suffixed != suffixed)
        {
            units.push(Unit {
                parent,
                suffixed,
                atoms,
                root,
                link,
            });
        }
    }
    let unit_atoms: AtomSet = units.iter().flat_map(|u| &u.atoms).copied().collect();
    if !(2..=4).contains(&units.len())
        || unit_atoms.len() != units.iter().map(|u| u.atoms.len()).sum::<usize>()
    {
        return None;
    }
    let linker: AtomSet = component
        .iter()
        .filter(|a| heavy(a) && !unit_atoms.contains(a))
        .copied()
        .collect();
    let cyclic = get_cyclic_atoms(mol, &AtomSet::new());
    let chain = |a: usize| {
        mol.atoms[&a].is_carbon()
            && !cyclic.contains(&a)
            && mol
                .get_neighbors(a)
                .iter()
                .all(|&n| !heavy(&n) || plain(a, n))
    };
    let oxo = |a: usize, b: usize| {
        mol.atoms[&a].symbol == "O"
            && mol.degree(a) == 1
            && mol.bond_order(a, b) == 2
            && matches!(mol.atoms[&b].symbol.as_str(), "S" | "Se" | "Te")
    };
    // Chains and aromatic rings are one component each; every other atom
    // stands alone, with its oxo groups.
    let joined: EdgeSet = mol
        .bonds
        .values()
        .filter(|b| linker.contains(&b.u) && linker.contains(&b.v))
        .filter(|b| {
            chain(b.u) && chain(b.v)
                || b.is_aromatic && mol.atoms[&b.u].is_carbon() && mol.atoms[&b.v].is_carbon()
                || oxo(b.u, b.v)
                || oxo(b.v, b.u)
        })
        .map(|b| normalize_edge(b.u, b.v))
        .collect();
    let mut owner = BTreeMap::new();
    let mut segments = Vec::new();
    for atoms in connected_components(&linker, &joined) {
        let ports: Vec<_> = atoms
            .iter()
            .flat_map(|&u| {
                mol.get_neighbors(u)
                    .iter()
                    .filter(|n| !atoms.contains(n) && heavy(n))
                    .map(move |&n| (u, n))
            })
            .collect();
        // P-51.3.1 (1): single bonds throughout. Rings other than benzene
        // are not named as linking groups, nor are heteroatom chains, which
        // have names of their own such as triselanediyl (P-68.4.1.2).
        if ports.iter().any(|&(u, n)| {
            !plain(u, n)
                || linker.contains(&n) && !mol.atoms[&u].is_carbon() && !mol.atoms[&n].is_carbon()
        }) || atoms
            .iter()
            .any(|a| cyclic.contains(a) && !mol.atoms[a].is_aromatic)
        {
            return None;
        }
        owner.extend(atoms.iter().map(|&a| (a, segments.len())));
        segments.push(Segment { atoms, ports });
    }
    // P-51.4.1.1: four heterounits call for skeletal replacement names.
    if segments
        .iter()
        .filter(|s| s.atoms.iter().any(|a| !mol.atoms[a].is_carbon()))
        .count()
        >= 4
    {
        return None;
    }
    // P-15.3.1.2.2.1: the central group is the component whose arms towards
    // the units are all alike; P-15.3.4.1.1 rules out any other linking group.
    let linker = segments.iter().find_map(|centre| {
        let mut arms = Vec::new();
        let mut used = 1;
        for &(mut from, mut at) in &centre.ports {
            let mut arm = Vec::new();
            while !units.iter().any(|u| (u.root, u.link) == (at, from)) {
                let segment = &segments[*owner.get(&at)?];
                let [a, b] = segment.ports[..] else {
                    return None;
                };
                let exit = if a == (at, from) { b } else { a };
                arm.push((
                    segment.name(mol, &[at, exit.0], false)?,
                    segment.atoms.len() > 1
                        && segment.atoms.iter().all(|a| mol.atoms[a].is_carbon()),
                ));
                used += 1;
                if used > segments.len() {
                    return None;
                }
                (from, at) = exit;
            }
            arms.push(arm);
        }
        if used != segments.len() || arms.iter().any(|arm| arm != &arms[0]) {
            return None;
        }
        let ports: Vec<_> = centre.ports.iter().map(|p| p.0).collect();
        let central = centre.name(mol, &ports, true)?;
        Some(match &arms[0][..] {
            [] if centre.atom_prefix(mol, &ports).is_some() => central,
            [] => format!("({central})"),
            // P-16.3.4 (a), (c): di(4,1-phenylene), di(tetradecane-14,1-diyl).
            [(arm, true)] => format!("({central}{}({arm}))", multipliers::basic(arms.len())),
            // P-16.3.5 (a), P-16.3.6 (a), (b): bis(oxy), bis(oxymethylene).
            arm => format!(
                "({central}{}({}))",
                multipliers::complex_(arms.len()),
                arm.iter().map(|a| a.0.as_str()).collect::<String>()
            ),
        })
    })?;
    let Unit {
        parent,
        suffixed,
        atoms,
        root,
        link,
    } = &units[0];
    // P-51.3.1 (3): every unit matches the named one from its attachment atom.
    if !units[1..].iter().all(|other| {
        super::identical_multiplicative_branches(
            mol,
            atoms,
            &other.atoms,
            suffixed,
            other.suffixed,
            *root,
            other.root,
        )
    }) {
        return None;
    }
    // P-15.3.2.4.1: substituted units are cited with their prefixes. Only
    // rings are built, and not where every position is substituted: those
    // names omit locants, as in bis(pentafluorodisilane).
    let substituted = atoms.len() > suffixed.len();
    let path = parent.primary_path();
    if substituted && (!parent.is_ring || path.iter().all(|a| mol.atoms[a].total_h_count == 0)) {
        return None;
    }
    // The unit keeps the locants it needs beside the linking group
    // (P-15.3.2.1), as in di(cyclohexane-1-carboxylic acid).
    let named = component_namer::name_component_with_atom_locants(
        mol,
        component,
        None,
        Some((parent, substituted.then_some((*root, *link)))),
    )
    .ok()?;
    // Without the linking group the engine may contract the unit to a
    // retained name such as anisole; the parent must stay the one it is in
    // the substitutive name, with every prefix cited before it.
    if substituted
        && (component_namer::name_component_with_atom_locants(
            mol,
            component,
            None,
            Some((parent, None)),
        )
        .ok()?
        .rendered_parent
            != named.rendered_parent
            || !named.text.ends_with(&named.rendered_parent)
            || !named
                .text
                .starts_with(|c: char| c.is_ascii_digit() || c == 'N'))
    {
        return None;
    }
    let locant = named.locants.get(root)?;
    let functionalized = suffixed.len() > path.len();
    // The locant 1 of a mononuclear parent is omitted (P-15.3.1.3 (a)), and
    // those names are not built.
    if path.len() < 2 || !locant.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    let text = named.rendered_parent;
    let stem = text.strip_suffix(" acid").unwrap_or(&text);
    let numerical = [
        "di", "tri", "tetr", "pent", "hex", "hept", "oct", "non", "dec", "undec", "do", "hen",
        "icos",
    ]
    .iter()
    .any(|term| stem.starts_with(term));
    // P-16.3.6 (c): Hantzsch-Widman and replacement names take 'bis' and are
    // not built; these retained names take 'di', as in 3,3'-[...]difuran
    // (P-51.3.2.1).
    if path.iter().any(|a| !mol.atoms[a].is_carbon())
        && ![
            "furan",
            "morpholine",
            "piperazine",
            "piperidine",
            "pyrazine",
            "pyridazine",
            "pyridine",
            "pyrimidine",
            "pyrrolidine",
            "selenophene",
            "tellurophene",
            "thiophene",
        ]
        .contains(&text.as_str())
    {
        return None;
    }
    // P-16.3.3 (d): diphenol, diacetic acid. P-16.3.4 (c), (d), (e):
    // functionalized parent hydrides that begin with a numerical term or cite
    // locants are enclosed, as in di(decanoic acid) and
    // di(cyclohexane-1-carboxylic acid). Other spellings are left alone.
    let unit = if stem.bytes().all(|b| b.is_ascii_lowercase()) && !numerical {
        text
    } else if functionalized
        && stem.starts_with(|c: char| c.is_ascii_lowercase())
        && stem
            .bytes()
            .all(|b| b.is_ascii_lowercase() || b.is_ascii_digit() || b == b',' || b == b'-')
    {
        format!("({text})")
    } else {
        return None;
    };
    let count = units.len();
    let locants: Vec<_> = (0..count)
        .map(|i| format!("{locant}{}", "'".repeat(i)))
        .collect();
    Some(if substituted {
        // P-16.3.5 (e): bis(2-bromobenzoic acid).
        format!(
            "{}-{linker}{}({})",
            locants.join(","),
            multipliers::complex_(count),
            named.text
        )
    } else {
        format!(
            "{}-{linker}{}{unit}",
            locants.join(","),
            multipliers::basic(count)
        )
    })
}

#[cfg(test)]
mod tests {
    use crate::test_support::{parse_strict, verify_identity};
    use nazuchem_ir::{NameReason, NameStatus};

    fn name(smiles: &str) -> crate::Name {
        crate::name_smiles(smiles, &Default::default()).unwrap()
    }

    /// The multiplicative name is preferred, the substitutive name is its
    /// alternative, and both parse back to the structure.
    fn assert_pin(smiles: &str, expected: &str, substitutive: &str) -> Result<(), String> {
        let name = name(smiles);
        assert_eq!(name.text, expected, "{smiles}");
        assert_eq!(name.status, NameStatus::Preferred, "{smiles}");
        let [alternative] = &name.alternatives[..] else {
            panic!("{smiles}: {:?}", name.alternatives);
        };
        assert_eq!(alternative.name, substitutive, "{smiles}");
        assert_eq!(
            alternative.status,
            NameStatus::NotPreferred(NameReason::Multiplicative),
            "{smiles}"
        );
        verify_identity(smiles, &parse_strict(expected)?)?;
        verify_identity(smiles, &parse_strict(substitutive)?)
    }

    #[test]
    fn identical_units_on_one_linking_atom_take_the_multiplicative_pin() -> Result<(), String> {
        for (smiles, expected, substitutive) in [
            // P-63.2.4.2, P-15.3.1.3: locants are always cited.
            ("c1ccccc1Oc1ccccc1", "1,1'-oxydibenzene", "phenoxybenzene"),
            // P-63.2.5.
            (
                "c1ccccc1Sc1ccccc1",
                "1,1'-sulfanediyldibenzene",
                "(phenylsulfanyl)benzene",
            ),
            // P-63.6, with the Se and Te prefixes of P-35.2.3.
            (
                "O=S(c1ccccc1)c1ccccc1",
                "1,1'-sulfinyldibenzene",
                "(benzenesulfinyl)benzene",
            ),
            (
                "O=[Se](c1ccccc1)c1ccccc1",
                "1,1'-seleninyldibenzene",
                "[oxo(phenyl)-lambda^4-selanyl]benzene",
            ),
            (
                "O=[Se](=O)(c1ccccc1)c1ccccc1",
                "1,1'-selenonyldibenzene",
                "[dioxo(phenyl)-lambda^6-selanyl]benzene",
            ),
            (
                "O=[Te](=O)(c1ccccc1)c1ccccc1",
                "1,1'-telluronyldibenzene",
                "[dioxo(phenyl)-lambda^6-tellanyl]benzene",
            ),
            (
                "O=S(=O)(c1ccccc1)c1ccccc1",
                "1,1'-sulfonyldibenzene",
                "(benzenesulfonyl)benzene",
            ),
            (
                "O=[Te](c1ccccc1)c1ccccc1",
                "1,1'-tellurinyldibenzene",
                "[oxo(phenyl)-lambda^4-tellanyl]benzene",
            ),
            // P-63.2.5, P-35.2.2.
            (
                "c1ccccc1[Se]c1ccccc1",
                "1,1'-selanediyldibenzene",
                "(phenylselanyl)benzene",
            ),
            (
                "c1ccccc1[Te]c1ccccc1",
                "1,1'-tellanediyldibenzene",
                "(phenyltellanyl)benzene",
            ),
            // P-16.3.4 (e): several suffix groups, and a ketone.
            (
                "O=C(O)c1cc(Oc2cc(C(=O)O)cc(C(=O)O)c2)cc(C(=O)O)c1",
                "5,5'-oxydi(benzene-1,3-dicarboxylic acid)",
                "5-(3,5-dicarboxyphenoxy)benzene-1,3-dicarboxylic acid",
            ),
            (
                "CC(=O)COCC(C)=O",
                "1,1'-oxydi(propan-2-one)",
                "1-(2-oxopropoxy)propan-2-one",
            ),
            // P-44.1.2.2: the ring is senior to the chain.
            (
                "c1ccccc1Cc1ccccc1",
                "1,1'-methylenedibenzene",
                "benzylbenzene",
            ),
            // P-51.3.1.
            (
                "S(CC(=O)O)CC(=O)O",
                "2,2'-sulfanediyldiacetic acid",
                "[(carboxymethyl)sulfanyl]acetic acid",
            ),
            (
                "P(CC(=O)O)(CC(=O)O)CC(=O)O",
                "2,2',2''-phosphanetriyltriacetic acid",
                "[bis(carboxymethyl)phosphanyl]acetic acid",
            ),
            (
                "OC(=O)C1CCC(CC1)OC1CCC(CC1)C(=O)O",
                "4,4'-oxydi(cyclohexane-1-carboxylic acid)",
                "4-[(4-carboxycyclohexyl)oxy]cyclohexane-1-carboxylic acid",
            ),
            // P-15.3.2.1.
            (
                "Oc1ccc(Oc2ccc(O)cc2)cc1",
                "4,4'-oxydiphenol",
                "4-(4-hydroxyphenoxy)phenol",
            ),
            (
                "Nc1ccc(Oc2ccc(N)cc2)cc1",
                "4,4'-oxydianiline",
                "4-(4-aminophenoxy)aniline",
            ),
            (
                "OCCSCCO",
                "2,2'-sulfanediyldi(ethan-1-ol)",
                "2-[(2-hydroxyethyl)sulfanyl]ethan-1-ol",
            ),
            (
                "NCCOCCN",
                "2,2'-oxydi(ethan-1-amine)",
                "2-(2-aminoethoxy)ethan-1-amine",
            ),
            (
                "OCCN(CCO)CCO",
                "2,2',2''-nitrilotri(ethan-1-ol)",
                "2-[bis(2-hydroxyethyl)amino]ethan-1-ol",
            ),
            (
                "OC(=O)COCC(=O)O",
                "2,2'-oxydiacetic acid",
                "(carboxymethoxy)acetic acid",
            ),
            (
                "OC(=O)CCOCCC(=O)O",
                "3,3'-oxydipropanoic acid",
                "3-(2-carboxyethoxy)propanoic acid",
            ),
            (
                "OS(=O)(=O)c1ccc(Oc2ccc(S(=O)(=O)O)cc2)cc1",
                "4,4'-oxydi(benzene-1-sulfonic acid)",
                "4-(4-sulfophenoxy)benzene-1-sulfonic acid",
            ),
            // P-16.3.4 (c), (d).
            (
                "OC(=O)CCCCCCCCCOCCCCCCCCCC(=O)O",
                "10,10'-oxydi(decanoic acid)",
                "10-[(9-carboxynonyl)oxy]decanoic acid",
            ),
            (
                "O=CCCCCOCCCCC=O",
                "5,5'-oxydi(pentanal)",
                "5-[(5-oxopentyl)oxy]pentanal",
            ),
            // P-13.6.1, P-62.2.5.1, P-63.1.5.
            (
                "N#Cc1ccccc1Cc1ccccc1C#N",
                "2,2'-methylenedibenzonitrile",
                "2-[(2-cyanophenyl)methyl]benzonitrile",
            ),
            (
                "N#Cc1ccc(Nc2ccc(C#N)cc2)cc1",
                "4,4'-azanediyldibenzonitrile",
                "4-[(4-cyanophenyl)amino]benzonitrile",
            ),
            (
                "OC(=O)CN(CC(=O)O)CC(=O)O",
                "2,2',2''-nitrilotriacetic acid",
                "[bis(carboxymethyl)amino]acetic acid",
            ),
            (
                "SC1CCCC1SC1CCCC1S",
                "2,2'-sulfanediyldi(cyclopentane-1-thiol)",
                "2-[(2-sulfanylcyclopentyl)sulfanyl]cyclopentane-1-thiol",
            ),
            (
                "Oc1ccccc1Sc1ccccc1O",
                "2,2'-sulfanediyldiphenol",
                "2-[(2-hydroxyphenyl)sulfanyl]phenol",
            ),
        ] {
            assert_pin(smiles, expected, substitutive)?;
        }
        Ok(())
    }

    #[test]
    fn concatenated_and_chain_linking_groups_take_the_multiplicative_pin() -> Result<(), String> {
        for (smiles, expected, substitutive) in [
            // P-52.2.5.3 Example 2, P-29.6.1.
            (
                "c1ccccc1Oc2ccc(Oc3ccccc3)cc2",
                "1,1'-[1,4-phenylenebis(oxy)]dibenzene",
                "1,4-diphenoxybenzene",
            ),
            (
                "c1ccccc1Cc2ccc(Cc3ccccc3)cc2",
                "1,1'-[1,4-phenylenebis(methylene)]dibenzene",
                "1,4-dibenzylbenzene",
            ),
            (
                "c1ccccc1Oc2cc(Oc3ccccc3)cc(Oc4ccccc4)c2",
                "1,1',1''-[benzene-1,3,5-triyltris(oxy)]tribenzene",
                "1,3,5-triphenoxybenzene",
            ),
            // P-15.3.1.2.2.1: benzene-1,2,4-triyltris(oxy).
            (
                "c1ccccc1Oc2ccc(Oc3ccccc3)c(Oc3ccccc3)c2",
                "1,1',1''-[benzene-1,2,4-triyltris(oxy)]tribenzene",
                "1,2,4-triphenoxybenzene",
            ),
            // P-51.3.1.
            (
                "O=C(O)c1ccc(OCCOc2ccc(C(=O)O)cc2)cc1",
                "4,4'-[ethane-1,2-diylbis(oxy)]dibenzoic acid",
                "4-[2-(4-carboxyphenoxy)ethoxy]benzoic acid",
            ),
            (
                "OC(=O)COCCOCC(=O)O",
                "2,2'-[ethane-1,2-diylbis(oxy)]diacetic acid",
                "[2-(carboxymethoxy)ethoxy]acetic acid",
            ),
            // P-15.3.1.2.2.4: the arm is numbered from the unit.
            (
                "OC(=O)COCCOCCOCC(=O)O",
                "2,2'-[oxybis(ethane-2,1-diyloxy)]diacetic acid",
                "{2-[2-(carboxymethoxy)ethoxy]ethoxy}acetic acid",
            ),
            (
                "OCCOCCOCCOCCO",
                "2,2'-[oxybis(ethane-2,1-diyloxy)]di(ethan-1-ol)",
                "2-{2-[2-(2-hydroxyethoxy)ethoxy]ethoxy}ethan-1-ol",
            ),
            // P-15.3.2.1, P-16.3.4 (c): 'di' before one chain or ring arm.
            (
                "Oc1ccc(cc1)CCCCCCCCCCCCCCOCCCCCCCCCCCCCCc1ccc(O)cc1",
                "4,4'-[oxydi(tetradecane-14,1-diyl)]diphenol",
                "4-(14-{[14-(4-hydroxyphenyl)tetradecyl]oxy}tetradecyl)phenol",
            ),
            // P-15.3.1.2.2.1, P-16.3.6 (a), (b): 'bis' before single atoms.
            (
                "c1ccccc1COCc1ccccc1",
                "1,1'-[oxybis(methylene)]dibenzene",
                "[(benzyloxy)methyl]benzene",
            ),
            (
                "c1ccccc1SCSc1ccccc1",
                "1,1'-[methylenebis(sulfanediyl)]dibenzene",
                "{[(phenylsulfanyl)methyl]sulfanyl}benzene",
            ),
            (
                "c1ccccc1S(=O)(=O)CS(=O)(=O)c1ccccc1",
                "1,1'-[methylenebis(sulfonyl)]dibenzene",
                "[(benzenesulfonyl)methanesulfonyl]benzene",
            ),
            // P-15.3.3.1, P-29.3.2.2: unsymmetrical central chains.
            (
                "CC(c1ccccc1)Cc1ccccc1",
                "1,1'-(propane-1,2-diyl)dibenzene",
                "(1-phenylpropan-2-yl)benzene",
            ),
            (
                "CC(C)(c1ccc(O)cc1)c1ccc(O)cc1",
                "4,4'-(propane-2,2-diyl)diphenol",
                "4-[2-(4-hydroxyphenyl)propan-2-yl]phenol",
            ),
            (
                "c1ccccc1CCc1ccccc1",
                "1,1'-(ethane-1,2-diyl)dibenzene",
                "(2-phenylethyl)benzene",
            ),
            // P-51.3.2.1: the most units are multiplied.
            (
                "OC(=O)c1ccc(cc1)CC(c1ccc(cc1)C(=O)O)c1ccc(cc1)C(=O)O",
                "4,4',4''-(ethane-1,1,2-triyl)tribenzoic acid",
                "4-[1,2-bis(4-carboxyphenyl)ethyl]benzoic acid",
            ),
            // P-15.3.3.2.1: methanetriyl.
            (
                "c1ccccc1C(c1ccccc1)c1ccccc1",
                "1,1',1''-methanetriyltribenzene",
                "(diphenylmethyl)benzene",
            ),
        ] {
            assert_pin(smiles, expected, substitutive)?;
        }
        Ok(())
    }

    #[test]
    fn heterocyclic_and_substituted_units_take_the_multiplicative_pin() -> Result<(), String> {
        for (smiles, expected, substitutive) in [
            // P-15.3.2.4.2, P-51.3.2.1: oxydipyridine, difuran.
            (
                "c1ccncc1Oc1cccnc1",
                "3,3'-oxydipyridine",
                "3-[(pyridin-3-yl)oxy]pyridine",
            ),
            (
                "c1ccsc1Sc1cccs1",
                "2,2'-sulfanediyldithiophene",
                "2-[(thiophen-2-yl)sulfanyl]thiophene",
            ),
            // P-15.3.2.4.1, P-16.3.5 (e).
            (
                "Brc1ccc(Oc2ccc(Br)cc2)cc1",
                "1,1'-oxybis(4-bromobenzene)",
                "1-bromo-4-(4-bromophenoxy)benzene",
            ),
            (
                "Brc1ccc(Sc2ccc(Br)cc2)cc1",
                "1,1'-sulfanediylbis(4-bromobenzene)",
                "1-bromo-4-[(4-bromophenyl)sulfanyl]benzene",
            ),
            (
                "O=S(c1ccc(Br)cc1)c1ccc(Br)cc1",
                "1,1'-sulfinylbis(4-bromobenzene)",
                "1-bromo-4-(4-bromobenzene-1-sulfinyl)benzene",
            ),
            (
                "O=S(=O)(c1ccc(Br)cc1)c1ccc(Br)cc1",
                "1,1'-sulfonylbis(4-bromobenzene)",
                "1-bromo-4-(4-bromobenzene-1-sulfonyl)benzene",
            ),
            (
                "OC(=O)c1ccc(Oc2ccc(C(=O)O)c(Br)c2)cc1Br",
                "4,4'-oxybis(2-bromobenzoic acid)",
                "2-bromo-4-(3-bromo-4-carboxyphenoxy)benzoic acid",
            ),
            (
                "CC(Cl)c1cncc(Oc2cncc(C(C)Cl)c2)c1",
                "3,3'-oxybis[5-(1-chloroethyl)pyridine]",
                "3-(1-chloroethyl)-5-{[5-(1-chloroethyl)pyridin-3-yl]oxy}pyridine",
            ),
            (
                "CN(C)C(=O)c1ccc(Oc2ccc(C(=O)N(C)C)cc2)cc1",
                "4,4'-oxybis(N,N-dimethylbenzamide)",
                "4-[4-(dimethylcarbamoyl)phenoxy]-N,N-dimethylbenzamide",
            ),
            (
                "Oc1ccc(cc1Cl)C(C)(C)c1ccc(O)c(Cl)c1",
                "4,4'-(propane-2,2-diyl)bis(2-chlorophenol)",
                "2-chloro-4-[2-(3-chloro-4-hydroxyphenyl)propan-2-yl]phenol",
            ),
            (
                "Brc1ccc(OCCOc2ccc(Br)cc2)cc1",
                "1,1'-[ethane-1,2-diylbis(oxy)]bis(4-bromobenzene)",
                "1-bromo-4-[2-(4-bromophenoxy)ethoxy]benzene",
            ),
        ] {
            assert_pin(smiles, expected, substitutive)?;
        }
        Ok(())
    }

    #[test]
    fn other_multiplicative_pins_keep_the_substitutive_name_and_its_status() {
        for (smiles, expected) in [
            // P-15.3.1.3: the linking group takes locant 1 here, against the
            // lowest locant set of the substitutive name; this unit numbering
            // is not built.
            (
                "ClCc1ccc(Cc2ccc(CCl)c(Br)c2)cc1Br",
                "2-bromo-4-{[3-bromo-4-(chloromethyl)phenyl]methyl}-1-(chloromethyl)benzene",
            ),
            (
                "Brc1ccc(Oc2ccc(Br)c(Br)c2)cc1Br",
                "1,2-dibromo-4-(3,4-dibromophenoxy)benzene",
            ),
            // The unit alone would be contracted to anisole.
            (
                "COc1ccc(Cc2ccc(OC)cc2)cc1",
                "1-methoxy-4-[(4-methoxyphenyl)methyl]benzene",
            ),
            // P-15.3.2.4.1: fully substituted units omit locants, and
            // substituted chain units are not built.
            (
                "Fc1c(F)c(F)c(F)c(F)c1Oc1c(F)c(F)c(F)c(F)c1F",
                "pentafluoro(2,3,4,5,6-pentafluorophenoxy)benzene",
            ),
            (
                "CN(C)C(=O)CCSCCC(=O)N(C)C",
                "3-{[3-(dimethylamino)-3-oxopropyl]sulfanyl}-N,N-dimethylpropanamide",
            ),
            // P-15.3.2.1: ring assembly, branched chain and other ring
            // linking groups are not built.
            (
                "O=C(O)CCCCCCCCCOc1ccc(-c2ccc(OCCCCCCCCCC(=O)O)cc2)cc1",
                "10-({4'-[(9-carboxynonyl)oxy][1,1'-biphenyl]-4-yl}oxy)decanoic acid",
            ),
            (
                "Oc1ccc(OCC(C)COc2ccc(O)cc2)cc1",
                "4-[3-(4-hydroxyphenoxy)-2-methylpropoxy]phenol",
            ),
            (
                "c1ccccc1C(c1ccccc1)(c1ccccc1)c1ccccc1",
                "(triphenylmethyl)benzene",
            ),
            // P-68.4.1.2: the PIN is 4,4'-triselanediyldibenzoic acid.
            (
                "O=C(O)c1ccc([Se][Se][Se]c2ccc(C(=O)O)cc2)cc1",
                "4-({[(4-carboxyphenyl)selanyl]selanyl}selanyl)benzoic acid",
            ),
            // P-16.3.6 (c): Hantzsch-Widman units take 'bis'.
            ("C1CCOC1OC1CCCO1", "2-[(oxolan-2-yl)oxy]oxolane"),
            // P-45.6.2: configuration needs primed stereodescriptors.
            (
                "C[C@H](O)CSC[C@H](C)O",
                "(2S)-1-{[(2S)-2-hydroxypropyl]sulfanyl}propan-2-ol",
            ),
            // A unit name outside the audited spellings.
            ("O=CCOCC=O", "2-(2-oxoethoxy)ethanal"),
            // P-16.3.4 (f), P-16.3.6 (d): names with brackets.
            (
                "C1CC2CC1CC2OC1CC2CCC1C2",
                "2-[(bicyclo[2.2.1]heptan-2-yl)oxy]bicyclo[2.2.1]heptane",
            ),
        ] {
            let name = name(smiles);
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(
                name.status,
                NameStatus::NotPreferred(NameReason::Multiplicative),
                "{smiles}"
            );
            assert!(name.alternatives.is_empty(), "{smiles}");
        }
    }

    #[test]
    fn senior_parents_and_unlike_units_are_not_multiplied() {
        for (smiles, expected) in [
            // P-62.2.2.1: not azanediyldibenzene.
            ("c1ccc(Nc2ccccc2)cc1", "N-phenylaniline"),
            // P-51.3.2.3: the senior class is the parent.
            ("O=C(c1ccccc1)c1ccccc1", "diphenylmethanone"),
            ("c1ccc(P(c2ccccc2)c2ccccc2)cc1", "triphenylphosphane"),
            // P-15.3.4.2, P-63.6: alkanes are not multiplied.
            ("CCOCC", "ethoxyethane"),
            ("CCS(=O)(=O)CC", "(ethanesulfonyl)ethane"),
            // P-51.3.3: unlike attachment locants.
            (
                "OC(=O)C1CCCCC1CC1CCC(CC1)C(=O)O",
                "2-[(4-carboxycyclohexyl)methyl]cyclohexane-1-carboxylic acid",
            ),
            ("OC(=O)CCCC(=O)O", "pentanedioic acid"),
        ] {
            let name = name(smiles);
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, NameStatus::Preferred, "{smiles}");
        }
    }
}
