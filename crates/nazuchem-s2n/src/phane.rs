// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: ReShiki | src/naming/openclatura/phane.rs | 5f559fd482e237307a589dfcb7ea333cee0a6e13
// Copyright (c) 2026 Ameyanagi and ReShiki contributors

//! Cyclic and linear phanes with disjoint monocyclic amplificants
//! (PhI-2/PhI-3, P-26).
use super::{
    graph_kernel::{Edge, adjacency_from_edges, connected_components, normalize_edge},
    molecule::{AtomSet, Molecule},
    rules::multipliers,
};
use std::collections::{BTreeMap, BTreeSet};

struct Amplificant {
    prefix: &'static str,
    // PhI-2.2.3: heterocycles precede benzene, with N > O > S.
    seniority: usize,
    numberings: Vec<BTreeMap<usize, usize>>,
}

/// P-52.2.5.3, P-51.3.2: a linear phane has at least four whole
/// ring systems and seven skeletal nodes along a ring-to-ring path.
pub(super) fn has_linear_skeleton(mol: &Molecule, splits: &BTreeSet<Edge>) -> bool {
    let systems = super::ring_assembly::cyclic_units(mol);
    let rings: AtomSet = systems.iter().filter_map(|s| s.first().copied()).collect();
    let mut collapsed: BTreeMap<_, _> = mol
        .atoms
        .values()
        .filter(|a| a.symbol != "H")
        .map(|a| (a.idx, a.idx))
        .collect();
    for system in &systems {
        if let Some(&root) = system.first() {
            for &atom in system {
                collapsed.insert(atom, root);
            }
        }
    }
    let nodes = collapsed.values().copied().collect();
    let edges = mol
        .bonds
        .values()
        .filter(|bond| !splits.contains(&normalize_edge(bond.u, bond.v)))
        .filter_map(|bond| {
            let (&u, &v) = (collapsed.get(&bond.u)?, collapsed.get(&bond.v)?);
            (u != v).then_some(normalize_edge(u, v))
        })
        .collect();
    let adjacency = adjacency_from_edges(&nodes, &edges);
    // Collapsing every cyclic component leaves a forest: avoiding the
    // preceding node visits exactly the simple paths from each terminal ring.
    rings.iter().any(|&start| {
        let mut pending = vec![(start, start, 1, 1)];
        while let Some((previous, current, length, ring_count)) = pending.pop() {
            if rings.contains(&current) && ring_count >= 4 && length >= 7 {
                return true;
            }
            for &next in &adjacency[&current] {
                if next != previous {
                    pending.push((
                        current,
                        next,
                        length + 1,
                        ring_count + usize::from(rings.contains(&next)),
                    ));
                }
            }
        }
        false
    })
}

/// The phane name and whether it is a linear phane. Collapsing the
/// amplificants leaves one cycle or one unbranched chain; every atom is
/// accounted for.
pub(super) fn name(mol: &Molecule) -> Option<(String, bool)> {
    if mol.bonds.len() <= mol.atoms.len()
        || mol.atoms.values().any(|a| {
            a.charge != 0 || a.isotope.is_some() || a.has_stereo() || a.radical_electrons != 0
        })
        || mol.bonds.values().any(|b| b.stereo.is_some())
    {
        return None;
    }
    let aromatic: AtomSet = mol
        .atoms
        .values()
        .filter(|a| a.is_aromatic)
        .map(|a| a.idx)
        .collect();
    let aromatic_edges = mol
        .bonds
        .values()
        .filter(|b| b.is_aromatic)
        .map(|b| normalize_edge(b.u, b.v))
        .collect();
    let mut ring_atoms = aromatic.clone();
    let mut ring_edges: BTreeSet<_> = aromatic_edges;
    // An isolated saturated six-carbon ring is a cyclohexana amplificant.
    for system in super::ring_assembly::cyclic_units(mol) {
        if system.len() == 6
            && system.iter().all(|i| {
                let a = &mol.atoms[i];
                !a.is_aromatic
                    && a.is_carbon()
                    && mol.degree(*i) + a.total_h_count == 4
                    && mol
                        .get_neighbors(*i)
                        .iter()
                        .filter(|n| system.contains(n))
                        .count()
                        == 2
            })
        {
            for &atom in &system {
                ring_atoms.insert(atom);
                for &neighbor in mol.get_neighbors(atom) {
                    if system.contains(&neighbor) {
                        ring_edges.insert(normalize_edge(atom, neighbor));
                    }
                }
            }
        }
    }
    let adjacency = adjacency_from_edges(&ring_atoms, &ring_edges);
    let mut amplificants = BTreeMap::new();
    let mut collapsed: BTreeMap<_, _> = mol.atoms.keys().map(|&i| (i, i)).collect();
    for atoms in connected_components(&ring_atoms, &ring_edges) {
        if !matches!(atoms.len(), 5 | 6) || atoms.iter().any(|i| adjacency[i].len() != 2) {
            return None;
        }
        let start = *atoms.first()?;
        let mut ring = vec![start];
        let mut previous = start;
        let mut current = *adjacency[&start].first()?;
        while current != start {
            ring.push(current);
            let next = *adjacency[&current].iter().find(|&&n| n != previous)?;
            previous = current;
            current = next;
        }
        let hetero: Vec<_> = ring
            .iter()
            .filter(|&&i| !mol.atoms[&i].is_carbon())
            .collect();
        let (prefix, seniority) = match (ring.len(), hetero.as_slice()) {
            (6, []) if !aromatic.contains(&start) => ("cyclohexana", 4),
            (6, []) => ("benzena", 3),
            (6, [i]) if mol.atoms[i].symbol == "N" => ("pyridina", 0),
            (5, [i]) if mol.atoms[i].symbol == "O" => ("furana", 1),
            (5, [i]) if mol.atoms[i].symbol == "S" => ("thiophena", 2),
            _ => return None,
        };
        if ring.iter().any(|i| {
            let a = &mol.atoms[i];
            if a.is_carbon() {
                mol.degree(*i) + a.total_h_count != if a.is_aromatic { 3 } else { 4 }
            } else {
                mol.degree(*i) != 2 || a.total_h_count != 0
            }
        }) {
            return None;
        }
        let mut numberings = vec![];
        for offset in 0..ring.len() {
            if !hetero.is_empty() && ring[offset] != *hetero[0] {
                continue;
            }
            for reverse in [false, true] {
                numberings.push(
                    (0..ring.len())
                        .map(|i| {
                            (
                                ring[(offset + if reverse { ring.len() - i } else { i })
                                    % ring.len()],
                                i + 1,
                            )
                        })
                        .collect(),
                );
            }
        }
        for i in atoms {
            collapsed.insert(i, start);
        }
        amplificants.insert(
            start,
            Amplificant {
                prefix,
                seniority,
                numberings,
            },
        );
    }
    if amplificants.len() < 2 {
        return None;
    }
    // Number the parent hydride first, then the substituent (P-26.4.3).
    let mut branches = mol.atoms.values().filter_map(|atom| {
        let group = match (atom.symbol.as_str(), atom.total_h_count) {
            ("O", 1) => "ol",
            ("C", 3) => "methyl",
            ("F", 0) => "fluoro",
            ("Cl", 0) => "chloro",
            ("Br", 0) => "bromo",
            ("I", 0) => "iodo",
            _ => return None,
        };
        let neighbors = mol.get_neighbors(atom.idx);
        let [ring_atom] = neighbors else { return None };
        (aromatic.contains(ring_atom) && mol.bond_order(atom.idx, *ring_atom) == 1)
            .then_some((atom.idx, *ring_atom, group))
    });
    let mut branch = branches.next();
    if branches.next().is_some() {
        return None;
    }
    if let Some((atom, _, _)) = branch {
        collapsed.remove(&atom);
    }
    // In the presence of a senior carboxylic acid, a bridging carbonate
    // can contribute an O-C-O path and an oxo prefix (P-41, P-65.6.3.2.3).
    let mut carbonyl = None;
    for atom in mol
        .atoms
        .values()
        .filter(|a| a.is_carbon() && !a.is_aromatic && mol.degree(a.idx) == 3)
    {
        let neighbors = mol.get_neighbors(atom.idx);
        let Some(&oxygen) = neighbors.iter().find(|&&n| {
            mol.atoms[&n].symbol == "O" && mol.degree(n) == 1 && mol.bond_order(atom.idx, n) == 2
        }) else {
            continue;
        };
        if atom.total_h_count != 0
            || mol.atoms[&oxygen].total_h_count != 0
            || neighbors
                .iter()
                .any(|&n| n != oxygen && mol.bond_order(atom.idx, n) != 1)
        {
            return None;
        }
        let others: Vec<_> = neighbors.iter().copied().filter(|&n| n != oxygen).collect();
        if let Some(&ring) = others.iter().find(|n| aromatic.contains(n)) {
            let hydroxyl = *others.iter().find(|&&n| n != ring)?;
            if mol.atoms[&hydroxyl].symbol != "O"
                || mol.atoms[&hydroxyl].total_h_count != 1
                || mol.degree(hydroxyl) != 1
                || branch.is_some()
            {
                return None;
            }
            branch = Some((atom.idx, ring, "carboxylic acid"));
            for i in [atom.idx, hydroxyl, oxygen] {
                collapsed.remove(&i);
            }
        } else if others
            .iter()
            .all(|n| mol.atoms[n].symbol == "O" && mol.degree(*n) == 2)
        {
            if carbonyl.replace(atom.idx).is_some() {
                return None;
            }
            collapsed.remove(&oxygen);
        }
    }
    if carbonyl.is_some() && !matches!(branch, Some((_, _, "carboxylic acid"))) {
        return None;
    }
    let nodes: AtomSet = collapsed.values().copied().collect();
    let mut edges = BTreeSet::new();
    // For each skeleton edge, retain the real atom at either end for attachments.
    let mut attachments = BTreeMap::new();
    for b in mol.bonds.values() {
        let (Some(&u), Some(&v)) = (collapsed.get(&b.u), collapsed.get(&b.v)) else {
            continue;
        };
        if u == v {
            if !ring_edges.contains(&normalize_edge(b.u, b.v)) {
                return None;
            }
            continue;
        }
        if b.order != 1 || b.is_aromatic || !edges.insert(normalize_edge(u, v)) {
            return None;
        }
        attachments.insert((u, v), b.u);
        attachments.insert((v, u), b.v);
    }
    let adjacency = adjacency_from_edges(&nodes, &edges);
    let terminals: Vec<_> = nodes
        .iter()
        .copied()
        .filter(|n| adjacency[n].len() == 1)
        .collect();
    let linear = terminals.len() == 2;
    // P-52.2.5.1: a cyclophane has at least six nodes; a linear phane has at
    // least seven, with four rings of which two are terminal.
    if adjacency.values().any(|neighbors| neighbors.len() > 2)
        || connected_components(&nodes, &edges).len() != 1
        || !(if linear {
            (7..=1000).contains(&nodes.len())
                && amplificants.len() >= 4
                && terminals.iter().all(|i| amplificants.contains_key(i))
        } else {
            branch.is_none() && ring_atoms == aromatic && terminals.is_empty() && nodes.len() >= 6
        })
    {
        return None;
    }
    for &i in &nodes {
        if amplificants.contains_key(&i) {
            let ends: Vec<_> = adjacency[&i]
                .iter()
                .map(|&n| attachments[&(i, n)])
                .collect();
            // P-52.2.5.2.1: adjacent attachments make a fused ring system
            // of a cyclophane.
            if !linear && (ends[0] == ends[1] || mol.get_neighbors(ends[0]).contains(&ends[1])) {
                return None;
            }
        } else {
            let a = &mol.atoms[&i];
            if Some(i) == carbonyl {
                continue;
            }
            let hydrogens = match a.symbol.as_str() {
                "C" => 2,
                "O" | "S" => 0,
                "N" => 1,
                _ => return None,
            };
            if a.is_aromatic || a.total_h_count != hydrogens || mol.degree(i) != 2 {
                return None;
            }
        }
    }
    let order_starts: Vec<_> = nodes.iter().copied().collect();
    // Each numbering is compared numerically, never by the rendered name.
    let mut best = None;

    // A chain is numbered from either end, a cycle from any node either way.
    for &start in if linear {
        &terminals[..]
    } else {
        &order_starts[..]
    } {
        for &second in &adjacency[&start] {
            let mut order = vec![start];
            let mut previous = start;
            let mut current = second;
            while current != start {
                order.push(current);
                let Some(&next) = adjacency[&current].iter().find(|&&n| n != previous) else {
                    break;
                };
                previous = current;
                current = next;
            }
            let locants: BTreeMap<_, _> =
                order.iter().enumerate().map(|(i, &n)| (n, i + 1)).collect();
            let superatoms: Vec<_> = order
                .iter()
                .filter(|n| amplificants.contains_key(n))
                .map(|n| locants[n])
                .collect();
            let hetero = |symbol: Option<&str>| -> Vec<usize> {
                order
                    .iter()
                    .filter(|n| {
                        !amplificants.contains_key(n)
                            && mol.atoms[n].symbol != "C"
                            && symbol.is_none_or(|s| mol.atoms[n].symbol == s)
                    })
                    .map(|n| locants[n])
                    .collect()
            };
            let seniority: Vec<Vec<usize>> = (0..5)
                .map(|rank| {
                    order
                        .iter()
                        .filter(|n| amplificants.get(n).is_some_and(|a| a.seniority == rank))
                        .map(|n| locants[n])
                        .collect()
                })
                .collect();
            let mut numbered = vec![];
            let mut attachment_sets = vec![];
            let mut attachment_order = vec![];
            let mut branch_locant = None;
            let mut center = None;
            for &node in &order {
                let Some(amplificant) = amplificants.get(&node) else {
                    continue;
                };
                let mut neighbors: Vec<_> = adjacency[&node].iter().copied().collect();
                neighbors.sort_by_key(|n| locants[n]);
                let (set, pair, position) = amplificant
                    .numberings
                    .iter()
                    .map(|numbering| {
                        let pair: Vec<_> = neighbors
                            .iter()
                            .map(|&n| numbering[&attachments[&(node, n)]])
                            .collect();
                        let mut set = pair.clone();
                        set.sort();
                        (
                            set,
                            pair,
                            branch.and_then(|(_, i, _)| numbering.get(&i).copied()),
                        )
                    })
                    .min()?;
                if let Some(position) = position {
                    branch_locant = Some((locants[&node], position));
                }
                // Only the middle node can remain unpaired when a linear
                // skeleton has two otherwise equivalent numberings.
                if linear && locants[&node] * 2 == nodes.len() + 1 && pair.len() == 2 {
                    let reversible = amplificant.numberings.iter().any(|numbering| {
                        numbering[&attachments[&(node, neighbors[0])]] == pair[1]
                            && numbering[&attachments[&(node, neighbors[1])]] == pair[0]
                    });
                    if !reversible {
                        center = Some((amplificant.seniority, set.clone(), pair.clone()));
                    }
                }
                attachment_sets.push((amplificant.seniority, locants[&node], set));
                attachment_order.push(pair.clone());
                numbered.push((amplificant.prefix, locants[&node], pair));
            }
            // PhI-3.2.2 compares attachment sets by amplificant seniority,
            // then by superatom locant within each cited amplification prefix.
            // P-26.4.2.3 precedes comparison of all attachment locants
            // only when a single unsymmetrical amplificant remains.
            let single_asymmetric = center.and_then(|(rank, set, pair)| {
                (attachment_sets
                    .iter()
                    .filter(|(r, _, s)| *r == rank && *s == set)
                    .count()
                    == 1)
                    .then_some(pair)
            });
            attachment_sets.sort_by_key(|&(rank, locant, _)| (rank, locant));
            let attachment_sets: Vec<_> =
                attachment_sets.into_iter().map(|(_, _, set)| set).collect();
            // PhI-4.1/4.5: heteroatoms are placed on the numbering of the
            // carbocyclic phane and decide only the choices left by PhI-3.
            let key = (
                superatoms,
                seniority,
                attachment_sets,
                // P-26.4.2.4 completes the hydrocarbon numbering before
                // P-26.5.1 places the heteroatoms on it.
                single_asymmetric,
                attachment_order,
                hetero(None),
                hetero(Some("O")),
                hetero(Some("S")),
                hetero(Some("N")),
                branch_locant,
                carbonyl.map(|c| locants[&c]),
            );
            if best.as_ref().is_none_or(|(old, _)| key < *old) {
                let mut parts = vec![];
                for (symbol, prefix) in [("O", "oxa"), ("S", "thia"), ("N", "aza")] {
                    let positions = hetero(Some(symbol));
                    if !positions.is_empty() {
                        let multiplier = if positions.len() == 1 {
                            String::new()
                        } else {
                            multipliers::basic(positions.len())
                        };
                        parts.push(format!(
                            "{}-{multiplier}{prefix}",
                            positions
                                .iter()
                                .map(usize::to_string)
                                .collect::<Vec<_>>()
                                .join(",")
                        ));
                    }
                }
                for rank in 0..5 {
                    let Some(amplificant) = amplificants.values().find(|a| a.seniority == rank)
                    else {
                        continue;
                    };
                    let members: Vec<_> = numbered
                        .iter()
                        .filter(|(prefix, _, _)| *prefix == amplificant.prefix)
                        .collect();
                    let mut groups: Vec<(&Vec<usize>, Vec<usize>)> = vec![];
                    for (_, locant, pair) in &members {
                        let locant = *locant;
                        if let Some((_, positions)) = groups.iter_mut().find(|(p, _)| *p == pair) {
                            positions.push(locant);
                        } else {
                            groups.push((pair, vec![locant]));
                        }
                    }
                    let descriptor = groups
                        .iter()
                        .map(|(pair, positions)| {
                            let list = |locants: &[usize]| {
                                locants
                                    .iter()
                                    .map(usize::to_string)
                                    .collect::<Vec<_>>()
                                    .join(",")
                            };
                            format!("{}({})", list(positions), list(pair))
                        })
                        .collect::<Vec<_>>()
                        .join(",");
                    let multiplier = if members.len() == 1 {
                        String::new()
                    } else {
                        multipliers::basic(members.len())
                    };
                    parts.push(format!("{descriptor}-{multiplier}{}", amplificant.prefix));
                }
                let mut name = format!(
                    "{}{}{}phane",
                    parts.join("-"),
                    if linear { "" } else { "cyclo" },
                    multipliers::numerical_term(nodes.len())
                );
                if let (Some((node, ring)), Some((_, _, group))) = (branch_locant, branch) {
                    let superscript = ['⁰', '¹', '²', '³', '⁴', '⁵', '⁶'];
                    let locant = format!("{node}{}", superscript[ring]);
                    if group == "ol" {
                        name.pop();
                        name.push_str(&format!("-{locant}-ol"));
                    } else if group == "carboxylic acid" {
                        name.push_str(&format!("-{locant}-carboxylic acid"));
                        if let Some(carbon) = carbonyl {
                            name = format!("{}-oxo-{name}", locants[&carbon]);
                        }
                    } else {
                        name = format!("{locant}-{group}-{name}");
                    }
                }
                // The parser also accepts flattened composite locants. If
                // no explicit superscript distinguishes the notation, a
                // skeletal replacement locant can collide with a ring atom.
                let ambiguous = linear
                    && branch.is_none()
                    && hetero(None).iter().any(|&i| {
                        let locant = i.to_string();
                        numbered.iter().any(|(_, n, _)| {
                            locant
                                .strip_prefix(&n.to_string())
                                .filter(|s| !s.starts_with('0'))
                                .and_then(|s| s.parse::<usize>().ok())
                                .is_some_and(|ring| {
                                    (1..=amplificants[&order[*n - 1]].numberings[0].len())
                                        .contains(&ring)
                                })
                        })
                    });
                best = Some((key, (!ambiguous).then_some(name)));
            }
        }
    }
    best.and_then(|(_, name)| name.map(|name| (name, linear)))
}

#[cfg(test)]
mod tests;
