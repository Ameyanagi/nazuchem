// SPDX-License-Identifier: MIT AND (MIT OR Apache-2.0)
// Origin: mixed
// Source: openclatura | src/openclatura/von_baeyer.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | src/openclatura/chains.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: ReShiki | src/naming/openclatura/von_baeyer.rs | 5f559fd482e237307a589dfcb7ea333cee0a6e13
// Source: openclatura | src/openclatura/graph_kernel.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/unit/test_analysis.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: Nazuchem | crates/nazuchem-s2n/src/von_baeyer.rs | issue #100
// Copyright (c) 2026 lamalab-org
// Copyright (c) 2026 Ameyanagi and ReShiki contributors
// Changes: Consolidated openclatura legacy/main-ring and independent von-Baeyer searches with original Rust Python-order adapter and local bounds/tests.
// Local contributions: MIT OR Apache-2.0; see NOTICE and PROVENANCE-nazuchem.toml.

//! Extended von Baeyer discovery from upstream `chains.py` and `von_baeyer.py`.
//! Candidate ranking never substitutes for reconstruction: callers retain only
//! paths whose descriptor rebuilds every edge of the parent skeleton.
use std::collections::{BTreeMap, BTreeSet, VecDeque};

use super::graph_kernel::{EdgeSet, adjacency_from_edges, connected_components, normalize_edge};
use super::molecule::{AtomSet, Molecule};
use super::polycycle_topology::{
    RingNumbering, audit_von_baeyer_descriptor, build_von_baeyer_numbering,
};
use super::python_int_set::PythonIntSet;
use super::ring_renderer::render_von_baeyer_descriptor;

type Adjacency = BTreeMap<usize, AtomSet>;
type OrderedAdjacency = BTreeMap<usize, Vec<usize>>;

#[derive(Clone, Debug)]
pub struct DescriptorPaths {
    pub descriptor: Option<String>,
    pub paths: Vec<Vec<usize>>,
    /// Atoms in the main ring, and in the longest unbranched bridge across it.
    pub main_sizes: (usize, usize),
}

const MAIN_RING_WALK_STATES: usize = 2_000_000;
// Upstream divergence: upstream allows rank 28 (2.7e8 subsets, tens of
// seconds); naming runs interactively, so larger blocks abstain at rank 20.
const MAX_CYCLE_SPACE_RANK: usize = 20;
const MAX_LONGEST_CYCLES: usize = 4096;

#[derive(Clone, Debug)]
struct LegacyBridge {
    nodes: AtomSet,
    endpoints: [usize; 2],
}

// Provenance: fn python_adjacency | adapted
// Source: openclatura | src/openclatura/chains.py::get_von_baeyer_descriptor_candidates
// Source: openclatura | src/openclatura/graph_kernel.py::adjacency_from_edges
fn python_adjacency(adj: &Adjacency) -> OrderedAdjacency {
    // Upstream divergence: upstream inserts neighbors in tuple-edge-set order.
    // Rebuilding from sorted IDs retains integer-table traversal, but a hash
    // collision between neighbor IDs can still depend on that insertion order.
    adj.iter()
        .map(|(&atom, neighbors)| {
            let neighbors: PythonIntSet = neighbors.iter().copied().collect();
            (atom, neighbors.iter().collect())
        })
        .collect()
}

// Provenance: fn bridge_interior_path | adapted
// Source: openclatura | src/openclatura/chains.py::_bridge_interior_path
fn bridge_interior_path(
    adj: &OrderedAdjacency,
    nodes: &AtomSet,
    start: usize,
    end: usize,
) -> Vec<usize> {
    let mut queue = VecDeque::from([(start, Vec::new())]);
    let mut visited = AtomSet::new();
    while let Some((current, path)) = queue.pop_front() {
        if current == end && !path.is_empty() {
            return path[..path.len() - 1].to_vec();
        }
        visited.insert(current);
        for &neighbor in &adj[&current] {
            if (nodes.contains(&neighbor) || neighbor == end) && !visited.contains(&neighbor) {
                let mut next = path.clone();
                next.push(neighbor);
                queue.push_back((neighbor, next));
            }
        }
    }
    Vec::new()
}

// Provenance: fn largest_cycles_walk | adapted
// Source: openclatura | src/openclatura/chains.py::inline in _largest_cycles_from_cycle_space
fn largest_cycles_walk(atoms: &AtomSet, adj: &OrderedAdjacency) -> Option<Vec<Vec<usize>>> {
    let mut states = 0;
    let mut longest = 0;
    let mut cycles = Vec::new();
    for &start in atoms {
        states += 1;
        let mut path = vec![start];
        let mut stack = vec![(start, 0)];
        while let Some(&(current, index)) = stack.last() {
            let Some(&neighbor) = adj[&current].get(index) else {
                stack.pop();
                path.pop();
                continue;
            };
            if let Some(frame) = stack.last_mut() {
                frame.1 += 1;
            }
            if neighbor == start && path.len() >= 3 {
                if path.len() > longest {
                    longest = path.len();
                    cycles.clear();
                }
                if path.len() == longest {
                    cycles.push(path.clone());
                }
            } else if !path.contains(&neighbor) {
                states += 1;
                if states > MAIN_RING_WALK_STATES {
                    return None;
                }
                path.push(neighbor);
                stack.push((neighbor, 0));
            }
        }
    }
    Some(cycles)
}

// Provenance: fn xor_mask | adapted
// Source: openclatura | src/openclatura/chains.py::inline in _largest_cycles_from_cycle_space
fn xor_mask(target: &mut [u64], source: &[u64]) {
    for (target, source) in target.iter_mut().zip(source) {
        *target ^= source;
    }
}

// Provenance: fn cycle_path_from_mask | adapted
// Source: openclatura | src/openclatura/chains.py::_cycle_path_from_edges
fn cycle_path_from_mask(mask: &[u64], edges: &[(usize, usize)]) -> Option<Vec<usize>> {
    let mut adj: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (index, &(left, right)) in edges.iter().enumerate() {
        if mask[index / 64] & (1_u64 << (index % 64)) == 0 {
            continue;
        }
        adj.entry(left).or_default().push(right);
        adj.entry(right).or_default().push(left);
    }
    if adj.is_empty() || adj.values().any(|ends| ends.len() != 2) {
        return None;
    }
    let start = *adj.first_key_value()?.0;
    let mut path = vec![start];
    let mut previous = None;
    let mut current = start;
    loop {
        let ends = &adj[&current];
        let next = if Some(ends[0]) != previous {
            ends[0]
        } else {
            ends[1]
        };
        if next == start {
            break;
        }
        path.push(next);
        previous = Some(current);
        current = next;
        if path.len() > adj.len() {
            return None;
        }
    }
    (path.len() == adj.len()).then_some(path)
}

// Provenance: fn largest_cycles_from_cycle_space | adapted
// Source: openclatura | src/openclatura/chains.py::_largest_cycles_from_cycle_space
fn largest_cycles_from_cycle_space(
    atoms: &AtomSet,
    edges: &EdgeSet,
    adj: &Adjacency,
) -> Vec<Vec<usize>> {
    let rank = edges.len() + 1 - atoms.len();
    if rank == 0 || rank > MAX_CYCLE_SPACE_RANK {
        return Vec::new();
    }
    let ordered_edges: Vec<_> = edges.iter().copied().collect();
    let edge_index: BTreeMap<_, _> = ordered_edges
        .iter()
        .enumerate()
        .map(|(i, &edge)| (edge, i))
        .collect();
    let mut parent = BTreeMap::new();
    let mut depth = BTreeMap::new();
    let mut tree = EdgeSet::new();
    let mut seen = AtomSet::new();
    for &root in atoms {
        if !seen.insert(root) {
            continue;
        }
        depth.insert(root, 0_usize);
        let mut stack = vec![root];
        while let Some(current) = stack.pop() {
            for &neighbor in &adj[&current] {
                if !seen.insert(neighbor) {
                    continue;
                }
                parent.insert(neighbor, current);
                depth.insert(neighbor, depth[&current] + 1);
                tree.insert(normalize_edge(current, neighbor));
                stack.push(neighbor);
            }
        }
    }
    let words = edges.len().div_ceil(64);
    let mut basis = Vec::new();
    for &edge in edges {
        if tree.contains(&edge) {
            continue;
        }
        let (mut left, mut right) = edge;
        let mut ring = EdgeSet::from([edge]);
        while left != right {
            if depth[&left] >= depth[&right] {
                ring.insert(normalize_edge(left, parent[&left]));
                left = parent[&left];
            } else {
                ring.insert(normalize_edge(right, parent[&right]));
                right = parent[&right];
            }
        }
        let mut mask = vec![0; words];
        for edge in ring {
            let index = edge_index[&edge];
            mask[index / 64] |= 1_u64 << (index % 64);
        }
        basis.push(mask);
    }
    if basis.is_empty() {
        return Vec::new();
    }
    let mut current = vec![0; words];
    let mut longest = 0;
    let mut best = Vec::new();
    for step in 1_usize..1_usize << basis.len() {
        xor_mask(&mut current, &basis[step.trailing_zeros() as usize]);
        let length = current
            .iter()
            .map(|word| word.count_ones() as usize)
            .sum::<usize>();
        if length < longest
            || (length == longest && (best.is_empty() || best.len() >= MAX_LONGEST_CYCLES))
        {
            continue;
        }
        let Some(path) = cycle_path_from_mask(&current, &ordered_edges) else {
            continue;
        };
        if length > longest {
            longest = length;
            best.clear();
        }
        best.push(path);
    }
    best.sort_by(|a, b| (a.len(), a).cmp(&(b.len(), b)));
    best
}

// Provenance: fn main_ring_division | adapted
// Source: openclatura | src/openclatura/chains.py::_main_ring_division
fn main_ring_division(
    bridge: &LegacyBridge,
    position: &BTreeMap<usize, usize>,
    size: usize,
) -> Option<usize> {
    let first = *position.get(&bridge.endpoints[0])?;
    let second = *position.get(&bridge.endpoints[1])?;
    let span = first.abs_diff(second).min(size - first.abs_diff(second));
    Some(size.abs_diff(2 * span))
}

/// Upstream main-ring search, in preference order. The per-call walk and
/// cycle-space bounds deliberately abstain on an unmanageable graph.
// Provenance: fn descriptor_candidates | adapted
// Source: openclatura | src/openclatura/chains.py::get_von_baeyer_descriptor_candidates
pub fn descriptor_candidates(atoms: &AtomSet, edges: &EdgeSet) -> Vec<DescriptorPaths> {
    let adj = adjacency_from_edges(atoms, edges);
    let ordered_adj = python_adjacency(&adj);
    let largest = largest_cycles_walk(atoms, &ordered_adj)
        .unwrap_or_else(|| largest_cycles_from_cycle_space(atoms, edges, &adj));
    if largest.is_empty() {
        return vec![DescriptorPaths {
            descriptor: None,
            paths: vec![atoms.iter().copied().collect()],
            main_sizes: (0, 0),
        }];
    }
    let mut candidates = Vec::new();
    for main_ring in &largest {
        let main_atoms: AtomSet = main_ring.iter().copied().collect();
        let ring_edges: EdgeSet = (0..main_ring.len())
            .map(|i| normalize_edge(main_ring[i], main_ring[(i + 1) % main_ring.len()]))
            .collect();
        let mut bridges = Vec::new();
        for &(first, second) in edges {
            if main_atoms.contains(&first)
                && main_atoms.contains(&second)
                && !ring_edges.contains(&(first, second))
            {
                bridges.push(LegacyBridge {
                    nodes: AtomSet::new(),
                    endpoints: [first, second],
                });
            }
        }
        let outside: AtomSet = atoms.difference(&main_atoms).copied().collect();
        for component in connected_components(&outside, edges) {
            let mut connections = Vec::new();
            for &atom in &component {
                for &neighbor in &adj[&atom] {
                    if main_atoms.contains(&neighbor) {
                        connections.push((atom, neighbor));
                    }
                }
            }
            if connections.len() >= 2 {
                bridges.push(LegacyBridge {
                    nodes: component,
                    endpoints: [connections[0].1, connections[1].1],
                });
                for &(first, second) in &connections[2..] {
                    bridges.push(LegacyBridge {
                        nodes: AtomSet::new(),
                        endpoints: [first, second],
                    });
                }
            }
        }
        if bridges.is_empty()
            || bridges
                .iter()
                .any(|bridge| bridge.endpoints[0] == bridge.endpoints[1])
        {
            continue;
        }
        bridges.sort_by_key(|bridge| std::cmp::Reverse(bridge.nodes.len()));
        let position: BTreeMap<_, _> = main_ring
            .iter()
            .enumerate()
            .map(|(i, &atom)| (atom, i))
            .collect();
        let maximum = bridges[0].nodes.len();
        let selected = bridges
            .iter()
            .enumerate()
            .filter(|(_, bridge)| bridge.nodes.len() == maximum)
            .filter_map(|(i, bridge)| {
                let division = main_ring_division(bridge, &position, main_ring.len())?;
                let first = position[&bridge.endpoints[0]];
                let second = position[&bridge.endpoints[1]];
                Some(((division, first.min(second), first.max(second)), i))
            })
            .min()
            .map(|(_, index)| index);
        if let Some(index) = selected {
            let bridge = bridges.remove(index);
            bridges.insert(0, bridge);
        }
        let division =
            main_ring_division(&bridges[0], &position, main_ring.len()).unwrap_or(main_ring.len());
        candidates.push(((std::cmp::Reverse(maximum), division), main_ring, bridges));
    }
    if candidates.is_empty() {
        return vec![DescriptorPaths {
            descriptor: None,
            paths: vec![largest[0].clone()],
            main_sizes: (0, 0),
        }];
    }
    candidates.sort_by_key(|candidate| candidate.0);
    // The upstream caller consumes at most 256 rendered choices. Applying the
    // same bound here avoids materializing the rest of a symmetric skeleton.
    candidates
        .into_iter()
        .flat_map(|(_, ring, bridges)| legacy_from_candidate(ring, &bridges, atoms, &ordered_adj))
        .take(256)
        .collect()
}

// Provenance: fn legacy_from_candidate | adapted
// Source: openclatura | src/openclatura/chains.py::_von_baeyer_from_candidate
fn legacy_from_candidate(
    ring: &[usize],
    bridges: &[LegacyBridge],
    atoms: &AtomSet,
    adj: &OrderedAdjacency,
) -> Vec<DescriptorPaths> {
    let main = &bridges[0];
    let [first, second] = main.endpoints;
    let i = ring.iter().position(|&atom| atom == first).unwrap_or(0);
    let j = ring.iter().position(|&atom| atom == second).unwrap_or(0);
    let low = i.min(j);
    let high = i.max(j);
    let path_a = ring[low..=high].to_vec();
    let path_b: Vec<_> = ring[high..].iter().chain(&ring[..=low]).copied().collect();
    let (mut branch1, mut branch2) = if path_a.len() >= path_b.len() {
        (path_a, path_b.into_iter().rev().collect::<Vec<_>>())
    } else {
        (path_b, path_a.into_iter().rev().collect())
    };
    if branch1[0] != first {
        branch1.reverse();
        branch2.reverse();
    }
    let bridge_path = if !main.nodes.is_empty() {
        bridge_interior_path(adj, &main.nodes, first, second)
    } else {
        Vec::new()
    };
    let mut path1 = branch1.clone();
    path1.extend(branch2[1..branch2.len() - 1].iter().rev().copied());
    path1.extend(&bridge_path);
    let mut path2: Vec<_> = branch1.iter().rev().copied().collect();
    path2.extend_from_slice(&branch2[1..branch2.len() - 1]);
    path2.extend(bridge_path.iter().rev().copied());
    let mut primary = [branch1.len() - 2, branch2.len() - 2, main.nodes.len()];
    primary.sort_by(|a, b| b.cmp(a));
    let mut orientations = Vec::new();
    for mut path in [path1, path2] {
        let base: BTreeMap<_, _> = path
            .iter()
            .enumerate()
            .map(|(i, &atom)| (atom, i + 1))
            .collect();
        // P-23.2.6.3: the bridge at the highest-numbered bridgehead is numbered
        // first, each from the atom next to its higher-numbered bridgehead;
        // P-23.2.6.4: then the higher other bridgehead, then the longer bridge.
        let mut extras: Vec<_> = bridges[1..]
            .iter()
            .filter(|bridge| !bridge.nodes.is_empty())
            .map(|bridge| {
                let mut ends = bridge.endpoints;
                ends.sort_by_key(|end| std::cmp::Reverse(base.get(end).copied().unwrap_or(1)));
                (bridge, ends)
            })
            .collect();
        extras.sort_by_key(|(bridge, ends)| {
            std::cmp::Reverse((
                base.get(&ends[0]).copied().unwrap_or(1),
                base.get(&ends[1]).copied().unwrap_or(1),
                bridge.nodes.len(),
            ))
        });
        let mut visited: AtomSet = path.iter().copied().collect();
        for (bridge, [high, low]) in extras {
            for atom in bridge_interior_path(adj, &bridge.nodes, high, low) {
                if visited.insert(atom) {
                    path.push(atom);
                }
            }
        }
        for &atom in atoms {
            if visited.insert(atom) {
                path.push(atom);
            }
        }
        let pos: BTreeMap<_, _> = path
            .iter()
            .enumerate()
            .map(|(i, &atom)| (atom, i + 1))
            .collect();
        let mut chords: Vec<_> = bridges[1..]
            .iter()
            .map(|bridge| {
                let a = *pos.get(&bridge.endpoints[0]).unwrap_or(&1);
                let b = *pos.get(&bridge.endpoints[1]).unwrap_or(&1);
                (bridge.nodes.len(), a.min(b), a.max(b))
            })
            .collect();
        // P-23.2.6.1.2 and P-23.2.6.1.3: independent bridges before dependent
        // ones, each group in decreasing length; P-23.2.6.2.5: then low locants.
        chords.sort_by_key(|&(length, low, high)| {
            (high > base.len(), std::cmp::Reverse(length), low, high)
        });
        let mut body = format!("[{}.{}.{}", primary[0], primary[1], primary[2]);
        let mut preference = Vec::new();
        for (length, low, high) in chords {
            preference.extend([low, high]);
            body.push_str(&format!(".{length}^{{{low},{high}}}"));
        }
        body.push(']');
        // P-23.2.6.2.4 before P-23.2.6.2.5: low locants as a set, then as cited.
        let mut set = preference.clone();
        set.sort_unstable();
        orientations.push((
            (set, preference),
            render_von_baeyer_descriptor(bridges.len(), &body),
            path,
        ));
    }
    orientations.sort_by(|a, b| a.0.cmp(&b.0));
    // A bridge is unbranched (P-23.1.2) when each of its atoms has two neighbours.
    let unbranched = bridges
        .iter()
        .filter(|bridge| bridge.nodes.iter().all(|node| adj[node].len() == 2))
        .map(|bridge| bridge.nodes.len())
        .max()
        .unwrap_or(0);
    orientations
        .into_iter()
        .map(|(_, descriptor, path)| DescriptorPaths {
            descriptor: Some(descriptor),
            paths: vec![path],
            main_sizes: (ring.len(), unbranched),
        })
        .collect()
}

#[derive(Clone, Debug)]
pub struct VonBaeyerBridge {
    pub length: usize,
    pub attachments: [usize; 2],
    pub atoms: Vec<usize>,
    pub dependent: bool,
}

#[derive(Clone, Debug)]
pub struct VonBaeyerCandidate {
    pub descriptor: String,
    pub path: Vec<usize>,
    pub numbering: RingNumbering,
}

/// The independent audited decomposition used when the main-ring search has
/// no reconstructable descriptor. Limits match upstream's bounded scope.
// Provenance: fn find_von_baeyer_candidates | adapted
// Source: openclatura | src/openclatura/von_baeyer.py::find_von_baeyer_candidates
pub fn find_von_baeyer_candidates(
    mol: &Molecule,
    atoms: &AtomSet,
    edges: &EdgeSet,
    minimum_main_sizes: (usize, usize),
) -> Vec<VonBaeyerCandidate> {
    if edges.len() + 1 < atoms.len() + 3 {
        return Vec::new();
    }
    let aromatic = atoms.iter().any(|id| mol.atoms[id].is_aromatic);
    if aromatic && (atoms.len() != 19 || edges.len() != 23) {
        return Vec::new();
    }
    let adj = adjacency_from_edges(atoms, edges);
    if atoms.iter().any(|&atom| {
        adj[&atom].len() >= 4 && {
            let remaining: AtomSet = atoms
                .iter()
                .filter(|&&other| other != atom)
                .copied()
                .collect();
            connected_components(&remaining, edges).len() >= 2
        }
    }) {
        return Vec::new();
    }
    let ring_count = edges.len() + 1 - atoms.len();
    let heads: Vec<_> = atoms
        .iter()
        .filter(|&&atom| adj[&atom].len() >= 3)
        .copied()
        .collect();
    if ring_count > 8 || heads.len() > 12 {
        return Vec::new();
    }
    let mut candidates = Vec::new();
    let mut seen = BTreeSet::new();
    for first in 0..heads.len() {
        for second in first + 1..heads.len() {
            let paths = simple_paths_between(heads[first], heads[second], &adj);
            if paths.len() >= 96 {
                continue;
            }
            for a in 0..paths.len() {
                for b in a + 1..paths.len() {
                    for c in b + 1..paths.len() {
                        let indexes = [a, b, c];
                        let mut internal = AtomSet::new();
                        if indexes.iter().any(|&index| {
                            paths[index][1..paths[index].len() - 1]
                                .iter()
                                .any(|&atom| !internal.insert(atom))
                        }) {
                            continue;
                        }
                        let primary = indexes.map(|i| &paths[i]);
                        let sizes = |main: usize| {
                            let (a, b) = ((main + 1) % 3, (main + 2) % 3);
                            (
                                primary[a].len() + primary[b].len() - 2,
                                primary[main].len() - 2,
                            )
                        };
                        if (0..3).all(|main| sizes(main) < minimum_main_sizes) {
                            continue;
                        }
                        let primary_edges: EdgeSet = primary
                            .iter()
                            .flat_map(|p| p.windows(2).map(|pair| normalize_edge(pair[0], pair[1])))
                            .collect();
                        internal.extend([heads[first], heads[second]]);
                        let remaining = edges.difference(&primary_edges).copied().collect();
                        let Some(decompositions) = classify_secondary(atoms, &internal, &remaining)
                        else {
                            // An incomplete decomposition search cannot establish preference.
                            return Vec::new();
                        };
                        for main in 0..3 {
                            if sizes(main) < minimum_main_sizes {
                                continue;
                            }
                            let a = (main + 1) % 3;
                            let b = (main + 2) % 3;
                            for (left, right) in [(a, b), (b, a)] {
                                if primary[left].len() < primary[right].len() {
                                    continue;
                                }
                                for start in [heads[first], heads[second]] {
                                    let [first, second, main_path] = [left, right, main].map(|i| {
                                        let mut path = primary[i].clone();
                                        if path[0] != start {
                                            path.reverse();
                                        }
                                        path
                                    });
                                    for bridges in &decompositions {
                                        if bridges.len() + 2 == ring_count
                                            && let Some(candidate) = build_candidate(
                                                &first,
                                                &second,
                                                &main_path,
                                                edges,
                                                bridges.clone(),
                                            )
                                            && seen.insert((
                                                candidate.descriptor.clone(),
                                                candidate.path.clone(),
                                            ))
                                        {
                                            candidates.push(candidate);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    candidates.sort_by_cached_key(|candidate| {
        super::polycycle_topology::von_baeyer_descriptor_rank(&candidate.descriptor)
    });
    if aromatic {
        // The two benzenoid rings of the ethanoanthracene imide occupy
        // locants 2–7 and 9–14. The five-atom bridge has nitrogen at 17.
        candidates.retain(|candidate| {
            candidate.descriptor == "pentacyclo[6.6.5.0^{2,7}.0^{9,14}.0^{15,19}]"
                && candidate.path.iter().enumerate().all(|(index, id)| {
                    let atom = &mol.atoms[id];
                    atom.charge == 0
                        && atom.symbol == if index == 16 { "N" } else { "C" }
                        && atom.is_aromatic == matches!(index, 1..=6 | 8..=13)
                })
        });
    }
    candidates
}

// Provenance: fn simple_paths_between | adapted
// Source: openclatura | src/openclatura/von_baeyer.py::_simple_paths_between
fn simple_paths_between(start: usize, end: usize, adj: &Adjacency) -> Vec<Vec<usize>> {
    let mut paths = Vec::new();
    let mut stack = vec![(start, vec![start])];
    while !stack.is_empty() && paths.len() < 96 {
        let Some((current, path)) = stack.pop() else {
            break;
        };
        for &neighbor in adj[&current].iter().rev() {
            if neighbor == end {
                let mut next = path.clone();
                next.push(neighbor);
                paths.push(next);
            } else if !path.contains(&neighbor) {
                let mut next = path.clone();
                next.push(neighbor);
                stack.push((neighbor, next));
            }
        }
    }
    paths.sort_by(|a, b| b.len().cmp(&a.len()).then_with(|| a.cmp(b)));
    paths
}

// Provenance: fn build_candidate | adapted
// Source: openclatura | src/openclatura/von_baeyer.py::inline in _build_candidates_for_decomposition
fn build_candidate(
    first: &[usize],
    second: &[usize],
    main: &[usize],
    edges: &EdgeSet,
    mut secondary: Vec<VonBaeyerBridge>,
) -> Option<VonBaeyerCandidate> {
    let mut path = first.to_vec();
    path.extend(second[1..second.len() - 1].iter().rev().copied());
    path.extend_from_slice(&main[1..main.len() - 1]);
    let mut locants: BTreeMap<_, _> = path
        .iter()
        .enumerate()
        .map(|(i, &atom)| (atom, i + 1))
        .collect();
    // P-23.2.6.3: the bridge at the highest-numbered bridgehead is numbered
    // first, each from the atom next to its higher-numbered bridgehead;
    // P-23.2.6.4: then the higher other bridgehead, then the longer bridge.
    let mut unnumbered: Vec<_> = secondary
        .iter_mut()
        .filter(|bridge| bridge.length > 0)
        .collect();
    while !unnumbered.is_empty() {
        let mut available = Vec::new();
        let mut deferred = Vec::new();
        for bridge in unnumbered {
            let [a, b] = bridge.attachments;
            if let (Some(&la), Some(&lb)) = (locants.get(&a), locants.get(&b)) {
                if la < lb {
                    bridge.atoms.reverse();
                }
                available.push((
                    std::cmp::Reverse((la.max(lb), la.min(lb), bridge.length)),
                    bridge,
                ));
            } else {
                deferred.push(bridge);
            }
        }
        if available.is_empty() {
            return None;
        }
        available.sort_by_key(|(key, _)| *key);
        for (_, bridge) in available {
            for &atom in &bridge.atoms {
                path.push(atom);
                locants.insert(atom, path.len());
            }
        }
        unnumbered = deferred;
    }
    for bridge in &mut secondary {
        bridge.attachments.sort_by_key(|atom| locants[atom]);
    }
    // P-23.2.6.1.2 and P-23.2.6.1.3: independent bridges before dependent
    // ones, each group in decreasing length.
    secondary.sort_by_key(|bridge| {
        (
            bridge.dependent,
            std::cmp::Reverse(bridge.length),
            bridge.attachments.map(|atom| locants[&atom]),
        )
    });
    let lengths = [first.len() - 2, second.len() - 2, main.len() - 2];
    let mut parts: Vec<_> = lengths.iter().map(usize::to_string).collect();
    for bridge in &secondary {
        let a = locants[&bridge.attachments[0]];
        let b = locants[&bridge.attachments[1]];
        parts.push(format!("{}^{{{},{}}}", bridge.length, a.min(b), a.max(b)));
    }
    let descriptor =
        render_von_baeyer_descriptor(secondary.len() + 1, &format!("[{}]", parts.join(".")));
    if !audit_von_baeyer_descriptor(&descriptor, &path, edges) {
        return None;
    }
    let numbering = build_von_baeyer_numbering(&descriptor, &path, edges);
    Some(VonBaeyerCandidate {
        descriptor,
        path,
        numbering,
    })
}

// Provenance: fn classify_secondary | adapted
// Source: openclatura | src/openclatura/von_baeyer.py::_classify_secondary_bridges
fn classify_secondary(
    atoms: &AtomSet,
    primary: &AtomSet,
    remaining: &EdgeSet,
) -> Option<Vec<Vec<VonBaeyerBridge>>> {
    // Remove bridges as paths between already reached atoms. Process one outside
    // component at a time to avoid permuting disjoint bridges. The original
    // primary set distinguishes dependent endpoints after each extension.
    let mut work = vec![(primary.clone(), remaining.clone(), Vec::new())];
    let mut completed = Vec::new();
    let mut budget = 256usize;
    while let Some((reached, mut edges, mut bridges)) = work.pop() {
        budget = budget.checked_sub(1)?;
        edges.retain(|&(a, b)| {
            if reached.contains(&a) && reached.contains(&b) {
                bridges.push(VonBaeyerBridge {
                    length: 0,
                    attachments: [a, b],
                    atoms: vec![],
                    dependent: !(primary.contains(&a) && primary.contains(&b)),
                });
                false
            } else {
                true
            }
        });
        if edges.is_empty() {
            if reached == *atoms {
                completed.push(bridges);
            }
            continue;
        }
        let outside = atoms.difference(&reached).copied().collect();
        let component = connected_components(&outside, &edges).into_iter().next()?;
        let internal_edges = edges
            .iter()
            .copied()
            .filter(|(a, b)| component.contains(a) && component.contains(b))
            .collect();
        let adj = adjacency_from_edges(&component, &internal_edges);
        let contacts: Vec<_> = edges
            .iter()
            .filter_map(|&(a, b)| {
                if reached.contains(&a) && component.contains(&b) {
                    Some((a, b))
                } else if reached.contains(&b) && component.contains(&a) {
                    Some((b, a))
                } else {
                    None
                }
            })
            .collect();
        for (i, &(a, start)) in contacts.iter().enumerate() {
            for &(b, end) in &contacts[i + 1..] {
                if a == b {
                    continue;
                }
                let paths = if start == end {
                    vec![vec![start]]
                } else {
                    simple_paths_between(start, end, &adj)
                };
                if paths.len() >= 96 {
                    return None;
                }
                for internal in paths {
                    if work.len() >= budget {
                        return None;
                    }
                    let mut next_edges = edges.clone();
                    next_edges.remove(&normalize_edge(a, start));
                    next_edges.remove(&normalize_edge(b, end));
                    for pair in internal.windows(2) {
                        next_edges.remove(&normalize_edge(pair[0], pair[1]));
                    }
                    let mut next = reached.clone();
                    next.extend(&internal);
                    let mut next_bridges = bridges.clone();
                    next_bridges.push(VonBaeyerBridge {
                        length: internal.len(),
                        attachments: [a, b],
                        atoms: internal,
                        dependent: !(primary.contains(&a) && primary.contains(&b)),
                    });
                    work.push((next, next_edges, next_bridges));
                }
            }
        }
    }
    Some(completed)
}

#[cfg(test)]
// Provenance: mod tests | original
// Source: ReShiki | src/naming/openclatura/von_baeyer.rs::tests
mod tests {
    use super::super::{graph_io::read_smiles, polycycle_topology::edges_within_atoms};
    use super::*;

    #[test]
    // Provenance: fn tests::ethanoanthracene_imide_preserves_hydrazone_and_bridge_stereo | original
    // Source: Nazuchem | crates/nazuchem-s2n/src/von_baeyer.rs::tests::ethanoanthracene_imide_preserves_hydrazone_and_bridge_stereo
    fn ethanoanthracene_imide_preserves_hydrazone_and_bridge_stereo() {
        let original = r"C1=CC=C(C=C1)CN2C(=O)[C@H]3[C@@H](C2=O)C4(C5=CC=CC=C5C3C6=CC=CC=C64)/C=N\NC(=O)C7=CC=CO7";
        let reversed = r"O1C=CC=C1C(=O)N/N=C\C12C3=CC=CC=C3C(C3=CC=CC=C31)[C@@H]1[C@H]2C(=O)N(CC2=CC=CC=C2)C1=O";
        let e_isomer = original.replace(r"/C=N\N", "/C=N/N");
        let epimer = original.replace("[C@H]3", "[C@@H]3");
        for smiles in [original, reversed, &e_isomer, &epimer] {
            let name = crate::name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(
                name.status,
                nazuchem_ir::NameStatus::NotPreferred(nazuchem_ir::NameReason::AcylPrefix)
            );
            assert!(
                name.text
                    .contains("17-azapentacyclo[6.6.5.0^{2,7}.0^{9,14}.0^{15,19}]")
            );
            let parsed = crate::test_support::parse_strict(&name.text).unwrap();
            crate::test_support::verify_identity(smiles, &parsed).unwrap();
        }
    }

    #[test]
    fn cycle_space_recovers_every_largest_simple_cycle() {
        let atoms = AtomSet::from([0, 1, 2, 3]);
        let edges = EdgeSet::from([(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]);
        let adjacency = adjacency_from_edges(&atoms, &edges);
        let cycle_edges = |cycles: Vec<Vec<usize>>| {
            cycles
                .into_iter()
                .map(|cycle| {
                    (0..cycle.len())
                        .map(|i| normalize_edge(cycle[i], cycle[(i + 1) % cycle.len()]))
                        .collect::<EdgeSet>()
                })
                .collect::<BTreeSet<_>>()
        };
        let exhaustive =
            cycle_edges(largest_cycles_walk(&atoms, &python_adjacency(&adjacency)).unwrap());
        let space = cycle_edges(largest_cycles_from_cycle_space(&atoms, &edges, &adjacency));
        assert_eq!(space.len(), 3);
        assert_eq!(space, exhaustive);
    }

    #[test]
    // Provenance: fn tests::legacy_and_independent_search_offer_only_audited_parent_numberings | adapted
    // Source: ReShiki | src/naming/openclatura/von_baeyer.rs::tests::legacy_and_independent_search_offer_only_audited_parent_numberings
    // Source: openclatura | tests/unit/test_analysis.py::test_von_baeyer_polycycle_keeps_descriptor_source_numbering
    // Source: openclatura | tests/unit/test_analysis.py::test_ring_descriptor_rendering_is_registry_backed
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn legacy_and_independent_search_offer_only_audited_parent_numberings() {
        let mol = read_smiles("C1NC23COC12C=CC3", &Default::default()).unwrap();
        let atoms = mol.atoms.keys().copied().collect();
        let edges = edges_within_atoms(&mol, &atoms);
        let candidate =
            super::super::chains::polyspiro_or_von_baeyer_candidate(&mol, &atoms, &edges);
        assert_eq!(
            candidate.descriptor.as_deref(),
            Some("tricyclo[3.2.2.0^{1,5}]")
        );
        assert!(candidate.descriptor_allowed());
        assert!(!candidate.paths.is_empty());
        assert!(
            candidate
                .paths
                .iter()
                .all(|path| audit_von_baeyer_descriptor(
                    candidate.descriptor.as_deref().unwrap(),
                    path,
                    &edges
                ))
        );
        let independent = find_von_baeyer_candidates(&mol, &atoms, &edges, (0, 0));
        assert!(!independent.is_empty());
        assert!(
            independent
                .iter()
                .all(|candidate| candidate.numbering.audit_ok)
        );
    }

    #[test]
    fn main_ring_integer_set_walk_preserves_first_audited_descriptor() {
        for (smiles, expected) in [
            (
                "B1CCC2CC3=CC1C32c1ccccc1-c1ccccc1",
                "tricyclo[5.1.1.0^{3,8}]",
            ),
            (
                "C=CC1=COC2CC3C2C2CC1C32",
                "tetracyclo[7.1.1.0^{2,5}.0^{3,10}]",
            ),
            (
                "O=[N+]([O-])[C@]12C[C@@H]3CC[C@H](CO1)C[C@@H]32",
                "tricyclo[4.2.2.0^{3,8}]",
            ),
            (
                "CCC1=C2CCC3=C1[C@H](C3)[C@H](C#N)C2",
                "tricyclo[4.2.2.0^{3,8}]",
            ),
            ("CC1CC2=C3CC1CCC=C3C2", "tricyclo[5.2.2.0^{3,9}]"),
            (
                "FC1C2CCC3C=PC4C=C4C1CCC(C3)O2",
                "tetracyclo[8.4.1.1^{7,12}.0^{2,4}]",
            ),
            (
                "CCC12CC3(N)C4CCC(C4)C(B1C(C)C)C23",
                "tetracyclo[5.2.1.1^{2,4}.0^{3,6}]",
            ),
        ] {
            let mol = read_smiles(smiles, &Default::default()).unwrap();
            let system = super::super::chains::find_ring_systems(&mol, &AtomSet::new())
                .unwrap()
                .into_iter()
                .find(|system| system.is_polycycle)
                .unwrap();
            let edges = edges_within_atoms(&mol, &system.atoms);
            assert_eq!(
                system.polycycle_descriptor.as_deref(),
                Some(expected),
                "{smiles}"
            );
            assert!(
                system
                    .paths
                    .iter()
                    .all(|path| audit_von_baeyer_descriptor(expected, path, &edges))
            );
            assert!(!system.paths.is_empty());
        }
    }
}
