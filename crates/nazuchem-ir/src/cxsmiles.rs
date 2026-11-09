// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | crates/nazuchem-ir/src/cxsmiles.rs | Phase B step 8a
// Copyright (c) 2026 Ameyanagi and Nazuchem contributors
// Changes: deterministic IR v1 SMILES traversal and conservative CXSMILES projection.

//! Deterministic, perception-free CXSMILES export.
//!
//! CXSMILES is a convenience view; versioned IR JSON is the lossless contract.
//! See `docs/cxsmiles.md` for standard extensions and the two locant properties.
//! Traversal and the compatibility projection follow the observable conventions
//! of the existing OPSIN writer, but this module is an original implementation:
//! no OPSIN implementation is included or copied here.

use crate::{
    BondOrder, Connectivity, Molecule, Ref, Relation, Stereo, StereoGroupKind, ValidationError,
    Winding,
};
use std::collections::{BTreeMap, VecDeque};
use std::fmt;

/// A structure validation or serialization error.
#[derive(Debug)]
#[non_exhaustive]
pub enum CxsmilesError {
    InvalidStructure(ValidationError),
    ContradictoryStereo,
    /// Standard CXSMILES cannot retain the stated stereo without changing it
    /// or asserting a configuration absent from the IR.
    StereoNotExpressible,
}

/// Compatibility spelling; new code should use [`CxsmilesError`].
pub type CxSmilesError = CxsmilesError;

impl fmt::Display for CxsmilesError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidStructure(error) => write!(f, "invalid structure: {error}"),
            Self::ContradictoryStereo => f.write_str("contradictory double bond stereochemistry"),
            Self::StereoNotExpressible => {
                f.write_str("stereochemistry cannot be expressed in standard CXSMILES")
            }
        }
    }
}

impl std::error::Error for CxsmilesError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::InvalidStructure(error) => Some(error),
            Self::ContradictoryStereo | Self::StereoNotExpressible => None,
        }
    }
}

impl From<ValidationError> for CxsmilesError {
    fn from(error: ValidationError) -> Self {
        Self::InvalidStructure(error)
    }
}

impl Molecule {
    /// Export a deterministic standard CXSMILES convenience view.
    ///
    /// First locants are atom values; custom properties retain all locants and
    /// parent membership. Use IR JSON to preserve every IR distinction.
    /// Unrepresentable stereo returns a typed error rather than changing it.
    pub fn to_cxsmiles(&self) -> Result<String, CxsmilesError> {
        Writer::new(self)?.write(false)
    }

    /// Export the legacy OPSIN semantic projection, without locants, hydrogen
    /// origin, explicit radical layers, or locant properties. Test support only.
    ///
    /// This is deterministic, not canonical. Atom order, ring numbering and
    /// globally equivalent slash directions can differ from other writers.
    #[cfg(any(test, feature = "test-support"))]
    #[doc(hidden)]
    pub fn semantic_cxsmiles_for_tests(&self) -> Result<String, CxsmilesError> {
        self.semantic_cxsmiles_with_atom_order_for_tests()
            .map(|(text, _)| text)
    }

    /// Test-only compatibility text and the original IR identities of its
    /// emitted atoms. Ordinary unannotated terminal H vertices may be folded
    /// into counts; the source molecule and public CXSMILES view are unchanged.
    #[cfg(any(test, feature = "test-support"))]
    #[doc(hidden)]
    pub fn semantic_cxsmiles_with_atom_order_for_tests(
        &self,
    ) -> Result<(String, Vec<u32>), CxsmilesError> {
        let (projection, identities) = compatibility_projection(self)?;
        let writer = Writer::new(&projection)?;
        let order = writer
            .order
            .iter()
            .map(|&atom| identities[atom as usize])
            .collect();
        Ok((writer.write(true)?, order))
    }

    /// The first locant of each principal-parent atom, in IR identity order.
    /// Atoms outside the principal parent have no displayed locant.
    pub fn displayed_locants(&self) -> Vec<Option<String>> {
        let mut parent = vec![false; self.atoms.len()];
        for &atom in &self.parent {
            if let Some(member) = parent.get_mut(atom as usize) {
                *member = true;
            }
        }
        self.atoms
            .iter()
            .zip(parent)
            .map(|(atom, is_parent)| is_parent.then(|| atom.locants.first().cloned()).flatten())
            .collect()
    }

    /// IR atom identities in the order used by the public CXSMILES projection.
    pub fn smiles_atom_order(&self) -> Result<Vec<u32>, CxsmilesError> {
        Ok(Writer::new(self)?.order)
    }
}

struct Writer<'a> {
    molecule: &'a Molecule,
    adjacent: Vec<Vec<usize>>,
    order: Vec<u32>,
    index: Vec<usize>,
    roots: Vec<u32>,
    parent_bond: Vec<Option<usize>>,
    children: Vec<Vec<(u32, usize)>>,
    rings: Vec<Vec<usize>>,
    bond_index: Vec<usize>,
    bond_to: Vec<u32>,
    directions: Vec<Option<bool>>,
    unknown_tetra_bonds: Vec<Option<usize>>,
}

impl<'a> Writer<'a> {
    fn new(molecule: &'a Molecule) -> Result<Self, CxsmilesError> {
        molecule.validate()?;
        let count = molecule.atoms.len();
        let mut writer = Self {
            molecule,
            adjacent: vec![Vec::new(); count],
            order: Vec::with_capacity(count),
            index: vec![usize::MAX; count],
            roots: Vec::new(),
            parent_bond: vec![None; count],
            children: vec![Vec::new(); count],
            rings: vec![Vec::new(); count],
            bond_index: vec![usize::MAX; molecule.bonds.len()],
            bond_to: vec![0; molecule.bonds.len()],
            directions: vec![None; molecule.bonds.len()],
            unknown_tetra_bonds: vec![None; count],
        };
        for (id, bond) in molecule.bonds.iter().enumerate() {
            for &atom in &bond.atoms {
                writer.adjacent[atom as usize].push(id);
            }
        }
        // Attachment points are component roots before ordinary atoms. Within
        // each component, bond insertion order is the deterministic tie-break.
        let roots = (0..count)
            .filter(|&atom| molecule.atoms[atom].element == 0)
            .chain((0..count).filter(|&atom| molecule.atoms[atom].element != 0));
        for root in roots {
            if writer.index[root] != usize::MAX {
                continue;
            }
            writer.roots.push(root as u32);
            writer.visit(root as u32);
        }
        // Ring markers themselves use original bond order, not discovery order.
        for rings in &mut writer.rings {
            rings.sort_unstable();
        }
        let mut next = 0;
        for &atom in &writer.order {
            if let Some(bond) = writer.parent_bond[atom as usize] {
                writer.bond_index[bond] = next;
                writer.bond_to[bond] = atom;
                next += 1;
            }
            for &bond in &writer.rings[atom as usize] {
                let other = writer.other(bond, atom);
                if writer.index[other as usize] < writer.index[atom as usize] {
                    writer.bond_index[bond] = next;
                    writer.bond_to[bond] = atom;
                    next += 1;
                }
            }
        }
        Ok(writer)
    }

    fn other(&self, bond: usize, atom: u32) -> u32 {
        let pair = self.molecule.bonds[bond].atoms;
        if pair[0] == atom { pair[1] } else { pair[0] }
    }

    fn visit(&mut self, root: u32) {
        self.index[root as usize] = self.order.len();
        self.order.push(root);
        let mut stack = vec![(root, 0)];
        while let Some((atom, cursor)) = stack.last_mut() {
            if *cursor == self.adjacent[*atom as usize].len() {
                stack.pop();
                continue;
            }
            let bond = self.adjacent[*atom as usize][*cursor];
            *cursor += 1;
            if self.parent_bond[*atom as usize] == Some(bond) {
                continue;
            }
            let atom = *atom;
            let other = self.other(bond, atom);
            if self.index[other as usize] == usize::MAX {
                self.parent_bond[other as usize] = Some(bond);
                self.children[atom as usize].push((other, bond));
                self.index[other as usize] = self.order.len();
                self.order.push(other);
                stack.push((other, 0));
            } else if self.index[other as usize] < self.index[atom as usize] {
                self.rings[atom as usize].push(bond);
                self.rings[other as usize].push(bond);
            }
        }
    }

    fn find_bond(&self, first: u32, second: u32) -> Option<usize> {
        self.adjacent[first as usize]
            .iter()
            .copied()
            .find(|&bond| self.other(bond, first) == second)
    }

    fn validate_projection(&self) -> Result<(), CxsmilesError> {
        for stereo in &self.molecule.stereo {
            match stereo {
                Stereo::Tetrahedral {
                    center, winding, ..
                } if *winding != Winding::Unknown && self.parity(*center).is_empty() => {
                    return Err(CxsmilesError::StereoNotExpressible);
                }
                // No portable standard CXSMILES encoding is implemented for
                // ordered axial references, including explicit unknown axes.
                // JSON retains the exact IR record.
                Stereo::Atropisomer { .. } => {
                    return Err(CxsmilesError::StereoNotExpressible);
                }
                _ => {}
            }
        }
        for group in &self.molecule.stereo_groups {
            if group.stereo.iter().any(|&entry| {
                !matches!(
                    self.molecule.stereo[entry as usize],
                    Stereo::Tetrahedral { winding, .. } if winding != Winding::Unknown
                )
            }) {
                // CX enhanced groups identify tetrahedral atom centres, not
                // arbitrary stereo records or unknown configurations.
                return Err(CxsmilesError::StereoNotExpressible);
            }
        }
        Ok(())
    }

    fn assign_unknown_tetra_bonds(&mut self) -> Result<(), CxsmilesError> {
        // A CX wavy bond can identify only one unknown tetrahedral centre.
        // Find an injective centre-to-bond matching instead of choosing each
        // centre's first bond independently. Alternating paths can move a
        // previously assigned centre to another carrier when needed.
        self.unknown_tetra_bonds.fill(None);
        let mut owner: Vec<Option<usize>> = vec![None; self.molecule.bonds.len()];
        for stereo in &self.molecule.stereo {
            let Stereo::Tetrahedral {
                center: root,
                winding: Winding::Unknown,
                ..
            } = stereo
            else {
                continue;
            };
            let root = *root as usize;
            let mut seen = vec![false; self.molecule.atoms.len()];
            let mut previous = vec![None; self.molecule.atoms.len()];
            let mut queue = VecDeque::from([root]);
            seen[root] = true;
            let mut free = None;
            'search: while let Some(center) = queue.pop_front() {
                for &bond in &self.adjacent[center] {
                    if self.molecule.bonds[bond].order != BondOrder::Single
                        || self.directions[bond].is_some()
                    {
                        continue;
                    }
                    if let Some(other) = owner[bond] {
                        if !seen[other] {
                            seen[other] = true;
                            previous[other] = Some((center, bond));
                            queue.push_back(other);
                        }
                    } else {
                        free = Some((center, bond));
                        break 'search;
                    }
                }
            }
            let (mut center, mut bond) = free.ok_or(CxsmilesError::StereoNotExpressible)?;
            loop {
                owner[bond] = Some(center);
                self.unknown_tetra_bonds[center] = Some(bond);
                match previous[center] {
                    Some((prior, contested)) => {
                        center = prior;
                        bond = contested;
                    }
                    None => break,
                }
            }
        }
        Ok(())
    }

    fn assign_directions(&mut self) -> Result<(), CxsmilesError> {
        // Slash directions form binary parity constraints. Solving the whole
        // connected constraint component also handles conjugated systems.
        let mut equations = vec![Vec::new(); self.molecule.bonds.len()];
        let mut seeds = Vec::new();
        let mut specified = vec![false; self.molecule.bonds.len()];
        let mut entries = self
            .molecule
            .stereo
            .iter()
            .filter_map(|stereo| match stereo {
                Stereo::DoubleBond {
                    bond,
                    refs,
                    relation,
                } if *relation != Relation::Unknown => Some((*bond as usize, *refs, *relation)),
                _ => None,
            })
            .collect::<Vec<_>>();
        entries.sort_by_key(|(bond, _, _)| self.bond_index[*bond]);
        for (central, refs, relation) in entries {
            specified[central] = true;
            let endpoints = self.molecule.bonds[central].atoms;
            let first = self.find_bond(refs[0], endpoints[0]).unwrap();
            let second = self.find_bond(endpoints[1], refs[1]).unwrap();
            if [first, second]
                .iter()
                .any(|&bond| self.molecule.bonds[bond].order != BondOrder::Single)
            {
                return Err(CxsmilesError::StereoNotExpressible);
            }
            let first_flip = self.bond_to[first] != endpoints[0];
            let second_flip = self.bond_to[second] != refs[1];
            let parity = (relation == Relation::Cis) ^ first_flip ^ second_flip;
            equations[first].push((second, parity));
            equations[second].push((first, parity));
            seeds.push((first, first_flip));
            for (endpoint, primary) in [(endpoints[0], first), (endpoints[1], second)] {
                let others = self.adjacent[endpoint as usize]
                    .iter()
                    .copied()
                    .filter(|&bond| bond != central && bond != primary)
                    .collect::<Vec<_>>();
                if others.len() > 1 {
                    return Err(CxsmilesError::StereoNotExpressible);
                }
                if let [other] = others.as_slice()
                    && self.molecule.bonds[*other].order == BondOrder::Single
                {
                    let opposite = true
                        ^ (self.bond_to[primary] == endpoint)
                        ^ (self.bond_to[*other] == endpoint);
                    equations[primary].push((*other, opposite));
                    equations[*other].push((primary, opposite));
                }
            }
        }
        for (seed, direction) in seeds {
            if self.directions[seed].is_some() {
                continue;
            }
            self.directions[seed] = Some(direction);
            let mut queue = VecDeque::from([seed]);
            while let Some(bond) = queue.pop_front() {
                let direction = self.directions[bond].unwrap();
                for &(other, parity) in &equations[bond] {
                    let expected = direction ^ parity;
                    match self.directions[other] {
                        Some(existing) if existing != expected => {
                            return Err(CxsmilesError::ContradictoryStereo);
                        }
                        Some(_) => {}
                        None => {
                            self.directions[other] = Some(expected);
                            queue.push_back(other);
                        }
                    }
                }
            }
        }
        let original = self.directions.clone();
        // Avoid accidentally specifying an otherwise unspecified double bond
        // using redundant marks introduced by neighbouring stereocentres.
        let mut bonds = (0..self.molecule.bonds.len()).collect::<Vec<_>>();
        bonds.sort_by_key(|&bond| self.bond_index[bond]);
        for central in bonds {
            if self.molecule.bonds[central].order != BondOrder::Double || specified[central] {
                continue;
            }
            let marked = self.molecule.bonds[central].atoms.map(|atom| {
                self.adjacent[atom as usize]
                    .iter()
                    .copied()
                    .filter(|&bond| self.directions[bond].is_some())
                    .collect::<Vec<_>>()
            });
            if marked.iter().any(Vec::is_empty) {
                continue;
            }
            for side in &marked {
                let removable = side.iter().all(|&bond| {
                    self.molecule.bonds[bond].atoms.iter().all(|&atom| {
                        !self.adjacent[atom as usize]
                            .iter()
                            .any(|&edge| specified[edge])
                            || self.adjacent[atom as usize].iter().any(|&edge| {
                                !side.contains(&edge) && self.directions[edge].is_some()
                            })
                    })
                });
                if removable {
                    for &bond in side {
                        self.directions[bond] = None;
                    }
                    break;
                }
            }
        }
        if self.slashes_match(&specified) {
            return Ok(());
        }
        // Greedy pruning keeps the established spelling when it works. When
        // its local choices fail, solve the complete mark-presence problem:
        // every known double bond needs a mark at each end, and every absent
        // or unknown double bond must have an unmarked end. Each known end
        // has at most two single bonds, so these are exactly 2-SAT clauses.
        self.directions = original;
        let mut clauses = Vec::new();
        for (central, bond) in self.molecule.bonds.iter().enumerate() {
            if bond.order != BondOrder::Double {
                continue;
            }
            let marked = bond.atoms.map(|atom| self.marked_at(atom));
            if specified[central] {
                for side in marked {
                    match side.as_slice() {
                        [first] => clauses.push(((*first, true), (*first, true))),
                        [first, second] => clauses.push(((*first, true), (*second, true))),
                        _ => return Err(CxsmilesError::StereoNotExpressible),
                    }
                }
            } else {
                for &first in &marked[0] {
                    for &second in &marked[1] {
                        clauses.push(((first, false), (second, false)));
                    }
                }
            }
        }
        let chosen = solve_marks(self.directions.len(), &clauses)
            .ok_or(CxsmilesError::StereoNotExpressible)?;
        for (direction, selected) in self.directions.iter_mut().zip(chosen) {
            if !selected {
                *direction = None;
            }
        }
        if !self.slashes_match(&specified) {
            return Err(CxsmilesError::StereoNotExpressible);
        }
        Ok(())
    }

    fn marked_at(&self, atom: u32) -> Vec<usize> {
        self.adjacent[atom as usize]
            .iter()
            .copied()
            .filter(|&bond| self.directions[bond].is_some())
            .collect()
    }

    fn slashes_match(&self, specified: &[bool]) -> bool {
        self.molecule.bonds.iter().enumerate().all(|(id, bond)| {
            bond.order != BondOrder::Double
                || bond
                    .atoms
                    .iter()
                    .all(|&atom| !self.marked_at(atom).is_empty())
                    == specified[id]
        })
    }

    fn bond_text(&self, bond: usize) -> &'static str {
        match self.molecule.bonds[bond].order {
            BondOrder::Double => "=",
            BondOrder::Triple => "#",
            BondOrder::Quadruple => "$",
            BondOrder::Single => match self.directions[bond] {
                Some(true) => "/",
                Some(false) => "\\",
                None => "",
            },
        }
    }

    fn write(mut self, semantic: bool) -> Result<String, CxsmilesError> {
        self.assign_directions()?;
        self.assign_unknown_tetra_bonds()?;
        self.validate_projection()?;
        let mut output = String::new();
        let mut closure_ids = vec![None; self.molecule.bonds.len()];
        let mut available = (1..=99).chain(std::iter::once(0)).collect::<VecDeque<_>>();
        let mut extra_closure = 100;
        enum Event {
            Atom(u32),
            Open,
            Close,
            Bond(usize),
        }
        for &root in &self.roots {
            if !output.is_empty() {
                output.push('.');
            }
            let mut events = vec![Event::Atom(root)];
            while let Some(event) = events.pop() {
                let atom = match event {
                    Event::Open => {
                        output.push('(');
                        continue;
                    }
                    Event::Close => {
                        output.push(')');
                        continue;
                    }
                    Event::Bond(bond) => {
                        output.push_str(self.bond_text(bond));
                        continue;
                    }
                    Event::Atom(atom) => atom,
                };
                output.push_str(&self.atom_text(atom, semantic));
                let mut released = Vec::new();
                for &bond in &self.rings[atom as usize] {
                    let other = self.other(bond, atom);
                    if self.index[other as usize] < self.index[atom as usize] {
                        let number = closure_ids[bond].unwrap();
                        output.push_str(&ring_number(number));
                        released.push(number);
                    }
                }
                for &bond in &self.rings[atom as usize] {
                    let other = self.other(bond, atom);
                    if self.index[other as usize] > self.index[atom as usize] {
                        let number = available.pop_front().unwrap_or_else(|| {
                            let number = extra_closure;
                            extra_closure += 1;
                            number
                        });
                        closure_ids[bond] = Some(number);
                        output.push_str(self.bond_text(bond));
                        output.push_str(&ring_number(number));
                    }
                }
                for number in released.into_iter().rev() {
                    available.push_front(number);
                }
                for (index, &(child, bond)) in self.children[atom as usize].iter().enumerate().rev()
                {
                    let branch = index + 1 < self.children[atom as usize].len();
                    if branch {
                        events.push(Event::Close);
                    }
                    events.push(Event::Atom(child));
                    events.push(Event::Bond(bond));
                    if branch {
                        events.push(Event::Open);
                    }
                }
            }
        }
        let layers = self.layers(semantic)?;
        if !layers.is_empty() {
            output.push_str(" |");
            output.push_str(&layers.join(","));
            output.push('|');
        }
        Ok(output)
    }

    fn tetrahedral(&self, atom: u32) -> Option<(&[Ref; 4], Winding)> {
        self.molecule.stereo.iter().find_map(|stereo| match stereo {
            Stereo::Tetrahedral {
                center,
                neighbors,
                winding,
            } if *center == atom => Some((neighbors, *winding)),
            _ => None,
        })
    }

    fn parity(&self, atom: u32) -> &'static str {
        let Some((stored, winding)) = self.tetrahedral(atom) else {
            return "";
        };
        if winding == Winding::Unknown {
            return "";
        }
        let mut current = Vec::with_capacity(4);
        if let Some(parent) = self.parent_bond[atom as usize] {
            current.push(Ref::Atom(self.other(parent, atom)));
        }
        // A virtual ligand is first at a root, or follows the incoming atom.
        // SMILES has only one implicit ligand slot; preflight rejects known
        // stereo records requiring more than one virtual ligand.
        let virtuals = stored
            .iter()
            .filter(|reference| !matches!(reference, Ref::Atom(_)));
        for reference in virtuals {
            current.push(*reference);
        }
        for &bond in &self.rings[atom as usize] {
            if self.index[self.other(bond, atom) as usize] < self.index[atom as usize] {
                current.push(Ref::Atom(self.other(bond, atom)));
            }
        }
        for &bond in &self.rings[atom as usize] {
            if self.index[self.other(bond, atom) as usize] > self.index[atom as usize] {
                current.push(Ref::Atom(self.other(bond, atom)));
            }
        }
        for &(child, _) in &self.children[atom as usize] {
            current.push(Ref::Atom(child));
        }
        if current.len() != 4
            || stored
                .iter()
                .filter(|reference| !matches!(reference, Ref::Atom(_)))
                .count()
                > 1
        {
            return "";
        }
        let permutation = stored.map(|reference| {
            current
                .iter()
                .position(|candidate| *candidate == reference)
                .unwrap()
        });
        let odd = (0..4)
            .flat_map(|first| (first + 1..4).map(move |second| (first, second)))
            .filter(|&(first, second)| permutation[first] > permutation[second])
            .count()
            % 2
            == 1;
        if (winding == Winding::Clockwise) ^ odd {
            "@@"
        } else {
            "@"
        }
    }

    fn atom_text(&self, atom: u32, semantic: bool) -> String {
        let a = &self.molecule.atoms[atom as usize];
        let bonds = self.adjacent[atom as usize]
            .iter()
            .map(|&bond| order_value(self.molecule.bonds[bond].order))
            .sum::<u32>();
        let total = bonds + u32::from(a.hydrogens);
        let expected: &[u32] = match a.element {
            5 => &[3],
            6 => &[4],
            7 | 15 => &[3, 5],
            8 => &[2],
            16 => &[2, 4, 6],
            9 | 17 | 35 | 53 => &[1],
            0 => &[1, 2, 3, 4, 5, 6, 7, 8, 9],
            _ => &[],
        };
        let implicit = expected
            .iter()
            .copied()
            .find(|&valence| valence >= bonds)
            .unwrap_or(bonds);
        let described =
            expected.contains(&total) || expected.last().is_some_and(|&last| total > last);
        let parity = self.parity(atom);
        let atom_map = if semantic {
            self.polymer_class(atom)
                .map(|value| value as u32)
                .or(a.atom_map)
        } else {
            a.atom_map
        };
        let brackets = atom_map.is_some()
            || expected.is_empty()
            || !described
            || total != implicit
            || a.charge != 0
            || a.isotope.is_some()
            || !parity.is_empty()
            || (!semantic && (a.hydrogens_stated || a.radical_electrons > 0));
        let mut text = String::new();
        if brackets {
            text.push('[');
        }
        if let Some(isotope) = a.isotope {
            text.push_str(&isotope.to_string());
        }
        text.push_str(ELEMENTS[a.element as usize]);
        text.push_str(parity);
        if brackets && a.hydrogens > 0 {
            text.push('H');
            if a.hydrogens != 1 {
                text.push_str(&a.hydrogens.to_string());
            }
        }
        match a.charge {
            0 => {}
            1 => text.push('+'),
            -1 => text.push('-'),
            charge => {
                if charge > 0 {
                    text.push('+');
                }
                text.push_str(&charge.to_string());
            }
        }
        if let Some(class) = atom_map {
            text.push(':');
            text.push_str(&class.to_string());
        }
        if brackets {
            text.push(']');
        }
        text
    }

    fn polymer_class(&self, atom: u32) -> Option<usize> {
        let mut unique = Vec::new();
        for &point in self
            .molecule
            .polymers
            .iter()
            .flat_map(|unit| &unit.attachment_points)
        {
            if !unique.contains(&point) {
                unique.push(point);
            }
        }
        unique
            .iter()
            .position(|&point| point == atom)
            .map(|index| index + 1)
    }

    fn label(&self, atom: u32) -> String {
        let atom_data = &self.molecule.atoms[atom as usize];
        if let Some(label) = &atom_data.label {
            escape_label(label)
        } else if self.polymer_class(atom).is_some() {
            "star_e".into()
        } else if let Some(point) = atom_data.attachment_point {
            format!("_AP{point}")
        } else {
            String::new()
        }
    }

    fn layers(&self, semantic: bool) -> Result<Vec<String>, CxsmilesError> {
        let mut layers = Vec::new();
        let labels = self
            .order
            .iter()
            .map(|&id| self.label(id))
            .collect::<Vec<_>>();
        if let Some(last) = labels.iter().rposition(|label| !label.is_empty()) {
            layers.push(format!("${}$", labels[..=last].join(";")));
        }
        if !semantic {
            let values = self
                .order
                .iter()
                .map(|&id| {
                    self.molecule.atoms[id as usize]
                        .locants
                        .first()
                        .map(|value| escape(value))
                        .unwrap_or_default()
                })
                .collect::<Vec<_>>();
            if let Some(last) = values.iter().rposition(|value| !value.is_empty()) {
                layers.push(format!("$_AV:{}$", values[..=last].join(";")));
            }
            let mut radicals = BTreeMap::<u8, Vec<usize>>::new();
            for (index, &atom) in self.order.iter().enumerate() {
                match self.molecule.atoms[atom as usize].radical_electrons {
                    1 => radicals.entry(1).or_default().push(index),
                    2 => radicals.entry(2).or_default().push(index),
                    3 => radicals.entry(5).or_default().push(index),
                    _ => {}
                }
            }
            for (kind, atoms) in radicals {
                layers.push(format!("^{kind}:{}", indices(atoms)));
            }
        }
        let mut groups = self
            .molecule
            .stereo_groups
            .iter()
            .map(|group| {
                let mut atoms = group
                    .stereo
                    .iter()
                    .filter_map(|&entry| match &self.molecule.stereo[entry as usize] {
                        Stereo::Tetrahedral {
                            center, winding, ..
                        } if *winding != Winding::Unknown => Some(self.index[*center as usize]),
                        _ => None,
                    })
                    .collect::<Vec<_>>();
                atoms.sort_unstable();
                (group.kind, atoms)
            })
            .filter(|(_, atoms)| !atoms.is_empty())
            .collect::<Vec<_>>();
        groups.sort_by(|first, second| first.1.cmp(&second.1));
        if semantic && groups.len() == 1 {
            match groups[0].0 {
                StereoGroupKind::And(1 | 2) => layers.push("r".into()),
                StereoGroupKind::And(number) => layers.push(format!(
                    "&{number}:{}",
                    indices(groups[0].1.iter().copied())
                )),
                StereoGroupKind::Or(number) => layers.push(format!(
                    "o{number}:{}",
                    indices(groups[0].1.iter().copied())
                )),
            }
        } else {
            let (mut and, mut or) = (1, 1);
            for (kind, atoms) in groups {
                let prefix = match kind {
                    StereoGroupKind::And(number) => {
                        let number = if semantic {
                            let value = and;
                            and += 1;
                            value
                        } else {
                            number
                        };
                        format!("&{number}")
                    }
                    StereoGroupKind::Or(number) => {
                        let number = if semantic {
                            let value = or;
                            or += 1;
                            value
                        } else {
                            number
                        };
                        format!("o{number}")
                    }
                };
                layers.push(format!("{prefix}:{}", indices(atoms)));
            }
        }
        let mut positions = Vec::new();
        for position in &self.molecule.position_variation {
            let bond = &self.molecule.bonds[position.bond as usize];
            // Bond orientation identifies the source explicitly, including
            // non-wildcard sources and bonds with two wildcard endpoints.
            let attachment = bond.atoms[0];
            positions.push(format!(
                "{}:{}",
                self.index[attachment as usize],
                position
                    .endpoints
                    .iter()
                    .map(|&atom| self.index[atom as usize].to_string())
                    .collect::<Vec<_>>()
                    .join(".")
            ));
        }
        if !positions.is_empty() {
            layers.push(format!("m:{}", positions.join(",")));
        }
        for polymer in &self.molecule.polymers {
            let atoms = polymer.atoms.iter().map(|&atom| self.index[atom as usize]);
            let connectivity = match polymer.connectivity {
                Connectivity::HeadToTail => "ht",
                Connectivity::HeadToHead => "hh",
                Connectivity::Unknown => "eu",
            };
            let label = polymer
                .label
                .as_deref()
                .map(escape_label)
                .unwrap_or_default();
            layers.push(format!("Sg:n:{}:{label}:{connectivity}", indices(atoms)));
        }
        if !semantic {
            let mut wiggly = Vec::new();
            let mut unknown_double = Vec::new();
            for stereo in &self.molecule.stereo {
                match stereo {
                    Stereo::Tetrahedral {
                        center,
                        winding: Winding::Unknown,
                        ..
                    } => {
                        let bond = self.unknown_tetra_bonds[*center as usize]
                            .expect("unknown tetrahedral carriers are assigned before writing");
                        wiggly.push(format!(
                            "{}.{}",
                            self.index[*center as usize], self.bond_index[bond]
                        ));
                    }
                    Stereo::DoubleBond {
                        bond,
                        relation: Relation::Unknown,
                        ..
                    } => {
                        let bond = *bond as usize;
                        wiggly.push(format!(
                            "{}.{}",
                            self.index[self.molecule.bonds[bond].atoms[0] as usize],
                            self.bond_index[bond]
                        ));
                        unknown_double.push(self.bond_index[bond]);
                    }
                    _ => {}
                }
            }
            if !wiggly.is_empty() {
                layers.push(format!("w:{}", wiggly.join(",")));
            }
            if !unknown_double.is_empty() {
                layers.push(format!("ctu:{}", indices(unknown_double)));
            }
            self.metadata(&mut layers);
        }
        Ok(layers)
    }

    fn metadata(&self, layers: &mut Vec<String>) {
        let mut properties = Vec::new();
        for (index, &id) in self.order.iter().enumerate() {
            let atom = &self.molecule.atoms[id as usize];
            if self.molecule.parent.contains(&id) {
                properties.push(format!("{index}.nazuchemParent.1"));
            }
            if !atom.locants.is_empty() {
                // Serializing a string vector cannot fail: it has no map keys,
                // floats, or custom serializers.
                let locants = serde_json::to_string(&atom.locants)
                    .expect("serializing a string vector is infallible");
                properties.push(format!("{index}.nazuchemLocants.{}", escape(&locants)));
            }
        }
        if !properties.is_empty() {
            layers.push(format!("atomProp:{}", properties.join(":")));
        }
    }
}

/// Prepare the legacy oracle's hydrogen-count view without changing the IR.
/// Annotation-bearing vertices and bonds keep their identities. Only a
/// tetrahedral atom reference can become an implicit-H reference; other stereo
/// kinds keep their explicit references because their wire shapes differ.
#[cfg(any(test, feature = "test-support"))]
fn compatibility_projection(molecule: &Molecule) -> Result<(Molecule, Vec<u32>), CxsmilesError> {
    let writer = Writer::new(molecule)?;
    let mut protected = vec![false; molecule.atoms.len()];
    for &atom in &molecule.parent {
        protected[atom as usize] = true;
    }
    for stereo in &molecule.stereo {
        match stereo {
            Stereo::Tetrahedral {
                center, neighbors, ..
            } => {
                protected[*center as usize] = true;
                if neighbors.iter().any(|r| !matches!(r, Ref::Atom(_))) {
                    for reference in neighbors {
                        if let Ref::Atom(atom) = reference {
                            protected[*atom as usize] = true;
                        }
                    }
                }
            }
            Stereo::DoubleBond { bond, refs, .. } => {
                for atom in molecule.bonds[*bond as usize]
                    .atoms
                    .into_iter()
                    .chain(*refs)
                {
                    protected[atom as usize] = true;
                }
            }
            Stereo::Atropisomer { bond, refs, .. } => {
                for atom in molecule.bonds[*bond as usize].atoms {
                    protected[atom as usize] = true;
                }
                for reference in refs {
                    if let Ref::Atom(atom) = reference {
                        protected[*atom as usize] = true;
                    }
                }
            }
        }
    }
    for position in &molecule.position_variation {
        for &atom in molecule.bonds[position.bond as usize]
            .atoms
            .iter()
            .chain(&position.endpoints)
        {
            protected[atom as usize] = true;
        }
    }
    for unit in &molecule.polymers {
        for &atom in unit.atoms.iter().chain(&unit.attachment_points) {
            protected[atom as usize] = true;
        }
        for &bond in &unit.crossing_bonds {
            for atom in molecule.bonds[bond as usize].atoms {
                protected[atom as usize] = true;
            }
        }
    }
    let mut hydrogen_counts = molecule
        .atoms
        .iter()
        .map(|atom| atom.hydrogens)
        .collect::<Vec<_>>();
    let mut folded = vec![false; molecule.atoms.len()];
    for (id, atom) in molecule.atoms.iter().enumerate() {
        if protected[id]
            || atom.element != 1
            || atom.isotope.is_some()
            || atom.charge != 0
            || atom.hydrogens != 0
            || atom.radical_electrons != 0
            || atom.label.is_some()
            || atom.atom_map.is_some()
            || atom.attachment_point.is_some()
            || !atom.locants.is_empty()
        {
            continue;
        }
        let [bond] = writer.adjacent[id].as_slice() else {
            continue;
        };
        if molecule.bonds[*bond].order != BondOrder::Single {
            continue;
        }
        let neighbor = writer.other(*bond, id as u32) as usize;
        if matches!(molecule.atoms[neighbor].element, 0 | 1)
            || hydrogen_counts[neighbor] == u8::MAX
            || (writer.tetrahedral(neighbor as u32).is_some() && hydrogen_counts[neighbor] != 0)
        {
            continue;
        }
        folded[id] = true;
        hydrogen_counts[neighbor] += 1;
    }
    let identities = (0..molecule.atoms.len())
        .filter(|&atom| !folded[atom])
        .map(|atom| atom as u32)
        .collect::<Vec<_>>();
    let mut atom_indices = vec![0; molecule.atoms.len()];
    let atoms = identities
        .iter()
        .enumerate()
        .map(|(index, &id)| {
            atom_indices[id as usize] = index as u32;
            let mut atom = molecule.atoms[id as usize].clone();
            atom.hydrogens = hydrogen_counts[id as usize];
            atom
        })
        .collect();
    let mut bond_indices = vec![0; molecule.bonds.len()];
    let bonds = molecule
        .bonds
        .iter()
        .enumerate()
        .filter_map(|(id, bond)| {
            if bond.atoms.iter().any(|&atom| folded[atom as usize]) {
                None
            } else {
                Some((id, bond))
            }
        })
        .enumerate()
        .map(|(index, (id, bond))| {
            bond_indices[id] = index as u32;
            crate::Bond {
                atoms: bond.atoms.map(|atom| atom_indices[atom as usize]),
                order: bond.order,
            }
        })
        .collect();
    let reference = |value: Ref| match value {
        Ref::Atom(atom) if folded[atom as usize] => Ref::ImplicitHydrogen,
        Ref::Atom(atom) => Ref::Atom(atom_indices[atom as usize]),
        other => other,
    };
    let stereo = molecule
        .stereo
        .iter()
        .map(|stereo| match stereo {
            Stereo::Tetrahedral {
                center,
                neighbors,
                winding,
            } => Stereo::Tetrahedral {
                center: atom_indices[*center as usize],
                neighbors: neighbors.map(reference),
                winding: *winding,
            },
            Stereo::DoubleBond {
                bond,
                refs,
                relation,
            } => Stereo::DoubleBond {
                bond: bond_indices[*bond as usize],
                refs: refs.map(|atom| atom_indices[atom as usize]),
                relation: *relation,
            },
            Stereo::Atropisomer {
                bond,
                refs,
                winding,
            } => Stereo::Atropisomer {
                bond: bond_indices[*bond as usize],
                refs: refs.map(reference),
                winding: *winding,
            },
        })
        .collect();
    let position_variation = molecule
        .position_variation
        .iter()
        .map(|position| crate::PositionVariation {
            bond: bond_indices[position.bond as usize],
            endpoints: position
                .endpoints
                .iter()
                .map(|&atom| atom_indices[atom as usize])
                .collect(),
        })
        .collect();
    let polymers = molecule
        .polymers
        .iter()
        .map(|unit| crate::RepeatUnit {
            atoms: unit
                .atoms
                .iter()
                .map(|&atom| atom_indices[atom as usize])
                .collect(),
            crossing_bonds: unit
                .crossing_bonds
                .iter()
                .map(|&bond| bond_indices[bond as usize])
                .collect(),
            attachment_points: unit
                .attachment_points
                .iter()
                .map(|&atom| atom_indices[atom as usize])
                .collect(),
            connectivity: unit.connectivity,
            label: unit.label.clone(),
        })
        .collect();
    Ok((
        Molecule {
            atoms,
            bonds,
            stereo,
            position_variation,
            polymers,
            stereo_groups: molecule.stereo_groups.clone(),
            parent: molecule
                .parent
                .iter()
                .map(|&atom| atom_indices[atom as usize])
                .collect(),
        },
        identities,
    ))
}

type MarkLiteral = (usize, bool);

/// Solve mark-presence clauses without recursion, including very large graphs.
fn solve_marks(count: usize, clauses: &[(MarkLiteral, MarkLiteral)]) -> Option<Vec<bool>> {
    let mut graph = vec![Vec::new(); 2 * count];
    let mut reverse = vec![Vec::new(); 2 * count];
    for &((first, first_value), (second, second_value)) in clauses {
        let first = 2 * first + usize::from(first_value);
        let second = 2 * second + usize::from(second_value);
        for (from, to) in [(first ^ 1, second), (second ^ 1, first)] {
            graph[from].push(to);
            reverse[to].push(from);
        }
    }
    let mut seen = vec![false; graph.len()];
    let mut order = Vec::with_capacity(graph.len());
    for root in 0..graph.len() {
        if seen[root] {
            continue;
        }
        seen[root] = true;
        let mut stack = vec![(root, 0)];
        while let Some((node, cursor)) = stack.last_mut() {
            if *cursor == graph[*node].len() {
                order.push(*node);
                stack.pop();
            } else {
                let next = graph[*node][*cursor];
                *cursor += 1;
                if !seen[next] {
                    seen[next] = true;
                    stack.push((next, 0));
                }
            }
        }
    }
    let mut component = vec![usize::MAX; graph.len()];
    let mut next_component = 0;
    for &root in order.iter().rev() {
        if component[root] != usize::MAX {
            continue;
        }
        component[root] = next_component;
        let mut stack = vec![root];
        while let Some(node) = stack.pop() {
            for &next in &reverse[node] {
                if component[next] == usize::MAX {
                    component[next] = next_component;
                    stack.push(next);
                }
            }
        }
        next_component += 1;
    }
    let mut selected = Vec::with_capacity(count);
    for bond in 0..count {
        let (no, yes) = (component[2 * bond], component[2 * bond + 1]);
        if no == yes {
            return None;
        }
        selected.push(yes > no);
    }
    Some(selected)
}

fn order_value(order: BondOrder) -> u32 {
    match order {
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
    }
}

fn ring_number(number: usize) -> String {
    match number {
        0..=9 => number.to_string(),
        10..=99 => format!("%{number}"),
        _ => format!("%({number})"),
    }
}

fn indices(values: impl IntoIterator<Item = usize>) -> String {
    values
        .into_iter()
        .map(|value| value.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

fn escape(value: &str) -> String {
    // UTF-16 numeric entities match the established OPSIN label convention,
    // including supplementary characters represented by surrogate pairs.
    value
        .encode_utf16()
        .map(|unit| {
            if unit <= 127 && (unit as u8).is_ascii_alphanumeric() {
                (unit as u8 as char).to_string()
            } else {
                format!("&#{unit};")
            }
        })
        .collect()
}

const ELEMENTS: [&str; 119] = [
    "*", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
    "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
    "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
    "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
    "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
    "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
];

fn escape_label(value: &str) -> String {
    if (value.strip_prefix("_AP").is_some_and(|suffix| {
        !suffix.is_empty() && suffix.bytes().all(|byte| byte.is_ascii_digit())
    })) || value == "star_e"
    {
        value.to_owned()
    } else if let Some(prefix) = value.strip_suffix("_p") {
        format!("{}_p", escape(prefix))
    } else {
        escape(value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Atom, AxialWinding, Bond, PositionVariation, RepeatUnit, StereoGroup};

    fn atom(element: u8, hydrogens: u8) -> Atom {
        Atom {
            element,
            hydrogens,
            ..Atom::default()
        }
    }

    fn bond(first: u32, second: u32, order: BondOrder) -> Bond {
        Bond {
            atoms: [first, second],
            order,
        }
    }

    fn polyene(states: &[u8], reversed: bool) -> Molecule {
        let count = 2 * states.len() + 2;
        let mut molecule = Molecule {
            atoms: (0..count)
                .map(|i| atom(6, if i == 0 || i + 1 == count { 3 } else { 1 }))
                .collect(),
            bonds: (0..count - 1)
                .map(|i| {
                    let order = if i % 2 == 0 {
                        BondOrder::Single
                    } else {
                        BondOrder::Double
                    };
                    if reversed {
                        bond(i as u32 + 1, i as u32, order)
                    } else {
                        bond(i as u32, i as u32 + 1, order)
                    }
                })
                .collect(),
            ..Molecule::default()
        };
        for (i, &state) in states.iter().enumerate() {
            let relation = match state {
                0 => continue,
                1 => Relation::Cis,
                2 => Relation::Trans,
                3 => Relation::Unknown,
                _ => unreachable!(),
            };
            let refs = [2 * i as u32, 2 * i as u32 + 3];
            molecule.stereo.push(Stereo::DoubleBond {
                bond: 2 * i as u32 + 1,
                refs: if reversed { [refs[1], refs[0]] } else { refs },
                relation,
            });
        }
        molecule
    }

    #[test]
    fn convenience_projection_has_only_the_two_locant_properties() {
        let mut molecule = Molecule {
            atoms: vec![atom(6, 3), atom(0, 0), atom(8, 1)],
            bonds: vec![bond(0, 2, BondOrder::Single), bond(0, 1, BondOrder::Single)],
            parent: vec![0],
            ..Molecule::default()
        };
        molecule.atoms[0].hydrogens = 2;
        molecule.atoms[0].locants = vec!["1".into(), "1′".into(), "".into()];
        molecule.atoms[2].locants = vec!["outside parent".into()];
        let output = molecule.to_cxsmiles().unwrap();
        assert!(output.contains("$_AV:;1;outside&#32;parent$"), "{output}");
        assert!(output.contains("1.nazuchemParent.1"));
        assert!(output.contains("1.nazuchemLocants."));
        assert!(output.contains("2.nazuchemLocants."));
        let properties = output
            .split_once("atomProp:")
            .unwrap()
            .1
            .trim_end_matches('|');
        for property in properties.split(':') {
            let key = property.split('.').nth(1).unwrap();
            assert!(
                ["nazuchemParent", "nazuchemLocants"].contains(&key),
                "{property}"
            );
        }
        assert_eq!(molecule.smiles_atom_order().unwrap(), [1, 0, 2]);
        assert_eq!(molecule.displayed_locants(), [Some("1".into()), None, None]);
        assert_eq!(output, molecule.to_cxsmiles().unwrap());
        assert_eq!(Molecule::default().to_cxsmiles().unwrap(), "");
    }

    #[test]
    fn atom_maps_attachment_roles_and_labels_are_separate() {
        let molecule = Molecule {
            atoms: vec![
                Atom {
                    atom_map: Some(7),
                    attachment_point: Some(3),
                    ..atom(0, 0)
                },
                atom(6, 2),
                Atom {
                    label: Some("R;🧪".into()),
                    ..atom(0, 0)
                },
            ],
            bonds: vec![bond(0, 1, BondOrder::Single), bond(1, 2, BondOrder::Single)],
            ..Molecule::default()
        };
        assert_eq!(
            molecule.to_cxsmiles().unwrap(),
            "[*:7]C* |$_AP3;;R&#59;&#55358;&#56810;$|"
        );
    }

    #[test]
    fn standard_radical_and_stated_hydrogen_layers_remain() {
        let molecule = Molecule {
            atoms: vec![
                Atom {
                    hydrogens_stated: true,
                    ..atom(6, 4)
                },
                Atom {
                    radical_electrons: 1,
                    ..atom(6, 3)
                },
                Atom {
                    radical_electrons: 2,
                    ..atom(6, 2)
                },
                Atom {
                    radical_electrons: 3,
                    ..atom(6, 1)
                },
                Atom {
                    isotope: Some(2),
                    ..atom(1, 0)
                },
            ],
            ..Molecule::default()
        };
        assert_eq!(
            molecule.to_cxsmiles().unwrap(),
            "[CH4].[CH3].[CH2].[CH].[2H] |^1:1,^2:2,^5:3|"
        );
    }

    #[test]
    fn tetrahedral_winding_and_unknown_are_distinct() {
        let mut molecule = Molecule {
            atoms: vec![atom(6, 1), atom(9, 0), atom(17, 0), atom(35, 0)],
            bonds: vec![
                bond(0, 1, BondOrder::Single),
                bond(0, 2, BondOrder::Single),
                bond(0, 3, BondOrder::Single),
            ],
            stereo: vec![Stereo::Tetrahedral {
                center: 0,
                neighbors: [
                    Ref::ImplicitHydrogen,
                    Ref::Atom(1),
                    Ref::Atom(2),
                    Ref::Atom(3),
                ],
                winding: Winding::Clockwise,
            }],
            ..Molecule::default()
        };
        assert_eq!(molecule.to_cxsmiles().unwrap(), "[C@@H](F)(Cl)Br");
        let Stereo::Tetrahedral { winding, .. } = &mut molecule.stereo[0] else {
            unreachable!()
        };
        *winding = Winding::Anticlockwise;
        assert_eq!(molecule.to_cxsmiles().unwrap(), "[C@H](F)(Cl)Br");
        molecule.stereo_groups.push(StereoGroup {
            kind: StereoGroupKind::And(1),
            stereo: vec![0],
        });
        assert_eq!(
            molecule.semantic_cxsmiles_for_tests().unwrap(),
            "[C@H](F)(Cl)Br |r|"
        );
        let Stereo::Tetrahedral { winding, .. } = &mut molecule.stereo[0] else {
            unreachable!()
        };
        *winding = Winding::Unknown;
        assert!(matches!(
            molecule.to_cxsmiles(),
            Err(CxsmilesError::StereoNotExpressible)
        ));
        molecule.stereo_groups.clear();
        assert_eq!(molecule.to_cxsmiles().unwrap(), "C(F)(Cl)Br |w:0.0|");
        assert_eq!(
            molecule.semantic_cxsmiles_for_tests().unwrap(),
            "C(F)(Cl)Br"
        );
    }

    #[test]
    fn compatibility_hydrogen_folding_retains_original_identity_mapping() {
        let molecule = Molecule {
            atoms: vec![atom(9, 0), atom(6, 0), atom(1, 0), atom(17, 0), atom(35, 0)],
            bonds: vec![
                bond(1, 0, BondOrder::Single),
                bond(1, 3, BondOrder::Single),
                bond(1, 4, BondOrder::Single),
                bond(1, 2, BondOrder::Single),
            ],
            stereo: vec![Stereo::Tetrahedral {
                center: 1,
                neighbors: [Ref::Atom(0), Ref::Atom(2), Ref::Atom(3), Ref::Atom(4)],
                winding: Winding::Clockwise,
            }],
            ..Molecule::default()
        };
        let original = molecule.clone();
        assert_eq!(molecule.to_cxsmiles().unwrap(), "F[C@@](Cl)(Br)[H]");
        let (semantic, order) = molecule
            .semantic_cxsmiles_with_atom_order_for_tests()
            .unwrap();
        assert_eq!(semantic, "F[C@@H](Cl)Br");
        assert_eq!(order, [0, 1, 3, 4]);
        assert_eq!(molecule.smiles_atom_order().unwrap(), [0, 1, 3, 4, 2]);
        assert_eq!(molecule, original);

        for annotation in 0..8 {
            let mut annotated = molecule.clone();
            match annotation {
                0 => annotated.parent.push(2),
                1 => annotated.atoms[2].locants.push("H".into()),
                2 => annotated.atoms[2].isotope = Some(2),
                3 => annotated.atoms[2].charge = 1,
                4 => annotated.atoms[2].radical_electrons = 1,
                5 => annotated.atoms[2].label = Some("hydrogen".into()),
                6 => annotated.atoms[2].atom_map = Some(7),
                7 => annotated.atoms[2].attachment_point = Some(1),
                _ => unreachable!(),
            }
            let (_, order) = annotated
                .semantic_cxsmiles_with_atom_order_for_tests()
                .unwrap();
            assert!(order.contains(&2), "annotation {annotation}");
        }
    }

    #[test]
    fn compatibility_folding_never_duplicates_a_virtual_hydrogen() {
        let molecule = Molecule {
            atoms: vec![atom(6, 1), atom(9, 0), atom(17, 0), atom(1, 0)],
            bonds: vec![
                bond(0, 1, BondOrder::Single),
                bond(0, 2, BondOrder::Single),
                bond(0, 3, BondOrder::Single),
            ],
            stereo: vec![Stereo::Tetrahedral {
                center: 0,
                neighbors: [
                    Ref::ImplicitHydrogen,
                    Ref::Atom(1),
                    Ref::Atom(2),
                    Ref::Atom(3),
                ],
                winding: Winding::Clockwise,
            }],
            ..Molecule::default()
        };
        let (projection, identities) = compatibility_projection(&molecule).unwrap();
        projection.validate().unwrap();
        assert_eq!(identities, [0, 1, 2, 3]);
        assert_eq!(projection, molecule);
    }

    #[test]
    fn unknown_tetrahedral_without_a_wiggly_single_bond_is_an_error() {
        let molecule = Molecule {
            atoms: vec![
                atom(6, 0),
                atom(9, 0),
                atom(17, 0),
                atom(35, 0),
                atom(53, 0),
            ],
            bonds: (1..5)
                .map(|atom| bond(0, atom, BondOrder::Double))
                .collect(),
            stereo: vec![Stereo::Tetrahedral {
                center: 0,
                neighbors: [Ref::Atom(1), Ref::Atom(2), Ref::Atom(3), Ref::Atom(4)],
                winding: Winding::Unknown,
            }],
            ..Molecule::default()
        };
        molecule.validate().unwrap();
        assert!(matches!(
            molecule.to_cxsmiles(),
            Err(CxsmilesError::StereoNotExpressible)
        ));
    }

    fn adjacent_unknown_centers() -> Molecule {
        Molecule {
            atoms: vec![
                atom(6, 1),
                atom(6, 1),
                atom(9, 0),
                atom(17, 0),
                atom(35, 0),
                atom(53, 0),
            ],
            bonds: vec![
                bond(0, 1, BondOrder::Single),
                bond(0, 2, BondOrder::Single),
                bond(0, 3, BondOrder::Single),
                bond(1, 4, BondOrder::Single),
                bond(1, 5, BondOrder::Single),
            ],
            stereo: vec![
                Stereo::Tetrahedral {
                    center: 0,
                    neighbors: [
                        Ref::ImplicitHydrogen,
                        Ref::Atom(1),
                        Ref::Atom(2),
                        Ref::Atom(3),
                    ],
                    winding: Winding::Unknown,
                },
                Stereo::Tetrahedral {
                    center: 1,
                    neighbors: [
                        Ref::ImplicitHydrogen,
                        Ref::Atom(0),
                        Ref::Atom(4),
                        Ref::Atom(5),
                    ],
                    winding: Winding::Unknown,
                },
            ],
            ..Molecule::default()
        }
    }

    #[test]
    fn adjacent_unknown_centers_get_distinct_wiggly_carriers() {
        let molecule = adjacent_unknown_centers();
        assert_eq!(
            molecule.to_cxsmiles().unwrap(),
            "C(C(Br)I)(F)Cl |w:0.0,1.1|"
        );
    }

    #[test]
    fn unknown_carrier_matching_can_reassign_an_earlier_center() {
        let mut molecule = adjacent_unknown_centers();
        // Centre 1 can use only the shared bond. Centre 0 initially chooses
        // it too, but can move to its unshared single bond instead.
        molecule.bonds[3].order = BondOrder::Double;
        molecule.bonds[4].order = BondOrder::Double;
        let mut writer = Writer::new(&molecule).unwrap();
        writer.assign_directions().unwrap();
        writer.assign_unknown_tetra_bonds().unwrap();
        assert_eq!(writer.unknown_tetra_bonds[0], Some(1));
        assert_eq!(writer.unknown_tetra_bonds[1], Some(0));
        molecule.to_cxsmiles().unwrap();
    }

    #[test]
    fn hall_deficient_unknown_carriers_return_a_typed_error() {
        let mut molecule = adjacent_unknown_centers();
        for bond in &mut molecule.bonds[1..] {
            bond.order = BondOrder::Double;
        }
        molecule.validate().unwrap();
        assert!(matches!(
            molecule.to_cxsmiles(),
            Err(CxsmilesError::StereoNotExpressible)
        ));
    }

    #[test]
    fn groups_never_drop_unknown_members_from_the_standard_projection() {
        let mut molecule = adjacent_unknown_centers();
        molecule.stereo_groups.push(StereoGroup {
            kind: StereoGroupKind::And(1),
            stereo: vec![0],
        });
        assert!(matches!(
            molecule.to_cxsmiles(),
            Err(CxsmilesError::StereoNotExpressible)
        ));
        let Stereo::Tetrahedral { winding, .. } = &mut molecule.stereo[0] else {
            unreachable!()
        };
        *winding = Winding::Clockwise;
        molecule.stereo_groups[0].stereo = vec![0, 1];
        assert!(matches!(
            molecule.to_cxsmiles(),
            Err(CxsmilesError::StereoNotExpressible)
        ));
        molecule.stereo_groups[0].stereo = vec![0];
        molecule.to_cxsmiles().unwrap();

        let mut double = polyene(&[3], false);
        double.stereo_groups.push(StereoGroup {
            kind: StereoGroupKind::Or(2),
            stereo: vec![0],
        });
        assert!(matches!(
            double.to_cxsmiles(),
            Err(CxsmilesError::StereoNotExpressible)
        ));
    }

    #[test]
    fn every_axis_state_is_a_typed_error_without_a_standard_encoding() {
        let mut molecule = Molecule {
            atoms: vec![
                atom(6, 0),
                atom(6, 0),
                atom(9, 0),
                atom(17, 0),
                atom(35, 0),
                atom(53, 0),
            ],
            bonds: vec![
                bond(0, 1, BondOrder::Single),
                bond(0, 2, BondOrder::Single),
                bond(0, 3, BondOrder::Single),
                bond(1, 4, BondOrder::Single),
                bond(1, 5, BondOrder::Single),
            ],
            stereo: vec![Stereo::Atropisomer {
                bond: 0,
                refs: [Ref::Atom(2), Ref::Atom(3), Ref::Atom(4), Ref::Atom(5)],
                winding: AxialWinding::M,
            }],
            ..Molecule::default()
        };
        for state in [AxialWinding::M, AxialWinding::P, AxialWinding::Unknown] {
            let Stereo::Atropisomer { winding, .. } = &mut molecule.stereo[0] else {
                unreachable!()
            };
            *winding = state;
            assert!(matches!(
                molecule.to_cxsmiles(),
                Err(CxsmilesError::StereoNotExpressible)
            ));
        }
    }

    #[test]
    fn partially_specified_triene_does_not_assert_the_middle_bond() {
        // (2E,6E)-octa-2,4,6-triene. The required slash on each side of
        // the unspecified middle double bond cannot be removed.
        for state in [0, 3] {
            let molecule = polyene(&[2, state, 2], false);
            assert!(matches!(
                molecule.to_cxsmiles(),
                Err(CxsmilesError::StereoNotExpressible)
            ));
            assert!(matches!(
                molecule.semantic_cxsmiles_for_tests(),
                Err(CxsmilesError::StereoNotExpressible)
            ));
        }
    }

    #[test]
    fn exhaustive_linear_polyene_slash_patterns_are_exact_or_fail_closed() {
        // All 4^5 combinations of absent/cis/trans/unknown, with both
        // endpoint orientations. A chain has one possible mark at either
        // side of each double bond; only an unknown/absent bond flanked by
        // known neighbours makes its standard projection impossible.
        for reversed in [false, true] {
            for code in 0_u32..4_u32.pow(5) {
                let states = (0..5)
                    .map(|i| ((code >> (2 * i)) & 3) as u8)
                    .collect::<Vec<_>>();
                let known = states
                    .iter()
                    .map(|&state| state == 1 || state == 2)
                    .collect::<Vec<_>>();
                let impossible = (1..4).any(|i| !known[i] && known[i - 1] && known[i + 1]);
                let molecule = polyene(&states, reversed);
                let mut writer = Writer::new(&molecule).unwrap();
                match writer.assign_directions() {
                    Err(CxsmilesError::StereoNotExpressible) => assert!(impossible, "{states:?}"),
                    Err(error) => panic!("{states:?}: {error}"),
                    Ok(()) => {
                        assert!(!impossible, "{states:?}");
                        let output = writer.write(false).unwrap();
                        let smiles = output
                            .split_once(" |")
                            .map_or(output.as_str(), |(smiles, _)| smiles);
                        // This independent chain tokenizer checks the standard
                        // output rather than the solver's assignment arrays.
                        let mut separators = Vec::new();
                        let mut separator = None;
                        let mut first = true;
                        for ch in smiles.chars() {
                            match ch {
                                '/' | '\\' | '=' => separator = Some(ch),
                                'C' => {
                                    if !first {
                                        separators.push(separator.take());
                                    }
                                    first = false;
                                }
                                _ => {}
                            }
                        }
                        for i in 0..5 {
                            let left = separators[2 * i];
                            let right = separators[2 * i + 2];
                            let asserted = matches!(left, Some('/' | '\\'))
                                && matches!(right, Some('/' | '\\'));
                            assert_eq!(asserted, known[i], "{states:?}: {output}");
                            if asserted {
                                assert_eq!(left == right, states[i] == 2, "{states:?}: {output}");
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn mark_solver_agrees_with_exhaustive_boolean_assignments() {
        let literals = [(0, false), (0, true), (1, false), (1, true)];
        let possibilities = (0..literals.len())
            .flat_map(|i| (i..literals.len()).map(move |j| (literals[i], literals[j])))
            .collect::<Vec<_>>();
        for mask in 0..1 << possibilities.len() {
            let clauses = possibilities
                .iter()
                .enumerate()
                .filter_map(|(i, &clause)| (mask & (1 << i) != 0).then_some(clause))
                .collect::<Vec<_>>();
            let satisfies = |assignment: &[bool]| {
                clauses
                    .iter()
                    .all(|&((a, av), (b, bv))| assignment[a] == av || assignment[b] == bv)
            };
            let expected =
                (0..4).any(|assignment| satisfies(&[assignment & 1 != 0, assignment & 2 != 0]));
            let actual = solve_marks(2, &clauses);
            assert_eq!(actual.is_some(), expected, "{clauses:?}");
            if let Some(assignment) = actual {
                assert!(satisfies(&assignment), "{clauses:?}: {assignment:?}");
            }
        }
    }

    #[test]
    fn polymer_endpoints_and_position_variation_use_standard_layers() {
        let molecule = Molecule {
            atoms: vec![atom(0, 0), atom(8, 0), atom(6, 2), atom(6, 2), atom(0, 0)],
            bonds: vec![
                bond(0, 1, BondOrder::Single),
                bond(1, 2, BondOrder::Single),
                bond(2, 3, BondOrder::Single),
                bond(3, 4, BondOrder::Single),
            ],
            polymers: vec![RepeatUnit {
                atoms: vec![1, 2, 3],
                crossing_bonds: vec![0, 3],
                attachment_points: vec![0, 4],
                connectivity: Connectivity::HeadToTail,
                label: None,
            }],
            position_variation: vec![PositionVariation {
                bond: 0,
                endpoints: vec![1, 2],
            }],
            ..Molecule::default()
        };
        assert_eq!(
            molecule.to_cxsmiles().unwrap(),
            "*OCC* |$star_e;;;;star_e$,m:0:1.2,Sg:n:1,2,3::ht|"
        );
        assert_eq!(
            molecule.semantic_cxsmiles_for_tests().unwrap(),
            "[*:1]OCC[*:2] |$star_e;;;;star_e$,m:0:1.2,Sg:n:1,2,3::ht|"
        );
    }

    #[test]
    fn dense_graphs_have_no_fixed_ring_closure_limit() {
        let mut molecule = Molecule {
            atoms: vec![atom(6, 0); 24],
            ..Molecule::default()
        };
        for first in 0..24 {
            for second in first + 1..24 {
                molecule.bonds.push(bond(first, second, BondOrder::Single));
            }
        }
        assert!(molecule.to_cxsmiles().unwrap().contains("%(100)"));
    }
}
