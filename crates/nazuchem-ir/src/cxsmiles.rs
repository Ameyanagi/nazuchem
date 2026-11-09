// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | crates/nazuchem-ir/src/cxsmiles.rs | Phase B step 8a
// Copyright (c) 2026 Ameyanagi and Nazuchem contributors
// Changes: deterministic IR v1 SMILES traversal and lossless CXSMILES projection.

//! Deterministic, perception-free CXSMILES export.
//!
//! See `docs/cxsmiles.md` for the standard extension mapping and the versioned
//! `nazuchem*` atom properties that preserve information beyond standard CXSMILES.
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
pub enum CxSmilesError {
    InvalidStructure(ValidationError),
    Metadata(serde_json::Error),
    ContradictoryStereo,
}

impl fmt::Display for CxSmilesError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidStructure(error) => write!(f, "invalid structure: {error}"),
            Self::Metadata(error) => write!(f, "cannot encode CXSMILES metadata: {error}"),
            Self::ContradictoryStereo => f.write_str("contradictory double bond stereochemistry"),
        }
    }
}

impl std::error::Error for CxSmilesError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::InvalidStructure(error) => Some(error),
            Self::Metadata(error) => Some(error),
            Self::ContradictoryStereo => None,
        }
    }
}

impl From<ValidationError> for CxSmilesError {
    fn from(error: ValidationError) -> Self {
        Self::InvalidStructure(error)
    }
}

impl From<serde_json::Error> for CxSmilesError {
    fn from(error: serde_json::Error) -> Self {
        Self::Metadata(error)
    }
}

impl Molecule {
    /// Export a lossless CXSMILES view of this IR, retaining atom identities,
    /// complete locants, hydrogen origin, and exact ordered graph annotations.
    ///
    /// Losslessness requires preserving the documented `nazuchem*` properties;
    /// consumers that discard arbitrary atom properties keep the standard
    /// chemical projection but cannot reconstruct every IR field.
    pub fn to_cxsmiles(&self) -> Result<String, CxSmilesError> {
        Writer::new(self)?.write(false)
    }

    /// Export the legacy OPSIN semantic projection, without locants, hydrogen
    /// origin, explicit radical layers, or IR-specific metadata.
    ///
    /// This is deterministic, not canonical. Atom order, ring numbering and
    /// globally equivalent slash directions can differ from other writers.
    pub fn semantic_cxsmiles(&self) -> Result<String, CxSmilesError> {
        Writer::new(self)?.write(true)
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

    /// IR atom identities in the order used by both CXSMILES projections.
    pub fn smiles_atom_order(&self) -> Result<Vec<u32>, CxSmilesError> {
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
}

impl<'a> Writer<'a> {
    fn new(molecule: &'a Molecule) -> Result<Self, CxSmilesError> {
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

    fn assign_directions(&mut self, semantic: bool) -> Result<(), CxSmilesError> {
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
                    known: true,
                } => Some((*bond as usize, *refs, *relation)),
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
                // Exact references remain in metadata; cumulene-like records
                // have no direct slash encoding in the v1 stereo vocabulary.
                continue;
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
                            if semantic {
                                return Err(CxSmilesError::ContradictoryStereo);
                            }
                            // The exact incompatible records are still lossless
                            // in metadata; avoid falsely asserting a stereo view.
                            self.directions.fill(None);
                            return Ok(());
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
        Ok(())
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

    fn write(mut self, semantic: bool) -> Result<String, CxSmilesError> {
        self.assign_directions(semantic)?;
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

    fn tetrahedral(&self, atom: u32) -> Option<(&[Ref; 4], Winding, bool)> {
        self.molecule.stereo.iter().find_map(|stereo| match stereo {
            Stereo::Tetrahedral {
                center,
                neighbors,
                winding,
                known,
            } if *center == atom => Some((neighbors, *winding, *known)),
            _ => None,
        })
    }

    fn parity(&self, atom: u32, semantic: bool) -> &'static str {
        let Some((stored, winding, known)) = self.tetrahedral(atom) else {
            return "";
        };
        if !semantic && !known {
            return "";
        }
        let mut current = Vec::with_capacity(4);
        if let Some(parent) = self.parent_bond[atom as usize] {
            current.push(Ref::Atom(self.other(parent, atom)));
        }
        // A virtual ligand is first at a root, or follows the incoming atom.
        // SMILES has one implicit ligand slot; exact exotic records containing
        // two virtual ligands are carried by the metadata instead.
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
        let parity = self.parity(atom, semantic);
        let polymer_class = if semantic {
            self.polymer_class(atom)
        } else {
            None
        };
        let brackets = polymer_class.is_some()
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
        if let Some(class) = polymer_class {
            text.push(':');
            text.push_str(&class.to_string());
        }
        if brackets {
            text.push(']');
        }
        text
    }

    fn polymer_class(&self, atom: u32) -> Option<usize> {
        if self.molecule.atoms[atom as usize].label.as_deref() != Some("star_e") {
            return None;
        }
        let points = self
            .molecule
            .polymers
            .iter()
            .flat_map(|unit| &unit.crossing_bonds)
            .filter_map(|&bond| {
                self.molecule.bonds[bond as usize]
                    .atoms
                    .into_iter()
                    .find(|&id| self.molecule.atoms[id as usize].label.as_deref() == Some("star_e"))
            })
            .collect::<Vec<_>>();
        let mut unique = Vec::new();
        for point in points {
            if !unique.contains(&point) {
                unique.push(point);
            }
        }
        unique
            .iter()
            .position(|&point| point == atom)
            .map(|index| index + 1)
    }

    fn layers(&self, semantic: bool) -> Result<Vec<String>, CxSmilesError> {
        let mut layers = Vec::new();
        let labels = self
            .order
            .iter()
            .map(|&id| {
                self.molecule.atoms[id as usize]
                    .label
                    .as_deref()
                    .map(escape_label)
                    .unwrap_or_default()
            })
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
                        Stereo::Tetrahedral { center, .. } => Some(self.index[*center as usize]),
                        Stereo::DoubleBond { .. } => None,
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
                StereoGroupKind::Abs => {}
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
                    StereoGroupKind::Abs if semantic => continue,
                    StereoGroupKind::Abs => "a".into(),
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
                        known: false,
                        ..
                    } => {
                        if let Some(&bond) = self.adjacent[*center as usize]
                            .iter()
                            .find(|&&bond| self.molecule.bonds[bond].order == BondOrder::Single)
                        {
                            wiggly.push(format!(
                                "{}.{}",
                                self.index[*center as usize], self.bond_index[bond]
                            ));
                        }
                    }
                    Stereo::DoubleBond {
                        bond, known: false, ..
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
            self.metadata(&mut layers)?;
        }
        Ok(layers)
    }

    fn metadata(&self, layers: &mut Vec<String>) -> Result<(), CxSmilesError> {
        if self.order.is_empty() {
            return Ok(());
        }
        let mut properties = vec!["0.nazuchemVersion.1".to_owned()];
        let reordered = self
            .order
            .iter()
            .enumerate()
            .any(|(index, &atom)| index != atom as usize);
        for (index, &id) in self.order.iter().enumerate() {
            let atom = &self.molecule.atoms[id as usize];
            if reordered {
                properties.push(format!("{index}.nazuchemIndex.{id}"));
            }
            if self.molecule.parent.contains(&id) {
                properties.push(format!("{index}.nazuchemParent.1"));
            }
            if !atom.locants.is_empty() {
                properties.push(format!(
                    "{index}.nazuchemLocants.{}",
                    escape(&serde_json::to_string(&atom.locants)?)
                ));
            }
            // Bracketed atoms otherwise look as though their H count was stated.
            if !atom.hydrogens_stated && self.atom_text(id, false).starts_with('[') {
                properties.push(format!("{index}.nazuchemHydrogensStated.0"));
            }
            // IR validation intentionally permits chemically unusual counts.
            // Never let a consumer infer these values from valence or radicals
            // suggested by bracket notation.
            properties.push(format!("{index}.nazuchemHydrogens.{}", atom.hydrogens));
            properties.push(format!(
                "{index}.nazuchemRadicalElectrons.{}",
                atom.radical_electrons
            ));
            if let Some(label) = &atom.label {
                properties.push(format!(
                    "{index}.nazuchemLabel.{}",
                    escape(&serde_json::to_string(label)?)
                ));
            }
        }
        // Standard CX does not retain array order, bond endpoint orientation,
        // stereo reference order, group references to bonds, or all SRU bonds.
        let mut topology = serde_json::Map::new();
        for (name, value, empty) in [
            (
                "bonds",
                serde_json::to_value(&self.molecule.bonds)?,
                self.molecule.bonds.is_empty(),
            ),
            (
                "stereo",
                serde_json::to_value(&self.molecule.stereo)?,
                self.molecule.stereo.is_empty(),
            ),
            (
                "stereo_groups",
                serde_json::to_value(&self.molecule.stereo_groups)?,
                self.molecule.stereo_groups.is_empty(),
            ),
            (
                "position_variation",
                serde_json::to_value(&self.molecule.position_variation)?,
                self.molecule.position_variation.is_empty(),
            ),
            (
                "polymers",
                serde_json::to_value(&self.molecule.polymers)?,
                self.molecule.polymers.is_empty(),
            ),
            (
                "parent",
                serde_json::to_value(&self.molecule.parent)?,
                self.molecule.parent.is_empty(),
            ),
        ] {
            if !empty {
                topology.insert(name.into(), value);
            }
        }
        if !topology.is_empty() {
            properties.push(format!(
                "0.nazuchemTopology.{}",
                escape(&serde_json::to_string(&topology)?)
            ));
        }
        layers.push(format!("atomProp:{}", properties.join(":")));
        Ok(())
    }
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
    use crate::{Atom, Bond, PositionVariation, RepeatUnit, StereoGroup};

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

    // A deliberately small test oracle, not a production CXSMILES importer.
    // It reads atom tokens and the documented exact metadata, without valence
    // perception or access to the source molecule.
    fn restore_export(text: &str) -> Molecule {
        if text.is_empty() {
            return Molecule::default();
        }
        let (smiles, cx) = text.split_once(" |").unwrap();
        let mut atoms = Vec::new();
        let bytes = smiles.as_bytes();
        let mut cursor = 0;
        while cursor < bytes.len() {
            let bracketed = bytes[cursor] == b'[';
            if !bracketed && bytes[cursor] != b'*' && !bytes[cursor].is_ascii_uppercase() {
                cursor += 1;
                continue;
            }
            let start = cursor + usize::from(bracketed);
            let end = if bracketed {
                start + smiles[start..].find(']').unwrap()
            } else {
                cursor + 1 + usize::from(bytes.get(cursor + 1).is_some_and(u8::is_ascii_lowercase))
            };
            let token = &smiles[start..end];
            let digit_count = token.bytes().take_while(u8::is_ascii_digit).count();
            let isotope = (digit_count > 0).then(|| token[..digit_count].parse().unwrap());
            let symbol_len = 1 + usize::from(
                token
                    .as_bytes()
                    .get(digit_count + 1)
                    .is_some_and(u8::is_ascii_lowercase),
            );
            let symbol = &token[digit_count..digit_count + symbol_len];
            let element = ELEMENTS.iter().position(|&value| value == symbol).unwrap() as u8;
            let charge = token
                .find(['+', '-'])
                .map(|position| {
                    let sign = if token.as_bytes()[position] == b'+' {
                        1_i16
                    } else {
                        -1_i16
                    };
                    let tail = &token[position + 1..];
                    let digits = tail.bytes().take_while(u8::is_ascii_digit).count();
                    (sign
                        * if digits == 0 {
                            1
                        } else {
                            tail[..digits].parse::<i16>().unwrap()
                        }) as i8
                })
                .unwrap_or_default();
            atoms.push(Atom {
                element,
                isotope,
                charge,
                hydrogens_stated: bracketed,
                ..Atom::default()
            });
            cursor = end + usize::from(bracketed);
        }
        let block = cx.rsplit_once("atomProp:").unwrap().1.trim_end_matches('|');
        let mut identities = (0..atoms.len()).collect::<Vec<_>>();
        let mut result = Molecule::default();
        for property in block.split(':') {
            let mut parts = property.splitn(3, '.');
            let index = parts.next().unwrap().parse::<usize>().unwrap();
            let key = parts.next().unwrap();
            let value = decode(parts.next().unwrap());
            match key {
                "nazuchemVersion" => assert_eq!(value, "1"),
                "nazuchemIndex" => identities[index] = value.parse().unwrap(),
                "nazuchemParent" => assert_eq!(value, "1"),
                "nazuchemHydrogens" => atoms[index].hydrogens = value.parse().unwrap(),
                "nazuchemHydrogensStated" => atoms[index].hydrogens_stated = value == "1",
                "nazuchemRadicalElectrons" => {
                    atoms[index].radical_electrons = value.parse().unwrap()
                }
                "nazuchemLabel" => atoms[index].label = Some(serde_json::from_str(&value).unwrap()),
                "nazuchemLocants" => atoms[index].locants = serde_json::from_str(&value).unwrap(),
                "nazuchemTopology" => {
                    let mut topology: serde_json::Map<String, serde_json::Value> =
                        serde_json::from_str(&value).unwrap();
                    macro_rules! restore {
                        ($field:ident) => {
                            if let Some(value) = topology.remove(stringify!($field)) {
                                result.$field = serde_json::from_value(value).unwrap();
                            }
                        };
                    }
                    restore!(bonds);
                    restore!(stereo);
                    restore!(stereo_groups);
                    restore!(position_variation);
                    restore!(polymers);
                    restore!(parent);
                    assert!(topology.is_empty());
                }
                other => panic!("unknown metadata property {other}"),
            }
        }
        result.atoms.resize(atoms.len(), Atom::default());
        for (identity, atom) in identities.into_iter().zip(atoms) {
            result.atoms[identity] = atom;
        }
        result.validate().unwrap();
        result
    }

    fn decode(value: &str) -> String {
        let mut units = Vec::new();
        let mut rest = value;
        while !rest.is_empty() {
            if let Some(entity) = rest.strip_prefix("&#") {
                let (number, tail) = entity.split_once(';').unwrap();
                units.push(number.parse().unwrap());
                rest = tail;
            } else {
                let character = rest.chars().next().unwrap();
                let mut encoded = [0; 2];
                units.extend_from_slice(character.encode_utf16(&mut encoded));
                rest = &rest[character.len_utf8()..];
            }
        }
        String::from_utf16(&units).unwrap()
    }

    #[test]
    fn full_export_round_trips_every_ir_field_without_chemical_inference() {
        let mut molecule = Molecule {
            atoms: vec![
                atom(6, 1),
                atom(9, 0),
                atom(17, 0),
                atom(35, 0),
                atom(6, 1),
                atom(6, 1),
                atom(9, 0),
                atom(17, 0),
                atom(16, 0),
                atom(8, 0),
                atom(7, 0),
                atom(6, 3),
                atom(0, 0),
                atom(6, 2),
                atom(0, 0),
                atom(1, 0),
                atom(118, 3),
            ],
            bonds: vec![
                bond(0, 1, BondOrder::Single),
                bond(0, 2, BondOrder::Single),
                bond(0, 3, BondOrder::Single),
                bond(5, 4, BondOrder::Double),
                bond(4, 6, BondOrder::Single),
                bond(7, 5, BondOrder::Single),
                bond(8, 9, BondOrder::Single),
                bond(10, 8, BondOrder::Single),
                bond(8, 11, BondOrder::Single),
                bond(12, 13, BondOrder::Single),
                bond(13, 14, BondOrder::Single),
                bond(15, 13, BondOrder::Single),
            ],
            stereo: vec![
                Stereo::Tetrahedral {
                    center: 0,
                    neighbors: [
                        Ref::Atom(2),
                        Ref::ImplicitHydrogen,
                        Ref::Atom(3),
                        Ref::Atom(1),
                    ],
                    winding: Winding::Anticlockwise,
                    known: false,
                },
                Stereo::DoubleBond {
                    bond: 3,
                    refs: [7, 6],
                    relation: Relation::Cis,
                    known: true,
                },
                Stereo::Tetrahedral {
                    center: 8,
                    neighbors: [Ref::LonePair, Ref::Atom(11), Ref::Atom(9), Ref::Atom(10)],
                    winding: Winding::Clockwise,
                    known: true,
                },
            ],
            stereo_groups: vec![
                StereoGroup {
                    kind: StereoGroupKind::Or(99),
                    stereo: vec![1],
                },
                StereoGroup {
                    kind: StereoGroupKind::And(17),
                    stereo: vec![0],
                },
                StereoGroup {
                    kind: StereoGroupKind::Abs,
                    stereo: vec![2],
                },
            ],
            position_variation: vec![PositionVariation {
                bond: 9,
                endpoints: vec![13, 0, 5],
            }],
            polymers: vec![
                RepeatUnit {
                    atoms: vec![13],
                    crossing_bonds: vec![11, 10, 9],
                    connectivity: Connectivity::HeadToHead,
                    label: Some(String::new()),
                },
                RepeatUnit {
                    atoms: vec![5, 4],
                    crossing_bonds: vec![5, 4],
                    connectivity: Connectivity::Unknown,
                    label: None,
                },
            ],
            parent: vec![13, 0, 5],
        };
        molecule.atoms[0].hydrogens_stated = true;
        molecule.atoms[0].locants = vec!["1".into(), "1′".into(), "".into()];
        molecule.atoms[5].locants = vec!["".into()];
        molecule.atoms[6].locants = vec!["outside parent".into()];
        molecule.atoms[12].label = Some("_x;|$:,&🧪".into());
        molecule.atoms[13].label = Some("".into());
        molecule.atoms[14].label = Some("star_e".into());
        molecule.atoms[15].isotope = Some(2);
        molecule.atoms[16].isotope = Some(u16::MAX);
        molecule.atoms[16].charge = i8::MIN;
        molecule.atoms[16].radical_electrons = u8::MAX;
        let output = molecule.to_cxsmiles().unwrap();
        assert!(!output.contains("_x;|$"));
        assert!(output.contains("$_AV:"));
        assert!(output.contains("m:"));
        assert!(output.contains("Sg:n:"));
        assert!(output.contains("w:"));
        assert_eq!(restore_export(&output), molecule);
        assert_eq!(molecule.smiles_atom_order().unwrap()[0], 12);
        assert_eq!(molecule.displayed_locants()[0], Some("1".into()));
        assert_eq!(molecule.displayed_locants()[6], None);
        assert_eq!(output, molecule.to_cxsmiles().unwrap());
    }

    #[test]
    fn simple_hydrogen_origin_and_radical_counts_are_exact() {
        for stated in [false, true] {
            for element in [0, 1, 6, 17, 118] {
                for hydrogens in [0, 1, 4, 255] {
                    for radical in [0, 1, 2, 3, 4, 255] {
                        let molecule = Molecule {
                            atoms: vec![Atom {
                                hydrogens_stated: stated,
                                radical_electrons: radical,
                                ..atom(element, hydrogens)
                            }],
                            ..Molecule::default()
                        };
                        assert_eq!(restore_export(&molecule.to_cxsmiles().unwrap()), molecule);
                    }
                }
            }
        }
    }

    #[test]
    fn semantic_projection_keeps_legacy_stereo_and_polymer_conventions() {
        let molecule = Molecule {
            atoms: vec![
                Atom {
                    label: Some("star_e".into()),
                    ..atom(0, 0)
                },
                atom(8, 0),
                atom(6, 2),
                atom(6, 2),
                Atom {
                    label: Some("star_e".into()),
                    ..atom(0, 0)
                },
            ],
            bonds: vec![
                bond(0, 1, BondOrder::Single),
                bond(1, 2, BondOrder::Single),
                bond(2, 3, BondOrder::Single),
                bond(3, 4, BondOrder::Single),
            ],
            polymers: vec![RepeatUnit {
                atoms: vec![1, 2, 3],
                crossing_bonds: vec![0, 3],
                connectivity: Connectivity::HeadToTail,
                label: None,
            }],
            ..Molecule::default()
        };
        assert_eq!(
            molecule.semantic_cxsmiles().unwrap(),
            "[*:1]OCC[*:2] |$star_e;;;;star_e$,Sg:n:1,2,3::ht|"
        );
        let mut stereo = Molecule {
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
                known: false,
            }],
            ..Molecule::default()
        };
        assert_eq!(stereo.semantic_cxsmiles().unwrap(), "[C@@H](F)(Cl)Br");
        assert!(
            !stereo
                .to_cxsmiles()
                .unwrap()
                .split_once(" |")
                .unwrap()
                .0
                .contains('@')
        );
        stereo.stereo_groups.push(StereoGroup {
            kind: StereoGroupKind::And(1),
            stereo: vec![0],
        });
        assert!(stereo.semantic_cxsmiles().unwrap().ends_with(" |r|"));
    }

    #[test]
    fn bonds_and_ring_closures_keep_exact_identities() {
        let molecule = Molecule {
            atoms: vec![atom(6, 0), atom(6, 0), atom(6, 0), atom(6, 0)],
            bonds: vec![
                bond(2, 0, BondOrder::Quadruple),
                bond(1, 2, BondOrder::Triple),
                bond(3, 1, BondOrder::Single),
                bond(0, 3, BondOrder::Double),
                bond(0, 1, BondOrder::Single),
            ],
            ..Molecule::default()
        };
        let output = molecule.to_cxsmiles().unwrap();
        assert!(output.contains('$'));
        assert!(output.contains('#'));
        assert_eq!(restore_export(&output), molecule);
        assert_eq!(
            restore_export(&Molecule::default().to_cxsmiles().unwrap()),
            Molecule::default()
        );
    }

    #[test]
    fn cis_trans_and_unknown_double_bond_projection_are_distinct() {
        let mut molecule = Molecule {
            atoms: vec![atom(6, 3), atom(6, 1), atom(6, 1), atom(6, 3)],
            bonds: vec![
                bond(0, 1, BondOrder::Single),
                bond(1, 2, BondOrder::Double),
                bond(2, 3, BondOrder::Single),
            ],
            stereo: vec![Stereo::DoubleBond {
                bond: 1,
                refs: [0, 3],
                relation: Relation::Trans,
                known: true,
            }],
            ..Molecule::default()
        };
        assert_eq!(molecule.semantic_cxsmiles().unwrap(), r"C\C=C\C");
        molecule.stereo[0] = Stereo::DoubleBond {
            bond: 1,
            refs: [0, 3],
            relation: Relation::Cis,
            known: true,
        };
        assert_eq!(molecule.semantic_cxsmiles().unwrap(), r"C\C=C/C");
        molecule.stereo[0] = Stereo::DoubleBond {
            bond: 1,
            refs: [0, 3],
            relation: Relation::Trans,
            known: false,
        };
        let output = molecule.to_cxsmiles().unwrap();
        assert!(output.starts_with("CC=CC |w:1.1,ctu:1,"));
        assert_eq!(restore_export(&output), molecule);
        assert_eq!(molecule.semantic_cxsmiles().unwrap(), "CC=CC");
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
        let output = molecule.to_cxsmiles().unwrap();
        assert!(output.contains("%(100)"));
        assert_eq!(restore_export(&output), molecule);
    }

    #[test]
    fn position_variation_source_is_the_first_bond_endpoint() {
        for (first_element, second_element, pair, expected) in [
            (6, 0, [0, 1], "*CC |m:1:1.2|"),
            (6, 0, [1, 0], "*CC |m:0:1.2|"),
            (0, 0, [1, 0], "**C |m:1:1.2|"),
            (0, 0, [0, 1], "**C |m:0:1.2|"),
        ] {
            let molecule = Molecule {
                atoms: vec![
                    atom(first_element, if first_element == 6 { 2 } else { 0 }),
                    atom(second_element, 0),
                    atom(6, 3),
                ],
                bonds: vec![
                    bond(pair[0], pair[1], BondOrder::Single),
                    bond(if first_element == 6 { 0 } else { 1 }, 2, BondOrder::Single),
                ],
                position_variation: vec![PositionVariation {
                    bond: 0,
                    endpoints: vec![if first_element == 6 { 0 } else { 1 }, 2],
                }],
                ..Molecule::default()
            };
            assert_eq!(molecule.semantic_cxsmiles().unwrap(), expected);
            let expected_layer = expected.split_once(" |").unwrap().1.trim_end_matches('|');
            let full = molecule.to_cxsmiles().unwrap();
            assert!(
                full.contains(&format!(" |{expected_layer},atomProp:")),
                "{full}"
            );
            assert_eq!(restore_export(&full), molecule);
        }
    }
}
