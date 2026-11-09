// SPDX-License-Identifier: MIT
// Origin: original
// Source: Nazuchem | crates/nazuchem-n2s/tests/support/ir_cxsmiles.rs | Phase B step 8a

//! Normalize only the standard semantic CX layers emitted by the legacy writer.
//! Every written atom index is mapped back to the adapter's stable IR identity.
//! Unknown layers fail closed: graph identity alone must never erase metadata.

use std::collections::BTreeMap;

use nazuchem_ir::{BondOrder, Molecule};
use nazuchem_kernel::smiles;

#[derive(Debug, Default, PartialEq, Eq)]
pub struct SemanticLayers {
    labels: BTreeMap<u32, String>,
    relative: bool,
    stereo: Vec<(char, Vec<u32>)>,
    positions: Vec<(u32, Vec<u32>)>,
    polymers: Vec<(Vec<u32>, String)>,
}

pub fn bare(text: &str) -> &str {
    text.split_once(" |").map_or(text, |(smiles, _)| smiles)
}

/// Legacy bracket valence deficits carry radical counts even when no OPSIN
/// out-atom records remain. Check that typed IR field independently of text.
pub fn legacy_radicals(molecule: &Molecule, legacy: &str) -> Result<(), String> {
    let written = smiles::parse(bare(legacy)).map_err(|error| error.to_string())?;
    let expected = written.graph.assign_radicals()?;
    let actual: Vec<_> = molecule
        .atoms
        .iter()
        .map(|atom| atom.radical_electrons)
        .collect();
    if expected != actual {
        return Err(format!(
            "IR radical counts differ from legacy SMILES: expected {expected:?}, actual {actual:?}"
        ));
    }
    Ok(())
}

/// Prove that the writer's index map describes the graph actually emitted.
/// Canonical identity alone could hide an incorrect map between equivalent
/// graphs and would then associate locants or CX layers with the wrong atoms.
pub fn written_graph(molecule: &Molecule, text: &str, order: &[u32]) -> Result<(), String> {
    let written = smiles::parse(bare(text)).map_err(|error| error.to_string())?;
    if written.graph.atoms.len() != order.len() {
        return Err("writer index map has the wrong atom count".into());
    }
    for (output, &identity) in written.graph.atoms.iter().zip(order) {
        let expected = &molecule.atoms[identity as usize];
        if (output.atomic_number, output.isotope, output.charge)
            != (
                expected.element,
                expected.isotope.unwrap_or(0),
                expected.charge,
            )
        {
            return Err(format!(
                "writer atom index map disagrees at identity {identity}"
            ));
        }
    }
    let edge = |a: u32, b: u32, order: u8| (a.min(b), a.max(b), order);
    let mut actual: Vec<_> = written
        .graph
        .bonds
        .iter()
        .map(|bond| edge(order[bond.a], order[bond.b], bond.order))
        .collect();
    let mut expected = Vec::new();
    for bond in &molecule.bonds {
        let order = match bond.order {
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Quadruple => 6,
            _ => return Err("unknown IR bond order".into()),
        };
        expected.push(edge(bond.atoms[0], bond.atoms[1], order));
    }
    expected.sort_unstable();
    actual.sort_unstable();
    if expected != actual {
        return Err("writer atom index map disagrees with emitted bonds".into());
    }
    Ok(())
}

pub fn semantic_layers(text: &str, order: &[u32]) -> Result<SemanticLayers, String> {
    let Some((_, extension)) = text.split_once(" |") else {
        return Ok(SemanticLayers::default());
    };
    let extension = extension
        .strip_suffix('|')
        .ok_or_else(|| "unterminated semantic CX extension".to_owned())?;
    let mut segments = Vec::new();
    let (mut start, mut in_labels) = (0, false);
    let bytes = extension.as_bytes();
    for (index, &byte) in bytes.iter().enumerate() {
        if byte == b'$' {
            in_labels = !in_labels;
        } else if byte == b','
            && !in_labels
            && bytes
                .get(index + 1)
                .is_some_and(|next| next.is_ascii_alphabetic() || matches!(*next, b'$' | b'&'))
        {
            segments.push(&extension[start..index]);
            start = index + 1;
        }
    }
    segments.push(&extension[start..]);
    let atom = |text: &str| -> Result<u32, String> {
        let index = text
            .parse::<usize>()
            .map_err(|_| format!("invalid CX atom index {text}"))?;
        order
            .get(index)
            .copied()
            .ok_or_else(|| format!("CX index {index} outside written atoms"))
    };
    let atoms = |text: &str, separator| -> Result<Vec<u32>, String> {
        let mut atoms = text
            .split(separator)
            .map(&atom)
            .collect::<Result<Vec<_>, _>>()?;
        atoms.sort_unstable();
        Ok(atoms)
    };
    let mut normalized = SemanticLayers::default();
    for segment in segments {
        if let Some(labels) = segment.strip_prefix('$').and_then(|s| s.strip_suffix('$')) {
            for (index, label) in split_labels(labels)?.into_iter().enumerate() {
                if !label.is_empty() {
                    let identity = *order.get(index).ok_or("CX label outside written atoms")?;
                    if normalized
                        .labels
                        .insert(identity, label.to_owned())
                        .is_some()
                    {
                        return Err("duplicate semantic atom label".into());
                    }
                }
            }
        } else if segment == "r" {
            normalized.relative = true;
        } else if let Some(body) = segment.strip_prefix("m:") {
            for position in body.split(',') {
                let (source, targets) = position
                    .split_once(':')
                    .ok_or("invalid position variation")?;
                normalized
                    .positions
                    .push((atom(source)?, atoms(targets, '.')?));
            }
        } else if let Some(body) = segment.strip_prefix("Sg:n:") {
            let (members, details) = body.split_once(':').ok_or("invalid polymer unit")?;
            // Retain the complete suffix, including connectivity and brackets.
            normalized
                .polymers
                .push((atoms(members, ',')?, details.to_owned()));
        } else if let Some((group, members)) = segment.split_once(':') {
            let kind = group.chars().next().ok_or("empty enhanced stereo group")?;
            if !matches!(kind, 'a' | 'o' | '&')
                || (kind == 'a' && group != "a")
                || (kind != 'a' && group[1..].parse::<u32>().is_err())
            {
                return Err(format!("unsupported semantic CX layer {segment}"));
            }
            // Group numbers are names; kind and atom membership are semantics.
            normalized.stereo.push((kind, atoms(members, ',')?));
        } else {
            return Err(format!("unsupported semantic CX layer {segment}"));
        }
    }
    normalized.stereo.sort();
    normalized.positions.sort();
    normalized.polymers.sort();
    Ok(normalized)
}

fn split_labels(text: &str) -> Result<Vec<&str>, String> {
    let (mut start, mut index) = (0, 0);
    let mut labels = Vec::new();
    let bytes = text.as_bytes();
    while index < bytes.len() {
        if bytes[index..].starts_with(b"&#") {
            let end = bytes[index + 2..]
                .iter()
                .position(|byte| *byte == b';')
                .map(|offset| index + 2 + offset)
                .ok_or("unterminated escaped CX label")?;
            if end == index + 2 || !bytes[index + 2..end].iter().all(u8::is_ascii_digit) {
                return Err("invalid escaped CX label".into());
            }
            index = end;
        } else if bytes[index] == b';' {
            labels.push(&text[start..index]);
            start = index + 1;
        }
        index += 1;
    }
    labels.push(&text[start..]);
    Ok(labels)
}
