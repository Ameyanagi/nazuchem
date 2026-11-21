// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: ReShiki | src/naming.rs | 4460ac34c3ad6bf79358c4e4ae07103f69e5f42d
// Source: ReShiki | src/naming/reverse.rs | 4460ac34c3ad6bf79358c4e4ae07103f69e5f42d
// Source: ReShiki | src/naming/local.rs | 4460ac34c3ad6bf79358c4e4ae07103f69e5f42d
// Source: ReShiki | src/naming/tests.rs | 4460ac34c3ad6bf79358c4e4ae07103f69e5f42d
// Source: ReShiki | crates/model/src/molecule_names.rs | 4460ac34c3ad6bf79358c4e4ae07103f69e5f42d
// Changes: move naming admission and verification onto kernel State with an injected Budget.

//! Admission and identity checks shared by naming and generated-name verification.
use nazuchem_ir::{NameReason, NameStatus};
use nazuchem_kernel::{
    self as chemistry,
    budget::Budget,
    identity::{self as molecule_names, IdentityRejection},
};

const MAX_NAMING_ATOMS: usize = 512;
const MAX_SMILES_BYTES: usize = 32_768;

fn is_tetrahedral(tag: u8) -> bool {
    matches!(tag, 1 | 2)
}

pub(crate) fn canonical_state(
    state: &chemistry::stereo::perception::State,
    budget: &Budget,
) -> Result<String, String> {
    canonical(state, &[], budget)
}

/// Odd-cumulene geometry is supplied separately by the SMILES reader.
fn canonical(
    state: &chemistry::stereo::perception::State,
    cumulenes: &[chemistry::stereo::CumuleneEz],
    budget: &Budget,
) -> Result<String, String> {
    if state.graph.atoms.is_empty() || state.graph.atoms.len() > MAX_NAMING_ATOMS {
        return Err("Naming supports molecular graphs with 1–512 atoms".into());
    }
    // The kernel checks chemical semantics; this adapter retains the naming
    // bounds, diagnostics and admission order.
    molecule_names::canonical_smiles_identity(state, cumulenes, budget).map_err(|error| match error {
        IdentityRejection::UnsupportedAtomsOrBonds => {
            "Naming does not support query atoms, radicals, exotic bonds or relative stereo groups"
                .into()
        }
        IdentityRejection::NonclassicalStereo => {
            const MESSAGE: &str = "Naming supports tetrahedral, ordinary alkene and biaryl axis stereochemistry; other stereo classes cannot be discarded";
            MESSAGE.into()
        }
        IdentityRejection::Write(message) => message,
    })
}

/// Canonical isomeric graph identity, including isotope, charge and specified
/// stereo. Reserved axis markers carry biaryl configuration; other maps and
/// atom order are not molecular identity.
pub fn canonical_smiles(text: &str, budget: &Budget) -> Result<String, String> {
    let prepared = prepare_identity(text, budget)?;
    canonical(&prepared.state, &prepared.cumulene_ez, budget)
}

fn prepare_identity(text: &str, budget: &Budget) -> Result<chemistry::smiles::Prepared, String> {
    if text.len() > MAX_SMILES_BYTES || text.contains('|') || text.contains('>') {
        return Err(
            "Naming requires ordinary molecular SMILES without reaction or CX extensions".into(),
        );
    }
    bound_atom_tokens(text)?;
    // Native import cleans unusable winding and can discard non-tetrahedral
    // classes. Check the raw graph first, before that cleanup changes identity.
    // Unlike canonical_identity, this raw check intentionally covers atoms only.
    let raw = chemistry::smiles::parse(text).map_err(|e| e.to_string())?;
    if raw.metadata.atoms.iter().any(|m| m.chiral_tag > 2) {
        return Err("Naming does not support allene or coordination stereo classes; no plain graph was substituted".into());
    }
    let specified_atoms = raw
        .metadata
        .atoms
        .iter()
        .filter(|m| is_tetrahedral(m.chiral_tag))
        .count();
    let mut imported = chemistry::smiles::read(text, budget).map_err(|e| e.to_string())?;
    if imported.dropped_cumulene_ez {
        return Err("Naming does not support this cumulene E/Z configuration yet; no less specific name was substituted".into());
    }
    if specified_atoms
        != imported
            .prepared
            .state
            .metadata
            .atoms
            .iter()
            .filter(|m| is_tetrahedral(m.chiral_tag))
            .count()
    {
        return Err("Native import could not retain every specified tetrahedral center; no name or graph was assigned".into());
    }
    imported.prepared.state = molecule_names::apply_axis_markers(&imported.prepared.state, budget)
        .map_err(|_| "Naming received invalid biaryl axis markers".to_owned())?;
    Ok(imported.prepared)
}

/// Reject oversized output before native sanitization/ranking can allocate or
/// expand a graph. Bracket atoms count once; Cl/Br are two-character atoms.
fn bound_atom_tokens(text: &str) -> Result<(), String> {
    let bytes = text.as_bytes();
    let mut i = 0;
    let mut atoms = 0usize;
    while let Some(&b) = bytes.get(i) {
        if b == b'[' {
            atoms += 1;
            let end = bytes
                .get(i + 1..)
                .and_then(|s| s.iter().position(|&b| b == b']'))
                .ok_or("Unclosed naming SMILES atom")?;
            i += end + 2;
        } else {
            if b.is_ascii_alphabetic() || b == b'*' {
                atoms += 1;
                if (b == b'C' && bytes.get(i + 1) == Some(&b'l'))
                    || (b == b'B' && bytes.get(i + 1) == Some(&b'r'))
                {
                    i += 1;
                }
            }
            i += 1;
        }
        if atoms > MAX_NAMING_ATOMS {
            return Err(
                "Naming supports at most 512 input atom tokens; no oversized graph was imported"
                    .into(),
            );
        }
    }
    Ok(())
}

/// Verify a nomenclature round trip. P-72.2.2.2.1.1 gives a thiocarboxylate
/// the same name for C(=O)S− and C(=S)O−. The graph, including isotope,
/// hydrogen placement, stereo and charges outside those groups, must match.
pub fn verify_identity(expected: &str, actual: &str, budget: &Budget) -> Result<(), String> {
    let mut expected = prepare_identity(expected, budget)?;
    let expected_exact = canonical(&expected.state, &expected.cumulene_ez, budget)?;
    let mut actual = prepare_identity(actual, budget)?;
    let actual_exact = canonical(&actual.state, &actual.cumulene_ez, budget)?;
    if expected_exact == actual_exact {
        return Ok(());
    }
    let expected_normalized = normalize_thiocarboxylates(&mut expected, budget)?;
    let actual_normalized = normalize_thiocarboxylates(&mut actual, budget)?;
    if (expected_normalized.is_some() || actual_normalized.is_some())
        && expected_normalized.as_ref().unwrap_or(&expected_exact)
            == actual_normalized.as_ref().unwrap_or(&actual_exact)
    {
        return Ok(());
    }
    Err("The source structure differs in connectivity, charge, isotope, tautomer or stereochemistry. No name was assigned".into())
}

fn normalize_thiocarboxylates(
    prepared: &mut chemistry::smiles::Prepared,
    budget: &Budget,
) -> Result<Option<String>, String> {
    let state = &mut prepared.state;
    let graph = &mut state.graph;
    if !graph
        .atoms
        .iter()
        .any(|atom| atom.atomic_number == 16 && atom.charge == -1)
    {
        return Ok(None);
    }
    let mut adjacent = vec![Vec::new(); graph.atoms.len()];
    for (id, bond) in graph.bonds.iter().enumerate() {
        adjacent[bond.a].push((bond.b, id));
        adjacent[bond.b].push((bond.a, id));
    }
    let mut groups = Vec::new();
    for (carbon, atom) in graph.atoms.iter().enumerate() {
        if atom.atomic_number != 6
            || atom.charge != 0
            || atom.aromatic
            || atom.radical_electrons != 0
        {
            continue;
        }
        let hydrogens = |id: usize| {
            u32::from(graph.atoms[id].explicit_hydrogens) + state.valences[id].implicit_hydrogens
        };
        // R is carbon or hydrogen in methanethioate. Isotopic hydrogen
        // survives preparation as a separate vertex.
        if !matches!((adjacent[carbon].len(), hydrogens(carbon)), (2, 1) | (3, 0)) {
            continue;
        }
        let terminal = |element, charge, order| {
            adjacent[carbon].iter().copied().find(|&(id, bond)| {
                let other = &graph.atoms[id];
                other.atomic_number == element
                    && other.charge == charge
                    && other.radical_electrons == 0
                    && !other.aromatic
                    && adjacent[id].len() == 1
                    && hydrogens(id) == 0
                    && graph.bonds[bond].order == order
                    && state.metadata.bonds[bond].stereo == 0
                    && state.directions[bond] == chemistry::kekulize::Direction::None
            })
        };
        if let (Some(oxygen), Some(sulfur)) = (terminal(8, 0, 2), terminal(16, -1, 1))
            && adjacent[carbon].iter().all(|&(id, bond)| {
                let other = &graph.atoms[id];
                id == oxygen.0
                    || id == sulfur.0
                    || graph.bonds[bond].order == 1
                        && (other.atomic_number == 6
                            || other.atomic_number == 1
                                && other.charge == 0
                                && other.radical_electrons == 0
                                && adjacent[id].len() == 1)
            })
        {
            groups.push((oxygen, sulfur));
        }
    }
    if groups.is_empty() {
        return Ok(None);
    }
    for ((oxygen, co), (sulfur, cs)) in groups {
        graph.atoms[oxygen].charge = -1;
        graph.atoms[sulfur].charge = 0;
        graph.bonds[co].order = 1;
        graph.bonds[cs].order = 2;
    }
    state.valences = graph.valences()?;
    canonical(state, &prepared.cumulene_ez, budget).map(Some)
}

/// Verify a generated alternative, allowing either enantiomer only for a
/// validated relative endo/exo name.
pub fn verify_alternative_identity(
    expected: &str,
    actual: &str,
    status: NameStatus,
    budget: &Budget,
) -> Result<(), String> {
    if status != NameStatus::NotPreferred(NameReason::EndoExoStereo) {
        return verify_identity(expected, actual, budget);
    }
    // These relative descriptors identify either enantiomer of the cage.
    // Invert every tetrahedral centre together.
    let mut mirror = String::with_capacity(expected.len());
    let mut chars = expected.chars().peekable();
    while let Some(c) = chars.next() {
        mirror.push(c);
        if c == '@' && chars.next_if_eq(&'@').is_none() {
            mirror.push('@');
        }
    }
    let relative = relative_stereo_smiles(actual);
    let chiral = canonical_smiles(expected, budget)? != canonical_smiles(&mirror, budget)?;
    if chiral != relative.is_some() {
        return Err("The reconstructed endo/exo name must carry a relative stereo group exactly when the source structure is chiral".into());
    }
    if relative.is_some() {
        let imported =
            chemistry::smiles::read(actual, budget).map_err(|error| error.to_string())?;
        let metadata = &imported.prepared.state.metadata;
        let valid_group = match metadata.groups.as_slice() {
            [group] if group.kind == 1 && group.bonds.is_empty() => {
                let mut atoms = group.atoms.clone();
                atoms.sort_unstable();
                atoms.into_iter().eq(metadata
                    .atoms
                    .iter()
                    .enumerate()
                    .filter_map(|(atom, meta)| is_tetrahedral(meta.chiral_tag).then_some(atom)))
            }
            _ => false,
        };
        if !valid_group {
            return Err("The reconstructed endo/exo name must carry one relative stereo group containing exactly every tetrahedral stereocentre and no bonds".into());
        }
    }
    let actual = relative.unwrap_or(actual);
    verify_identity(expected, actual, budget).or_else(|_| verify_identity(&mirror, actual, budget))
}

/// Check only the syntax of a single relative CX group in a generation reply.
/// Reverse verification checks its parsed membership before accepting identity.
pub(crate) fn relative_stereo_smiles(smiles: &str) -> Option<&str> {
    let (graph, group) = smiles.split_once(" |o")?;
    let (number, atoms) = group.strip_suffix('|')?.split_once(':')?;
    let digits = |text: &str| !text.is_empty() && text.bytes().all(|b| b.is_ascii_digit());
    (digits(number) && atoms.split(',').all(digits)).then_some(graph)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn methyl_radical_with_hydrogen_vertices(isotopes: &[Option<u16>]) -> nazuchem_ir::Molecule {
        use nazuchem_ir::{Atom, Bond, BondOrder, Molecule};

        let carbon = isotopes.len() as u32;
        let mut atoms = isotopes
            .iter()
            .map(|&isotope| Atom {
                element: 1,
                isotope,
                ..Atom::default()
            })
            .collect::<Vec<_>>();
        atoms.push(Atom {
            element: 6,
            hydrogens: 3 - isotopes.len() as u8,
            hydrogens_stated: true,
            radical_electrons: 1,
            ..Atom::default()
        });
        Molecule {
            atoms,
            bonds: (0..carbon)
                .map(|hydrogen| Bond {
                    atoms: [hydrogen, carbon],
                    order: BondOrder::Single,
                })
                .collect(),
            ..Molecule::default()
        }
    }

    #[test]
    fn radical_identity_folds_ordinary_h_vertices_without_changing_caller_identities() {
        let budget = Budget::default();
        let expected = crate::name_smiles("[CH3]", &budget).unwrap();
        assert_eq!(expected.text, "methyl");
        assert_eq!(
            serde_json::to_value(crate::name_smiles("[H][CH2]", &budget).unwrap()).unwrap(),
            serde_json::to_value(&expected).unwrap()
        );
        // Construct IR directly: the SMILES fixture bridge has already removed
        // ordinary H vertices and cannot exercise this admission boundary.
        for count in 1..=3 {
            let source = methyl_radical_with_hydrogen_vertices(&vec![None; count]);
            for reversed in [false, true] {
                let mut molecule = source.clone();
                if reversed {
                    molecule.atoms.reverse();
                    for bond in &mut molecule.bonds {
                        bond.atoms = bond.atoms.map(|atom| count as u32 - atom);
                    }
                }
                let original = molecule.clone();
                let state = crate::ir_input::prepare(&molecule, &budget).unwrap();
                let state_before = serde_json::to_value(&state).unwrap();
                let canonical = canonical_state(&state, &budget).unwrap();
                assert_eq!(canonical_smiles(&canonical, &budget).unwrap(), "[CH3]");
                assert_eq!(serde_json::to_value(&state).unwrap(), state_before);
                let carbon = if reversed { 0 } else { count };
                let graph = crate::graph_io::read_molecule(&molecule, &budget).unwrap();
                assert_eq!(graph.atoms.keys().copied().collect::<Vec<_>>(), [carbon]);
                assert_eq!(graph.atoms[&carbon].total_h_count, 3);
                assert_eq!(graph.atoms[&carbon].radical_electrons, 1);
                let actual =
                    crate::name_molecule(&molecule, nazuchem_ir::IsotopeStyle::default(), &budget)
                        .unwrap();
                assert_eq!(
                    serde_json::to_value(actual).unwrap(),
                    serde_json::to_value(&expected).unwrap(),
                    "{count} H vertices, reversed={reversed}"
                );
                assert_eq!(molecule, original);
            }
        }
    }

    #[test]
    fn radical_identity_retains_isotopic_h_while_folding_ordinary_h() {
        let budget = Budget::default();
        // Issue #27, P-82.2.1: the radical cites its label as the prefix does.
        let expected = crate::name_smiles("[2H][CH2]", &budget).unwrap();
        assert_eq!(expected.text, "(2H1)methyl");
        for isotopes in [vec![Some(2)], vec![None, Some(2)]] {
            let molecule = methyl_radical_with_hydrogen_vertices(&isotopes);
            let original = molecule.clone();
            let carbon = isotopes.len();
            let isotope = carbon - 1;
            let state = crate::ir_input::prepare(&molecule, &budget).unwrap();
            let state_before = serde_json::to_value(&state).unwrap();
            let canonical = canonical_state(&state, &budget).unwrap();
            assert_eq!(
                canonical_smiles(&canonical, &budget).unwrap(),
                canonical_smiles("[2H][CH2]", &budget).unwrap()
            );
            assert_eq!(serde_json::to_value(&state).unwrap(), state_before);
            let graph = crate::graph_io::read_molecule(&molecule, &budget).unwrap();
            assert_eq!(
                graph.atoms.keys().copied().collect::<Vec<_>>(),
                [isotope, carbon]
            );
            assert_eq!(graph.atoms[&isotope].isotope, Some(2));
            assert_eq!(graph.atoms[&carbon].total_h_count, 2);
            assert_eq!(graph.atoms[&carbon].radical_electrons, 1);
            let actual =
                crate::name_molecule(&molecule, nazuchem_ir::IsotopeStyle::default(), &budget)
                    .unwrap();
            assert_eq!(
                serde_json::to_value(actual).unwrap(),
                serde_json::to_value(&expected).unwrap()
            );
            assert_eq!(molecule, original);
        }
    }

    #[test]
    fn oversized_smiles_are_rejected_before_graph_preparation() {
        assert!(bound_atom_tokens(&"C".repeat(513)).is_err());
        assert!(bound_atom_tokens(&"[13CH3]".repeat(513)).is_err());
        assert!(bound_atom_tokens(&"Cl".repeat(512)).is_ok());
        assert!(bound_atom_tokens(&"Cc".repeat(257)).is_err());
        assert!(bound_atom_tokens("C[CH3").is_err());
    }
}
