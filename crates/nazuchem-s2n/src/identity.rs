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
    canonical(state, false, budget)
}

/// `cumulenes` keeps the odd-cumulene E/Z markers of a state read from SMILES.
fn canonical(
    state: &chemistry::stereo::perception::State,
    cumulenes: bool,
    budget: &Budget,
) -> Result<String, String> {
    if state.graph.atoms.is_empty() || state.graph.atoms.len() > MAX_NAMING_ATOMS {
        return Err("Naming supports molecular graphs with 1–512 atoms".into());
    }
    // The kernel checks chemical semantics; this adapter retains the naming
    // bounds, diagnostics and admission order.
    if cumulenes {
        molecule_names::canonical_smiles_identity(state, budget)
    } else {
        molecule_names::canonical_identity(state, budget)
    }
    .map_err(|error| match error {
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
    if imported.dropped_cumulene_ez
        || molecule_names::mark_cumulenes(&mut imported.prepared).is_err()
    {
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
    let state = molecule_names::apply_axis_markers(&imported.prepared.state, budget)
        .map_err(|_| "Naming received invalid biaryl axis markers".to_owned())?;
    canonical(&state, true, budget)
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

/// Require exact connectivity, charge, isotope, tautomer and stereo identity.
pub fn verify_identity(expected: &str, actual: &str, budget: &Budget) -> Result<(), String> {
    if canonical_smiles(expected, budget)? != canonical_smiles(actual, budget)? {
        return Err("The source structure differs in connectivity, charge, isotope, tautomer or stereochemistry. No name was assigned".into());
    }
    Ok(())
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
