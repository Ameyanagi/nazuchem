// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | crates/nazuchem-s2n/src/multiplicative.rs | issue #4
// Copyright (c) 2026 Ameyanagi and Nazuchem contributors

//! Multiplicative names for identical units on one single-atom linking group
//! (P-15.3.1.3, P-51.3.1).
use super::{
    chains::get_cyclic_atoms,
    component_namer,
    molecule::{AtomSet, Molecule},
    subgraph_tools::subgraph_component,
};
use nazuchem_ir::{NameReason, NameStatus};

/// Two or three identical units on a simple linking atom, where each unit is
/// an all-carbon ring or functionalized chain bearing nothing but suffix
/// groups. Everything else keeps its substitutive name.
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
    let unit = unit_name(mol, component);
    let clean = mol.name_status() == NameStatus::Unverified(NameReason::UnauditedPath);
    mol.set_name_status(status);
    let (locant, linker, count, unit) = unit.filter(|_| clean)?;
    let locants: Vec<_> = (0..count)
        .map(|i| format!("{locant}{}", "'".repeat(i)))
        .collect();
    Some(format!(
        "{}-{linker}{}{unit}",
        locants.join(","),
        ["di", "tri"][count - 2]
    ))
}

/// The attachment locant, the linking prefix, the number of units and the
/// name of one unit, enclosed where needed, as the substitutive name numbers it.
fn unit_name(mol: &Molecule, component: &AtomSet) -> Option<(String, &'static str, usize, String)> {
    // The unit keeps the locants it needs beside the linking group
    // (P-15.3.2.1), as in di(cyclohexane-1-carboxylic acid).
    let parent = component_namer::name_component_with_atom_locants(mol, component, None).ok()?;
    let component_namer::Derivation::Plan {
        parent_atoms,
        principal_group,
        ..
    } = &parent.derivation
    else {
        return None;
    };
    // The senior parent and its suffix groups carry one substituent: the
    // linking atom with the other units.
    let mut unit = parent_atoms.clone();
    unit.extend(principal_group.iter().flat_map(|(_, atoms)| atoms));
    let links: Vec<_> = unit
        .iter()
        .flat_map(|&u| {
            mol.get_neighbors(u)
                .iter()
                .filter(|n| !unit.contains(n) && mol.atoms[n].symbol != "H")
                .map(move |&n| (u, n))
        })
        .collect();
    let [(root, centre)] = links[..] else {
        return None;
    };
    let locant = parent.locants.get(&root)?;
    let cyclic = get_cyclic_atoms(mol, &AtomSet::new());
    // P-15.3.1.1, P-15.3.4.2: rings and functionalized chains, not alkanes.
    // Heteroatom parents need the 'bis' rules of P-16.3.6 and are left alone.
    if parent_atoms.len() < 2
        || parent_atoms.iter().any(|i| !mol.atoms[i].is_carbon())
        || !(cyclic.contains(&root) || principal_group.is_some())
        || !locant.bytes().all(|b| b.is_ascii_digit())
    {
        return None;
    }
    let atom = &mol.atoms[&centre];
    let (oxo, arms): (Vec<usize>, Vec<usize>) = mol
        .get_neighbors(centre)
        .iter()
        .filter(|n| mol.atoms[n].symbol != "H")
        .partition(|&&n| {
            mol.atoms[&n].symbol == "O" && mol.degree(n) == 1 && mol.bond_order(centre, n) == 2
        });
    // P-15.3.1.2.1.1 (methylene, oxy, sulfanediyl, nitrilo), P-15.3.1.2.1.2
    // (azanediyl), P-51.3.1 (phosphanetriyl), P-63.2.5 (selanediyl,
    // tellanediyl), P-63.6 and P-35.2.3 (sulfinyl to telluronyl).
    let linker = match (
        atom.symbol.as_str(),
        oxo.len(),
        arms.len(),
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
        ("N", 0, 2, 1) => "azanediyl",
        ("N", 0, 3, 0) => "nitrilo",
        ("P", 0, 3, 0) => "phosphanetriyl",
        _ => return None,
    };
    if cyclic.contains(&centre) || arms.iter().any(|&arm| mol.bond_order(centre, arm) != 1) {
        return None;
    }
    // P-51.3.1 (3): every unit matches the named one from its attachment atom.
    let exclude = AtomSet::from([centre]);
    let unit = subgraph_component(mol, root, &exclude);
    if !arms.iter().filter(|&&arm| arm != root).all(|&arm| {
        let other = subgraph_component(mol, arm, &exclude);
        super::identical_multiplicative_branches(mol, &unit, &other, &unit, &other, root, arm)
    }) {
        return None;
    }
    let text = parent.rendered_parent;
    let stem = text.strip_suffix(" acid").unwrap_or(&text);
    let functionalized = principal_group.is_some();
    let numerical = [
        "di", "tri", "tetr", "pent", "hex", "hept", "oct", "non", "dec", "undec", "do", "hen",
        "icos",
    ]
    .iter()
    .any(|term| stem.starts_with(term));
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
    Some((locant.clone(), linker, arms.len(), unit))
}

#[cfg(test)]
mod tests {
    use crate::test_support::{parse_strict, verify_identity};
    use nazuchem_ir::{NameReason, NameStatus};

    fn name(smiles: &str) -> crate::Name {
        crate::name_smiles(smiles, &Default::default()).unwrap()
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
                "(carboxymethylsulfanyl)acetic acid",
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
            verify_identity(smiles, &parse_strict(substitutive)?)?;
        }
        Ok(())
    }

    #[test]
    fn other_multiplicative_pins_keep_the_substitutive_name_and_its_status() {
        for (smiles, expected) in [
            // P-15.3.2.4.1: substituted units take 'bis' and their own numbering.
            (
                "Brc1ccc(Oc2ccc(Br)cc2)cc1",
                "1-bromo-4-(4-bromophenoxy)benzene",
            ),
            // P-15.3.1.2.2, P-52.2.5.3 Example 2: concatenated linking groups.
            ("c1ccccc1Oc2ccc(Oc3ccccc3)cc2", "1,4-diphenoxybenzene"),
            (
                "O=C(O)c1ccc(OCCOc2ccc(C(=O)O)cc2)cc1",
                "4-[2-(4-carboxyphenoxy)ethoxy]benzoic acid",
            ),
            // P-15.3.3.1: chains as linking groups.
            ("CC(c1ccccc1)Cc1ccccc1", "(2-phenylpropyl)benzene"),
            // P-16.3.6 (c): heterocyclic units.
            ("c1ccncc1Oc1cccnc1", "3-[(pyridin-3-yl)oxy]pyridine"),
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
