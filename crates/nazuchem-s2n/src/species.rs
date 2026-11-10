// SPDX-License-Identifier: MIT AND (MIT OR Apache-2.0)
// Origin: mixed
// Source: Nazuchem | crates/nazuchem-s2n/src/species.rs | issue #25
// Source: ReShiki | src/naming/openclatura/species.rs | 5f559fd482e237307a589dfcb7ea333cee0a6e13
// Copyright (c) 2026 lamalab-org
// Source: openclatura | tests/exploratory/test_charge_representation_gaps.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/data/parser_xml_resources/simpleGroups.json | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/unit/test_chalcogen_structural_roles.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/unit/test_chalcogen_functional_groups.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/unit/test_analysis.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/data/parser_xml_resources/arylGroups.json | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/data/parser_xml_resources/functionalTerms.json | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/data/parser_xml_resources/nonCarboxylicAcids.json | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Copyright (c) 2026 Ameyanagi and ReShiki contributors

//! Binary salt names and monatomic metal cations, outside the organic parent engine.
use super::{Name, component_namer::Derivation, molecule::Molecule, name_smiles_with_style};
use nazuchem_ir::{IsotopeStyle, NameReason, NameStatus};
use std::collections::BTreeMap;

// Provenance: fn monatomic_cation | original
// Source: ReShiki | src/naming/openclatura/species.rs::monatomic_cation
pub(super) fn monatomic_cation(
    smiles: &str,
    in_salt: bool,
    budget: &nazuchem_kernel::budget::Budget,
) -> Result<Option<Name>, String> {
    let parsed = nazuchem_kernel::smiles::prepare(smiles, budget).map_err(|e| e.to_string())?;
    let [atom] = parsed.state.graph.atoms.as_slice() else {
        return Ok(None);
    };
    if atom.charge <= 0
        || atom.isotope != 0
        || atom.radical_electrons != 0
        || atom.explicit_hydrogens != 0
    {
        return Ok(None);
    }
    let (name, common_charge) = match atom.atomic_number {
        3 => ("lithium", Some(1)),
        4 => ("beryllium", Some(2)),
        11 => ("sodium", Some(1)),
        12 => ("magnesium", Some(2)),
        19 => ("potassium", Some(1)),
        20 => ("calcium", Some(2)),
        26 => ("iron", None),
        27 => ("cobalt", None),
        28 => ("nickel", None),
        29 => ("copper", None),
        30 => ("zinc", Some(2)),
        37 => ("rubidium", Some(1)),
        38 => ("strontium", Some(2)),
        47 => ("silver", Some(1)),
        48 => ("cadmium", Some(2)),
        55 => ("caesium", Some(1)),
        56 => ("barium", Some(2)),
        80 => ("mercury", None),
        82 => ("lead", None),
        87 => ("francium", Some(1)),
        88 => ("radium", Some(2)),
        _ => return Ok(None),
    };
    Ok(Some(Name {
        text: if in_salt && common_charge == Some(atom.charge) {
            name.into()
        } else {
            format!("{name}({}+)", atom.charge)
        },
        status: if in_salt {
            NameStatus::Preferred
        } else {
            NameStatus::Unverified(NameReason::UnauditedPath)
        },
        alternatives: vec![],
        omitted_stereo: vec![],
        omitted_units: 0,
        verification: None,
    }))
}

/// P-72.7: assess the anionic parent actually cited by the component pipeline.
// Provenance: fn assess_anion_parent | original
// Source: ReShiki | src/naming/openclatura/species.rs::assess_anion_parent
pub(super) fn assess_anion_parent(mol: &Molecule, derivation: &Derivation) {
    let Derivation::Plan {
        parent_atoms,
        principal_group,
        retained_anion_name,
        ..
    } = derivation
    else {
        return;
    };
    let anionic_principal = principal_group.as_ref().is_some_and(|(key, atoms)| {
        super::perception::is_anionic_group_key(key)
            && atoms.iter().map(|i| mol.atoms[i].charge).sum::<i32>() < 0
    });
    // Carbon-free anions take retained or additive names, not a P-72.7 parent.
    if mol.atoms.values().all(|atom| atom.symbol != "C")
        || !mol.atoms.values().any(|atom| atom.charge < 0)
    {
        return;
    }
    let neutral_group_atoms = neutral_group_atoms(mol);
    let anionic_centre = |i: &usize| mol.atoms[i].charge < 0 && !neutral_group_atoms.contains(i);
    if !mol.atoms.keys().any(anionic_centre) {
        return;
    }
    let principal_atoms = principal_group.as_ref().map(|(_, atoms)| atoms);
    // P-70.2 / P-74.2: only an anionic suffix group cites the charge; a
    // charge-separated or deprotonated neutral class, as in acetamide for
    // acetylazanide, does not. With no characteristic group, the parent cites -ide.
    let anionic_parent = principal_group.as_ref().map_or_else(
        || parent_atoms.iter().any(anionic_centre),
        |(key, _)| {
            anionic_principal
                // A cationic characteristic suffix can accompany skeletal -ide.
                || (matches!(key.as_str(), "aminium" | "iminium" | "diazonio")
                    && parent_atoms.iter().any(anionic_centre))
        },
    );
    let unresolved_alkoxide = matches!(
        retained_anion_name,
        super::assembly_parts::RetainedAnionName::Required { .. }
            | super::assembly_parts::RetainedAnionName::SubstitutedAlkyl
    );
    let external_centre = mol.atoms.keys().any(|i| {
        anionic_centre(i)
            && !parent_atoms.contains(i)
            && principal_atoms.is_none_or(|atoms| !atoms.contains(i))
    });
    let unresolved_parent = external_centre
        && !super::parent_selection::anion_parent_is_senior(mol, parent_atoms, principal_atoms);
    if parent_atoms.iter().any(anionic_centre)
        && principal_group
            .as_ref()
            .is_some_and(|(key, _)| matches!(key.as_str(), "aminium" | "iminium" | "diazonio"))
    {
        mol.mark_name(NameStatus::Unverified(NameReason::IonicSuffixOrder));
    }
    if !anionic_parent || unresolved_alkoxide || unresolved_parent {
        let status = NameStatus::Unverified(NameReason::AnionParent);
        // The unresolved principal anion spelling remains the reason when a
        // separately staged cationic prefix also needs verification.
        if unresolved_alkoxide
            && mol.name_status() == NameStatus::Unverified(NameReason::CationPrefix)
        {
            mol.set_name_status(status);
        } else {
            mol.mark_name(status);
        }
    }
}

/// P-41, P-73.1: a cation is senior to every neutral class, so its parent must
/// carry a cationic centre. A name that cites every cationic centre in a prefix,
/// as in (1-methylpyridin-1-ium-4-yl)methanol, is not certified.
// Provenance: fn assess_cation_parent | original
// Source: ReShiki | src/naming/openclatura/species.rs::assess_cation_parent
pub(super) fn assess_cation_parent(mol: &Molecule, derivation: &Derivation) {
    let Derivation::Plan {
        parent_atoms,
        principal_group,
        ..
    } = derivation
    else {
        return;
    };
    let neutral_group_atoms = neutral_group_atoms(mol);
    let centres: Vec<_> = mol
        .atoms
        .iter()
        .filter(|(i, atom)| atom.charge != 0 && !neutral_group_atoms.contains(i))
        .collect();
    // Zwitterions and anions take the anionic parent (P-72, P-74).
    if centres.is_empty() || centres.iter().any(|(_, atom)| atom.charge < 0) {
        return;
    }
    let cited = |i: &usize| {
        parent_atoms.contains(i)
            || principal_group
                .as_ref()
                .is_some_and(|(_, atoms)| atoms.contains(i))
    };
    if !centres.iter().any(|(i, _)| cited(i)) {
        mol.mark_name(NameStatus::Unverified(NameReason::CationParent));
    }
}

/// Nitro, N-oxide, azido, diazo and isocyano charges express neutral groups
/// rather than ionic centres. A separate counter-ion does not cancel a real centre.
// Provenance: fn neutral_group_atoms | original
// Source: ReShiki | src/naming/openclatura/species.rs::neutral_group_atoms
fn neutral_group_atoms(mol: &Molecule) -> std::collections::BTreeSet<usize> {
    super::perception::perceive_groups(mol)
        .into_iter()
        .filter(|group| {
            (matches!(
                group.key.as_str(),
                "nitro"
                    | "azido"
                    | "diazo"
                    | "acid_azide"
                    | "ring_acid_azide"
                    | "isocyano"
                    | "acid_isocyanide"
                    | "ring_acid_isocyanide"
            ) || group
                .descriptor
                .as_ref()
                .is_some_and(|descriptor| descriptor.derivative == "zwitterion"))
                && group
                    .atoms_involved
                    .iter()
                    .map(|i| mol.atoms[i].charge)
                    .sum::<i32>()
                    == 0
        })
        .flat_map(|group| group.atoms_involved)
        .collect()
}

/// Nonitalic Roman letters, then locants (italic element locants before numerals).
type SaltIonKey = (String, Option<Vec<(bool, usize, String)>>);

/// P-14.5 key of a plainly spelled ion name: Roman letters first, then locants
/// at the first difference. `None` for any italic descriptor other than element
/// locants and indicated hydrogen (tert, sec, the s/as of indacene, fusion
/// letters and other descriptor brackets), lambda, superscripts, marks other
/// than hyphens, commas, primes and enclosing marks, and composite locants.
/// Locants are compared only in names without enclosing marks, where every
/// locant is a primary one (P-14.5).
// Provenance: fn salt_ion_sort_key | original
// Source: ReShiki | src/naming/openclatura/species.rs::salt_ion_sort_key
fn salt_ion_sort_key(name: &str) -> Option<SaltIonKey> {
    // Enclosing brackets always hold parentheses; other brackets are descriptors.
    let descriptor_bracket = name
        .split('[')
        .skip(1)
        .any(|rest| !rest.split(']').next().unwrap_or(rest).contains('('));
    if !name
        .chars()
        .all(|c| c.is_ascii_alphanumeric() || " ,-'()[]{}".contains(c))
        || name.contains("lambda")
        || name.contains("indacen")
        || descriptor_bracket
    {
        return None;
    }
    let (mut letters, mut locants) = (String::new(), Vec::new());
    let primary_locants = !name.contains('(');
    for token in name.split(|c: char| !c.is_ascii_alphanumeric() && c != '\'') {
        let digits = token
            .find(|c: char| !c.is_ascii_digit())
            .unwrap_or(token.len());
        if digits > 0 {
            // 1H and 4aH cite indicated hydrogen at locants 1 and 4a.
            let rest = token[digits..]
                .strip_suffix('H')
                .unwrap_or(&token[digits..]);
            let letter = rest.trim_end_matches('\'');
            if letter.len() > 1 || !letter.bytes().all(|b| b.is_ascii_lowercase()) {
                return None;
            }
            locants.push((true, token[..digits].parse().ok()?, rest.to_owned()));
        } else if matches!(
            token.trim_end_matches('\''),
            "N" | "O" | "S" | "Se" | "Te" | "P"
        ) {
            locants.push((false, 0, token.to_owned()));
        } else if token.starts_with(|c: char| c.is_ascii_uppercase())
            || token.contains('\'')
            || token.len() == 1
            || matches!(token, "tert" | "sec")
        {
            return None;
        } else {
            letters.push_str(token);
        }
    }
    Some((letters, primary_locants.then_some(locants)))
}

// Provenance: fn name_salt | original
// Source: ReShiki | src/naming/openclatura/species.rs::name_salt
pub(super) fn name_salt(
    smiles: &str,
    style: IsotopeStyle,
    budget: &nazuchem_kernel::budget::Budget,
) -> Result<Name, String> {
    let mut ions = Vec::new();
    let mut charge = 0i32;
    for component in smiles.split('.') {
        let state = nazuchem_kernel::smiles::prepare(component, budget)
            .map_err(|e| e.to_string())?
            .state;
        let positive = state.graph.atoms.iter().any(|a| a.charge > 0);
        let negative = state.graph.atoms.iter().any(|a| a.charge < 0);
        let q: i32 = state.graph.atoms.iter().map(|a| i32::from(a.charge)).sum();
        if q == 0
            || positive && negative
            || state.graph.atoms.iter().any(|a| a.radical_electrons != 0)
        {
            return Err("Salt naming requires charged components without zwitterions or radicals; neutral components and mixtures are unsupported".into());
        }
        charge += q;
        let name = match monatomic_cation(component, true, budget)? {
            Some(name) => name,
            None => name_smiles_with_style(component, style, budget)?,
        };
        // P-14.5 leaves out stereodescriptors and isotope descriptors, which
        // the cited names cannot omit yet.
        let key = (state.graph.atoms.iter().all(|a| a.isotope == 0)
            && state.metadata.atoms.iter().all(|a| a.chiral_tag == 0)
            && state.metadata.bonds.iter().all(|b| b.stereo == 0))
        .then(|| salt_ion_sort_key(&name.text))
        .flatten();
        ions.push((q < 0, name, component, key));
    }
    // Ions of one sign are ordered by the key; when any of them has no key or
    // two different names are not told apart by it, main's character order
    // stays and is unverified.
    let ordered = ions.iter().all(|(negative, name, _, key)| {
        ions.iter().all(|other| {
            other.0 != *negative
                || other.1.text == name.text
                || match (key, &other.3) {
                    (Some((a, a_locants)), Some((b, b_locants))) => {
                        a != b
                            || a_locants.is_some() && b_locants.is_some() && a_locants != b_locants
                    }
                    _ => false,
                }
        })
    });
    let mut grouped = BTreeMap::new();
    for (negative, name, component, key) in ions {
        let plain: String = name
            .text
            .chars()
            .filter(char::is_ascii_alphanumeric)
            .map(|c| c.to_ascii_lowercase())
            .collect();
        let entry = grouped
            .entry((
                negative,
                key.filter(|_| ordered),
                plain,
                name.text.clone(),
                component,
            ))
            .or_insert((name, 0));
        entry.1 += 1;
    }
    if charge != 0 {
        return Err(
            "Salt naming requires a charge-balanced set of ions; the total charge must be zero"
                .into(),
        );
    }
    let mut status = NameStatus::Preferred;
    let mut words = Vec::new();
    for (_, (name, count)) in grouped {
        // The strongest ion downgrade wins; ties keep the first ion's reason.
        if name.status.severity() > status.severity() {
            status = name.status;
        }
        let complex = name.text.chars().any(|c| !c.is_ascii_alphabetic())
            || super::formatting::starts_with_multiplier(&name.text)
            || name.text.ends_with("aminium");
        let word = if count == 1 {
            name.text
        } else {
            let prefix = if complex {
                super::rules::multipliers::complex_(count)
            } else {
                super::rules::multipliers::basic(count)
            };
            if complex {
                format!("{prefix}({})", name.text)
            } else {
                format!("{prefix}{}", name.text)
            }
        };
        words.push(word);
    }
    if !ordered && status == NameStatus::Preferred {
        status = NameStatus::Unverified(NameReason::SaltIonOrder);
    }
    Ok(Name {
        text: super::nest_enclosing_marks(&words.join(" ")),
        status,
        alternatives: vec![],
        omitted_stereo: vec![],
        omitted_units: 0,
        verification: None,
    })
}

/// P-14.8, P-77.1.3 (1): a set with neutral components, such as a base drawn
/// with HCl or a hydrate, is named as an adduct: the component names joined by
/// em dashes, then the ratio. Organic components come first, by seniority of
/// classes (P-14.8.1, P-41; a salt is ionic, so it leads), then the inorganic
/// component, then water (P-14.8.2). The hydrohalide form of P-77.1.3 (3) and
/// the hydrate form of P-14.8.2 are alternatives.
// Provenance: fn name_adduct | original
// Source: Nazuchem | crates/nazuchem-s2n/src/species.rs::name_adduct
pub(super) fn name_adduct(
    smiles: &str,
    style: IsotopeStyle,
    budget: &nazuchem_kernel::budget::Budget,
) -> Result<Name, String> {
    const ORGANIC: u8 = 0;
    const INORGANIC: u8 = 1;
    const WATER: u8 = 2;
    struct Part {
        // Formula order: kind, then class seniority, where a component
        // without a characteristic group is last.
        rank: (u8, Vec<i64>),
        hydride: u8,
        labelled: bool,
        base: bool,
        salt: bool,
        text: String,
        status: NameStatus,
        count: usize,
    }
    let mut ions = Vec::new();
    let mut organic_ion = false;
    // Per neutral component: whether it is organic and has nitrogen, the
    // element of a mononuclear hydride (0 otherwise), whether it is
    // isotopically labelled, and its count.
    let mut neutral = BTreeMap::<&str, (bool, bool, u8, bool, usize)>::new();
    let mut charged_neutral = false;
    for component in smiles.split('.') {
        let atoms = nazuchem_kernel::smiles::prepare(component, budget)
            .map_err(|e| e.to_string())?
            .state
            .graph
            .atoms;
        let organic = atoms.iter().any(|a| a.atomic_number == 6);
        if atoms.iter().map(|a| i32::from(a.charge)).sum::<i32>() == 0 {
            let nitrogen = atoms.iter().any(|a| a.atomic_number == 7);
            let mut heavy = atoms.iter().filter(|a| a.atomic_number != 1);
            let hydride = match (heavy.next(), heavy.next()) {
                (Some(atom), None) => atom.atomic_number,
                _ => 0,
            };
            // Radicals are the most senior class (P-41) and are not ordered here.
            if atoms.iter().any(|a| a.radical_electrons != 0) {
                return Err("Adducts with radical components are unsupported".into());
            }
            let labelled = atoms.iter().any(|a| a.isotope != 0);
            charged_neutral |= atoms.iter().any(|a| a.charge != 0);
            neutral
                .entry(component)
                .or_insert((organic, nitrogen, hydride, labelled, 0))
                .4 += 1;
        } else {
            ions.push(component);
            organic_ion |= organic;
        }
    }
    if neutral.is_empty() {
        return name_salt(smiles, style, budget);
    }
    let mut parts = Vec::new();
    if !ions.is_empty() {
        // How a zwitterion pairs with separate ions is not decided here.
        if charged_neutral {
            return Err(
                "Zwitterions and charge-separated components mixed with ions are unsupported"
                    .into(),
            );
        }
        let salt = name_salt(&ions.join("."), style, budget)?;
        parts.push(Part {
            rank: (if organic_ion { ORGANIC } else { INORGANIC }, vec![]),
            hydride: 0,
            labelled: false,
            base: false,
            salt: true,
            text: salt.text,
            status: salt.status,
            count: 1,
        });
    }
    for (component, (organic, nitrogen, hydride, labelled, count)) in neutral {
        let name = name_smiles_with_style(component, style, budget)?;
        let rank = if hydride == 8 {
            (WATER, vec![])
        } else if organic {
            let mol = super::graph_io::read_smiles(component, budget)?;
            // A zwitterion is senior to every neutral class (P-41), unlike a
            // charge-separated neutral group such as nitro.
            let neutral_groups = neutral_group_atoms(&mol);
            if mol
                .atoms
                .iter()
                .any(|(i, atom)| atom.charge != 0 && !neutral_groups.contains(i))
            {
                return Err("Adducts with zwitterionic components are unsupported".into());
            }
            let groups = super::perception::perceive_groups(&mol);
            let key = super::principal_groups::component_principal_key(&mol, &groups, false);
            let seniority = groups
                .iter()
                .find(|g| Some(&g.key) == key.as_ref())
                .and_then(|g| g.resolved_rule.as_ref()?.seniority.clone());
            (ORGANIC, seniority.unwrap_or(vec![i64::MAX]))
        } else {
            (INORGANIC, vec![])
        };
        parts.push(Part {
            rank,
            hydride,
            labelled,
            base: organic && nitrogen,
            salt: false,
            text: name.text,
            status: name.status,
            count,
        });
    }
    parts.sort_by(|a, b| a.rank.cmp(&b.rank));
    if parts.len() < 2 {
        return Err("Adduct naming requires at least two different components".into());
    }
    // Organic components of one class, or without a characteristic group, are
    // ordered by the further criteria of P-41 and P-44, which are not covered.
    if parts
        .windows(2)
        .any(|pair| pair[0].rank.0 == ORGANIC && pair[0].rank == pair[1].rank)
    {
        return Err("Adducts of organic components of the same class are unsupported".into());
    }
    // Inorganic components are ordered by the Red Book, which is not covered.
    if parts.iter().filter(|p| p.rank.0 == INORGANIC).count() > 1 {
        return Err(
            "Adducts with more than one inorganic component besides water are unsupported".into(),
        );
    }
    // A mixed adduct has no preferred IUPAC name (P-14.8.2); the class order
    // of an organic adduct is not certified yet.
    let mut status = if parts.iter().any(|p| p.rank.0 != ORGANIC) {
        NameStatus::NotPreferred(NameReason::Adduct)
    } else {
        NameStatus::Unverified(NameReason::Adduct)
    };
    for part in &parts {
        if part.status.severity() > status.severity() {
            status = part.status;
        }
    }
    let multiplied = |n: usize, word: &str| match n {
        1 => word.to_owned(),
        n => format!("{}{word}", super::rules::multipliers::basic(n)),
    };
    let word = match parts.as_slice() {
        // P-14.8.2: hydrate with a numerical prefix; hemi and sesqui are also used.
        [first, water] if water.rank.0 == WATER && !water.labelled => {
            match (first.count, water.count) {
                (1, 1) => Some("monohydrate".to_owned()),
                (1, n) => Some(multiplied(n, "hydrate")),
                (2, 1) => Some("hemihydrate".to_owned()),
                (2, 3) => Some("sesquihydrate".to_owned()),
                _ => None,
            }
        }
        // P-77.1.3 (3): the unaltered name of the organic base, then the hydrohalide.
        [first, acid] if first.count == 1 && first.base && !acid.labelled => [
            (9, "hydrofluoride"),
            (17, "hydrochloride"),
            (35, "hydrobromide"),
            (53, "hydroiodide"),
        ]
        .iter()
        .find(|(halogen, _)| *halogen == acid.hydride)
        .map(|(_, word)| multiplied(acid.count, word)),
        _ => None,
    };
    let alternatives = word
        .map(|word| nazuchem_ir::AlternativeName {
            name: format!("{} {word}", parts[0].text),
            status: NameStatus::NotPreferred(NameReason::Adduct),
        })
        .into_iter()
        .collect();
    // A binary salt name is enclosed, as in {piperazin-1-ium ...carboxylate}—ethyl
    // acetate—dichloromethane (1/1/1) (P-14.8.1).
    let names: Vec<_> = parts
        .iter()
        .map(|p| match p.salt {
            true => super::nest_enclosing_marks(&format!("({})", p.text)),
            false => p.text.clone(),
        })
        .collect();
    let ratio: Vec<_> = parts.iter().map(|p| p.count.to_string()).collect();
    Ok(Name {
        text: format!("{} ({})", names.join("\u{2014}"), ratio.join("/")),
        status,
        alternatives,
        omitted_stereo: vec![],
        omitted_units: 0,
        verification: None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    // Provenance: fn tests::salt_ion_key_orders_letters_then_locants_and_declines_uncovered_spellings | adapted
    // Source: ReShiki | src/naming/openclatura/species.rs::tests::salt_ion_key_orders_letters_then_locants_and_declines_uncovered_spellings
    // Source: openclatura | tests/exploratory/test_charge_representation_gaps.py::test_wrong_charge_single_atoms_fail_honestly_instead_of_mislabeling
    // Source: openclatura | tests/data/parser_xml_resources/simpleGroups.json::<fixture data>
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn salt_ion_key_orders_letters_then_locants_and_declines_uncovered_spellings() {
        for (first, second) in [
            ("acetate", "4-methylbenzoate"),
            ("2-methylundecanoate", "10-methylundecanoate"),
            ("N-methylanilinium", "2-methylanilinium"),
            ("sodium", "N,N,N-trimethylmethanaminium"),
            ("1H-inden-1-ide", "2H-inden-2-ide"),
        ] {
            assert!(
                salt_ion_sort_key(first).unwrap() < salt_ion_sort_key(second).unwrap(),
                "{first} before {second}"
            );
        }
        for name in [
            "[1,1'-biphenyl]-4-carboxylate",
            "benzo[a]anthracen-7-ide",
            "(2H3)acetate",
            "1lambda^5-phosphinin-1-ylium",
            "N1-methylbutanediamide",
            "3-chloros-indacene-1-carboxylate",
            "3-chloros-indacen-1-olate",
            "as-indacen-1-ide",
            "tert-butoxide",
            "N-(tert-butyl)anilinium",
            "s-triazin-2-olate",
            "(2R)-butan-2-olate",
        ] {
            assert_eq!(salt_ion_sort_key(name), None, "{name}");
        }
    }

    #[test]
    // Provenance: fn tests::salt_ion_order_falls_back_for_uncovered_ions | original
    // Source: ReShiki | src/naming/openclatura/species.rs::tests::salt_ion_order_falls_back_for_uncovered_ions
    fn salt_ion_order_falls_back_for_uncovered_ions() {
        for (smiles, expected) in [
            (
                "[Na+].[Na+].CC[C@@H](C)[C@@H](O)C([O-])=O.CC[C@H](O)[C@H](C)C([O-])=O",
                "disodium (2R,3R)-2-hydroxy-3-methylpentanoate (2S,3S)-3-hydroxy-2-methylpentanoate",
            ),
            (
                "[Na+].[Na+].[O-]C(=O)c1ccc(cc1)-c1ccccc1.CC(=O)[O-]",
                "disodium [1,1'-biphenyl]-4-carboxylate acetate",
            ),
        ] {
            let name = name_salt(smiles, IsotopeStyle::Iupac, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(
                name.status,
                NameStatus::Unverified(NameReason::SaltIonOrder),
                "{smiles}"
            );
        }
        // A letter tie that needs internal locants is not ordered.
        let (first, second) = (
            salt_ion_sort_key("2-(3-chloropropyl)octanoate").unwrap(),
            salt_ion_sort_key("2-chloro-2-propyloctanoate").unwrap(),
        );
        assert_eq!(first.0, second.0);
        assert!(first.1.is_none() || second.1.is_none());
        // One stereospecified ion per sign leaves nothing to order.
        let name = name_salt(
            "C[C@@H](CC)[NH3+].CC(=O)[O-]",
            IsotopeStyle::Iupac,
            &Default::default(),
        )
        .unwrap();
        assert_eq!(name.status, NameStatus::Preferred, "{}", name.text);
    }

    #[test]
    // Provenance: fn tests::mixed_salt_anion_spellings_preserve_status | original
    // Source: ReShiki | src/naming/openclatura/species.rs::tests::mixed_salt_anion_spellings_preserve_status
    fn mixed_salt_anion_spellings_preserve_status() {
        for (smiles, expected, status) in [
            (
                "[Na+].[Na+].Cc1ccc(cc1)C(=O)[O-].CC(=O)[O-]",
                "disodium acetate 4-methylbenzoate",
                NameStatus::Preferred,
            ),
            (
                "[Na+].[Na+].[Na+].CC(=O)CC[O-].O=C([O-])c1ccc(cc1)Oc1ccc(cc1)C(=O)[O-]",
                "trisodium 4-(4-carboxylatophenoxy)benzoate 3-oxobutan-1-olate",
                NameStatus::NotPreferred(NameReason::Multiplicative),
            ),
        ] {
            let name = name_salt(smiles, IsotopeStyle::Iupac, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, status, "{smiles}");
        }
    }

    // Provenance: fn tests::assess_component | original
    // Source: ReShiki | src/naming/openclatura/species.rs::tests::assess_component
    fn assess_component(smiles: &str) -> NameStatus {
        let mol = super::super::graph_io::read_smiles(smiles, &Default::default()).unwrap();
        let atoms = mol.atoms.keys().copied().collect();
        let name =
            super::super::component_namer::name_component_with_citation(&mol, &atoms, None, None)
                .unwrap();
        assert!(
            matches!(&name.derivation, Derivation::Plan { .. }),
            "{smiles}"
        );
        // Isolate this assessor from marks made while deriving the name.
        mol.set_name_status(NameStatus::Unverified(NameReason::UnauditedPath));
        assess_anion_parent(&mol, &name.derivation);
        mol.certify();
        mol.name_status()
    }

    #[test]
    // Provenance: fn tests::anion_parent_assessment_uses_the_cited_derivation | original
    // Source: ReShiki | src/naming/openclatura/species.rs::tests::anion_parent_assessment_uses_the_cited_derivation
    fn anion_parent_assessment_uses_the_cited_derivation() {
        for smiles in ["C[N+](C)(C)CC[O-]", "[NH3+]CC([NH3+])C[O-]"] {
            assert_eq!(
                assess_component(smiles),
                NameStatus::Unverified(NameReason::AnionParent),
                "{smiles}"
            );
        }
    }

    #[test]
    // Provenance: fn tests::neutral_charge_separated_groups_are_not_anionic_centres | adapted
    // Source: ReShiki | src/naming/openclatura/species.rs::tests::neutral_charge_separated_groups_are_not_anionic_centres
    // Source: openclatura | tests/unit/test_chalcogen_structural_roles.py::test_acyl_pseudohalides_use_the_complete_leaving_group_route
    // Source: openclatura | tests/unit/test_chalcogen_functional_groups.py::test_nitrogen_chalcogenides_are_classified_before_rendering
    // Source: openclatura | tests/unit/test_analysis.py::test_azinic_acid_is_not_named_as_a_nitro_compound
    // Source: openclatura | tests/unit/test_chalcogen_functional_groups.py::test_nitrile_chalcogenides_use_zwitterion_priority_and_class_names
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn neutral_charge_separated_groups_are_not_anionic_centres() {
        for smiles in [
            "[O-][N+](=O)c1ccccc1",
            "C[N+](C)(C)[O-]",
            "CN=[N+]=[N-]",
            "C=[N+]=[N-]",
            "CC#[N+][O-]",
            "CC(=O)N=[N+]=[N-]",
            "C[N+]#[C-]",
            "CC(=O)[N+]#[C-]",
        ] {
            assert_eq!(assess_component(smiles), NameStatus::Preferred, "{smiles}");
        }
    }

    #[test]
    // Provenance: fn tests::chosen_anionic_parents_keep_preferred_status | adapted
    // Source: ReShiki | src/naming/openclatura/species.rs::tests::chosen_anionic_parents_keep_preferred_status
    // Source: openclatura | tests/unit/test_chalcogen_functional_groups.py::test_shared_chalcogen_ligand_classifier
    // Source: openclatura | tests/data/parser_xml_resources/arylGroups.json::<fixture data>
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn chosen_anionic_parents_keep_preferred_status() {
        for smiles in [
            "[n-]1nnnc1C(=O)O",
            "[O-]C(=O)c1ccccc1Cc1cc[n-]c1",
            "[O-]C(=O)c1ccccc1[CH2-]",
            "O=C([O-])c1cc[n-]c1",
            "[O-]C(=O)c1ccccc1-c1[n-]cc(C(=O)[O-])c1",
            "C[O-]",
            "[O-]c1ccccc1",
            "Oc1ccccc1[CH2-]",
            "[O-]c1ccccc1Cc1ccccc1C(=O)[O-]",
            "Oc1ccccc1-c1[n-]ccc1C",
            "OCCc1cc[n-]c1",
            "OC(=O)c1ccccc1-c1[n-]ccc1",
        ] {
            assert_eq!(assess_component(smiles), NameStatus::Preferred, "{smiles}");
        }
    }

    #[test]
    // Provenance: fn tests::skeletal_anion_with_characteristic_cation_has_unsettled_suffix_order | original
    // Source: ReShiki | src/naming/openclatura/species.rs::tests::skeletal_anion_with_characteristic_cation_has_unsettled_suffix_order
    fn skeletal_anion_with_characteristic_cation_has_unsettled_suffix_order() {
        for smiles in ["[NH3+]CC[CH-]C#N", "[NH3+][C-]1C=CC2=C1N=NO2"] {
            assert_eq!(
                assess_component(smiles),
                NameStatus::Unverified(NameReason::IonicSuffixOrder),
                "{smiles}"
            );
        }
    }

    #[test]
    // Provenance: fn tests::tied_external_anionic_parents_remain_unverified | original
    // Source: ReShiki | src/naming/openclatura/species.rs::tests::tied_external_anionic_parents_remain_unverified
    fn tied_external_anionic_parents_remain_unverified() {
        assert_eq!(
            assess_component("[n-]1ccc(c1)CCc1cc[n-]c1"),
            NameStatus::Unverified(NameReason::AnionParent),
        );
    }

    #[test]
    // Provenance: fn tests::retained_anion_assessment_requires_application_evidence | adapted
    // Source: ReShiki | src/naming/openclatura/species.rs::tests::retained_anion_assessment_requires_application_evidence
    // Source: openclatura | tests/data/parser_xml_resources/simpleGroups.json::<fixture data>
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn retained_anion_assessment_requires_application_evidence() {
        use super::super::assembly_parts::RetainedAnionName;

        let mol = super::super::graph_io::read_smiles("C[O-]", &Default::default()).unwrap();
        let atoms = mol.atoms.keys().copied().collect();
        let mut component =
            super::super::component_namer::name_component_with_citation(&mol, &atoms, None, None)
                .unwrap();
        let Derivation::Plan {
            retained_anion_name,
            ..
        } = &mut component.derivation
        else {
            panic!("expected an assembled retained anion parent");
        };
        assert_eq!(
            retained_anion_name,
            &RetainedAnionName::Applied("methoxide".into())
        );
        // The same displayed text cannot certify a selected but unapplied name.
        *retained_anion_name = RetainedAnionName::Required {
            name: "methoxide".into(),
            absorbed_substituent: None,
        };
        assess_anion_parent(&mol, &component.derivation);
        assert_eq!(
            mol.name_status(),
            NameStatus::Unverified(NameReason::AnionParent)
        );
    }

    #[test]
    // Provenance: fn tests::sulfonium_ylide_status_follows_substituent_anions | original
    // Source: ReShiki | src/naming/openclatura/species.rs::tests::sulfonium_ylide_status_follows_substituent_anions
    fn sulfonium_ylide_status_follows_substituent_anions() {
        for (smiles, expected, status) in [
            (
                "C[S+](C)[CH-]c1cc([O-])ccc1[O-]",
                "(dimethylsulfaniumyl)2,5-dioxidophenylmethanide",
                NameStatus::Unverified(NameReason::AnionParent),
            ),
            (
                "C[S+](C)[CH-]c1ccc([N+](=O)[O-])cc1",
                "(dimethylsulfaniumyl)4-nitrophenylmethanide",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
        ] {
            let name = super::super::name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected);
            assert_eq!(name.status, status, "{smiles}");
        }
    }

    #[test]
    // Provenance: fn tests::unsupported_anionic_parent_examples_still_decline | original
    // Source: ReShiki | src/naming/openclatura/species.rs::tests::unsupported_anionic_parent_examples_still_decline
    fn unsupported_anionic_parent_examples_still_decline() {
        for smiles in ["[Na+].Oc1ccccc1-[n-]1cccc1", "[Na+].Oc1ccccc1-n1cccc1"] {
            assert!(
                super::super::name_smiles(smiles, &Default::default()).is_err(),
                "{smiles}"
            );
        }
    }

    #[test]
    // Provenance: fn tests::anionic_parents_precede_neutral_groups_and_roundtrip_complete_salts | original
    // Source: ReShiki | src/naming/openclatura/species.rs::tests::anionic_parents_precede_neutral_groups_and_roundtrip_complete_salts
    fn anionic_parents_precede_neutral_groups_and_roundtrip_complete_salts() {
        let parser = nazuchem_n2s::Parser::new().unwrap();
        for (smiles, expected, status) in [
            (
                "[Na+].Oc1ccccc1[C-]1C=CC=C1",
                "sodium 1-(2-hydroxyphenyl)cyclopenta-2,4-dien-1-ide",
                NameStatus::Preferred,
            ),
            (
                "[Na+].OCC[CH2-]",
                "sodium 3-hydroxypropan-1-ide",
                NameStatus::Preferred,
            ),
            (
                "[Na+].Oc1ccccc1-c1[n-]ccc1C",
                "sodium 2-(2-hydroxyphenyl)-3-methyl-1H-pyrrol-1-ide",
                NameStatus::Preferred,
            ),
            (
                "Oc1ccccc1-c1[n-]ccc1C",
                "2-(2-hydroxyphenyl)-3-methyl-1H-pyrrol-1-ide",
                NameStatus::Preferred,
            ),
            ("[n-]1cccc1", "1H-pyrrol-1-ide", NameStatus::Preferred),
            (
                "[Na+].Cc1cc[n-]c1",
                "sodium 3-methyl-1H-pyrrol-1-ide",
                NameStatus::Preferred,
            ),
            (
                "[Na+].Oc1ccccc1Cc1cc[n-]c1",
                "sodium 3-[(2-hydroxyphenyl)methyl]-1H-pyrrol-1-ide",
                NameStatus::Preferred,
            ),
            (
                "[Na+].OC(=O)c1ccccc1-c1[n-]ccc1",
                "sodium 2-(2-carboxyphenyl)-1H-pyrrol-1-ide",
                NameStatus::Preferred,
            ),
            (
                "[Na+].Oc1ccccc1[CH2-]",
                "sodium (2-hydroxyphenyl)methanide",
                NameStatus::Preferred,
            ),
            (
                "[Na+].Oc1ccc(cc1)[CH2-]",
                "sodium (4-hydroxyphenyl)methanide",
                NameStatus::Preferred,
            ),
            (
                "[Na+].OC(=O)c1ccccc1Cc1ccccc1[O-]",
                "sodium 2-[(2-carboxyphenyl)methyl]phenoxide",
                NameStatus::Preferred,
            ),
            (
                "[Na+].OC(=O)c1ccccc1Cc1ccccc1C(=O)[O-]",
                "sodium 2-[(2-carboxyphenyl)methyl]benzoate",
                NameStatus::Preferred,
            ),
            (
                "[Na+].[O-]C(=O)c1ccccc1Cc1ccccc1O",
                "sodium 2-[(2-hydroxyphenyl)methyl]benzoate",
                NameStatus::Preferred,
            ),
            (
                "[Na+].[O-]C(=O)Cc1ccccc1O",
                "sodium (2-hydroxyphenyl)acetate",
                NameStatus::Preferred,
            ),
            (
                "[Na+].Oc1ccccc1Cc1ccccc1[O-]",
                "sodium 2-[(2-hydroxyphenyl)methyl]phenoxide",
                NameStatus::Preferred,
            ),
            (
                "[Na+].NC(=O)c1ccccc1Cc1ccccc1[O-]",
                "sodium 2-[(2-carbamoylphenyl)methyl]phenoxide",
                NameStatus::Preferred,
            ),
            (
                "[Na+].O=Cc1ccccc1Cc1ccccc1[O-]",
                "sodium 2-[(2-formylphenyl)methyl]phenoxide",
                NameStatus::Preferred,
            ),
            (
                "[Na+].CC(=O)c1ccccc1Cc1ccccc1[O-]",
                "sodium 2-[(2-acetylphenyl)methyl]phenoxide",
                NameStatus::Preferred,
            ),
            (
                "[Na+].[Na+].[O-]c1ccccc1Cc1ccccc1C(=O)[O-]",
                "disodium 2-[(2-oxidophenyl)methyl]benzoate",
                NameStatus::Preferred,
            ),
            (
                "[Na+].OCCc1cc[n-]c1",
                "sodium 3-(2-hydroxyethyl)-1H-pyrrol-1-ide",
                NameStatus::Preferred,
            ),
            (
                "[Na+].Oc1ccccc1-c1nc[n-]c1",
                "sodium 4-(2-hydroxyphenyl)-1H-imidazol-1-ide",
                NameStatus::Preferred,
            ),
            (
                "[Na+].[n-]1nnnc1-c1ccc(O)cc1",
                "sodium 5-(4-hydroxyphenyl)-1H-tetrazol-1-ide",
                NameStatus::Preferred,
            ),
            (
                "[Na+].[n-]1nnnc1C(=O)O",
                "sodium 5-carboxy-1H-tetrazol-1-ide",
                NameStatus::Preferred,
            ),
            (
                "[Na+].c1ccc2[n-]cc(C(=O)O)c2c1",
                "sodium 3-carboxy-1H-indol-1-ide",
                NameStatus::Preferred,
            ),
            (
                "[Na+].OC(=O)c1nc[n-]c1",
                "sodium 4-carboxy-1H-imidazol-1-ide",
                NameStatus::Preferred,
            ),
            (
                "[Na+].[n-]1ncc(C(=O)O)c1",
                "sodium 4-carboxy-1H-pyrazol-1-ide",
                NameStatus::Preferred,
            ),
            (
                "[O-]C(=O)c1ccccc1-c1[n-]cc(C(=O)[O-])c1",
                "5-(2-carboxylatophenyl)-1H-pyrrol-1-ide-3-carboxylate",
                NameStatus::Preferred,
            ),
            (
                "[NH3+]CC1[N-]C(=O)NC1=O",
                "5-(azaniumylmethyl)-2,4-dioxoimidazolidin-1-ide",
                NameStatus::Preferred,
            ),
            ("CC(=O)[O-]", "acetate", NameStatus::Preferred),
            (
                "O=C([O-])CCC(=O)O",
                "3-carboxypropanoate",
                NameStatus::Preferred,
            ),
            (
                "[Na+].O=C([O-])CCC(=O)O",
                "sodium 3-carboxypropanoate",
                NameStatus::Preferred,
            ),
        ] {
            let name = super::super::name_smiles(smiles, &Default::default())
                .unwrap_or_else(|error| panic!("{smiles}: {error}"));
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, status, "{smiles}");
            let parsed = parser.parse(&name.text, &nazuchem_n2s::ParseOptions::strict());
            assert_eq!(
                parsed.status,
                nazuchem_n2s::Status::Success,
                "{}: {}",
                name.text,
                parsed.message
            );
            let reconstructed = parsed.structure.unwrap().semantic_cxsmiles().unwrap();
            crate::test_support::verify_identity(smiles, &reconstructed).unwrap_or_else(|error| {
                panic!("{smiles} -> {} -> {reconstructed}: {error}", name.text)
            });
        }
    }

    #[test]
    // Provenance: fn tests::retained_alkoxide_parent_is_preferred | adapted
    // Source: ReShiki | src/naming/openclatura/species.rs::tests::retained_alkoxide_parent_is_preferred
    // Source: openclatura | tests/data/parser_xml_resources/simpleGroups.json::<fixture data>
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn retained_alkoxide_parent_is_preferred() {
        let name = super::super::name_smiles("C[O-]", &Default::default()).unwrap();
        assert_eq!(name.text, "methoxide");
        assert_eq!(name.status, NameStatus::Preferred);
    }

    #[test]
    // Provenance: fn tests::acid_anion_parent_precedes_the_neutral_acid | original
    // Source: ReShiki | src/naming/openclatura/species.rs::tests::acid_anion_parent_precedes_the_neutral_acid
    fn acid_anion_parent_precedes_the_neutral_acid() {
        let name = super::super::name_smiles("O=C([O-])CCC(=O)O", &Default::default()).unwrap();
        assert_eq!(name.text, "3-carboxypropanoate");
        assert_eq!(name.status, NameStatus::Preferred);
    }

    #[test]
    // Provenance: fn tests::salt_inherits_retained_alkoxide_parent_preference | original
    // Source: ReShiki | src/naming/openclatura/species.rs::tests::salt_inherits_retained_alkoxide_parent_preference
    fn salt_inherits_retained_alkoxide_parent_preference() {
        let name = super::super::name_smiles("[Na+].C[O-]", &Default::default()).unwrap();
        assert_eq!(name.text, "sodium methoxide");
        assert_eq!(name.status, NameStatus::Preferred);
    }

    #[test]
    // Provenance: fn tests::salt_inherits_the_preferred_acid_anion_parent | original
    // Source: ReShiki | src/naming/openclatura/species.rs::tests::salt_inherits_the_preferred_acid_anion_parent
    fn salt_inherits_the_preferred_acid_anion_parent() {
        let name =
            super::super::name_smiles("[Na+].O=C([O-])CCC(=O)O", &Default::default()).unwrap();
        assert_eq!(name.text, "sodium 3-carboxypropanoate");
        assert_eq!(name.status, NameStatus::Preferred);
    }

    #[test]
    // Provenance: fn tests::organic_ions_and_unaudited_inorganic_names_keep_their_status | adapted
    // Source: ReShiki | src/naming/openclatura/species.rs::tests::organic_ions_and_unaudited_inorganic_names_keep_their_status
    // Source: openclatura | tests/data/parser_xml_resources/functionalTerms.json::<fixture data>
    // Source: openclatura | tests/unit/test_analysis.py::test_public_api_golden_names
    // Source: openclatura | tests/data/parser_xml_resources/nonCarboxylicAcids.json::<fixture data>
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn organic_ions_and_unaudited_inorganic_names_keep_their_status() {
        for (smiles, expected, status) in [
            ("CC(=O)[O-]", "acetate", NameStatus::Preferred),
            (
                "O=S(=O)([O-])[O-]",
                "sulfate",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "[Cl-]",
                "chloride",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "[N-]=[N+]=[N-]",
                "azide",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "[O-]N=O",
                "nitrite",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            ("C[NH3+]", "methanaminium", NameStatus::Preferred),
            ("CC(=O)O", "acetic acid", NameStatus::Preferred),
        ] {
            let name = super::super::name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, status, "{smiles}");
        }
    }

    #[test]
    // Provenance: fn tests::balanced_salts_roundtrip_as_complete_graphs | adapted
    // Source: ReShiki | src/naming/openclatura/species.rs::tests::balanced_salts_roundtrip_as_complete_graphs
    // Source: openclatura | tests/unit/test_analysis.py::test_component_shortcuts_are_final_binding_audited
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn balanced_salts_roundtrip_as_complete_graphs() {
        let parser = nazuchem_n2s::Parser::new().unwrap();
        for (smiles, expected, status) in [
            ("CC(=O)[O-].[Na+]", "sodium acetate", NameStatus::Preferred),
            (
                "C[NH3+].[Cl-]",
                "methanaminium chloride",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "[Na+].[Na+].O=S(=O)([O-])[O-]",
                "disodium sulfate",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "[Na+].[K+].O=S(=O)([O-])[O-]",
                "potassium sodium sulfate",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "[Ca+2].[Cl-].[Cl-]",
                "calcium dichloride",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "[Ca+2].CC(=O)[O-].CC(=O)[O-]",
                "calcium diacetate",
                NameStatus::Preferred,
            ),
            (
                "C[NH3+].C[NH3+].O=S(=O)([O-])[O-]",
                "bis(methanaminium) sulfate",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "[Fe+2]",
                "iron(2+)",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "[Na+]",
                "sodium(1+)",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "[Ca+2]",
                "calcium(2+)",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
        ] {
            let name = super::super::name_smiles(smiles, &Default::default())
                .unwrap_or_else(|e| panic!("{smiles}: {e}"));
            assert_eq!(name.text, expected);
            assert_eq!(name.status, status, "{smiles}");
            let parsed = parser.parse(&name.text, &nazuchem_n2s::ParseOptions::strict());
            assert_eq!(
                parsed.status,
                nazuchem_n2s::Status::Success,
                "{}: {}",
                name.text,
                parsed.message
            );
            let reconstructed = parsed.structure.unwrap().semantic_cxsmiles().unwrap();
            crate::test_support::verify_identity(smiles, &reconstructed)
                .unwrap_or_else(|e| panic!("{smiles} -> {} -> {reconstructed}: {e}", name.text));
        }
    }

    /// Issue #25. P-14.8.1: adduct names join the components with em dashes
    /// and cite the ratio; organic components are in order of class seniority
    /// (P-41) and an enclosed salt leads. P-14.8.2: inorganic components
    /// follow, water is last, a mixed adduct has no PIN, and a hydrate may be
    /// cited with a numerical prefix. P-77.1.3 (1), (3): a salt of an organic
    /// base drawn with the neutral acid is an adduct or a hydrohalide.
    #[test]
    // Provenance: fn tests::neutral_components_are_named_as_adducts | original
    // Source: Nazuchem | crates/nazuchem-s2n/src/species.rs::tests::neutral_components_are_named_as_adducts
    fn neutral_components_are_named_as_adducts() {
        let parser = nazuchem_n2s::Parser::new().unwrap();
        let mixed = NameStatus::NotPreferred(NameReason::Adduct);
        let organic = NameStatus::Unverified(NameReason::Adduct);
        for (smiles, expected, status, alternative) in [
            // P-14.8.2 examples.
            (
                "CN1CCC[C@H]1c1cccnc1.Cl",
                "3-[(2S)-1-methylpyrrolidin-2-yl]pyridine\u{2014}hydrogen chloride (1/1)",
                mixed,
                Some("3-[(2S)-1-methylpyrrolidin-2-yl]pyridine hydrochloride"),
            ),
            (
                "NCC(=O)Nc1ccc(OCC)cc1.O",
                "2-amino-N-(4-ethoxyphenyl)acetamide\u{2014}water (1/1)",
                mixed,
                Some("2-amino-N-(4-ethoxyphenyl)acetamide monohydrate"),
            ),
            (
                "OC(=O)C(O)=O.NCCN.OC(=O)C(O)=O.NCCN.O.O.O",
                "oxalic acid\u{2014}ethane-1,2-diamine\u{2014}water (2/2/3)",
                mixed,
                None,
            ),
            (
                "NCc1ccc(cc1)S(=O)(=O)c1ccc(N)cc1.Cl.O",
                "4-[4-(aminomethyl)benzene-1-sulfonyl]aniline\u{2014}hydrogen chloride\u{2014}water (1/1/1)",
                mixed,
                None,
            ),
            // P-77.1.3 example (1).
            (
                "CN(C)C1SCCN1.CN(C)C1SCCN1.OS(=O)(=O)O",
                "N,N-dimethyl-1,3-thiazolidin-2-amine\u{2014}sulfuric acid (2/1)",
                mixed,
                None,
            ),
            // P-77.1.3 (3): hydrohalides.
            (
                "Nc1ccccc1.Cl",
                "aniline\u{2014}hydrogen chloride (1/1)",
                mixed,
                Some("aniline hydrochloride"),
            ),
            (
                "CCN.Br",
                "ethanamine\u{2014}hydrogen bromide (1/1)",
                mixed,
                Some("ethanamine hydrobromide"),
            ),
            (
                "C1C=CCC2C1CNC2.I",
                "2,3,3a,4,7,7a-hexahydro-1H-isoindole\u{2014}hydrogen iodide (1/1)",
                mixed,
                Some("2,3,3a,4,7,7a-hexahydro-1H-isoindole hydroiodide"),
            ),
            (
                "NCCN.Cl.Cl",
                "ethane-1,2-diamine\u{2014}hydrogen chloride (1/2)",
                mixed,
                Some("ethane-1,2-diamine dihydrochloride"),
            ),
            // The ratio counts the drawn components and is not reduced.
            (
                "Nc1ccccc1.Nc1ccccc1.Cl.Cl",
                "aniline\u{2014}hydrogen chloride (2/2)",
                mixed,
                None,
            ),
            (
                "[H]Cl.NCC",
                "ethanamine\u{2014}hydrogen chloride (1/1)",
                mixed,
                Some("ethanamine hydrochloride"),
            ),
            // A hydrocarbon is not a base, so it takes no hydrohalide form.
            (
                "c1ccccc1.Cl",
                "benzene\u{2014}hydrogen chloride (1/1)",
                mixed,
                None,
            ),
            // P-14.8.2: organic components precede inorganic ones.
            (
                "CCO.[Na+].[Cl-]",
                "ethanol\u{2014}(sodium chloride) (1/1)",
                mixed,
                None,
            ),
            // P-14.8.1: a salt is enclosed; P-14.8.2: hydrates, hemi included.
            (
                "CC(=O)[O-].[Na+].O.O.O",
                "(sodium acetate)\u{2014}water (1/3)",
                mixed,
                Some("sodium acetate trihydrate"),
            ),
            (
                "C[NH3+].C[NH3+].O=S(=O)([O-])[O-].O.O",
                "[bis(methanaminium) sulfate]\u{2014}water (1/2)",
                mixed,
                Some("bis(methanaminium) sulfate dihydrate"),
            ),
            (
                "Nc1ccccc1.Nc1ccccc1.O",
                "aniline\u{2014}water (2/1)",
                mixed,
                Some("aniline hemihydrate"),
            ),
            // P-14.8.2: labelled water is still last; no hydrate form is given.
            (
                "CN.Cl.[18OH2]",
                "methanamine\u{2014}hydrogen chloride\u{2014}(18O)water (1/1/1)",
                mixed,
                None,
            ),
            // P-14.8.1 examples: organic adducts in order of class seniority.
            (
                "CCO.c1ccncc1",
                "ethanol\u{2014}pyridine (1/1)",
                organic,
                None,
            ),
            (
                "Oc1cccc2cccnc12.OC(=O)c1ccc(cc1)[N+]([O-])=O",
                "4-nitrobenzoic acid\u{2014}quinolin-8-ol (1/1)",
                organic,
                None,
            ),
            // P-77.1.3 (1): organic acids are senior to amines (P-41).
            (
                "Nc1ccccc1.OC=O",
                "formic acid\u{2014}aniline (1/1)",
                organic,
                None,
            ),
            (
                "Nc1ccccc1.Cc1ccc(S(=O)(=O)O)cc1",
                "4-methylbenzene-1-sulfonic acid\u{2014}aniline (1/1)",
                organic,
                None,
            ),
            (
                "Nc1ccccc1.OC(=O)C(F)(F)F",
                "trifluoroacetic acid\u{2014}aniline (1/1)",
                organic,
                None,
            ),
        ] {
            let name = super::super::name_smiles(smiles, &Default::default())
                .unwrap_or_else(|e| panic!("{smiles}: {e}"));
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, status, "{smiles}");
            let alternatives: Vec<_> = name.alternatives.iter().map(|a| &*a.name).collect();
            assert_eq!(alternatives, Vec::from_iter(alternative), "{smiles}");
            assert!(name.alternatives.iter().all(|a| a.status == mixed));
            for text in std::iter::once(expected).chain(alternative) {
                let parsed = parser.parse(text, &nazuchem_n2s::ParseOptions::strict());
                assert_eq!(
                    parsed.status,
                    nazuchem_n2s::Status::Success,
                    "{text}: {}",
                    parsed.message
                );
                let reconstructed = parsed.structure.unwrap().semantic_cxsmiles().unwrap();
                crate::test_support::verify_identity(smiles, &reconstructed)
                    .unwrap_or_else(|e| panic!("{smiles} -> {text} -> {reconstructed}: {e}"));
            }
        }
        // Identical molecules are no adduct, and the Red Book order of several
        // inorganic components is not covered.
        for smiles in [
            "O.O",
            "Nc1ccccc1.Cl.Br",
            "[NH3+]CC(=O)[O-].[Na+].[Cl-]",
            "[NH3+]CC(=O)[O-].CC(=O)O",
            "[CH3].Cl",
            // Same class: the further P-41 and P-44 criteria are not covered.
            "ClCCl.c1ccncc1",
            "Oc1cccc(O)c1O.Oc1cccc2cccnc12",
        ] {
            assert!(super::super::name_smiles(smiles, &Default::default()).is_err());
        }
    }

    /// A salt takes its strongest ion assessment (Preferred < Unverified <
    /// NotPreferred; ties keep the first ion's reason in citation order), so a
    /// known departure is never masked by an unchecked ion (P-77.1.1, P-77.2.1).
    /// Retained anion spellings and salt order must preserve stronger departures.
    #[test]
    // Provenance: fn tests::salt_status_is_the_strongest_ion_downgrade | adapted
    // Source: ReShiki | src/naming/openclatura/species.rs::tests::salt_status_is_the_strongest_ion_downgrade
    // Source: openclatura | tests/unit/test_analysis.py::test_public_api_golden_names
    // Source: openclatura | tests/data/parser_xml_resources/arylGroups.json::<fixture data>
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn salt_status_is_the_strongest_ion_downgrade() {
        use NameReason::{BridgedFusedParent, Multiplicative, RingAssembly, UnauditedParent};
        use NameStatus::{NotPreferred, Preferred, Unverified};
        const NTA: &str = "O=C([O-])CN(CC(=O)[O-])CC(=O)[O-]";
        const OXYDIBENZOATE: &str = "O=C([O-])c1ccc(cc1)Oc1ccc(cc1)C(=O)[O-]";
        const EPOXY_AMINIUM: &str = "C[N+](C)(C)c1ccc2c(c1)C1CCC2O1";
        const NONANYL_AMINIUM: &str = "C[N+](C)(C)CC1C2CC3CC1CC2C3";
        const BIS_AMINIUM: &str = "C[N+](C)(C)Cc1ccc(cc1)C[N+](C)(C)C";
        const PENTALENE_CARBOXYLATE: &str = "C1C2CC3CC1CC2C3C(=O)[O-]";
        for (smiles, expected, status) in [
            (
                format!("[Na+].[Na+].[Na+].[Na+].[O-]c1ccccc1.{NTA}"),
                "tetrasodium [bis(carboxylatomethyl)amino]acetate phenoxide",
                NotPreferred(Multiplicative),
            ),
            (
                format!("{NTA}.[O-]c1ccccc1.[Na+].[Na+].[Na+].[Na+]"),
                "tetrasodium [bis(carboxylatomethyl)amino]acetate phenoxide",
                NotPreferred(Multiplicative),
            ),
            (
                format!("{EPOXY_AMINIUM}.[Na+].[Na+].[Na+].C1=CCCCC1-C1=CCCC(C1)C(=O)[O-].{NTA}"),
                "trisodium N,N,N-trimethyl-1,2,3,4-tetrahydro-1,4-epoxynaphthalen-6-aminium [bis(carboxylatomethyl)amino]acetate 3-(cyclohex-2-en-1-yl)cyclohex-3-ene-1-carboxylate",
                NotPreferred(Multiplicative),
            ),
            (
                "C[N+](C)(C)C1=CCCCC1-C1=CCCCC1.[Na+].[O-]c1ccccc1.[O-]C(=O)c1ccccc1".into(),
                "6-(cyclohex-1-en-1-yl)-N,N,N-trimethylcyclohex-1-en-1-aminium sodium benzoate phenoxide",
                NotPreferred(RingAssembly),
            ),
            (
                format!("[Na+].[Na+].[Na+].CC(=O)CC[O-].{OXYDIBENZOATE}"),
                "trisodium 4-(4-carboxylatophenoxy)benzoate 3-oxobutan-1-olate",
                NotPreferred(Multiplicative),
            ),
            (
                format!("[K+].[S-]c1ccccc1.[K+].[K+].[K+].{NTA}"),
                "tetrapotassium benzenethiolate [bis(carboxylatomethyl)amino]acetate",
                NotPreferred(Multiplicative),
            ),
            (
                format!("{EPOXY_AMINIUM}.{EPOXY_AMINIUM}.{OXYDIBENZOATE}"),
                "bis(N,N,N-trimethyl-1,2,3,4-tetrahydro-1,4-epoxynaphthalen-6-aminium) 4-(4-carboxylatophenoxy)benzoate",
                NotPreferred(Multiplicative),
            ),
            (
                "[Na+].[O-]C(=O)c1ccccc1".into(),
                "sodium benzoate",
                Preferred,
            ),
            ("[Na+].[O-]c1ccccc1".into(), "sodium phenoxide", Preferred),
            (
                "[Na+].[O-]c1ccccc1.[Na+].[O-]C(=O)c1ccccc1".into(),
                "disodium benzoate phenoxide",
                Preferred,
            ),
            (
                format!("[Na+].[Na+].{OXYDIBENZOATE}"),
                "disodium 4-(4-carboxylatophenoxy)benzoate",
                NotPreferred(Multiplicative),
            ),
            (
                format!("[O-]c1ccccc1.[O-]c1ccccc1.{BIS_AMINIUM}"),
                "N,N,N-trimethyl-1-{4-[(trimethylammonio)methyl]phenyl}methanaminium diphenoxide",
                NotPreferred(Multiplicative),
            ),
            (
                format!("{PENTALENE_CARBOXYLATE}.{BIS_AMINIUM}.{PENTALENE_CARBOXYLATE}"),
                "N,N,N-trimethyl-1-{4-[(trimethylammonio)methyl]phenyl}methanaminium bis(octahydro-2,5-methanopentalene-1-carboxylate)",
                NotPreferred(Multiplicative),
            ),
            (
                format!("{NONANYL_AMINIUM}.[O-]C(=O)c1ccccc1"),
                "N,N,N-trimethyl-1-(tricyclo[3.3.1.0^{3,7}]nonan-2-yl)methanaminium benzoate",
                Unverified(UnauditedParent),
            ),
            (
                format!("{EPOXY_AMINIUM}.[O-]c1ccccc1"),
                "N,N,N-trimethyl-1,2,3,4-tetrahydro-1,4-epoxynaphthalen-6-aminium phenoxide",
                Unverified(BridgedFusedParent),
            ),
            (
                format!("{NONANYL_AMINIUM}.{NONANYL_AMINIUM}.{OXYDIBENZOATE}"),
                "bis[N,N,N-trimethyl-1-(tricyclo[3.3.1.0^{3,7}]nonan-2-yl)methanaminium] 4-(4-carboxylatophenoxy)benzoate",
                NotPreferred(Multiplicative),
            ),
            ("[O-]c1ccccc1".into(), "phenoxide", Preferred),
            (
                NTA.into(),
                "[bis(carboxylatomethyl)amino]acetate",
                NotPreferred(Multiplicative),
            ),
            ("CC(=O)O".into(), "acetic acid", Preferred),
        ] {
            let name = super::super::name_smiles(&smiles, &Default::default())
                .unwrap_or_else(|e| panic!("{smiles}: {e}"));
            assert_eq!(
                (name.text.as_str(), name.status),
                (expected, status),
                "{smiles}"
            );
        }
    }

    #[test]
    // Provenance: fn tests::salt_inherits_the_strongest_ion_status | original
    // Source: ReShiki | src/naming/openclatura/species.rs::tests::salt_inherits_the_strongest_ion_status
    fn salt_inherits_the_strongest_ion_status() {
        let name =
            super::super::name_smiles("[Na+].[Na+].[Cl-].C[O-]", &Default::default()).unwrap();
        assert_eq!(name.text, "disodium chloride methoxide");
        assert_eq!(
            name.status,
            NameStatus::Unverified(NameReason::UnauditedPath)
        );
    }

    #[test]
    // Provenance: fn tests::monatomic_cations_are_certified_only_in_salts | adapted
    // Source: ReShiki | src/naming/openclatura/species.rs::tests::monatomic_cations_are_certified_only_in_salts
    // Source: openclatura | tests/unit/test_analysis.py::test_component_shortcuts_are_final_binding_audited
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn monatomic_cations_are_certified_only_in_salts() {
        for smiles in ["[Na+]", "[Fe+2]", "[Ca+2]"] {
            assert_eq!(
                monatomic_cation(smiles, false, &Default::default())
                    .unwrap()
                    .unwrap()
                    .status,
                NameStatus::Unverified(NameReason::UnauditedPath),
                "{smiles}"
            );
            assert_eq!(
                monatomic_cation(smiles, true, &Default::default())
                    .unwrap()
                    .unwrap()
                    .status,
                NameStatus::Preferred,
                "{smiles}"
            );
        }
    }

    #[test]
    // Provenance: fn tests::mixtures_unbalanced_sets_and_zwitterions_are_refused | original
    // Source: ReShiki | src/naming/openclatura/species.rs::tests::mixtures_unbalanced_sets_and_zwitterions_are_refused
    fn mixtures_unbalanced_sets_and_zwitterions_are_refused() {
        for smiles in ["C.C", "[Na+].[Na+].[Cl-]", "[NH3+]CC(=O)[O-].[Na+].[Cl-]"] {
            assert!(
                super::super::name_smiles(smiles, &Default::default()).is_err(),
                "{smiles}"
            );
        }
    }
}
