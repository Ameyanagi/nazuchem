// SPDX-License-Identifier: MIT AND (MIT OR Apache-2.0)
// Origin: mixed
// Source: Nazuchem | crates/nazuchem-s2n/src/species.rs | issue #25
// Source: ReShiki | src/naming/openclatura/species.rs | 5f559fd482e237307a589dfcb7ea333cee0a6e13
// Source: Nazuchem | crates/nazuchem-s2n/src/species.rs | issue #32
// Source: Nazuchem | crates/nazuchem-s2n/src/species.rs | issue #11
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
use super::{
    Name, component_namer::Derivation, molecule::Molecule, name_smiles_with_style,
    special_cases::alphanumerical_order_key,
};
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
    let centres = cationic_centres(mol);
    let cited = |i: &usize| {
        parent_atoms.contains(i)
            || principal_group
                .as_ref()
                .is_some_and(|(_, atoms)| atoms.contains(i))
    };
    if !centres.is_empty() && !centres.iter().any(cited) {
        mol.mark_name(NameStatus::Unverified(NameReason::CationParent));
    }
}

/// P-70.2, P-73.7(a): the cationic centres of a cation, leaving out
/// charge-separated neutral groups. Empty for neutral species and for anions
/// and zwitterions, which take the anionic parent (P-72, P-74).
// Provenance: fn cationic_centres | original
// Source: Nazuchem | crates/nazuchem-s2n/src/species.rs::cationic_centres
pub(super) fn cationic_centres(mol: &Molecule) -> std::collections::BTreeSet<usize> {
    let neutral_group_atoms = neutral_group_atoms(mol);
    let centres: std::collections::BTreeSet<usize> = mol
        .atoms
        .iter()
        .filter(|(i, atom)| atom.charge != 0 && !neutral_group_atoms.contains(i))
        .map(|(i, _)| *i)
        .collect();
    if centres.iter().any(|i| mol.atoms[i].charge < 0) {
        return Default::default();
    }
    centres
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
            ) || group.descriptor.as_ref().is_some_and(|descriptor| {
                descriptor.derivative == "zwitterion"
                    // A sulfonyl azide (P-65.5.2.1).
                    || descriptor.family == "central_acid"
                        && descriptor
                            .leaving_group
                            .as_ref()
                            .is_some_and(|l| l.kind == "N3")
            })) && group
                .atoms_involved
                .iter()
                .map(|i| mol.atoms[i].charge)
                .sum::<i32>()
                == 0
        })
        .flat_map(|group| group.atoms_involved)
        .collect()
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
        let labelled = state.graph.atoms.iter().any(|a| a.isotope != 0)
            || state.metadata.atoms.iter().any(|a| a.chiral_tag != 0)
            || state.metadata.bonds.iter().any(|b| b.stereo != 0);
        // A mononuclear oxoanion with two or more oxido or hydroxy groups on
        // its central atom, such as carbonate, sulfate or phosphate. Its acid
        // condenses to a polynuclear one (dicarbonic acid, P-65.2.3;
        // disulfuric acid, P-67.2.1), so `dicarbonate` would name that anion;
        // nitrite has a single such group and keeps `di`.
        let (atoms, bonds) = (&state.graph.atoms, &state.graph.bonds);
        let central: Vec<_> = (0..atoms.len())
            .filter(|&i| atoms[i].atomic_number != 8)
            .collect();
        let oxoanion = q < 0
            && matches!(central[..], [c] if atoms[c].explicit_hydrogens == 0
                && state.valences[c].implicit_hydrogens == 0
                && bonds
                    .iter()
                    .filter(|b| {
                        let o = if b.a == c { b.b } else { b.a };
                        b.order == 1
                            && (b.a == c || b.b == c)
                            && bonds.iter().filter(|x| x.a == o || x.b == o).count() == 1
                    })
                    .count()
                    >= 2);
        // P-16.3.6 (a): a monatomic anion is multiplied with `bis` when the
        // `di` form names a homopolyatomic anion: dioxide (O2 2-) and
        // disulfide (S2 2-) for dianions, triiodide (I3-) for halides.
        let monatomic = q < 0
            && atoms.len() == 1
            && atoms[0].explicit_hydrogens == 0
            && state.valences[0].implicit_hydrogens == 0;
        // P-16.3.6 (a): two NC-SO2-O(-) or NC-S(O)-O(-) ions are
        // `bis(sulfurocyanidate)`; `disulfuro...` would be read as a
        // derivative of disulfuric acid (P-67.2.1), and `diphosphoro...` as
        // one of diphosphoric acid.
        let mut elements: Vec<_> = central.iter().map(|&i| atoms[i].atomic_number).collect();
        elements.sort_unstable();
        let oxoanion = oxoanion
            || q < 0
                && matches!(elements[..], [6, 7, 15 | 16])
                && atoms.len() >= 5
                && bonds.iter().any(|b| b.order == 3);
        let bis_from = if oxoanion || monatomic && q <= -2 {
            2
        } else if monatomic {
            3
        } else {
            usize::MAX
        };
        ions.push((q < 0, name, component, labelled, bis_from));
    }
    // P-14.5: alphanumerical order leaves out stereo and isotope descriptors,
    // so a labelled ion is keyed by its cited name rendered without them: the
    // same parent, prefixes and locants. Removing the labels from the
    // structure instead can change the parent, as for a gonane ion.
    let keys: Vec<_> = ions
        .iter()
        .map(|(negative, name, component, labelled, _)| {
            if !labelled {
                return Some((alphanumerical_order_key(&name.text), name.text.clone()));
            }
            // Nothing to order against.
            if !ions
                .iter()
                .any(|o| o.0 == *negative && o.1.text != name.text)
            {
                return None;
            }
            let bare = super::formatting::without_descriptors(|| {
                name_smiles_with_style(component, style, budget)
            })
            .ok()?
            .text;
            // The descriptors are the only text left out, so the Roman
            // letters are those of the cited name: prefixes that differ only
            // in their descriptors are multiplied once these are gone.
            let key = alphanumerical_order_key(&bare);
            (key.0 == alphanumerical_order_key(&name.text).0).then_some((key, bare))
        })
        .collect();
    // Ions of one sign are ordered by the key when it tells every two
    // different names apart; otherwise main's character order stays and is
    // unverified. Nonitalic Roman letters decide first (P-14.5), then italic
    // letters (P-14.5.3), then the locants of the primary substituents, where
    // no locant precedes any locant (P-14.5, P-14.5.4), then the internal
    // locants of names that differ in nothing but locants (P-14.5.4).
    // Descriptor brackets hold no such locants, and the order of enclosing
    // marks (P-14.6) is not ranked here. Nor are composite locants with a
    // superscript, as 1^2 (P-14.3.5): the key reads only the bonding number
    // of lambda^5, which takes no part.
    let superscript = |text: &str| text.replace("lambda^", "").contains('^');
    let locant_free = |text: &str| {
        let descriptor_bracket = text
            .split('[')
            .skip(1)
            .any(|rest| !rest.split(']').next().unwrap_or(rest).contains('('));
        (!descriptor_bracket).then(|| text.replace(|c: char| c.is_ascii_digit(), ""))
    };
    let ordered = ions.iter().zip(&keys).all(|((negative, name, ..), key)| {
        ions.iter().zip(&keys).all(|(other, other_key)| {
            other.0 != *negative
                || other.1.text == name.text
                || match (key, other_key) {
                    (Some((a, a_text)), Some((b, b_text))) => {
                        (&a.0, &a.1) != (&b.0, &b.1)
                            || !superscript(a_text)
                                && !superscript(b_text)
                                && (a.2 != b.2
                                    || a.3 != b.3
                                        && locant_free(a_text)
                                            .is_some_and(|x| Some(x) == locant_free(b_text)))
                    }
                    _ => false,
                }
        })
    });
    let mut grouped = BTreeMap::new();
    for ((negative, name, component, _, bis_from), key) in ions.into_iter().zip(keys) {
        let plain: String = name
            .text
            .chars()
            .filter(char::is_ascii_alphanumeric)
            .map(|c| c.to_ascii_lowercase())
            .collect();
        let entry = grouped
            .entry((
                negative,
                key.filter(|_| ordered).map(|(key, _)| key),
                plain,
                name.text.clone(),
                component,
            ))
            .or_insert((name, bis_from, 0));
        entry.2 += 1;
    }
    if charge != 0 {
        return Err(
            "Salt naming requires a charge-balanced set of ions; the total charge must be zero"
                .into(),
        );
    }
    let mut status = NameStatus::Preferred;
    let mut words = Vec::new();
    for (_, (name, bis_from, count)) in grouped {
        // The strongest ion downgrade wins; ties keep the first ion's reason.
        if name.status.severity() > status.severity() {
            status = name.status;
        }
        // P-16.3.6 (a): `bis(carbonate)` for two carbonate ions, since
        // `dicarbonate` (like `disulfate` and `diphosphate`) names the
        // condensed dinuclear anion.
        let complex = count >= bis_from
            || name.text.chars().any(|c| !c.is_ascii_alphabetic())
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
    const UNORDERED: u16 = u16::MAX;
    struct Part {
        // Formula order: kind, then class seniority (P-41).
        rank: (u8, u16),
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
            rank: (if organic_ion { ORGANIC } else { INORGANIC }, 0),
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
            (WATER, 0)
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
            // P-41, Table 4.1: the class of the principal characteristic
            // group, for the classes whose perceived key is that class alone.
            let class = match key.as_deref() {
                Some("carboxylic_acid" | "ring_carboxylic_acid") => 70,
                Some("sulfonic_acid") => 71,
                Some("amide" | "ring_amide") => 110,
                Some("nitrile" | "ring_nitrile") => 140,
                Some("aldehyde" | "ring_aldehyde") => 150,
                Some("ketone") => 160,
                Some("alcohol") => 170,
                Some("amine") => 190,
                Some("imine") => 200,
                // Classes 21 to 43 have no suffix.
                None => 210,
                Some(_) => UNORDERED,
            };
            (ORGANIC, class)
        } else {
            (INORGANIC, 0)
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
    parts.sort_by_key(|a| a.rank);
    if parts.len() < 2 {
        return Err("Adduct naming requires at least two different components".into());
    }
    // Organic components of one class, or without a characteristic group, are
    // ordered by the further criteria of P-41 and P-44, which are not covered.
    let neutral_organic = parts.iter().filter(|p| p.rank.0 == ORGANIC && !p.salt);
    if neutral_organic.clone().count() > 1 && neutral_organic.clone().any(|p| p.rank.1 == UNORDERED)
        || parts
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
    // Provenance: fn tests::descriptor_free_rendering_keeps_the_cited_name | original
    // Source: Nazuchem | crates/nazuchem-s2n/src/species.rs::tests::descriptor_free_rendering_keeps_the_cited_name
    fn descriptor_free_rendering_keeps_the_cited_name() {
        // P-14.5: alphanumerical order leaves out stereochemical and isotopic
        // descriptors. One case per descriptor writer: the parent, prefixes
        // and locants are those of the cited name.
        for (smiles, cited, bare) in [
            (
                "C[C@H](O)C(=O)[O-]",
                "(2S)-2-hydroxypropanoate",
                "2-hydroxypropanoate",
            ),
            (
                "C[C@H]1CC[C@@H](C)CC1",
                "(1s,4s)-1,4-dimethylcyclohexane",
                "1,4-dimethylcyclohexane",
            ),
            // Without its stereo the structure is named from another parent.
            (
                "C[C@]12CC[C@H]3[C@@H](CCC4CC(CC[C@]34C)[NH3+])[C@@H]1CCC2",
                "(8S,9S,10S,13S,14S)-10,13-dimethylgonan-3-aminium",
                "10,13-dimethylgonan-3-aminium",
            ),
            (
                "C/C=C/SS/C=C/C",
                "1,2-bis[(1E)-prop-1-en-1-yl]disulfane",
                "1,2-bis(prop-1-en-1-yl)disulfane",
            ),
            (
                "C[S@](=O)CC",
                "(S)-ethyl methyl sulfoxide",
                "ethyl methyl sulfoxide",
            ),
            (
                "C[S@](=O)c1ccc(C(=O)O)cc1",
                "4-[(S)-methanesulfinyl]benzoic acid",
                "4-(methanesulfinyl)benzoic acid",
            ),
            (
                "OC(=O)c1ccc(/N=N/c2ccccc2)cc1",
                "4-[(E)-phenyldiazenyl]benzoic acid",
                "4-(phenyldiazenyl)benzoic acid",
            ),
            (
                "OC(=O)c1ccc(N/C(N)=N/C)cc1",
                "4-[(E)-2-methylguanidino]benzoic acid",
                "4-(2-methylguanidino)benzoic acid",
            ),
            (
                "OC(=O)c1ccc(cc1)/C(N)=N/C",
                "4-[(Z)-N'-methylcarbamimidoyl]benzoic acid",
                "4-(N'-methylcarbamimidoyl)benzoic acid",
            ),
            (
                "C/C(CC)=N/N",
                "(2Z)-butan-2-one hydrazone",
                "butan-2-one hydrazone",
            ),
            (
                "C/C(/CC)=N\\N=C(/C)\\CC",
                "(1E,2E)-di(butan-2-ylidene)hydrazine",
                "di(butan-2-ylidene)hydrazine",
            ),
            ("[2H]C([2H])([2H])C(=O)[O-]", "(2H3)acetate", "acetate"),
            (
                "[2H]C([2H])([2H])[C@H](O)C(=O)[O-]",
                "(2S)-2-hydroxy(3,3,3-2H3)propanoate",
                "2-hydroxypropanoate",
            ),
        ] {
            let name =
                |smiles| name_smiles_with_style(smiles, IsotopeStyle::Iupac, &Default::default());
            assert_eq!(name(smiles).unwrap().text, cited, "{smiles}");
            let rendered = crate::formatting::without_descriptors(|| name(smiles));
            assert_eq!(rendered.unwrap().text, bare, "{smiles}");
            // The mode ends with the render.
            assert_eq!(name(smiles).unwrap().text, cited, "{smiles}");
        }
    }

    #[test]
    // Provenance: fn tests::salt_ions_follow_alphanumerical_order | original
    // Source: Nazuchem | crates/nazuchem-s2n/src/species.rs::tests::salt_ions_follow_alphanumerical_order
    fn salt_ions_follow_alphanumerical_order() {
        for (smiles, expected) in [
            // P-14.5: nonitalic Roman letters first; the stereodescriptors of
            // the gonane ion take no part, and its cited parent stays gonane.
            (
                "[Cl-].[Cl-].C[C@]12CC[C@H]3[C@@H](CCC4CC(CC[C@]34C)[NH3+])[C@@H]1CCC2.CC[NH3+]",
                "(8S,9S,10S,13S,14S)-10,13-dimethylgonan-3-aminium ethanaminium dichloride",
            ),
            (
                "[Cl-].[Cl-].C[C@]12CC[C@H]3[C@@H](CCC4CC(CC[C@]34C)[NH3+])[C@@H]1CCC2.CCCC[NH3+]",
                "butan-1-aminium (8S,9S,10S,13S,14S)-10,13-dimethylgonan-3-aminium dichloride",
            ),
            // P-14.5: locants of the primary substituents; no locant precedes
            // a locant, so the single prefix at 2 precedes the two at 2,2.
            (
                "[Na+].[Na+].O=C([O-])C(CCCCl)CCCCCC.O=C([O-])C(Cl)(CCC)CCCCCC",
                "disodium 2-(3-chloropropyl)octanoate 2-chloro-2-propyloctanoate",
            ),
            // P-14.5.4: equal primary locants, then the lowest internal
            // locant at the first point of difference (1 before 2).
            (
                "[Na+].[Na+].[O-]C(=O)c1ccc(cc1)[C@H](C)Cl.[O-]C(=O)c1ccc(cc1)CCCl",
                "disodium 4-[(1S)-1-chloroethyl]benzoate 4-(2-chloroethyl)benzoate",
            ),
            // P-14.5.3: italic letters are considered once the Roman letters
            // tie, ahead of the locants of P-14.5.4.
            (
                "CC(=O)[O-].CC(=O)[O-].CC(C)(C)[NH2+]c1ccccc1.[NH3+]c1ccccc1CCCC",
                "2-butylanilinium N-(tert-butyl)anilinium diacetate",
            ),
            // P-14.5.4: the numerals of a von Baeyer descriptor are no locants.
            (
                "[Na+].[Na+].[O-]C(=O)C12CCC(CC1)C2.[O-]C(=O)C1CC2CCC1C2",
                "disodium bicyclo[2.2.1]heptane-1-carboxylate bicyclo[2.2.1]heptane-2-carboxylate",
            ),
            // P-14.5: stereodescriptors take no part; the locants decide.
            (
                "[Na+].[Na+].CC[C@@H](C)[C@@H](O)C([O-])=O.CC[C@H](O)[C@H](C)C([O-])=O",
                "disodium (2R,3R)-2-hydroxy-3-methylpentanoate (2S,3S)-3-hydroxy-2-methylpentanoate",
            ),
        ] {
            let name = name_salt(smiles, IsotopeStyle::Iupac, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_ne!(
                name.status,
                NameStatus::Unverified(NameReason::SaltIonOrder),
                "{smiles}"
            );
        }
        // P-14.5: isotopic descriptors take no part either, so the labelled
        // acetate keeps its place before butanoate.
        let name = name_salt(
            "[Na+].[Na+].[2H]C([2H])([2H])C(=O)[O-].CCCC(=O)[O-]",
            IsotopeStyle::Iupac,
            &Default::default(),
        )
        .unwrap();
        assert_eq!(name.text, "disodium (2H3)acetate butanoate");
        assert_ne!(
            name.status,
            NameStatus::Unverified(NameReason::SaltIonOrder)
        );
    }

    #[test]
    // Provenance: fn tests::salt_ion_order_falls_back_for_uncovered_ions | original
    // Source: ReShiki | src/naming/openclatura/species.rs::tests::salt_ion_order_falls_back_for_uncovered_ions
    fn salt_ion_order_falls_back_for_uncovered_ions() {
        // Identical Roman letters (P-14.5) that the key does not rank.
        for (smiles, expected) in [
            // Ions that differ only in stereodescriptors, which P-14.5 leaves out.
            (
                "[Na+].[Na+].C/C=C/C(=O)[O-].C/C=C\\C(=O)[O-]",
                "disodium (2E)-but-2-enoate (2Z)-but-2-enoate",
            ),
            // Composite locants with superscripts are not ranked (P-14.3.5).
            (
                "[Na+].[Na+].Clc1ccccc1-c1c([O-])c(C)c(-c2ccccc2)cc1.[O-]c1c(Cl)c(C)ccc1-c1ccc(-c2ccccc2)cc1",
                "disodium 1^2-chloro-2^3-methyl[1^1,2^1:2^4,3^1-terphenyl]-2^2-olate 1^3-chloro-1^4-methyl[1^1,2^1:2^4,3^1-terphenyl]-1^2-olate",
            ),
            // Without descriptors the two ligands of the first ion are
            // rendered as one multiplied prefix, which is not the cited name.
            (
                "[Na+].[Na+].O=C([O-])CN([C@H](Br)CN=C=O)[C@@H](Br)CN=C=O.BrCC(=O)[O-]",
                "disodium {[(1R)-1-bromo-2-isocyanatoethyl][(1S)-1-bromo-2-isocyanatoethyl]amino}acetate bromoacetate",
            ),
            // Only the numerals of the von Baeyer descriptors differ.
            (
                "[Na+].[Na+].[O-]C(=O)C1CC2CCC1CC2.[O-]C(=O)C1CCC2CCC1C2",
                "disodium bicyclo[2.2.2]octane-2-carboxylate bicyclo[3.2.1]octane-2-carboxylate",
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
            // P-14.5: stereodescriptors take no part; the ions are compared as
            // 2-hydroxy-3-methylpentanoate and 3-hydroxy-2-methylpentanoate,
            // identical Roman letters, then locants 2,3 before 3,2 (P-14.5.4).
            (
                "[Na+].[Na+].CC[C@@H](C)[C@@H](O)C([O-])=O.CC[C@H](O)[C@H](C)C([O-])=O",
                "disodium (2R,3R)-2-hydroxy-3-methylpentanoate (2S,3S)-3-hydroxy-2-methylpentanoate",
                NameStatus::Preferred,
            ),
            // P-14.5: acetate before hydroxypropanoate, butenoate before
            // propanoate and ethenyl before propenyl, whatever the descriptors.
            (
                "[Na+].[Na+].C[C@H](O)C(=O)[O-].CC(=O)[O-]",
                "disodium acetate (2S)-2-hydroxypropanoate",
                NameStatus::Preferred,
            ),
            (
                "[Na+].[Na+].C/C=C/C(=O)[O-].CCC(=O)[O-]",
                "disodium (2E)-but-2-enoate propanoate",
                NameStatus::Preferred,
            ),
            (
                "[Na+].[Na+].[O-]C(=O)c1ccc(cc1)/C=C/C.[O-]C(=O)c1ccc(cc1)C=C",
                "disodium 4-ethenylbenzoate 4-[(1E)-prop-1-en-1-yl]benzoate",
                NameStatus::Preferred,
            ),
            (
                "[Na+].[Na+].C[C@@H]1CC[C@H](CC1)C(=O)[O-].CC(=O)[O-]",
                "disodium acetate (1r,4r)-4-methylcyclohexane-1-carboxylate",
                NameStatus::Preferred,
            ),
            // P-14.5.4: 2-chloropropanoate before 3-chloropropanoate.
            (
                "[Na+].[Na+].C[C@H](Cl)C(=O)[O-].ClCCC(=O)[O-]",
                "disodium (2S)-2-chloropropanoate 3-chloropropanoate",
                NameStatus::Preferred,
            ),
            // P-14.5.4 with compound locants in the descriptors: 1 before 2.
            (
                "[Cl-].[Cl-].[NH3+][C@@H]1CCC[C@H]2CCCC[C@@H]12.[NH3+][C@@H]1CC[C@@H]2CCCC[C@H]2C1",
                "(1R,4aR,8aR)-decahydronaphthalen-1-aminium (2R,4aS,8aS)-decahydronaphthalen-2-aminium dichloride",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            // P-14.5: isotopic descriptors take no part; acetate before
            // chloroacetate and before hydroxypropanoate.
            (
                "[Na+].[Na+].C[13C](=O)[O-].[O-]C(=O)CCl",
                "disodium (1-13C)acetate chloroacetate",
                NameStatus::Unverified(NameReason::IsotopePlacement),
            ),
            (
                "[Na+].[Na+].[2H]C([2H])([2H])[C@H](O)C(=O)[O-].CC(=O)[O-]",
                "disodium acetate (2S)-2-hydroxy(3,3,3-2H3)propanoate",
                NameStatus::Unverified(NameReason::IsotopePlacement),
            ),
            // P-14.5.3: the italic tert is not a Roman letter (butoxide before
            // methoxide) and breaks a tie after the plain name.
            (
                "[Na+].CC(C)(C)[O-].[Na+].C[O-]",
                "disodium tert-butoxide methoxide",
                NameStatus::Preferred,
            ),
            (
                "[Na+].[Na+].CC(C)(C)[O-].CCCC[O-]",
                "disodium butoxide tert-butoxide",
                NameStatus::Preferred,
            ),
            (
                "[Na+].[Na+].CC(C)(C)c1ccc(cc1)C(=O)[O-].CCCCc1ccc(cc1)C(=O)[O-]",
                "disodium 4-butylbenzoate 4-(tert-butyl)benzoate",
                NameStatus::Preferred,
            ),
            // P-14.5: enclosing marks, fusion letters and von Baeyer
            // descriptors hold no Roman letters.
            (
                "[Na+].[Na+].[O-]C(=O)c1ccc(cc1)-c1ccccc1.CC(=O)[O-]",
                "disodium acetate [1,1'-biphenyl]-4-carboxylate",
                NameStatus::Preferred,
            ),
            (
                "[Na+].[Na+].[O-]C(=O)c1ccc2ccc3ccccc3c2c1.[O-]C(=O)c1ccc2c(c1)ccc1ccc3ccccc3c12",
                "disodium benzo[c]phenanthrene-3-carboxylate phenanthrene-3-carboxylate",
                NameStatus::Preferred,
            ),
            (
                "[Na+].[Na+].[O-]C(=O)[C@@H]1C[C@H]2CC[C@@H]1C2.CC(=O)[O-]",
                "disodium acetate (1R,2R,4S)-bicyclo[2.2.1]heptane-2-carboxylate",
                NameStatus::Preferred,
            ),
            // P-14.5.3, P-16.2.1: the italic ghi of a peri-fusion descriptor is
            // not Roman either; benzoindole before benzoperylene.
            (
                "[Na+].[Na+].C1(=CC=2C=3C4=C1C=CC1=CC=CC(C5=CC=CC(=CC2)C35)=C14)C(=O)[O-].C1C(=NC=2C=CC3=C(C12)C=CC=C3)C(=O)[O-]",
                "disodium 1H-benzo[e]indole-2-carboxylate benzo[ghi]perylene-1-carboxylate",
                NameStatus::Preferred,
            ),
            // P-14.5.4: primary locants outside the enclosing marks, 2 before 3.
            (
                "[Na+].[Na+].[O-]C(=O)C(CCl)CCCCC.[O-]C(=O)CC(CCl)CCCC",
                "disodium 2-(chloromethyl)heptanoate 3-(chloromethyl)heptanoate",
                NameStatus::Preferred,
            ),
            // P-14.5.3, P-16.2.4.1(d): the italic s of s-indacene is set off by
            // a hyphen and is no Roman letter; chloroindacenolate first.
            (
                "[Na+].[Na+].ClC1=c2cc3c(cc2C([O-])=C1)=CC=C3.[O-]c1ccc(Cl)c2ccccc12",
                "disodium 3-chloro-s-indacen-1-olate 4-chloronaphthalen-1-olate",
                NameStatus::Preferred,
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
        for smiles in [
            "C[N+](C)(C)CC[O-]",
            "[NH3+]CC([NH3+])C[O-]",
            // No prefix keeps these cations, so they stay suffixes beside -ide.
            "[CH2-]CC(=O)[NH3+]",
            "[CH2-]CC#[NH+]",
        ] {
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
            // P-74.0: the aminium cation becomes a prefix of the -ide parent.
            "[NH3+]CC[CH-]C#N",
            "[NH3+][C-]1C=CC2=C1N=NO2",
            "C[N+](C)(C)[C-](C)C",
        ] {
            assert_eq!(assess_component(smiles), NameStatus::Preferred, "{smiles}");
        }
    }

    #[test]
    // Provenance: fn tests::skeletal_anion_with_characteristic_cation_has_unsettled_suffix_order | original
    // Source: ReShiki | src/naming/openclatura/species.rs::tests::skeletal_anion_with_characteristic_cation_has_unsettled_suffix_order
    fn skeletal_anion_with_characteristic_cation_has_unsettled_suffix_order() {
        // Aminium becomes a prefix (P-74.0), but the diazonio prefix is not
        // the preselected diazyn-1-ium-1-yl (P-73.6), so diazonium stays.
        assert_eq!(
            assess_component("[CH2-]CC[N+]#N"),
            NameStatus::Unverified(NameReason::IonicSuffixOrder),
        );
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
            // Issue #92, P-16.3.6 (a): two mononuclear oxoanions take `bis`,
            // because `dicarbonate`, `disulfate`, `disulfite`, `diphosphate`
            // and `diborate` name the condensed dinuclear anions.
            (
                "[K+].[K+].[K+].[K+].[O-]C([O-])=O.[O-]C([O-])=O",
                "tetrapotassium bis(carbonate)",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "[Ca+2].[Ca+2].[O-]C([O-])=O.[O-]C([O-])=O",
                "dicalcium bis(carbonate)",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "[Na+].[Na+].[Na+].[Na+].[O-]S([O-])(=O)=O.[O-]S([O-])(=O)=O",
                "tetrasodium bis(sulfate)",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "[Na+].[Na+].[Na+].[Na+].[O-]S([O-])=O.[O-]S([O-])=O",
                "tetrasodium bis(sulfite)",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "[Na+].[Na+].[Na+].[Na+].[Na+].[Na+].[O-]P([O-])([O-])=O.[O-]P([O-])([O-])=O",
                "hexasodium bis(phosphate)",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "[Na+].[Na+].[Na+].[Na+].[Na+].[Na+].[O-]B([O-])[O-].[O-]B([O-])[O-]",
                "hexasodium bis(borate)",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "[Na+].[Na+].[Na+].[Na+].[Na+].[Na+].[O-]C([O-])=O.[O-]C([O-])=O.[O-]C([O-])=O",
                "hexasodium tris(carbonate)",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            // Issue #109, P-65.6.2.3.2, P-72.2.2.2.1.2: acid hydrogen of an
            // inorganic oxoanion is cited as `hydrogen`; P-16.3.6 (a) gives
            // `bis(hydrogen sulfate)` for two such anions.
            (
                "[Na+].OC(=O)[O-]",
                "sodium hydrogen carbonate",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "[Ca+2].OC(=O)[O-].OC(=O)[O-]",
                "calcium bis(hydrogen carbonate)",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "[Ca+2].OS(=O)(=O)[O-].OS(=O)(=O)[O-]",
                "calcium bis(hydrogen sulfate)",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "[Na+].[Na+].OP(=O)([O-])[O-]",
                "disodium hydrogen phosphate",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "[Na+].OP(=O)(O)[O-]",
                "sodium dihydrogen phosphate",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            // P-77.1.1 example.
            (
                "CC[NH+](CC)CC.OS(=O)(=O)[O-]",
                "N,N-diethylethanaminium hydrogen sulfate",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            // Issue #109, P-16.3.6 (a): `disulfide` and `dioxide` name the
            // S2(2-) and O2(2-) anions and `triiodide` names I3-, so separate
            // monatomic anions take `bis`/`tris`; `diiodide` stays.
            (
                "[Ca+2].[Ca+2].[S-2].[S-2]",
                "dicalcium bis(sulfide)",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "[Na+].[Na+].[Na+].[Na+].[O-2].[O-2]",
                "tetrasodium bis(oxide)",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "[Na+].[Na+].[S-2]",
                "disodium sulfide",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "[Na+].[Ca+2].[I-].[I-].[I-]",
                "calcium sodium tris(iodide)",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "[Na+].[Na+].[Na+].[Na+].[Na+].[Na+].[S-2].[S-2].[S-2]",
                "hexasodium tris(sulfide)",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "[Ca+2].[Ca+2].[Cl-].[Cl-].[Cl-].[Cl-]",
                "dicalcium tetrakis(chloride)",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "[Ca+2].[I-].[I-]",
                "calcium diiodide",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            // P-16.3.6 (a) applies only to ambiguous names: nitrite has no
            // condensed `dinitrite` anion and `diformate` names two formates.
            (
                "[Ca+2].[O-]N=O.[O-]N=O",
                "calcium dinitrite",
                NameStatus::Unverified(NameReason::UnauditedPath),
            ),
            (
                "[Ca+2].[O-]C=O.[O-]C=O",
                "calcium diformate",
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
            // Issue #92, P-16.3.6 (a): the hydrate keeps `bis(carbonate)`.
            (
                "[K+].[K+].[K+].[K+].[O-]C([O-])=O.[O-]C([O-])=O.O",
                "[tetrapotassium bis(carbonate)]\u{2014}water (1/1)",
                mixed,
                Some("tetrapotassium bis(carbonate) monohydrate"),
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
            (
                "CC(N)=O.CC#N",
                "acetamide\u{2014}acetonitrile (1/1)",
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
            // Classes whose order is not covered: an ester and an anhydride.
            "CS(=O)(=O)OC.CC(=O)OC(C)=O",
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
                "1-{4-[(N,N-dimethylmethanaminiumyl)methyl]phenyl}-N,N,N-trimethylmethanaminium diphenoxide",
                NotPreferred(Multiplicative),
            ),
            (
                format!("{PENTALENE_CARBOXYLATE}.{BIS_AMINIUM}.{PENTALENE_CARBOXYLATE}"),
                "1-{4-[(N,N-dimethylmethanaminiumyl)methyl]phenyl}-N,N,N-trimethylmethanaminium bis(octahydro-2,5-methanopentalene-1-carboxylate)",
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

    #[test]
    // Provenance: fn tests::cationic_ring_is_the_parent | original
    // Source: Nazuchem | crates/nazuchem-s2n/src/species.rs::tests::cationic_ring_is_the_parent
    // P-41 and P-70.2: cations are senior to every neutral class, so neutral
    // groups become prefixes on the cationic ring (P-73.1.2.1:
    // 2,2-dimethyl-1,3-dioxo-2,3-dihydro-1H-isoindol-2-ium). P-73.7(a): the
    // parent has the most cationic centres, including cationic suffixes.
    fn cationic_ring_is_the_parent() {
        for (smiles, expected) in [
            // Issue #32.
            (
                "OCc1cc[n+](C)cc1",
                "4-(hydroxymethyl)-1-methylpyridin-1-ium",
            ),
            ("OCC[n+]1ccccc1", "1-(2-hydroxyethyl)pyridin-1-ium"),
            ("C[n+]1ccc(C(=O)O)cc1", "4-carboxy-1-methylpyridin-1-ium"),
            (
                "C[n+]1ccc(CCN)cc1",
                "4-(2-aminoethyl)-1-methylpyridin-1-ium",
            ),
            ("O=C1CC[NH+](C)CC1", "1-methyl-4-oxopiperidin-1-ium"),
            // P-73.7(a) precedes ring seniority (P-44.2).
            (
                "C1CC[NH+](C1)CC2=C(N=C3N2C=CS3)C4=CC=CC=C4",
                "1-[(6-phenylimidazo[2,1-b][1,3]thiazol-5-yl)methyl]pyrrolidin-1-ium",
            ),
            // Cationic suffixes keep their place (P-73.7(a), (d)).
            (
                "[NH3+]CCc1cc[n+](C)cc1",
                "2-(1-methylpyridin-1-ium-4-yl)ethan-1-aminium",
            ),
            ("OCC[NH3+]", "2-hydroxyethan-1-aminium"),
            ("OCC1CC[NH2+]C1", "3-(hydroxymethyl)pyrrolidin-1-ium"),
            // Nitro is a neutral group, not a cationic centre.
            (
                "OCc1cc[n+](cc1)c1ccc([N+](=O)[O-])cc1",
                "4-(hydroxymethyl)-1-(4-nitrophenyl)pyridin-1-ium",
            ),
            // P-73.5.3.2, P-75.2: a skeletal cationic centre takes the low
            // locant before detachable prefixes (P-14.4(f)).
            ("C[NH+]1CCN(CC)CC1", "4-ethyl-1-methylpiperazin-1-ium"),
            (
                "CCN1CC[NH+](CC1)CCNC2=NC3=CC(=NN3C4=CC=CC=C42)C",
                "4-ethyl-1-{2-[(2-methylpyrazolo[1,5-a]quinazolin-5-yl)amino]ethyl}piperazin-1-ium",
            ),
            // P-73.7(c): a cationic N precedes a cationic O before ring
            // seniority, which would prefer morpholine (P-44.2.1).
            (
                "C[NH+]1CCC(CC2C[O+](C)CCN2)CC1",
                "1-methyl-4-[(1-methylmorpholin-1-ium-3-yl)methyl]piperidin-1-ium",
            ),
        ] {
            let name = super::super::name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, NameStatus::Preferred, "{smiles}");
        }
        // A charge-separated N-oxide is not a cationic centre (P-74).
        let name = super::super::name_smiles("[O-][n+]1ccc(CO)cc1", &Default::default()).unwrap();
        assert_eq!(name.text, "(1-oxidopyridin-1-ium-4-yl)methanol");
    }
}
