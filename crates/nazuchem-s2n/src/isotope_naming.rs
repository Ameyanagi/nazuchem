// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: ReShiki | src/naming/openclatura/isotope_naming.rs | 327d32fed63a5b3df15b7b8dc5ab84a7c36e0865
// Copyright (c) 2026 Ameyanagi and ReShiki contributors

//! Isotopically substituted structures (P-82), never sample-label notation.
//!
//! The component engine supplies the parent locants. Only those positions,
//! unique heteroatoms and completely substituted elements are supported. A
//! bounded automorphism proof chooses the lowest equivalent isotope locants;
//! the finished name must independently reconstruct the original graph.
use super::{
    Name, component_namer,
    locants::retained_locant_sort_key,
    molecule::{AtomSet, Molecule},
};
use nazuchem_ir::{IsotopeStyle, NameReason, NameStatus};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone, Debug)]
pub(super) struct Site {
    pub(super) anchor: usize,
    symbol: String,
    mass: u16,
}

impl Site {
    /// P-82.5.2: the higher atomic number, then the higher mass number.
    pub(super) fn priority(&self) -> (usize, u16) {
        let element = super::rules::elements::get(&self.symbol);
        (element.map_or(0, |e| e.atomic_number), self.mass)
    }
}

/// The descriptor that a prefix name starts with, as in (2H3)methyl.
pub(super) fn leading_descriptor<'a>(name: &str, descriptors: &'a [String]) -> &'a str {
    descriptors
        .iter()
        .find(|descriptor| name.starts_with(descriptor.as_str()))
        .map_or("", String::as_str)
}

/// P-82.5.1: nuclides are considered last, so alphanumerical order and
/// numbering use the prefix name without its descriptors.
pub(super) fn citation_name(name: &str, descriptors: &[String]) -> String {
    descriptors
        .iter()
        .fold(name.to_owned(), |name, descriptor| {
            name.replace(&format!("{descriptor}-"), "")
                .replace(descriptor.as_str(), "")
        })
}

/// Keep CIP labels calculated with the isotopes present. Recalculating them
/// after stripping can invert an existing centre or erase an isotope-created
/// centre. The final strict parser check refuses any unsupported stereo path.
fn stripped_parent(source: &Molecule) -> Result<(Molecule, Vec<Site>), String> {
    let mut parent = Molecule::new();
    parent.isotope_labels_stripped = true;
    let mut sites = Vec::new();
    for (&idx, atom) in &source.atoms {
        if atom.symbol == "H" {
            let neighbors = source.get_neighbors(idx);
            if atom.charge != 0
                || neighbors.len() != 1
                || source.atoms[&neighbors[0]].symbol == "H"
                || source.bond_order(idx, neighbors[0]) != 1
            {
                return Err(
                    "openclatura: this isotopic hydrogen attachment is not supported".into(),
                );
            }
            if let Some(mass) = atom.isotope {
                sites.push(Site {
                    anchor: neighbors[0],
                    symbol: "H".into(),
                    mass,
                });
            }
            continue;
        }
        let mut stripped = atom.clone();
        if let Some(mass) = stripped.isotope.take() {
            sites.push(Site {
                anchor: idx,
                symbol: atom.symbol.clone(),
                mass,
            });
        }
        let attached_h = source
            .get_neighbors(idx)
            .iter()
            .filter(|n| source.atoms[n].symbol == "H")
            .count();
        stripped.explicit_h_count += attached_h;
        stripped.total_h_count += attached_h;
        parent.add_atom(stripped)?;
    }
    for bond in source.bonds.values() {
        if parent.atoms.contains_key(&bond.u) && parent.atoms.contains_key(&bond.v) {
            parent.add_bond(bond.clone())?;
        }
    }
    parent.accurate_cip = source.accurate_cip.clone();
    parent.legacy_cip = source.legacy_cip.clone();
    // Like the CIP labels, the possible stereo units keep isotope-created centres.
    parent.stereogenic_atoms = source
        .stereogenic_atoms
        .iter()
        .copied()
        .filter(|idx| parent.atoms.contains_key(idx))
        .collect();
    parent.stereogenic_bonds = source
        .stereogenic_bonds
        .iter()
        .copied()
        .filter(|idx| parent.bonds.contains_key(idx))
        .collect();
    if parent.atoms.is_empty() {
        return Err("openclatura: isolated isotope atoms are not supported".into());
    }
    Ok((parent, sites))
}

fn same_atom(mol: &Molecule, a: usize, b: usize) -> bool {
    let (a, b) = (&mol.atoms[&a], &mol.atoms[&b]);
    a.symbol == b.symbol
        && a.charge == b.charge
        && a.total_h_count == b.total_h_count
        && a.is_aromatic == b.is_aromatic
        && a.stereo == b.stereo
        && (matches!(a.stereo.as_deref(), Some("R" | "S")) || a.raw_stereo == b.raw_stereo)
        && mol.degree(a.idx) == mol.degree(b.idx)
}

fn same_bond(mol: &Molecule, a: usize, b: usize, x: usize, y: usize) -> bool {
    match (mol.get_bond(a, b), mol.get_bond(x, y)) {
        (None, None) => true,
        (Some(a), Some(b)) => {
            a.is_aromatic == b.is_aromatic
                && (a.is_aromatic || a.order == b.order)
                && a.stereo == b.stereo
        }
        _ => false,
    }
}

/// Complete labeled-graph automorphisms, rather than relying on source atom
/// order or guessing a chain direction from a SMILES string. Refuse if proving
/// completeness would exceed the request budget. Pinned atoms stay in place, so
/// that the descriptors of separate name parts describe one and the same
/// labelled structure.
fn lowest_equivalent_map(
    mol: &Molecule,
    sites: &[Site],
    locants: &BTreeMap<usize, String>,
    pinned: &AtomSet,
) -> Result<BTreeMap<usize, usize>, String> {
    mol.spend_naming_work(mol.atoms.len().saturating_mul(mol.atoms.len()))?;
    let candidates: BTreeMap<_, Vec<_>> = mol
        .atoms
        .keys()
        .map(|&a| {
            (
                a,
                mol.atoms
                    .keys()
                    .copied()
                    .filter(|&b| a == b || !pinned.contains(&a) && same_atom(mol, a, b))
                    .collect(),
            )
        })
        .collect();
    let mut order: Vec<_> = mol.atoms.keys().copied().collect();
    order.sort_by_key(|a| (candidates[a].len(), std::cmp::Reverse(mol.degree(*a))));
    struct Search<'a> {
        mol: &'a Molecule,
        candidates: &'a BTreeMap<usize, Vec<usize>>,
        order: &'a [usize],
        sites: &'a [Site],
        locants: &'a BTreeMap<usize, String>,
        best: Option<(IsotopeKey, BTreeMap<usize, usize>)>,
    }
    impl Search<'_> {
        fn visit(
            &mut self,
            mapping: &mut BTreeMap<usize, usize>,
            used: &mut BTreeSet<usize>,
        ) -> Result<(), String> {
            self.mol.spend_naming_work(1)?;
            if mapping.len() == self.order.len() {
                // A mapping that moves a label onto an atom without a locant
                // is not a usable numbering; keep searching the others.
                let Ok(key) = numbering_key(self.mol, self.sites, self.locants, mapping) else {
                    return Ok(());
                };
                if self
                    .best
                    .as_ref()
                    .is_none_or(|(previous, _)| &key < previous)
                {
                    self.best = Some((key, mapping.clone()));
                }
                return Ok(());
            }
            let atom = self.order[mapping.len()];
            for &target in &self.candidates[&atom] {
                if used.contains(&target)
                    || !mapping
                        .iter()
                        .all(|(&a, &b)| same_bond(self.mol, atom, a, target, b))
                {
                    continue;
                }
                mapping.insert(atom, target);
                used.insert(target);
                self.visit(mapping, used)?;
                used.remove(&target);
                mapping.remove(&atom);
            }
            Ok(())
        }
    }
    let mut search = Search {
        mol,
        candidates: &candidates,
        order: &order,
        sites,
        locants,
        best: None,
    };
    search.visit(&mut BTreeMap::new(), &mut BTreeSet::new())?;
    search.best.map(|(_, map)| map).ok_or_else(uncited)
}

fn uncited() -> String {
    "openclatura: an isotope label on a characteristic group or on this kind of substituent prefix is not supported".into()
}

fn site_locant(
    mol: &Molecule,
    map: &BTreeMap<usize, String>,
    anchor: usize,
) -> Result<String, String> {
    if mol.atoms.len() == 1 {
        return Ok("1".into());
    }
    if let Some(locant) = map.get(&anchor) {
        return Ok(locant.clone());
    }
    let symbol = &mol.atoms[&anchor].symbol;
    // Letter locants can identify a unique characteristic heteroatom, but
    // neither branch atom numbers nor repeated suffix atoms can be inferred.
    if ["N", "O", "S", "Se", "Te", "P", "As"].contains(&symbol.as_str())
        && mol.atoms.values().filter(|a| &a.symbol == symbol).count() == 1
    {
        return Ok(symbol.clone());
    }
    Err("openclatura: isotope position has no proven parent or unique heteroatom locant".into())
}

type IsotopeKey = Vec<Vec<super::locants::LocantKey>>;
fn numbering_key(
    mol: &Molecule,
    sites: &[Site],
    locants: &BTreeMap<usize, String>,
    symmetry: &BTreeMap<usize, usize>,
) -> Result<IsotopeKey, String> {
    let mut all = Vec::new();
    let mut priority = BTreeMap::<(usize, u16), Vec<_>>::new();
    for site in sites {
        let locant = site_locant(mol, locants, symmetry[&site.anchor])?;
        let key = retained_locant_sort_key(&locant);
        all.push(key.clone());
        let atomic_number = super::rules::elements::get(&site.symbol)
            .ok_or("openclatura: unknown isotope element")?
            .atomic_number;
        priority
            .entry((atomic_number, site.mass))
            .or_default()
            .push(key);
    }
    all.sort();
    let mut key = vec![all];
    // P-82.5.2: first all modified positions together, then the higher atomic
    // number and, for one element, the higher mass number.
    for (_, mut locants) in priority.into_iter().rev() {
        locants.sort();
        key.push(locants);
    }
    Ok(key)
}

/// The name part that a descriptor precedes, for the omission of locants.
struct Unit<'a> {
    /// The part with its substituents: no locants when all of it is modified
    /// in the same way (P-82.6.1.3),
    atoms: &'a AtomSet,
    /// counting the labels that the descriptors of its prefixes cite.
    cited: &'a [Site],
    /// Its parent: a single atom of an element needs no locant (P-82.6.1.2).
    parent: &'a AtomSet,
}

fn isotope_descriptors(
    mol: &Molecule,
    sites: &[Site],
    locants: &BTreeMap<usize, String>,
    symmetry: &BTreeMap<usize, usize>,
    unit: &Unit,
    style: IsotopeStyle,
) -> Result<String, String> {
    let mut groups = BTreeMap::<(&str, u16), Vec<usize>>::new();
    for site in sites {
        groups
            .entry((&site.symbol, site.mass))
            .or_default()
            .push(symmetry[&site.anchor]);
    }
    if style == IsotopeStyle::Boughton && groups.len() != 1 {
        return Err("openclatura: mixed-nuclide Boughton notation is not supported by the local parser; use IUPAC isotope style".into());
    }
    let mut descriptors = Vec::new();
    for ((symbol, mass), anchors) in groups {
        let atoms = unit.atoms.iter().map(|atom| &mol.atoms[atom]);
        let available = if symbol == "H" {
            atoms.map(|a| a.total_h_count).sum()
        } else {
            atoms.filter(|a| a.symbol == symbol).count()
        } - unit
            .cited
            .iter()
            .filter(|site| site.symbol == symbol && site.mass == mass)
            .count();
        let single = symbol != "H"
            && unit
                .parent
                .iter()
                .filter(|atom| mol.atoms[*atom].symbol == symbol)
                .count()
                == 1;
        let only_carbon_parent = symbol == "H"
            && locants.len() == 1
            && anchors
                .iter()
                .all(|anchor| locants.contains_key(anchor) && mol.atoms[anchor].symbol == "C");
        let omit = available == anchors.len() || single || only_carbon_parent;
        let mut positions = if omit {
            vec![]
        } else {
            anchors
                .iter()
                .map(|&a| site_locant(mol, locants, a))
                .collect::<Result<Vec<_>, _>>()?
        };
        positions.sort_by_key(|s| retained_locant_sort_key(s));
        let position = if positions.is_empty() {
            String::new()
        } else {
            format!("{}-", positions.join(","))
        };
        let nuclide = if style == IsotopeStyle::Boughton && symbol == "H" && mass == 2 {
            "d".into()
        } else {
            format!("{mass}{symbol}")
        };
        let count = if symbol == "H" || anchors.len() > 1 {
            anchors.len().to_string()
        } else {
            String::new()
        };
        descriptors.push(format!("{position}{nuclide}{count}"));
    }
    Ok(match style {
        IsotopeStyle::Iupac => format!("({})", descriptors.join(",")),
        IsotopeStyle::Boughton => format!("-{}", descriptors.join("-")),
    })
}

fn with_descriptor(
    base: &str,
    rendered_parent: &str,
    descriptor: &str,
    style: IsotopeStyle,
) -> Result<String, String> {
    if style == IsotopeStyle::Boughton {
        return Ok(format!("{base}{descriptor}"));
    }
    // P-82.2.1: substituent prefixes precede the descriptor; the engine's
    // rendered parent identifies the boundary without interpreting the name.
    let prefix = base.strip_suffix(rendered_parent).ok_or_else(|| {
        "openclatura: isotope parent boundary is not preserved by ionic naming".to_owned()
    })?;
    let hyphen = if rendered_parent.starts_with(|c: char| c.is_ascii_digit()) {
        "-"
    } else {
        ""
    };
    Ok(format!("{prefix}{descriptor}{hyphen}{rendered_parent}"))
}

/// P-82.2.1: a substituent prefix cites its own labels before its parent, as
/// in 4-(2-14C)ethylbenzoic acid and 4-[(3-14C)thiolan-2-yl]pyridine. Returns
/// the descriptor and the labelled atoms of the prefix parent numbered by
/// `locants`; `branch` is the whole prefix, the only part a numbering may
/// move, and `nested` holds the labelled atoms that its own prefixes cite.
pub(super) fn prefix_descriptor(
    mol: &Molecule,
    locants: &BTreeMap<usize, String>,
    branch: &AtomSet,
    nested: &AtomSet,
) -> Result<Option<(String, AtomSet)>, String> {
    let (sites, other): (Vec<_>, Vec<_>) = mol
        .isotope_sites()
        .into_iter()
        .partition(|site| locants.contains_key(&site.anchor));
    if sites.is_empty() {
        return Ok(None);
    }
    let pinned = mol
        .atoms
        .keys()
        .copied()
        .filter(|atom| !branch.contains(atom))
        .chain(other.iter().map(|site| site.anchor))
        .collect();
    let symmetry = lowest_equivalent_map(mol, &sites, locants, &pinned)?;
    let cited: Vec<_> = other
        .into_iter()
        .filter(|site| nested.contains(&site.anchor))
        .collect();
    let unit = Unit {
        atoms: branch,
        cited: &cited,
        parent: &locants.keys().copied().collect(),
    };
    let descriptor =
        isotope_descriptors(mol, &sites, locants, &symmetry, &unit, IsotopeStyle::Iupac)?;
    Ok(Some((
        descriptor,
        sites.iter().map(|site| site.anchor).collect(),
    )))
}

#[cfg(test)]
fn verify(smiles: &str, name: &str) -> Result<(), String> {
    let parser = nazuchem_n2s::Parser::new().map_err(|e| e.to_string())?;
    let result = parser.parse(name, &nazuchem_n2s::ParseOptions::strict());
    if result
        .warnings
        .iter()
        .any(|warning| warning.kind == nazuchem_n2s::WarningKind::StereochemistryIgnored)
    {
        return Err(
            "openclatura: isotope naming could not preserve the specified CIP stereochemistry"
                .into(),
        );
    }
    let structure = result.structure.ok_or_else(|| {
        format!(
            "openclatura: isotope descriptor placement is not supported by the local parser: {}",
            result.message
        )
    })?;
    let reconstructed = structure.semantic_cxsmiles().map_err(|e| e.to_string())?;
    crate::test_support::verify_identity(smiles, &reconstructed)
        .map_err(|_| "openclatura: isotope descriptor does not exactly reconstruct the labelled structure and stereochemistry".into())
}

/// Name the already admitted graph without serializing its atom identities.
pub(super) fn name_molecule(source: &Molecule, style: IsotopeStyle) -> Result<Name, String> {
    let (parent, sites) = stripped_parent(source)?;
    if sites.is_empty() {
        return Err("openclatura: isotope naming requires an isotope".into());
    }
    if parent
        .atoms
        .values()
        .any(|atom| atom.radical_electrons != 0)
    {
        // P-71.2.1 names the radical like a substituent prefix, which cites
        // its labels by P-82.2.1: (2H3)methyl.
        parent.set_isotope_sites(sites);
        return super::radical_naming::name_radical(&parent)?
            .filter(|_| style == IsotopeStyle::Iupac && parent.isotope_sites().is_empty())
            .ok_or_else(|| {
                "openclatura: combined isotope and radical nomenclature is not supported here"
                    .into()
            });
    }
    if parent.atoms.values().any(|atom| atom.symbol == "Po") {
        return Err("openclatura: polonium nomenclature is not validated".into());
    }
    super::validate_parent(&parent)?;
    parent.mark_isotope_parent();
    parent.set_isotope_sites(sites.clone());
    let atoms: AtomSet = parent.atoms.keys().copied().collect();
    let base = component_namer::name_component_with_atom_locants(&parent, &atoms, Some(source))?;
    let locants = &base.locants;
    super::assess_derivation(&parent, &base.derivation);
    // P-82.2.1: the substituent prefixes have cited their own labels, and a
    // label left in a prefix has no place in the parent's descriptor.
    let (cited, sites): (Vec<_>, Vec<_>) = sites
        .into_iter()
        .partition(|site| base.prefix_isotopes.contains(&site.anchor));
    if sites
        .iter()
        .any(|site| base.prefix_atoms.contains(&site.anchor))
    {
        return Err(uncited());
    }
    let text = if sites.is_empty() {
        base.text
    } else {
        let symmetry = lowest_equivalent_map(&parent, &sites, locants, &base.prefix_isotopes)?;
        let unit = Unit {
            atoms: &atoms,
            cited: &cited,
            parent: &atoms,
        };
        let descriptor = isotope_descriptors(&parent, &sites, locants, &symmetry, &unit, style)?;
        with_descriptor(&base.text, &base.rendered_parent, &descriptor, style)?
    };
    if style == IsotopeStyle::Boughton && !base.prefix_isotopes.is_empty() {
        return Err("openclatura: Boughton notation cannot place a label on a substituent prefix; use IUPAC isotope style".into());
    }
    let text = super::nest_enclosing_marks(&text);
    let status = if style == IsotopeStyle::Boughton {
        NameStatus::NotPreferred(NameReason::BoughtonIsotope)
    } else {
        parent.mark_name(NameStatus::Unverified(NameReason::IsotopePlacement));
        parent.name_status()
    };
    Ok(Name {
        text,
        status,
        alternatives: vec![],
        omitted_stereo: vec![],
        omitted_units: 0,
        verification: None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::name_smiles_with_style as name_smiles;

    #[test]
    fn stripped_parents_keep_isotope_provenance_and_original_cip() {
        for smiles in [
            "[13CH3][C@H](F)C([C@H](F)C)CCCC",
            "[2H][C@H](F)C([C@H](F)C)CCCC",
        ] {
            let source = crate::graph_io::read_smiles(smiles, &Default::default()).unwrap();
            let (parent, _) = stripped_parent(&source).unwrap();
            assert!(parent.atoms.values().all(|a| a.isotope.is_none()));
            assert!(parent.isotope_labels_stripped);
            assert_eq!(parent.accurate_cip, source.accurate_cip);
            assert_eq!(parent.legacy_cip, source.legacy_cip);
            for (idx, atom) in &parent.atoms {
                assert_eq!(atom.stereo, source.atoms[idx].stereo);
            }
            assert!(parent.clone().isotope_labels_stripped);
            assert!(
                parent
                    .subgraph(&parent.atoms.keys().copied().collect())
                    .isotope_labels_stripped
            );
        }
    }

    #[test]
    fn isotopic_names_reconstruct_exact_graphs() {
        for smiles in [
            "[13CH4]",
            "[2H]C",
            "[2H]C([2H])([2H])O",
            "CO[2H]",
            "[2H]C([2H])(C)O",
            "[2H]C([2H])([2H])O[2H]",
            "C[13CH2]CO",
            "CC[13CH2]C",
            "[15NH2]C",
            "C[18OH]",
            "[2H][C@H](O)C",
            "[2H][C@@H](O)C",
        ] {
            for style in [IsotopeStyle::Iupac, IsotopeStyle::Boughton] {
                let name = name_smiles(smiles, style, &Default::default())
                    .unwrap_or_else(|e| panic!("{smiles} {style:?}: {e}"));
                verify(smiles, &name.text).unwrap();
                assert_eq!(
                    name.status,
                    if style == IsotopeStyle::Boughton {
                        NameStatus::NotPreferred(NameReason::BoughtonIsotope)
                    } else {
                        NameStatus::Unverified(NameReason::IsotopePlacement)
                    }
                );
            }
        }
    }

    #[test]
    fn nuclides_are_ordered_by_symbol_then_increasing_mass() {
        for (smiles, expected) in [
            ("[13CH3]C([2H])([2H])C", "(1-13C,2,2-2H2)propane"),
            ("[2H]C([3H])CO", "(2-2H1,2-3H1)ethan-1-ol"),
            ("C[13CH2]C([2H])C", "(2-13C,3-2H1)butane"),
            ("CC([2H])C([3H])C", "(3-2H1,2-3H1)butane"),
        ] {
            let name = name_smiles(smiles, IsotopeStyle::Iupac, &Default::default()).unwrap();
            assert_eq!(name.text, expected);
            verify(smiles, &name.text).unwrap();
            assert!(name_smiles(smiles, IsotopeStyle::Boughton, &Default::default()).is_err());
        }
    }

    #[test]
    fn engine_locants_and_equivalent_numbering_ignore_atom_order() {
        for (a, b) in [
            ("CC[13CH2]C", "C[13CH2]CC"),
            ("[2H]CCC", "CCC[2H]"),
            ("[13CH3]CO", "OC[13CH3]"),
        ] {
            let a = name_smiles(a, IsotopeStyle::Iupac, &Default::default())
                .unwrap()
                .text;
            let b = name_smiles(b, IsotopeStyle::Iupac, &Default::default())
                .unwrap()
                .text;
            assert_eq!(a, b);
        }
        assert_eq!(
            name_smiles("CC[13CH2]C", IsotopeStyle::Iupac, &Default::default())
                .unwrap()
                .text,
            "(2-13C)butane"
        );
        assert_eq!(
            name_smiles("[13CH3]CO", IsotopeStyle::Iupac, &Default::default())
                .unwrap()
                .text,
            "(2-13C)ethan-1-ol"
        );
        assert_eq!(
            name_smiles(
                "[2H]C([2H])([2H])O",
                IsotopeStyle::Iupac,
                &Default::default()
            )
            .unwrap()
            .text,
            "(2H3)methanol"
        );
        assert_eq!(
            name_smiles("CO[2H]", IsotopeStyle::Iupac, &Default::default())
                .unwrap()
                .text,
            "(O-2H1)methanol"
        );
        assert_eq!(
            name_smiles(
                "[2H]C([2H])(C)O",
                IsotopeStyle::Boughton,
                &Default::default()
            )
            .unwrap()
            .text,
            "ethan-1-ol-1,1-d2"
        );
    }

    #[test]
    fn isotope_descriptors_follow_substituents_and_stereodescriptors() {
        for (smiles, expected) in [
            ("[13CH3]CCl", "1-chloro(2-13C)ethane"),
            ("ClC([2H])([2H])Cl", "dichloro(2H2)methane"),
            ("[13CH3][C@@H](F)Cl", "(1S)-1-chloro-1-fluoro(2-13C)ethane"),
            // P-91.3: the only stereocentre of an unlocanted name.
            ("F[13C@H](Cl)Br", "(R)-bromochlorofluoro(13C)methane"),
        ] {
            let name = name_smiles(smiles, IsotopeStyle::Iupac, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            verify(smiles, &name.text).unwrap();
        }
    }

    #[test]
    fn labelled_amide_keeps_the_labelled_parent() {
        // The cross-key amide parent would move the label into a substituent.
        for (smiles, expected) in [
            ("N[13C](=O)CNC(=O)c1ccccc1", "2-benzamido(1-13C)acetamide"),
            ("NC(=O)[13CH2]NC(=O)c1ccccc1", "2-benzamido(2-13C)acetamide"),
        ] {
            let name = name_smiles(smiles, IsotopeStyle::Iupac, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            verify(smiles, &name.text).unwrap();
        }
    }

    #[test]
    fn ambiguous_branch_and_suffix_positions_are_refused() {
        for smiles in ["C(=[18O])O", "[2H][2H]", "[13CH3]"] {
            // Generation keeps only names that the parser reconstructs exactly.
            assert!(
                name_smiles(smiles, IsotopeStyle::Iupac, &Default::default())
                    .and_then(|name| verify(smiles, &name.text))
                    .is_err(),
                "{smiles}"
            );
        }
    }

    #[test]
    fn substituent_prefixes_cite_their_own_labels() {
        // Issue #27. P-82.2.1: the descriptor precedes the part of the name
        // that is isotopically substituted.
        for (smiles, expected) in [
            // P-82.2.1 examples.
            ("[2H]C([2H])([2H])Oc1ccccc1", "(2H3)methoxybenzene"),
            ("[13CH3]c1ccccc1[13CH3]", "1,2-di[(13C)methyl]benzene"),
            ("[13CH3]c1[13cH]cccc1", "1-(13C)methyl(2-13C)benzene"),
            ("OC1(CCCC1)[14CH2]N", "1-[amino(14C)methyl]cyclopentan-1-ol"),
            ("c1cnccc1C1SCC[14CH2]1", "4-[(3-14C)thiolan-2-yl]pyridine"),
            // P-82.2.4 example.
            ("OC(=O)c1ccc(C[14CH3])cc1", "4-(2-14C)ethylbenzoic acid"),
            // P-82.2.2.1 examples: the modified prefix is cited first.
            ("[13CH3]c1ncccc1C", "2-(13C)methyl-3-methylpyridine"),
            (
                "OCC(CC([2H])[2H])C(CC)CCC",
                "2-(2,2-2H2)ethyl-3-ethylhexan-1-ol",
            ),
            // P-82.6.1.1 example.
            (
                "[2H]COC([2H])([2H])SCOO",
                "{[(2H1)methoxy(2H2)methyl]sulfanyl}methaneperoxol",
            ),
            // P-82.6.1.2, P-82.6.1.3: no locants for a single atom of the
            // element or for a completely substituted prefix.
            ("[11CH3]Oc1ccccc1", "(11C)methoxybenzene"),
            ("[13CH3]C1CCCCC1", "(13C)methylcyclohexane"),
            (
                "[2H]c1c([2H])c([2H])c([2H])c([2H])c1-c1ncncn1",
                "2-(2H5)phenyl-1,3,5-triazine",
            ),
            (
                "CC(=O)[14c]1[14cH][14cH][14cH][14cH][14cH]1",
                "1-(14C6)phenylethan-1-one",
            ),
            // P-82.6.3.2: a labelled position keeps its numbered parent.
            (
                "[13CH3]C(C)(C)c1ccccc1",
                "[2-methyl(1-13C)propan-2-yl]benzene",
            ),
            ("c1ccc(cc1)[13CH2]Oc1ccccc1", "[phenoxy(13C)methyl]benzene"),
            // The PubChem rows of the issue (CID 155606068, CID 10475160).
            (
                "[2H]C([2H])(C)N1C2=C(C=C(C=C2)C(=O)OC)C=C1CC3=C(C=C(C=C3)Cl)C(F)(F)F",
                "methyl 2-{[4-chloro-2-(trifluoromethyl)phenyl]methyl}-1-(1,1-2H2)ethyl-1H-indole-5-carboxylate",
            ),
            (
                "CCOC(=O)NC1=NC(=C(C(=C1)N[C@@H](C)C(=O)[14C]2=[14CH][14CH]=[14CH][14CH]=[14CH]2)[N+](=O)[O-])N",
                "ethyl (6-amino-5-nitro-4-{[(2S)-1-oxo-1-(14C6)phenylpropan-2-yl]amino}pyridin-2-yl)carbamate",
            ),
        ] {
            let name = name_smiles(smiles, IsotopeStyle::Iupac, &Default::default())
                .unwrap_or_else(|e| panic!("{smiles}: {e}"));
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(
                name.status,
                NameStatus::Unverified(NameReason::IsotopePlacement),
                "{smiles}"
            );
            verify(smiles, &name.text).unwrap();
            // Boughton notation has no place inside a prefix.
            assert!(
                name_smiles(smiles, IsotopeStyle::Boughton, &Default::default()).is_err(),
                "{smiles}"
            );
        }
    }

    #[test]
    fn labelled_prefixes_take_the_lowest_equivalent_locants() {
        // P-82.5.2 within the prefix, whatever the atom order of the input.
        for (smiles, expected) in [
            (
                "OC1CCC(CC1)c1c([2H])cccc1",
                "4-(2-2H1)phenylcyclohexan-1-ol",
            ),
            ("OC1CCC(CC1)c1ccccc1[2H]", "4-(2-2H1)phenylcyclohexan-1-ol"),
            ("c1ccccc1C([13CH3])C", "[(1-13C)propan-2-yl]benzene"),
            ("c1ccccc1C(C)[13CH3]", "[(1-13C)propan-2-yl]benzene"),
            (
                "OC1CCC(CC1)c1cc([2H])cc(C)c1",
                "4-[3-methyl(5-2H1)phenyl]cyclohexan-1-ol",
            ),
        ] {
            let name = name_smiles(smiles, IsotopeStyle::Iupac, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            verify(smiles, &name.text).unwrap();
        }
    }

    #[test]
    fn labelled_prefixes_take_the_lower_parent_locant() {
        // P-82.5.2: lowest locants to the modified groups, in either atom order.
        for (smiles, expected) in [
            ("[13CH3]c1ccccc1C", "1-(13C)methyl-2-methylbenzene"),
            ("Cc1ccccc1[13CH3]", "1-(13C)methyl-2-methylbenzene"),
            ("[13CH3]c1cccc(C)c1", "1-(13C)methyl-3-methylbenzene"),
            ("Cc1cccc([13CH3])c1", "1-(13C)methyl-3-methylbenzene"),
            ("CCC1CCC(C[13CH3])CC1", "1-(2-13C)ethyl-4-ethylcyclohexane"),
            ("[13CH3]CC1CCC(CC)CC1", "1-(2-13C)ethyl-4-ethylcyclohexane"),
            // Then the higher atomic number, then the higher mass number.
            (
                "[13CH3]c1ccccc1[14CH3]",
                "1-(14C)methyl-2-(13C)methylbenzene",
            ),
            (
                "[14CH3]c1ccccc1[13CH3]",
                "1-(14C)methyl-2-(13C)methylbenzene",
            ),
            ("[2H]Cc1ccccc1[13CH3]", "1-(13C)methyl-2-(2H1)methylbenzene"),
            ("[13CH3]c1ccccc1C[2H]", "1-(13C)methyl-2-(2H1)methylbenzene"),
            // P-82.5.1: a descriptor inside a prefix changes neither the
            // alphanumerical order nor the numbering of the unmodified compound.
            (
                "OC1CC([13CH2]N)CC(CCN)C1",
                "3-(2-aminoethyl)-5-[amino(13C)methyl]cyclohexan-1-ol",
            ),
            (
                "OC1CC(CCN)CC([13CH2]N)C1",
                "3-(2-aminoethyl)-5-[amino(13C)methyl]cyclohexan-1-ol",
            ),
        ] {
            let name = name_smiles(smiles, IsotopeStyle::Iupac, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            verify(smiles, &name.text).unwrap();
        }
    }

    #[test]
    fn equal_parents_prefer_the_isotopically_modified_one() {
        // P-44.4.1.11.1: the parent with more modified atoms, in either atom order.
        for (smiles, expected) in [
            ("OC(=O)CC(C)(C)[13CH3]", "3,3-dimethyl(4-13C)butanoic acid"),
            ("OC(=O)CC([13CH3])(C)C", "3,3-dimethyl(4-13C)butanoic acid"),
            ("CC(C)C(C)[13CH3]", "2,3-dimethyl(1-13C)butane"),
            ("[13CH3]C(C)C(C)C", "2,3-dimethyl(1-13C)butane"),
            ("CC(C)(C[2H])C(=O)O", "2,2-dimethyl(3-2H1)propanoic acid"),
            (
                "CC(C)([13CH3])c1ccccc1",
                "[2-methyl(1-13C)propan-2-yl]benzene",
            ),
            ("CCC(CC)(CC[2H])O", "3-ethyl(1-2H1)pentan-3-ol"),
        ] {
            let name = name_smiles(smiles, IsotopeStyle::Iupac, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            verify(smiles, &name.text).unwrap();
        }
    }

    #[test]
    fn completely_substituted_parts_omit_their_locants() {
        // P-82.6.1.3: the labels a prefix cites count towards the completely
        // substituted compound; an unlabelled hydrogen anywhere keeps locants.
        for (smiles, expected) in [
            (
                "[13CH3][13CH]([13CH3])[13CH3]",
                "2-(13C)methyl(13C3)propane",
            ),
            (
                "[2H]C([2H])([2H])c1c([2H])c([2H])c([2H])c([2H])c1[2H]",
                "1-(2H3)methyl(2H5)benzene",
            ),
            (
                "[2H]Cc1c([2H])c([2H])c([2H])c([2H])c1[2H]",
                "1-(2H1)methyl(2,3,4,5,6-2H5)benzene",
            ),
            (
                "Cc1c([2H])c([2H])c([2H])c([2H])c1[2H]",
                "1-methyl(2,3,4,5,6-2H5)benzene",
            ),
            // The same inside a prefix: its own substituents count.
            (
                "OC1CCC(CC1)c1c([2H])c([2H])c(C)c([2H])c1[2H]",
                "4-[4-methyl(2,3,5,6-2H4)phenyl]cyclohexan-1-ol",
            ),
            (
                "OC1CCC(CC1)c1c([2H])c([2H])c(C([2H])([2H])[2H])c([2H])c1[2H]",
                "4-[4-(2H3)methyl(2H4)phenyl]cyclohexan-1-ol",
            ),
            // P-82.6.1.2: one atom of the element in the prefix's parent.
            ("c1ccccc1[13CH2]c1ccncc1", "4-[phenyl(13C)methyl]pyridine"),
        ] {
            let name = name_smiles(smiles, IsotopeStyle::Iupac, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            verify(smiles, &name.text).unwrap();
        }
    }

    #[test]
    fn labels_that_no_name_part_cites_still_decline() {
        // Characteristic groups and the substituents of a suffix nitrogen
        // do not compose a descriptor. Nor does a heteroatom inside a prefix:
        // the parent's descriptor cannot reach it (P-82.2.1).
        for smiles in [
            "[13CH3][18O]c1ccccc1",
            "[13CH3][18O][13CH3]",
            "CC([35Cl])CC",
            "C[13CH2]N(C)c1ccccc1",
            "c1ccccc1[13CH2]NN",
            "CC(C)C(=O)O[13CH2]C",
        ] {
            let error = name_smiles(smiles, IsotopeStyle::Iupac, &Default::default()).unwrap_err();
            assert!(error.contains("isotope label"), "{smiles}: {error}");
        }
    }

    #[test]
    fn isotopically_substituted_radicals_cite_their_labels() {
        // Issue #27. P-71.2.1 names the radical like a prefix, and P-82.2.1
        // places the descriptor before it.
        let parser = nazuchem_n2s::Parser::new().unwrap();
        for (smiles, expected) in [
            ("[2H][C]([2H])[2H]", "(2H3)methyl"),
            ("[13CH3]", "(13C)methyl"),
            ("[2H][CH]C", "(1-2H1)ethyl"),
            ("[13CH3][CH]C", "(1-13C)propan-2-yl"),
            ("[13CH3][O]", "(13C)methoxyl"),
            ("[2H]c1cc[c]cc1", "(4-2H1)phenyl"),
        ] {
            let name = name_smiles(smiles, IsotopeStyle::Iupac, &Default::default())
                .unwrap_or_else(|e| panic!("{smiles}: {e}"));
            assert_eq!(name.text, expected, "{smiles}");
            let options = nazuchem_n2s::ParseOptions {
                allow_radicals: true,
                ..Default::default()
            };
            let parsed = parser.parse(&name.text, &options).structure.unwrap();
            crate::test_support::verify_identity(smiles, &parsed.semantic_cxsmiles().unwrap())
                .unwrap();
            assert!(
                name_smiles(smiles, IsotopeStyle::Boughton, &Default::default()).is_err(),
                "{smiles}"
            );
        }
        // Radical names that are not rendered as a prefix cite no label.
        for smiles in ["[13CH3][NH]", "[2H][N][2H]", "[13CH2][13CH2]"] {
            assert!(
                name_smiles(smiles, IsotopeStyle::Iupac, &Default::default()).is_err(),
                "{smiles}"
            );
        }
    }

    #[test]
    fn isotopic_onium_keeps_the_ordinary_parent_plan() {
        // The sulfanium shortcut has no atom locants for the isotope descriptor.
        let smiles = "[13CH3][S+](C)C";
        let name = name_smiles(smiles, IsotopeStyle::Iupac, &Default::default()).unwrap();
        assert_eq!(name.text, "(dimethylsulfaniumyl)(1-13C)methane");
        verify(smiles, &name.text).unwrap();
    }

    #[test]
    fn isotopic_hydrazinecarboxylates_keep_the_formohydrazide_plan() {
        // The contracted hydrazinecarboxylate parent has no carbon locant.
        for (smiles, expected) in [
            ("NN[13C](=O)OC", "1-methoxy(1-13C)formohydrazide"),
            ("NN[13C](=O)O", "1-hydroxy(13C)formohydrazide"),
        ] {
            let name = name_smiles(smiles, IsotopeStyle::Iupac, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            verify(smiles, &name.text).unwrap();
        }
    }

    #[test]
    fn isotopic_spiro_parent_uses_the_engine_locants_and_exact_graph() {
        let smiles = "[13CH2]1CCC2(CC1)CCCC2";
        let name = name_smiles(smiles, IsotopeStyle::Iupac, &Default::default()).unwrap();
        assert!(name.text.contains("13C"));
        verify(smiles, &name.text).unwrap();
        assert!(!matches!(name.status, NameStatus::Preferred));
    }
}
