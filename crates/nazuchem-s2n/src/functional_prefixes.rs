// SPDX-License-Identifier: MIT AND (MIT OR Apache-2.0)
// Origin: mixed
// Source: openclatura | src/openclatura/functional_prefixes.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | src/openclatura/chalcogen_vocabulary.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: ReShiki | src/naming/openclatura/functional_prefixes.rs | 5f559fd482e237307a589dfcb7ea333cee0a6e13
// Source: openclatura | tests/data/parser_xml_resources/simpleSubstituents.json | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Copyright (c) 2026 lamalab-org
// Copyright (c) 2026 Ameyanagi and ReShiki contributors
// Changes: Openclatura characteristic-prefix handlers consolidated into a Rust dispatcher.
// Local contributions: MIT OR Apache-2.0; see NOTICE and PROVENANCE-nazuchem.toml.

//! Characteristic groups cited as prefixes, in upstream registry order.
use super::{
    assembly_parts::{RenderedSubstituentName, SubstituentItem},
    assembly_prefixes::substituent_sort_key,
    component_group_rules::{amide_nitrogen, ester_or_peroxy_single_oxygen, is_carboxylate_group},
    formatting::{format_center_ligands_by, format_multiplier, strip_outer_parentheses},
    molecule::{AtomSet, Molecule, bond_ids_within},
    namer::{name_subgraph_rendered, name_subgraph_with_upstream},
    nomenclature::{ChalcogenElement, NAMER_RULES, RULES},
    perception::{LigandRole, PerceivedGroup},
    subgraph_tools::subgraph_component,
};
use nazuchem_ir::{NameReason, NameStatus};
use std::collections::BTreeMap;

pub struct PrefixContext<'a> {
    pub mol: &'a Molecule,
    pub parent_path: &'a [usize],
    pub sub_exclude: &'a AtomSet,
}

// Provenance: fn chalcogen_term | adapted
// Source: openclatura | src/openclatura/chalcogen_vocabulary.py::_element_data
fn chalcogen_term(
    element: &str,
    field: &str,
    column: fn(&'static ChalcogenElement) -> &'static Option<String>,
) -> Result<&'static str, String> {
    NAMER_RULES
        .chalcogen_elements
        .get(element)
        .and_then(|e| column(e).as_deref())
        .ok_or_else(|| format!("openclatura: {element} {field} not ported"))
}

// Provenance: fn element_ligand | adapted
// Source: openclatura | src/openclatura/functional_prefixes.py::inline in ester_prefix_handler
fn element_ligand(branch: RenderedSubstituentName, element: &str) -> Result<String, String> {
    // Upstream divergence: amidooxy contractions change the attachment of N-acyl aminooxy ligands.
    let plain = strip_outer_parentheses(&branch.text);
    if plain.ends_with("amido") && !super::formatting::is_complex_prefix(&plain) {
        return Err("openclatura: amido chalcogen-ligand attachment is not supported".into());
    }
    if element == "O" {
        return Ok(branch.into_oxy_prefix());
    }
    Ok(super::formatting::format_element_substituent(
        "",
        &branch.text,
        chalcogen_term(element, "substituent_suffix", |e| &e.substituent_suffix)?,
        false,
    ))
}

// Provenance: fn central_linkage_prefix | adapted
// Source: openclatura | src/openclatura/chalcogen_vocabulary.py::central_linkage_prefix
fn central_linkage_prefix(group: &PerceivedGroup) -> Option<&'static str> {
    let descriptor = group.descriptor.as_ref()?;
    let element = descriptor.central_element.as_deref()?;
    let count = descriptor.ligands_with_role(LigandRole::DoubleBonded).len();
    NAMER_RULES
        .chalcogen_central_acid_origins
        .iter()
        .find(|r| r.element == element && r.double_bond_count == count)
        .map(|r| r.linkage_prefix.as_str())
}

// Provenance: fn descriptor_ester_prefix | adapted
// Source: openclatura | src/openclatura/functional_prefixes.py::ester_prefix_handler
// Source: openclatura | src/openclatura/functional_prefixes.py::central_ester_prefix_handler
fn descriptor_ester_prefix(
    context: &PrefixContext<'_>,
    group: &PerceivedGroup,
) -> Result<Option<String>, String> {
    let Some(d) = &group.descriptor else {
        return Ok(None);
    };
    let doubles = d.ligands_with_role(LigandRole::DoubleBonded);
    let linkers: Vec<_> = d
        .ligands
        .iter()
        .filter(|l| l.role != LigandRole::DoubleBonded)
        .collect();
    let Some(terminal) = linkers.last() else {
        return Ok(None);
    };
    let Some(root) = terminal.attachment_atom else {
        return Ok(None);
    };
    if doubles.len() != 1 && d.family != "central_acid" {
        return Ok(None);
    }
    let mut blocked = context.sub_exclude.clone();
    blocked.insert(terminal.atom);
    let branch = name_subgraph_rendered(context.mol, root, &blocked, Some(terminal.atom), false)?;
    let mut ligand = element_ligand(branch, &terminal.element)?;
    for linker in linkers[..linkers.len() - 1].iter().rev() {
        ligand = element_ligand(
            RenderedSubstituentName {
                text: ligand,
                ..Default::default()
            },
            &linker.element,
        )?;
    }
    let base = if d.family == "central_acid" {
        central_linkage_prefix(group).unwrap_or("")
    } else {
        chalcogen_term(&doubles[0].element, "acyl_center_prefix", |e| {
            &e.acyl_center_prefix
        })?
    };
    Ok(Some(format!(
        "({}{base})",
        strip_outer_parentheses(&ligand)
    )))
}

// Provenance: fn central_chalcogen_linkage_prefix | adapted
// Source: openclatura | src/openclatura/functional_prefixes.py::central_chalcogen_linkage_prefix_handler
fn central_chalcogen_linkage_prefix(
    context: &PrefixContext<'_>,
    group: &PerceivedGroup,
) -> Result<String, String> {
    let d = group
        .descriptor
        .as_ref()
        .ok_or("openclatura: missing central linkage descriptor")?;
    let bridge: Vec<_> = d
        .ligands
        .iter()
        .filter(|l| l.role != LigandRole::DoubleBonded)
        .collect();
    let Some(base) = central_linkage_prefix(group) else {
        return Ok(group.prefix().unwrap_or("").into());
    };
    if bridge.len() != 2 {
        return Ok(group.prefix().unwrap_or("").into());
    }
    let terminal = bridge[1];
    let terminal_text = if terminal.role == LigandRole::HydrogenBearing {
        let key = NAMER_RULES
            .chalcogen_simple_group_keys
            .iter()
            .find(|r| {
                r.family == "hydroxy" && r.derivative == "alcohol" && r.element == terminal.element
            })
            .map(|r| r.key.as_str());
        key.and_then(|k| RULES.functional_groups.cited_prefix_for(k))
            .unwrap_or_default()
    } else if matches!(
        terminal.role,
        LigandRole::CarbonLink | LigandRole::HeteroatomLink
    ) {
        if let Some(root) = terminal.attachment_atom {
            let mut blocked = context.sub_exclude.clone();
            blocked.insert(terminal.atom);
            element_ligand(
                name_subgraph_rendered(context.mol, root, &blocked, Some(terminal.atom), false)?,
                &terminal.element,
            )?
        } else {
            String::new()
        }
    } else {
        String::new()
    };
    if terminal_text.is_empty() {
        return Ok(group.prefix().unwrap_or("").into());
    }
    Ok(format!(
        "({}{base})",
        element_ligand(
            RenderedSubstituentName {
                text: terminal_text,
                ..Default::default()
            },
            &bridge[0].element
        )?
    ))
}

// Provenance: fn nitrogen_group_prefix | adapted
// Source: openclatura | src/openclatura/functional_prefixes.py::amidine_prefix_handler
// Source: openclatura | src/openclatura/functional_prefixes.py::hydrazide_prefix_handler
// Source: openclatura | src/openclatura/functional_prefixes.py::hydrazine_prefix_handler
fn nitrogen_group_prefix(
    context: &PrefixContext<'_>,
    group: &PerceivedGroup,
) -> Result<String, String> {
    let mol = context.mol;
    let hydrazine = group.key == "hydrazine";
    let amidine = matches!(group.key.as_str(), "amidine" | "ring_amidine");
    let hydrazide = !hydrazine && !amidine;
    let mut ns: Vec<_> = group
        .atoms_involved
        .iter()
        .copied()
        .filter(|i| mol.atoms[i].symbol == "N")
        .collect();
    // The acyl centre of a hydrazide: carbon, or sulfur in hydrazinesulfonyl.
    let core = group
        .atoms_involved
        .iter()
        .copied()
        .find(|i| {
            mol.atoms[i].symbol != "N"
                && (hydrazide || mol.atoms[i].is_carbon())
                && ns.iter().any(|&n| mol.get_bond(*i, n).is_some())
        })
        .unwrap_or(group.attachment_carbon);
    ns.sort_by_key(|&n| {
        if amidine {
            mol.bond_order(n, core)
        } else {
            usize::from(mol.get_bond(n, core).is_none())
        }
    });
    let stereo = if amidine {
        ns.last()
            .and_then(|&n| mol.get_bond(n, core))
            .and_then(|b| b.stereo.as_ref())
            .filter(|_| !super::formatting::descriptors_omitted())
            .map(|s| format!("({s})-"))
            .unwrap_or_default()
    } else {
        String::new()
    };
    // P-66.3.2.2: ending the parent chain, the group is hydrazinyl with oxo
    // (cited by collect_component_prefix_substituents); its substituents take
    // the locants 1 and 2 of hydrazinyl (P-66.3.2.3, P-68.3.1.2.1).
    let chain_end = hydrazide && context.parent_path.contains(&core);
    // P-66.3.2.1: otherwise it is the acyl prefix of the hydrazine acid, as
    // hydrazinecarbonyl, hydrazinecarbothioyl or hydrazinesulfonyl.
    let acyl = group.prefix().unwrap_or("hydrazinecarbonyl");
    let base = if hydrazine || chain_end {
        "hydrazinyl".into()
    } else if amidine {
        "carbamimidoyl".into()
    } else {
        // P-16.5.1.4: it names the parent hydride hydrazine, so it is
        // enclosed although simple.
        format!("({acyl})")
    };
    let mut names: Vec<(String, Vec<String>)> = Vec::new();
    for (index, &n) in ns.iter().take(2).enumerate() {
        let loc = if hydrazine || amidine {
            if index == 0 { "N" } else { "N'" }
        } else if index == 0 {
            "1"
        } else {
            "2"
        };
        for &root in mol.get_neighbors(n) {
            if root == core
                || group.atoms_involved.contains(&root)
                || mol.atoms[&root].symbol == "H"
            {
                continue;
            }
            let mut blocked = context.sub_exclude.clone();
            blocked.insert(n);
            let text = name_subgraph_with_upstream(mol, root, &blocked, Some(n))?;
            // P-66.3.2.3, P-66.3.3.3: with a second acyl group the preferred
            // prefix may be the hydrazido form of the other acid, as
            // acetohydrazido; that choice is not made here. R-CO-N(NH2)- is
            // always a substituted hydrazinyl (1-acetylhydrazin-1-yl).
            if hydrazide
                && !(chain_end && index == 0)
                && mol.get_neighbors(root).iter().any(|&x| {
                    mol.bond_order(root, x) == 2
                        && matches!(mol.atoms[&x].symbol.as_str(), "O" | "S" | "Se" | "Te")
                })
            {
                mol.mark_name(NameStatus::NotPreferred(NameReason::AminoPrefix));
            }
            if let Some((_, locants)) = names.iter_mut().find(|(name, _)| name == &text) {
                locants.push(loc.into());
            } else {
                names.push((text, vec![loc.into()]));
            }
        }
    }
    if names.is_empty() {
        return Ok(if stereo.is_empty() {
            base
        } else {
            format!("({stereo}{base})")
        });
    }
    if !hydrazine {
        names.sort_by_key(|(name, _)| substituent_sort_key(name));
    }
    let prefixes = names
        .into_iter()
        .map(|(name, locs)| {
            format!(
                "{}-{}",
                locs.join(","),
                format_multiplier(
                    &name,
                    locs.len(),
                    super::formatting::is_complex_prefix(&name)
                )
            )
        })
        .collect::<Vec<_>>()
        .join("-");
    let base = if chain_end {
        "hydrazin-1-yl".into()
    } else if hydrazide {
        acyl.replacen("hydrazine", "hydrazine-1-", 1)
    } else {
        base
    };
    Ok(format!("({stereo}{prefixes}{base})"))
}

// Provenance: fn ester_prefix_from_group | adapted
// Source: openclatura | src/openclatura/functional_prefixes.py::ester_prefix_from_group
pub fn ester_prefix_from_group(
    mol: &Molecule,
    group: &PerceivedGroup,
    exclude: &AtomSet,
    suffix: &str,
) -> Result<String, String> {
    let Some(oxygen) = ester_or_peroxy_single_oxygen(mol, group) else {
        return Ok(String::new());
    };
    let Some(&root) = mol
        .get_neighbors(oxygen)
        .iter()
        .find(|n| !group.atoms_involved.contains(n))
    else {
        return Ok(String::new());
    };
    let mut blocked = exclude.clone();
    blocked.insert(oxygen);
    let branch = name_subgraph_rendered(mol, root, &blocked, Some(oxygen), false)?;
    if RULES.prefixes.peroxy_ester_groups.contains(&group.key) {
        return Ok(format!(
            "({}peroxycarbonyl)",
            format_multiplier(&strip_outer_parentheses(&branch.text), 1, false)
        ));
    }
    Ok(format!("({}{suffix})", branch.into_oxy_prefix()))
}

// Provenance: fn amide_prefix_from_group | adapted
// Source: openclatura | src/openclatura/functional_prefixes.py::amide_prefix_from_group
pub fn amide_prefix_from_group(
    mol: &Molecule,
    group: &PerceivedGroup,
    exclude: &AtomSet,
) -> Result<String, String> {
    let Some(nitrogen) = amide_nitrogen(mol, group) else {
        return Ok(String::new());
    };
    let base = group
        .resolved_rule
        .as_ref()
        .and_then(|r| r.prefix.as_deref())
        .unwrap_or("");
    let mut blocked = exclude.clone();
    blocked.insert(nitrogen);
    let subs: Vec<_> = mol
        .get_neighbors(nitrogen)
        .iter()
        .filter(|n| !group.atoms_involved.contains(n) && mol.atoms[n].symbol != "H")
        .copied()
        .collect();
    if subs.is_empty() {
        return Ok(base.into());
    }
    let names = subs
        .into_iter()
        .map(|n| name_subgraph_with_upstream(mol, n, &blocked, Some(nitrogen)))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(format!(
        "({}{base})",
        format_center_ligands_by(&names, substituent_sort_key)
    ))
}

// Provenance: fn prefix_from_group | adapted
// Source: openclatura | src/openclatura/functional_prefixes.py::prefix_from_group
pub fn prefix_from_group(
    context: &PrefixContext<'_>,
    group: &PerceivedGroup,
) -> Result<String, String> {
    let Some(rule) = group
        .resolved_rule
        .as_ref()
        .or_else(|| RULES.functional_groups.get(&group.key))
    else {
        return Ok(String::new());
    };
    let families = &rule.families;
    // P-72.6.1: an acid anion cited as a prefix remains carboxylato. The
    // ester-like registry family is shared with neutral carboxylate esters.
    if is_carboxylate_group(group) {
        return Ok("carboxylato".into());
    }
    if RULES.functional_groups.get(&group.key).is_none()
        && let Some(d) = &group.descriptor
        && d.family == "central_acid"
    {
        if d.ligands
            .iter()
            .filter(|l| l.role != LigandRole::DoubleBonded)
            .count()
            == 2
        {
            return central_chalcogen_linkage_prefix(context, group);
        }
        if rule.has_capability("ester_like")
            && let Some(prefix) = descriptor_ester_prefix(context, group)?
        {
            return Ok(prefix);
        }
    }
    if families.contains("direct_prefix") {
        return Ok(RULES
            .functional_groups
            .cited_prefix_for(&group.key)
            .unwrap_or_default());
    }
    if rule.has_capability("hydrazide")
        || matches!(group.key.as_str(), "amidine" | "ring_amidine" | "hydrazine")
    {
        return nitrogen_group_prefix(context, group);
    }
    if group.key == "imino_prefix" || group.key == "iminium" {
        let base = if group.key == "iminium" {
            "iminio"
        } else {
            "imino"
        };
        let Some(n) = amide_nitrogen(context.mol, group) else {
            return Ok(base.into());
        };
        let roots: Vec<_> = context
            .mol
            .get_neighbors(n)
            .iter()
            .copied()
            .filter(|r| {
                *r != group.attachment_carbon
                    && !group.atoms_involved.contains(r)
                    && context.mol.atoms[r].symbol != "H"
            })
            .collect();
        if group.key == "iminium" {
            if let Some(prefix) = super::heteroatom_subgraphs::aminium_prefix(
                context.mol,
                n,
                context.sub_exclude,
                Some(group.attachment_carbon),
            )? {
                return Ok(prefix);
            }
            super::amine_oxides::assess_cationic_prefix(context.mol, n);
            if context.mol.get_neighbors(n).iter().all(|&neighbor| {
                neighbor == group.attachment_carbon || context.mol.atoms[&neighbor].symbol == "H"
            }) {
                return Ok("azaniumylidene".into());
            }
            context
                .mol
                .mark_name(NameStatus::Unverified(NameReason::CationPrefix));
        }
        if roots.is_empty() {
            return Ok(base.into());
        }
        let mut blocked = context.sub_exclude.clone();
        blocked.insert(n);
        let names = roots
            .into_iter()
            .map(|r| name_subgraph_with_upstream(context.mol, r, &blocked, Some(n)))
            .collect::<Result<Vec<_>, _>>()?;
        return Ok(format!(
            "({}{base})",
            format_center_ligands_by(&names, substituent_sort_key)
        ));
    }
    if families.contains("sulfonyl") {
        let name = ester_prefix_from_group(context.mol, group, context.sub_exclude, "sulfonyl")?;
        // P-72.6.1: without an ester group the registry prefix applies,
        // `sulfo` for the acid and `sulfonato` for its anion.
        return Ok(if name.is_empty() {
            rule.prefix.clone().unwrap_or_default()
        } else {
            name
        });
    }
    if families.contains("peroxy_acid") {
        // Upstream divergence: hydroperoxycarbonyl duplicates a center already in the parent.
        if group.descriptor.as_ref().is_none_or(|descriptor| {
            descriptor
                .centers
                .iter()
                .any(|center| context.parent_path.contains(center))
        }) {
            return Err("openclatura: peroxy-acid prefix overlaps its parent carbon".into());
        }
        // P-65.1.5.3: a chalcogen peroxy acid is cited from its own atoms,
        // as in `4-[(hydroxysulfanyl)carbonyl]cyclohexane-1-carboxylic acid`.
        if let Some([double, link, terminal]) = group
            .descriptor
            .as_ref()
            .map(|descriptor| descriptor.ligands.as_slice())
            .filter(|ligands| ligands.iter().any(|ligand| ligand.element != "O"))
        {
            let hydride = NAMER_RULES
                .chalcogen_simple_group_keys
                .iter()
                .find(|r| {
                    r.family == "hydroxy"
                        && r.derivative == "alcohol"
                        && r.element == terminal.element
                })
                .and_then(|r| RULES.functional_groups.cited_prefix_for(&r.key))
                .ok_or("openclatura: peroxy-acid prefix has no terminal chalcogen name")?;
            // P-63.4.2.2: -OOH is `hydroperoxy`, -SSH `disulfanyl`, -OSH
            // `sulfanyloxy`.
            let ligand = if link.element == "O" && terminal.element == "O" {
                "hydroperoxy".into()
            } else if link.element == terminal.element {
                format!("di{hydride}")
            } else {
                element_ligand(
                    RenderedSubstituentName {
                        text: hydride,
                        ..Default::default()
                    },
                    &link.element,
                )?
            };
            return Ok(format!(
                "({ligand}{})",
                chalcogen_term(&double.element, "acyl_center_prefix", |e| {
                    &e.acyl_center_prefix
                })?
            ));
        }
        return Ok("hydroperoxycarbonyl".into());
    }
    if rule.has_capability("acyl_leaving_group") {
        return Ok(rule
            .prefix
            .as_ref()
            .map(|s| format!("({s})"))
            .unwrap_or_default());
    }
    if families.contains("cyano_prefix") {
        return Ok("cyano".into());
    }
    if families.contains("carboxy_prefix") {
        // P-65.2.1.6: a chalcogen analogue of carboxy names its own atoms,
        // `sulfanylcarbonyl` (HS-CO-) and `hydroxycarbonothioyl` (HO-CS-);
        // `dithiocarboxy` (HS-CS-) needs no position.
        if let Some([double, single]) = group
            .descriptor
            .as_ref()
            .filter(|d| d.derivative == "acid")
            .map(|d| d.ligands.as_slice())
            && (double.element != "O" || single.element != "O")
        {
            if double.element == single.element {
                return Ok(format!("({})", rule.prefix.as_deref().unwrap_or_default()));
            }
            let hydride = if single.element == "O" {
                "hydroxy"
            } else {
                chalcogen_term(&single.element, "substituent_suffix", |e| {
                    &e.substituent_suffix
                })?
            };
            return Ok(format!(
                "({hydride}{})",
                chalcogen_term(&double.element, "acyl_center_prefix", |e| {
                    &e.acyl_center_prefix
                })?
            ));
        }
        // `carboxy` would drop the charge of a chalcogen analogue of
        // carboxylato (all-oxygen carboxylato is handled above).
        if group
            .descriptor
            .as_ref()
            .is_some_and(|d| d.derivative == "anion")
        {
            return Err("openclatura: chalcogen carboxylate anion prefix not ported".into());
        }
        return Ok("carboxy".into());
    }
    if families.contains("amide_like") || rule.has_capability("amide_like") {
        return amide_prefix_from_group(context.mol, group, context.sub_exclude);
    }
    if families.contains("ester_like") || rule.has_capability("ester_like") {
        if !RULES.prefixes.peroxy_ester_groups.contains(&group.key)
            && let Some(prefix) = descriptor_ester_prefix(context, group)?
        {
            return Ok(prefix);
        }
        return ester_prefix_from_group(context.mol, group, context.sub_exclude, "carbonyl");
    }
    Ok(if context.parent_path.contains(&group.attachment_carbon) {
        if group.key == "aminium"
            && !amide_nitrogen(context.mol, group).is_some_and(|n| {
                context.mol.get_neighbors(n).iter().all(|&neighbor| {
                    neighbor == group.attachment_carbon
                        || context.mol.atoms[&neighbor].symbol == "H"
                })
            })
        {
            if let Some(n) = amide_nitrogen(context.mol, group) {
                super::amine_oxides::assess_cationic_prefix(context.mol, n);
            }
            context
                .mol
                .mark_name(NameStatus::Unverified(NameReason::CationPrefix));
            return Ok("ammonio".into());
        }
        rule.prefix.clone().unwrap_or_default()
    } else {
        String::new()
    })
}

// Provenance: fn collect_component_prefix_substituents | adapted
// Source: openclatura | src/openclatura/functional_prefixes.py::collect_component_prefix_substituents
pub fn collect_component_prefix_substituents(
    mol: &Molecule,
    groups: &[PerceivedGroup],
    path: &[usize],
    exclude: &AtomSet,
) -> Result<(BTreeMap<usize, Vec<SubstituentItem>>, AtomSet), String> {
    let parent: AtomSet = path.iter().copied().collect();
    let mut mapping: BTreeMap<usize, Vec<SubstituentItem>> = BTreeMap::new();
    let mut handled = AtomSet::new();
    let context = PrefixContext {
        mol,
        parent_path: path,
        sub_exclude: exclude,
    };
    for group in groups {
        if group
            .resolved_rule
            .as_ref()
            .is_some_and(|r| r.families.contains("prefix_skip"))
            || !parent.contains(&group.attachment_carbon)
            // P-66.1.1.4.1.1: an amide carbonyl carbon ending a parent chain
            // is cited as amino with oxo, both named as ordinary substituents
            // (4-amino-4-oxobutanoic acid). An acyl carbon that is the whole
            // parent has no carbon left for a carbamoyl or carbonyl prefix,
            // so its ligands are cited the same way; a hydrazide is handled
            // below. P-65.1.4.2 (2), P-65.1.5.1, P-65.1.5.3: a peroxy acid
            // or a mixed chalcogen acid ending the chain is cited the same
            // way (6-hydroperoxy-6-oxohexanoic acid,
            // 4-oxo-4-sulfanylbutanoic acid).
            || group.descriptor.as_ref().is_some_and(|d| {
                d.family == "acyl"
                    && (d.derivative == "amide"
                        || d.names_its_chalcogens()
                        || path.len() == 1 && d.derivative != "hydrazide")
                    && d.centers.contains(&group.attachment_carbon)
            })
        {
            continue;
        }
        // P-66.5.1.1.4: cyano is the whole -CN group as a substituent, so a
        // nitrile carbon that is itself the whole parent cannot be cited
        // that way; cyanomethane would count the carbon twice.
        if group.key == "nitrile" && path.len() == 1 && group.atoms_involved.contains(&path[0]) {
            return Err(
                "openclatura: a nitrile carbon that is the whole parent has no cyano prefix".into(),
            );
        }
        let name = prefix_from_group(&context, group)?;
        if name.is_empty() {
            continue;
        }
        let mut atoms = group.atoms_involved.clone();
        for &idx in &group.atoms_involved {
            for &n in mol.get_neighbors(idx) {
                if n != group.attachment_carbon
                    && !group.atoms_involved.contains(&n)
                    && !parent.contains(&n)
                    && !exclude.contains(&n)
                    && mol.atoms[&n].symbol != "H"
                {
                    let mut blocked = exclude.clone();
                    blocked.insert(idx);
                    atoms.extend(subgraph_component(mol, n, &blocked));
                }
            }
        }
        let mut bond_atoms = atoms.clone();
        bond_atoms.insert(group.attachment_carbon);
        let mut item = SubstituentItem::new(name, vec![]);
        // Cyano is a simple characteristic prefix even when its parent
        // locant is omitted and another substituent is cited beside it.
        item.outer_parentheses_optional = matches!(group.key.as_str(), "nitrile" | "ring_nitrile");
        item.bond_ids = bond_ids_within(mol, &bond_atoms);
        item.charge_atom_ids = atoms
            .iter()
            .copied()
            .filter(|i| mol.atoms[i].charge != 0)
            .collect();
        item.atom_ids = atoms;
        mapping
            .entry(group.attachment_carbon)
            .or_default()
            .push(item);
        handled.extend(&group.atoms_involved);
        // P-66.3.2.2: a hydrazide carbonyl carbon ending the parent chain
        // cites its chalcogen beside hydrazinyl, as in
        // 3-hydrazinyl-3-oxopropanoic acid and
        // 2-hydrazinyl-2-sulfanylideneacetohydrazide (P-66.3.4).
        if let Some(double) = group
            .descriptor
            .as_ref()
            .filter(|d| d.derivative == "hydrazide" && d.centers.contains(&group.attachment_carbon))
            .and_then(|d| {
                d.ligands
                    .iter()
                    .find(|l| l.role == LigandRole::DoubleBonded)
            })
        {
            let name = if double.element == "O" {
                "oxo".into()
            } else {
                format!(
                    "{}idene",
                    chalcogen_term(&double.element, "substituent_suffix", |e| {
                        &e.substituent_suffix
                    })?
                )
            };
            let mut item = SubstituentItem::new(name, vec![]);
            item.atom_ids = AtomSet::from([double.atom]);
            item.bond_ids =
                bond_ids_within(mol, &AtomSet::from([group.attachment_carbon, double.atom]));
            mapping
                .entry(group.attachment_carbon)
                .or_default()
                .push(item);
        }
    }
    Ok((mapping, handled))
}

#[cfg(test)]
// Provenance: mod tests | original
// Source: ReShiki | src/naming/openclatura/functional_prefixes.rs::tests
mod tests {
    use super::super::name_smiles;

    #[test]
    // Provenance: fn tests::acid_anion_prefix_uses_its_charge_descriptor_not_the_shared_ester_key | adapted
    // Source: ReShiki | src/naming/openclatura/functional_prefixes.rs::tests::acid_anion_prefix_uses_its_charge_descriptor_not_the_shared_ester_key
    // Source: openclatura | tests/data/parser_xml_resources/simpleSubstituents.json::<fixture data>
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn acid_anion_prefix_uses_its_charge_descriptor_not_the_shared_ester_key() {
        for (smiles, derivative, expected) in [
            ("CC(=O)[O-]", "anion", "carboxylato"),
            ("COC(=O)c1ccccc1", "ester", "(methoxycarbonyl)"),
        ] {
            let mol = super::super::graph_io::read_smiles(smiles, &Default::default()).unwrap();
            let group = super::super::perception::perceive_groups(&mol)
                .into_iter()
                .find(|g| {
                    g.descriptor
                        .as_ref()
                        .is_some_and(|d| d.family == "acyl" && d.derivative == derivative)
                })
                .unwrap();
            let context = super::PrefixContext {
                mol: &mol,
                parent_path: &[group.attachment_carbon],
                sub_exclude: &super::AtomSet::from([group.attachment_carbon]),
            };
            assert_eq!(
                super::prefix_from_group(&context, &group).unwrap(),
                expected
            );
        }
    }

    /// Issue #148: an amide cited as a prefix.
    #[test]
    fn amide_prefixes_follow_p_66_1_1_4_1() {
        use nazuchem_ir::NameStatus;

        for (smiles, expected) in [
            // P-66.1.1.4.1.1 (1), P-65.1.6.1: -CO-NH2 ending a chain of two
            // or more carbons is amino + oxo (4-amino-4-oxobutanoic acid,
            // 3-amino-3-oxopropanoic acid).
            ("OC(=O)CCC(N)=O", "4-amino-4-oxobutanoic acid"),
            ("OC(=O)CC(N)=O", "3-amino-3-oxopropanoic acid"),
            ("NC(=O)CCC(=O)OC", "methyl 4-amino-4-oxobutanoate"),
            ("OS(=O)(=O)CC(N)=O", "2-amino-2-oxoethane-1-sulfonic acid"),
            // P-103.1.1.1: asparagine and glutamine.
            ("NC(=O)CC(N)C(O)=O", "2,4-diamino-4-oxobutanoic acid"),
            ("NC(=O)CCC(N)C(O)=O", "2,5-diamino-5-oxopentanoic acid"),
            // P-66.1.1.4.1.2: -CO-CO-NH2 (4-amino-3,4-dioxobutanoic acid).
            ("OC(=O)CC(=O)C(N)=O", "4-amino-3,4-dioxobutanoic acid"),
            // P-66.1.4.4: chalcogen analogues
            // (3-amino-3-sulfanylidenepropanoic acid).
            ("OC(=O)CC(N)=S", "3-amino-3-sulfanylidenepropanoic acid"),
            ("OC(=O)CC(N)=[Se]", "3-amino-3-selanylidenepropanoic acid"),
            ("NC(=O)CCC(N)=S", "4-amino-4-sulfanylidenebutanamide"),
            // P-66.1.1.4.1.1, P-65.1.6.2: the amino group is substituted in
            // the normal way (5-anilino-5-oxopentanoic acid,
            // anilino(oxo)acetic acid).
            ("OC(=O)CCCC(=O)Nc1ccccc1", "5-anilino-5-oxopentanoic acid"),
            ("OC(=O)C(=O)Nc1ccccc1", "anilino(oxo)acetic acid"),
            ("OC(=O)C(=O)NC1=CC=CC=C1", "anilino(oxo)acetic acid"),
            ("OC(=O)CCC(=O)NC", "4-(methylamino)-4-oxobutanoic acid"),
            ("OC(=O)CCC(=O)N(C)C", "4-(dimethylamino)-4-oxobutanoic acid"),
            // P-45.2.1: of two chains of equal length, the one through the
            // amide carbon has more substituents cited as prefixes, as in
            // 2-methyl-3-(methylamino)-3-oxopropyl.
            ("OC(=O)CC(C)C(N)=O", "4-amino-3-methyl-4-oxobutanoic acid"),
            (
                "OC(=O)c1ccc(cc1)CC(C)C(N)=O",
                "4-(3-amino-2-methyl-3-oxopropyl)benzoic acid",
            ),
            (
                "OC(=O)CC(CC(N)=O)CC(N)=O",
                "5-amino-3-(2-amino-2-oxoethyl)-5-oxopentanoic acid",
            ),
            // P-66.1.1.4.1.1 (2): carbamoyl on a ring and on a nonterminal
            // chain atom (2-carbamoylbenzoic acid,
            // 3-(dimethylcarbamoyl)pentanedioic acid).
            ("OC(=O)c1ccccc1C(N)=O", "2-carbamoylbenzoic acid"),
            (
                "OC(=O)CC(CC(O)=O)C(=O)N(C)C",
                "3-(dimethylcarbamoyl)pentanedioic acid",
            ),
            ("OC(=O)CC(CC)C(N)=O", "3-carbamoylpentanoic acid"),
            (
                "OC(=O)CC(C(N)=O)C(N)=O",
                "4-amino-3-carbamoyl-4-oxobutanoic acid",
            ),
            (
                "OC(=O)c1ccc(o1)CC(N)=O",
                "5-(2-amino-2-oxoethyl)furan-2-carboxylic acid",
            ),
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, NameStatus::Preferred, "{smiles}");
            let parsed = crate::test_support::parse_strict(expected).unwrap();
            assert_eq!(
                crate::test_support::canonical_smiles(&parsed).unwrap(),
                crate::test_support::canonical_smiles(smiles).unwrap(),
                "{expected}"
            );
        }
        // P-65.1.1.1, P-65.1.6.1: oxamic acid is retained, may be
        // substituted and is not generated; its systematic names are not
        // labelled preferred.
        for (smiles, expected) in [
            ("OC(=O)C(N)=O", "amino(oxo)acetic acid"),
            ("COC(=O)C(N)=O", "methyl amino(oxo)acetate"),
            ("OC(=O)C(=O)NC", "(methylamino)(oxo)acetic acid"),
            ("OC(=O)C(=O)N(C)C", "(dimethylamino)(oxo)acetic acid"),
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(
                name.status,
                NameStatus::NotPreferred(nazuchem_ir::NameReason::UnretainedParent),
                "{smiles}"
            );
            let parsed = crate::test_support::parse_strict(expected).unwrap();
            assert_eq!(
                crate::test_support::canonical_smiles(&parsed).unwrap(),
                crate::test_support::canonical_smiles(smiles).unwrap(),
                "{expected}"
            );
        }
    }

    /// Issue #153: an amide, ester or acid halide carbon that is the whole
    /// parent of a sulfur acid.
    #[test]
    fn one_carbon_acyl_parents_cite_their_ligands() {
        use nazuchem_ir::{NameReason, NameStatus};

        for (smiles, expected) in [
            // P-65.3.1: the sulfonic or sulfinic acid suffix sits on methane,
            // whose carbon is the amide carbon, so carbamoyl (P-65.2.1.5)
            // would cite that carbon twice. Amino with oxo is the remaining
            // substitutive form; P-66.1.1.4.1.1 (1) states it only for chains
            // of two or more carbons, so the names are not labelled preferred.
            ("NC(=O)S(O)(=O)=O", "amino(oxo)methanesulfonic acid"),
            ("NC(=O)S(=O)(=O)[O-]", "amino(oxo)methanesulfonate"),
            ("NC(=O)S(=O)(=O)OC", "methyl amino(oxo)methanesulfonate"),
            ("NC(=O)S(Cl)(=O)=O", "amino(oxo)methanesulfonyl chloride"),
            ("NC(=O)S(O)=O", "amino(oxo)methanesulfinic acid"),
            ("NC(=O)[Se](O)(=O)=O", "amino(oxo)methaneselenonic acid"),
            // P-65.3.1.1: functional replacement in the suffix.
            ("NC(=O)S(O)(=O)=S", "amino(oxo)methanesulfonothioic O-acid"),
            // P-66.1.4.4: chalcogen analogues of the amide.
            (
                "NC(=S)S(O)(=O)=O",
                "amino(sulfanylidene)methanesulfonic acid",
            ),
            // P-66.1.1.4.1.1: the amino group is substituted in the normal
            // way.
            (
                "CNC(=O)S(O)(=O)=O",
                "(methylamino)(oxo)methanesulfonic acid",
            ),
            (
                "CN(C)C(=O)S(O)(=O)=O",
                "(dimethylamino)(oxo)methanesulfonic acid",
            ),
            (
                "O=C(Nc1ccccc1)S(O)(=O)=O",
                "anilino(oxo)methanesulfonic acid",
            ),
            // Esters and acid halides share the cause: alkoxy or halo with
            // oxo (P-65.6.3.2.3, P-65.5.4 (3)). Neither rule says whether a
            // single carbon counts as a chain end, so these are not labelled
            // preferred either.
            ("COC(=O)S(O)(=O)=O", "methoxy(oxo)methanesulfonic acid"),
            ("ClC(=O)S(O)(=O)=O", "chloro(oxo)methanesulfonic acid"),
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(
                name.status,
                NameStatus::Unverified(NameReason::UnauditedParent),
                "{smiles}"
            );
            let parsed = crate::test_support::parse_strict(expected).unwrap();
            assert_eq!(
                crate::test_support::canonical_smiles(&parsed).unwrap(),
                crate::test_support::canonical_smiles(smiles).unwrap(),
                "{expected}"
            );
        }
        // Unchanged neighbours: a chain of two carbons (P-66.1.1.4.1.1 (1))
        // and a formamide parent.
        for (smiles, expected) in [
            ("NC(=O)CS(O)(=O)=O", "2-amino-2-oxoethane-1-sulfonic acid"),
            ("NC(=O)S(N)(=O)=O", "1-sulfamoylformamide"),
            ("NC(=O)S(=O)(=O)c1ccccc1", "1-(benzenesulfonyl)formamide"),
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
        }
    }

    /// Issue #141: a hydrazide cited as a prefix.
    #[test]
    fn hydrazide_prefixes_follow_p_66_3_2() {
        use nazuchem_ir::NameStatus;

        for (smiles, expected) in [
            // P-66.3.2.2: -CO-NHNH2 ending a chain is hydrazinyl + oxo
            // (3-hydrazinyl-3-oxopropanoic acid); in a substituent chain as
            // the amide in 2-amino-2-oxoethyl (P-66.1.1.4.1.1).
            ("OC(=O)CC(=O)NN", "3-hydrazinyl-3-oxopropanoic acid"),
            ("OC(=O)CCC(=O)NN", "4-hydrazinyl-4-oxobutanoic acid"),
            ("NC(=N)CCC(=O)NN", "4-hydrazinyl-4-oxobutanimidamide"),
            ("NNC(=O)CCC(N)=O", "4-hydrazinyl-4-oxobutanamide"),
            (
                "OC(=O)c1ccc(cc1)CCC(=O)NN",
                "4-(3-hydrazinyl-3-oxopropyl)benzoic acid",
            ),
            ("NC(N)=NCCC(=O)NN", "2-(3-hydrazinyl-3-oxopropyl)guanidine"),
            (
                "OC(=O)CC(CC(=O)NN)CC(=O)NN",
                "5-hydrazinyl-3-(2-hydrazinyl-2-oxoethyl)-5-oxopentanoic acid",
            ),
            // P-45.2.1: of two chains of equal length, the one through the
            // hydrazide carbon has more substituents cited as prefixes.
            (
                "OC(=O)CC(C)C(=O)NN",
                "4-hydrazinyl-3-methyl-4-oxobutanoic acid",
            ),
            (
                "OC(=O)c1ccc(cc1)CC(C)C(=O)NN",
                "4-(3-hydrazinyl-2-methyl-3-oxopropyl)benzoic acid",
            ),
            // P-66.3.4: 2-hydrazinyl-2-sulfanylideneacetohydrazide.
            (
                "NNC(=O)C(=S)NN",
                "2-hydrazinyl-2-sulfanylideneacetohydrazide",
            ),
            (
                "OC(=O)CCC(=S)NN",
                "4-hydrazinyl-4-sulfanylidenebutanoic acid",
            ),
            // P-66.3.2.3 (2), P-66.3.6: a substituted hydrazinyl takes the
            // locants 1 and 2, as in 2-(benzylidenehydrazinyl)-2-oxoacetamide.
            (
                "OC(=O)CCC(=O)NNC",
                "4-(2-methylhydrazin-1-yl)-4-oxobutanoic acid",
            ),
            (
                "OC(=O)CCC(=O)N(C)N",
                "4-(1-methylhydrazin-1-yl)-4-oxobutanoic acid",
            ),
            // P-66.3.2.3: R-CO-N(NH2)- is a substituted hydrazinyl, as in
            // 4-(1-acetyl-2-ethylhydrazin-1-yl)benzoic acid.
            (
                "OC(=O)CCC(=O)N(N)C(C)=O",
                "4-(1-acetylhydrazin-1-yl)-4-oxobutanoic acid",
            ),
            // P-66.3.2.1 (1), P-16.5.1.4: on a ring, a heteroatom or a
            // nonterminal chain atom the acyl prefix of the hydrazine acid,
            // enclosed because it names the parent hydride hydrazine, as in
            // 2-(hydrazinecarbonyl)benzene-1-sulfonic acid,
            // (hydrazinesulfonyl)acetic acid and
            // 3-(hydrazinesulfinyl)naphthalene-2-carboxylic acid.
            (
                "OS(=O)(=O)c1ccccc1C(=O)NN",
                "2-(hydrazinecarbonyl)benzene-1-sulfonic acid",
            ),
            ("OC(=O)CS(=O)(=O)NN", "(hydrazinesulfonyl)acetic acid"),
            (
                "OC(=O)c1cc2ccccc2cc1S(=O)NN",
                "3-(hydrazinesulfinyl)naphthalene-2-carboxylic acid",
            ),
            (
                "NC(=N)c1ccc(cc1)C(=O)NN",
                "4-(hydrazinecarbonyl)benzene-1-carboximidamide",
            ),
            (
                "OC(=O)c1ccc(cc1)C(=S)NN",
                "4-(hydrazinecarbothioyl)benzoic acid",
            ),
            (
                "OC(=O)CC(CC(=O)O)C(=O)NN",
                "3-(hydrazinecarbonyl)pentanedioic acid",
            ),
            (
                "OC(=O)c1ccc(cc1)CC(C(=O)NN)CC",
                "4-[2-(hydrazinecarbonyl)butyl]benzoic acid",
            ),
            (
                "OC(=O)c1cc(C(=O)NN)cc(C(=O)NN)c1",
                "3,5-bis(hydrazinecarbonyl)benzoic acid",
            ),
            // Substituted acyl prefixes keep every substituent and its locant
            // (2-methylhydrazine-1-carboxylic acid, P-68.3.1.2.1).
            (
                "OC(=O)c1ccc(cc1)C(=O)NNC",
                "4-(2-methylhydrazine-1-carbonyl)benzoic acid",
            ),
            (
                "OC(=O)c1ccc(cc1)C(=S)N(C)N",
                "4-(1-methylhydrazine-1-carbothioyl)benzoic acid",
            ),
            (
                "OC(=O)c1ccc(cc1)S(=O)(=O)NNC",
                "4-(2-methylhydrazine-1-sulfonyl)benzoic acid",
            ),
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            assert_eq!(name.status, NameStatus::Preferred, "{smiles}");
            let parsed = crate::test_support::parse_strict(expected).unwrap();
            assert_eq!(
                crate::test_support::canonical_smiles(&parsed).unwrap(),
                crate::test_support::canonical_smiles(smiles).unwrap(),
                "{expected}"
            );
        }
        // Not labelled preferred: the hydrazinecarboxylic ester and amide
        // prefixes (P-68.3.1.2) are not generated, and a diacylhydrazine
        // prefix is not chosen between P-66.3.2.3 (1) and P-66.3.3.3.
        for smiles in [
            "OC(=O)CCOC(=O)NN",
            "OC(=O)CCNC(=O)NN",
            "OC(=O)CCC(=O)NNC(C)=O",
            "OC(=O)c1ccc(cc1)C(=O)NNC(C)=O",
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_ne!(
                name.status,
                NameStatus::Preferred,
                "{smiles}: {}",
                name.text
            );
            let parsed = crate::test_support::parse_strict(&name.text).unwrap();
            assert_eq!(
                crate::test_support::canonical_smiles(&parsed).unwrap(),
                crate::test_support::canonical_smiles(smiles).unwrap(),
                "{}",
                name.text
            );
        }
    }

    #[test]
    fn external_acid_groups_are_cited_before_recursive_parent_assembly() {
        for (smiles, expected) in [
            (
                "[O-]C(=O)c1ccccc1-c1[n-]cc(C(=O)[O-])c1",
                "5-(2-carboxylatophenyl)-1H-pyrrol-1-ide-3-carboxylate",
            ),
            (
                "[O-]C(=O)c1ccccc1CC(=O)[O-]",
                "2-(carboxylatomethyl)benzoate",
            ),
            (
                "O=C([O-])c1ccc(cc1)Oc1ccc(cc1)C(=O)[O-]",
                "4-(4-carboxylatophenoxy)benzoate",
            ),
            ("[O-]C(=O)c1ccccc1CC(=O)O", "2-(carboxymethyl)benzoate"),
            (
                "[O-]C(=O)c1cc([NH3+])ccc1CC(=O)O",
                "5-azaniumyl-2-(carboxymethyl)benzoate",
            ),
        ] {
            // Exercise the same canonical input boundary as the CLI.
            let canonical = crate::test_support::canonical_smiles(smiles).unwrap();
            assert_eq!(
                name_smiles(&canonical, &Default::default()).unwrap().text,
                expected,
                "{smiles}"
            );
        }
    }

    #[test]
    fn unsubstituted_aminium_prefix_uses_azaniumyl() {
        let mol = super::super::graph_io::read_smiles("C[NH3+]", &Default::default()).unwrap();
        let group = super::super::perception::perceive_groups(&mol)
            .into_iter()
            .find(|group| group.key == "aminium")
            .unwrap();
        let context = super::PrefixContext {
            mol: &mol,
            parent_path: &[group.attachment_carbon],
            sub_exclude: &super::AtomSet::new(),
        };
        assert_eq!(
            super::prefix_from_group(&context, &group).unwrap(),
            "azaniumyl"
        );
        assert_eq!(
            mol.name_status(),
            nazuchem_ir::NameStatus::Unverified(nazuchem_ir::NameReason::UnauditedPath)
        );
    }

    #[test]
    fn peroxy_ester_prefix_encloses_only_complex_ligands() {
        for (smiles, expected) in [
            ("COOC(=O)CC(=O)O", "(methylperoxycarbonyl)acetic acid"),
            (
                "CC(C)OOC(=O)CC(=O)O",
                "[(propan-2-yl)peroxycarbonyl]acetic acid",
            ),
        ] {
            let name = name_smiles(smiles, &Default::default()).unwrap();
            assert_eq!(name.text, expected, "{smiles}");
            let parsed = crate::test_support::parse_strict(&name.text).unwrap();
            assert_eq!(
                crate::test_support::canonical_smiles(smiles).unwrap(),
                crate::test_support::canonical_smiles(&parsed).unwrap(),
                "{}",
                name.text
            );
        }
    }

    #[test]
    fn peroxy_carbamate_keeps_competing_ester_prefix_route() {
        assert_eq!(
            name_smiles("COOC(=O)N1[C@H](C=O)CC[C@H]1C(=O)O", &Default::default())
                .map(|name| name.text)
                .as_deref(),
            Ok("(2S,5S)-5-formyl-1-[(methoxyoxy)carbonyl]pyrrolidine-2-carboxylic acid")
        );
    }
}
