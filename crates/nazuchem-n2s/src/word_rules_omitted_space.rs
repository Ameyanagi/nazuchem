// SPDX-License-Identifier: MIT
// Origin: adapted
// Source: OPSIN | opsin-core/src/main/java/uk/ac/cam/ch/wwmm/opsin/WordRulesOmittedSpaceCorrector.java | 42e87ac6efb1bd9f578986f5ea7ee1faa8e93c91
// Source: OPSIN | opsin-core/src/main/java/uk/ac/cam/ch/wwmm/opsin/WordRulesOmittedSpaceCorrector.java | b91b610af5ab07560fedb20730d7aef46bb2bca0
// Copyright Daniel Lowe and OPSIN contributors; MIT (see LICENSE).
// Copyright 2017 Daniel Lowe
// Changes: Preserve signed integer count ranges and check multiplied valency arithmetic.

//! Post-component omitted-space correction from `WordRulesOmittedSpaceCorrector`.
//! OPSIN 2.9.0, b91b610af5ab07560fedb20730d7aef46bb2bca0.
//! Copyright Daniel Lowe and OPSIN contributors; MIT (see LICENSE).
//!
//! This pass runs after group fragments, locants and multipliers have been
//! assigned. Its context supplies graph/stereo facts; it never infers these
//! facts from token text or invokes a Java runtime.

use crate::{
    ParsingError,
    graph::{Element as ChemEl, FragmentId},
    parse_tree::{Arena, NodeId, ParseTree},
};
use std::collections::HashSet;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OutAtomFacts {
    pub element: ChemEl,
    pub valency: u32,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FragmentFacts {
    pub out_atoms: Vec<OutAtomFacts>,
    pub functional_atom_count: usize,
    pub has_default_in_atom: bool,
    /// One environment per substitutable hydrogen, in fragment atom order.
    /// Characteristic atoms are excluded; counts use determineValency(true)
    /// minus incoming, spare and outgoing valency exactly as upstream.
    pub substitutable_hydrogen_environments: Vec<String>,
}

pub trait OmittedSpaceContext {
    fn fragment_facts(&self, fragment: FragmentId) -> Result<FragmentFacts, ParsingError>;
    /// Equivalent to FragmentManager.cloneElement: clone token fragments as
    /// well as tree elements, retaining the build state's fragment mappings.
    fn clone_element(&mut self, arena: &mut Arena, element: NodeId)
    -> Result<NodeId, ParsingError>;
}

pub fn correct_omitted_spaces(
    tree: &mut ParseTree,
    context: &mut impl OmittedSpaceContext,
) -> Result<(), ParsingError> {
    let rules = tree.arena.descendants_named(tree.root, "wordRule");
    let name = tree.arena[tree.root]
        .attribute("name")
        .unwrap_or("")
        .to_owned();
    for rule in rules {
        match tree.arena[rule].attribute("wordRule") {
            Some("divalentFunctionalGroup") => correct_divalent(&mut tree.arena, rule, context)?,
            Some("simple") => correct_ester(&mut tree.arena, rule, &name, context)?,
            _ => {}
        }
    }
    Ok(())
}

fn error(message: impl Into<String>) -> ParsingError {
    ParsingError(message.into())
}

pub fn rightmost_group(arena: &Arena, mut element: NodeId) -> Option<NodeId> {
    while arena[element].name == "bracket" {
        element =
            arena[element].children.iter().rev().copied().find(|&id| {
                matches!(arena[id].name.as_str(), "bracket" | "substituent" | "root")
            })?;
    }
    arena[element]
        .children
        .iter()
        .rev()
        .copied()
        .find(|&id| arena[id].name == "group")
}

fn facts(
    arena: &Arena,
    element: NodeId,
    context: &impl OmittedSpaceContext,
) -> Result<FragmentFacts, ParsingError> {
    let group = rightmost_group(arena, element).ok_or_else(|| {
        error("OPSIN bug: Unable to find rightmost group for omitted-space correction")
    })?;
    group_facts(arena, group, context)
}
fn group_facts(
    arena: &Arena,
    group: NodeId,
    context: &impl OmittedSpaceContext,
) -> Result<FragmentFacts, ParsingError> {
    context.fragment_facts(
        arena[group].fragment.ok_or_else(|| {
            error("OPSIN bug: Group has no fragment during omitted-space correction")
        })?,
    )
}

fn single_carbon_or_silicon_radical(facts: &FragmentFacts) -> bool {
    facts.out_atoms.len() == 1
        && facts.out_atoms[0].valency == 1
        && matches!(facts.out_atoms[0].element, ChemEl::C | ChemEl::Si)
}

fn relevant_children(arena: &Arena, element: NodeId, include_root: bool) -> Vec<NodeId> {
    arena[element]
        .children
        .iter()
        .copied()
        .filter(|&id| {
            matches!(arena[id].name.as_str(), "substituent" | "bracket")
                || include_root && arena[id].name == "root"
        })
        .collect()
}

fn correct_divalent(
    arena: &mut Arena,
    rule: NodeId,
    context: &impl OmittedSpaceContext,
) -> Result<(), ParsingError> {
    let words: Vec<_> = arena
        .children_named(rule, "word")
        .into_iter()
        .filter(|&id| arena[id].attribute("type") == Some("substituent"))
        .collect();
    if words.len() == 1 {
        let children = relevant_children(arena, words[0], false);
        if children.len() == 2 {
            let first = children[0];
            if arena[first].attribute("locant").is_none()
                && arena[first].attribute("multiplier").is_none()
                && single_carbon_or_silicon_radical(&facts(arena, first, context)?)
            {
                let moved = children[1];
                arena.detach(moved);
                let new_word = arena.grouping("word");
                arena[new_word].add_attribute("type", "substituent");
                arena.add_child(new_word, moved);
                arena.insert_after(words[0], new_word);
            }
        }
    }
    Ok(())
}

fn multiplier(arena: &Arena, element: NodeId) -> Result<usize, ParsingError> {
    // Changes: Keep user-derived counts within the source signed-int range.
    arena[element]
        .attribute("multiplier")
        .map_or(Ok(1), |value| {
            value
                .parse::<i32>()
                .ok()
                .and_then(|value| usize::try_from(value).ok())
                .ok_or_else(|| error(format!("OPSIN bug: Invalid multiplier {value}")))
        })
}

fn suitable_for_ester(
    arena: &Arena,
    element: NodeId,
    functional_count: usize,
    context: &impl OmittedSpaceContext,
) -> Result<bool, ParsingError> {
    if arena[element].attribute("locant").is_some() {
        return Ok(false);
    }
    if !single_carbon_or_silicon_radical(&facts(arena, element, context)?) {
        return Ok(false);
    }
    Ok(multiplier(arena, element)? <= functional_count)
}

fn ate_or_ite_ending(value: &str) -> bool {
    let lower = value.to_ascii_lowercase();
    let lower = lower.trim_end_matches([']', ')', '}']);
    ["at", "ate", "it", "ite"]
        .iter()
        .any(|ending| lower.ends_with(ending))
}

fn correct_ester(
    arena: &mut Arena,
    rule: NodeId,
    name: &str,
    context: &mut impl OmittedSpaceContext,
) -> Result<(), ParsingError> {
    let words = arena.children_named(rule, "word");
    let value = arena[rule].attribute("value").unwrap_or("").to_owned();
    if words.is_empty() || !ate_or_ite_ending(&value) {
        return Ok(());
    }
    // Changes: A multiplied anion, as in calcium bis(methyl hydrogen
    // phosphate), is one word per copy here, each a bracket around the
    // anion (Blue Book P-16.3.6, tris(dihydrogen phosphate)). Each copy is
    // read as the word without the multiplier would be: the organyl group
    // before the acid can be its ester group (P-67.1.3.2, sodium methyl
    // hydrogen phosphate), and the ester becomes a component of its own.
    // It is the ester group where the name has a space inside the brackets,
    // which the tokenizer removed: bis(methyl phosphonate) is not
    // bis(methylphosphonate). The ester must leave an anion (issue #170).
    let multiplied = words.len() > 1;
    let spaced = multiplied && !name.contains(&value);
    correct_ester_word(arena, rule, words[0], multiplied, spaced, context)?;
    if arena[rule].attribute("wordRule") == Some("ester") {
        let mut last = rule;
        for &word in &words[1..] {
            let copy = arena.grouping("wordRule");
            arena[copy].add_attribute("type", "full");
            arena[copy].add_attribute("wordRule", "simple");
            arena[copy].add_attribute("value", &value);
            arena.detach(word);
            arena.add_child(copy, word);
            arena.insert_after(last, copy);
            correct_ester_word(arena, copy, word, multiplied, spaced, context)?;
            last = copy;
        }
    }
    Ok(())
}

fn correct_ester_word(
    arena: &mut Arena,
    rule: NodeId,
    word: NodeId,
    multiplied: bool,
    spaced: bool,
    context: &mut impl OmittedSpaceContext,
) -> Result<(), ParsingError> {
    let mut children = relevant_children(arena, word, true);
    if multiplied {
        match children[..] {
            [bracket] if arena[bracket].name == "bracket" => {
                children = relevant_children(arena, bracket, true);
            }
            _ => return Ok(()),
        }
    }
    if children.len() < 2 {
        return Ok(());
    }
    let root = *children.last().unwrap();
    let root_group = rightmost_group(arena, root)
        .ok_or_else(|| error("OPSIN bug: Root group missing for omitted ester correction"))?;
    let root_facts = group_facts(arena, root_group, context)?;
    let root_multiplier = multiplier(arena, root)?;
    let functional_count = root_facts.functional_atom_count * root_multiplier;
    if functional_count == 0 {
        return Ok(());
    }
    let subs = &children[..children.len() - 1];
    if subs.len() == 1 && root_multiplier > 1 {
        return Ok(());
    }
    let first = subs[0];
    if !suitable_for_ester(arena, first, functional_count, context)? {
        if arena[first].attribute("locant").is_some() {
            let mut last_without_locant = None;
            for (index, &sub) in subs.iter().enumerate().skip(1) {
                if arena[sub].attribute("locant").is_none() {
                    if !suitable_for_ester(arena, sub, 1, context)? {
                        return Ok(());
                    }
                    last_without_locant = Some(index);
                    break;
                }
            }
            if let Some(index) = last_without_locant
                && substitution_ambiguous(&root_facts, 1)
                && (!multiplied || functional_count > 1)
            {
                transform_multiple(arena, rule, &subs[..=index]);
            }
        }
        return Ok(());
    }
    let first_multiplier = multiplier(arena, first)?;
    let first_group = rightmost_group(arena, first).ok_or_else(|| {
        error("OPSIN bug: Substituent group missing for omitted ester correction")
    })?;
    if (spaced
        || ester_preferred(arena, first_group, first, root_group, subs.len())?
        || subs.len() > 1
            && (all_other_subs_locanted(arena, subs)
                || insufficient_hydrogens(arena, subs, &root_facts, root_multiplier, context)?)
        || (subs.len() == 1 || root_multiplier > 1)
            && substitution_ambiguous(&root_facts, first_multiplier))
        && (!multiplied || functional_count > first_multiplier)
    {
        transform_single(arena, rule, first, context)?;
    }
    Ok(())
}

fn all_other_subs_locanted(arena: &Arena, subs: &[NodeId]) -> bool {
    subs.len() > 1
        && subs
            .iter()
            .skip(1)
            .all(|&id| arena[id].attribute("locant").is_some())
}

fn total_out_valency(facts: &FragmentFacts) -> i64 {
    facts
        .out_atoms
        .iter()
        .map(|out| i64::from(out.valency))
        .sum()
}

fn insufficient_hydrogens(
    arena: &Arena,
    subs: &[NodeId],
    root: &FragmentFacts,
    root_multiplier: usize,
    context: &impl OmittedSpaceContext,
) -> Result<bool, ParsingError> {
    // Changes: Reject multiplication/subtraction overflow instead of allowing
    // a large multiplier to change the omitted-space decision through wrapping.
    let overflow = || error("Omitted-space multiplier arithmetic overflow");
    let mut hydrogens = root
        .substitutable_hydrogen_environments
        .len()
        .checked_mul(root_multiplier)
        .and_then(|count| i64::try_from(count).ok())
        .ok_or_else(overflow)?;
    for &sub in subs.iter().skip(1) {
        let count = i64::try_from(multiplier(arena, sub)?).map_err(|_| overflow())?;
        let valency = total_out_valency(&facts(arena, sub, context)?)
            .checked_mul(count)
            .ok_or_else(overflow)?;
        hydrogens = hydrogens.checked_sub(valency).ok_or_else(overflow)?;
    }
    let count = i64::try_from(multiplier(arena, subs[0])?).map_err(|_| overflow())?;
    let required = total_out_valency(&facts(arena, subs[0], context)?)
        .checked_mul(count)
        .ok_or_else(overflow)?;
    Ok(hydrogens >= 0 && hydrogens < required)
}

fn ester_preferred(
    arena: &Arena,
    substituent_group: NodeId,
    substituent: NodeId,
    root_group: NodeId,
    count: usize,
) -> Result<bool, ParsingError> {
    if arena[substituent].attribute("multiplier").is_some() && multiplier(arena, substituent)? == 1
    {
        return Ok(true);
    }
    let root_name = arena.value(arena[root_group].parent.expect("Root group has no parent"));
    if arena[substituent_group].attribute("type") == Some("chain")
        && arena[substituent_group].attribute("subType") == Some("alkaneStem")
    {
        let group_value = arena.value(substituent_group);
        let sub_value = arena.value(
            arena[substituent_group]
                .parent
                .expect("Substituent group has no parent"),
        );
        // Equivalent to the pinned group-value + "yl-?" full-string match.
        if (sub_value == format!("{group_value}yl") || sub_value == format!("{group_value}yl-"))
            && [
                "format",
                "formate",
                "formoat",
                "formoate",
                "methanat",
                "methanate",
                "methanoat",
                "methanoate",
                "acetat",
                "acetate",
                "acetoat",
                "acetoate",
                "ethanat",
                "ethanate",
                "ethanoat",
                "ethanoate",
            ]
            .iter()
            .any(|ending| root_name.ends_with(ending))
        {
            return Ok(true);
        }
    }
    if (root_name.ends_with("carbamate") || root_name.ends_with("carbamat")) && count >= 2 {
        let mut root = arena[substituent_group]
            .parent
            .expect("Substituent group has no parent");
        while let Some(parent) = arena[root].parent {
            root = parent;
        }
        if arena.children_named(root, "wordRule").len() == 1 {
            return Ok(true);
        }
    }
    Ok(false)
}

fn substitution_ambiguous(facts: &FragmentFacts, multiplier: usize) -> bool {
    if multiplier == 1 && facts.has_default_in_atom {
        return false;
    }
    let hydrogens = &facts.substitutable_hydrogen_environments;
    if hydrogens.len() == multiplier {
        return false;
    }
    let environments: HashSet<_> = hydrogens.iter().collect();
    if environments.len() == 1
        && (multiplier == 1 || Some(multiplier) == hydrogens.len().checked_sub(1))
    {
        return false;
    }
    true
}

fn strip_final_hyphen(arena: &mut Arena, element: NodeId) {
    if let Some(last) = arena.last_child(element)
        && arena[last].name == "hyphen"
    {
        arena.detach(last);
    }
}

fn transform_single(
    arena: &mut Arena,
    rule: NodeId,
    substituent: NodeId,
    context: &mut impl OmittedSpaceContext,
) -> Result<(), ParsingError> {
    // Changes: Validate the source integer range before multiplier-driven cloning.
    arena[rule].set_attribute("wordRule", "ester");
    strip_final_hyphen(arena, substituent);
    arena.detach(substituent);
    let word = arena.grouping("word");
    arena[word].add_attribute("type", "substituent");
    arena.add_child(word, substituent);
    arena.insert_child(rule, word, 0);
    if let Some(attribute) = arena[substituent].remove_attribute("multiplier") {
        let multiplier = attribute
            .value
            .parse::<i32>()
            .ok()
            .and_then(|value| usize::try_from(value).ok())
            .ok_or_else(|| error("OPSIN bug: Invalid ester substituent multiplier"))?;
        for _ in 1..multiplier {
            let clone = context.clone_element(arena, word)?;
            arena.insert_after(word, clone);
        }
    }
    Ok(())
}

fn transform_multiple(arena: &mut Arena, rule: NodeId, substituents: &[NodeId]) {
    arena[rule].set_attribute("wordRule", "ester");
    strip_final_hyphen(arena, *substituents.last().expect("No ester substituents"));
    let word = arena.grouping("word");
    arena[word].add_attribute("type", "substituent");
    for &substituent in substituents {
        arena.detach(substituent);
        arena.add_child(word, substituent);
    }
    arena.insert_child(rule, word, 0);
}
