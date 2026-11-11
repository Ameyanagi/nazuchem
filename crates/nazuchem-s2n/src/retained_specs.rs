// SPDX-License-Identifier: MIT
// Origin: port
// Source: openclatura | src/openclatura/retained_specs.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Copyright (c) 2026 lamalab-org

//! Retained parent attachment spelling and locant policy.
use super::nomenclature::NAMER_RULES;
#[derive(Clone, Debug)]
pub struct RetainedParentSpec {
    pub substituent_stem: Option<String>,
    pub substituent_terminal: Option<String>,
    pub print_substituent_locant: bool,
}
pub fn retained_parent_spec(name: Option<&str>) -> Option<RetainedParentSpec> {
    let name = name.filter(|n| !n.is_empty())?;
    let data = NAMER_RULES.retained_substituent_stems.get(name);
    let stem = data.map(|d| d.0.clone());
    let terminal = data.map(|d| d.1.clone());
    let normalized = stem.as_deref().unwrap_or(name).to_lowercase();
    let normalized = normalized.trim_end_matches('e');
    let needs_locant = NAMER_RULES
        .ambiguous_connection_substituent_stems
        .iter()
        .any(|s| normalized.ends_with(s.as_str()));
    Some(RetainedParentSpec {
        substituent_stem: stem,
        substituent_terminal: terminal,
        print_substituent_locant: needs_locant,
    })
}
