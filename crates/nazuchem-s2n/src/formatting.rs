// SPDX-License-Identifier: MIT AND (MIT OR Apache-2.0)
// Origin: mixed
// Source: Nazuchem | crates/nazuchem-s2n/src/formatting.rs | issue #45
// Source: Nazuchem | crates/nazuchem-s2n/src/formatting.rs | issue #134
// Source: Nazuchem | crates/nazuchem-s2n/src/formatting.rs | issue #11
// Source: openclatura | src/openclatura/formatting.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: ReShiki | src/naming/openclatura/formatting.rs | 5f559fd482e237307a589dfcb7ea333cee0a6e13
// Source: openclatura | tests/exploratory/test_parenthesis_functions_direct.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/unit/test_stems.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | tests/unit/test_functional_rendering_regressions.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Copyright (c) 2026 lamalab-org
// Copyright (c) 2026 Ameyanagi and ReShiki contributors
// Changes: Openclatura fragment/multiplier formatting with local preferred-prefix and enclosure rules.
// Local contributions: MIT OR Apache-2.0; see NOTICE and PROVENANCE-nazuchem.toml.

//! Name-fragment formatting helpers.
use super::{
    assembly_parts::{AssemblyParts, PrefixKind},
    assembly_utils::is_fully_enclosed,
    rules::{multipliers, stems},
};
use std::cell::Cell;
use std::collections::BTreeMap;
use std::sync::LazyLock;

/// The phenylamino spelling of the retained anilino prefixes (P-62.2.1.1.1), for
/// contractions that read the ligands of an amino prefix, as in phenylcarbamoyl.
// Provenance: fn phenylamino_spelling | original
// Source: ReShiki | src/naming/openclatura/formatting.rs::phenylamino_spelling
pub fn phenylamino_spelling(name: &str) -> &str {
    match name {
        "anilino" => "phenylamino",
        "N-phenylanilino" => "diphenylamino",
        _ => name,
    }
}
// Provenance: fn strip_outer_parentheses | port
// Source: openclatura | src/openclatura/formatting.py::strip_outer_parentheses
pub fn strip_outer_parentheses(name: &str) -> String {
    if name.starts_with('(') && name.ends_with(')') {
        name[1..name.len() - 1].into()
    } else {
        name.into()
    }
}
thread_local! {
    static DESCRIPTORS_OMITTED: Cell<bool> = const { Cell::new(false) };
}

/// Whether names are being rendered without their stereo and isotope
/// descriptors. Every descriptor writer checks this; nothing else does, so the
/// parent, prefixes, locants and numbering are those of the cited name.
// Provenance: fn descriptors_omitted | original
// Source: Nazuchem | crates/nazuchem-s2n/src/formatting.rs::descriptors_omitted
pub(super) fn descriptors_omitted() -> bool {
    DESCRIPTORS_OMITTED.get()
}

/// Run `render` with stereo and isotope descriptors left out of every name,
/// the form that alphanumerical order compares (P-14.5).
// Provenance: fn without_descriptors | original
// Source: Nazuchem | crates/nazuchem-s2n/src/formatting.rs::without_descriptors
pub(super) fn without_descriptors<T>(render: impl FnOnce() -> T) -> T {
    let _restore = RestoreDescriptors(DESCRIPTORS_OMITTED.replace(true));
    render()
}

struct RestoreDescriptors(bool);

// Provenance: impl <RestoreDescriptors as Drop> | original
// Source: Nazuchem | crates/nazuchem-s2n/src/formatting.rs::<RestoreDescriptors as Drop>
impl Drop for RestoreDescriptors {
    fn drop(&mut self) {
        DESCRIPTORS_OMITTED.set(self.0);
    }
}

static BASIC_MULTIPLIERS: LazyLock<Vec<String>> =
    LazyLock::new(|| (2..=1000).map(multipliers::basic).collect());
// Provenance: fn starts_with_multiplier | adapted
// Source: openclatura | src/openclatura/formatting.py::_starts_with_multiplier
pub fn starts_with_multiplier(name: &str) -> bool {
    BASIC_MULTIPLIERS.iter().any(|m| name.starts_with(m))
}
// Provenance: fn is_complex_prefix | adapted
// Source: openclatura | src/openclatura/formatting.py::is_complex_prefix
pub fn is_complex_prefix(name: &str) -> bool {
    name.contains('(')
        || name.starts_with(|c: char| c.is_ascii_digit())
        || name.contains('-')
        || name.contains(' ')
        || starts_with_multiplier(name)
}
// Provenance: fn is_composite_prefix | port
// Source: openclatura | src/openclatura/formatting.py::is_composite_prefix
pub fn is_composite_prefix(name: &str) -> bool {
    for ending in ["ylidyne", "ylidene", "yl"] {
        if let Some(before) = name.strip_suffix(ending) {
            return before.contains("yl");
        }
    }
    name.contains("yl")
}
/// P-16.5.1.1: a compound prefix (P-29.1.2), a simple prefix that carries a
/// substituent prefix such as chloromethyl, cyclopropylmethyl or methoxyamino,
/// is enclosed when it substitutes another prefix, as in [(chloromethyl)sulfanyl]
/// or (methoxyamino)phosphanyl (P-35.4.1). Substituted diazenyl, phosphanyl,
/// arsanyl and stibanyl prefixes such as methyldiazenyl and the compound
/// heterone prefix oxophosphanyl (P-61.6) are compound in the same way.
// Provenance: fn is_compound_prefix | original
// Source: Nazuchem | crates/nazuchem-s2n/src/formatting.rs::is_compound_prefix
pub fn is_compound_prefix(name: &str) -> bool {
    is_composite_prefix(name)
        || ["amino", "diazenyl", "phosphanyl", "arsanyl", "stibanyl"]
            .iter()
            .any(|parent| {
                name.strip_suffix(parent)
                    .is_some_and(|stem| !stem.is_empty())
            })
        || stems::terminal_stem(name)
            .is_some_and(|s| name != format!("{s}yl") && name != format!("cyclo{s}yl"))
}
// Provenance: fn acyl_prefix_stem | adapted
// Source: Nazuchem | crates/nazuchem-s2n/src/formatting.rs::format_element_substituent
fn acyl_prefix_stem(name: &str) -> Option<&str> {
    [
        "sulfonyl",
        "sulfinyl",
        "selenonyl",
        "seleninyl",
        "telluronyl",
        "tellurinyl",
        "carbonyl",
        "thioyl",
        "selenoyl",
        "telluroyl",
    ]
    .iter()
    .find_map(|acyl| name.strip_suffix(acyl))
    .map(|stem| stem.strip_suffix("carbo").unwrap_or(stem))
}

/// A spelled-out parent hydride in an acyl prefix requires enclosing marks
/// even when the prefix is simple (P-16.5.1.4).
// Provenance: fn names_parent_hydride_acyl | adapted
// Source: Nazuchem | crates/nazuchem-s2n/src/formatting.rs::format_element_substituent
fn names_parent_hydride_acyl(name: &str) -> bool {
    acyl_prefix_stem(name).is_some_and(|stem| {
        ["ane", "ene", "yne", "ine"]
            .iter()
            .any(|end| stem.ends_with(end))
    })
}

// Provenance: fn has_numeric_locants | original
// Source: Nazuchem | crates/nazuchem-s2n/src/formatting.rs::has_numeric_locants
fn has_numeric_locants(list: &str) -> bool {
    list.split(',').all(|locant| {
        let after = locant.trim_start_matches(|c: char| c.is_ascii_digit());
        let suffix = after.trim_end_matches('\'');
        after.len() < locant.len()
            && (suffix.is_empty() || suffix.len() == 1 && suffix.as_bytes()[0].is_ascii_lowercase())
    })
}

/// Indicated hydrogen and additive/subtractive hydro terms modify the
/// parent hydride. Read the same base for substitution and ambiguity.
/// P-14.7, P-31.2 and P-16.3.6(c).
// Provenance: fn without_hydride_modifiers | original
// Source: Nazuchem | crates/nazuchem-s2n/src/formatting.rs::without_hydride_modifiers
fn without_hydride_modifiers(mut name: &str) -> &str {
    while let Some((locants, tail)) = name.split_once('-') {
        if locants
            .split(',')
            .all(|locant| locant.strip_suffix('H').is_some_and(has_numeric_locants))
        {
            name = tail;
        } else if has_numeric_locants(locants)
            && let Some((quantity, rest)) = tail.split_once("hydro")
            && BASIC_MULTIPLIERS
                .iter()
                .any(|m| m == quantity.strip_suffix("de").unwrap_or(quantity))
            && !rest.starts_with("xy")
        {
            name = rest.trim_start_matches('-');
        } else {
            break;
        }
    }
    name
}

/// Read the unsubstituted stem and the locanted suffixes of a flattened
/// prefix. Ordinary branch prefixes carry their substitution state directly.
// Provenance: fn has_only_parent_locants | original
// Source: Nazuchem | crates/nazuchem-s2n/src/formatting.rs::has_only_parent_locants
fn has_only_parent_locants(name: &str) -> bool {
    let Some((stem, rest)) = name.split_once('-') else {
        return false;
    };
    if !stem.starts_with(|c: char| c.is_ascii_lowercase()) {
        return false;
    }
    let mut parts = rest.split('-');
    while let Some(locant) = parts.next() {
        if !has_numeric_locants(locant) {
            return false;
        }
        let Some(ending) = parts.next() else {
            return false;
        };
        if matches!(
            ending,
            "yl" | "diyl" | "ylidene" | "ylidyne" | "carbonyl" | "sulfonyl" | "sulfinyl"
        ) {
            return parts.next().is_none();
        }
        if !matches!(ending, "ium" | "ide" | "uide")
            && !["en", "yn"].iter().any(|unsaturation| {
                ending
                    .strip_suffix('e')
                    .unwrap_or(ending)
                    .strip_suffix(unsaturation)
                    .is_some_and(|quantity| {
                        quantity.is_empty() || BASIC_MULTIPLIERS.iter().any(|m| m == quantity)
                    })
            })
        {
            return false;
        }
    }
    false
}

// Provenance: fn retained_italic_prefix | original
// Source: Nazuchem | crates/nazuchem-s2n/src/formatting.rs::retained_italic_prefix
fn retained_italic_prefix(name: &str) -> bool {
    matches!(name, "tert-butyl" | "tert-butoxy")
}

/// P-16.3.6 requires a derived multiplier for several simple prefixes whose
/// multiplication could be read as part of a chain, replacement or fusion.
// Provenance: fn ambiguous_basic_multiplier | original
// Source: Nazuchem | crates/nazuchem-s2n/src/formatting.rs::ambiguous_basic_multiplier
fn ambiguous_basic_multiplier(name: &str) -> bool {
    let unmodified = name
        .strip_prefix("lambda^")
        .or_else(|| name.strip_prefix("λ^"))
        .and_then(|lambda| lambda.split_once('-'))
        .map_or(name, |(_, prefix)| prefix);
    matches!(
        unmodified,
        "sulfanyl"
            | "selanyl"
            | "tellanyl"
            | "phosphanyl"
            | "arsanyl"
            | "stibanyl"
            | "sulfanylidene"
            | "selanylidene"
            | "tellanylidene"
            | "azaniumyl"
            | "azaniumylidene"
            | "oxy"
            | "thio"
            | "methylene"
            | "azo"
            | "sulfonyl"
            | "sulfinyl"
            | "selenonyl"
            | "seleninyl"
            | "telluronyl"
            | "tellurinyl"
            | "diazenyl"
            | "diazonium"
            // P-16.3.6(e): the locanted diazonium prefix also begins with
            // the numerical term di for the two nitrogen atoms.
            | "diazyn-1-ium-1-yl"
            | "disulfanyl"
            | "diselanyl"
            | "ditellanyl"
    ) || {
        let stem = name.trim_start_matches(|c: char| c.is_ascii_digit() || "[,]-".contains(c));
        super::nomenclature::RULES
            .assembly
            .replacement_prefix_order
            .keys()
            .any(|element| {
                stem.starts_with(element.as_str())
                    || ["di", "tri", "tetra"].iter().any(|multiplier| {
                        stem.strip_prefix(multiplier)
                            .is_some_and(|base| base.starts_with(element.as_str()))
                    })
            })
    } || name.split_once('[').is_some_and(|(stem, tail)| {
        // A fusion descriptor ends with a letter inside the square brackets.
        // The brackets of von Baeyer and spiro parents have other meanings.
        !stem.is_empty()
            && !["bicyclo", "tricyclo", "spiro"]
                .iter()
                .any(|p| stem.ends_with(p))
            && tail.split_once(']').is_some_and(|(descriptor, _)| {
                descriptor.ends_with(|c: char| c.is_ascii_lowercase())
            })
    })
}

/// Read a generated IUPAC nuclide descriptor from a flattened prefix.
/// A locant list may precede a nuclide, as in (2,2-2H2)ethyl.
// Provenance: fn without_leading_isotope_descriptor | original
// Source: Nazuchem | crates/nazuchem-s2n/src/formatting.rs::without_leading_isotope_descriptor
fn without_leading_isotope_descriptor(name: &str) -> &str {
    let Some((descriptor, rest)) = name.strip_prefix('(').and_then(|s| s.split_once(')')) else {
        return name;
    };
    // The renderer puts a hyphen after a nuclide only before a numeric locant.
    let prefix = rest.strip_prefix('-').unwrap_or(rest);
    if prefix.is_empty()
        || rest.starts_with('-') && !prefix.starts_with(|c: char| c.is_ascii_digit())
    {
        return name;
    }
    let is_locant =
        |s: &str| has_numeric_locants(s) || matches!(s, "N" | "O" | "S" | "Se" | "Te" | "P" | "As");
    let is_nuclide = |s: &str| {
        let symbol_and_count = s.trim_start_matches(|c: char| c.is_ascii_digit());
        let element = symbol_and_count.trim_end_matches(|c: char| c.is_ascii_digit());
        let count = &symbol_and_count[element.len()..];
        s.len() > symbol_and_count.len()
            && !s.starts_with('0')
            && !count.starts_with('0')
            && super::rules::elements::is_known(element)
            // Generated hydrogen nuclides always give a count. This also
            // keeps indicated hydrogen such as (1H) out of the classifier.
            && (element != "H" || !count.is_empty())
    };
    let mut pending_locant = false;
    for term in descriptor.split(',') {
        if let Some((locant, nuclide)) = term.split_once('-') {
            if !is_locant(locant) || !is_nuclide(nuclide) {
                return name;
            }
            pending_locant = false;
        } else if !pending_locant && is_nuclide(term) {
            continue;
        } else if is_locant(term) {
            pending_locant = true;
        } else {
            return name;
        }
    }
    if pending_locant { name } else { prefix }
}

/// P-16.3.5(a) uses derived multipliers for substituted prefixes. Prefer
/// the classification retained from assembly, when the caller has it.
// Provenance: fn needs_complex_multiplier | adapted
// Source: openclatura | src/openclatura/formatting.py::needs_complex_multiplier
fn needs_complex_multiplier(name: &str, kind: Option<PrefixKind>) -> bool {
    let name = if is_fully_enclosed(name) {
        &name[1..name.len() - 1]
    } else {
        name
    };
    // Stereodescriptors identify the same simple or substituted prefix.
    let name = name
        .strip_prefix('(')
        .and_then(|name| name.split_once(")-"))
        .filter(|(descriptor, _)| {
            descriptor.split(',').all(|term| {
                let symbol = term.trim_start_matches(|c: char| c.is_ascii_digit() || c == '\'');
                matches!(
                    symbol,
                    "R" | "S" | "E" | "Z" | "M" | "P" | "r" | "s" | "m" | "p"
                ) || symbol.len() == 2
                    && symbol.as_bytes()[0].is_ascii_lowercase()
                    && symbol.ends_with(['R', 'S', 'E', 'Z', 'M', 'P'])
            })
        })
        .map_or(name, |(_, prefix)| prefix);
    // P-16.3.4(f): a simple ligand with a nuclide descriptor keeps its basic
    // multiplier and takes enclosing marks. Keep the descriptor in the name.
    let name = without_leading_isotope_descriptor(name);
    let name = without_hydride_modifiers(name);
    if ambiguous_basic_multiplier(name) {
        return true;
    }
    if retained_italic_prefix(name) {
        return false;
    }
    if let Some(kind) = kind {
        return kind == PrefixKind::Substituted;
    }
    // Some functional-class paths assemble strings directly. Their parent
    // locants, retained contractions and enclosed simple acyls still use di.
    if has_only_parent_locants(name) {
        return false;
    }
    if let Some(stem) = acyl_prefix_stem(name) {
        // P-35.3.2: concatenating an oxy prefix with carbonyl is compound.
        if stem.ends_with("oxy") {
            return true;
        }
        // P-16.5.1.4: enclosure alone does not make a simple acyl prefix
        // complex. Substitution on its hydride still takes bis.
        if names_parent_hydride_acyl(name) {
            let is_retained = stem == "hydrazine"
                || super::nomenclature::RULES
                    .retained
                    .substituent_stems
                    .contains_key(stem);
            let is_hydrocarbon = ["ane", "ene", "yne"]
                .iter()
                .find_map(|ending| stem.strip_suffix(ending))
                .is_some_and(|hydride| {
                    let hydride = hydride.strip_prefix("cyclo").unwrap_or(hydride);
                    (1..=1000).any(|length| stems::stem_for(length) == hydride)
                });
            return !(is_retained || is_hydrocarbon);
        }
    }
    // P-63.2.2.1: substituted sulfanyl and uncontracted oxy groups are compound.
    is_compound_prefix(name)
        || [
            "sulfanyl",
            "selanyl",
            "tellanyl",
            "sulfanylidene",
            "selanylidene",
            "tellanylidene",
        ]
        .iter()
        .filter_map(|ending| name.strip_suffix(ending))
        .any(|branch| !branch.is_empty())
        || ["oxy", "carbamoyl"]
            .iter()
            .filter_map(|ending| name.strip_suffix(ending))
            .any(|branch| branch.ends_with("yl") || branch.ends_with(')'))
        || name.contains(['(', ' '])
        || name.contains('-')
        || name.starts_with(|c: char| c.is_ascii_digit())
}

// Provenance: fn needs_multiplier_parentheses | original
// Source: ReShiki | src/naming/openclatura/formatting.rs::needs_simple_multiplier_parentheses
fn needs_multiplier_parentheses(name: &str) -> bool {
    // P-16.3.4: locants, other enclosing marks and leading numerical terms
    // require parentheses even when the multiplier is basic, e.g. di(tridecyl).
    !retained_italic_prefix(name)
        && (is_complex_prefix(name)
            || name.contains('[')
            || name.starts_with("dec")
            || names_parent_hydride_acyl(name))
}

/// Keep existing enclosing marks and enclose a multiplied isotope descriptor
/// with its prefix. Classify the prefix without the descriptor (P-82.2.1).
// Provenance: fn multiply_prefix | original
// Source: Nazuchem | crates/nazuchem-s2n/src/formatting.rs::multiply_prefix
pub fn multiply_prefix(
    name: &str,
    plain: &str,
    count: usize,
    kind: Option<PrefixKind>,
    enclose: bool,
) -> String {
    let derived = needs_complex_multiplier(plain, kind);
    let multiplier = if derived {
        multipliers::complex_(count)
    } else {
        multipliers::basic(count)
    };
    let enclosed = is_fully_enclosed(name);
    let bare = if enclosed {
        &name[1..name.len() - 1]
    } else {
        name
    };
    if !enclose && name == plain && retained_italic_prefix(bare) {
        format!("{multiplier}-{bare}")
    } else if enclosed {
        format!("{multiplier}{name}")
    } else if derived || enclose || name != plain || needs_multiplier_parentheses(name) {
        format!("{multiplier}({name})")
    } else {
        format!("{multiplier}{name}")
    }
}

// Provenance: fn format_multiplier | adapted
// Source: openclatura | src/openclatura/formatting.py::format_multiplier
pub fn format_multiplier(name: &str, count: usize, safe_enclose: bool) -> String {
    if count == 1 {
        return if (safe_enclose || is_complex_prefix(name)) && !is_fully_enclosed(name) {
            format!("({name})")
        } else {
            name.into()
        };
    }
    multiply_prefix(name, name, count, None, false)
}

// Provenance: fn count_names | port
// Source: openclatura | src/openclatura/formatting.py::count_names
pub fn count_names(names: &[String]) -> BTreeMap<String, usize> {
    let mut result = BTreeMap::new();
    for name in names {
        *result.entry(name.clone()).or_insert(0) += 1;
    }
    result
}
// Provenance: fn format_counted_prefixes | port
// Source: openclatura | src/openclatura/formatting.py::format_counted_prefixes
pub fn format_counted_prefixes(names: &[String]) -> String {
    let counts = count_names(names);
    let safe = counts.len() > 1 || counts.keys().any(|n| is_complex_prefix(n));
    counts
        .iter()
        .map(|(n, c)| format_multiplier(n, *c, safe))
        .collect()
}
// Provenance: fn format_center_ligands_by | adapted
// Source: openclatura | src/openclatura/formatting.py::format_center_ligands
pub fn format_center_ligands_by<K: Ord>(names: &[String], sort_key: impl Fn(&str) -> K) -> String {
    let counts = count_names(names);
    let mut ordered: Vec<_> = counts.keys().collect();
    ordered.sort_by_key(|s| sort_key(s));
    let mut rendered = String::new();
    for (i, name) in ordered.iter().enumerate() {
        let count = counts[*name];
        if count > 1 {
            rendered.push_str(&multiply_prefix(name, name, count, None, i > 0));
            continue;
        }
        let merges = ordered[i + 1..]
            .iter()
            .any(|later| counts[*later] > 1 && name.ends_with(later.as_str()));
        // P-16.5.1.3.1: later distinct ligands are enclosed, with a simple
        // multiplying prefix outside the parentheses: ethyldi(methyl).
        if (is_complex_prefix(name) || merges || i > 0) && !is_fully_enclosed(name) {
            rendered.push_str(&format!("({name})"));
        } else {
            rendered.push_str(name);
        }
    }
    rendered
}
// Provenance: fn substituted_alkoxy_prefix | adapted
// Source: openclatura | src/openclatura/formatting.py::substituted_alkoxy_prefix
pub fn substituted_alkoxy_prefix(branch: &str, parts: &AssemblyParts) -> Option<String> {
    // P-63.2.2.2: the retained parent, not a suffix of a nested prefix, must
    // be methyl, ethyl, propyl, butyl or phenyl. tert-Butoxy is not substitutable.
    if parts.attachment_order != 1
        || parts.attachment_locant != "1"
        || parts.is_acyl_substituent
        || parts.principal_group.is_some()
        || parts.retained_substituent_name.is_some()
        || !parts.a_prefixes.is_empty()
        || !parts.parent_charges.is_empty()
        || !parts.hydro_operations.is_empty()
    {
        return None;
    }
    let stem = if parts.is_ring {
        (parts.retained_name.as_deref() == Some("benzene") && !parts.is_ring_assembly())
            .then_some("phen")?
    } else {
        if parts.retained_name.is_some() || !parts.unsaturations.is_empty() {
            return None;
        }
        match parts.parent_length {
            1 => "meth",
            2 => "eth",
            3 => "prop",
            4 => "but",
            _ => return None,
        }
    };
    let branch = strip_outer_parentheses(branch);
    let terminal = format!("{stem}yl");
    let prefix = branch.strip_suffix(&terminal)?;
    Some(if prefix.is_empty() {
        format!("{stem}oxy")
    } else {
        format!("({prefix}{stem}oxy)")
    })
}
// Provenance: fn oxy_prefix_from_branch | adapted
// Source: openclatura | src/openclatura/formatting.py::oxy_prefix_from_branch
pub fn oxy_prefix_from_branch(branch: &str) -> String {
    let branch = strip_outer_parentheses(branch);
    if let Some(name) = super::nomenclature::RULES
        .heteroatoms
        .alkyl_oxy_prefixes
        .get(&branch)
    {
        return name.clone();
    }
    format_element_substituent("", &branch, "oxy", false)
}
// Provenance: fn format_element_substituent | adapted
// Source: openclatura | src/openclatura/formatting.py::format_element_substituent
pub fn format_element_substituent(
    stereo_prefix: &str,
    branch: &str,
    suffix: &str,
    is_double: bool,
) -> String {
    let branch = strip_outer_parentheses(branch);
    let suffix = format!("{suffix}{}", if is_double { "idene" } else { "" });
    // An acyl prefix that names a parent hydride, as methanesulfonyl or
    // ethanethioyl, is enclosed although simple (P-16.5.1.4).
    if is_complex_prefix(&branch)
        || is_compound_prefix(&branch)
        || names_parent_hydride_acyl(&branch)
    {
        format!("({stereo_prefix}({branch}){suffix})")
    } else {
        format!("({stereo_prefix}{branch}{suffix})")
    }
}

#[cfg(test)]
// Provenance: mod tests | original
// Source: ReShiki | src/naming/openclatura/formatting.rs::tests
mod tests {
    use super::format_multiplier;

    #[test]
    // Provenance: fn tests::simple_dec_prefixes_keep_basic_multipliers_and_parentheses | adapted
    // Source: ReShiki | src/naming/openclatura/formatting.rs::tests::simple_dec_prefixes_keep_basic_multipliers_and_parentheses
    // Source: openclatura | tests/exploratory/test_parenthesis_functions_direct.py::test_starts_with_multiplier_true_positives
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn simple_dec_prefixes_keep_basic_multipliers_and_parentheses() {
        for (name, count, expected) in [
            ("decyl", 1, "decyl"),
            ("decyl", 2, "di(decyl)"),
            ("decylidene", 2, "di(decylidene)"),
            ("decyl", 3, "tri(decyl)"),
            ("ethylidene", 2, "diethylidene"),
            // Issue #120, P-16.3.5 (a): bis(bromomethyl).
            ("chloromethyl", 2, "bis(chloromethyl)"),
            ("methyl", 2, "dimethyl"),
            ("ethenyl", 2, "diethenyl"),
            // P-16.3.6(a): the lambda convention leaves the ambiguity.
            ("lambda^4-sulfanyl", 2, "bis(lambda^4-sulfanyl)"),
        ] {
            assert_eq!(format_multiplier(name, count, false), expected);
        }
    }

    #[test]
    // Provenance: fn tests::oxy_prefix_keeps_parent_hydride_acyl_enclosed | adapted
    // Source: ReShiki | src/naming/openclatura/formatting.rs::tests::oxy_prefix_keeps_parent_hydride_acyl_enclosed
    // Source: openclatura | tests/unit/test_stems.py::test_substituted_alkoxy_prefix_rejects_unrelated_branches
    // Source: openclatura | tests/unit/test_functional_rendering_regressions.py::test_multiplier_preserves_ligand_boundaries
    // Source: openclatura | tests/exploratory/test_parenthesis_functions_direct.py::test_starts_with_multiplier_true_positives
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn oxy_prefix_keeps_parent_hydride_acyl_enclosed() {
        // P-16.5.1.4: (methanethioyl)oxy, as in (cyclohexanecarbonyl)sulfanyl (P-57.4).
        assert_eq!(
            super::oxy_prefix_from_branch("(methanethioyl)"),
            "((methanethioyl)oxy)"
        );
        assert_eq!(super::oxy_prefix_from_branch("pentyl"), "(pentyloxy)");
        assert_eq!(super::oxy_prefix_from_branch("ethyl"), "ethoxy");
    }
}
