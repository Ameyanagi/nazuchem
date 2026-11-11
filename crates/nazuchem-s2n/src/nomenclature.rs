// SPDX-License-Identifier: MIT
// Origin: adapted
// Source: openclatura | src/openclatura/nomenclature.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | src/openclatura/naming_data.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: openclatura | src/openclatura/chalcogen_vocabulary.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Copyright (c) 2026 lamalab-org
// Changes: Openclatura registry types, family/capability derivation and immutable JSON lookup adapted to Rust.
// Unresolved: Small Rust JSON conversion helpers have no one-to-one upstream function boundary; adapted records their role in the copied registry/schema rather than asserting literal function translation.

//! Nomenclature tables transcribed from Openclatura's checked-in rule registry.
use super::numbering::NumberingCriterion;
use super::rules::retained::MonocycleSpec;
use serde::Deserialize;
use std::collections::{BTreeMap, BTreeSet};
use std::ops::Deref;
use std::sync::LazyLock;

type StringMap = BTreeMap<String, String>;
type StringSet = BTreeSet<String>;

/// One registry section; it derefs to its typed values.
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Section<T> {
    #[serde(rename = "bluebook_rule")]
    _bluebook_rule: String,
    values: T,
}
impl<T> Deref for Section<T> {
    type Target = T;
    fn deref(&self) -> &T {
        &self.values
    }
}

/// The whole of `namer_rules.json`. Loading rejects a missing section, a
/// missing required column and any key that is not declared here.
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NamerRules {
    pub numbering_criteria: Section<NumberingCriteria>,
    pub single_atom_cations: Section<StringSet>,
    pub single_atom_anions: Section<StringMap>,
    pub mononuclear_parent_hydrides: Section<StringMap>,
    pub retained_homonuclear_chain_names: Section<StringMap>,
    pub replacement_parent_oxoacid_specs: Section<Vec<OxoacidSpec>>,
    pub retained_ionic_n_parents: Section<StringMap>,
    pub saturated_n_ring_ionic_parents: Section<BTreeMap<usize, String>>,
    pub parent_charge_suffixes: Section<BTreeMap<String, ParentChargeSuffixRule>>,
    pub replacement_charge_prefixes: Section<StringMap>,
    pub replacement_charge_states: Section<Vec<ReplacementChargeState>>,
    pub heteroatom_charge_prefixes: Section<StringMap>,
    pub alkyl_oxy_prefixes: Section<StringMap>,
    pub simple_sulfanyl_prefixes: Section<StringSet>,
    pub simple_selanyl_prefixes: Section<StringSet>,
    pub retained_sulfonyl_ligand_contractions: Section<Vec<SulfonylLigandContraction>>,
    pub retained_sulfonyl_group_names: Section<Vec<SulfonylGroupName>>,
    pub halogen_prefixes: Section<StringMap>,
    pub halogen_lambda_suffixes: Section<StringMap>,
    pub retained_ring_elements: Section<StringSet>,
    pub ring_descriptor_templates: Section<StringMap>,
    pub polycycle_prefixes: Section<BTreeMap<usize, String>>,
    pub direct_group_prefixes: Section<StringMap>,
    pub chain_external_carbonyl_groups: Section<StringSet>,
    pub prefix_groups_to_skip: Section<StringSet>,
    pub ester_like_prefix_groups: Section<StringSet>,
    pub peroxy_ester_groups: Section<StringSet>,
    pub amide_like_prefix_groups: Section<StringSet>,
    pub carboxy_prefix_groups: Section<StringSet>,
    pub cyano_prefix_groups: Section<StringSet>,
    pub acid_halide_prefixes: Section<StringMap>,
    pub peroxy_acid_prefix_groups: Section<StringSet>,
    pub sulfonyl_prefix_groups: Section<StringSet>,
    pub direct_prefix_groups: Section<StringMap>,
    pub front_modifier_principal_groups: Section<StringSet>,
    pub hydrazone_principal_groups: Section<StringSet>,
    pub chalcogen_elements: Section<BTreeMap<String, ChalcogenElement>>,
    pub chalcogen_simple_group_keys: Section<Vec<SimpleGroupKey>>,
    pub chalcogenide_class_names: Section<StringMap>,
    pub chalcogen_peroxide_class_names: Section<Vec<PeroxideClassName>>,
    pub chalcogen_peroxol_rules: Section<Vec<PeroxolRule>>,
    pub chalcogen_standard_acyl_keys: Section<Vec<StandardAcylKey>>,
    pub chalcogen_standard_peroxy_acyl_keys: Section<Vec<StandardPeroxyAcylKey>>,
    pub chalcogen_acyl_leaving_groups: Section<BTreeMap<String, AcylLeavingGroup>>,
    pub chalcogen_central_acid_origins: Section<Vec<CentralAcidOrigin>>,
    pub chalcogen_rule_policies: Section<ChalcogenRulePolicies>,
    pub chalcogen_templates: Section<StringMap>,
    pub functional_groups: Section<BTreeMap<String, FunctionalGroupRow>>,
    pub retained_substituent_stems: Section<BTreeMap<String, (String, String)>>,
    pub ambiguous_connection_substituent_stems: Section<StringSet>,
    pub suffix_nitrogen_markers: Section<Vec<String>>,
    pub replacement_prefix_order: Section<BTreeMap<String, usize>>,
    pub unsaturation_order: Section<BTreeMap<String, usize>>,
    pub substituent_attachment_suffixes: Section<StringMap>,
    pub postprocess_literal_replacements: Section<Vec<LiteralReplacement>>,
    pub retained_chain_functional_parents: Section<Vec<ChainFunctionalParent>>,
    pub retained_principal_citations: Section<Vec<PrincipalCitation>>,
    pub retained_monocycle_specs: Section<Vec<MonocycleSpec>>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NumberingCriteria {
    pub chain: Vec<NumberingCriterion>,
    pub ring: Vec<NumberingCriterion>,
}
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OxoacidSpec {
    pub name: String,
    pub ester_suffix: String,
    pub central: String,
    pub single_o: usize,
    pub double_o: usize,
    pub charge: Option<i32>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ReplacementChargeState {
    pub symbol: String,
    pub charge: i32,
    pub valence: usize,
    pub prefix: String,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SulfonylLigandContraction {
    pub ligand_ending: String,
    pub suffix: String,
    pub name: String,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SulfonylGroupName {
    pub branch: String,
    pub suffix: String,
    pub name: String,
}
/// Polonium has a rank only; the vocabulary columns are absent for it.
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ChalcogenElement {
    pub rank: i64,
    pub additive_prefix: Option<String>,
    pub substituent_suffix: Option<String>,
    pub replacement_infix: Option<String>,
    pub acyl_center_prefix: Option<String>,
    pub nitrogen_chalcogenide_citation: Option<String>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SimpleGroupKey {
    pub family: String,
    pub derivative: String,
    pub element: String,
    pub key: String,
    #[serde(default)]
    pub external: bool,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PeroxideClassName {
    pub elements: Vec<String>,
    pub name: String,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PeroxolRule {
    pub first: String,
    pub terminal: String,
    pub suffix: String,
    pub rank: i64,
    pub capabilities: StringSet,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StandardAcylKey {
    pub derivative: String,
    pub double_element: String,
    pub single_element: Option<String>,
    pub key: String,
    #[serde(default)]
    pub external: bool,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StandardPeroxyAcylKey {
    pub derivative: String,
    pub key: String,
    #[serde(default)]
    pub external: bool,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AcylLeavingGroup {
    pub word: String,
    pub prefix: String,
    pub rank: i64,
    pub chain_key: String,
    pub external_key: String,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CentralAcidOrigin {
    pub element: String,
    pub double_bond_count: usize,
    pub stem: String,
    pub rank: i64,
    pub linkage_prefix: String,
    pub acid_key: Option<String>,
    pub ester_key: Option<String>,
    pub amide_key: Option<String>,
    pub acid_prefix: Option<String>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ChalcogenRulePolicies {
    pub nitrile_chalcogenide: NitrilePolicy,
    pub peroxol: RulePolicy,
    pub acyl: BTreeMap<String, RulePolicy>,
    pub peroxy_acyl: BTreeMap<String, RulePolicy>,
    pub anhydride: RulePolicy,
    pub imide_route: ImideRoutePolicy,
    pub central_peroxy: BTreeMap<String, RulePolicy>,
    pub central: BTreeMap<String, RulePolicy>,
}
/// The default is the policy of a derivative that has no row: no templates,
/// no seniority, no families.
#[derive(Clone, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RulePolicy {
    pub prefix: Option<String>,
    pub suffix: Option<String>,
    pub citation_suffix: Option<String>,
    pub seniority: Option<i64>,
    #[serde(default)]
    pub ranked_seniority: bool,
    pub families: StringSet,
    #[serde(default)]
    pub capabilities: StringSet,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NitrilePolicy {
    pub prefix: String,
    pub suffix: String,
    pub seniority: i64,
    pub families: StringSet,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ImideRoutePolicy {
    pub key: String,
    pub seniority: i64,
    pub families: StringSet,
    pub capabilities: StringSet,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FunctionalGroupRow {
    role: String,
    prefix: Option<String>,
    suffix: Option<String>,
    multi_suffix: Option<MultiSuffix>,
    seniority: Option<Seniority>,
    #[serde(default)]
    suffix_with_locant: bool,
    needs_locant: bool,
    #[serde(default)]
    citation_suffixes: StringMap,
    families: Option<StringSet>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct MultiSuffix {
    multiplier_positions: Vec<usize>,
}
#[derive(Deserialize)]
#[serde(untagged)]
enum Seniority {
    Rank(i64),
    Ranks(Vec<i64>),
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LiteralReplacement {
    pub pattern: String,
    pub replacement: String,
    #[serde(rename = "category")]
    _category: String,
    #[serde(rename = "reason")]
    _reason: String,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ChainFunctionalParent {
    length: usize,
    group: String,
    count: usize,
    name: String,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PrincipalCitation {
    mode: String,
    name: String,
    citation: String,
}

// The immutable embedded registry is checked in with the port and verified by the oracle tests.
#[allow(clippy::expect_used)]
pub static NAMER_RULES: LazyLock<NamerRules> = LazyLock::new(|| {
    serde_json::from_str(include_str!("data/namer_rules.json"))
        .expect("checked-in nomenclature data")
});
pub fn chalcogen_replacement_infix(element: &str) -> Option<&'static str> {
    NAMER_RULES
        .chalcogen_elements
        .get(element)?
        .replacement_infix
        .as_deref()
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct FunctionalGroupRule {
    pub key: String,
    pub role: String,
    pub prefix: Option<String>,
    pub suffix: Option<String>,
    pub suffix_multiplier_positions: Vec<usize>,
    pub seniority: Option<Vec<i64>>,
    pub suffix_with_locant: bool,
    pub needs_locant: bool,
    pub families: BTreeSet<String>,
    pub capabilities: BTreeSet<String>,
    pub citation_suffixes: BTreeMap<String, String>,
}
impl FunctionalGroupRule {
    pub fn has_capability(&self, capability: &str) -> bool {
        self.capabilities.contains(capability)
    }
    pub fn suffix_for_citation(&self, mode: Option<&str>) -> Option<&str> {
        mode.and_then(|m| self.citation_suffixes.get(m).map(String::as_str))
            .or(self.suffix.as_deref())
    }
}
#[derive(Clone, Debug, Default)]
pub struct FunctionalGroupRules {
    pub by_key: BTreeMap<String, FunctionalGroupRule>,
}
impl FunctionalGroupRules {
    pub fn get(&self, key: &str) -> Option<&FunctionalGroupRule> {
        self.by_key.get(key)
    }
    pub fn prefix_for(&self, key: &str) -> Option<&str> {
        self.get(key)?.prefix.as_deref()
    }
    pub fn cited_prefix_for(&self, key: &str) -> Option<String> {
        let r = self.get(key)?;
        let p = r.prefix.as_ref()?;
        Some(if r.has_capability("acyl_leaving_group") {
            format!("({p})")
        } else {
            p.clone()
        })
    }
    pub fn direct_subgraph_prefix_for(&self, key: &str) -> Option<&str> {
        let r = self.get(key)?;
        (r.role == "prefix" || r.families.contains("direct_group"))
            .then_some(r.prefix.as_deref())
            .flatten()
    }
    pub fn keys_with_family(&self, family: &str) -> BTreeSet<String> {
        self.by_key
            .iter()
            .filter(|(_, r)| r.families.contains(family))
            .map(|(k, _)| k.clone())
            .collect()
    }
}
const CAPABILITIES: &[&str] = &[
    "acid_halide",
    "acyl_leaving_group",
    "amide_like",
    "anhydride",
    "blocks_structural_parent",
    "central_acid",
    "chain_external_carbonyl",
    "element_locanted_suffix",
    "ester_like",
    "front_modifier",
    "hydrazide",
    "promotes_member_groups",
    "urea",
];
fn derived_families(key: &str) -> BTreeSet<String> {
    let rules = &*NAMER_RULES;
    let mut out = BTreeSet::new();
    for (family, table) in [
        (
            "chain_external_carbonyl",
            &*rules.chain_external_carbonyl_groups,
        ),
        ("prefix_skip", &*rules.prefix_groups_to_skip),
        ("ester_like", &*rules.ester_like_prefix_groups),
        ("peroxy_ester", &*rules.peroxy_ester_groups),
        ("amide_like", &*rules.amide_like_prefix_groups),
        ("carboxy_prefix", &*rules.carboxy_prefix_groups),
        ("cyano_prefix", &*rules.cyano_prefix_groups),
        ("peroxy_acid", &*rules.peroxy_acid_prefix_groups),
        ("sulfonyl", &*rules.sulfonyl_prefix_groups),
        ("front_modifier", &*rules.front_modifier_principal_groups),
        ("hydrazone", &*rules.hydrazone_principal_groups),
    ] {
        if table.contains(key) {
            out.insert(family.to_owned());
        }
    }
    for (family, table) in [
        ("acid_halide", &*rules.acid_halide_prefixes),
        ("direct_prefix", &*rules.direct_prefix_groups),
        ("direct_group", &*rules.direct_group_prefixes),
    ] {
        if table.contains_key(key) {
            out.insert(family.to_owned());
        }
    }
    out
}
fn family_capabilities(families: &BTreeSet<String>) -> impl Iterator<Item = String> {
    families
        .iter()
        .filter(|s| CAPABILITIES.contains(&s.as_str()))
        .cloned()
}
/// A principal-group rule that carries a chalcogen policy's families and
/// capabilities; the caller fills in the affixes and the seniority.
pub fn policy_rule(
    key: &str,
    families: &BTreeSet<String>,
    capabilities: &BTreeSet<String>,
) -> FunctionalGroupRule {
    let mut capabilities = capabilities.clone();
    capabilities.extend(family_capabilities(families));
    FunctionalGroupRule {
        key: key.to_owned(),
        role: "principal".into(),
        suffix_multiplier_positions: vec![0],
        needs_locant: true,
        families: families.clone(),
        capabilities,
        ..Default::default()
    }
}
// The acid-halide policy belongs to the immutable checked-in registry schema.
#[allow(clippy::expect_used)]
fn functional_groups() -> FunctionalGroupRules {
    let rules = &*NAMER_RULES;
    let mut by_key = rules
        .functional_groups
        .iter()
        .map(|(key, row)| {
            let families = row
                .families
                .clone()
                .unwrap_or_else(|| derived_families(key));
            let rule = FunctionalGroupRule {
                key: key.clone(),
                role: row.role.clone(),
                prefix: row.prefix.clone(),
                suffix: row.suffix.clone(),
                suffix_multiplier_positions: row
                    .multi_suffix
                    .as_ref()
                    .map_or_else(|| vec![0], |m| m.multiplier_positions.clone()),
                seniority: row.seniority.as_ref().map(|s| match s {
                    Seniority::Rank(rank) => vec![*rank],
                    Seniority::Ranks(ranks) => ranks.clone(),
                }),
                suffix_with_locant: row.suffix_with_locant,
                needs_locant: row.needs_locant,
                capabilities: family_capabilities(&families).collect(),
                families,
                citation_suffixes: row.citation_suffixes.clone(),
            };
            (key.clone(), rule)
        })
        .collect::<BTreeMap<_, _>>();
    let policy = &rules.chalcogen_rule_policies.acyl["acid_halide"];
    let seniority = policy.seniority.expect("checked-in acyl seniority");
    for item in rules.chalcogen_acyl_leaving_groups.values() {
        for (external, key, suffix) in [
            (false, &item.chain_key, "oyl"),
            (true, &item.external_key, "carbonyl"),
        ] {
            let mut r = policy_rule(key, &policy.families, &policy.capabilities);
            r.prefix = Some(format!("{}carbonyl", item.prefix));
            r.suffix = Some(format!("{suffix} {}", item.word));
            r.seniority = Some(vec![seniority, item.rank]);
            r.suffix_multiplier_positions = vec![0, 1];
            r.suffix_with_locant = external;
            by_key.insert(key.clone(), r);
        }
    }
    FunctionalGroupRules { by_key }
}
#[derive(Debug, Default)]
pub struct RetainedNameRules {
    pub substituent_stems: BTreeMap<String, (String, String)>,
    pub chain_functional_parents: BTreeMap<(usize, String, usize), String>,
    pub principal_citations: BTreeMap<(String, String), String>,
}
#[derive(Debug, Default)]
pub struct HeteroatomRules {
    pub alkyl_oxy_prefixes: BTreeMap<String, String>,
    pub simple_sulfanyl_prefixes: BTreeSet<String>,
    pub simple_selanyl_prefixes: BTreeSet<String>,
}
#[derive(Debug, Default)]
pub struct PrefixRules {
    pub peroxy_ester_groups: BTreeSet<String>,
    pub amide_like_groups: BTreeSet<String>,
}
#[derive(Debug, Default)]
pub struct ComponentRules {
    pub mononuclear_parent_hydrides: BTreeMap<String, String>,
    pub retained_homonuclear_chain_names: BTreeMap<String, String>,
    pub replacement_parent_oxoacid_specs: Vec<OxoacidSpec>,
}
#[derive(Debug, Default)]
pub struct IonRules {
    pub single_atom_cations: BTreeSet<String>,
    pub single_atom_anions: BTreeMap<String, String>,
}
#[derive(Clone, Debug, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ParentChargeSuffixRule {
    pub suffix: String,
    #[serde(rename = "reason")]
    _reason: String,
}
#[derive(Debug, Default)]
pub struct ChargeRules {
    pub retained_ionic_n_parents: BTreeMap<String, String>,
    pub saturated_n_ring_ionic_parents: BTreeMap<usize, String>,
    pub parent_charge_suffixes: BTreeMap<String, ParentChargeSuffixRule>,
    pub replacement_charge_prefixes: BTreeMap<String, String>,
    pub heteroatom_charge_prefixes: BTreeMap<String, String>,
}
#[derive(Debug, Default)]
pub struct AssemblyRules {
    pub replacement_prefix_order: BTreeMap<String, usize>,
    pub unsaturation_order: BTreeMap<String, usize>,
    pub substituent_attachment_suffixes: BTreeMap<String, String>,
    pub ambiguous_connection_substituent_stems: BTreeSet<String>,
    pub suffix_nitrogen_markers: Vec<String>,
}
#[derive(Debug)]
pub struct NomenclatureRegistry {
    pub retained: RetainedNameRules,
    pub heteroatoms: HeteroatomRules,
    pub prefixes: PrefixRules,
    pub components: ComponentRules,
    pub ions: IonRules,
    pub charges: ChargeRules,
    pub assembly: AssemblyRules,
    pub functional_groups: FunctionalGroupRules,
}
pub static RULES: LazyLock<NomenclatureRegistry> = LazyLock::new(|| {
    let rules = &*NAMER_RULES;
    NomenclatureRegistry {
        functional_groups: functional_groups(),
        retained: RetainedNameRules {
            substituent_stems: rules.retained_substituent_stems.clone(),
            chain_functional_parents: rules
                .retained_chain_functional_parents
                .iter()
                .map(|v| ((v.length, v.group.clone(), v.count), v.name.clone()))
                .collect(),
            principal_citations: rules
                .retained_principal_citations
                .iter()
                .map(|v| ((v.mode.clone(), v.name.clone()), v.citation.clone()))
                .collect(),
        },
        heteroatoms: HeteroatomRules {
            alkyl_oxy_prefixes: rules.alkyl_oxy_prefixes.clone(),
            simple_sulfanyl_prefixes: rules.simple_sulfanyl_prefixes.clone(),
            simple_selanyl_prefixes: rules.simple_selanyl_prefixes.clone(),
        },
        prefixes: PrefixRules {
            peroxy_ester_groups: rules.peroxy_ester_groups.clone(),
            amide_like_groups: rules.amide_like_prefix_groups.clone(),
        },
        components: ComponentRules {
            mononuclear_parent_hydrides: rules.mononuclear_parent_hydrides.clone(),
            retained_homonuclear_chain_names: rules.retained_homonuclear_chain_names.clone(),
            replacement_parent_oxoacid_specs: rules.replacement_parent_oxoacid_specs.clone(),
        },
        ions: IonRules {
            single_atom_cations: rules.single_atom_cations.clone(),
            single_atom_anions: rules.single_atom_anions.clone(),
        },
        charges: ChargeRules {
            retained_ionic_n_parents: rules.retained_ionic_n_parents.clone(),
            saturated_n_ring_ionic_parents: rules.saturated_n_ring_ionic_parents.clone(),
            parent_charge_suffixes: rules.parent_charge_suffixes.clone(),
            replacement_charge_prefixes: rules.replacement_charge_prefixes.clone(),
            heteroatom_charge_prefixes: rules.heteroatom_charge_prefixes.clone(),
        },
        assembly: AssemblyRules {
            replacement_prefix_order: rules.replacement_prefix_order.clone(),
            unsaturation_order: rules.unsaturation_order.clone(),
            substituent_attachment_suffixes: rules.substituent_attachment_suffixes.clone(),
            ambiguous_connection_substituent_stems: rules
                .ambiguous_connection_substituent_stems
                .clone(),
            suffix_nitrogen_markers: rules.suffix_nitrogen_markers.clone(),
        },
    }
});

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn registry_rows_reject_misspelled_and_missing_columns() {
        let row = |json: &str| serde_json::from_str::<SimpleGroupKey>(json);
        let complete =
            r#"{"family":"hydroxy","derivative":"alcohol","element":"O","key":"alcohol"}"#;
        assert!(row(complete).is_ok());
        // A misspelled optional column no longer reads as its default.
        assert!(row(&complete.replace('}', r#","extrenal":true}"#)).is_err());
        assert!(row(&complete.replace(r#","key":"alcohol""#, "")).is_err());
    }
}
