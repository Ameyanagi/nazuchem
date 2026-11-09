// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | crates/nazuchem-ir/src/molecule.rs | Phase B step 8a

//! IR v1 stores atom identity and stated structure without chemistry perception.

#[cfg(any(test, feature = "schema"))]
use std::borrow::Cow;
use std::collections::{BTreeMap, BTreeSet};

#[cfg(any(test, feature = "schema"))]
use schemars::{JsonSchema, Schema, SchemaGenerator, json_schema};
use serde::ser::SerializeStruct;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use thiserror::Error;

/// Version written in every serialized molecule.
pub const IR_VERSION: u64 = 1;

/// Published JSON Schema for the serialized IR v1 document.
///
/// The schema checks the JSON representation. Cross-reference and graph
/// constraints additionally require [`Molecule::validate`].
pub const JSON_SCHEMA: &str = include_str!("../schema/nazuchem-ir-1.schema.json");

/// Generate the published schema from the same Rust types used for JSON input.
#[cfg(any(test, feature = "schema"))]
pub fn json_schema() -> Schema {
    let mut schema = schemars::schema_for!(Molecule);
    if let Some(object) = schema.as_object_mut() {
        for value in object.values_mut() {
            add_index_bounds(value);
        }
    }
    schema
}

// Schemars emits uint32 as a format annotation. JSON Schema treats formats as
// annotations by default, so publish explicit bounds for interoperable validators.
#[cfg(any(test, feature = "schema"))]
fn add_index_bounds(value: &mut serde_json::Value) {
    match value {
        serde_json::Value::Object(object) => {
            if object.get("format").and_then(serde_json::Value::as_str) == Some("uint32") {
                object.insert("maximum".into(), u32::MAX.into());
            }
            for value in object.values_mut() {
                add_index_bounds(value);
            }
        }
        serde_json::Value::Array(array) => {
            for value in array {
                add_index_bounds(value);
            }
        }
        _ => {}
    }
}

/// A structure whose atom and bond indices are stable identities.
///
/// Bonds carry localized Kekulé orders; aromaticity and coordinates are absent.
/// Fields remain public for adapters and callers constructing structures. Call
/// [`Self::validate`] after mutation; JSON input and output validate automatically.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct Molecule {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
    pub stereo: Vec<Stereo>,
    pub stereo_groups: Vec<StereoGroup>,
    pub position_variation: Vec<PositionVariation>,
    pub polymers: Vec<RepeatUnit>,
    pub parent: Vec<u32>,
}

/// An atom with its stated hydrogen origin, locants, and attachment label.
#[derive(Debug, Default, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(any(test, feature = "schema"), derive(JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Atom {
    /// Atomic number; zero denotes a wildcard or attachment point.
    #[cfg_attr(any(test, feature = "schema"), schemars(range(max = 118)))]
    pub element: u8,
    /// Isotope mass number. Isotopic hydrogen is an atom, never a count.
    #[cfg_attr(any(test, feature = "schema"), schemars(range(min = 1)))]
    pub isotope: Option<u16>,
    pub charge: i8,
    /// Hydrogen atoms not represented as graph vertices.
    pub hydrogens: u8,
    /// Distinguishes explicitly stated counts such as `[CH4]` from inferred `C`.
    pub hydrogens_stated: bool,
    pub radical_electrons: u8,
    pub label: Option<String>,
    /// All name locants, in name order. Parent membership controls display.
    pub locants: Vec<String>,
}

/// A localized bond; endpoint order defines double-bond references and the position-variation attachment end.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(any(test, feature = "schema"), derive(JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Bond {
    #[cfg_attr(any(test, feature = "schema"), schemars(extend("uniqueItems" = true)))]
    pub atoms: [u32; 2],
    pub order: BondOrder,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(any(test, feature = "schema"), derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
#[non_exhaustive]
pub enum BondOrder {
    Single,
    Double,
    Triple,
    Quadruple,
}

/// An ordered tetrahedral reference, including non-vertex ligands.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[cfg_attr(any(test, feature = "schema"), derive(JsonSchema))]
#[serde(
    tag = "kind",
    content = "atom",
    rename_all = "snake_case",
    deny_unknown_fields
)]
#[non_exhaustive]
pub enum Ref {
    Atom(u32),
    ImplicitHydrogen,
    LonePair,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(any(test, feature = "schema"), derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
#[non_exhaustive]
pub enum Winding {
    Clockwise,
    Anticlockwise,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(any(test, feature = "schema"), derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
#[non_exhaustive]
pub enum Relation {
    Cis,
    Trans,
}

/// Ordered structural stereochemistry, independent of CIP assignment.
///
/// `known: false` records unknown (wavy/ξ) stereochemistry. An absent entry means
/// unspecified stereochemistry. In a double-bond entry, `refs[i]` is attached
/// to `bonds[bond].atoms[i]`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(any(test, feature = "schema"), derive(JsonSchema))]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
#[non_exhaustive]
pub enum Stereo {
    Tetrahedral {
        center: u32,
        #[cfg_attr(any(test, feature = "schema"), schemars(extend("uniqueItems" = true)))]
        neighbors: [Ref; 4],
        winding: Winding,
        known: bool,
    },
    DoubleBond {
        bond: u32,
        #[cfg_attr(any(test, feature = "schema"), schemars(extend("uniqueItems" = true)))]
        refs: [u32; 2],
        relation: Relation,
        known: bool,
    },
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(any(test, feature = "schema"), derive(JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct StereoGroup {
    pub kind: StereoGroupKind,
    /// Indices into [`Molecule::stereo`], rather than atom indices.
    #[cfg_attr(any(test, feature = "schema"), schemars(length(min = 1), extend("uniqueItems" = true)))]
    pub stereo: Vec<u32>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[cfg_attr(any(test, feature = "schema"), derive(JsonSchema))]
#[serde(
    tag = "kind",
    content = "id",
    rename_all = "snake_case",
    deny_unknown_fields
)]
#[non_exhaustive]
pub enum StereoGroupKind {
    Abs,
    And(u32),
    Or(u32),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(any(test, feature = "schema"), derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
#[non_exhaustive]
pub enum Connectivity {
    HeadToTail,
    HeadToHead,
    Unknown,
}

/// A repeating unit, including its bonds to atoms outside the unit.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(any(test, feature = "schema"), derive(JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct RepeatUnit {
    #[cfg_attr(any(test, feature = "schema"), schemars(length(min = 1), extend("uniqueItems" = true)))]
    pub atoms: Vec<u32>,
    #[cfg_attr(any(test, feature = "schema"), schemars(extend("uniqueItems" = true)))]
    pub crossing_bonds: Vec<u32>,
    pub connectivity: Connectivity,
    pub label: Option<String>,
}

/// Possible attachment endpoints for a position-variation bond.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(any(test, feature = "schema"), derive(JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct PositionVariation {
    pub bond: u32,
    #[cfg_attr(any(test, feature = "schema"), schemars(length(min = 1), extend("uniqueItems" = true)))]
    pub endpoints: Vec<u32>,
}

/// A structural consistency error. Validation does not infer valence or aromaticity.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
#[non_exhaustive]
pub enum ValidationError {
    #[error("{collection} has {len} entries, exceeding the u32 identity space")]
    TooManyEntries {
        collection: &'static str,
        len: usize,
    },
    #[error("{path} references index {index}, but the collection has {len} entries")]
    IndexOutOfBounds {
        path: String,
        index: u32,
        len: usize,
    },
    #[error("atom {atom} has unsupported atomic number {element}")]
    InvalidElement { atom: usize, element: u8 },
    #[error("atom {atom} has invalid isotope mass number {isotope}")]
    InvalidIsotope { atom: usize, isotope: u16 },
    #[error("bond {bond} joins atom {atom} to itself")]
    SelfBond { bond: usize, atom: u32 },
    #[error("bond {bond} duplicates bond {previous}")]
    DuplicateBond { bond: usize, previous: usize },
    #[error("{path} contains duplicate index {index}")]
    DuplicateIndex { path: String, index: u32 },
    #[error("stereo entry {stereo} is invalid: {reason}")]
    InvalidStereo { stereo: usize, reason: String },
    #[error("stereo group {group} is invalid: {reason}")]
    InvalidStereoGroup { group: usize, reason: String },
    #[error("position variation {variation} is invalid: {reason}")]
    InvalidPositionVariation { variation: usize, reason: String },
    #[error("repeat unit {unit} is invalid: {reason}")]
    InvalidRepeatUnit { unit: usize, reason: String },
}

/// A JSON syntax/shape, unsupported-version, or structural validation error.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum JsonError {
    #[error("invalid IR JSON: {0}")]
    Json(#[from] serde_json::Error),
    #[error("unsupported nazuchem IR version {version}; supported version is {IR_VERSION}")]
    UnsupportedVersion { version: u64 },
    #[error(transparent)]
    Validation(#[from] ValidationError),
}

impl Molecule {
    /// Construct and validate a graph with empty optional annotations.
    pub fn new(atoms: Vec<Atom>, bonds: Vec<Bond>) -> Result<Self, ValidationError> {
        let molecule = Self {
            atoms,
            bonds,
            ..Self::default()
        };
        molecule.validate()?;
        Ok(molecule)
    }

    /// Validate and write a complete, versioned JSON document.
    pub fn to_json(&self) -> Result<String, JsonError> {
        self.validate()?;
        Ok(serde_json::to_string(self)?)
    }

    /// Read a complete IR document, rejecting unknown fields and variants.
    pub fn from_json(json: &str) -> Result<Self, JsonError> {
        serde_json::from_str::<MoleculeDocument>(json)?.into_molecule()
    }

    /// Check index bounds and metadata consistency without chemistry perception.
    pub fn validate(&self) -> Result<(), ValidationError> {
        for (collection, len) in [
            ("atoms", self.atoms.len()),
            ("bonds", self.bonds.len()),
            ("stereo", self.stereo.len()),
        ] {
            if len > u32::MAX as usize {
                return Err(ValidationError::TooManyEntries { collection, len });
            }
        }
        for (atom, data) in self.atoms.iter().enumerate() {
            if data.element > 118 {
                return Err(ValidationError::InvalidElement {
                    atom,
                    element: data.element,
                });
            }
            if data.isotope == Some(0) {
                return Err(ValidationError::InvalidIsotope { atom, isotope: 0 });
            }
        }
        let mut edges = BTreeMap::new();
        let mut adjacency = vec![BTreeSet::new(); self.atoms.len()];
        for (bond, data) in self.bonds.iter().enumerate() {
            for atom in data.atoms {
                check_index(atom, self.atoms.len(), format!("bonds[{bond}].atoms"))?;
            }
            let [a, b] = data.atoms;
            if a == b {
                return Err(ValidationError::SelfBond { bond, atom: a });
            }
            if let Some(previous) = edges.insert([a.min(b), a.max(b)], bond) {
                return Err(ValidationError::DuplicateBond { bond, previous });
            }
            adjacency[a as usize].insert(b);
            adjacency[b as usize].insert(a);
        }
        let mut tetra_centers = BTreeSet::new();
        let mut double_bonds = BTreeSet::new();
        for (stereo, data) in self.stereo.iter().enumerate() {
            let invalid = |reason: &str| ValidationError::InvalidStereo {
                stereo,
                reason: reason.to_owned(),
            };
            match data {
                Stereo::Tetrahedral {
                    center, neighbors, ..
                } => {
                    check_index(
                        *center,
                        self.atoms.len(),
                        format!("stereo[{stereo}].center"),
                    )?;
                    if !tetra_centers.insert(center) {
                        return Err(invalid("duplicate tetrahedral center"));
                    }
                    let mut refs = BTreeSet::new();
                    let mut atom_refs = BTreeSet::new();
                    for neighbor in neighbors {
                        if !refs.insert(neighbor) {
                            return Err(invalid("duplicate tetrahedral reference"));
                        }
                        match neighbor {
                            Ref::Atom(atom) => {
                                check_index(
                                    *atom,
                                    self.atoms.len(),
                                    format!("stereo[{stereo}].neighbors"),
                                )?;
                                if !adjacency[*center as usize].contains(atom) {
                                    return Err(invalid(
                                        "atom reference is not bonded to the center",
                                    ));
                                }
                                atom_refs.insert(*atom);
                            }
                            Ref::ImplicitHydrogen
                                if self.atoms[*center as usize].hydrogens == 0 =>
                            {
                                return Err(invalid(
                                    "implicit hydrogen reference has no hydrogen count",
                                ));
                            }
                            Ref::ImplicitHydrogen | Ref::LonePair => {}
                        }
                    }
                    if atom_refs != adjacency[*center as usize] {
                        return Err(invalid("tetrahedral references omit a bonded atom"));
                    }
                }
                Stereo::DoubleBond { bond, refs, .. } => {
                    check_index(*bond, self.bonds.len(), format!("stereo[{stereo}].bond"))?;
                    if !double_bonds.insert(bond) {
                        return Err(invalid("duplicate double-bond entry"));
                    }
                    let data = &self.bonds[*bond as usize];
                    if data.order != BondOrder::Double {
                        return Err(invalid("double-bond stereo requires a double bond"));
                    }
                    if refs[0] == refs[1] {
                        return Err(invalid("double-bond references must be distinct"));
                    }
                    for (end, atom) in refs.iter().enumerate() {
                        check_index(
                            *atom,
                            self.atoms.len(),
                            format!("stereo[{stereo}].refs[{end}]"),
                        )?;
                        if data.atoms.contains(atom)
                            || !adjacency[data.atoms[end] as usize].contains(atom)
                        {
                            return Err(invalid(
                                "double-bond reference is not a substituent of the corresponding end",
                            ));
                        }
                    }
                }
            }
        }
        let mut group_kinds = BTreeSet::new();
        let mut grouped_stereo = BTreeSet::new();
        for (group, data) in self.stereo_groups.iter().enumerate() {
            let invalid = |reason: &str| ValidationError::InvalidStereoGroup {
                group,
                reason: reason.to_owned(),
            };
            if data.stereo.is_empty() {
                return Err(invalid("group is empty"));
            }
            if !group_kinds.insert(data.kind) {
                return Err(invalid("duplicate group kind and identifier"));
            }
            check_indices(
                &data.stereo,
                self.stereo.len(),
                &format!("stereo_groups[{group}].stereo"),
            )?;
            for stereo in &data.stereo {
                if !grouped_stereo.insert(*stereo) {
                    return Err(invalid("stereo entry belongs to more than one group"));
                }
            }
        }
        let mut variation_bonds = BTreeSet::new();
        for (variation, data) in self.position_variation.iter().enumerate() {
            let invalid = |reason: &str| ValidationError::InvalidPositionVariation {
                variation,
                reason: reason.to_owned(),
            };
            check_index(
                data.bond,
                self.bonds.len(),
                format!("position_variation[{variation}].bond"),
            )?;
            if data.endpoints.is_empty() {
                return Err(invalid("endpoint list is empty"));
            }
            if !variation_bonds.insert(data.bond) {
                return Err(invalid("duplicate position-variation bond"));
            }
            check_indices(
                &data.endpoints,
                self.atoms.len(),
                &format!("position_variation[{variation}].endpoints"),
            )?;
        }
        for (unit, data) in self.polymers.iter().enumerate() {
            let invalid = |reason: &str| ValidationError::InvalidRepeatUnit {
                unit,
                reason: reason.to_owned(),
            };
            if data.atoms.is_empty() {
                return Err(invalid("atom list is empty"));
            }
            check_indices(
                &data.atoms,
                self.atoms.len(),
                &format!("polymers[{unit}].atoms"),
            )?;
            check_indices(
                &data.crossing_bonds,
                self.bonds.len(),
                &format!("polymers[{unit}].crossing_bonds"),
            )?;
            let atoms: BTreeSet<_> = data.atoms.iter().copied().collect();
            for bond in &data.crossing_bonds {
                let [a, b] = self.bonds[*bond as usize].atoms;
                if atoms.contains(&a) == atoms.contains(&b) {
                    return Err(invalid(
                        "crossing bond must have exactly one atom in the unit",
                    ));
                }
            }
        }
        check_indices(&self.parent, self.atoms.len(), "parent")?;
        Ok(())
    }
}

fn check_index(index: u32, len: usize, path: String) -> Result<(), ValidationError> {
    if index as usize >= len {
        Err(ValidationError::IndexOutOfBounds { path, index, len })
    } else {
        Ok(())
    }
}

fn check_indices(indices: &[u32], len: usize, path: &str) -> Result<(), ValidationError> {
    let mut seen = BTreeSet::new();
    for index in indices {
        check_index(*index, len, path.to_owned())?;
        if !seen.insert(index) {
            return Err(ValidationError::DuplicateIndex {
                path: path.to_owned(),
                index: *index,
            });
        }
    }
    Ok(())
}

// One owned document type provides deserialization and schema generation. The
// version stays outside the public graph so callers cannot accidentally mutate it.
#[derive(Deserialize)]
#[cfg_attr(any(test, feature = "schema"), derive(JsonSchema))]
#[serde(deny_unknown_fields)]
#[cfg_attr(any(test, feature = "schema"), schemars(rename = "Molecule"))]
struct MoleculeDocument {
    #[cfg_attr(any(test, feature = "schema"), schemars(with = "VersionOne"))]
    nazuchem_ir: u64,
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
    stereo: Vec<Stereo>,
    stereo_groups: Vec<StereoGroup>,
    position_variation: Vec<PositionVariation>,
    polymers: Vec<RepeatUnit>,
    #[cfg_attr(any(test, feature = "schema"), schemars(extend("uniqueItems" = true)))]
    parent: Vec<u32>,
}

impl MoleculeDocument {
    fn into_molecule(self) -> Result<Molecule, JsonError> {
        if self.nazuchem_ir != IR_VERSION {
            return Err(JsonError::UnsupportedVersion {
                version: self.nazuchem_ir,
            });
        }
        let molecule = Molecule {
            atoms: self.atoms,
            bonds: self.bonds,
            stereo: self.stereo,
            stereo_groups: self.stereo_groups,
            position_variation: self.position_variation,
            polymers: self.polymers,
            parent: self.parent,
        };
        molecule.validate()?;
        Ok(molecule)
    }
}

#[cfg(any(test, feature = "schema"))]
struct VersionOne;

#[cfg(any(test, feature = "schema"))]
impl JsonSchema for VersionOne {
    fn schema_name() -> Cow<'static, str> {
        "VersionOne".into()
    }
    fn inline_schema() -> bool {
        true
    }
    fn json_schema(_: &mut SchemaGenerator) -> Schema {
        json_schema!({ "type": "integer", "const": IR_VERSION })
    }
}

#[cfg(any(test, feature = "schema"))]
impl JsonSchema for Molecule {
    fn schema_name() -> Cow<'static, str> {
        "Molecule".into()
    }
    fn json_schema(generator: &mut SchemaGenerator) -> Schema {
        MoleculeDocument::json_schema(generator)
    }
}

impl Serialize for Molecule {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        self.validate().map_err(serde::ser::Error::custom)?;
        let mut document = serializer.serialize_struct("Molecule", 8)?;
        document.serialize_field("nazuchem_ir", &IR_VERSION)?;
        document.serialize_field("atoms", &self.atoms)?;
        document.serialize_field("bonds", &self.bonds)?;
        document.serialize_field("stereo", &self.stereo)?;
        document.serialize_field("stereo_groups", &self.stereo_groups)?;
        document.serialize_field("position_variation", &self.position_variation)?;
        document.serialize_field("polymers", &self.polymers)?;
        document.serialize_field("parent", &self.parent)?;
        document.end()
    }
}

impl<'de> Deserialize<'de> for Molecule {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        MoleculeDocument::deserialize(deserializer)?
            .into_molecule()
            .map_err(serde::de::Error::custom)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn published_schema_matches_the_derived_wire_document() {
        let published: serde_json::Value = serde_json::from_str(JSON_SCHEMA).unwrap();
        let generated = serde_json::to_value(json_schema()).unwrap();
        assert_eq!(published, generated);
    }
}
