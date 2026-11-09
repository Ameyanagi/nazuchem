// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | crates/nazuchem-ir/src/molecule.rs | Phase B step 8a

//! IR v1 stores atom identity and stated structure without chemistry perception.

#[cfg(any(test, feature = "schema"))]
use std::borrow::Cow;
use std::collections::{BTreeMap, BTreeSet};

#[cfg(any(test, feature = "schema"))]
use schemars::{JsonSchema, Schema, SchemaGenerator, json_schema};
use serde::de::{DeserializeOwned, MapAccess, Visitor, value::MapAccessDeserializer};
use serde::ser::SerializeStruct;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use thiserror::Error;

/// Version written in every serialized molecule.
///
/// IR v1 freezes at the first published release. Before that release, its shape
/// may change without a version bump; afterward, new fields or variants require
/// a new IR version.
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
/// Equality compares stored fields and reference order. It does not test graph
/// isomorphism or stereochemical equivalence under reordered references.
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
///
/// JSON omits zero charge, zero radical electrons, and false
/// `hydrogens_stated`, restoring those defaults when the fields are absent.
#[derive(Debug, Default, Clone, PartialEq, Eq, Serialize)]
#[cfg_attr(any(test, feature = "schema"), derive(JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Atom {
    /// Atomic number; zero denotes a wildcard or attachment point.
    #[cfg_attr(any(test, feature = "schema"), schemars(range(max = 118)))]
    pub element: u8,
    /// Isotope mass number. Isotopic hydrogen is an atom, never a count.
    #[cfg_attr(any(test, feature = "schema"), schemars(range(min = 1)))]
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub isotope: Option<u16>,
    /// Formal charge; defaults to zero when omitted from JSON.
    #[serde(default, skip_serializing_if = "is_zero")]
    pub charge: i8,
    /// Hydrogen atoms not represented as graph vertices.
    pub hydrogens: u8,
    /// Distinguishes explicitly stated counts such as `[CH4]` from inferred `C`.
    /// Defaults to false when omitted from JSON.
    #[serde(default, skip_serializing_if = "is_false")]
    pub hydrogens_stated: bool,
    /// Stated unpaired electron count; defaults to zero when omitted from JSON.
    #[serde(default, skip_serializing_if = "is_zero")]
    pub radical_electrons: u8,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
    /// An explicit SMILES atom-map number, independent of text labels.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub atom_map: Option<u32>,
    /// A numbered wildcard attachment, corresponding to OPSIN `_AP{n}`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub attachment_point: Option<u32>,
    /// All name locants, in name order. Parent membership controls display.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub locants: Vec<String>,
}

/// A localized bond; endpoint order defines double-bond references and the position-variation attachment end.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
#[cfg_attr(any(test, feature = "schema"), derive(JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Bond {
    #[cfg_attr(any(test, feature = "schema"), schemars(extend("uniqueItems" = true)))]
    pub atoms: [u32; 2],
    pub order: BondOrder,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[cfg_attr(any(test, feature = "schema"), derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
#[non_exhaustive]
pub enum BondOrder {
    Single,
    Double,
    Triple,
    Quadruple,
}

/// An ordered structural reference, including non-vertex ligands.
///
/// JSON represents atom indices as integers and the two virtual ligands as the
/// strings `"implicit_hydrogen"` and `"lone_pair"`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[non_exhaustive]
pub enum Ref {
    Atom(u32),
    ImplicitHydrogen,
    LonePair,
}

/// Tetrahedral winding, viewed from `neighbors[0]` toward the center.
///
/// `Clockwise` means `neighbors[1]`, `neighbors[2]`, `neighbors[3]` occur in
/// clockwise order in that view; `Anticlockwise` means the reverse order.
/// `Unknown` records explicitly unknown configuration. No stereo entry means
/// unspecified configuration.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[cfg_attr(any(test, feature = "schema"), derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
#[non_exhaustive]
pub enum Winding {
    Clockwise,
    Anticlockwise,
    Unknown,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[cfg_attr(any(test, feature = "schema"), derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
#[non_exhaustive]
pub enum Relation {
    Cis,
    Trans,
    Unknown,
}

/// Axis winding relative to its ordered references, independent of CIP ranking.
///
/// Look along `bond.atoms[0]` toward `bond.atoms[1]`. The first reference pair is
/// attached to the near end, the second to the far end. `Clockwise` means the
/// shortest turn from the projection of `refs[0]` to `refs[2]` is clockwise;
/// `Anticlockwise` means anticlockwise. Swapping either pair reverses winding.
/// Reversing the bond and exchanging the two pairs preserves winding. These
/// values describe the stored reference order, not CIP M/P descriptors; no
/// ligand priorities are inferred. `Unknown` records stated unknown winding.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[cfg_attr(any(test, feature = "schema"), derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
#[non_exhaustive]
pub enum AxialWinding {
    Clockwise,
    Anticlockwise,
    Unknown,
}

/// Ordered structural stereochemistry, independent of CIP assignment.
///
/// Unknown winding/relation records stated unknown (wavy/ξ) stereochemistry;
/// an absent entry means unspecified stereochemistry. For a double bond,
/// `refs[i]` is attached to `bonds[bond].atoms[i]`. For an atropisomer,
/// `refs[0..2]` are attached to the first bond end and `refs[2..4]` to the second.
/// Axial chirality of allenes and other cumulenes is not representable in v1.
/// Equality is structural: equivalent configurations with permuted references
/// and adjusted winding can compare unequal.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
#[cfg_attr(any(test, feature = "schema"), derive(JsonSchema))]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
#[non_exhaustive]
pub enum Stereo {
    Tetrahedral {
        center: u32,
        #[cfg_attr(any(test, feature = "schema"), schemars(extend("uniqueItems" = true)))]
        neighbors: [Ref; 4],
        winding: Winding,
    },
    DoubleBond {
        bond: u32,
        #[cfg_attr(any(test, feature = "schema"), schemars(extend("uniqueItems" = true)))]
        refs: [u32; 2],
        relation: Relation,
    },
    Atropisomer {
        bond: u32,
        refs: [Ref; 4],
        winding: AxialWinding,
    },
}

/// A relative or racemic group. Known ungrouped stereo is absolute.
/// Every member must have known winding or relation; stated unknown stereo
/// cannot participate in a relative or racemic group.
///
/// JSON flattens `kind` and its identifier:
/// `{"kind":"and","id":7,"stereo":[0]}`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StereoGroup {
    pub kind: StereoGroupKind,
    /// Indices into [`Molecule::stereo`], rather than atom indices.
    pub stereo: Vec<u32>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize)]
#[cfg_attr(any(test, feature = "schema"), derive(JsonSchema))]
#[serde(
    tag = "kind",
    content = "id",
    rename_all = "snake_case",
    deny_unknown_fields
)]
#[non_exhaustive]
pub enum StereoGroupKind {
    And(u32),
    Or(u32),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[cfg_attr(any(test, feature = "schema"), derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
#[non_exhaustive]
pub enum Connectivity {
    HeadToTail,
    HeadToHead,
    Unknown,
}

/// A repeating unit, including its bonds to atoms outside the unit.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
#[cfg_attr(any(test, feature = "schema"), derive(JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct RepeatUnit {
    #[cfg_attr(any(test, feature = "schema"), schemars(length(min = 1), extend("uniqueItems" = true)))]
    pub atoms: Vec<u32>,
    #[cfg_attr(any(test, feature = "schema"), schemars(extend("uniqueItems" = true)))]
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub crossing_bonds: Vec<u32>,
    /// External wildcard atom identities that represent polymer end roles.
    #[cfg_attr(any(test, feature = "schema"), schemars(extend("uniqueItems" = true)))]
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub attachment_points: Vec<u32>,
    pub connectivity: Connectivity,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
}

/// Possible attachment endpoints for a position-variation bond.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
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
    #[error("invalid IR JSON: {message}")]
    InvalidJson { message: String },
    #[error("unknown IR field {field:?} at {path}")]
    UnknownField { path: String, field: String },
    #[error("unknown IR variant {variant:?} at {path}")]
    UnknownVariant { path: String, variant: String },
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
        serde_json::to_string(&ValidatedDocument(self)).map_err(json_error)
    }

    /// Read a complete IR document, rejecting unknown fields and variants.
    pub fn from_json(json: &str) -> Result<Self, JsonError> {
        molecule_from_value(serde_json::from_str(json).map_err(json_error)?)
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
        let mut axial_bonds = BTreeSet::new();
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
                Stereo::Atropisomer { bond, refs, .. } => {
                    check_index(*bond, self.bonds.len(), format!("stereo[{stereo}].bond"))?;
                    if !axial_bonds.insert(*bond) {
                        return Err(invalid("duplicate atropisomer axis"));
                    }
                    let axis = &self.bonds[*bond as usize];
                    if axis.order != BondOrder::Single {
                        return Err(invalid("atropisomer stereo requires a single bond"));
                    }
                    for end in 0..2 {
                        let center = axis.atoms[end];
                        let pair = &refs[2 * end..2 * end + 2];
                        if pair[0] == pair[1] {
                            return Err(invalid("duplicate reference at an atropisomer axis end"));
                        }
                        let mut atom_refs = BTreeSet::new();
                        for reference in pair {
                            match reference {
                                Ref::Atom(atom) => {
                                    check_index(
                                        *atom,
                                        self.atoms.len(),
                                        format!("stereo[{stereo}].refs"),
                                    )?;
                                    if axis.atoms.contains(atom)
                                        || !adjacency[center as usize].contains(atom)
                                    {
                                        return Err(invalid(
                                            "axis reference is not a substituent of the corresponding end",
                                        ));
                                    }
                                    atom_refs.insert(*atom);
                                }
                                Ref::ImplicitHydrogen
                                    if self.atoms[center as usize].hydrogens == 0 =>
                                {
                                    return Err(invalid(
                                        "axis implicit hydrogen reference has no hydrogen count",
                                    ));
                                }
                                Ref::ImplicitHydrogen | Ref::LonePair => {}
                            }
                        }
                        let mut substituents = adjacency[center as usize].clone();
                        substituents.remove(&axis.atoms[1 - end]);
                        if atom_refs != substituents {
                            return Err(invalid(
                                "atropisomer references omit a bonded substituent",
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
                if matches!(
                    self.stereo[*stereo as usize],
                    Stereo::Tetrahedral {
                        winding: Winding::Unknown,
                        ..
                    } | Stereo::DoubleBond {
                        relation: Relation::Unknown,
                        ..
                    } | Stereo::Atropisomer {
                        winding: AxialWinding::Unknown,
                        ..
                    }
                ) {
                    return Err(invalid("group member has unknown configuration"));
                }
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
            let boundary: BTreeSet<_> = self
                .bonds
                .iter()
                .enumerate()
                .filter_map(|(bond, data)| {
                    let [a, b] = data.atoms;
                    (atoms.contains(&a) != atoms.contains(&b)).then_some(bond as u32)
                })
                .collect();
            if data.crossing_bonds.iter().copied().collect::<BTreeSet<_>>() != boundary {
                return Err(invalid(
                    "crossing bonds must be the full boundary of the unit atom set",
                ));
            }
            check_indices(
                &data.attachment_points,
                self.atoms.len(),
                &format!("polymers[{unit}].attachment_points"),
            )?;
            for atom in &data.attachment_points {
                if atoms.contains(atom)
                    || self.atoms[*atom as usize].element != 0
                    || !adjacency[*atom as usize]
                        .iter()
                        .any(|neighbor| atoms.contains(neighbor))
                {
                    return Err(invalid(
                        "attachment point must be an external wildcard bonded to the unit",
                    ));
                }
            }
        }
        check_indices(&self.parent, self.atoms.len(), "parent")?;
        Ok(())
    }
}

fn is_zero<T: Default + PartialEq>(value: &T) -> bool {
    *value == T::default()
}

fn is_false(value: &bool) -> bool {
    !*value
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

// A map visitor deliberately omits visit_seq: derived structs otherwise accept
// positional arrays, which are not part of the published object wire format.
fn deserialize_object<'de, T: DeserializeOwned, D: Deserializer<'de>>(
    deserializer: D,
) -> Result<T, D::Error> {
    struct ObjectVisitor<T>(std::marker::PhantomData<T>);
    impl<'de, T: DeserializeOwned> Visitor<'de> for ObjectVisitor<T> {
        type Value = T;
        fn expecting(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            formatter.write_str("an IR object with named fields")
        }
        fn visit_map<A: MapAccess<'de>>(self, map: A) -> Result<T, A::Error> {
            let mut value = serde_json::Value::deserialize(MapAccessDeserializer::new(map))?;
            normalize_integer_numbers(&mut value);
            T::deserialize(value).map_err(serde::de::Error::custom)
        }
    }
    deserializer.deserialize_map(ObjectVisitor(std::marker::PhantomData))
}

// JSON Schema defines integer by mathematical value, so `1.0` and `1e0`
// have the same wire meaning as `1`. Normalize only exact integral floats in
// range; retain nonintegral/out-of-range floats for the typed decoder to reject.
fn normalize_integer_numbers(value: &mut serde_json::Value) {
    match value {
        serde_json::Value::Object(object) => {
            for value in object.values_mut() {
                normalize_integer_numbers(value);
            }
        }
        serde_json::Value::Array(array) => {
            for value in array {
                normalize_integer_numbers(value);
            }
        }
        serde_json::Value::Number(number) if number.is_f64() => {
            let value = number.as_f64().expect("JSON floating point number");
            if value.fract() == 0.0 {
                if value >= 0.0 && value < u64::MAX as f64 {
                    *number = (value as u64).into();
                } else if value >= i64::MIN as f64 && value < 0.0 {
                    *number = (value as i64).into();
                }
            }
        }
        _ => {}
    }
}

macro_rules! string_enum_deserialize {
    ($name:ident { $($variant:ident => $wire:literal),* $(,)? }) => {
        impl $name {
            const WIRE_VARIANTS: &'static [&'static str] = &[$($wire),*];
        }
        impl<'de> Deserialize<'de> for $name {
            fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
                let value = String::deserialize(deserializer)?;
                match value.as_str() {
                    $($wire => Ok(Self::$variant),)*
                    _ => Err(serde::de::Error::unknown_variant(&value, Self::WIRE_VARIANTS)),
                }
            }
        }
    };
}
string_enum_deserialize!(BondOrder { Single => "single", Double => "double", Triple => "triple", Quadruple => "quadruple" });
string_enum_deserialize!(Winding { Clockwise => "clockwise", Anticlockwise => "anticlockwise", Unknown => "unknown" });
string_enum_deserialize!(Relation { Cis => "cis", Trans => "trans", Unknown => "unknown" });
string_enum_deserialize!(AxialWinding { Clockwise => "clockwise", Anticlockwise => "anticlockwise", Unknown => "unknown" });
string_enum_deserialize!(Connectivity { HeadToTail => "head_to_tail", HeadToHead => "head_to_head", Unknown => "unknown" });

macro_rules! object_deserialize {
    ($name:ident { $($(#[$attr:meta])* $field:ident: $ty:ty),* $(,)? }) => {
        impl $name {
            const WIRE_FIELDS: &'static [&'static str] = &[$(stringify!($field)),*];
        }
        impl<'de> Deserialize<'de> for $name {
            fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
                #[derive(Deserialize)]
                #[serde(deny_unknown_fields)]
                struct Wire { $($(#[$attr])* $field: $ty),* }
                let Wire { $($field),* } = deserialize_object(deserializer)?;
                Ok(Self { $($field),* })
            }
        }
    };
}

object_deserialize!(Atom {
    element: u8,
    #[serde(default)] isotope: Option<u16>,
    #[serde(default)] charge: i8,
    hydrogens: u8,
    #[serde(default)] hydrogens_stated: bool,
    #[serde(default)] radical_electrons: u8,
    #[serde(default)] label: Option<String>,
    #[serde(default)] atom_map: Option<u32>,
    #[serde(default)] attachment_point: Option<u32>,
    #[serde(default)] locants: Vec<String>,
});
object_deserialize!(Bond {
    atoms: [u32; 2],
    order: BondOrder
});
object_deserialize!(RepeatUnit {
    atoms: Vec<u32>,
    #[serde(default)] crossing_bonds: Vec<u32>,
    #[serde(default)] attachment_points: Vec<u32>,
    connectivity: Connectivity,
    #[serde(default)] label: Option<String>,
});
object_deserialize!(PositionVariation { bond: u32, endpoints: Vec<u32> });

impl Serialize for Ref {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        match self {
            Self::Atom(atom) => serializer.serialize_u32(*atom),
            Self::ImplicitHydrogen => serializer.serialize_str("implicit_hydrogen"),
            Self::LonePair => serializer.serialize_str("lone_pair"),
        }
    }
}
impl Ref {
    const WIRE_VARIANTS: &'static [&'static str] = &["implicit_hydrogen", "lone_pair"];
}
impl<'de> Deserialize<'de> for Ref {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        #[derive(Deserialize)]
        #[serde(untagged)]
        enum Wire {
            Atom(u32),
            Virtual(String),
        }
        let mut value = serde_json::Value::deserialize(deserializer)?;
        normalize_integer_numbers(&mut value);
        match Wire::deserialize(value).map_err(serde::de::Error::custom)? {
            Wire::Atom(atom) => Ok(Self::Atom(atom)),
            Wire::Virtual(value) => match value.as_str() {
                "implicit_hydrogen" => Ok(Self::ImplicitHydrogen),
                "lone_pair" => Ok(Self::LonePair),
                _ => Err(serde::de::Error::unknown_variant(
                    &value,
                    Self::WIRE_VARIANTS,
                )),
            },
        }
    }
}
#[cfg(any(test, feature = "schema"))]
impl JsonSchema for Ref {
    fn schema_name() -> Cow<'static, str> {
        "Ref".into()
    }
    fn json_schema(generator: &mut SchemaGenerator) -> Schema {
        json_schema!({"anyOf": [generator.subschema_for::<u32>(), {
            "type": "string", "enum": Self::WIRE_VARIANTS
        }]})
    }
}

impl Stereo {
    const WIRE_VARIANTS: &'static [&'static str] = &["tetrahedral", "double_bond", "atropisomer"];

    fn wire_fields(kind: &str) -> &'static [&'static str] {
        match kind {
            "tetrahedral" => &["kind", "center", "neighbors", "winding"],
            "double_bond" => &["kind", "bond", "refs", "relation"],
            "atropisomer" => &["kind", "bond", "refs", "winding"],
            _ => &["kind"],
        }
    }
}

impl<'de> Deserialize<'de> for Stereo {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        #[derive(Deserialize)]
        #[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
        enum Wire {
            Tetrahedral {
                center: u32,
                neighbors: [Ref; 4],
                winding: Winding,
            },
            DoubleBond {
                bond: u32,
                refs: [u32; 2],
                relation: Relation,
            },
            Atropisomer {
                bond: u32,
                refs: [Ref; 4],
                winding: AxialWinding,
            },
        }
        Ok(match deserialize_object(deserializer)? {
            Wire::Tetrahedral {
                center,
                neighbors,
                winding,
            } => Self::Tetrahedral {
                center,
                neighbors,
                winding,
            },
            Wire::DoubleBond {
                bond,
                refs,
                relation,
            } => Self::DoubleBond {
                bond,
                refs,
                relation,
            },
            Wire::Atropisomer {
                bond,
                refs,
                winding,
            } => Self::Atropisomer {
                bond,
                refs,
                winding,
            },
        })
    }
}

impl StereoGroupKind {
    const WIRE_VARIANTS: &'static [&'static str] = &["and", "or"];
    #[cfg(test)]
    const WIRE_FIELDS: &'static [&'static str] = &["kind", "id"];
}

impl<'de> Deserialize<'de> for StereoGroupKind {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        #[derive(Deserialize)]
        #[serde(
            tag = "kind",
            content = "id",
            rename_all = "snake_case",
            deny_unknown_fields
        )]
        enum Wire {
            And(u32),
            Or(u32),
        }
        Ok(match deserialize_object(deserializer)? {
            Wire::And(id) => Self::And(id),
            Wire::Or(id) => Self::Or(id),
        })
    }
}

#[derive(Deserialize)]
#[cfg_attr(any(test, feature = "schema"), derive(JsonSchema))]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
#[cfg_attr(any(test, feature = "schema"), schemars(rename = "StereoGroup"))]
enum StereoGroupDocument {
    And {
        id: u32,
        #[cfg_attr(any(test, feature = "schema"), schemars(length(min = 1), extend("uniqueItems" = true)))]
        stereo: Vec<u32>,
    },
    Or {
        id: u32,
        #[cfg_attr(any(test, feature = "schema"), schemars(length(min = 1), extend("uniqueItems" = true)))]
        stereo: Vec<u32>,
    },
}
impl<'de> Deserialize<'de> for StereoGroup {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        Ok(match deserialize_object(deserializer)? {
            StereoGroupDocument::And { id, stereo } => Self {
                kind: StereoGroupKind::And(id),
                stereo,
            },
            StereoGroupDocument::Or { id, stereo } => Self {
                kind: StereoGroupKind::Or(id),
                stereo,
            },
        })
    }
}
impl StereoGroup {
    const WIRE_FIELDS: &'static [&'static str] = &["kind", "id", "stereo"];
}
impl Serialize for StereoGroup {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let (kind, id) = match self.kind {
            StereoGroupKind::And(id) => ("and", id),
            StereoGroupKind::Or(id) => ("or", id),
        };
        let mut object = serializer.serialize_struct("StereoGroup", 3)?;
        object.serialize_field("kind", kind)?;
        object.serialize_field("id", &id)?;
        object.serialize_field("stereo", &self.stereo)?;
        object.end()
    }
}
#[cfg(any(test, feature = "schema"))]
impl JsonSchema for StereoGroup {
    fn schema_name() -> Cow<'static, str> {
        "StereoGroup".into()
    }
    fn json_schema(generator: &mut SchemaGenerator) -> Schema {
        StereoGroupDocument::json_schema(generator)
    }
}

// Version is checked on the raw envelope before inspecting version-specific
// fields. This remains true when a newer document introduces unfamiliar shapes.
#[derive(Deserialize)]
#[cfg_attr(any(test, feature = "schema"), derive(JsonSchema))]
#[serde(deny_unknown_fields)]
#[cfg_attr(any(test, feature = "schema"), schemars(rename = "Molecule"))]
struct MoleculeDocument {
    #[cfg_attr(any(test, feature = "schema"), schemars(with = "VersionOne"))]
    nazuchem_ir: u64,
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    stereo: Vec<Stereo>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    stereo_groups: Vec<StereoGroup>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    position_variation: Vec<PositionVariation>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    polymers: Vec<RepeatUnit>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    #[cfg_attr(any(test, feature = "schema"), schemars(extend("uniqueItems" = true)))]
    parent: Vec<u32>,
}

impl MoleculeDocument {
    const WIRE_FIELDS: &'static [&'static str] = &[
        "nazuchem_ir",
        "atoms",
        "bonds",
        "stereo",
        "stereo_groups",
        "position_variation",
        "polymers",
        "parent",
    ];

    fn into_molecule(self) -> Result<Molecule, JsonError> {
        debug_assert_eq!(self.nazuchem_ir, IR_VERSION);
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

fn json_error(error: serde_json::Error) -> JsonError {
    JsonError::InvalidJson {
        message: error.to_string(),
    }
}

fn molecule_from_value(mut value: serde_json::Value) -> Result<Molecule, JsonError> {
    normalize_integer_numbers(&mut value);
    let object = value.as_object().ok_or_else(|| JsonError::InvalidJson {
        message: "molecule must be an object with named fields".into(),
    })?;
    let version = object
        .get("nazuchem_ir")
        .and_then(serde_json::Value::as_u64)
        .ok_or_else(|| JsonError::InvalidJson {
            message: "nazuchem_ir must be an unsigned integer".into(),
        })?;
    if version != IR_VERSION {
        return Err(JsonError::UnsupportedVersion { version });
    }
    check_document_fields(&value)?;
    serde_json::from_value::<MoleculeDocument>(value)
        .map_err(json_error)?
        .into_molecule()
}

fn check_fields(value: &serde_json::Value, allowed: &[&str], path: &str) -> Result<(), JsonError> {
    let object = value.as_object().ok_or_else(|| JsonError::InvalidJson {
        message: format!("{path} must be an object with named fields"),
    })?;
    for field in object.keys() {
        if !allowed.contains(&field.as_str()) {
            return Err(JsonError::UnknownField {
                path: path.into(),
                field: field.clone(),
            });
        }
    }
    Ok(())
}

fn check_variant(
    value: Option<&serde_json::Value>,
    allowed: &[&str],
    path: &str,
) -> Result<(), JsonError> {
    if let Some(variant) = value.and_then(serde_json::Value::as_str)
        && !allowed.contains(&variant)
    {
        return Err(JsonError::UnknownVariant {
            path: path.into(),
            variant: variant.into(),
        });
    }
    Ok(())
}

fn check_document_fields(value: &serde_json::Value) -> Result<(), JsonError> {
    check_fields(value, MoleculeDocument::WIRE_FIELDS, "$")?;
    for collection in [
        "atoms",
        "bonds",
        "stereo",
        "stereo_groups",
        "position_variation",
        "polymers",
    ] {
        for (index, item) in value
            .get(collection)
            .and_then(serde_json::Value::as_array)
            .into_iter()
            .flatten()
            .enumerate()
        {
            let path = format!("$.{collection}[{index}]");
            let allowed: &[&str] = match collection {
                "atoms" => Atom::WIRE_FIELDS,
                "bonds" => Bond::WIRE_FIELDS,
                "stereo" => {
                    check_variant(
                        item.get("kind"),
                        Stereo::WIRE_VARIANTS,
                        &format!("{path}.kind"),
                    )?;
                    Stereo::wire_fields(
                        item.get("kind")
                            .and_then(serde_json::Value::as_str)
                            .unwrap_or(""),
                    )
                }
                "stereo_groups" => StereoGroup::WIRE_FIELDS,
                "position_variation" => PositionVariation::WIRE_FIELDS,
                "polymers" => RepeatUnit::WIRE_FIELDS,
                _ => unreachable!(),
            };
            check_fields(item, allowed, &path)?;
            match collection {
                "bonds" => check_variant(
                    item.get("order"),
                    BondOrder::WIRE_VARIANTS,
                    &format!("{path}.order"),
                )?,
                "stereo" => {
                    match item.get("kind").and_then(serde_json::Value::as_str) {
                        Some("tetrahedral") => check_variant(
                            item.get("winding"),
                            Winding::WIRE_VARIANTS,
                            &format!("{path}.winding"),
                        )?,
                        Some("double_bond") => check_variant(
                            item.get("relation"),
                            Relation::WIRE_VARIANTS,
                            &format!("{path}.relation"),
                        )?,
                        Some("atropisomer") => check_variant(
                            item.get("winding"),
                            AxialWinding::WIRE_VARIANTS,
                            &format!("{path}.winding"),
                        )?,
                        _ => {}
                    }
                    for field in ["neighbors", "refs"] {
                        if item.get("kind").and_then(serde_json::Value::as_str)
                            == Some("double_bond")
                        {
                            continue;
                        }
                        for (reference, value) in item
                            .get(field)
                            .and_then(serde_json::Value::as_array)
                            .into_iter()
                            .flatten()
                            .enumerate()
                        {
                            check_variant(
                                Some(value),
                                Ref::WIRE_VARIANTS,
                                &format!("{path}.{field}[{reference}]"),
                            )?;
                        }
                    }
                }
                "stereo_groups" => check_variant(
                    item.get("kind"),
                    StereoGroupKind::WIRE_VARIANTS,
                    &format!("{path}.kind"),
                )?,
                "polymers" => check_variant(
                    item.get("connectivity"),
                    Connectivity::WIRE_VARIANTS,
                    &format!("{path}.connectivity"),
                )?,
                _ => {}
            }
        }
    }
    Ok(())
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

// This wrapper is only created after validation. Both public serialization
// paths validate exactly once while to_json retains its typed validation error.
struct ValidatedDocument<'a>(&'a Molecule);
impl Serialize for ValidatedDocument<'_> {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let molecule = self.0;
        let optional_count = [
            !molecule.stereo.is_empty(),
            !molecule.stereo_groups.is_empty(),
            !molecule.position_variation.is_empty(),
            !molecule.polymers.is_empty(),
            !molecule.parent.is_empty(),
        ]
        .into_iter()
        .filter(|present| *present)
        .count();
        let mut document = serializer.serialize_struct("Molecule", 3 + optional_count)?;
        document.serialize_field("nazuchem_ir", &IR_VERSION)?;
        document.serialize_field("atoms", &molecule.atoms)?;
        document.serialize_field("bonds", &molecule.bonds)?;
        if !molecule.stereo.is_empty() {
            document.serialize_field("stereo", &molecule.stereo)?;
        }
        if !molecule.stereo_groups.is_empty() {
            document.serialize_field("stereo_groups", &molecule.stereo_groups)?;
        }
        if !molecule.position_variation.is_empty() {
            document.serialize_field("position_variation", &molecule.position_variation)?;
        }
        if !molecule.polymers.is_empty() {
            document.serialize_field("polymers", &molecule.polymers)?;
        }
        if !molecule.parent.is_empty() {
            document.serialize_field("parent", &molecule.parent)?;
        }
        document.end()
    }
}
impl Serialize for Molecule {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        self.validate().map_err(serde::ser::Error::custom)?;
        ValidatedDocument(self).serialize(serializer)
    }
}
impl<'de> Deserialize<'de> for Molecule {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        molecule_from_value(serde_json::Value::deserialize(deserializer)?)
            .map_err(serde::de::Error::custom)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_round_trip<T>(value: &T) -> serde_json::Value
    where
        T: Serialize + DeserializeOwned + PartialEq + std::fmt::Debug,
    {
        let json = serde_json::to_value(value).unwrap();
        assert_eq!(&serde_json::from_value::<T>(json.clone()).unwrap(), value);
        json
    }

    fn assert_keys(value: &serde_json::Value, fields: &[&str]) {
        let serialized: BTreeSet<_> = value
            .as_object()
            .unwrap()
            .keys()
            .map(String::as_str)
            .collect();
        assert_eq!(serialized, fields.iter().copied().collect());
    }

    #[test]
    fn every_wire_type_and_variant_agrees_with_decoder_lists() {
        // Explicit fields keep this fixture fully populated when a record grows.
        // Exhaustive enum matches below likewise require every new variant to
        // have a serialization fixture before its decoder lists can drift.
        let atom = Atom {
            element: 6,
            isotope: Some(13),
            charge: -1,
            hydrogens: 1,
            hydrogens_stated: true,
            radical_electrons: 1,
            label: Some("C".into()),
            atom_map: Some(12),
            attachment_point: Some(2),
            locants: vec!["1".into()],
        };
        assert_keys(&assert_round_trip(&atom), Atom::WIRE_FIELDS);

        let bond = Bond {
            atoms: [0, 1],
            order: BondOrder::Single,
        };
        assert_keys(&assert_round_trip(&bond), Bond::WIRE_FIELDS);
        let stereo = [
            Stereo::Tetrahedral {
                center: 0,
                neighbors: [
                    Ref::Atom(1),
                    Ref::Atom(2),
                    Ref::Atom(3),
                    Ref::ImplicitHydrogen,
                ],
                winding: Winding::Clockwise,
            },
            Stereo::DoubleBond {
                bond: 0,
                refs: [2, 3],
                relation: Relation::Cis,
            },
            Stereo::Atropisomer {
                bond: 0,
                refs: [
                    Ref::Atom(2),
                    Ref::ImplicitHydrogen,
                    Ref::Atom(3),
                    Ref::LonePair,
                ],
                winding: AxialWinding::Anticlockwise,
            },
        ];
        let mut stereo_variants = BTreeSet::new();
        for entry in &stereo {
            match entry {
                Stereo::Tetrahedral { .. }
                | Stereo::DoubleBond { .. }
                | Stereo::Atropisomer { .. } => {}
            }
            let value = assert_round_trip(entry);
            let kind = value["kind"].as_str().unwrap();
            assert_keys(&value, Stereo::wire_fields(kind));
            stereo_variants.insert(kind.to_owned());
        }
        assert_eq!(
            stereo_variants,
            Stereo::WIRE_VARIANTS
                .iter()
                .map(|v| (*v).to_owned())
                .collect()
        );

        let groups = [StereoGroupKind::And(1), StereoGroupKind::Or(2)].map(|kind| StereoGroup {
            kind,
            stereo: vec![0],
        });
        let mut group_variants = BTreeSet::new();
        for group in &groups {
            match group.kind {
                StereoGroupKind::And(_) | StereoGroupKind::Or(_) => {}
            }
            let value = assert_round_trip(group);
            assert_keys(&value, StereoGroup::WIRE_FIELDS);
            let kind = assert_round_trip(&group.kind);
            assert_keys(&kind, StereoGroupKind::WIRE_FIELDS);
            group_variants.insert(value["kind"].as_str().unwrap().to_owned());
        }
        assert_eq!(
            group_variants,
            StereoGroupKind::WIRE_VARIANTS
                .iter()
                .map(|v| (*v).to_owned())
                .collect()
        );

        let polymer = RepeatUnit {
            atoms: vec![0, 2, 3],
            crossing_bonds: vec![0],
            attachment_points: vec![1],
            connectivity: Connectivity::HeadToTail,
            label: Some("n".into()),
        };
        assert_keys(&assert_round_trip(&polymer), RepeatUnit::WIRE_FIELDS);
        let variation = PositionVariation {
            bond: 0,
            endpoints: vec![2, 3],
        };
        assert_keys(
            &assert_round_trip(&variation),
            PositionVariation::WIRE_FIELDS,
        );
        let molecule = Molecule {
            atoms: vec![
                atom,
                Atom::default(),
                Atom {
                    element: 9,
                    ..Atom::default()
                },
                Atom {
                    element: 17,
                    ..Atom::default()
                },
            ],
            bonds: (1..=3)
                .map(|atom| Bond {
                    atoms: [0, atom],
                    order: BondOrder::Single,
                })
                .collect(),
            stereo: vec![stereo[0].clone()],
            stereo_groups: vec![groups[0].clone()],
            position_variation: vec![variation],
            polymers: vec![polymer],
            parent: vec![0],
        };
        let document = assert_round_trip(&molecule);
        assert_keys(&document, MoleculeDocument::WIRE_FIELDS);
        check_document_fields(&document).unwrap();

        macro_rules! assert_string_variants {
            ($name:ident { $($variant:ident),+ $(,)? }) => {{
                let mut serialized = BTreeSet::new();
                for value in [$($name::$variant),+] {
                    match value { $($name::$variant => {}),+ }
                    serialized.insert(assert_round_trip(&value).as_str().unwrap().to_owned());
                }
                assert_eq!(serialized, $name::WIRE_VARIANTS.iter().map(|v| (*v).to_owned()).collect());
            }};
        }
        assert_string_variants!(BondOrder {
            Single,
            Double,
            Triple,
            Quadruple
        });
        assert_string_variants!(Winding {
            Clockwise,
            Anticlockwise,
            Unknown
        });
        assert_string_variants!(Relation {
            Cis,
            Trans,
            Unknown
        });
        assert_string_variants!(AxialWinding {
            Clockwise,
            Anticlockwise,
            Unknown
        });
        assert_string_variants!(Connectivity {
            HeadToTail,
            HeadToHead,
            Unknown
        });
        let mut virtual_refs = BTreeSet::new();
        for reference in [Ref::Atom(1), Ref::ImplicitHydrogen, Ref::LonePair] {
            match reference {
                Ref::Atom(_) | Ref::ImplicitHydrogen | Ref::LonePair => {}
            }
            if let Some(name) = assert_round_trip(&reference).as_str() {
                virtual_refs.insert(name.to_owned());
            }
        }
        assert_eq!(
            virtual_refs,
            Ref::WIRE_VARIANTS.iter().map(|v| (*v).to_owned()).collect()
        );
    }

    #[test]
    fn published_schema_matches_the_derived_wire_document() {
        let generated = serde_json::to_value(json_schema()).unwrap();
        if std::env::var_os("UPDATE_NAZUCHEM_IR_SCHEMA").is_some() {
            std::fs::write(
                concat!(
                    env!("CARGO_MANIFEST_DIR"),
                    "/schema/nazuchem-ir-1.schema.json"
                ),
                format!("{}\n", serde_json::to_string_pretty(&generated).unwrap()),
            )
            .unwrap();
        }
        let published: serde_json::Value = serde_json::from_str(
            &std::fs::read_to_string(concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/schema/nazuchem-ir-1.schema.json"
            ))
            .unwrap(),
        )
        .unwrap();
        assert_eq!(published, generated);
    }
}
