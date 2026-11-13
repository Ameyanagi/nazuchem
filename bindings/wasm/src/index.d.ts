// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | bindings/wasm/src/index.d.ts | Phase B step 9

/** Rust serde vocabulary; preference is independent of identity verification. */
export type NameReason =
  | "sulfonyl_prefix" | "amino_prefix" | "cation_prefix" | "functional_class"
  | "hydride_parent" | "unretained_parent" | "diamide_locants" | "fusion_component"
  | "unregistered_retained_component" | "higher_order_fusion" | "benzo_fusion"
  | "multiplicative" | "acyl_prefix" | "ring_assembly" | "von_baeyer_fusion_fallback"
  | "linear_phane" | "opsin_fusion_layout" | "fusion_multiparent" | "hydro_prefix_locants"
  | "bridged_fused_parent" | "unaudited_parent" | "unaudited_path"
  | "substituted_retained_parent" | "atom_order" | "amine_chalcogenide"
  | "amine_oxide_zwitterion" | "principal_group_tie" | "lambda_parent_hydrogen"
  | "spiro_hydrogen" | "isotope_placement" | "boughton_isotope" | "radical_ion"
  | "anion_parent" | "cation_parent" | "fusion_nitrogen_hydron" | "ionic_suffix_order"
  | "salt_ion_order"
  | "adduct" | "phane" | "relative_ring_stereo" | "endo_exo_stereo" | "unknown_configuration"
  | "axial_alternative";

export type NameStatus = "preferred" | "not_preferred" | "unverified";

export interface AlternativeName {
  name: string;
  status: NameStatus;
  reason: NameReason | null;
}

export interface Verification {
  name: string;
  smiles: string;
}

/** Complete naming result; Rust's tagged status is flattened as status/reason. */
export interface Name {
  text: string;
  status: NameStatus;
  reason: NameReason | null;
  alternatives: AlternativeName[];
  omittedStereo: string[];
  omittedUnits: number;
  verification: Verification | null;
}

/** Integer ranges and graph constraints are enforced by the Rust IR reader. */
export interface Atom {
  element: number;
  hydrogens: number;
  isotope?: number | null;
  charge?: number;
  hydrogens_stated?: boolean;
  radical_electrons?: number;
  label?: string | null;
  atom_map?: number | null;
  attachment_point?: number | null;
  locants?: string[];
}

export type BondOrder = "single" | "double" | "triple" | "quadruple";
export interface Bond { atoms: [number, number]; order: BondOrder }
export type Ref = number | "implicit_hydrogen" | "lone_pair";
export type Winding = "clockwise" | "anticlockwise" | "unknown";
export type AxialWinding = "clockwise" | "anticlockwise" | "unknown";
export type Relation = "cis" | "trans" | "unknown";

/** Ordered structural stereo; these windings do not encode CIP labels. */
export type Stereo =
  | { kind: "tetrahedral"; center: number; neighbors: [Ref, Ref, Ref, Ref]; winding: Winding }
  | { kind: "double_bond"; bond: number; refs: [number, number]; relation: Relation }
  | { kind: "atropisomer"; bond: number; refs: [Ref, Ref, Ref, Ref]; winding: AxialWinding };

export type StereoGroup =
  | { kind: "and"; id: number; stereo: [number, ...number[]] }
  | { kind: "or"; id: number; stereo: [number, ...number[]] };

export type Connectivity = "head_to_tail" | "head_to_head" | "unknown";
export interface RepeatUnit {
  atoms: [number, ...number[]];
  connectivity: Connectivity;
  crossing_bonds?: number[];
  attachment_points?: number[];
  label?: string | null;
}
export interface PositionVariation { bond: number; endpoints: [number, ...number[]] }

/** Exact IR v1 JSON shape. Optional fields take the Rust defaults when absent. */
export interface Molecule {
  nazuchem_ir: 1;
  atoms: Atom[];
  bonds: Bond[];
  stereo?: Stereo[];
  stereo_groups?: StereoGroup[];
  position_variation?: PositionVariation[];
  polymers?: RepeatUnit[];
  parent?: number[];
}

export type ErrorKind = "init" | "parse" | "warning" | "serialize" | "naming" | "internal";
export class NazuchemError extends Error {
  readonly kind: ErrorKind;
  constructor(message: string, kind: ErrorKind);
}

export const IR_VERSION: 1;
/** Published draft-2020-12 schema as a JSON string, identical to Rust JSON_SCHEMA. */
export const JSON_SCHEMA: string;
export function nameToSmiles(name: string): string;
export function smilesToName(smiles: string): Name;
export function nameToMolecule(name: string): Molecule;
export function moleculeToName(molecule: Molecule | string): Name;
