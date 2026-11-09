<!-- SPDX-License-Identifier: MIT OR Apache-2.0 -->
<!-- Origin: original Nazuchem code -->

# IR v1 CXSMILES export

`Molecule::to_cxsmiles()` is a deterministic, lossless projection of a validated
IR molecule. It performs no aromaticity, valence, CIP, or other chemical
perception. Atom and bond insertion order determine the traversal; this is not
a canonical identifier. Disconnected components containing wildcards precede
other components. Isotopic and ordinary hydrogen vertices remain vertices.

Losslessness means that the standard SMILES graph together with the documented
`nazuchem*` properties contains every IR field, including vector order and atom
identity. A reader must preserve and apply these properties. A consumer that
discards arbitrary atom properties retains only the standard chemical view.
This phase does not provide a CXSMILES reader. Empty IR molecules export as the
empty string; validation requires their annotations to be empty too.

The extension syntax follows [ChemAxon's CXSMILES documentation](https://docs.chemaxon.com/latest/formats_chemaxon-extended-smiles-and-smarts-cxsmiles-and-cxsmarts.html).
Quadruple bonds use `$`, and more than 100 simultaneous ring closures use
`%(number)`, as described by the [RDKit Book](https://www.rdkit.org/docs/RDKit_Book.html).
The writer uses numeric entities for non-alphanumeric UTF-16 code units in
extension values, consistent with the existing OPSIN label writer. Reserved
attachment labels `_AP` followed by digits and `star_e` retain their standard
spelling; homology labels retain their terminal `_p` marker.

| IR information | Standard chemical projection |
| --- | --- |
| Element, isotope, charge, non-vertex hydrogen count | SMILES atom tokens; localized bonds retain their Kekulé order. Brackets retain a stated hydrogen count. |
| Atom label | `$label;...$`, at the emitted atom's index. The exact optional label is also stored as a JSON property. |
| First locant | `$_AV:value;...$`, for every atom with a locant, including atoms outside the principal parent. |
| Parent membership | `atomProp:i.nazuchemParent.1`. |
| Radical electrons 1, 2, 3 | `^1:`, `^2:`, `^5:` respectively. Spin multiplicity is not in IR v1 and is not inferred. The exact count, including zero, also has a property to prevent inference from valence. |
| Known tetrahedral stereo | `@`/`@@`, adjusted for emitted neighbour order. A lone pair occupies the virtual ligand slot. Records with two virtual ligands have only the metadata representation. |
| Known double-bond stereo | Adjacent `/` and `\` single bonds, solved as parity constraints across conjugated systems. Records without single reference bonds have only the metadata representation. |
| Unknown tetrahedral stereo | `w:atom.bond` on an incident single bond when available; no `@`/`@@` assertion. |
| Unknown double-bond stereo | `w:atom.bond` and `ctu:bond`; no slash assertion for that record. |
| Enhanced stereo | `a:`, `&n:`, `on:` for tetrahedral entries. Standard atom-based groups cannot encode group references to double-bond records; metadata retains them. |
| Position variation | `m:attachment:endpoint.endpoint`. The attachment source is always `bonds[bond].atoms[0]`, independent of either endpoint's element; adapters must orient the bond with its source first. The exact designated bond and ordered endpoints remain in metadata. |
| Repeat units | `Sg:n:atoms:label:ht`, `hh`, or `eu`. The full list of crossing bonds remains in metadata: IR v1 does not identify a head/tail bracket partition, so that partition is not invented. |

All standard indices refer to emitted SMILES atom order and CX bond order. A CX
bond receives its index when the later endpoint is written, including a ring
bond at its closing marker. The exact IR bond identities are independent of
these output indices.

The standard projection of structurally valid but chemically contradictory
stereo records cannot assert those records simultaneously. In that case the
full exporter suppresses slash assertions and retains the exact records in the
metadata. This does not change the IR or infer a replacement configuration.

## Versioned atom properties

Properties use standard `atomProp:index.key.value` syntax. Their values are
strings. JSON-valued properties are compact JSON escaped with the same numeric
entity convention. The following property contract is version 1:

| Key | Location and meaning |
| --- | --- |
| `nazuchemVersion` | Emitted atom 0, value `1`. Marks this metadata contract for every nonempty molecule. |
| `nazuchemIndex` | Every emitted atom when the traversal changes atom order; the value is the original IR atom index. If absent throughout, emitted and IR atom orders are identical. |
| `nazuchemParent` | Value `1` on principal-parent atoms; absence means false. This exposes parent membership to ordinary property-aware CX consumers. |
| `nazuchemLocants` | JSON string array containing **all** locants, in name order, on every atom with a nonempty locant vector. It preserves empty strings and distinguishes an empty vector from a first empty locant. |
| `nazuchemHydrogensStated` | Value `0` on atoms whose necessary brackets would otherwise falsely imply a stated hydrogen count. Otherwise the bracket state gives this flag: bracketed means true, unbracketed means false. |
| `nazuchemHydrogens` | Exact non-vertex hydrogen count on every atom, including hydrogen and wildcard atoms, so restoration never requires valence inference. |
| `nazuchemRadicalElectrons` | Exact radical-electron count on every atom, including zero, so restoration never infers radicals from under-valent bracket notation. |
| `nazuchemLabel` | JSON string for every present label, including the empty string; absence means `None`. This preserves the original label independently of reserved standard aliases. |
| `nazuchemTopology` | Emitted atom 0. A JSON object containing each nonempty graph annotation array listed below; omitted arrays are empty. |

`nazuchemTopology` contains `bonds`, `stereo`, `stereo_groups`,
`position_variation`, `polymers`, and `parent`, serialized in the IR v1 JSON
representation. All indices in this object are original IR identities. It
retains bond vector order and endpoint orientation, tetrahedral reference order
and virtual-reference types, double-bond reference order and relation, unknown
stereo winding, stereo-entry identities and group numbers/order, exact
position-variation bond identity, SRU crossing-bond identities and order,
optional SRU labels, and principal-parent vector order. Standard CXSMILES does
not encode these distinctions consistently. Atom data already represented by
the standard graph or per-atom properties is not duplicated here.

A future lossless importer must restore atom order using `nazuchemIndex`, apply
per-atom metadata, restore the annotation vectors from `nazuchemTopology`, and
validate the resulting molecule. It must treat a contradiction between the
standard projection and the exact metadata as an import error rather than
silently choosing one representation. Metadata with an unknown version must be
rejected with a typed error.

## Compatibility projection and locant display

`Molecule::semantic_cxsmiles()` exports the existing OPSIN semantic view. It
excludes locants, IR properties, hydrogen-origin brackets, and explicit radical
layers. As in the legacy writer, an unknown tetrahedral record retains its
stored winding in this compatibility projection; callers needing explicit
unknown stereo must use `to_cxsmiles()`. A sole AND group 1 or 2 uses legacy `r`;
multiple OR/AND groups are numbered sequentially in output order. Absolute
groups are implicit. Polymer `star_e` atoms receive sequential atom maps from
the crossing-bond order to reproduce the legacy polymer writer.

`displayed_locants()` returns one entry per IR atom: its first locant if it is
in `parent`, otherwise `None`. `smiles_atom_order()` exposes the emitted atom
identities; it can project this vector into serialized order without losing
identity. The full export also writes non-parent locants, so downstream uses do
not have to recover them from display labels.

The semantic projection is not an exact serialization of the IR contract.
Non-polymer SMILES atom-map provenance is absent from IR v1; a label `_AP1` does
not distinguish an explicit `[*:1]` from an automatically numbered bare `*`.
Graph-identity comparisons may therefore be needed when gating legacy text.
