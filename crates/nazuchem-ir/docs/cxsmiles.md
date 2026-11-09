<!-- SPDX-License-Identifier: MIT OR Apache-2.0 -->
<!-- Origin: original Nazuchem code -->

# IR v1 CXSMILES export

`Molecule::to_cxsmiles()` produces a deterministic convenience view of a
validated IR molecule. **Versioned IR JSON is the lossless contract.** CXSMILES
retains the chemical graph and supported annotations, together with complete
locants and principal-parent membership. It does not preserve every IR field,
identity, vector order, or distinction between absent and empty annotations.
There is no CXSMILES importer in this phase.

The writer performs no aromaticity, valence, CIP, or other chemical perception.
Atom and bond insertion order determine traversal, so the result is not a
canonical identifier. Components containing wildcard atoms precede ordinary
components. Explicit hydrogen vertices, including isotopic hydrogen, remain
vertices. An empty molecule exports as the empty string.

The extensions follow the [ChemAxon CXSMILES format](https://docs.chemaxon.com/latest/formats_chemaxon-extended-smiles-and-smarts-cxsmiles-and-cxsmarts.html).
The writer uses `$` for quadruple bonds and `%(number)` for ring closure numbers
above 99, as supported by the [RDKit SMILES extensions](https://www.rdkit.org/docs/RDKit_Book.html#smiles-support-and-extensions).

| IR information | Exported representation |
| --- | --- |
| Elements, isotopes, charge, hydrogen counts | SMILES atom tokens with localized Kekulé bond orders. Stated hydrogen counts retain brackets. |
| Atom maps | SMILES `:number` inside brackets, from `Atom.atom_map`. |
| Actual labels | `$label;...$`, at emitted atom indices. A homology label may include its terminal `_p` marker. |
| Attachment roles | `Atom.attachment_point` supplies a reserved `_APnumber` label; `RepeatUnit.attachment_points` supplies `star_e` labels. An actual label takes precedence. These generated display labels do not alter `Atom.label`. |
| First locants | `$_AV:value;...$` for every atom with a first locant, including atoms outside the parent. |
| Parent membership and complete locants | The two custom properties described below. |
| Radical electrons 1, 2, 3 | `^1:`, `^2:`, `^5:` respectively. Spin multiplicity is not inferred. Other counts have no dedicated exported radical layer. |
| Known tetrahedral stereo | `@` or `@@`, adjusted for the emitted neighbour order. An implicit hydrogen or lone pair occupies the virtual ligand slot. |
| Known double-bond stereo | `/` and `\` on adjacent single bonds, with parity constraints solved across conjugated systems. |
| Unknown tetrahedral stereo | No `@` or `@@`; `w:atom.bond` on a distinct incident single bond without a slash for each centre. A record without an assignable carrier returns a typed error. |
| Unknown double-bond stereo | No asserted slash configuration; `w:atom.bond` and `ctu:bond`. |
| Known axial stereo | A typed `CxsmilesError::StereoNotExpressible` error. No portable axial encoding is implemented. |
| Unknown axial stereo | The same typed error; omitting the record would lose its stated unknown state. Its record remains available in IR JSON. |
| Enhanced stereo | `&number:` or `onumber:` for groups consisting entirely of known tetrahedral centres. Other group members return a typed error. Known ungrouped entries are absolute. |
| Position variation | `m:attachment:endpoint.endpoint`, with the attachment source taken from `bonds[bond].atoms[0]`. |
| Repeat units | `Sg:n:atoms:label:ht`, `hh`, or `eu`. The complete crossing-bond boundary and ordered attachment-point roles remain in JSON. |

All standard indices refer to emitted atom order and CX bond order. A bond
receives its CX index when the later endpoint is written, including a ring bond
at its closing marker. `smiles_atom_order()` returns the corresponding IR atom
identities so callers can align external annotations without guessing.
Extension values escape non-alphanumeric UTF-16 code units with numeric entities,
following the established OPSIN convention. Reserved attachment labels and the
terminal homology `_p` marker retain their standard spelling.

## Locant properties

The writer emits exactly two custom property names using standard
`atomProp:index.key.value` syntax:

| Key | Value |
| --- | --- |
| `nazuchemParent` | `1` on each principal-parent atom; absence means false. |
| `nazuchemLocants` | A compact JSON string array containing all locants in name order, escaped with numeric entities. Emitted for each nonempty locant vector, including a vector whose first locant is an empty string. |

No topology document, atom-index reconstruction layer, or additional `nazuchem*`
properties are embedded. Consumers needing stable identities, exact ordered
stereo references, hydrogen origin, unsupported radicals, or full annotation
vectors must use IR JSON.

`displayed_locants()` returns the first locant for atoms in `parent`, and `None`
for other atoms, in IR identity order. The CX atom-value layer also includes
non-parent first locants; it is a general annotation layer rather than the
parent-only display projection.

## Stereo that cannot be expressed

The writer must never assert a configuration absent from the IR. Slash pruning
retains a mark at each end of every known double bond and leaves at least one
end unmarked for each absent or unknown double-bond record. If local pruning
cannot meet those conditions, a Boolean constraint solver checks whether any
subset of the available marks can do so. If no subset works, export returns
`CxsmilesError::StereoNotExpressible`.

For example, `(2E,6E)-octa-2,4,6-triene` has specified outer double bonds and an
unspecified middle bond. The slashes required by the outer configurations would
also assert a middle configuration, so this molecule returns the typed error.
Adding custom metadata would not undo that extra standard SMILES assertion.

The same error covers every axial state, known tetrahedral records requiring
more than one virtual ligand, unknown tetrahedral records without a usable
wiggly single bond, double-bond reference bonds that cannot carry slashes, and
enhanced groups containing any member other than a known tetrahedral record.
Unknown centres receive globally distinct wiggly carriers through augmenting
paths in a centre-to-bond matching. If their eligible bonds cannot supply one
distinct carrier per centre, export returns the typed error. Group membership
is never silently dropped or reduced.
Inconsistent slash parity constraints return
`CxsmilesError::ContradictoryStereo`. Invalid graphs return
`CxsmilesError::InvalidStructure`. The older spelling `CxSmilesError` is a type
alias for source compatibility. No JSON implementation error type appears in
this error API.

## Internal compatibility support

The opt-in `test-support` feature exposes the hidden
`semantic_cxsmiles_for_tests()` method for the independent legacy oracle. This
helper omits locant properties, stated-hydrogen brackets, explicit radical
layers, and unknown-stereo layers. It retains the same fail-closed stereo
checks as the public writer. A sole AND group 1 or 2 uses legacy `r`; multiple
AND and OR groups are numbered sequentially. Polymer attachment points receive
sequential atom maps in their recorded attachment-point order to match the
legacy writer. Other atom maps are retained.

The helper folds ordinary terminal hydrogen vertices into counts only when they
have no locants, parent membership, isotope, charge, radical, label, map, or
attachment role. A tetrahedral atom reference may become the single implicit-H
slot; H vertices referenced by double-bond or axial stereo, position variation,
or repeat-unit annotations remain explicit. This temporary projection changes
neither the source molecule nor its full CXSMILES export.
`semantic_cxsmiles_with_atom_order_for_tests()` also returns the original IR atom
identities in its emitted order; this is a subset when hydrogens are folded.
The public `smiles_atom_order()` continues to describe the full projection.

This test helper is not the application API or a lossless serialization format.
