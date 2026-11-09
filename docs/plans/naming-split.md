# Plan: move structure↔name out of ReShiki into an open-source library

Status: accepted 2026-11-09 (decisions D1–D6 made). Name: **nazuchem** (名付け + chem).

## 1. Where things stand (main @ e73e0fdd)

| Part                            | Location                                                                                                                                                                                                                                                      | Size                                                                           | Origin                                                                                                                                                                                                                                          |
| ------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| name→structure (n2s)            | `crates/opsin`                                                                                                                                                                                                                                                | ~50k lines Rust (src + tests), 5.7 MB resources, 5.7 MB oracle fixtures        | Port of OPSIN main `42e87ac6` (MIT, Daniel Lowe). Already a standalone crate with `PORTING.md`, `NOTICE`, `provenance/`                                                                                                                         |
| structure→name (s2n)            | `src/naming/openclatura`                                                                                                                                                                                                                                      | 103 Rust files / 54k lines (nested `fusion/`, `phane/`, `rules/`), 960 KB data | Began as a port of openclatura 0.4.3 (MIT, lamalab-org), heavily extended in ReShiki. **Also contains OPSIN-derived data** (`data/fusion_entry_geometry.json`, `retained_macrocycle_templates.json`). Many files have no per-file origin marker |
| App glue                        | `src/naming.rs`, `src/naming/*.rs`                                                                                                                                                                                                                            | ~3.8k lines                                                                    | ReShiki original; `catalog` data is OPSIN MIT + ChEBI CC BY 4.0                                                                                                                                                                                 |
| Chemistry the naming code needs | `crates/chemistry`: `smiles`, `kekulize`, `aromaticity`, `rings`, `ranking`, `stereo` (perception, CIP, and the spatial/drawing types they use), `sanitize`, `electronic`, `normalize`, `cx`, `graph`, `hydrogens`, `atomic_data`, `isotopes`, `native_order` | ~27.5k lines (exact transitive closure to be inventoried in step 1)            | RDKit-derived (BSD-3); extended-tetrahedral CIP from SiMolecule **centres (BSD-2)**                                                                                                                                                             |

Couplings that matter:

- `opsin` → `reshiki-chemistry` only in `modern_cip.rs`. Its result reaches ReShiki as semantic CXSMILES + a parent-only locant projection (`principal_parent_atoms`, first locant, writer atom order, folded H excluded).
- `openclatura` → `reshiki-chemistry`, `crate::naming` (status/reason enums, `IsotopeStyle`, `canonical_smiles`), and through `canonical_smiles` → `molecule_names::canonical_identity` in `crates/model` (admission checks, clearing unknown stereo, radical-electron preservation).
- `aromaticity` (via `sanitize`) calls `Capabilities::detect()` and draws on the process-wide reservation pool in `crates/process-heap`. `fusion` keeps thread-local caches.
- Tests that cross the boundary: catalog tests embed 17 OPSIN XML resources; phane tests use the InChI reader; `s2n_harness` verifies through the worker adapters. `crates/python` `drawing::from_name` is **production** code calling `opsin` directly.
- The ReShiki notice collector (`scripts/license_notices.py`) reads the repo-root `licenses/` tree (where the RDKit, centres and openclatura terms live today), `licenses/rust/manifest.json` supplements, and each separately licensed dependency's crate directory; workspace crates declared `MIT OR Apache-2.0` (such as `reshiki-chemistry`) are skipped and rely on the root tree.
- `canonical_smiles` in `src/naming.rs` holds s2n admission checks (SMILES size cap, CX/reaction rejection, atom-token bound, allene/coordination stereo rejection, kept tetrahedral centres). openclatura reaches it from `graph_io.rs`, `mod.rs`, `radical_ions.rs` and `unknown_stereo.rs`.
- openclatura's own tests call the n2s worker in-process through `crate::naming::{local, worker}` (27 call sites, e.g. `parent_pipeline.rs`, `unknown_stereo.rs`, `mod.rs`); `radical_ions.rs` tests use `crate::document`.
- Churn: 87 commits since 2026-09-01 touch naming paths; 29 of 77 open issues are naming bugs.

## 2. Goals and non-goals

Goals

1. ReShiki depends on the naming library like any other dependency; naming bugs, corpora and oracle gates live in the library repo.
2. A typed, versioned structure IR becomes the boundary in both directions (Phase B).
3. General-purpose: Rust, Python (PyPI), JS/WASM (npm), no ReShiki dependency.
4. Provenance explicit per file and per ported function; every distribution channel ships the complete third-party notices.
5. The move itself changes no behaviour: every existing gate is byte-identical before and after.

Non-goals for the move: naming-quality work, algorithm rewrites, switching chemistry kernel. Process supervision, browser allocator ownership, wire/heap limits and the PubChem lookup stay in ReShiki.

## 3. Two phases

**Phase A — move, same interface.** Extract the code into crates with today's interface: SMILES in for s2n, semantic CXSMILES + parent locants out for n2s, worker protocol v4 unchanged. Proven byte-identical, then lifted into the new repo. This alone delivers the development-speed goal.

**Phase B — introduce the IR.** Designed from a complete inventory of what the two engines carry today, added as adapters next to the string interface, gated separately, then made the primary API. Mixing the two phases would bundle five compatibility changes (extraction, IR, parser path, result types, protocol v5) into a move that must stay byte-identical.

## 4. Library layout and crate graph

```
nazuchem/
  crates/
    nazuchem-ir/    original; MIT OR Apache-2.0. Shared vocabulary: NameStatus, NameReason, IsotopeStyle,
                  warnings; (Phase B) the structure IR + JSON Schema + CXSMILES read/write
    nazuchem-kernel/  perception kernel + canonical identity (the pure part of canonical_identity);
                  RDKit BSD-3 + centres BSD-2 + original; published, documented as internal
    nazuchem-n2s/   name → structure; today's crates/opsin renamed; MIT (OPSIN)
    nazuchem-s2n/   structure → name; today's src/naming/openclatura; openclatura MIT + OPSIN MIT data + original
    nazuchem/       facade: Parser, Namer, verify(); features "n2s", "s2n"; re-exports nazuchem-ir
  bindings/python  PyO3 + maturin abi3 → PyPI `nazuchem`
  bindings/wasm    wasm-bindgen → npm `nazuchem` (or `@ameyanagi/nazuchem`)
  corpora/         s2n bake-off/baseline/rules TSVs, openclatura graph oracle
  PROVENANCE.toml  per-file origin ledger (§6)
```

Dependency direction: `ir` ← `kernel` ← {`n2s`, `s2n`} ← `nazuchem`. Shared types and identity sit **below** both engines, so the facade never has to be imported by an engine (no facade↔s2n cycle).

Two adjustments needed to keep the library host-neutral without changing behaviour:

- **Budget injection.** `nazuchem-kernel` must not detect machine capabilities or touch a process-wide reservation pool. Today `aromaticity::perceive()` calls `Capabilities::detect()` itself; in the kernel every perception entry point takes a `Budget` argument instead, and so do the s2n entry points (`Namer::new(budget)`; the n2s `Parser` likewise for `modern_cip`). ReShiki's worker passes `Capabilities::detect().resolve()` and its pool, exactly as today, so ReShiki output stays identical; the snapshot gate runs with that budget. Library users get `Budget::default()`, documented as unbounded; it differs from ReShiki only where ReShiki's budget would have run out. The host resolves its limits once per naming request (as process-heap's own `Budget` documents), where the copied code used to re-detect them at every perception call; the two differ only if free memory changes during a request that runs longer than one second.
- **Admission policy split.** `canonical_smiles` (from `src/naming.rs`) moves into `nazuchem-s2n` together with all its checks: size cap, CX/reaction rejection, atom-token bound, allene/coordination stereo rejection and the kept-tetrahedral-centre comparison. `canonical_identity` (from `crates/model`) splits into a pure identity function (to `nazuchem-kernel`) and ReShiki's document-level admission rules (stay in `crates/model`). The engine keeps its other rejections (e.g. N/ring lone-pair centres), so every input rejected today is still rejected.

Not called `opsin` on crates.io: that would read as the upstream project. Docs say "a Rust port of OPSIN by Daniel Lowe".

Library invariants (stated in the README): no I/O, no threads, no process-wide mutable state; bounded thread-local caches exist in fusion. Worst-case time on adversarial input is unbounded; hosts needing hard limits run it supervised, as ReShiki does.

### Chemistry kernel (decision D2)

- **A (recommended): copy the exact transitive subset of `reshiki-chemistry`** → byte-identical, every corpus stays green. Cost: ~27.5k lines exist twice; fixes are mirrored deliberately.
- B: COSMolKit (`cosmolkit-core/-smiles/-stereo`, 0.5.0-rc) → no duplication but canonical SMILES, Kekulé choice and ranking ties may differ, forcing a re-baseline mid-move. Evaluate later as its own project behind the same corpora.

## 5. The structure IR (Phase B)

nazuchem's typed graph and versioned JSON (`{"nazuchem_ir":1}`) are the
lossless public contract, with a published JSON Schema. CXSMILES is a convenience
export; Molfile V3000 may follow later. Standard formats were rejected as the
contract because they depend on external interpretation, reorder atom identities,
or cannot carry every name annotation.

Designer decisions for v1 (A–H, including final PR #35 revisions):

- **A — Stereo unknowns and convention.** `Tetrahedral { center, neighbors:
  [Ref; 4], winding: Clockwise | Anticlockwise | Unknown }` and
  `DoubleBond { bond, refs: [u32; 2], relation: Cis | Trans | Unknown }` have no
  `known` boolean. An absent entry means unspecified; `Unknown` explicitly
  records an unknown configuration. Looking from `neighbors[0]` toward the
  tetrahedral center, `neighbors[1]`, `[2]`, `[3]` run in the stated winding.
  Double-bond `refs[i]` is attached to `bonds[bond].atoms[i]`; cis/trans compares
  those two references. OPSIN parity conversion is checked against known R/S CIP
  labels. The pinned OPSIN parser **erases ξ/unknown name descriptors before
  `to_ir`**; retaining them is a separately tracked parser limitation, not an IR
  encoding restriction.
- **B — Atropisomers.** `Atropisomer { bond, refs: [Ref; 4], winding:
  Clockwise | Anticlockwise | Unknown }` is included in v1; its JSON values are
  `"clockwise"`, `"anticlockwise"`, and `"unknown"`. The first pair belongs to
  `bond.atoms[0]`, the second to `bond.atoms[1]`. Looking from the first endpoint
  to the second, the shortest turn from `refs[0]` to `refs[2]` has the stated
  winding. These are reference-relative structural windings, not
  cached CIP labels; switching one pair reverses it, while reversing the bond
  and exchanging both pairs preserves it. OPSIN axes, including
  the binaphthol M/P and Ra/Sa corpus, are retained.
- **C — Enhanced stereo.** A known center without a group is absolute.
  `StereoGroup` contains only `And(id)` and `Or(id)`, serialized flat as
  `{"kind":"and","id":7,"stereo":[0]}`; members index `Molecule.stereo`.
  Validation rejects any group member whose winding or relation is unknown.
- **D — JSON.** `Ref` is an untagged atom index or the string
  `"implicit_hydrogen"` / `"lone_pair"`. Optional fields default on input and
  empty/None annotations are omitted on output. Atom defaults `charge: 0`,
  `radical_electrons: 0`, and `hydrogens_stated: false` are also omitted on write
  and restored on read. Documents and nested records
  are objects, never positional arrays; each variant accepts only its own
  fields. Readers inspect `nazuchem_ir` before the document shape and report a
  typed `UnsupportedVersion`; unknown fields and variants also have typed
  errors. A fully populated serialization test covers every wire type/variant
  and compares decoder field/variant lists with serialized keys to prevent drift.
  **IR v1 freezes at the first published release; until then it may change
  without a version bump. After freezing, every new field or variant increments
  `nazuchem_ir`; readers accept
  supported versions up to their own; writers emit the lowest version able to
  represent the document.** Only v1 exists today; zero and later versions are
  unsupported. Graph constraints are validated in addition to schema shape.
- **E — Typed identities.** Atoms carry `atom_map: Option<u32>` and
  `attachment_point: Option<u32>` (OPSIN `_AP{n}`); ordinary maps are accepted.
  `label: Option<String>` is genuine text/R-group content. Polymer star/end roles
  are `RepeatUnit.attachment_points`, not atom label strings. Locants retain all
  name-order values, and `Molecule.parent` records the complete principal parent.
  A hydrogen carrying locants, parent membership, isotope, charge, label/map, or
  any stereo role stays a vertex; only unannotated ordinary H can become a count.
- **F — Scope.** No coordinates, zero-order or dative bonds in v1. Layout belongs
  to hosts; coordination chemistry is outside both engines' supported contract.
  Allene axial chirality is not representable in v1. Localized Kekulé orders and
  integer annotations preserve structural `Eq`; it does not test stereo-equivalence
  under reordered atoms or references. The step 8b corpus proof below confirms
  that v1 remains Kekulé-only; no aromatic atom or bond field is needed.
- **G — CXSMILES.** JSON is lossless. CX output uses standard fields plus
  `$_AV:` first-locant atom values and just `nazuchemParent` and
  `nazuchemLocants` properties (the latter stores every locant). No topology or
  index re-embedding is permitted. Export never asserts configuration absent
  from the IR: unrepresentable stereo, including the outer E bonds in
  `(2E,6E)-octa-2,4,6-triene` with its unspecified middle bond or
  `(2E,8E)-deca-2,4,5,6,8-pentaene` with its unspecified odd cumulene, returns
  `CxsmilesError::StereoNotExpressible`. The legacy comparison projection is
  hidden test support. Every emitted corpus structure is reparsed with the
  kernel and checked for specified/unspecified state and configuration, including
  maximal odd cumulene chains. Use simple pruning and a final audit, plus greedy
  distinct wavy-bond carriers. Remove 2-SAT and augmenting matching: comparing
  both implementations gives zero CX refusals in all three IR corpora, with
  unchanged exact/canonical counts. Ordinary-chain properties remain 940
  exports / 424 refusals; cumulene properties remain 336 / 48. Across 64 synthetic
  unknown-centre carrier cases, matching gives 48 / 16 and greedy assignment
  gives 42 / 22. These six additional conservative refusals retain their complete
  configuration in JSON.
- **H — Facade.** Re-export `AxialWinding`, `IR_VERSION`, and `JSON_SCHEMA`.
  Use `CxsmilesError` throughout the public API and docs, without a spelling alias.

**Step 8b aromaticity proof (2026-11-09).** The regular s2n corpus gate
compares all 1,650 distinct SMILES in the bakeoff (703), openclatura (894),
rules (94), and golden inputs (the golden set adds no new strings). A test-only
SMILES → kernel → IR bridge retains prepared atom/bond order and discards all
aromatic flags, leaving only localized Kekulé orders. Direct IR preparation
reproduces the aromatic atom/bond snapshot for **all 1,650** inputs before
naming admission. The 1,637 connected, admitted graph snapshots additionally
have identical canonically localized orders; the other 13 graphs are rejected
on both paths. These snapshots include 987 aromatic graphs, 11,212 aromatic
atoms, 11,760 aromatic bonds, and 149 localized nonaromatic bonds between
aromatic endpoints. Thus the fusion audit retains its localized-pi rejections.
Full `Name` comparison (text, status, ordered alternatives, omitted stereo and
units, and verification) gives **1,616 identical / 0 differing / 34 rejected on
both / 0 rejected on one**. All source strings convert to IR. **The proof passes;
IR v1 stays Kekulé-only, with no schema change.**

The complementary n2s proof runs every `Structure::to_ir` result from the
structure oracle, configuration oracle and generated-name golden corpus
through both `name_molecule` and legacy semantic CXSMILES → `name_smiles`.
It covers 4,630 structures from 4,914 parse records; the remaining 284 records
are parser failures and produce no IR. There are no IR conversion failures or
missing semantic-CXSMILES baselines.

| Corpus | Complete names identical | Differing | Rejected on both | Rejected on one |
| --- | ---: | ---: | ---: | ---: |
| Distinct SMILES inputs | 1,616 | 0 | 34 | 0 |
| n2s structure oracle (2,605 structures) | 2,293 | 0 | 312 | 0 |
| n2s configuration oracle (416 structures) | 336 | 0 | 80 | 0 |
| n2s generated-name golden (1,609 structures) | 1,608 | 0 | 1 | 0 |
| **Total (6,280 comparisons)** | **5,853** | **0** | **427** | **0** |

All 34 source-SMILES rejection messages also match byte for byte. The regular
test pins these counts and compares every `Name` field; the largest n2s oracle
uses four bounded test-only shards without skipping or caching any record.
Focused tests exercise R/S, E/Z and M/P through the naming API, reference
permutations, reversed axes, all three Unknown forms, enhanced groups, isotope
styles, hydrogen vertices and sparse original atom IDs. Admission regressions
cover Unknown with symmetric ligands/cumulene bonds, nonbiaryl unknown axes,
and disconnected unknown salts (including rejection precedence).

The direct adapter preserves atom/bond indices, hydrogen origin/counts, charge,
isotopes, radicals, maps and reference-relative stereo. The ordinary naming
projection folds nonisotopic H vertices with the same kernel removal/winding
rules as SMILES import, retaining original (possibly sparse) surviving IDs;
the caller's IR, including stereo-reference and annotated H vertices, is
unchanged. This fixes the initial n2s comparison differences caused by the
legacy engine encountering those H vertices as substituents. Isotopic H remains
a vertex. Enhanced groups and unsupported CX-only annotations retain the
existing rejection policy, and unknown axes receive the same biaryl eligibility
check as specified axes.

Atoms additionally retain element, isotope, charge, hydrogen count and whether
that count was stated, and radical electrons. Atom index is identity. Bonds
carry endpoint identities and localized orders. Position variation retains the
attachment bond and all possible endpoints. Repeat units carry atom membership,
connectivity (`HeadToTail`, `HeadToHead`, `Unknown`), optional text, attachment
roles and `crossing_bonds` equal to the **full boundary** of the member atom set.
The displayed locant is the first locant of a principal-parent atom.

Result types must keep today's semantics:

- `GeneratedName`: `text, status, alternatives, omitted_stereo, omitted_units` (unknown units ≠ omitted locants), `atom_locants` with a flag for when they do not apply (radical-ion names suppress the recovered parent's locants), and `verification: { target: Input | RecoveredNeutral { name, smiles }, outcome }`.
- Two parse modes stay separate: strict import parsing and generated-name verification parsing (radicals allowed), as worker v4's `Parse` and `Verify` do today.
- Before switching the worker to IR (protocol v5), measure IR JSON size against the 32 KiB request / 64 KiB response caps and the browser worker's memory limits; raise caps only with evidence.

## 6. Provenance, licensing and IUPAC citations

1. **Rust source files** start with an origin header:
   ```rust
   // SPDX-License-Identifier: MIT
   // Origin: port of OPSIN FusedRingNumberer.java (github.com/dan2097/opsin @ 42e87ac6), © 2017 Daniel Lowe
   // Changes: FR-5.4(b) per-element heteroatom seniority (Nazuchem; diverges from upstream) — see PORTING.md
   ```
   or `// SPDX-License-Identifier: MIT OR Apache-2.0` + `// Origin: original Nazuchem code`. Mixed files mark function-level boundaries: `// Ported from openclatura chains.py::_longest_chain` vs `// Nazuchem extension (not in openclatura 0.4.3): …`.
2. **Data, byte-preserved upstream XML and binary DFA tables** get no inline header (that would break byte equality); they are listed in `PROVENANCE-nazuchem.toml` with origin, upstream path, commit, licence and SHA-256.
3. `PROVENANCE-nazuchem.toml` covers every file recursively (source, data, fixtures, scripts). `uv run --locked python scripts/check_provenance.py` fails CI on a missing entry, a header/ledger mismatch, a missing mixed-function marker or a non-Rust hash change. The release notice collector gathers the crate-local terms.
4. **Notices travel inside each publishable crate**: every crate directory contains its own `LICENSE*`, `NOTICE` and a `licenses/` folder with the complete upstream texts that apply to it (n2s: OPSIN; s2n: openclatura + OPSIN + CPython; kernel: RDKit + centres + LLVM/MSVC; facade: all). Wheels and npm packages include the same. Release CI asserts the **contents**, not just filenames, of the `.crate`, wheel and npm tarball, and ReShiki's CI asserts its aggregated `THIRD-PARTY-NOTICES.txt` still contains OPSIN, openclatura, RDKit and centres after the switch. During Phase A (while the new crates still live in the ReShiki workspace), those crates are not `MIT OR Apache-2.0` as a whole: declare their real licence expressions, copy the matching root `licenses/` entries (rdkit, centres, stdlib, openclatura, cpython) into each crate's own `licenses/` folder, retain the root copies while ReShiki contains those implementations, and make sure the collector's "original workspace crate" shortcut does not skip them.
5. **Wording**: README, crate docs and package metadata say "port of OPSIN" / "began as a port of openclatura 0.4.3, extended" — never "original" for those parts.
6. **IUPAC references**: rule functions carry a doc tag such as `/// IUPAC 2013: P-44.1.2.2`; fused `FR-`, spiro `SP-`, von Baeyer `VB-`, ring assemblies `RB-`. Rule numbers and short quotes only (Blue Book text is IUPAC/RSC copyright). A script turns the tags into a rule-coverage page. Note that OPSIN-compatibility data (e.g. `fusion_entry_geometry.json`) already says it is "parser compatibility, not IUPAC PIN certification"; keep that distinction.
7. **Licences**: original code MIT OR Apache-2.0; OPSIN and openclatura MIT; RDKit BSD-3; centres BSD-2; CPython Python-2.0 (complete licence/history); LLVM libc++ and Microsoft STL Apache-2.0 WITH LLVM-exception; ChEBI CC BY 4.0 (only if the catalog moves).

**Provenance audit** (step 2): all 103 openclatura Rust files and their data against openclatura 0.4.3 (`90dc5bcb`) plus ReShiki git history, separating openclatura-ported, OPSIN-derived and original code, function by function; the chemistry subset against RDKit/centres. The output _is_ the first `PROVENANCE-nazuchem.toml`.

## 7. Migration steps

Baseline first: freeze the outputs of the OPSIN gates (8,100 lexical, 2,624 structure, 640 configuration), all s2n corpora, the openclatura graph oracle, and ReShiki's naming/app/Python tests. Every Phase A step must reproduce them byte-for-byte on macOS, Arch Linux x64, Ubuntu arm64 and Windows x64 (local test hosts first, CI only as the final check).

**Phase A (one ReShiki PR per step, each independently reviewed)**

| Step | Work                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| ---- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 0    | Land or close in-flight naming PRs                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| 1    | Inventory the chemistry transitive closure; create `crates/nazuchem-kernel` (exact copy + `Budget` injection); point `opsin` and openclatura at it                                                                                                                                                                                                                                                                                                                                                                                    |
| 2    | Provenance audit → `PROVENANCE-nazuchem.toml` (upstream version recorded per file; e.g. `fusion_entry_geometry.json` comes from OPSIN 2.9.0, not `42e87ac6`), headers, per-crate `licenses/`, collector fix                                                                                                                                                                                                                                                                                                                           |
| 3    | `crates/nazuchem-ir` with the shared vocabulary only; split `canonical_identity` and move it, `apply_axis_markers` and `stereo_lint::is_pseudoasymmetric` into `nazuchem-kernel`, removing the step 1 JSON bridge (`src/naming/kernel_bridge.rs`)                                                                                                                                                                                                                                                                                     |
| 4    | Move openclatura → `crates/nazuchem-s2n` (SMILES in, unchanged `Name` out) together with `canonical_smiles`, `verify_identity` and `reverse::verify_alternative_identity` (used by openclatura tests), switched to the kernel and the injected `Budget`; move its tests and fixtures; its tests that call n2s through `crate::naming::{local, worker}` (27 sites) use a dev-dependency on the n2s crate instead; replace the worker-adapter harness, `crate::document` and InChI-reader phane test helpers with library-local support |
| 5    | Rename `crates/opsin` → `crates/nazuchem-n2s` (update `nazuchem-s2n`'s dev-dependency in the same PR); facade crate `nazuchem` with today's string API; switch ReShiki callers incl. `crates/python` `from_name`; keep catalog XML evidence pinned                                                                                                                                                                                                                                                                                    |
| 6    | `git filter-repo` the crates, corpora and licences into the new repo (history kept); ReShiki depends by git rev and deletes its copies. Real-app check: name import + locant hover, generated name panel, browser naming worker                                                                                                                                                                                                                                                                                                       |
| 7    | New repo CI: fmt, clippy, tests, oracle gates, provenance, 3 OS + wasm32. Transfer naming issues (`gh issue transfer`)                                                                                                                                                                                                                                                                                                                                                                                                                |

**Phase B (in the new repo)**

| Step | Work                                                                                                              |
| ---- | ----------------------------------------------------------------------------------------------------------------- |
| 8a   | IR v1 to §5's requirements; n2s → IR adapter, convenience CXSMILES export, lossless JSON/schema and corpus gates next to the string API |
| 8b   | IR → structure-to-name adapter and the §5 aromaticity/localized-pi proof on the full corpora |
| 9    | Python and npm bindings on the IR API                                                                             |
| 10   | ReShiki: `IR → chemistry State` converter, worker protocol v5 after the size measurement                          |
| 11   | 0.1.0 via tag-triggered trusted publishing (crates.io, PyPI, npm); ReShiki moves from git rev to version          |

Before 0.1.0, decide whether to add `#[non_exhaustive]` to `NameStatus` and `NameReason` so the shared vocabulary's public API evolution policy is explicit.

Phase B step 8a decision: both `NameStatus` and `NameReason` are
`#[non_exhaustive]`, as are the IR's public enums. New naming reasons and stereo
variants may be added without requiring exhaustive downstream matches; unknown
JSON fields and variants still fail explicitly. The schema is generated with
the dev-only schemars dependency and checked against the published file. Corpus
tests validate each emitted document against that file with a strict test-only
validator (unknown schema keywords fail), avoiding a runtime schema dependency.
Step 8a adds the n2s adapter and CXSMILES exporter. Step 8b adds
`nazuchem_s2n::name_molecule(&Molecule, IsotopeStyle, &Budget)` and facade
`nazuchem::molecule_to_name(&Molecule)`, plus the regular corpus and focused
stereo/admission gates described in §5.

Step 8a comparison policy: exact legacy semantic CXSMILES is checked first;
the test-only projection may fold metadata-free stereo-reference H into a
virtual H for that legacy string comparison. Full IR and public CX output keep
and independently check those vertices. Otherwise the kernel canonical
isomeric identity must agree. The fallback fails
closed if preparation drops any specified stereo unit. Independently, source
OPSIN tetrahedral winding, double-bond relation and atropisomer configuration are
compared directly against the IR under the source identity map, alongside full
locants and parent membership before H folding. Every successful CX export is
reparsed without stereo cleanup to detect lost or newly asserted stereo; a typed
`StereoNotExpressible` is counted explicitly while its JSON/source graph checks
still run. Standard labels/groups/position variation/SRU layers are compared
separately under the writer identity map. The original golden gate and fixtures
remain unchanged. See `crates/nazuchem-ir/docs/cxsmiles.md` for export limits.

Risks: step 4 is the hardest (openclatura uses `crate::naming` in 112 places plus `super::super` paths; only about 32 are the status/reason enums, the rest are `canonical_smiles`, `IsotopeStyle` and test calls into the worker). The two chemistry copies can drift; mirror fixes deliberately. The npm package carries ~5.7 MB of OPSIN resources; features allow n2s-only or s2n-only builds.

## 8. Name

**nazuchem** (decided 2026-11-09). Free on crates.io, PyPI, npm (`nazuchem`, `@ameyanagi/nazuchem`) and `github.com/Ameyanagi/nazuchem` when checked on 2026-11-09. Crates: `nazuchem` (facade), `nazuchem-ir`, `nazuchem-kernel` (chemistry; avoids `nazuchem-chem`), `nazuchem-n2s`, `nazuchem-s2n`.

## 9. Decisions (made 2026-11-09)

- D1 Name: nazuchem.
- D2 Chemistry kernel: copy the reshiki-chemistry subset now; evaluate COSMolKit later as its own project.
- D3 Repository: github.com/Ameyanagi/nazuchem.
- D4 Licence for original code: MIT OR Apache-2.0.
- D5 Naming work during Phase A: freeze naming features and non-urgent fixes until step 6.
- D6 IR: nazuchem's own typed graph with versioned JSON (§5).
- Optional courtesy: tell OPSIN's author and lamalab-org before the first public release.
