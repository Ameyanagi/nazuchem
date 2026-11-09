# nazuchem

Rust library for IUPAC-style chemical name ↔ structure conversion. It accepts
names and molecular SMILES through string APIs and exposes a versioned structure
IR for name-to-structure results, with Python and JavaScript/WebAssembly bindings.

Extracted from [ReShiki](https://github.com/Ameyanagi/ReShiki), a molecular
drawing app. The name-to-structure engine is a Rust port of OPSIN by Daniel Lowe
(MIT). The structure-to-name engine began as a port of openclatura 0.4.3 by
lamalab-org (MIT), extended in ReShiki. The chemistry kernel adapts RDKit and
SiMolecule centres and uses Rust standard library stable sorting; the naming
engine also includes a CPython integer-set adaptation.

The 0.x API is not yet stable.

```rust
fn main() -> Result<(), nazuchem::Error> {
    // Equivalent to CCO; atom order follows the parser.
    assert_eq!(nazuchem::name_to_smiles("ethanol")?, "C(C)O");
    let name = nazuchem::smiles_to_name("CCO")?;
    assert_eq!(name.text, "ethanol");
    Ok(())
}
```

The facade enables `n2s` and `s2n` by default; either feature can be selected
independently with default features disabled. `name_to_smiles` applies strict
OPSIN options, rejects warnings and returns semantic CXSMILES.
`smiles_to_name` returns the complete `Name`, including its preference
assessment, alternatives and omitted stereochemistry.

| Crate | Purpose |
| --- | --- |
| [`nazuchem`](crates/nazuchem) | Public string-API facade and engine re-exports |
| [`nazuchem-ir`](crates/nazuchem-ir) | Shared naming vocabulary and versioned structure IR |
| [`nazuchem-kernel`](crates/nazuchem-kernel) | Internal molecular perception and canonical identity |
| [`nazuchem-n2s`](crates/nazuchem-n2s) | OPSIN name-to-structure port |
| [`nazuchem-s2n`](crates/nazuchem-s2n) | Extended openclatura structure-to-name port |
| [`nazuchem-python`](bindings/python) | Python 3.10+ abi3 extension and typed Python API |
| [`nazuchem-wasm`](bindings/wasm) | JavaScript package with bundler, web and node builds |

## Python

Install the `nazuchem` distribution with uv:

```sh
uv add nazuchem
```

For source builds and local tests with uv and maturin, see
[the binding instructions](bindings/python/README.md).

```python
import nazuchem

smiles = nazuchem.name_to_smiles("ethanol")
name = nazuchem.smiles_to_name(smiles)
assert name.text == "ethanol"

molecule = nazuchem.name_to_molecule("(2R)-butan-2-ol")
assert nazuchem.molecule_to_name(molecule).text == "(2R)-butan-2-ol"
```

`Name` is a frozen dataclass with `text`, `status`, `reason`, `alternatives`,
`omitted_stereo`, `omitted_units`, and `verification`. Its collections are tuples;
alternatives and verification are frozen dataclasses too. `name_to_molecule`
returns a plain IR dictionary, and `molecule_to_name` accepts that dictionary or
its JSON string. `IR_VERSION` is an integer, `JSON_SCHEMA` is a JSON string, and
`NazuchemError` subclasses `ValueError`; its `kind` identifies the error stage
(`init`, `parse`, `warning`, `serialize`, `naming`, or `internal`). The package
carries stubs and `py.typed`.

## JavaScript

Install the `nazuchem` npm package with bun:

```sh
bun add nazuchem
```

For source builds, see [the package instructions](bindings/wasm/README.md).
The root entry selects node in Node.js and the bundler build otherwise;
`nazuchem/node` and `nazuchem/web` select a target explicitly. The bundler build
requires WebAssembly module support, such as
[`vite-plugin-wasm`](https://github.com/Menci/vite-plugin-wasm) for Vite or
[`experiments.asyncWebAssembly`](https://webpack.js.org/configuration/experiments/#experiments)
for webpack.

```typescript
import { nameToSmiles, smilesToName, nameToMolecule, moleculeToName } from "nazuchem";

console.log(smilesToName(nameToSmiles("ethanol")).text); // ethanol
const molecule = nameToMolecule("(2R)-butan-2-ol");
console.log(moleculeToName(molecule).text); // (2R)-butan-2-ol
```

For direct browser use, initialize the web build before calling the functions:

```typescript
import init, { nameToMolecule } from "nazuchem/web";

await init();
const molecule = nameToMolecule("ethanol");
```

The web entry loads the shared `dist/nazuchem_bg.wasm` asset, currently about
18.6 MB uncompressed. See the [Wasm size measurements](docs/wasm-size.md) for
build comparisons.

JavaScript returns plain `Name` and IR objects, accepts an IR object or JSON
string in `moleculeToName`, and throws `NazuchemError` with the same `kind` values
as Python. `Name` fields use camelCase (`omittedStereo` and `omittedUnits`);
IR objects preserve the schema's snake_case field names. Naming results share
Python's status and reason vocabulary. The package includes `IR_VERSION`,
`JSON_SCHEMA`, hand-written TypeScript declarations, and the schema asset at
`nazuchem/schema.json`. CI checks the IR declarations against the Rust JSON
schema.

Both packages include `LICENSE-MIT`, `LICENSE-APACHE`, `NOTICE`, and generated
`THIRD-PARTY-NOTICES.txt`. The shared collector follows normal and build Cargo
dependencies for the selected Rust target and includes complete crate-local
licences and notices, including embedded engine/kernel terms. Wheel checks
compare the packaged bytes with freshly generated files; the npm package check
tests a real bun-produced archive and compares its notice bytes with freshly
generated files.

## Scope and limits

The pinned OPSIN revision and supported s2n rules do not cover every chemical
name or structure. A successful generated name is not necessarily a preferred
IUPAC name. Unsupported chemistry returns an error; OPSIN compatibility data
does not certify IUPAC PIN compliance. See the [facade documentation](crates/nazuchem/README.md)
for admission limits and unsupported stereochemistry.

`name_to_molecule` returns IR v1 with stable identities, full locants, parent
membership, and tetrahedral, double-bond and atropisomer stereo. JSON is the
lossless interchange format. The facade also exports `AxialWinding`,
`CxsmilesError`, `IR_VERSION`, and `JSON_SCHEMA` without requiring engine
features. CXSMILES is a convenience export and returns
`CxsmilesError::StereoNotExpressible` when standard notation cannot preserve
specified/unspecified stereo, including some partially specified polyenes and
atropisomers. Allene axial chirality is not representable in v1.
Layout coordinates belong to hosts; v1 has no dative or
zero-order bonds and coordination chemistry is outside both engines' contract.

IR v1 freezes at the first published release; until then it may change without
a version bump. After freezing, new fields or variants require a new IR version.
`Eq` is structural equality of the ordered graph and annotations, not
stereochemical equivalence. Atom JSON omits zero charge, zero radical electrons,
and false `hydrogens_stated`, restoring them on read; for example, an
unannotated methane atom is `{"element":6,"hydrogens":4}`. Axial JSON windings
are `"clockwise"`, `"anticlockwise"`, or `"unknown"` relative to the ordered
references, as defined in the [facade documentation](crates/nazuchem/README.md).

The IR distinguishes explicit `Unknown` stereo from an absent assertion, but
the pinned OPSIN parser currently erases ξ/unknown descriptors in names before
`Structure::to_ir`. Preserving those source annotations is a separate parser
limitation tracked in the [IR plan](docs/plans/naming-split.md#5-the-structure-ir-phase-b).

The library performs no I/O, creates no threads and has no process-wide mutable
naming state. Embedded parser resources are immutable; fusion uses bounded
thread-local caches. Default budgets are unbounded. Hosts can inject a `Budget`
into the engine APIs; its analysis timeout covers kernel perception, not an
entire naming request. Worst-case time on adversarial input is unbounded; hosts
requiring hard limits must supervise requests externally.

## Development

For registry setup and the release workflow, see
[Releasing nazuchem](docs/releasing.md).

Use the pinned Rust toolchain and checked lockfile:

```sh
export CARGO_TARGET_DIR="$PWD/target"
# On hosts whose default Python is older than 3.10:
export PYO3_PYTHON="$(uv python find 3.12)"
cargo fmt --all --check
cargo clippy --workspace --all-targets --locked --offline -- -D warnings
cargo test --workspace --locked --offline
uv run --no-project python scripts/check_provenance.py
uv run --no-project python -m unittest tests.test_check_provenance tests.test_binding_notices
```

The golden gate reads `bakeoff_corpus.tsv`, `openclatura_tests.tsv` and
`rules_cases.tsv` from `tests/fixtures/s2n/`, in that order, keeping each SMILES
at its first occurrence. It captures the complete s2n `Name` or error and
strict n2s results for the main name, alternatives and any recovered neutral
verification name, deduplicated in first-occurrence order within each row, in
[`tests/fixtures/golden/s2n_n2s.jsonl`](tests/fixtures/golden/s2n_n2s.jsonl).
Review intentional behavior changes before regenerating it; do not edit rows
by hand:

```sh
NAZUCHEM_UPDATE_GOLDEN=1 cargo test -p nazuchem --test golden --locked --offline
```

Frozen OPSIN fixtures run without Java or a live oracle. The
[extraction plan](docs/plans/naming-split.md) records the migration decisions
and later IR/bindings work.

## Licensing and provenance

Original code uses [MIT](LICENSE-MIT) OR [Apache-2.0](LICENSE-APACHE), except
original additions inside `nazuchem-n2s`, which are MIT only like the OPSIN
port they extend.
Incorporated material keeps its own terms: OPSIN/openclatura MIT, RDKit
BSD-3-Clause, centres BSD-2-Clause, and CPython Python-2.0. See [NOTICE](NOTICE),
the per-file [`PROVENANCE-nazuchem.toml`](PROVENANCE-nazuchem.toml) ledger and
each crate's `NOTICE` and `licenses/` texts for the applicable scope and pins.

Sorts with equal keys keep their input order (Rust's stable sort); nazuchem
does not reproduce C++ standard-library tie order.

`crates/nazuchem-s2n/src/data/bluebook/p25_retained.tsv` is retained as factual
data transcribed from IUPAC 2013 Recommendations P-25.2.1 Table 2.8.
It is excluded from the MIT/Apache-2.0 grant; its ledger licence remains
`NOASSERTION` (owner decision 2026-11-10).

`tests/fixtures/s2n/bakeoff_corpus.tsv` is a repository-only test fixture whose
CIDs and IUPAC names come from PubChem (NCBI). It is excluded from all published
crates, wheels and npm packages (owner decision 2026-11-10).
