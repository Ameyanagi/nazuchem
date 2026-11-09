<!-- SPDX-License-Identifier: MIT OR Apache-2.0 -->
<!-- Origin: original -->
<!-- Source: Nazuchem | bindings/wasm/README.md | Phase B step 9 -->

# nazuchem npm/WebAssembly bindings

Private, unpublished bindings for the Rust facade's strict chemical naming APIs.
The package enables both engines and exposes four synchronous operations:

```ts
import { nameToSmiles, smilesToName, nameToMolecule, moleculeToName } from "nazuchem";

nameToSmiles("ethanol"); // "C(C)O"
smilesToName("CCO").text; // "ethanol"
const molecule = nameToMolecule("ethanol"); // plain IR v1 object
moleculeToName(molecule).text; // "ethanol"
moleculeToName(JSON.stringify(molecule)).text; // "ethanol"
```

`smilesToName` and `moleculeToName` return the complete `Name`: `text`,
`status`, `reason`, `alternatives`, `omitted_stereo`, `omitted_units`, and
`verification`. Rust's tagged preference assessment is exposed as `status`
(`preferred`, `not_preferred`, or `unverified`) and nullable `reason`; each
alternative has its own `name`, `status`, and `reason`. `verification` is either
`null` or the recovered neutral `{ name, smiles }` used for a radical-ion name.
These assessments preserve the Rust engine's values. A successful name does
not certify preferred IUPAC name compliance.

IR objects use the exact Rust JSON field names and optional/default semantics.
`IR_VERSION` is the numeric version and `JSON_SCHEMA` is the published JSON
Schema string. TypeScript declarations cover every IR field and stereo variant;
the checks derive types from the complete schema and require exact type equality.
Integer ranges, uniqueness, references, graph validity and naming admission are
enforced at runtime by Rust. The IR is not mutated during naming.

Errors are instances of `NazuchemError`, with a `kind` retaining the facade
variant (`Init`, `Parse`, `Warning`, `Serialize`, `Naming`) and the engine's
message. Malformed JSON or rejected IR input uses `Serialize`; an admitted IR
that the naming engine cannot name uses `Naming`.

The root export supports bundlers and selects CommonJS on Node. Explicit Node
imports use `nazuchem/node`. The browser-native module requires initialization:

```ts
import init, { nameToMolecule } from "nazuchem/web";
await init(); // loads the adjacent .wasm resource
console.log(nameToMolecule("ethanol"));
```

The package's release build emits `dist/bundler`, `dist/web`, and `dist/node`.
Use Bun for JavaScript tools and uv for Python tooling. From this directory:

```sh
bun install --frozen-lockfile
rustup target add wasm32-unknown-unknown
cargo install wasm-bindgen-cli --version 0.2.128 --locked
bun run build
bun run test
bun run check:package
```

The build checks that the installed CLI, its exact pin, and `Cargo.lock` all use
wasm-bindgen 0.2.128. `WASM_BINDGEN` can point to a matching binary. It uses
`CARGO_TARGET_DIR` when set and otherwise the repository's `target` directory.
`bun scripts/build.ts --skip-cargo` packages an already compiled release.
Every build reports raw and gzip (level 9, Bun's `node:zlib`) bytes for the default
n2s+s2n Wasm; no `wasm-opt` pass is applied. The test command includes strict TypeScript
checking and schema equivalence before running Bun tests against the Node build.
The package check creates a real local tarball and verifies its notices.

Measured with the pinned Rust toolchain and wasm-bindgen 0.2.128, the default
`n2s` + `s2n` release build is **24,634,099 raw bytes** and **5,720,815 gzip
bytes** (Bun 1.4.4 `node:zlib`, level 9, no `wasm-opt`). Each target uses
that same Wasm payload. These are Wasm-only transfer sizes; JavaScript wrappers,
TypeScript declarations and licence notices add package overhead.

Builds collect the complete normal/build Cargo dependency licence texts for
the wasm32 target through `uv run --no-project python` and the shared notice
collector. npm packages include `LICENSE-MIT`, `LICENSE-APACHE`, `NOTICE`, and
`THIRD-PARTY-NOTICES.txt`. Those generated files retain third-party scope and
attribution. Outstanding provenance decisions described in the root `NOTICE`
still apply; `private: true` and Cargo `publish = false` remain in force.

IR v1 is not yet frozen. Supported chemistry, stereochemistry limits, unbounded
default budgets, and externally supervised time limits follow the
[Rust facade documentation](../../crates/nazuchem/README.md).
