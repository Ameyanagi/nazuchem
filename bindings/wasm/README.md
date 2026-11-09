<!-- SPDX-License-Identifier: MIT OR Apache-2.0 -->
<!-- Origin: original -->
<!-- Source: Nazuchem | bindings/wasm/README.md | Phase B step 9 -->

# nazuchem npm/WebAssembly bindings

The package enables both Rust naming engines and exposes four synchronous operations:

```sh
bun add nazuchem
```

```ts
import { nameToSmiles, smilesToName, nameToMolecule, moleculeToName } from "nazuchem";

nameToSmiles("ethanol"); // "C(C)O"
smilesToName("CCO").text; // "ethanol"
const molecule = nameToMolecule("ethanol"); // plain IR v1 object
moleculeToName(molecule).text; // "ethanol"
moleculeToName(JSON.stringify(molecule)).text; // "ethanol"
```

`smilesToName` and `moleculeToName` return the complete `Name`: `text`,
`status`, `reason`, `alternatives`, `omittedStereo`, `omittedUnits`, and
`verification`. The preference assessment is exposed as `status`
(`preferred`, `not_preferred`, or `unverified`) and nullable `reason`; each
alternative has its own `name`, `status`, and `reason`. `verification` is either
`null` or the recovered neutral `{ name, smiles }` used for a radical-ion name.
These assessments preserve the Rust engine's values. A successful name does
not certify preferred IUPAC name compliance.

IR objects retain the Rust JSON field names and optional/default semantics.
`IR_VERSION` is the numeric version and `JSON_SCHEMA` is the published JSON
Schema string; the schema is also exported as `nazuchem/schema.json`.
TypeScript declarations cover every IR field and stereo variant. Integer ranges,
uniqueness, references, graph validity and naming admission are enforced at
runtime by Rust. Object input rejects nonfinite numbers before JSON encoding.
The IR is not mutated during naming.

Errors are instances of `NazuchemError`, with stable lowercase `kind` values
`init`, `parse`, `warning`, `serialize`, `naming`, and `internal`, and the engine's
message. Malformed JSON or rejected IR uses `serialize`; an admitted IR that the
naming engine cannot name uses `naming`. WebAssembly runtime traps are mapped to
`internal`. Caught errors do not disable the instance: subsequent calls remain
usable, including naming ethanol after a rejected extreme spiro descriptor.

The root entry supports bundlers and selects CommonJS on Node. Explicit Node
imports use `nazuchem/node`. Bundlers must support asynchronous WebAssembly
imports, for example Vite with `vite-plugin-wasm` or webpack with
`experiments.asyncWebAssembly`. The web path loads an approximately **18.6 MB**
Wasm resource before transfer compression. The browser-native module requires
initialization:

```ts
import init, { nameToMolecule } from "nazuchem/web";
await init(); // loads the shared ../nazuchem_bg.wasm resource
console.log(nameToMolecule("ethanol"));
```

The build emits `dist/bundler`, `dist/web`, and `dist/node`, referencing one
`dist/nazuchem_bg.wasm`. The bundler's generated import module,
`dist/nazuchem_bg.js`, lives beside it so imports inside the Wasm resolve
correctly. Serve the complete `dist` tree so the web loader can resolve its parent
asset. CommonJS consumers of `nazuchem/web` receive the same explicit loader,
including async initialization and `initSync`. Generated payloads must be byte-identical before the builder shares them.

Use Bun for JavaScript tools and uv for Python tooling. From this directory:

```sh
bun install --frozen-lockfile
rustup target add wasm32-unknown-unknown
# Install matching official wasm-bindgen 0.2.128 and Binaryen 133 release binaries.
# Set WASM_BINDGEN and WASM_OPT if they are outside PATH.
bun run build
bun run test
bun run check:package
```

The build verifies the wasm-bindgen CLI against its exact pin and `Cargo.lock`,
and requires Binaryen 133. `CARGO_TARGET_DIR` selects the build directory;
otherwise it uses the repository's `target` directory. The maintained V4 build
selects `wasm-release`: fat LTO, one codegen unit, Rust optimization level 3, and
panic abort, followed by `wasm-opt -O3`. Native and Python release settings are
independent. `bun scripts/build.ts --skip-cargo` requires an existing artifact
built with that profile. Every build reports raw, gzip-9, and Brotli-11 sizes.
The historical size/performance experiment, method, and results are documented in
[docs/wasm-size.md](../../docs/wasm-size.md); its full harness is preserved at the
`wasm-size-experiment` git tag.

The tests include exact schema/type equivalence, the exported Rust naming-reason
vocabulary, strict bare-package TypeScript consumers under Node16 and bundler
resolution, Node API tests, web initialization, and static bundler instantiation
in a fresh Node process. The package check runs the same runtime suite against an
extracted real tarball and checks its four notices against an independent,
nonempty wasm32 dependency-notice collection.

Builds collect complete normal/build Cargo dependency licence texts for wasm32
through `uv run --no-project python` and the shared notice collector. Packages
include `LICENSE-MIT`, `LICENSE-APACHE`, `NOTICE`, and `THIRD-PARTY-NOTICES.txt`.
Those generated files retain third-party scope and attribution. Outstanding
provenance decisions described in the root `NOTICE` still apply.

IR v1 is not yet frozen. Supported chemistry, stereochemistry limits, unbounded
default budgets, and externally supervised time limits follow the
[Rust facade documentation](../../crates/nazuchem/README.md).
