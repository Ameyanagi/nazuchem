<!-- SPDX-License-Identifier: MIT OR Apache-2.0 -->
<!-- Origin: original -->
<!-- Source: Nazuchem | docs/wasm-size.md | Wasm size experiment -->

# Wasm size experiment

The maintained builder hard-codes V4. The historical harness, frozen corpus,
complete samples, and original variant builder are preserved by the
`wasm-size-experiment` git tag. Run the commands below in a checkout of that tag;
variant flags and the benchmark harness are intentionally absent from current
package sources. Relative `bench/` paths in the method refer to
`bindings/wasm/bench/` at that tag.

The benchmark runs with Bun as an orchestrator and real Node processes for all
timed work. From `bindings/wasm` in that checkout, preserve a V0 measurement
before trying other variants:

```sh
bun scripts/build.ts --variant=V0
bun run bench --label V0 --warmup 5 --output .generated/V0.json
bun scripts/build.ts --variant=V4
bun run bench --label V4 --warmup 5 --compare .generated/V0.json --output .generated/V4.json
bun run build # restore the selected default after experiments
```

Repeat the V4 commands with each of V1–V7. `--variant=V0` restores the original
plain release build with three Wasm copies; V1 only shares the binary; V2 also
removes `name`/`producers`; V3 adds the custom profile; V4/V5 add Binaryen -O3/-O4;
V6/V7 override the custom profile's Rust optimization level to s/z and use -Os/-Oz.
The optional Binaryen passes require pinned version 133; set `WASM_OPT` to its
`wasm-opt` executable. Other versions fail the build rather than silently change
the experiment. Official release assets and checksums are available from
[Binaryen version_133](https://github.com/WebAssembly/binaryen/releases/tag/version_133).

Each reported metric is the median of seven samples, with the full min–max range.
Cold samples use seven fresh Node processes and time `require()` (including Wasm
compilation/instantiation) through the first `nameToSmiles("ethanol")` result.
Process startup is excluded; the OS filesystem cache is uncontrolled and the
harness fingerprints build files before timing. Steady samples run in one
separate process after five untimed complete passes per direction. Each sample
calls all 300 names or all 300 SMILES, alternating direction order between runs.
The corpus is frozen before variants, stratified across golden/s2n fixture
categories, with fixed structural caps documented in `bench/corpus.ts`.
The benchmark rejects failed cases, unstable results, and cross-variant output,
runtime, corpus, or warmup/sample-count differences. Complete outputs are hashed
outside timed regions. The benchmark data retains the fixture provenance and is
excluded from npm packages. This is a bounded corpus measurement on the recorded
host/runtime, not a performance claim about every input or JavaScript engine.

Measured 2026-11-09 on Apple M1 / macOS arm64 (Darwin 27), Rust 1.99.1,
wasm-bindgen 0.2.128, Node 26.12.1 / V8 14.6.202.34-node.40, Bun 1.4.4,
and Binaryen 133 for V4–V7. Sizes below are **per Wasm payload**, except the
explicit shipped-total column; V0 includes three identical copies. gzip/Brotli
measure external transfer compression, not an embedded resource-format change.
All timing cells are milliseconds, median of seven with full min–max spread;
lower batch time means higher throughput. Builds and tests were excluded from
the serial timing windows. Initial V0 and every initial variant measurement were
retained; the complete matrix was repeated after correcting the bundler layout,
with unchanged timed Node and Wasm hashes.

| Variant | Recipe | Copies | Raw bytes / Wasm | gzip-9 bytes / Wasm | Brotli-11 bytes / Wasm | Shipped Wasm raw bytes |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| V0 | current release, three copies | 3 | 24,634,099 | 5,720,815 | 3,237,045 | 73,902,297 |
| V1 | shared Wasm only | 1 | 24,634,099 | 5,720,815 | 3,237,045 | 24,634,099 |
| V2 | V1 + strip name/producers | 1 | 20,346,641 | 5,322,977 | 3,055,049 | 20,346,641 |
| V3 | V2 + fat LTO / 1 CGU / Rust O3 / abort | 1 | 19,270,586 | 5,149,696 | 3,087,360 | 19,270,586 |
| V4 | V3 + wasm-opt -O3 | 1 | 18,621,450 | 5,147,532 | 3,163,236 | 18,621,450 |
| V5 | V3 + wasm-opt -O4 | 1 | 18,639,704 | 5,155,284 | 3,173,561 | 18,639,704 |
| V6 | V3 with Rust s + wasm-opt -Os | 1 | 14,564,831 | 4,148,804 | 2,549,090 | 14,564,831 |
| V7 | V3 with Rust z + wasm-opt -Oz | 1 | 13,352,715 | 3,645,153 | 2,218,378 | 13,352,715 |

| Variant | Cold ms, median [min–max] | 300 names ms, median [min–max] | 300 SMILES ms, median [min–max] | Eligible vs V0 |
| --- | ---: | ---: | ---: | --- |
| V0 | 303.45 [299.76–306.97] | 87.79 [86.66–89.14] | 526.42 [520.72–530.06] | yes |
| V1 | 301.03 [299.77–305.71] | 87.83 [87.16–89.33] | 526.83 [521.81–533.05] | yes |
| V2 | 303.84 [298.28–305.45] | 88.20 [86.86–89.84] | 527.08 [521.48–538.35] | yes |
| V3 | 301.57 [298.12–305.95] | 86.85 [85.75–88.44] | 534.96 [533.58–541.01] | yes |
| V4 | 312.73 [311.08–314.17] | 83.76 [82.72–84.70] | 499.20 [498.75–506.03] | no: cold |
| V5 | 312.67 [312.16–313.67] | 83.94 [82.89–84.93] | 499.85 [498.50–504.03] | no: cold |
| V6 | 143.80 [143.41–145.25] | 98.40 [97.63–99.49] | 666.14 [664.41–673.18] | no: nameToSmiles, smilesToName |
| V7 | 203.46 [202.71–204.00] | 135.22 [134.04–144.19] | 874.40 [868.01–943.66] | no: nameToSmiles, smilesToName |

**Decision: V4 is the default** (owner's choice, 2026-11-09). The rule was to
shrink the package without making naming slower. V4 is V3 plus Binaryen
`wasm-opt -O3`: one 18,621,450-byte Wasm (5,147,532 bytes gzip), down from three
24.6 MB copies (73,902,297 bytes, 17,162,445 gzip). Naming is faster: −4.6% for
300 names and −5.2% for 300 SMILES against final V0. The one cost is the first
load (instantiate plus first call), about 3% slower (+3.06%, ≈10 ms), at the
edge of the measured noise; the automatic threshold rule would have picked V3
for that reason. The noise thresholds were fixed from the initial V0 run before
optimization: `max(3%, (max − min) / median)` = **3.00% cold, 3.83% names,
3.08% SMILES**. The build requires Binaryen 133 (`WASM_OPT` may point at its
`wasm-opt`); CI installs the pinned release.

V3 (no wasm-opt) keeps first-load time unchanged at 19,270,586 bytes; V4/V5's
first load is +3.06%/+3.04% in the final matrix and +4.14%/+4.22% in the
independent initial matrix. V6/V7 materially
regress both steady directions despite their smaller payloads and faster cold
start. V6/V7 were experiment options; the maintained build uses only V4. `panic = "abort"` is scoped to the explicitly
selected Wasm profile and passes the binding/error-handling gates.
No native/Python release profile or embedded data format changes.

Individual samples, complete output checksums, build fingerprints, initial-run
metrics, tooling pins, and the decision are recorded in
[`bindings/wasm/bench/results.json` at `wasm-size-experiment`](https://github.com/Ameyanagi/nazuchem/blob/wasm-size-experiment/bindings/wasm/bench/results.json). The corpus checksum is
`a71ec3cb59ef70a8b4a3f5e4278041ac603687dcef75f53f986cee98aa6c5fa8`.

