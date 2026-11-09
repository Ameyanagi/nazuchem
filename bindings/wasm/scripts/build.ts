// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | bindings/wasm/scripts/build.ts | Phase B step 9

import { copyFile, mkdir, readFile, rm, writeFile } from "node:fs/promises";
import { resolve } from "node:path";
import { brotliCompressSync, constants, gzipSync } from "node:zlib";

const root = resolve(import.meta.dir, "../../..");
const packageRoot = resolve(import.meta.dir, "..");
const targetDir = resolve(root, process.env.CARGO_TARGET_DIR ?? "target");
const cli = process.env.WASM_BINDGEN ?? "wasm-bindgen";
const cliVersion = "0.2.128";
const wasmOpt = process.env.WASM_OPT ?? "wasm-opt";
const binaryenVersion = "133";
// Keep the full experiment reproducible. V4 is the owner-chosen default (2026-11-09):
// smallest variant that is not slower at naming (5% faster), accepting a ~3% slower first load.
const variants = {
  V0: { shared: false, strip: false, profile: "release" },
  V1: { shared: true, strip: false, profile: "release" },
  V2: { shared: true, strip: true, profile: "release" },
  V3: { shared: true, strip: true, profile: "wasm-release" },
  V4: { shared: true, strip: true, profile: "wasm-release", optimize: "-O3" },
  V5: { shared: true, strip: true, profile: "wasm-release", optimize: "-O4" },
  V6: { shared: true, strip: true, profile: "wasm-release", optLevel: "s", optimize: "-Os" },
  V7: { shared: true, strip: true, profile: "wasm-release", optLevel: "z", optimize: "-Oz" },
} satisfies Record<string, { shared: boolean; strip: boolean; profile: string; optLevel?: string; optimize?: string }>;
type Variant = keyof typeof variants;
let variantName = "V4";
for (let index = 2; index < process.argv.length; index++) {
  const arg = process.argv[index]!;
  if (arg === "--skip-cargo") continue;
  if (arg === "--variant") variantName = process.argv[++index] ?? "";
  else if (arg.startsWith("--variant=")) variantName = arg.slice("--variant=".length);
  else throw new Error(`Unknown build option: ${arg}`);
}
if (!Object.hasOwn(variants, variantName)) throw new Error(`Unknown size experiment variant: ${variantName}`);
const variant: { shared: boolean; strip: boolean; profile: string; optLevel?: string; optimize?: string } = variants[variantName as Variant];
if (variant.optimize) {
  const version = Bun.spawnSync([wasmOpt, "--version"], { stdout: "pipe", stderr: "pipe" });
  if (version.exitCode !== 0 || !new RegExp(`^wasm-opt version ${binaryenVersion}(?:\\s|$)`).test(version.stdout.toString().trim())) {
    throw new Error(`Expected Binaryen ${binaryenVersion}; set WASM_OPT to its pinned wasm-opt binary`);
  }
}

function run(cmd: string[], cwd = root): void {
  const result = Bun.spawnSync(cmd, {
    cwd, env: { ...process.env, CARGO_TARGET_DIR: targetDir }, stdout: "inherit", stderr: "inherit",
  });
  if (result.exitCode !== 0) throw new Error(`Command failed: ${cmd.join(" ")}`);
}

const lock = await readFile(resolve(root, "Cargo.lock"), "utf8");
const lockedVersion = lock.match(/name = "wasm-bindgen"\nversion = "([^"]+)"/)?.[1];
if (lockedVersion !== cliVersion) throw new Error(`Update the exact CLI pin: Cargo.lock has ${lockedVersion}`);
let version;
try {
  version = Bun.spawnSync([cli, "--version"], { stdout: "pipe", stderr: "pipe" });
} catch {
  throw new Error(`Install the pinned CLI: cargo install wasm-bindgen-cli --version ${cliVersion} --locked`);
}
if (version.exitCode !== 0 || version.stdout.toString().trim() !== `wasm-bindgen ${cliVersion}`) {
  throw new Error(`Expected wasm-bindgen ${cliVersion}; install with cargo install wasm-bindgen-cli --version ${cliVersion} --locked`);
}

if (!process.argv.includes("--skip-cargo")) {
  run(["cargo", "build", "-p", "nazuchem-wasm", "--target", "wasm32-unknown-unknown", "--profile", variant.profile, "--locked",
    ...(variant.optLevel ? ["--config", `profile.wasm-release.opt-level="${variant.optLevel}"`] : []),
  ]);
}
run(["uv", "run", "--no-project", "python", "scripts/collect_binding_notices.py", "--package", "nazuchem-wasm", "--target", "wasm32-unknown-unknown", "--output-dir", "bindings/wasm"]);

const dist = resolve(packageRoot, "dist");
await rm(dist, { recursive: true, force: true });
const header = "// SPDX-License-Identifier: MIT OR Apache-2.0\n// Origin: original\n// Source: Nazuchem | bindings/wasm/scripts/build.ts | Phase B step 9\n";
const declarations = await readFile(resolve(packageRoot, "src/index.d.ts"), "utf8");
const schemaText = await readFile(resolve(root, "crates/nazuchem-ir/schema/nazuchem-ir-1.schema.json"), "utf8");
const schema: { properties: { nazuchem_ir: { const: number } } } = JSON.parse(schemaText);
const wrapper = (await readFile(resolve(packageRoot, "src/wrapper.js"), "utf8"))
  .replace("__NAZUCHEM_SCHEMA__", JSON.stringify(schemaText))
  .replace("__NAZUCHEM_IR_VERSION__", String(schema.properties.nazuchem_ir.const));
const wasmPath = resolve(targetDir, `wasm32-unknown-unknown/${variant.profile}/nazuchem_wasm.wasm`);
const sharedPath = resolve(dist, "nazuchem_bg.wasm");
let sharedBytes: Buffer | undefined;

function redirect(content: string, from: string, to: string, count = 1): string {
  if (content.split(from).length !== count + 1) throw new Error(`Expected ${count} generated references: ${from}`);
  return content.replaceAll(from, to);
}

for (const [directory, target] of [["bundler", "bundler"], ["web", "web"], ["node", "nodejs"]] as const) {
  const out = resolve(dist, directory);
  await mkdir(out, { recursive: true });
  run([cli, wasmPath, "--out-dir", out, "--out-name", "nazuchem", "--target", target, "--no-typescript",
    ...(variant.strip ? ["--remove-name-section", "--remove-producers-section"] : []),
  ]);
  if (variant.shared) {
    const localWasm = resolve(out, "nazuchem_bg.wasm");
    const bytes = await readFile(localWasm);
    if (sharedBytes && !bytes.equals(sharedBytes)) throw new Error(`Cannot share different Wasm payloads: ${directory}`);
    if (!sharedBytes) {
      sharedBytes = bytes;
      await writeFile(sharedPath, bytes);
    }
    const gluePath = resolve(out, "nazuchem.js");
    let glue = await readFile(gluePath, "utf8");
    if (target === "nodejs") glue = redirect(glue, '${__dirname}/nazuchem_bg.wasm', '${__dirname}/../nazuchem_bg.wasm');
    else if (target === "web") glue = redirect(glue, "new URL('nazuchem_bg.wasm', import.meta.url)", "new URL('../nazuchem_bg.wasm', import.meta.url)");
    else {
      glue = redirect(glue, '"./nazuchem_bg.wasm"', '"../nazuchem_bg.wasm"');
      // Static Wasm imports resolve relative to the Wasm file, so its bundler
      // import module must move alongside it. Both JS imports use this instance.
      glue = redirect(glue, '"./nazuchem_bg.js"', '"../nazuchem_bg.js"', 2);
      await copyFile(resolve(out, "nazuchem_bg.js"), resolve(dist, "nazuchem_bg.js"));
      await rm(resolve(out, "nazuchem_bg.js"));
    }
    await writeFile(gluePath, glue);
    await rm(localWasm);
  }
  const isNode = target === "nodejs";
  const load = isNode ? 'const raw = require("./nazuchem.js");' : 'import * as raw from "./nazuchem.js";';
  let content = header + load + "\n" + wrapper;
  if (isNode) {
    content = content.replaceAll("export ", "") + "\nmodule.exports = { NazuchemError, IR_VERSION, JSON_SCHEMA, nameToSmiles, smilesToName, nameToMolecule, moleculeToName };\n";
    await writeFile(resolve(out, "package.json"), '{"type":"commonjs"}\n');
  } else {
    // Browser initialization is explicit. Bundlers instantiate their Wasm import.
    if (target === "web") content += '\nexport { default, initSync } from "./nazuchem.js";\n';
    await writeFile(resolve(out, "package.json"), '{"type":"module"}\n');
  }
  await writeFile(resolve(out, isNode ? "index.cjs" : "index.js"), content);
  let types = declarations;
  if (target === "web") {
    types += await readFile(resolve(packageRoot, "src/web.d.ts"), "utf8");
  }
  await writeFile(resolve(out, isNode ? "index.d.cts" : "index.d.ts"), types);
}
await copyFile(resolve(root, "crates/nazuchem-ir/schema/nazuchem-ir-1.schema.json"), resolve(dist, "nazuchem-ir-1.schema.json"));

if (variant.optimize) {
  run([wasmOpt, sharedPath, variant.optimize, "--strip-debug", "--strip-producers", "-o", `${sharedPath}.optimized`]);
  await copyFile(`${sharedPath}.optimized`, sharedPath);
  await rm(`${sharedPath}.optimized`);
}
const bytes = await readFile(variant.shared ? sharedPath : resolve(dist, "node/nazuchem_bg.wasm"));
console.log(`${variantName}: ${variant.shared ? "one shared Wasm" : "three Wasm copies"}; ${bytes.length} raw bytes; ${gzipSync(bytes, { level: 9 }).length} gzip bytes; ${brotliCompressSync(bytes, { params: { [constants.BROTLI_PARAM_QUALITY]: 11 } }).length} brotli bytes per payload`);
