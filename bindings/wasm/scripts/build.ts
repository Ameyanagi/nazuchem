// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | bindings/wasm/scripts/build.ts | Phase B step 9

import { copyFile, mkdir, readFile, rm, writeFile } from "node:fs/promises";
import { resolve } from "node:path";
import { gzipSync } from "node:zlib";

const root = resolve(import.meta.dir, "../../..");
const packageRoot = resolve(import.meta.dir, "..");
const targetDir = process.env.CARGO_TARGET_DIR ?? resolve(root, "target");
const cli = process.env.WASM_BINDGEN ?? "wasm-bindgen";
const cliVersion = "0.2.128";

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
  run(["cargo", "build", "-p", "nazuchem-wasm", "--target", "wasm32-unknown-unknown", "--release", "--locked"]);
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
const wasmPath = resolve(targetDir, "wasm32-unknown-unknown/release/nazuchem_wasm.wasm");

for (const [directory, target] of [["bundler", "bundler"], ["web", "web"], ["node", "nodejs"]] as const) {
  const out = resolve(dist, directory);
  await mkdir(out, { recursive: true });
  run([cli, wasmPath, "--out-dir", out, "--out-name", "nazuchem", "--target", target, "--no-typescript"]);
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
  const bytes = await readFile(resolve(out, "nazuchem_bg.wasm"));
  console.log(`${directory}: ${bytes.length} raw bytes; ${gzipSync(bytes, { level: 9 }).length} gzip bytes (default n2s+s2n, release, no wasm-opt)`);
}
await copyFile(resolve(root, "crates/nazuchem-ir/schema/nazuchem-ir-1.schema.json"), resolve(dist, "nazuchem-ir-1.schema.json"));
