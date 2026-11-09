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
// V4: one shared, stripped Wasm; fat LTO / one CGU / Rust O3 / abort,
// followed by Binaryen -O3. Historical measurements live in docs/wasm-size.md.
for (const arg of process.argv.slice(2)) {
  if (arg !== "--skip-cargo") throw new Error(`Unknown build option: ${arg}`);
}
const optimizerVersion = Bun.spawnSync([wasmOpt, "--version"], { stdout: "pipe", stderr: "pipe" });
if (optimizerVersion.exitCode !== 0 || !new RegExp(`^wasm-opt version ${binaryenVersion}(?:\\s|$)`).test(optimizerVersion.stdout.toString().trim())) {
  throw new Error(`Expected Binaryen ${binaryenVersion}; set WASM_OPT to its pinned wasm-opt binary`);
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
  throw new Error(`Install the pinned wasm-bindgen ${cliVersion} release binary and set WASM_BINDGEN`);
}
if (version.exitCode !== 0 || version.stdout.toString().trim() !== `wasm-bindgen ${cliVersion}`) {
  throw new Error(`Expected wasm-bindgen ${cliVersion}; set WASM_BINDGEN to the matching release binary`);
}

if (!process.argv.includes("--skip-cargo")) {
  run(["cargo", "build", "-p", "nazuchem-wasm", "--target", "wasm32-unknown-unknown", "--profile", "wasm-release", "--locked",
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
const wasmPath = resolve(targetDir, "wasm32-unknown-unknown/wasm-release/nazuchem_wasm.wasm");
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
    "--remove-name-section", "--remove-producers-section",
  ]);
  {
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
  if (target === "web") {
    await writeFile(resolve(out, "index.d.cts"), types);
    // CJS consumers get the same explicit browser loader. Preserve a runtime
    // asset URL: Bun otherwise substitutes the builder's absolute source path.
    const cjs = await Bun.build({
      entrypoints: [resolve(out, "index.js")], target: "node", format: "cjs",
      outdir: out, naming: "index.cjs",
      banner: header + 'const __nazuchem_web_url = require("node:url").pathToFileURL(__filename).href;',
      plugins: [{ name: "relative-web-wasm", setup(build) {
        build.onLoad({ filter: /nazuchem\.js$/ }, async ({ path }) => ({
          contents: redirect(await readFile(path, "utf8"), "import.meta.url", "__nazuchem_web_url"),
          loader: "js",
        }));
      } }],
    });
    if (!cjs.success) throw new Error(cjs.logs.map(String).join("\n"));
  }
}
await copyFile(resolve(root, "crates/nazuchem-ir/schema/nazuchem-ir-1.schema.json"), resolve(dist, "nazuchem-ir-1.schema.json"));

run([wasmOpt, sharedPath, "-O3", "--strip-debug", "--strip-producers", "-o", `${sharedPath}.optimized`]);
await copyFile(`${sharedPath}.optimized`, sharedPath);
await rm(`${sharedPath}.optimized`);
const bytes = await readFile(sharedPath);
console.log(`V4: one shared Wasm; ${bytes.length} raw bytes; ${gzipSync(bytes, { level: 9 }).length} gzip bytes; ${brotliCompressSync(bytes, { params: { [constants.BROTLI_PARAM_QUALITY]: 11 } }).length} brotli bytes`);
