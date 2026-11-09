// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | bindings/wasm/scripts/check-package.ts | Phase B step 9

import { mkdir, readFile } from "node:fs/promises";
import { posix, resolve } from "node:path";

const packageRoot = resolve(import.meta.dir, "..");
const archive = resolve(packageRoot, ".generated/nazuchem-package.tgz");
await mkdir(resolve(packageRoot, ".generated"), { recursive: true });
function run(command: string[], maxBuffer?: number): Buffer {
  const result = Bun.spawnSync(command, { cwd: packageRoot, stdout: "pipe", stderr: "pipe", maxBuffer });
  if (result.exitCode !== 0) throw new Error(result.stderr.toString());
  return result.stdout;
}
run(["bun", "pm", "pack", "--ignore-scripts", "--filename", archive]);
const entries = run(["tar", "-tzf", archive]).toString().split("\n").filter(Boolean);
for (const name of ["LICENSE-MIT", "LICENSE-APACHE", "NOTICE", "THIRD-PARTY-NOTICES.txt"]) {
  const path = `package/${name}`;
  if (!entries.includes(path)) throw new Error(`Package is missing ${path}`);
  const packed = run(["tar", "-xOf", archive, path]);
  if (!packed.equals(await readFile(resolve(packageRoot, name)))) throw new Error(`Package changed ${name}`);
}
const sharedWasm = "package/dist/nazuchem_bg.wasm";
const sharedBundlerGlue = "package/dist/nazuchem_bg.js";
const wasmEntries = entries.filter((entry) => entry.endsWith(".wasm"));
if (wasmEntries.length !== 1 || wasmEntries[0] !== sharedWasm) {
  throw new Error(`Package must contain exactly one shared Wasm at ${sharedWasm}; found ${JSON.stringify(wasmEntries)}`);
}
if (!entries.includes(sharedBundlerGlue)) throw new Error(`Package is missing ${sharedBundlerGlue}`);
// Inspect only module metadata. Compiling the packed bytes neither evaluates
// JavaScript imports nor instantiates or executes the Wasm module.
const packedWasm = run(["tar", "-xOf", archive, sharedWasm], 128 * 1024 * 1024);
const wasmModule = new WebAssembly.Module(new Uint8Array(packedWasm));
const importedModules = [...new Set(WebAssembly.Module.imports(wasmModule).map((entry) => entry.module))];
if (importedModules.length !== 1) {
  throw new Error(`Shared Wasm must import one relative bundler glue module; found ${JSON.stringify(importedModules)}`);
}
for (const specifier of importedModules) {
  if (!specifier.startsWith("./") && !specifier.startsWith("../")) {
    throw new Error(`Shared Wasm has a nonrelative dependency: ${specifier}`);
  }
  const importedPath = posix.normalize(posix.join(posix.dirname(sharedWasm), specifier));
  if (importedPath !== sharedBundlerGlue || !entries.includes(importedPath)) {
    throw new Error(`Shared Wasm dependency resolves to ${importedPath}, expected packed ${sharedBundlerGlue}`);
  }
}
for (const target of ["bundler", "web", "node"] as const) {
  const gluePath = `package/dist/${target}/nazuchem.js`;
  const wrapperPath = `package/dist/${target}/${target === "node" ? "index.cjs" : "index.js"}`;
  for (const path of [gluePath, wrapperPath]) {
    if (!entries.includes(path)) throw new Error(`Package is missing ${path}`);
  }
  const glue = run(["tar", "-xOf", archive, gluePath]).toString();
  const pattern = target === "node"
    ? /`\$\{__dirname\}\/([^`\r\n]+\.wasm)`/g
    : target === "web"
      ? /new URL\(\s*["']([^"'\r\n]+\.wasm)["']\s*,\s*import\.meta\.url\s*\)/g
      : /\bfrom\s+["']([^"'\r\n]+\.wasm)["']/g;
  const references = [...glue.matchAll(pattern)].map((match) => match[1]);
  const wasmLiterals = [...glue.matchAll(/(["'`])([^"'`\r\n]+\.wasm)\1/g)];
  if (references.length !== 1 || wasmLiterals.length !== references.length) {
    throw new Error(`${gluePath}: expected exactly one recognizable generated Wasm reference`);
  }
  const referencedPath = posix.normalize(posix.join(posix.dirname(gluePath), references[0]));
  if (referencedPath !== sharedWasm || !entries.includes(referencedPath)) {
    throw new Error(`${gluePath}: Wasm reference resolves to ${referencedPath}, expected packed ${sharedWasm}`);
  }
  if (target === "bundler") {
    const jsReferences = [...glue.matchAll(/\bfrom\s+["']([^"'\r\n]+\.js)["']/g)].map((match) => match[1]);
    if (jsReferences.length !== 2) throw new Error(`${gluePath}: expected both generated bundler glue references`);
    for (const reference of jsReferences) {
      const referencedJs = posix.normalize(posix.join(posix.dirname(gluePath), reference));
      if (referencedJs !== sharedBundlerGlue || !entries.includes(referencedJs)) {
        throw new Error(`${gluePath}: JavaScript reference resolves to ${referencedJs}, expected packed ${sharedBundlerGlue}`);
      }
    }
  }
}
if (entries.some((entry) => /\/(node_modules|tests|scripts|src|bench)\//.test(entry))) throw new Error("Development source leaked into package");
const manifest: { private: boolean } = JSON.parse(run(["tar", "-xOf", archive, "package/package.json"]).toString());
if (manifest.private !== true) throw new Error("npm publishing must remain disabled");
console.log(`Verified real bun pm pack archive: ${entries.length} files, ${(await Bun.file(archive).arrayBuffer()).byteLength} bytes, one shared Wasm reached by all three targets, all four notices intact.`);
