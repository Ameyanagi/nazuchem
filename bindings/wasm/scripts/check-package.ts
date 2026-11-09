// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | bindings/wasm/scripts/check-package.ts | Phase B step 9

import { mkdir, readFile } from "node:fs/promises";
import { resolve } from "node:path";

const packageRoot = resolve(import.meta.dir, "..");
const archive = resolve(packageRoot, ".generated/nazuchem-package.tgz");
await mkdir(resolve(packageRoot, ".generated"), { recursive: true });
function run(command: string[]): string {
  const result = Bun.spawnSync(command, { cwd: packageRoot, stdout: "pipe", stderr: "pipe" });
  if (result.exitCode !== 0) throw new Error(result.stderr.toString());
  return result.stdout.toString();
}
run(["bun", "pm", "pack", "--ignore-scripts", "--filename", archive]);
const entries = run(["tar", "-tzf", archive]).split("\n").filter(Boolean);
for (const name of ["LICENSE-MIT", "LICENSE-APACHE", "NOTICE", "THIRD-PARTY-NOTICES.txt"]) {
  const path = `package/${name}`;
  if (!entries.includes(path)) throw new Error(`Package is missing ${path}`);
  const packed = run(["tar", "-xOf", archive, path]);
  if (packed !== await readFile(resolve(packageRoot, name), "utf8")) throw new Error(`Package changed ${name}`);
}
for (const target of ["bundler", "web", "node"]) {
  if (!entries.includes(`package/dist/${target}/nazuchem_bg.wasm`)) throw new Error(`Package missing ${target} Wasm`);
}
if (entries.some((entry) => /\/(node_modules|tests|scripts|src)\//.test(entry))) throw new Error("Development source leaked into package");
const manifest: { private: boolean } = JSON.parse(run(["tar", "-xOf", archive, "package/package.json"]));
if (manifest.private !== true) throw new Error("npm publishing must remain disabled");
console.log(`Verified real bun pm pack archive: ${entries.length} files, ${(await Bun.file(archive).arrayBuffer()).byteLength} bytes, all four notices intact.`);
