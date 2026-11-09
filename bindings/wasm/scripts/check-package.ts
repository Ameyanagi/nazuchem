// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | bindings/wasm/scripts/check-package.ts | Phase B step 9

import { mkdir, mkdtemp, readFile, rm } from "node:fs/promises";
import { tmpdir } from "node:os";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../../..");
const packageRoot = resolve(import.meta.dir, "..");
const archive = resolve(packageRoot, ".generated/nazuchem-package.tgz");
await mkdir(resolve(packageRoot, ".generated"), { recursive: true });
function run(command: string[], cwd = packageRoot): Buffer {
  const result = Bun.spawnSync(command, { cwd, stdout: "pipe", stderr: "pipe" });
  if (result.exitCode !== 0) throw new Error(`${command.join(" ")}\n${result.stdout}\n${result.stderr}`);
  return result.stdout;
}
run(["bun", "pm", "pack", "--ignore-scripts", "--filename", archive]);
const scratch = await mkdtemp(resolve(tmpdir(), "nazuchem-wasm-package-"));
try {
  // The archive must match an independent, nonempty target-specific collection,
  // never merely the ignored local files that were packed moments earlier.
  const notices = resolve(scratch, "notices");
  run(["uv", "run", "--no-project", "python", "scripts/collect_binding_notices.py", "--package", "nazuchem-wasm", "--target", "wasm32-unknown-unknown", "--output-dir", notices], root);
  for (const name of ["LICENSE-MIT", "LICENSE-APACHE", "NOTICE", "THIRD-PARTY-NOTICES.txt"]) {
    const fresh = await readFile(resolve(notices, name));
    if (fresh.length === 0) throw new Error(`Fresh ${name} is empty`);
    const packed = run(["tar", "-xOf", archive, `package/${name}`]);
    if (!packed.equals(fresh)) throw new Error(`Packed ${name} differs from fresh wasm32 notices`);
  }
  const installed = resolve(scratch, "node_modules/nazuchem");
  await mkdir(installed, { recursive: true });
  run(["tar", "-xzf", archive, "--strip-components=1", "-C", installed]);
  const tests = Bun.spawnSync(["bun", "test", "tests"], {
    cwd: packageRoot, env: { ...process.env, NAZUCHEM_PACKAGE_ROOT: installed },
    stdout: "inherit", stderr: "inherit",
  });
  if (tests.exitCode !== 0) throw new Error("Packed runtime tests failed");
  console.log(`Verified packed runtime APIs and four fresh, nonempty notices: ${archive}`);
} finally {
  await rm(scratch, { recursive: true, force: true });
}
