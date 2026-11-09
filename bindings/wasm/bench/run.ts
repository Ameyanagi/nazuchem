// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | bindings/wasm/bench/run.ts | npm/wasm size benchmark

import { createHash } from "node:crypto";
import { mkdir, readFile, readdir, stat, writeFile } from "node:fs/promises";
import { dirname, join, relative, resolve } from "node:path";
import { fileURLToPath } from "node:url";
import { arch, cpus, platform, release } from "node:os";
import { corpus } from "./corpus.ts";

// Bun orchestrates only. Node's CommonJS loader and V8/WebAssembly runtime do
// every timed operation. Keeping the worker as plain JS avoids a Node TS loader.
const worker = String.raw`
const { createHash } = require("node:crypto");
const { readFileSync } = require("node:fs");
const { performance } = require("node:perf_hooks");
if (process.versions.bun) throw new Error("Benchmark worker must run in Node");
const [entry, mode, runCount, warmupCount] = process.argv.slice(1);
const hash = (value) => createHash("sha256").update(JSON.stringify(value)).digest("hex");
const engine = { node: process.version, v8: process.versions.v8, arch: process.arch, platform: process.platform };
if (mode === "cold") {
  const start = performance.now();
  const api = require(entry);
  const value = api.nameToSmiles("ethanol");
  const ms = performance.now() - start;
  if (value !== "C(C)O") throw new Error("Unexpected ethanol result: " + value);
  process.stdout.write(JSON.stringify({ engine, ms, checksum: hash(value) }));
} else {
  const inputs = JSON.parse(readFileSync(0, "utf8"));
  const api = require(entry);
  const methods = ["nameToSmiles", "smilesToName"];
  const outputs = {};
  const samples = { nameToSmiles: [], smilesToName: [] };
  const checksums = {};
  function pass(method, measured) {
    const values = inputs[method];
    const results = new Array(values.length);
    const start = performance.now();
    for (let i = 0; i < values.length; i++) {
      try { results[i] = api[method](values[i]); }
      catch (error) { throw new Error(method + " case " + i + " (" + values[i] + "): " + error.message, { cause: error }); }
    }
    const ms = performance.now() - start;
    const checksum = hash(results);
    if (checksums[method] && checksums[method] !== checksum) throw new Error(method + " output changed between passes");
    checksums[method] = checksum;
    if (measured) samples[method].push(ms);
    outputs[method] = results.length;
  }
  for (let i = 0; i < Number(warmupCount); i++) for (const method of methods) pass(method, false);
  for (let i = 0; i < Number(runCount); i++) {
    // Alternate which direction goes first, while preserving fixed corpus order.
    for (const method of i % 2 ? [...methods].reverse() : methods) pass(method, true);
  }
  process.stdout.write(JSON.stringify({ engine, samples, checksums, cases: outputs }));
}
`;

function usage(): never {
  console.log(`Usage: bun bindings/wasm/bench/run.ts [options]
  --dist PATH       Distribution containing node/index.cjs (default: bindings/wasm/dist)
  --output PATH     Save all samples and statistics as JSON (default: stdout only)
  --compare PATH    Compare with an earlier benchmark JSON; require identical outputs
  --label NAME      Optional variant label
  --runs N          Recorded samples per metric (default: 7)
  --warmup N        Full untimed corpus passes per direction (default: 2)
  --node PATH       Node executable (default: node found on PATH)
  --help

Cold: seven fresh Node processes, timer immediately before require through the
first nameToSmiles("ethanol") call. Steady: one separate Node process, warm both
300-case directions, then record seven complete passes per direction. Output
hashing and JSON encoding occur outside timers. No variant-dependent filtering.
Do not run builds or other benchmarks concurrently with recorded measurements.`);
  process.exit(0);
}

const options: Record<string, string> = {};
for (let i = 2; i < process.argv.length; i++) {
  const arg = process.argv[i]!;
  if (arg === "--help" || arg === "-h") usage();
  if (!["--dist", "--output", "--compare", "--label", "--runs", "--warmup", "--node"].includes(arg)) throw new Error(`Unknown option ${arg}`);
  const value = process.argv[++i];
  if (!value || value.startsWith("--")) throw new Error(`Missing value for ${arg}`);
  options[arg.slice(2)] = value;
}
const bindingRoot = resolve(dirname(fileURLToPath(import.meta.url)), "..");
const dist = resolve(options.dist ?? join(bindingRoot, "dist"));
const entry = join(dist, "node", "index.cjs");
const runs = Number(options.runs ?? 7);
const warmup = Number(options.warmup ?? 2);
if (!Number.isInteger(runs) || runs < 1 || !Number.isInteger(warmup) || warmup < 1) throw new Error("runs and warmup must be positive integers");
await stat(entry);
if (corpus.length !== 300 || new Set(corpus.map((row) => row.name)).size !== 300 || new Set(corpus.map((row) => row.smiles)).size !== 300) throw new Error("Expected 300 unique cases per direction");
const inputs = { nameToSmiles: corpus.map((row) => row.name), smilesToName: corpus.map((row) => row.smiles) };
const sha256 = (value: string | Uint8Array) => createHash("sha256").update(value).digest("hex");
const corpusSha256 = sha256(JSON.stringify(inputs));

interface Engine { node: string; v8: string; arch: string; platform: string }
interface ColdResult { engine: Engine; ms: number; checksum: string }
interface SteadyResult { engine: Engine; samples: Record<string, number[]>; checksums: Record<string, string>; cases: Record<string, number> }
async function runWorker<T>(mode: string): Promise<T> {
  const child = Bun.spawn([options.node ?? "node", "--input-type=commonjs", "-e", worker, entry, mode, String(runs), String(warmup)], {
    stdin: "pipe", stdout: "pipe", stderr: "pipe",
  });
  if (mode !== "cold") child.stdin.write(JSON.stringify(inputs));
  child.stdin.end();
  const [stdout, stderr, exitCode] = await Promise.all([new Response(child.stdout).text(), new Response(child.stderr).text(), child.exited]);
  if (exitCode !== 0) throw new Error(`Node ${mode} worker exited ${exitCode}:\n${stderr}\n${stdout}`);
  if (stderr) process.stderr.write(stderr);
  return JSON.parse(stdout) as T;
}

function statistics(samplesMs: number[]) {
  if (samplesMs.length !== runs || samplesMs.some((sample) => !Number.isFinite(sample) || sample <= 0)) throw new Error("Invalid timing samples");
  const sorted = [...samplesMs].sort((a, b) => a - b);
  const medianMs = sorted.length % 2 ? sorted[(sorted.length - 1) / 2]! : (sorted[sorted.length / 2 - 1]! + sorted[sorted.length / 2]!) / 2;
  const minMs = sorted[0]!;
  const maxMs = sorted[sorted.length - 1]!;
  return { samplesMs, medianMs, minMs, maxMs, spreadMs: maxMs - minMs, spreadPercent: 100 * (maxMs - minMs) / medianMs };
}

// Fingerprint the complete supplied build before timing; includes wrapper/glue
// changes as well as Wasm. This is intentionally outside every timed interval.
async function buildFiles(directory: string): Promise<{ path: string; bytes: number; sha256: string }[]> {
  const files = [];
  for (const child of (await readdir(directory, { withFileTypes: true })).sort((a, b) => a.name.localeCompare(b.name))) {
    const path = join(directory, child.name);
    if (child.isDirectory()) files.push(...await buildFiles(path));
    else if (child.isFile()) { const bytes = await readFile(path); files.push({ path: relative(dist, path), bytes: bytes.length, sha256: sha256(bytes) }); }
    else throw new Error(`Unsupported non-file build entry: ${path}`);
  }
  return files;
}
const files = await buildFiles(dist);
const cold: ColdResult[] = [];
for (let i = 0; i < runs; i++) {
  process.stderr.write(`cold ${i + 1}/${runs}\n`);
  cold.push(await runWorker<ColdResult>("cold"));
}
process.stderr.write(`steady: ${warmup} warmup passes + ${runs} measured passes, 300 cases per direction\n`);
const steady = await runWorker<SteadyResult>("steady");
if (cold.some((sample) => sample.checksum !== cold[0]!.checksum || JSON.stringify(sample.engine) !== JSON.stringify(steady.engine))) throw new Error("Inconsistent cold worker result/runtime");
const metrics = {
  cold: statistics(cold.map((sample) => sample.ms)),
  nameToSmiles: { ...statistics(steady.samples.nameToSmiles!), casesPerSample: corpus.length },
  smilesToName: { ...statistics(steady.samples.smilesToName!), casesPerSample: corpus.length },
};
const checksums = { cold: cold[0]!.checksum, ...steady.checksums };
const result = {
  schema: "nazuchem-wasm-benchmark-v1", label: options.label ?? null, recordedAt: new Date().toISOString(),
  dist, engine: steady.engine,
  host: { platform: platform(), release: release(), arch: arch(), cpu: cpus()[0]?.model ?? "unknown" },
  orchestrator: { bun: Bun.version },
  corpus: { version: 1, sha256: corpusSha256, casesPerDirection: corpus.length, categories: new Set(corpus.map((row) => `${row.source.split(":")[0]}:${row.category}`)).size },
  methodology: { runs, warmupPassesPerDirection: warmup, coldFreshProcesses: runs, steadyFreshProcesses: 1, coldIncludesProcessStartup: false, coldIncludesRequireAndFirstEthanolCall: true, coldFilesystemCache: "uncontrolled; build fingerprinting reads distribution before samples", steadyIncludesOutputHashing: false, steadyOrder: "fixed corpus order; alternate direction order between samples" },
  checksums, metrics, files,
};
if (options.compare) {
  const baseline = JSON.parse(await readFile(resolve(options.compare), "utf8"));
  if (baseline.schema !== result.schema || baseline.corpus?.sha256 !== corpusSha256 || JSON.stringify(baseline.checksums) !== JSON.stringify(checksums)) throw new Error("Comparison requires identical corpus and complete API output checksums");
  if (JSON.stringify(baseline.engine) !== JSON.stringify(result.engine)) throw new Error("Comparison requires identical Node/V8 runtime");
  if (baseline.methodology?.runs !== runs || baseline.methodology?.warmupPassesPerDirection !== warmup) throw new Error("Comparison requires identical runs and warmup counts");
  const comparisons = Object.fromEntries(Object.entries(metrics).map(([metric, current]) => {
    const prior = baseline.metrics[metric];
    return [metric, { medianChangePercent: 100 * (current.medianMs / prior.medianMs - 1), baselineMedianMs: prior.medianMs, currentMedianMs: current.medianMs }];
  }));
  Object.assign(result, { comparison: { baseline: resolve(options.compare), metrics: comparisons } });
}
const json = JSON.stringify(result, null, 2) + "\n";
if (options.output) { const path = resolve(options.output); await mkdir(dirname(path), { recursive: true }); await writeFile(path, json); }
process.stdout.write(json);
