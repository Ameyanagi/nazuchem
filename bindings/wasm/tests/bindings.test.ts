// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | bindings/wasm/tests/bindings.test.ts | Phase B step 9

import { expect, test } from "bun:test";
import Ajv2020 from "ajv/dist/2020.js";
import { createRequire } from "node:module";
import { readFile, realpath } from "node:fs/promises";
import { resolve } from "node:path";
import { fileURLToPath, pathToFileURL } from "node:url";
import type { Molecule, Name, Stereo } from "../src/index.d.ts";

// The identical suite runs against local dist and an extracted real tarball.
const packageRoot = await realpath(resolve(process.env.NAZUCHEM_PACKAGE_ROOT ?? fileURLToPath(new URL("..", import.meta.url))));
const packageRequire = createRequire(resolve(packageRoot, "package.json"));
const {
  IR_VERSION, JSON_SCHEMA, NazuchemError, moleculeToName,
  nameToMolecule, nameToSmiles, smilesToName,
}: typeof import("../dist/node/index.cjs") = packageRequire("nazuchem");

const ajv = new Ajv2020({ strict: true, validateFormats: false });
const validate = ajv.compile<Molecule>(JSON.parse(JSON_SCHEMA));

function assertValid(molecule: Molecule): void {
  expect(validate(molecule), JSON.stringify(validate.errors)).toBe(true);
  expect(Object.getPrototypeOf(molecule)).toBe(Object.prototype);
  for (const atom of molecule.atoms) expect(Object.getPrototypeOf(atom)).toBe(Object.prototype);
}

test("ethanol in both directions retains complete naming assessment", () => {
  expect(nameToSmiles("ethanol")).toBe("C(C)O");
  const name = smilesToName("CCO");
  expect(name).toEqual({ text: "ethanol", status: "preferred", reason: null, alternatives: [], omittedStereo: [], omittedUnits: 0, verification: null });
  const molecule = nameToMolecule("ethanol");
  assertValid(molecule);
  expect(molecule.atoms).toHaveLength(3);
  expect(moleculeToName(molecule)).toEqual(name);
  expect(moleculeToName(JSON.stringify(molecule))).toEqual(name);
  expect(smilesToName("CCCCO").text).toBe("butan-1-ol");
});

test("IR constants match the exported package schema byte for byte", async () => {
  expect(IR_VERSION).toBe(1);
  const schemaPath = packageRequire.resolve("nazuchem/schema.json");
  expect(JSON_SCHEMA).toBe(await readFile(schemaPath, "utf8"));
  expect(packageRequire("nazuchem/schema.json").properties.nazuchem_ir.const).toBe(IR_VERSION);
  expect(packageRequire("nazuchem/node").nameToSmiles("ethanol")).toBe("C(C)O");
});

test("bare root, Node subpath, and schema exports resolve in fresh Node consumers", () => {
  for (const [module, program] of [
    ["commonjs", `const assert = require("node:assert/strict");
const api = require("nazuchem");
assert.equal(api.nameToSmiles("ethanol"), "C(C)O");
assert.equal(require("nazuchem/node").smilesToName("CCO").text, "ethanol");
assert.equal(require("nazuchem/schema.json").properties.nazuchem_ir.const, 1);`],
    ["module", `import assert from "node:assert/strict";
import { nameToSmiles } from "nazuchem";
import { smilesToName } from "nazuchem/node";
assert.equal(nameToSmiles("ethanol"), "C(C)O");
assert.equal(smilesToName("CCO").text, "ethanol");`],
  ]) {
    const result = Bun.spawnSync(["node", `--input-type=${module}`, "-e", program!], {
      cwd: packageRoot, stdout: "pipe", stderr: "pipe", timeout: 30_000,
    });
    expect(result.exitCode, result.stderr.toString()).toBe(0);
  }
}, 30_000);

test("browser default initialization fetches the shared Wasm and runs the naming API", async () => {
  const sharedWasm = pathToFileURL(resolve(packageRoot, "dist/nazuchem_bg.wasm"));
  const fetchedUrls: string[] = [];
  const originalFetch = globalThis.fetch;
  globalThis.fetch = (async (input) => {
    const url = new URL(input instanceof Request ? input.url : input.toString());
    fetchedUrls.push(url.href);
    expect(url.href).toBe(sharedWasm.href);
    const bytes = await readFile(url);
    const body = new Uint8Array(bytes.length);
    body.set(bytes);
    return new Response(body, { headers: { "Content-Type": "application/wasm" } });
  }) as typeof fetch;
  try {
    const web: typeof import("../dist/web/index.js") = await import(pathToFileURL(resolve(packageRoot, "dist/web/index.js")).href);
    await web.default();
    expect(fetchedUrls).toEqual([sharedWasm.href]);
    expect(web.nameToSmiles("ethanol")).toBe("C(C)O");
    expect(web.nameToMolecule("ethanol").atoms).toHaveLength(3);
  } finally {
    globalThis.fetch = originalFetch;
  }
});

test("CommonJS web initialization supports default URL, supplied bytes, and initSync", () => {
  const program = `const assert = require("node:assert/strict");
const { readFileSync } = require("node:fs");
const { dirname, resolve } = require("node:path");
const { pathToFileURL } = require("node:url");
const schemaPath = require.resolve("nazuchem/schema.json");
const wasmPath = resolve(dirname(schemaPath), "nazuchem_bg.wasm");
const bytes = readFileSync(wasmPath);
(async () => {
  for (const mode of ["default", "async", "sync"]) {
    delete require.cache[require.resolve("nazuchem/web")];
    const web = require("nazuchem/web");
    let fetches = 0;
    globalThis.fetch = async (input) => {
      assert.equal(input.toString(), pathToFileURL(wasmPath).href);
      fetches++;
      return new Response(bytes, { headers: { "Content-Type": "application/wasm" } });
    };
    if (mode === "default") await web.default();
    else if (mode === "async") await web.default({ module_or_path: bytes });
    else web.initSync({ module: bytes });
    assert.equal(fetches, mode === "default" ? 1 : 0);
    assert.equal(web.nameToSmiles("ethanol"), "C(C)O");
    assert.equal(web.smilesToName("CCO").text, "ethanol");
    assert.equal(web.moleculeToName(web.nameToMolecule("ethanol")).text, "ethanol");
    assert.equal(web.IR_VERSION, 1);
    assert.equal(web.JSON_SCHEMA, readFileSync(schemaPath, "utf8"));
    assert.throws(() => web.nameToSmiles("not a chemical name"), (error) => error instanceof web.NazuchemError && error.kind === "parse");
  }
})().catch((error) => { console.error(error); process.exitCode = 1; });`;
  const result = Bun.spawnSync(["node", "--input-type=commonjs", "-e", program], {
    cwd: packageRoot, stdout: "pipe", stderr: "pipe", timeout: 30_000,
  });
  expect(result.exitCode, result.stderr.toString()).toBe(0);
}, 30_000);

test("bundler static Wasm imports run both naming directions in a fresh Node process", () => {
  const entry = pathToFileURL(resolve(packageRoot, "dist/bundler/index.js")).href;
  const program = [
    'import { strict as assert } from "node:assert";',
    `import { nameToSmiles, smilesToName } from ${JSON.stringify(entry)};`,
    'assert.equal(nameToSmiles("ethanol"), "C(C)O");',
    'assert.equal(smilesToName("CCO").text, "ethanol");',
  ].join("\n");
  const result = Bun.spawnSync(["node", "--experimental-wasm-modules", "--input-type=module", "-e", program], {
    stdout: "pipe", stderr: "pipe", timeout: 30_000,
  });
  expect(result.exitCode, result.stderr.toString()).toBe(0);
}, 30_000);

test("parse errors expose stable lowercase kind and engine message", () => {
  try {
    nameToMolecule("not a chemical name");
    throw new Error("Expected name parsing to fail");
  } catch (error: unknown) {
    expect(error).toBeInstanceOf(NazuchemError);
    if (!(error instanceof NazuchemError)) throw error;
    expect(error.name).toBe("NazuchemError");
    expect(error.kind).toBe("parse");
    expect(error.message.length).toBeGreaterThan(0);
  }
});

test("malformed IR and naming admission errors have distinct kinds", () => {
  for (const value of ["{", '{"nazuchem_ir":2,"atoms":[],"bonds":[]}', '{"nazuchem_ir":1,"atoms":[],"bonds":[],"surprise":true}', '{"nazuchem_ir":1,"atoms":[{"element":119,"hydrogens":0}],"bonds":[]}']) {
    try { moleculeToName(value); throw new Error("Expected IR rejection"); }
    catch (error: unknown) {
      expect(error).toBeInstanceOf(NazuchemError);
      if (error instanceof NazuchemError) expect(error.kind).toBe("serialize");
    }
  }
  expect(() => smilesToName("not a smiles")).toThrow(NazuchemError);
  try { smilesToName("not a smiles"); }
  catch (error: unknown) { if (error instanceof NazuchemError) expect(error.kind).toBe("naming"); }
});

test("tetrahedral, double-bond and atropisomer stereo survive IR naming", () => {
  const fixtures: [string, Stereo][] = [
    ["(2S)-2-hydroxypropanoic acid", { kind: "tetrahedral", center: 1, neighbors: [0, 6, 5, 2], winding: "clockwise" }],
    ["(2E)-but-2-ene", { kind: "double_bond", bond: 1, refs: [0, 3], relation: "trans" }],
    ["(1M)-[1,1'-binaphthalene]-2,2'-diol", { kind: "atropisomer", bond: 22, refs: [12, 20, 1, 9], winding: "anticlockwise" }],
  ];
  for (const [input, expectedStereo] of fixtures) {
    const molecule = nameToMolecule(input);
    assertValid(molecule);
    expect(molecule.stereo).toEqual([expectedStereo]);
    const before = JSON.stringify(molecule);
    const result = moleculeToName(molecule);
    expect(result.text).toBe(input);
    expect(JSON.stringify(molecule)).toBe(before);
    const roundtrip = nameToMolecule(result.text);
    assertValid(roundtrip);
    expect(roundtrip.stereo).toEqual(molecule.stereo);
  }
});

test("unknown stereo retains omission counts and independently assessed alternatives", () => {
  const molecule = nameToMolecule("(2S)-2-hydroxypropanoic acid");
  const stereo = molecule.stereo?.find((entry) => entry.kind === "tetrahedral");
  if (!stereo || stereo.kind !== "tetrahedral") throw new Error("Missing tetrahedral fixture");
  stereo.winding = "unknown";
  assertValid(molecule);
  const result = moleculeToName(molecule);
  expect(result.status).toBe("not_preferred");
  expect(result.reason).toBe("unknown_configuration");
  expect(result.omittedUnits).toBe(1);
  expect(result.omittedStereo).toEqual(["2"]);
  expect(result.alternatives.length).toBeGreaterThan(0);
  expect(result.alternatives[0]).toEqual({ name: "(2Ξ)-2-hydroxypropanoic acid", status: "preferred", reason: null });
});

test("radical-ion verification preserves the recovered neutral target", () => {
  const result: Name = smilesToName("C[N+](C)C");
  expect(result.status).toBe("not_preferred");
  expect(result.reason).toBe("radical_ion");
  expect(result.verification).not.toBeNull();
  expect(result.verification?.name).toBe("N,N-dimethylmethanamine");
  expect(result.verification?.smiles.length).toBeGreaterThan(0);
});

test("extreme spiro bridge lengths fail without trapping either name parser", () => {
  const inputs = ["2147483647", "4294967295", "18446744073709551615", "9".repeat(39)].flatMap((number) => [
    `spiro[${number}.5]pentane`, `spiro[${number}.${number}.5]pentane`,
    `bicyclo[${number}.2.1]heptane`, `[${number}]annulene`,
    `${number}methane`, `${number}-methylpentane`, `[1,${number}'-biphenyl]`,
    `methanol compound with ethanol (1:${number})`, `lambda^${number}-phosphane`,
  ]);
  for (const input of inputs) {
    for (const parse of [nameToSmiles, nameToMolecule]) {
      try {
        parse(input);
        throw new Error(`Expected parser rejection for ${input}`);
      } catch (error: unknown) {
        expect(error).toBeInstanceOf(NazuchemError);
        if (error instanceof NazuchemError) expect(error.kind).toBe("parse");
      }
      expect(nameToSmiles("ethanol")).toBe("C(C)O");
    }
  }
});

test("object input rejects nonfinite numbers before nullable fields can discard them", () => {
  for (const field of ["isotope", "atom_map", "charge"] as const) {
    for (const value of [Infinity, -Infinity, NaN]) {
      const molecule: Molecule = { nazuchem_ir: 1, atoms: [{ element: 6, hydrogens: 4, [field]: value }], bonds: [] };
      try {
        moleculeToName(molecule);
        throw new Error(`Expected nonfinite ${field} to fail`);
      } catch (error: unknown) {
        expect(error).toBeInstanceOf(NazuchemError);
        if (error instanceof NazuchemError) expect(error.kind).toBe("serialize");
      }
    }
  }
  expect(moleculeToName('{"nazuchem_ir":1,"atoms":[{"element":6,"hydrogens":4,"isotope":null}],"bonds":[]}').text).toBe("methane");
});

test("a WebAssembly trap maps to internal and subsequent naming still works", () => {
  const raw = packageRequire("./dist/node/nazuchem.js") as { nameToSmiles: (name: string) => string };
  const original = raw.nameToSmiles;
  // A real unreachable instruction exercises RuntimeError handling without a
  // production panic hook or a test-only export in the chemical engine.
  const trapModule = new WebAssembly.Module(new Uint8Array([
    0, 97, 115, 109, 1, 0, 0, 0,
    1, 4, 1, 96, 0, 0,
    3, 2, 1, 0,
    7, 8, 1, 4, 116, 114, 97, 112, 0, 0,
    10, 5, 1, 3, 0, 0, 11,
  ]));
  const trap = new WebAssembly.Instance(trapModule).exports.trap as () => never;
  raw.nameToSmiles = () => trap();
  try {
    expect(() => nameToSmiles("ethanol")).toThrow(NazuchemError);
    try { nameToSmiles("ethanol"); }
    catch (error: unknown) {
      expect(error).toBeInstanceOf(NazuchemError);
      if (error instanceof NazuchemError) expect(error.kind).toBe("internal");
    }
  } finally {
    raw.nameToSmiles = original;
  }
  expect(nameToSmiles("ethanol")).toBe("C(C)O");
  expect(smilesToName("CCO").text).toBe("ethanol");
});
