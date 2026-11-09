// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | bindings/wasm/tests/bindings.test.ts | Phase B step 9

import { expect, test } from "bun:test";
import Ajv2020 from "ajv/dist/2020.js";
import { createRequire } from "node:module";
import { readFile } from "node:fs/promises";
import {
  IR_VERSION, JSON_SCHEMA, NazuchemError, moleculeToName,
  nameToMolecule, nameToSmiles, smilesToName,
} from "../dist/node/index.cjs";
import type { Molecule, Name, Stereo } from "../src/index.d.ts";

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
  expect(name).toEqual({ text: "ethanol", status: "preferred", reason: null, alternatives: [], omitted_stereo: [], omitted_units: 0, verification: null });
  const molecule = nameToMolecule("ethanol");
  assertValid(molecule);
  expect(molecule.atoms).toHaveLength(3);
  expect(moleculeToName(molecule)).toEqual(name);
  expect(moleculeToName(JSON.stringify(molecule))).toEqual(name);
  expect(smilesToName("CCCCO").text).toBe("butan-1-ol");
});

test("IR constants match Rust and the checked-in schema byte for byte", async () => {
  expect(IR_VERSION).toBe(1);
  expect(JSON_SCHEMA).toBe(await readFile(new URL("../../../crates/nazuchem-ir/schema/nazuchem-ir-1.schema.json", import.meta.url), "utf8"));
  const raw: { irVersion(): number; jsonSchema(): string } = createRequire(import.meta.url)("../dist/node/nazuchem.js");
  expect(raw.irVersion()).toBe(IR_VERSION);
  expect(raw.jsonSchema()).toBe(JSON_SCHEMA);
});

test("browser default initialization fetches the shared Wasm and runs the naming API", async () => {
  const sharedWasm = new URL("../dist/nazuchem_bg.wasm", import.meta.url);
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
    const web = await import("../dist/web/index.js");
    await web.default();
    expect(fetchedUrls).toEqual([sharedWasm.href]);
    expect(web.nameToSmiles("ethanol")).toBe("C(C)O");
    expect(web.nameToMolecule("ethanol").atoms).toHaveLength(3);
  } finally {
    globalThis.fetch = originalFetch;
  }
});

test("bundler static Wasm imports run both naming directions in a fresh Node process", () => {
  const entry = new URL("../dist/bundler/index.js", import.meta.url).href;
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

test("parse errors preserve facade variant and message", () => {
  try {
    nameToMolecule("not a chemical name");
    throw new Error("Expected name parsing to fail");
  } catch (error: unknown) {
    expect(error).toBeInstanceOf(NazuchemError);
    if (!(error instanceof NazuchemError)) throw error;
    expect(error.name).toBe("NazuchemError");
    expect(error.kind).toBe("Parse");
    expect(error.message.length).toBeGreaterThan(0);
  }
});

test("malformed IR and naming admission errors have distinct kinds", () => {
  for (const value of ["{", '{"nazuchem_ir":2,"atoms":[],"bonds":[]}', '{"nazuchem_ir":1,"atoms":[],"bonds":[],"surprise":true}', '{"nazuchem_ir":1,"atoms":[{"element":119,"hydrogens":0}],"bonds":[]}']) {
    try { moleculeToName(value); throw new Error("Expected IR rejection"); }
    catch (error: unknown) {
      expect(error).toBeInstanceOf(NazuchemError);
      if (error instanceof NazuchemError) expect(error.kind).toBe("Serialize");
    }
  }
  expect(() => smilesToName("not a smiles")).toThrow(NazuchemError);
  try { smilesToName("not a smiles"); }
  catch (error: unknown) { if (error instanceof NazuchemError) expect(error.kind).toBe("Naming"); }
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
  expect(result.omitted_units).toBe(1);
  expect(result.omitted_stereo).toEqual(["2"]);
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
