// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | bindings/wasm/src/wrapper.js | Phase B step 9

export class NazuchemError extends Error {
  constructor(message, kind) {
    super(message);
    this.name = "NazuchemError";
    this.kind = kind;
  }
}

// The web target initializes Wasm after importing the module. Keep the constants
// in JS so they are available before initialization, while the build derives
// their values from the same Rust schema bundled in the release.
export const IR_VERSION = __NAZUCHEM_IR_VERSION__;
export const JSON_SCHEMA = __NAZUCHEM_SCHEMA__;

function call(operation, value) {
  try {
    return operation(value);
  } catch (source) {
    if (source instanceof Error && source.name === "NazuchemError") {
      throw new NazuchemError(source.message, source.kind);
    }
    if (source instanceof WebAssembly.RuntimeError) {
      throw new NazuchemError(source.message, "internal");
    }
    throw source;
  }
}

export function nameToSmiles(name) { return call(raw.nameToSmiles, name); }
function flattenName(result) {
  return {
    text: result.text,
    omittedStereo: result.omitted_stereo,
    omittedUnits: result.omitted_units,
    verification: result.verification,
    status: result.status.status,
    reason: result.status.reason ?? null,
    alternatives: result.alternatives.map((alternative) => ({
      name: alternative.name,
      status: alternative.status.status,
      reason: alternative.status.reason ?? null,
    })),
  };
}
export function smilesToName(smiles) { return flattenName(call(raw.smilesToName, smiles)); }
export function nameToMolecule(name) { return call(raw.nameToMolecule, name); }
// JSON.stringify would silently turn nonfinite numbers into null. Validate the
// supplied object first so nullable IR fields cannot discard invalid input.
function assertFiniteNumbers(value, seen = new WeakSet()) {
  if (typeof value === "number" && !Number.isFinite(value)) {
    throw new NazuchemError("molecule numbers must be finite", "serialize");
  }
  if (value === null || typeof value !== "object" || seen.has(value)) return;
  seen.add(value);
  for (const child of Object.values(value)) assertFiniteNumbers(child, seen);
}
export function moleculeToName(molecule) {
  assertFiniteNumbers(molecule);
  return flattenName(call(raw.moleculeToName, molecule));
}
