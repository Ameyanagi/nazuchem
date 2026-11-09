// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | bindings/wasm/scripts/check-types.ts | Phase B step 9

import { mkdir, readFile, writeFile } from "node:fs/promises";
import { resolve } from "node:path";

interface Schema {
  $ref?: string;
  const?: string | number;
  enum?: (string | number)[];
  type?: string | string[];
  anyOf?: Schema[];
  oneOf?: Schema[];
  properties?: Record<string, Schema>;
  required?: string[];
  items?: Schema;
  minItems?: number;
  maxItems?: number;
  $defs?: Record<string, Schema>;
}

const root = resolve(import.meta.dir, "../../..");
const packageRoot = resolve(import.meta.dir, "..");
const schema: Schema = JSON.parse(await readFile(resolve(root, "crates/nazuchem-ir/schema/nazuchem-ir-1.schema.json"), "utf8"));

// Derive the entire structural TypeScript model from the published schema. The
// compiler proves exact equality, including optionality, every union branch,
// nullable field and tuple length. JSON numeric/uniqueness constraints are checked
// separately by Ajv because TypeScript numbers cannot express those predicates.
function typeOf(node: Schema): string {
  if (node.$ref) {
    if (!node.$ref.startsWith("#/$defs/")) throw new Error(`Unsupported reference ${node.$ref}`);
    return node.$ref.slice("#/$defs/".length);
  }
  if (node.const !== undefined) return JSON.stringify(node.const);
  if (node.enum) return node.enum.map((value) => JSON.stringify(value)).join(" | ");
  const union = node.oneOf ?? node.anyOf;
  if (union) return union.map(typeOf).join(" | ");
  if (Array.isArray(node.type)) return node.type.map((type) => typeOf({ ...node, type })).join(" | ");
  switch (node.type) {
    case "null": return "null";
    case "integer":
    case "number": return "number";
    case "string": return "string";
    case "boolean": return "boolean";
    case "object": {
      if (!node.properties) throw new Error("Object without declared properties");
      return `{ ${Object.entries(node.properties).map(([key, value]) => `${JSON.stringify(key)}${node.required?.includes(key) ? "" : "?"}: ${typeOf(value)}`).join("; ")} }`;
    }
    case "array": {
      if (!node.items) throw new Error("Array without items");
      const item = typeOf(node.items);
      const min = node.minItems ?? 0;
      if (node.maxItems === min) return `[${Array.from({ length: min }, () => item).join(", ")}]`;
      if (node.maxItems !== undefined) throw new Error("Unsupported bounded variable-length array");
      if (min > 0) return `[${[...Array.from({ length: min }, () => item), `...(${item})[]`].join(", ")}]`;
      return `(${item})[]`;
    }
    default: throw new Error(`Unsupported schema shape: ${JSON.stringify(node)}`);
  }
}

const definitions = schema.$defs;
if (!definitions) throw new Error("Missing schema definitions");
const header = "// SPDX-License-Identifier: MIT OR Apache-2.0\n// Origin: original\n// Source: Nazuchem | bindings/wasm/scripts/check-types.ts | Phase B step 9\n";
const names = [...Object.keys(definitions), "Molecule"];
const derived = [...Object.entries(definitions), ["Molecule", schema] as const]
  .map(([name, node]) => `export type ${name} = ${typeOf(node)};`).join("\n");
await mkdir(resolve(packageRoot, ".generated"), { recursive: true });
await writeFile(resolve(packageRoot, ".generated/schema-types.ts"), header + derived + "\n");
const checks = names.map((name) => `type Check${name} = Assert<Equal<Handwritten.${name}, Generated.${name}>>;`).join("\n");
await writeFile(resolve(packageRoot, ".generated/equivalence.ts"), header + `
import type * as Handwritten from "../src/index.d.ts";
import type * as Generated from "./schema-types.ts";
type Equal<A, B> = (<T>() => T extends A ? 1 : 2) extends (<T>() => T extends B ? 1 : 2) ? true : false;
type Assert<T extends true> = T;
${checks}
`);

// Naming reasons are not in the Molecule schema. Check their exact serde names
// against Rust so new engine reasons cannot silently disappear from declarations.
const rustVocabulary = await readFile(resolve(root, "crates/nazuchem-ir/src/lib.rs"), "utf8");
const reasons = rustVocabulary.match(/pub enum NameReason \{([^}]+)\}/)?.[1];
if (!reasons) throw new Error("Cannot locate Rust NameReason");
const reasonTypes = reasons.split(",").map((name) => name.trim()).filter(Boolean)
  .map((name) => JSON.stringify(name.replace(/[A-Z]/g, (letter, index) => `${index ? "_" : ""}${letter.toLowerCase()}`))).join(" | ");
await writeFile(resolve(packageRoot, ".generated/name-reasons.ts"), header + `
import type { NameReason } from "../src/index.d.ts";
type Equal<A, B> = (<T>() => T extends A ? 1 : 2) extends (<T>() => T extends B ? 1 : 2) ? true : false;
type Assert<T extends true> = T;
type CheckNameReason = Assert<Equal<NameReason, ${reasonTypes}>>;
`);

const checked = Bun.spawnSync(["bun", "x", "--no-install", "tsc", "--project", "tsconfig.json"], { cwd: packageRoot, stdout: "inherit", stderr: "inherit" });
if (checked.exitCode !== 0) process.exit(checked.exitCode);
console.log(`TypeScript checked; all ${names.length} IR types exactly match schema-derived types; NameReason matches Rust.`);
