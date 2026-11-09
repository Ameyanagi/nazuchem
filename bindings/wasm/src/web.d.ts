// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | bindings/wasm/src/web.d.ts | Phase B step 9

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;
export interface InitOutput {
  readonly memory: WebAssembly.Memory;
}
export function initSync(module: { module: BufferSource | WebAssembly.Module }): InitOutput;
export default function init(input?: { module_or_path: InitInput | Promise<InitInput> }): Promise<InitOutput>;
