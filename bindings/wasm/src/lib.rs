// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | bindings/wasm/src/lib.rs | Phase B step 9

//! Private WebAssembly bridge. Public JavaScript declarations and errors live in
//! the package wrapper; all chemical data retains the facade's JSON shape.
#![forbid(unsafe_code)]

#[cfg(target_arch = "wasm32")]
mod bindings {
    use serde::Serialize;
    use wasm_bindgen::prelude::*;

    fn error(kind: &str, message: &str) -> JsValue {
        let error = js_sys::Error::new(message);
        error.set_name("NazuchemError");
        // A fresh Error is extensible, so defining this own property succeeds.
        let _ = js_sys::Reflect::set(&error, &"kind".into(), &kind.into());
        error.into()
    }

    fn facade_error(source: nazuchem::Error) -> JsValue {
        error(&format!("{:?}", source.kind()), &source.to_string())
    }

    fn to_object(value: &impl Serialize) -> Result<JsValue, JsValue> {
        let json = serde_json::to_string(value)
            .map_err(|source| error("Serialize", &source.to_string()))?;
        js_sys::JSON::parse(&json).map_err(|_| error("Serialize", "invalid serialized JSON"))
    }

    /// Interpret a name using the facade's strict name-to-structure policy.
    #[wasm_bindgen(js_name = nameToSmiles)]
    pub fn name_to_smiles(name: &str) -> Result<String, JsValue> {
        nazuchem::name_to_smiles(name).map_err(facade_error)
    }

    /// Return the complete Rust Name, including nested preference assessment.
    #[wasm_bindgen(js_name = smilesToName)]
    pub fn smiles_to_name(smiles: &str) -> Result<JsValue, JsValue> {
        to_object(&nazuchem::smiles_to_name(smiles).map_err(facade_error)?)
    }

    /// Return a validated IR v1 plain object with the Rust JSON field names.
    #[wasm_bindgen(js_name = nameToMolecule)]
    pub fn name_to_molecule(name: &str) -> Result<JsValue, JsValue> {
        to_object(&nazuchem::name_to_molecule(name).map_err(facade_error)?)
    }

    /// Accept the same IR object or JSON string as the facade's JSON reader.
    #[wasm_bindgen(js_name = moleculeToName)]
    pub fn molecule_to_name(input: JsValue) -> Result<JsValue, JsValue> {
        let json = if let Some(json) = input.as_string() {
            json
        } else {
            js_sys::JSON::stringify(&input)
                .ok()
                .and_then(|value| value.as_string())
                .ok_or_else(|| {
                    error("Serialize", "molecule must be a JSON object or JSON string")
                })?
        };
        let molecule = nazuchem::Molecule::from_json(&json)
            .map_err(|source| error("Serialize", &source.to_string()))?;
        to_object(&nazuchem::molecule_to_name(&molecule).map_err(facade_error)?)
    }

    /// The wrapper exposes this Rust constant as the JavaScript JSON_SCHEMA.
    #[wasm_bindgen(js_name = jsonSchema)]
    pub fn json_schema() -> String {
        nazuchem::JSON_SCHEMA.to_owned()
    }

    /// The wrapper exposes this Rust constant as the JavaScript IR_VERSION.
    #[wasm_bindgen(js_name = irVersion)]
    pub fn ir_version() -> u32 {
        // IR v1 uses a small integer; avoid a JavaScript BigInt for this constant.
        u32::try_from(nazuchem::IR_VERSION).expect("supported IR version fits a JavaScript number")
    }
}
