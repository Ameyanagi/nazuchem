// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | bindings/python/src/lib.rs | Phase B step 9

//! The Python facade owns ergonomic types; this module transports lossless JSON.

use pyo3::create_exception;
use pyo3::exceptions::PyException;
use pyo3::prelude::*;

create_exception!(nazuchem, NazuchemError, PyException);

fn binding_error(py: Python<'_>, kind: &str, message: String) -> PyErr {
    let error = NazuchemError::new_err(message);
    // The exception instance has a dict, so this attribute is available to callers.
    if let Err(attribute_error) = error.value(py).setattr("kind", kind) {
        return attribute_error;
    }
    error
}

fn engine_error(py: Python<'_>, error: nazuchem::Error) -> PyErr {
    binding_error(py, &format!("{:?}", error.kind()), error.to_string())
}

fn name_json(py: Python<'_>, name: &nazuchem::s2n::Name) -> PyResult<String> {
    serde_json::to_string(name).map_err(|error| binding_error(py, "Serialize", error.to_string()))
}

#[pyfunction]
fn name_to_smiles(py: Python<'_>, name: &str) -> PyResult<String> {
    py.detach(|| nazuchem::name_to_smiles(name))
        .map_err(|error| engine_error(py, error))
}

#[pyfunction]
fn smiles_to_name(py: Python<'_>, smiles: &str) -> PyResult<String> {
    let name = py
        .detach(|| nazuchem::smiles_to_name(smiles))
        .map_err(|error| engine_error(py, error))?;
    name_json(py, &name)
}

#[pyfunction]
fn name_to_molecule(py: Python<'_>, name: &str) -> PyResult<String> {
    let molecule = py
        .detach(|| nazuchem::name_to_molecule(name))
        .map_err(|error| engine_error(py, error))?;
    molecule
        .to_json()
        .map_err(|error| binding_error(py, "Serialize", error.to_string()))
}

#[pyfunction]
fn molecule_to_name(py: Python<'_>, molecule_json: &str) -> PyResult<String> {
    let molecule = nazuchem::Molecule::from_json(molecule_json)
        .map_err(|error| binding_error(py, "Serialize", error.to_string()))?;
    let name = py
        .detach(|| nazuchem::molecule_to_name(&molecule))
        .map_err(|error| engine_error(py, error))?;
    name_json(py, &name)
}

#[pymodule]
fn _native(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("NazuchemError", module.py().get_type::<NazuchemError>())?;
    module.add("IR_VERSION", nazuchem::IR_VERSION)?;
    module.add("JSON_SCHEMA", nazuchem::JSON_SCHEMA)?;
    module.add_function(wrap_pyfunction!(name_to_smiles, module)?)?;
    module.add_function(wrap_pyfunction!(smiles_to_name, module)?)?;
    module.add_function(wrap_pyfunction!(name_to_molecule, module)?)?;
    module.add_function(wrap_pyfunction!(molecule_to_name, module)?)?;
    Ok(())
}
