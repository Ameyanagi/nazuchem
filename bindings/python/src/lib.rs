// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | bindings/python/src/lib.rs | Phase B step 9

//! The Python facade owns ergonomic types; this module transports lossless JSON.

use std::panic::{UnwindSafe, catch_unwind};

use pyo3::create_exception;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

create_exception!(nazuchem, NazuchemError, PyValueError);

struct BindingFailure {
    kind: &'static str,
    message: String,
}

fn binding_error(py: Python<'_>, failure: BindingFailure) -> PyErr {
    let error = NazuchemError::new_err(failure.message);
    // The exception instance has a dict, so this attribute is available to callers.
    if let Err(attribute_error) = error.value(py).setattr("kind", failure.kind) {
        return attribute_error;
    }
    error
}

fn engine_error(error: nazuchem::Error) -> BindingFailure {
    let kind = match error.kind() {
        nazuchem::ErrorKind::Init => "init",
        nazuchem::ErrorKind::Parse => "parse",
        nazuchem::ErrorKind::Warning => "warning",
        nazuchem::ErrorKind::Serialize => "serialize",
        nazuchem::ErrorKind::Naming => "naming",
        _ => "internal",
    };
    BindingFailure {
        kind,
        message: error.to_string(),
    }
}

fn serialization_error(error: impl std::fmt::Display) -> BindingFailure {
    BindingFailure {
        kind: "serialize",
        message: error.to_string(),
    }
}

fn run_engine<T, F>(py: Python<'_>, operation: F) -> PyResult<T>
where
    T: Send,
    F: FnOnce() -> Result<T, BindingFailure> + Send + UnwindSafe,
{
    let result = py
        .detach(|| catch_unwind(operation))
        .unwrap_or_else(|panic| {
            let detail = panic
                .downcast_ref::<String>()
                .map(String::as_str)
                .or_else(|| panic.downcast_ref::<&str>().copied())
                .unwrap_or("unknown panic payload");
            Err(BindingFailure {
                kind: "internal",
                message: format!("chemical engine panicked: {detail}"),
            })
        });
    result.map_err(|failure| binding_error(py, failure))
}

#[pyfunction]
fn name_to_smiles(py: Python<'_>, name: &str) -> PyResult<String> {
    run_engine(py, || nazuchem::name_to_smiles(name).map_err(engine_error))
}

#[pyfunction]
fn smiles_to_name(py: Python<'_>, smiles: &str) -> PyResult<String> {
    run_engine(py, || {
        let name = nazuchem::smiles_to_name(smiles).map_err(engine_error)?;
        serde_json::to_string(&name).map_err(serialization_error)
    })
}

#[pyfunction]
fn name_to_molecule(py: Python<'_>, name: &str) -> PyResult<String> {
    run_engine(py, || {
        let molecule = nazuchem::name_to_molecule(name).map_err(engine_error)?;
        molecule.to_json().map_err(serialization_error)
    })
}

#[pyfunction]
fn molecule_to_name(py: Python<'_>, molecule_json: &str) -> PyResult<String> {
    run_engine(py, || {
        let molecule = nazuchem::Molecule::from_json(molecule_json).map_err(serialization_error)?;
        let name = nazuchem::molecule_to_name(&molecule).map_err(engine_error)?;
        serde_json::to_string(&name).map_err(serialization_error)
    })
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn panic_is_an_internal_value_error_and_interpreter_stays_usable() {
        Python::initialize();
        Python::attach(|py| {
            let failure = run_engine::<(), _>(py, || panic!("induced engine panic")).unwrap_err();
            assert!(failure.is_instance_of::<NazuchemError>(py));
            assert!(failure.is_instance_of::<PyValueError>(py));
            assert_eq!(
                failure
                    .value(py)
                    .getattr("kind")
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                "internal"
            );
            assert!(failure.to_string().contains("induced engine panic"));
            assert_eq!(name_to_smiles(py, "ethanol").unwrap(), "C(C)O");
        });
    }
}
