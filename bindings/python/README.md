<!-- SPDX-License-Identifier: MIT OR Apache-2.0 -->
<!-- Origin: original -->
<!-- Source: Nazuchem | bindings/python/README.md | Phase B step 9 -->

# Nazuchem for Python

Local Python bindings for chemical names, semantic CXSMILES, and the versioned
Nazuchem molecule IR. Python 3.10 or newer is required. The PyO3 extension uses
`abi3-py310`, with a thin typed Python facade and no runtime Python dependencies.
Distribution remains blocked by the repository's unresolved provenance decisions;
the Rust binding crate has `publish = false`. See the generated `NOTICE` and the
root `NOTICE` for the established licence scope and unresolved materials.

```python
import json
import nazuchem

smiles = nazuchem.name_to_smiles("ethanol")
name = nazuchem.smiles_to_name(smiles)
assert name.text == "ethanol"

molecule = nazuchem.name_to_molecule("(2R)-2-hydroxypropanoic acid")
assert molecule["nazuchem_ir"] == nazuchem.IR_VERSION
roundtrip = nazuchem.molecule_to_name(molecule)
from_json = nazuchem.molecule_to_name(json.dumps(molecule))
schema = json.loads(nazuchem.JSON_SCHEMA)
```

The four functions are `name_to_smiles(str) -> str`, `smiles_to_name(str) -> Name`,
`name_to_molecule(str) -> dict`, and `molecule_to_name(dict | str) -> Name`.
The IR object includes its version and all supported structure metadata. Rust
validates every IR input against the versioned format and graph invariants;
unknown fields, unsupported versions, and invalid atom/bond/stereo references
fail explicitly. The Python facade rejects non-JSON values, non-string object
keys, and non-finite numbers in dict inputs. It never changes the input dict.

`Name` is a frozen dataclass containing `text`, `status`, `reason`, `alternatives`,
`omitted_stereo`, `omitted_units`, and `verification`. Status and reason strings
retain Rust's snake_case JSON values. `alternatives` contains frozen
`AlternativeName(name, status, reason)` objects. `verification`, when present,
is a frozen `Verification(name, smiles)` for a recovered neutral radical-ion
target. Collection fields are tuples. No Rust result fields are discarded.

Engine failures raise `NazuchemError`. Its message is the original Rust message,
and `.kind` is the Rust `ErrorKind` variant name: `Init`, `Parse`, `Warning`,
`Serialize`, or `Naming`. IR decode/validation failures use `Serialize`.
Python argument type errors raise `TypeError`; non-finite numbers and circular
dict/list references raise `ValueError`.

From this directory, build and test the actual wheel with uv:

```sh
uv sync --no-install-project --python 3.12
uv run --no-sync python scripts/build_wheel.py
uv run --no-project --python 3.12 --with ./dist/<wheel>.whl --with pytest --with jsonschema python -m pytest tests
uv run --no-sync ruff check .
uv run --no-sync ruff format --check .
uv run --no-sync ty check
```

`build_wheel.py` generates all four legal files using the root collector, builds
with `uv` and `maturin --locked`, then unzips the wheel to byte-compare each
PEP 639 `.dist-info/licenses/` member against freshly generated files. It fails
on missing, stale, or changed notice bytes. The generated files, virtualenv,
native libraries, and wheels are ignored by Git. Run the checker separately with
`uv run --no-sync python scripts/check_wheel_notices.py dist/<wheel>.whl`.
Add `--release` for optimization, or `--target <Rust triple>` for a target build;
the wrapper forwards that target to the collector and wheel checker as well.
For manual maturin builds, run the collector first and run the checker afterward.
For plain Cargo checks, set `PYO3_PYTHON` to a Python 3.10+ interpreter if the
system default Python is older. These commands build local artifacts only.
