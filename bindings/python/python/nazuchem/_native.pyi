# SPDX-License-Identifier: MIT OR Apache-2.0
# Origin: original
# Source: Nazuchem | bindings/python/python/nazuchem/_native.pyi | Phase B step 9

from . import ErrorKind

IR_VERSION: int
JSON_SCHEMA: str

class NazuchemError(ValueError):
    kind: ErrorKind

def name_to_smiles(name: str) -> str: ...
def smiles_to_name(smiles: str) -> str: ...
def name_to_molecule(name: str) -> str: ...
def molecule_to_name(molecule_json: str) -> str: ...
