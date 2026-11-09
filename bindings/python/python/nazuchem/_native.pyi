# SPDX-License-Identifier: MIT OR Apache-2.0
# Origin: original
# Source: Nazuchem | bindings/python/python/nazuchem/_native.pyi | Phase B step 9

IR_VERSION: int
JSON_SCHEMA: str

class NazuchemError(Exception):
    kind: str

def name_to_smiles(name: str) -> str: ...
def smiles_to_name(smiles: str) -> str: ...
def name_to_molecule(name: str) -> str: ...
def molecule_to_name(molecule_json: str) -> str: ...
