# SPDX-License-Identifier: MIT OR Apache-2.0
# Origin: original
# Source: Nazuchem | bindings/python/python/nazuchem/__init__.pyi | Phase B step 9

from dataclasses import dataclass
from typing import TypeAlias

from ._native import IR_VERSION as IR_VERSION
from ._native import JSON_SCHEMA as JSON_SCHEMA
from ._native import NazuchemError as NazuchemError

JsonValue: TypeAlias = bool | int | float | str | None | list[JsonValue] | dict[str, JsonValue]
Molecule: TypeAlias = dict[str, JsonValue]

@dataclass(frozen=True)
class AlternativeName:
    name: str
    status: str
    reason: str | None

@dataclass(frozen=True)
class Verification:
    name: str
    smiles: str

@dataclass(frozen=True)
class Name:
    text: str
    status: str
    reason: str | None
    alternatives: tuple[AlternativeName, ...]
    omitted_stereo: tuple[str, ...]
    omitted_units: int
    verification: Verification | None

def name_to_smiles(name: str) -> str: ...
def smiles_to_name(smiles: str) -> Name: ...
def name_to_molecule(name: str) -> Molecule: ...
def molecule_to_name(molecule: Molecule | str) -> Name: ...
