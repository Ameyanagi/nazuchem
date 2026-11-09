# SPDX-License-Identifier: MIT OR Apache-2.0
# Origin: original
# Source: Nazuchem | bindings/python/python/nazuchem/__init__.py | Phase B step 9

"""Chemical naming with lossless, versioned JSON structure documents."""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Literal, TypeAlias, TypedDict, cast

from . import _native
from ._native import IR_VERSION, JSON_SCHEMA, NazuchemError

__all__ = [
    "IR_VERSION",
    "JSON_SCHEMA",
    "AlternativeName",
    "ErrorKind",
    "JsonValue",
    "Molecule",
    "Name",
    "NameStatus",
    "NazuchemError",
    "Verification",
    "molecule_to_name",
    "name_to_molecule",
    "name_to_smiles",
    "smiles_to_name",
]

ErrorKind: TypeAlias = Literal["init", "parse", "warning", "serialize", "naming", "internal"]
NameStatus: TypeAlias = Literal["preferred", "not_preferred", "unverified"]
JsonValue: TypeAlias = bool | int | float | str | None | list["JsonValue"] | dict[str, "JsonValue"]
Molecule: TypeAlias = dict[str, JsonValue]


@dataclass(frozen=True)
class AlternativeName:
    """An alternative name and its independent preference assessment."""

    name: str
    status: NameStatus
    reason: str | None


@dataclass(frozen=True)
class Verification:
    """The recovered neutral target of a descriptive radical-ion name."""

    name: str
    smiles: str


@dataclass(frozen=True)
class Name:
    """A generated name, retaining all of the Rust engine's result metadata.

    Status and reason use the Rust JSON vocabulary, such as ``preferred`` and
    ``unaudited_path``. Collection fields are tuples to keep the result immutable.
    """

    text: str
    status: NameStatus
    reason: str | None
    alternatives: tuple[AlternativeName, ...]
    omitted_stereo: tuple[str, ...]
    omitted_units: int
    verification: Verification | None


class _Status(TypedDict, total=False):
    status: NameStatus
    reason: str


class _Alternative(TypedDict):
    name: str
    status: _Status


class _Verification(TypedDict):
    name: str
    smiles: str


class _Name(TypedDict):
    text: str
    status: _Status
    alternatives: list[_Alternative]
    omitted_stereo: list[str]
    omitted_units: int
    verification: _Verification | None


def _name_from_json(document: str) -> Name:
    result = cast(_Name, json.loads(document))
    status = result["status"]
    verification = result["verification"]
    return Name(
        text=result["text"],
        status=status["status"],
        reason=status.get("reason"),
        alternatives=tuple(
            AlternativeName(item["name"], item["status"]["status"], item["status"].get("reason"))
            for item in result["alternatives"]
        ),
        omitted_stereo=tuple(result["omitted_stereo"]),
        omitted_units=result["omitted_units"],
        verification=(
            Verification(verification["name"], verification["smiles"])
            if verification is not None
            else None
        ),
    )


def name_to_smiles(name: str) -> str:
    """Parse a name strictly and return semantic CXSMILES."""
    return _native.name_to_smiles(name)


def smiles_to_name(smiles: str) -> Name:
    """Generate a name and its preference, alternatives, and omission metadata."""
    return _name_from_json(_native.smiles_to_name(smiles))


def name_to_molecule(name: str) -> Molecule:
    """Parse a name strictly and return a validated IR v1 JSON object."""
    return cast(Molecule, json.loads(_native.name_to_molecule(name)))


def molecule_to_name(molecule: Molecule | str) -> Name:
    """Name an IR dict or JSON string without modifying the caller's object.

    Python values outside JSON raise TypeError or ValueError. Invalid IR raises
    NazuchemError with kind ``serialize``; naming failures have kind ``naming``.
    """
    if isinstance(molecule, str):
        return _name_from_json(_native.molecule_to_name(molecule))
    if not isinstance(molecule, dict):
        raise TypeError("molecule must be an IR JSON object (dict) or JSON string")
    document = json.dumps(molecule, allow_nan=False, separators=(",", ":"))
    return _name_from_json(_native.molecule_to_name(document))
