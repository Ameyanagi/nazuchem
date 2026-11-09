# SPDX-License-Identifier: MIT OR Apache-2.0
# Origin: original
# Source: Nazuchem | bindings/python/tests/test_api.py | Phase B step 9

import copy
import json
from dataclasses import FrozenInstanceError, asdict, is_dataclass
from typing import cast

import jsonschema
import nazuchem
import pytest


def test_ethanol_roundtrip_and_frozen_result() -> None:
    smiles = nazuchem.name_to_smiles("ethanol")
    result = nazuchem.smiles_to_name(smiles)
    assert result.text == "ethanol"
    assert is_dataclass(result)
    assert result.status == "preferred"
    assert result.reason is None
    assert result.alternatives == ()
    assert result.omitted_stereo == ()
    assert result.omitted_units == 0
    assert result.verification is None
    assert set(asdict(result)) == {
        "text",
        "status",
        "reason",
        "alternatives",
        "omitted_stereo",
        "omitted_units",
        "verification",
    }
    with pytest.raises(FrozenInstanceError):
        field = "text"
        setattr(result, field, "changed")
    assert nazuchem.name_to_smiles(result.text) == smiles


def test_ir_schema_and_dict_and_json_roundtrip() -> None:
    molecule = nazuchem.name_to_molecule("ethanol")
    before = copy.deepcopy(molecule)
    assert type(nazuchem.IR_VERSION) is int
    assert isinstance(nazuchem.JSON_SCHEMA, str)
    assert molecule["nazuchem_ir"] == nazuchem.IR_VERSION
    jsonschema.validate(molecule, json.loads(nazuchem.JSON_SCHEMA))
    assert nazuchem.molecule_to_name(molecule).text == "ethanol"
    assert nazuchem.molecule_to_name(json.dumps(molecule)).text == "ethanol"
    assert molecule == before


@pytest.mark.parametrize("function", [nazuchem.name_to_smiles, nazuchem.name_to_molecule])
def test_invalid_name_preserves_error_kind(function) -> None:
    with pytest.raises(nazuchem.NazuchemError) as raised:
        function("not a chemical name")
    assert raised.value.kind == "Parse"
    assert str(raised.value)


def test_naming_error_kind() -> None:
    with pytest.raises(nazuchem.NazuchemError) as raised:
        nazuchem.smiles_to_name("C>>O")
    assert raised.value.kind == "Naming"
    assert str(raised.value)


def test_radical_ion_verification_is_preserved() -> None:
    result = nazuchem.smiles_to_name("[CH2+][CH2]")
    assert result.text == "ethene radical cation"
    assert (result.status, result.reason) == ("not_preferred", "radical_ion")
    assert result.verification == nazuchem.Verification(name="ethene", smiles="C=C")


def test_axial_alternative_keeps_its_own_status_and_reason() -> None:
    result = nazuchem.smiles_to_name("Oc1ccc2ccccc2[c:900004]1-[c:900004]1c(O)ccc2ccccc12")
    assert result.status == "preferred"
    assert result.reason is None
    assert len(result.alternatives) == 1
    alternative = result.alternatives[0]
    assert alternative.name == result.text.replace("1M", "1Ra")
    assert (alternative.status, alternative.reason) == ("not_preferred", "axial_alternative")


def test_unknown_stereo_omissions_are_preserved() -> None:
    result = nazuchem.smiles_to_name("O=C(O)[CH:900001](F)Cl")
    assert result.omitted_stereo == ("2",)
    assert result.omitted_units == 1


@pytest.mark.parametrize(
    ("name", "kind", "descriptor"),
    [
        ("(2R)-butan-2-ol", "tetrahedral", "2R"),
        ("(2R)-2-hydroxypropanoic acid", "tetrahedral", "2R"),
        ("(2E)-but-2-ene", "double_bond", "E"),
    ],
)
def test_stereo_survives_ir(name: str, kind: str, descriptor: str) -> None:
    molecule = nazuchem.name_to_molecule(name)
    stereo = cast(list[dict[str, nazuchem.JsonValue]], molecule["stereo"])
    assert any(unit["kind"] == kind for unit in stereo)
    generated = nazuchem.molecule_to_name(molecule)
    assert descriptor in generated.text
    reconstructed = nazuchem.name_to_molecule(generated.text)
    jsonschema.validate(reconstructed, json.loads(nazuchem.JSON_SCHEMA))
    assert reconstructed["stereo"] == molecule["stereo"]


@pytest.mark.parametrize(
    "document",
    [
        {},
        {"nazuchem_ir": 999},
        {"nazuchem_ir": 1, "atoms": [], "bonds": [], "unknown": True},
        "not json",
        "[]",
    ],
)
def test_invalid_ir_is_a_serialization_error(document) -> None:
    with pytest.raises(nazuchem.NazuchemError) as raised:
        nazuchem.molecule_to_name(document)
    assert raised.value.kind == "Serialize"
    assert str(raised.value)


def test_bad_bond_reference_is_rejected() -> None:
    molecule = nazuchem.name_to_molecule("ethanol")
    bonds = cast(list[dict[str, nazuchem.JsonValue]], molecule["bonds"])
    bonds[0]["atoms"] = [0, 999]
    with pytest.raises(nazuchem.NazuchemError) as raised:
        nazuchem.molecule_to_name(molecule)
    assert raised.value.kind == "Serialize"


@pytest.mark.parametrize(
    "function",
    [nazuchem.name_to_smiles, nazuchem.smiles_to_name, nazuchem.name_to_molecule],
)
def test_string_arguments_are_not_coerced(function) -> None:
    with pytest.raises(TypeError):
        function(42)


@pytest.mark.parametrize("document", [[], 42, None, {1: "bad key"}, {"atoms": ()}])
def test_python_only_ir_values_are_rejected(document) -> None:
    with pytest.raises(TypeError):
        nazuchem.molecule_to_name(document)


@pytest.mark.parametrize("value", [float("nan"), float("inf"), float("-inf")])
def test_nonfinite_json_values_are_rejected(value: float) -> None:
    with pytest.raises(ValueError):
        nazuchem.molecule_to_name({"nazuchem_ir": value})


def test_circular_json_is_rejected() -> None:
    molecule: nazuchem.Molecule = {}
    molecule["cycle"] = molecule
    with pytest.raises(ValueError, match="Circular reference"):
        nazuchem.molecule_to_name(molecule)
