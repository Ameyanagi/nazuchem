# SPDX-License-Identifier: MIT OR Apache-2.0
# Origin: original
# Source: Nazuchem | scripts/check_release.py

"""Check release versions/gates, or package the five libraries in a scratch copy."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory

import tomllib

ROOT = Path(__file__).resolve().parents[1]
LIBRARIES = (
    "nazuchem-ir",
    "nazuchem-kernel",
    "nazuchem-n2s",
    "nazuchem-s2n",
    "nazuchem",
)
BINDINGS = ("nazuchem-python", "nazuchem-wasm")


def check(root: Path, *, dry_run: bool, tag: str | None, output: Path | None) -> None:
    workspace = tomllib.loads((root / "Cargo.toml").read_text())["workspace"]
    versions = {"workspace.package": workspace["package"]["version"]}
    packages = {}
    for member in workspace["members"]:
        manifest = root / member / "Cargo.toml"
        package = tomllib.loads(manifest.read_text())["package"]
        name = package["name"]
        if name in packages:
            raise ValueError(f"duplicate workspace package: {name}")
        packages[name] = package
        version = package["version"]
        versions[str(manifest.relative_to(root))] = (
            workspace["package"]["version"] if version == {"workspace": True} else version
        )
    if set(packages) != set(LIBRARIES + BINDINGS):
        raise ValueError(f"update release package inventory: {sorted(packages)}")
    versions["bindings/python/pyproject.toml"] = tomllib.loads(
        (root / "bindings/python/pyproject.toml").read_text()
    )["project"]["version"]
    npm = json.loads((root / "bindings/wasm/package.json").read_text())
    versions["bindings/wasm/package.json"] = npm["version"]
    version = versions["workspace.package"]
    if not isinstance(version, str) or not re.fullmatch(
        r"\d+\.\d+\.\d+(?:-[0-9A-Za-z.-]+)?", version
    ):
        raise ValueError(f"invalid release version: {version!r}")
    for path, actual in versions.items():
        if actual != version:
            raise ValueError(f"{path}: version {actual!r} differs from {version}")
    expected_tag = f"v{version}"
    if tag is not None and tag != expected_tag:
        raise ValueError(f"tag {tag!r} differs from manifest tag {expected_tag}")
    if not dry_run:
        if tag is None:
            raise ValueError("uploads require running on the version tag, not a branch")
        gated = [name for name in LIBRARIES if packages[name].get("publish") is False]
        if gated:
            raise ValueError(f"owner decisions still gate publication: {', '.join(gated)}")
    if any(packages[name].get("publish") is not False for name in BINDINGS):
        raise ValueError("binding crates must keep publish = false")
    if npm["name"] != "nazuchem" or npm.get("private", False):
        raise ValueError("the npm package must be the public, unscoped nazuchem package")
    if (
        tomllib.loads((root / "bindings/python/pyproject.toml").read_text())["project"]["name"]
        != "nazuchem"
    ):
        raise ValueError("the PyPI project must be named nazuchem")
    print(f"All {len(versions)} manifest versions match {expected_tag}; dry_run={dry_run}")
    if output:
        with output.open("a") as stream:
            stream.write(f"tag={expected_tag}\ndry_run={str(dry_run).lower()}\n")


def package(root: Path) -> None:
    # Include uncommitted release changes in local rehearsals, excluding ignored build output.
    files = (
        subprocess.check_output(
            ["git", "ls-files", "-z", "--cached", "--others", "--exclude-standard"],
            cwd=root,
        )
        .decode()
        .split("\0")
    )
    with TemporaryDirectory(prefix="nazuchem-release-") as temporary:
        scratch = Path(temporary)
        for name in filter(None, files):
            source = root / name
            if source.is_file():
                destination = scratch / name
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, destination)
        for name in LIBRARIES:
            manifest = scratch / "crates" / name / "Cargo.toml"
            manifest.write_text(
                re.sub(r"(?m)^publish = false$", "publish = true", manifest.read_text())
            )
        # Cargo's multi-package mode resolves first-release dependencies together.
        # Archive verification currently hits rust-lang/cargo#14396 ("no hash listed").
        # The release job separately compiles/tests the full, unmodified workspace.
        subprocess.run(
            [
                "cargo",
                "package",
                "--workspace",
                "--exclude",
                "nazuchem-python",
                "--exclude",
                "nazuchem-wasm",
                "--no-verify",
                "--allow-dirty",
                "--locked",
            ],
            cwd=scratch,
            env={
                **os.environ,
                "CARGO_TARGET_DIR": os.environ.get("CARGO_TARGET_DIR", str(root / "target")),
            },
            check=True,
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("check", "package"))
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--tag")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.command == "check":
        check(args.root.resolve(), dry_run=args.dry_run, tag=args.tag, output=args.output)
    else:
        package(args.root.resolve())


if __name__ == "__main__":
    main()
