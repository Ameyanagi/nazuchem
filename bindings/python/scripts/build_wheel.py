# SPDX-License-Identifier: MIT OR Apache-2.0
# Origin: original
# Source: Nazuchem | bindings/python/scripts/build_wheel.py | Phase B step 9

"""Generate the target's notices, build an abi3 wheel, and verify its contents."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

PACKAGE_DIR = Path(__file__).resolve().parents[1]
ROOT = PACKAGE_DIR.parents[1]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", help="Rust target triple (defaults to rustc host)")
    parser.add_argument("--release", action="store_true", help="enable Rust release optimization")
    parser.add_argument("--out", type=Path, default=PACKAGE_DIR / "dist")
    args = parser.parse_args()
    environment = os.environ.copy()
    environment.setdefault("CARGO_TARGET_DIR", str(ROOT / "target"))
    collect = [
        sys.executable,
        str(ROOT / "scripts/collect_binding_notices.py"),
        "--package",
        "nazuchem-python",
        "--output-dir",
        str(PACKAGE_DIR),
    ]
    if args.target:
        collect.extend(["--target", args.target])
    subprocess.run(collect, cwd=ROOT, env=environment, check=True)
    args.out.mkdir(parents=True, exist_ok=True)
    # A fresh output directory prevents a stale wheel from passing verification.
    with TemporaryDirectory(prefix="nazuchem-wheels-") as temporary:
        build = [
            "uv",
            "run",
            "--no-project",
            "--with",
            "maturin>=1.8,<2",
            "maturin",
            "build",
            "--locked",
            "--manifest-path",
            str(PACKAGE_DIR / "Cargo.toml"),
            "--interpreter",
            sys.executable,
            "--out",
            temporary,
        ]
        if args.target:
            build.extend(["--target", args.target])
        if args.release:
            build.append("--release")
        subprocess.run(build, cwd=PACKAGE_DIR, env=environment, check=True)
        wheels = sorted(Path(temporary).glob("*.whl"))
        if not wheels:
            raise SystemExit("maturin did not produce a wheel")
        for wheel in wheels:
            check = [
                sys.executable,
                str(PACKAGE_DIR / "scripts/check_wheel_notices.py"),
                str(wheel),
            ]
            if args.target:
                check.extend(["--target", args.target])
            subprocess.run(check, cwd=ROOT, env=environment, check=True)
            destination = args.out / wheel.name
            shutil.copyfile(wheel, destination)
            print(f"Built and checked {destination}")


if __name__ == "__main__":
    main()
