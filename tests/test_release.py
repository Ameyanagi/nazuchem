# SPDX-License-Identifier: MIT OR Apache-2.0
# Origin: original
# Source: Nazuchem | tests/test_release.py

"""Exercise release admission before any registry job can run."""

import shutil
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import tomllib

from scripts.check_release import LIBRARIES, ROOT, check

VERSION = tomllib.loads((ROOT / "Cargo.toml").read_text())["workspace"]["package"]["version"]
TAG = f"v{VERSION}"


class ReleaseGateTests(unittest.TestCase):
    def setUp(self):
        self.temporary = TemporaryDirectory(prefix="nazuchem-release-test-")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        manifests = [
            ROOT / "Cargo.toml",
            *ROOT.glob("crates/*/Cargo.toml"),
            *ROOT.glob("bindings/*/Cargo.toml"),
        ]
        manifests += [ROOT / "bindings/python/pyproject.toml", ROOT / "bindings/wasm/package.json"]
        for source in manifests:
            target = self.root / source.relative_to(ROOT)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, target)

    def admit(self, *, dry_run=True, tag=None):
        check(self.root, dry_run=dry_run, tag=tag, output=None)

    def release_gates(self):
        for name in LIBRARIES:
            path = self.root / "crates" / name / "Cargo.toml"
            path.write_text(path.read_text().replace("[package]\n", "[package]\npublish = true\n", 1))

    def test_versions_require_three_ascii_numeric_parts(self):
        path = self.root / "Cargo.toml"
        original = path.read_text()
        for version in ("1.2", "1.2.3-rc.1", "1.2.3+build", "١.٢.٣"):
            path.write_text(original.replace(f'version = "{VERSION}"', f'version = "{version}"'))
            with self.subTest(version), self.assertRaisesRegex(ValueError, "invalid release"):
                self.admit()

    def test_tag_mismatch_is_rejected_even_in_rehearsal(self):
        with self.assertRaisesRegex(ValueError, "tag .* differs"):
            self.admit(tag="v9.9.9")

    def test_each_manifest_version_is_checked(self):
        for relative in (
            "bindings/python/pyproject.toml",
            "bindings/wasm/package.json",
            "crates/nazuchem-ir/Cargo.toml",
        ):
            path = self.root / relative
            original = path.read_text()
            if relative.endswith("Cargo.toml"):
                changed = original.replace("version.workspace = true", 'version = "9.9.9"')
            else:
                changed = original.replace(VERSION, "9.9.9")
            path.write_text(changed)
            with (
                self.subTest(manifest=relative),
                self.assertRaisesRegex(ValueError, "version .* differs"),
            ):
                self.admit(tag=TAG)
            path.write_text(original)

    def test_uploads_require_a_tag_and_every_library_gate_removed(self):
        self.release_gates()
        with self.assertRaisesRegex(ValueError, "uploads require"):
            self.admit(dry_run=False)
        for name in LIBRARIES:
            path = self.root / "crates" / name / "Cargo.toml"
            original = path.read_text()
            path.write_text(original.replace("publish = true", "publish = false"))
            with self.subTest(library=name), self.assertRaisesRegex(ValueError, "owner decisions"):
                self.admit(dry_run=False, tag=TAG)
            path.write_text(original)

    def test_rehearsal_and_ready_tag_succeed(self):
        self.admit()
        self.release_gates()
        output = self.root / "output"
        check(self.root, dry_run=False, tag=TAG, output=output)
        self.assertEqual(output.read_text(), f"tag={TAG}\ndry_run=false\n")


if __name__ == "__main__":
    unittest.main()
