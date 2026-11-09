# SPDX-License-Identifier: MIT OR Apache-2.0
# Origin: original
# Source: Nazuchem | tests/test_binding_notices.py | Phase B step 9
"""Notice collection must follow dependency edges and fail on missing texts."""

import tempfile
import unittest
from pathlib import Path

from scripts.collect_binding_notices import dependency_closure, render_notices


class BindingNoticeTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        root = Path(temporary.name)
        self.packages = []
        for name in ("binding", "runtime", "build", "dev"):
            directory = root / name
            directory.mkdir()
            (directory / "LICENSE-MIT").write_text(f"{name} licence\n", encoding="utf-8")
            self.packages.append(
                {
                    "id": name,
                    "name": name,
                    "version": "1.0.0",
                    "license": "MIT",
                    "manifest_path": str(directory / "Cargo.toml"),
                }
            )
        self.metadata: dict = {
            "packages": self.packages,
            "resolve": {
                "nodes": [
                    {
                        "id": "binding",
                        "deps": [
                            {"pkg": name, "dep_kinds": [{"kind": kind}]}
                            for name, kind in (
                                ("runtime", None),
                                ("build", "build"),
                                ("dev", "dev"),
                            )
                        ],
                    },
                    *({"id": name, "deps": []} for name in ("runtime", "build", "dev")),
                ]
            },
        }

    def test_normal_and_build_closure_excludes_test_only_dependencies(self):
        self.assertEqual(
            [item["name"] for item in dependency_closure(self.metadata, "binding")],
            ["build", "runtime"],
        )
        data = render_notices(self.metadata, "binding", "test-target")
        self.assertIn(b"runtime licence\n", data)
        self.assertIn(b"build licence\n", data)
        self.assertNotIn(b"dev licence", data)

    def test_nested_terms_and_explicit_license_file_are_preserved(self):
        directory = Path(self.packages[1]["manifest_path"]).parent
        (directory / "licenses/upstream").mkdir(parents=True)
        (directory / "licenses/upstream/NOTICE").write_text("upstream notice", encoding="utf-8")
        (directory / "custom-terms.txt").write_text("custom licence", encoding="utf-8")
        self.packages[1]["license_file"] = "custom-terms.txt"
        data = render_notices(self.metadata, "binding", "test-target")
        self.assertIn(b"upstream notice", data)
        self.assertIn(b"custom licence", data)

    def test_empty_closure_and_missing_licences_fail(self):
        dependencies = self.metadata["resolve"]["nodes"][0]["deps"]
        self.metadata["resolve"]["nodes"][0]["deps"] = []
        with self.assertRaisesRegex(ValueError, "no dependency notices"):
            render_notices(self.metadata, "binding", "test-target")
        self.metadata["resolve"]["nodes"][0]["deps"] = dependencies
        (Path(self.packages[1]["manifest_path"]).parent / "LICENSE-MIT").unlink()
        with self.assertRaisesRegex(ValueError, "no licence/notice files"):
            render_notices(self.metadata, "binding", "test-target")


if __name__ == "__main__":
    unittest.main()
