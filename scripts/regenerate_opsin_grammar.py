"""Regenerate the n2s DFA table and every hash that pins it from the checked grammar.

Run after editing the grammar XML in crates/nazuchem-n2s/resources/xml (normally
regexTokens.xml or regexes.xml), or after a merge conflict in the generated files:
merge the XML and any prose by hand, then rerun. Conflict hunks that differ only in
hashes are resolved here, and git keeps one side of automata.bin. --check verifies
without writing. See crates/nazuchem-n2s/provenance/README.md.
"""

import argparse
import hashlib
import os
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CRATE = ROOT / "crates/nazuchem-n2s"
MANIFEST = CRATE / "provenance/manifest.json"
LEDGER = ROOT / "PROVENANCE-nazuchem.toml"
AUTOMATA = CRATE / "resources/automata.bin"
GRAMMAR = CRATE / "resources/xml"
RESOURCE_DIR = "uk/ac/cam/ch/wwmm/opsin/resources"
# ExportDfa builds every automaton from these files alone (token files such as
# multipliers.xml feed token dictionaries, not automata). Each one that differs from the
# jar's upstream copy is overridden, and nothing else: freshly built automata (token-69
# here) get a different but isomorphic state numbering when more files load from disk.
DFA_GRAMMAR = ("regexes.xml", "regexes.dtd", "regexTokens.xml", "regexTokenList.dtd")
# The original oracle build. Maven jars embed build timestamps, so a rebuild of the
# pinned commit is accepted by JAR_CONTENT_SHA256 (see jar_content_sha256) instead.
JAR_SHA256 = "931be2952b290dea62548a3bb4129703078c763576a836ac99d8a54d30083f03"
JAR_CONTENT_SHA256 = "84949d143459a3de91b9d0483910832e70ee26c7a550380a012250a878c09de5"
HEX = "[0-9a-f]{64}"
CONFLICT = re.compile(
    r"<<<<<<< .*?\n(.*?)(?:\|{7} .*?\n.*?)?=======\n(.*?)>>>>>>> .*?\n", re.DOTALL
)


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def jar_content_sha256(jar):
    """SHA256 over sorted entry names and uncompressed bytes, ignoring zip metadata."""
    digest = hashlib.sha256()
    with zipfile.ZipFile(jar) as archive:
        for name in sorted(archive.namelist()):
            digest.update(name.encode() + b"\0" + archive.read(name) + b"\0")
    return digest.hexdigest()


def export_dfa(jar, javac, java, out):
    """Run the recorded exporter with each changed DFA grammar file as a resource override."""
    with zipfile.ZipFile(jar) as archive:
        overrides = [
            name
            for name in DFA_GRAMMAR
            if archive.read(f"{RESOURCE_DIR}/{name}") != (GRAMMAR / name).read_bytes()
        ]
    print(f"resource overrides: {', '.join(overrides) or 'none'}")
    with tempfile.TemporaryDirectory() as tmp:
        classes, work = Path(tmp, "classes"), Path(tmp, "work")
        classes.mkdir()
        (work / "resources" / RESOURCE_DIR).mkdir(parents=True)
        for name in overrides:
            shutil.copy(GRAMMAR / name, work / "resources" / RESOURCE_DIR / name)
        shutil.copy(CRATE / "provenance/ExportDfa.java.txt", classes / "ExportDfa.java")
        subprocess.run(
            [javac, "-cp", jar, "-d", classes, classes / "ExportDfa.java"], check=True
        )
        subprocess.run(
            [
                java,
                "-cp",
                f"{jar}{os.pathsep}{classes}",
                "uk.ac.cam.ch.wwmm.opsin.ExportDfa",
                out,
            ],
            cwd=work,
            check=True,
        )


def read(path):
    """Read `path`, taking our side of every merge conflict that differs only in hashes."""

    def resolve(hunk):
        ours, theirs = hunk.groups()
        return ours if re.sub(HEX, "", ours) == re.sub(HEX, "", theirs) else hunk[0]

    text = CONFLICT.sub(resolve, path.read_bytes().decode())
    if "<<<<<<< " in text:
        raise SystemExit(
            f"Resolve the remaining merge conflicts in {path.relative_to(ROOT)}."
        )
    return text


def pin(text, pattern, value, path):
    """Replace the one hash that follows `pattern` (a regex prefix) in `text`."""
    text, count = re.subn(f"({pattern}){HEX}", rf"\g<1>{value}", text)
    if count != 1:
        raise SystemExit(
            f"{path.relative_to(ROOT)}: expected one hash after {pattern!r}."
        )
    return text


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--jar",
        type=Path,
        default=os.environ.get("OPSIN_JAR"),
        help="pinned opsin-core-3.0-SNAPSHOT-jar-with-dependencies.jar (default: $OPSIN_JAR)",
    )
    parser.add_argument(
        "--java-home",
        type=Path,
        default=os.environ.get("JAVA_HOME"),
        help="JDK 21 home containing bin/javac (default: $JAVA_HOME)",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="fail if any generated file is stale instead of writing",
    )
    args = parser.parse_args()

    if args.jar is None or not args.jar.is_file():
        raise SystemExit(
            "Pass --jar or set OPSIN_JAR to opsin-core-3.0-SNAPSHOT-jar-with-dependencies.jar "
            "built from the pinned OPSIN commit; see crates/nazuchem-n2s/provenance/README.md."
        )
    if (
        sha256(args.jar.read_bytes()) != JAR_SHA256
        and (digest := jar_content_sha256(args.jar)) != JAR_CONTENT_SHA256
    ):
        raise SystemExit(
            f"{args.jar} is neither the pinned oracle (SHA256 {JAR_SHA256}) nor a rebuild "
            f"with jar content SHA256 {JAR_CONTENT_SHA256}; its content SHA256 is {digest}."
        )
    bin_dir = args.java_home and args.java_home.resolve() / "bin"
    javac, java = (
        bin_dir and shutil.which(tool, path=bin_dir) for tool in ("javac", "java")
    )
    if not (javac and java):
        raise SystemExit(
            "Pass --java-home or set JAVA_HOME to a JDK 21 home with bin/javac."
        )

    # Read every text file first, so unresolved conflicts fail before the slow export.
    for path in GRAMMAR.iterdir():
        read(path)
    manifest = read(MANIFEST)
    ledger = read(LEDGER)
    readme = CRATE / "provenance/README.md"
    pinned = [
        (path, read(path), pattern)
        for path, pattern in [
            (CRATE / "src/lib.rs", r'RESOURCE_FINGERPRINT: &str =\s*"'),
            (readme, r"exported crate constant is\s*`"),
        ]
    ]
    with tempfile.TemporaryDirectory() as tmp:
        table = Path(tmp, "automata.bin")
        export_dfa(args.jar.resolve(), javac, java, table)
        generated = {AUTOMATA: table.read_bytes()}

    # Hash every resource listed in the manifest, then the fingerprint over them in
    # sorted path order: path, NUL, raw bytes, NUL.
    fingerprint = hashlib.sha256()
    resources = sorted(re.findall(rf'"(resources/[^"]+)": "{HEX}"', manifest))
    for path in resources:
        data = generated.get(CRATE / path) or (CRATE / path).read_bytes()
        manifest = pin(manifest, f'"{re.escape(path)}": "', sha256(data), MANIFEST)
        fingerprint.update(path.encode() + b"\0" + data + b"\0")
    fingerprint = fingerprint.hexdigest()
    manifest = pin(
        manifest, r'"resource_fingerprint": \{[^}]*"sha256": "', fingerprint, MANIFEST
    )
    generated[MANIFEST] = manifest.encode()
    for path, text, pattern in pinned:
        generated[path] = pin(text, pattern, fingerprint, path).encode()

    # The repository ledger records the file hash of every resource, the manifest and
    # the README; Rust sources such as lib.rs are covered by headers, not hashes.
    for path in [*(CRATE / p for p in resources), MANIFEST, readme]:
        data = generated.get(path) or path.read_bytes()
        entry = rf'(?s:path = "{re.escape(path.relative_to(ROOT).as_posix())}"\n(?:(?!\[\[file\]\]).)*?sha256 = ")'
        ledger = pin(ledger, entry, sha256(data), LEDGER)
    generated[LEDGER] = ledger.encode()

    stale = [
        path
        for path, content in generated.items()
        if not path.exists() or path.read_bytes() != content
    ]
    for path in stale:
        print(f"{'stale' if args.check else 'updated'}: {path.relative_to(ROOT)}")
        if not args.check:
            path.write_bytes(generated[path])
    print(f"resource fingerprint {fingerprint}")
    if args.check and stale:
        sys.exit("Generated n2s grammar files are stale; rerun without --check.")


if __name__ == "__main__":
    main()
