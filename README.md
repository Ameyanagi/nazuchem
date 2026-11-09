# nazuchem

Rust library for IUPAC-style chemical name ↔ structure conversion. It accepts
names and molecular SMILES through the current string APIs. Python and npm/WASM
bindings and a versioned structure IR are planned.

Extracted from [ReShiki](https://github.com/Ameyanagi/ReShiki), a molecular
drawing app. The name-to-structure engine is a Rust port of OPSIN by Daniel Lowe
(MIT). The structure-to-name engine began as a port of openclatura 0.4.3 by
lamalab-org (MIT), extended in ReShiki. The chemistry kernel adapts RDKit and
SiMolecule centres, with LLVM libc++ and Microsoft STL sorting adaptations;
the naming engine also includes a CPython integer-set adaptation.

This is a pre-release extraction: every crate has `publish = false`, and the
0.x API is not yet stable.

```rust
fn main() -> Result<(), nazuchem::Error> {
    // Equivalent to CCO; atom order follows the parser.
    assert_eq!(nazuchem::name_to_smiles("ethanol")?, "C(C)O");
    let name = nazuchem::smiles_to_name("CCO")?;
    assert_eq!(name.text, "ethanol");
    Ok(())
}
```

The facade enables `n2s` and `s2n` by default; either feature can be selected
independently with default features disabled. `name_to_smiles` applies strict
OPSIN options, rejects warnings and returns semantic CXSMILES.
`smiles_to_name` returns the complete `Name`, including its preference
assessment, alternatives and omitted stereochemistry.

| Crate | Purpose |
| --- | --- |
| [`nazuchem`](crates/nazuchem) | Public string-API facade and engine re-exports |
| [`nazuchem-ir`](crates/nazuchem-ir) | Shared naming vocabulary; structure IR planned |
| [`nazuchem-kernel`](crates/nazuchem-kernel) | Internal molecular perception and canonical identity |
| [`nazuchem-n2s`](crates/nazuchem-n2s) | OPSIN name-to-structure port |
| [`nazuchem-s2n`](crates/nazuchem-s2n) | Extended openclatura structure-to-name port |

## Scope and limits

The pinned OPSIN revision and supported s2n rules do not cover every chemical
name or structure. A successful generated name is not necessarily a preferred
IUPAC name. Unsupported chemistry returns an error; OPSIN compatibility data
does not certify IUPAC PIN compliance. See the [facade documentation](crates/nazuchem/README.md)
for admission limits and unsupported stereochemistry.

The library performs no I/O, creates no threads and has no process-wide mutable
naming state. Embedded parser resources are immutable; fusion uses bounded
thread-local caches. Default budgets are unbounded. Hosts can inject a `Budget`
into the engine APIs; its analysis timeout covers kernel perception, not an
entire naming request. Worst-case time on adversarial input is unbounded; hosts
requiring hard limits must supervise requests externally.

## Development

Use the pinned Rust toolchain and checked lockfile:

```sh
export CARGO_TARGET_DIR="$PWD/target"
cargo fmt --all --check
cargo clippy --workspace --all-targets --locked --offline -- -D warnings
cargo test --workspace --locked --offline
python3 scripts/check_provenance.py
python3 -m unittest tests.test_check_provenance
```

The golden gate reads `bakeoff_corpus.tsv`, `openclatura_tests.tsv` and
`rules_cases.tsv` from `tests/fixtures/s2n/`, in that order, keeping each SMILES
at its first occurrence. It captures the complete s2n `Name` or error and
strict n2s results for the main name, alternatives and any recovered neutral
verification name, deduplicated in first-occurrence order within each row, in
[`tests/fixtures/golden/s2n_n2s.jsonl`](tests/fixtures/golden/s2n_n2s.jsonl).
Review intentional behavior changes before regenerating it; do not edit rows
by hand:

```sh
NAZUCHEM_UPDATE_GOLDEN=1 cargo test -p nazuchem --test golden --locked --offline
```

Frozen OPSIN fixtures run without Java or a live oracle. The
[extraction plan](docs/plans/naming-split.md) records the migration decisions
and later IR/bindings work.

## Licensing and provenance

Original code uses [MIT](LICENSE-MIT) OR [Apache-2.0](LICENSE-APACHE), except
original additions inside `nazuchem-n2s`, which are MIT only like the OPSIN
port they extend.
Incorporated material keeps its own terms: OPSIN/openclatura MIT, RDKit
BSD-3-Clause, centres BSD-2-Clause, LLVM libc++/Microsoft STL
Apache-2.0 WITH LLVM-exception, and CPython Python-2.0. See [NOTICE](NOTICE),
the per-file [`PROVENANCE-nazuchem.toml`](PROVENANCE-nazuchem.toml) ledger and
each crate's `NOTICE` and `licenses/` texts for the applicable scope and pins.

Owner decisions remain pending for native sort-policy provenance, IUPAC Blue
Book table reuse permission and the test-only PubChem-derived corpus's source
terms. These entries retain `NOASSERTION`; the established license expressions
do not grant rights for that material. Publishing remains disabled.
