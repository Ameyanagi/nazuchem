# nazuchem-kernel

Internal molecular perception kernel for nazuchem. The standalone library copies
the transitive chemistry modules used by OPSIN and openclatura from ReShiki at
`be63a799aa03a83dc6ec74bb294f4d42b89298d3`, retaining the algorithms,
attribution comments, and relevant unit tests except for the sorting change
described below.

The initial extraction inventory recorded 88 copied Rust module files and
27,886 copied source lines, plus the required `AtomFacts`, `Element`, atomic
table exports, module declarations, and crate attributes from the original root. The unused
generated-wedging module and the host allocation-measurement test stay outside
the extraction. That inventory predates removal of the six sort module files
and three sort fingerprint fixtures. The remaining ring-order fixture retains
its original bytes and recorded provenance.

Owner decision 2026-11-10, **Rust stable sort**: ordering uses Rust standard
library stable sort calls. Equal keys preserve input order. Reproducing the
per-OS C++ standard libraries' equal-key permutations is no longer a target.
The native sort reproductions, shared helpers and C++ fingerprints were removed;
no Rust standard library implementation was copied.

It also includes ReShiki's pure canonical identity, private M/P axis
marker codec, and pseudoasymmetric-centre classifier. These operate on kernel
states; identity and axis helpers take an explicit `Budget`. Document selection,
component admission and caption policy belong to the host application and are
outside this library. ReShiki's integration parity tests are maintained there.

Molecular perception is adapted from RDKit under BSD-3-Clause. Extended
tetrahedral CIP code is adapted from SiMolecule centres under BSD-2-Clause.
Original ReShiki/nazuchem code uses MIT OR Apache-2.0. Complete applicable
third-party texts and notices travel with this crate in `licenses/rdkit/`
and `licenses/centres/`; `LICENSE-MIT` and `LICENSE-APACHE` cover the local
contributions. Historical ReShiki paths in copied notices identify the source
repository.

The manifest combines these established terms as
`(MIT OR Apache-2.0) AND BSD-3-Clause AND BSD-2-Clause`. Repository
provenance work remains incomplete, and `publish = false` remains in force.
`PROVENANCE-nazuchem.toml` at the repository root records per-file ancestry,
local function boundaries and the unresolved facts; the crate's `NOTICE`
explains its incorporated scope.

`budget::Budget` is explicit at entry points that perform aromaticity
perception, including SMILES import, sanitization, and stereo axis validation.
`Budget::default()` is unbounded. A host can supply heap and work limits, an
optional analysis timeout, and a reservation callback whose returned guard is
held for the operation. The kernel does not detect machine capabilities or
maintain a process-wide reservation pool. `budget::Cancellation` supports
cooperative cancellation through `aromaticity::perceive_with_budget`.

```rust
use nazuchem_kernel::{budget::Budget, smiles};

let budget = Budget::default();
let molecule = smiles::read("C[C@H](O)C(=O)O", &budget)?;
let canonical = smiles::write::write(&molecule.prepared.state, Default::default())?;
# Ok::<(), Box<dyn std::error::Error>>(())
```

This crate performs no I/O and creates no threads. Its resource policy has no
global mutable state. Existing structural and algorithmic limits remain in
place; worst-case time on adversarial input is unbounded. Hosts needing hard
limits should supervise the operation externally.
