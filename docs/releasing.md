# Releasing nazuchem

The release workflow publishes the five Rust libraries, the `nazuchem` PyPI
wheels, the unscoped `nazuchem` npm package, and a GitHub Release with generated
notes, the wheels and the npm package. All versions are currently **0.1.0**.
Releases accept only numeric **X.Y.Z** versions and matching **vX.Y.Z** tags; prerelease
and build-metadata suffixes are not supported.

The owner gates were lifted on 2026-11-10 after the sort-policy, Blue Book and PubChem decisions (D7–D9).

## One-time setup

In `Ameyanagi/nazuchem`, create GitHub environments **`crates-io`**, **`pypi`**,
and **`npm`**, each with a required reviewer and deployment rules allowing
release tags (`v*`). Build/rehearsal jobs do not use these environments. Use a
public GitHub repository for npm provenance.

- **crates.io:** [current trusted publishing documentation](https://crates.io/docs/trusted-publishing)
  requires each crate to exist already; its first publication needs an API
  token. For **0.1.0**, put a short-lived token authorized to create these crates
  in the `crates-io` environment secret **`CRATES_IO_BOOTSTRAP_TOKEN`**. The
  workflow publishes `nazuchem-ir`, `nazuchem-kernel`, `nazuchem-n2s`,
  `nazuchem-s2n`, then `nazuchem`. Afterward configure a trusted publisher in
  **each crate's** Settings → Trusted Publishing: owner `Ameyanagi`, repository
  `nazuchem`, workflow **`release.yml`**, environment **`crates-io`**. Delete the
  secret and revoke its token before the next release. Future runs use
  [`crates-io-auth-action`](https://github.com/rust-lang/crates-io-auth-action),
  which obtains a short-lived token and revokes it after the job.
- **PyPI:** add an account-level
  [pending publisher](https://docs.pypi.org/trusted-publishers/creating-a-project-through-oidc/)
  for project **`nazuchem`**, owner **`Ameyanagi`**, repository **`nazuchem`**,
  workflow **`release.yml`**, environment **`pypi`**. It creates the project on
  first publication; it does not reserve the name. No token is needed.
  [`gh-action-pypi-publish`](https://github.com/pypa/gh-action-pypi-publish)
  uses OIDC and creates/uploads attestations by default. If the project already
  exists under this owner's control, add a normal trusted publisher instead.
- **npm:** the package must exist before a
  [trusted publisher can be configured](https://docs.npmjs.com/cli/v12/commands/npm-trust/#prerequisites).
  For **0.1.0**, put a short-lived granular token permitted to create/publish
  the public package, with automation-compatible 2FA bypass, in the `npm`
  environment secret **`NPM_BOOTSTRAP_TOKEN`**. After publication configure
  [trusted publishing](https://docs.npmjs.com/trusted-publishers/) for package
  **`nazuchem`**: GitHub owner **`Ameyanagi`**, repository **`nazuchem`**,
  workflow **`release.yml`**, environment **`npm`**.
  Delete the bootstrap secret and revoke its token. OIDC needs npm ≥11.5.1 and
  Node ≥22.14.0; the workflow's pinned Node includes a suitable npm. The package
  repository URL must match this repository. npm generates provenance for the
  public package/repository; the workflow also explicitly requests it.

Only the two first-publication bootstrap secrets above are needed. Existing
crates/packages can use OIDC immediately after configuration. GitHub Releases
use the job's automatic `GITHUB_TOKEN` with `contents: write`.

## Release procedure

1. Bump `[workspace.package].version`, versions in `[workspace.dependencies]`,
   `bindings/python/pyproject.toml`, `bindings/python/uv.lock`, and
   `bindings/wasm/package.json` together. Refresh `Cargo.lock` as needed. Review
   package metadata and the provenance ledger; update hashes only for edited
   recorded files. Keep all licence/NOTICE texts and source pins intact.
2. Push the reviewed release preparation branch. From Actions → Release → Run
   workflow, select that branch and leave **`dry_run=true`**. This checks every
   manifest, builds/tests the Rust workspace, packages all five libraries in a
   scratch copy, runs the existing bindings matrix and notice checks, and
   validates wheel metadata. Inspect the retained workflow artifacts. It skips
   registry uploads and GitHub Release creation.
3. After the setup above is complete, merge the reviewed
   release commit and push its tag, e.g. `git tag v0.1.0` then
   `git push origin v0.1.0`. Strict `vX.Y.Z` tags must match every manifest
   version. Approve the three protected environment jobs. A manual run with
   `dry_run=false` is also allowed, but only when selecting the matching tag.
4. Verify all five crates, PyPI wheels/attestations, npm package/provenance, and
   the GitHub Release. Remove/revoke bootstrap credentials and configure the
   trusted publishers as described above.

All build and package checks finish before the first registry upload. The
workflow does not require a prior main-CI result for the commit. Cargo already
[waits for index availability](https://doc.rust-lang.org/cargo/commands/cargo-publish.html#description)
after each publication, so dependent crates need no sleeps. A timeout may follow
a successful upload. On partial crates.io failure, inspect what was uploaded,
then use **Re-run failed jobs** in the same workflow run. The crates.io job
checks each exact version in the sparse index and skips versions already
published. Do not move a published tag or reuse a published
version.

## Build choices

`release.yml` calls `bindings.yml`, building abi3 wheels for Linux x86_64,
macOS arm64, and Windows x64, with wheel-only scope (no sdist). Intel macOS and
Linux aarch64 wheels are not built. Linux uses Maturin's
[Zig/manylinux support](https://www.maturin.rs/distribution.html#use-zig) to build
PyPI-compatible `manylinux_2_28` wheels; other platforms enforce `pypi`
compatibility. The npm artifact is the exact Bun-packed, tested V4 package,
using checksum-verified wasm-bindgen **0.2.128** and Binaryen **133**. Bun handles
install/build/pack; only upload uses npm because current
[`bun publish`](https://bun.sh/docs/pm/cli/publish) cannot generate provenance
and npm trusted publishing is not supported there yet.

Scratch Rust packaging uses Cargo's multi-package mode to resolve unpublished
dependencies in a temporary copy. Its `rust-crates-preflight` artifact is for rehearsal; the GitHub
Release attaches only the wheels and the npm package. Preflight
uses `--no-verify` because Cargo 1.99.1 archive verification
hits the existing temporary-registry checksum failure
([cargo#14396](https://github.com/rust-lang/cargo/issues/14396)). The original
workspace is compiled/tested separately. Per-crate `cargo publish` verification
against crates.io is first exercised on the real tag run.
