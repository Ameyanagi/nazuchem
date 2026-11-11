# nazuchem

Chemical names to a typed structure IR or semantic CXSMILES, and molecular
IR or SMILES to names with preference assessments.

The name-to-structure engine is a Rust port of OPSIN by Daniel Lowe (MIT).
The structure-to-name engine began as a port of openclatura 0.4.3, extended
in ReShiki. The chemistry kernel contains RDKit and SiMolecule centres
adaptations; an internal copy preserves their upstream provenance.

Extracted from [ReShiki](https://github.com/Ameyanagi/ReShiki), a molecular drawing app.

```rust
fn main() -> Result<(), nazuchem::Error> {
    #[cfg(feature = "n2s")]
    {
        // Equivalent to CCO; atom order follows the parser.
        assert_eq!(nazuchem::name_to_smiles("ethanol")?, "C(C)O");
        let molecule = nazuchem::name_to_molecule("ethanol")?;
        let json = molecule.to_json().expect("validated IR JSON");
        assert_eq!(nazuchem::Molecule::from_json(&json).unwrap(), molecule);
        #[cfg(feature = "s2n")]
        assert_eq!(nazuchem::molecule_to_name(&molecule)?.text, "ethanol");
    }
    #[cfg(feature = "s2n")]
    {
        let name = nazuchem::smiles_to_name("CCO")?;
        assert_eq!(name.text, "ethanol");
    }
    Ok(())
}
```

The default features are `n2s` and `s2n`. Disable default features and select
either independently. `nazuchem::n2s` and `nazuchem::s2n` expose the engine
APIs; IR types, shared vocabulary and `Budget` are re-exported at the facade root.
The IR exports include `AxialWinding`, `CxsmilesError`, `IR_VERSION`, and
`JSON_SCHEMA`, and are available even with both engine features disabled.
`name_to_smiles` uses strict OPSIN options, rejects parser warnings and returns
semantic CXSMILES. Axial configurations require the engine API that returns
axes separately. `smiles_to_name` uses IUPAC isotope notation and returns the
engine's complete `Name` result. `molecule_to_name` applies the same defaults
directly to an IR graph; it is available with `s2n` alone. The lower-level
`s2n::name_molecule` accepts an isotope style and a caller-supplied `Budget`.

`name_to_molecule` applies the same strict parsing policy and returns validated
IR v1, with stable atom identities, all locants and the principal-parent set.
`Molecule::to_json`/`from_json` use a strict versioned envelope;
`Molecule::to_cxsmiles` provides a convenience export using the
[documented CX extensions](../nazuchem-ir/docs/cxsmiles.md). JSON is lossless;
unrepresentable CX stereo returns `CxsmilesError::StereoNotExpressible`. IR v1
includes atropisomer axes; allene axial chirality is not representable in v1.
The pinned OPSIN parser erases explicit ξ/unknown
name descriptors before conversion, a separately tracked parser limitation.
The IR-to-name adapter validates the graph and chemical admission constraints
before naming, retaining the caller's atom indices. Locants, principal-parent
membership and ordinary atom maps do not constrain the generated name.
Annotations the naming engine cannot represent, including text labels,
attachment points, position variation, polymers and enhanced stereo groups,
return errors instead of being silently discarded. There is no CXSMILES reader.

IR v1 freezes at the first published release. Until that release, its fields and
variants may change without a version bump; after freezing, new fields or
variants require a new IR version. `Eq` compares the complete ordered structure
and annotations, not stereochemical equivalence under atom or reference
reordering.

Atom JSON omits `charge: 0`, `radical_electrons: 0`, and
`hydrogens_stated: false`; readers restore those defaults. Existing optional
fields also default on input and omit empty or `None` annotations on output.
For example, an otherwise unannotated methane atom is
`{"element":6,"hydrogens":4}`. `AxialWinding` uses `"clockwise"`,
`"anticlockwise"`, or `"unknown"`, relative to its ordered references rather
than CIP priorities: looking from `bond.atoms[0]` toward `bond.atoms[1]`, the
shortest turn from `refs[0]` to `refs[2]` has the stated winding. Swapping one
reference pair reverses it; reversing the axis and exchanging both pairs
preserves it. Stereo groups may contain only known configurations.

## Scope and limits

`n2s` interprets supported IUPAC/CAS-style systematic and retained names
according to the pinned OPSIN 3.0-SNAPSHOT revision; it does not cover every
chemical name or promise newer OPSIN behavior. `s2n` names chains,
monocycles, retained and systematic fused parents, spiro and von Baeyer ring
systems, supported ions and charge-balanced salts, and adducts of these with
neutral components such as acids and water. It also names supported
neutral mono-/diradicals and radical ions when a closed-shell neutral parent
can be recovered from one charge and one unpaired electron on the same atom
or the two ends of a broken pi bond.

`s2n` SMILES input is limited to 32,768 bytes and 512 atom tokens, with 1–512
atoms after import; IR input has the same 1–512 atom limit. Naming also retains
the native atom-map admission ceiling of 2,147,483,639, although IR JSON can
store the full `u32` range. It rejects CXSMILES extensions, reaction SMILES, allene and
coordination stereo classes, zwitterions mixed with ions,
and isotope-labelled neutral radicals. Unsupported chemistry or naming paths
return `Err`; comprehensive IUPAC coverage is not claimed. A successful name
includes a preference assessment and is not necessarily the preferred IUPAC
name. This is a pre-release 0.x API; compatibility is not yet guaranteed.

The convenience functions use `Budget::default()`, which is unbounded.
Hosts can pass explicit budgets to `n2s::Parser::with_budget` and the s2n
entry points. A budget's optional analysis timeout covers kernel molecular
perception only; it is not a deadline for the complete naming operation.
The library performs no I/O, creates no threads and has no
process-wide mutable naming state. The parser shares immutable embedded
resources; fusion retains bounded thread-local caches. Worst-case time on
adversarial input is unbounded; hosts needing hard limits should supervise
requests externally.

Original facade code is MIT OR Apache-2.0. The package expression also ANDs
the engines' and kernel's established terms: MIT (OPSIN and openclatura),
BSD-3-Clause (RDKit), BSD-2-Clause (centres), and Python-2.0 (a CPython adaptation).
Complete license texts, copyright notices and dependency scope notices are
included in `licenses/`; `LICENSE-MIT` and `LICENSE-APACHE` cover local code.

Sorts with equal keys keep their input order (Rust's stable sort); nazuchem
does not reproduce C++ standard-library tie order.

`nazuchem-s2n`'s `src/data/bluebook/p25_retained.tsv` lists the retained heterocycle
names of IUPAC Nomenclature of Organic Chemistry: Recommendations and Preferred
Names 2013, P-25.2.1 Table 2.8, transcribed from
<https://iupac.qmul.ac.uk/BlueBook/P2.html>; IUPAC holds copyright in the
Recommendations. The names are included as factual data with this attribution;
this package's licences do not cover the file, and no licence is asserted for it
(SPDX NOASSERTION).
The crate's `NOTICE` and copied dependency notices describe the exact scope.
