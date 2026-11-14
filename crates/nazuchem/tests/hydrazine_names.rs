// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | crates/nazuchem/tests/hydrazine_names.rs | issue #168

#![cfg(all(feature = "n2s", feature = "s2n"))]

//! Issue #168: a nitro group on an amine nitrogen no longer costs the parent
//! chain a carbon, acylhydrazones are named as N'-ylidene hydrazides, and the
//! prefix of H2N-N= is hydrazinylidene (IUPAC 2013 recommendations). Every
//! generated name reads back as the structure it was generated from.

use nazuchem::{
    Budget, NameReason, NameStatus, name_to_smiles, s2n::canonical_smiles, smiles_to_name,
};

fn check(smiles: &str, expected: &str, status: NameStatus) {
    let name = smiles_to_name(smiles).unwrap_or_else(|e| panic!("{smiles}: {e:?}"));
    assert_eq!(name.text, expected, "{smiles}");
    assert_eq!(name.status, status, "{smiles}");
    let parsed = name_to_smiles(expected).unwrap_or_else(|e| panic!("{expected}: {e:?}"));
    let budget = Budget::default();
    assert_eq!(
        canonical_smiles(&parsed, &budget).unwrap(),
        canonical_smiles(smiles, &budget).unwrap(),
        "{expected}"
    );
}

#[test]
fn a_nitro_group_on_nitrogen_keeps_the_parent_chain() {
    // P-67.1.4.3.2: -NH-NO2 is the preselected prefix nitramido, other
    // derivatives substitute nitro into amino, and a nitro group on another
    // nitrogen does not lengthen a chain. The acid or ester stays the
    // principal characteristic group (P-41) with its whole chain. Without a
    // senior group the compound is a nitramide (P-67.1.2.6.3, issue #176).
    for (smiles, expected) in [
        ("OC(=O)CCN[N+](=O)[O-]", "3-nitramidopropanoic acid"),
        ("OC(=O)CCNN(=O)=O", "3-nitramidopropanoic acid"),
        ("OC(=O)CN[N+](=O)[O-]", "nitramidoacetic acid"),
        ("OC(=O)CCCCN[N+](=O)[O-]", "5-nitramidopentanoic acid"),
        ("OC(=O)CC(N[N+](=O)[O-])C", "3-nitramidobutanoic acid"),
        (
            "OC(=O)C(C)CN[N+](=O)[O-]",
            "2-methyl-3-nitramidopropanoic acid",
        ),
        (
            "OC(=O)CCN(C)[N+](=O)[O-]",
            "3-[methyl(nitro)amino]propanoic acid",
        ),
        (
            "OC(=O)CCN([N+](=O)[O-])[N+](=O)[O-]",
            "3-(dinitroamino)propanoic acid",
        ),
        ("CCOC(=O)CCN[N+](=O)[O-]", "ethyl 3-nitramidopropanoate"),
        (
            "OC(=O)CCON[N+](=O)[O-]",
            "3-[(nitroamino)oxy]propanoic acid",
        ),
        // Names that were already right.
        ("NC(=O)CCN[N+](=O)[O-]", "3-nitramidopropanamide"),
        ("OC(=O)CCNN=O", "3-(nitrosoamino)propanoic acid"),
        // N(OH)2 is not a nitro group.
        (
            "OC(=O)CCNN(O)O",
            "3-(2,2-dihydroxyhydrazin-1-yl)propanoic acid",
        ),
    ] {
        check(smiles, expected, NameStatus::Preferred);
    }
}

#[test]
fn acylhydrazones_are_ylidene_hydrazides() {
    for (smiles, expected) in [
        // P-66.3.1.1: hydrazides are named with the hydrazide suffixes, not
        // as acyl derivatives of hydrazine. P-66.3.3.1: substituents on the
        // nitrogen atoms take N and N'; an ylidene on N' as in
        // N'-[(phenyldiazenyl)methylidene]acetohydrazide (P-68.3.1.3.5.1).
        ("CC(=O)NN=C(C)C", "N'-(propan-2-ylidene)acetohydrazide"),
        ("CC(C)=NNC(C)=O", "N'-(propan-2-ylidene)acetohydrazide"),
        ("CC(=O)NN=CC", "N'-ethylideneacetohydrazide"),
        ("CC(=O)NN=C1CCCCC1", "N'-cyclohexylideneacetohydrazide"),
        (
            "CC(=O)N(C)N=C(C)C",
            "N-methyl-N'-(propan-2-ylidene)acetohydrazide",
        ),
        ("CCCC(=O)NN=C(C)C", "N'-(propan-2-ylidene)butanehydrazide"),
        (
            "O=C(NN=C(C)C)c1ccccc1",
            "N'-(propan-2-ylidene)benzohydrazide",
        ),
        (
            "O=C(NN=C(C)C)C1CCCCC1",
            "N'-(propan-2-ylidene)cyclohexanecarbohydrazide",
        ),
        (
            "CS(=O)(=O)NN=C(C)C",
            "N'-(propan-2-ylidene)methanesulfonohydrazide",
        ),
        // P-66.3.4: chalcogen analogues.
        ("CC(=S)NN=C(C)C", "N'-(propan-2-ylidene)ethanethiohydrazide"),
        // P-29.6.2.1: benzylidene is the preferred prefix when
        // unsubstituted, as in N'-benzylideneoxamic hydrazide (P-66.3.6).
        ("CC(=O)NN=Cc1ccccc1", "N'-benzylideneacetohydrazide"),
        (
            "CS(=O)(=O)NN=Cc1ccccc1",
            "N'-benzylidenemethanesulfonohydrazide",
        ),
        (
            "CC(=O)NN=Cc1ccc(Cl)cc1",
            "N'-[(4-chlorophenyl)methylidene]acetohydrazide",
        ),
        ("CC(=O)N=Cc1ccccc1", "N-benzylideneacetamide"),
        (
            "OC(=O)C1CCC(CC1)=Cc1ccccc1",
            "4-benzylidenecyclohexane-1-carboxylic acid",
        ),
        // A cumulated carbon on N' is not an ylidene; the double bond is
        // kept by the isothiocyanato form.
        ("CC(=O)NN=C=S", "1-(isothiocyanatoamino)ethan-1-one"),
        // P-66.3.2.1 (1): under a senior acid the group is the preferred
        // acyl prefix hydrazinecarbonyl.
        (
            "OC(=O)c1ccc(cc1)C(=O)NN=C(C)C",
            "4-[2-(propan-2-ylidene)hydrazine-1-carbonyl]benzoic acid",
        ),
        // P-68.3.1.2.5: a semicarbazone stays a hydrazinecarboxamide.
        (
            "NC(=O)NN=C(C)C",
            "2-(propan-2-ylidene)hydrazine-1-carboxamide",
        ),
    ] {
        check(smiles, expected, NameStatus::Preferred);
    }
    // No Blue Book example places the descriptor of a configured N'=C bond
    // in a hydrazide name, so the amino form is kept and not preferred.
    check(
        "CC(=O)N/N=C/c1ccccc1",
        "1-({[(E)-phenylmethylidene]amino}amino)ethan-1-one",
        NameStatus::NotPreferred(NameReason::AminoPrefix),
    );
}

#[test]
fn the_prefix_of_a_hydrazone_is_hydrazinylidene() {
    for (smiles, expected) in [
        // P-68.3.1.2.1: H2N-N= is hydrazinylidene (not hydrazono), a simple
        // prefix; only compound and complex prefixes are enclosed (P-16.5.1.1).
        ("OC(=O)CC=NN", "3-hydrazinylidenepropanoic acid"),
        ("OC(=O)C(=NN)C", "2-hydrazinylidenepropanoic acid"),
        ("OC(=O)CC(=NN)CC(O)=O", "3-hydrazinylidenepentanedioic acid"),
        (
            "OC(=O)C1CCC(CC1)=NN",
            "4-hydrazinylidenecyclohexane-1-carboxylic acid",
        ),
        ("CC(=O)CC=NN", "4-hydrazinylidenebutan-2-one"),
        ("OC(=O)C/C=N/N", "(3E)-3-hydrazinylidenepropanoic acid"),
        // P-66.4.2.3.2: 3-amino-3-hydrazinylidenepropanoic acid.
        ("OC(=O)CC(N)=NN", "3-amino-3-hydrazinylidenepropanoic acid"),
        // P-65.1.3.2.2: 5-hydrazinylidene-5-hydroxypentanoic acid.
        (
            "OC(=O)CCCC(O)=NN",
            "5-hydrazinylidene-5-hydroxypentanoic acid",
        ),
        // P-68.3.1.2.2: 4-(phenylhydrazinylidene)cyclohexane-1-carboxylic
        // acid. The substitutable hydrogens are all on nitrogen 2, so no
        // locants are cited (P-14.3.4.6).
        (
            "OC(=O)C1CCC(CC1)=NNc1ccccc1",
            "4-(phenylhydrazinylidene)cyclohexane-1-carboxylic acid",
        ),
        ("OC(=O)CC=NNC", "3-(methylhydrazinylidene)propanoic acid"),
        (
            "OC(=O)CC=NN(C)C",
            "3-(dimethylhydrazinylidene)propanoic acid",
        ),
        (
            "OC(=O)CC=NN(C)CC",
            "3-[ethyl(methyl)hydrazinylidene]propanoic acid",
        ),
        // P-16.5.1.1, P-16.5.1.3.1: as (chloromethyl)(methyl)silane.
        (
            "OC(=O)CC=NN(C)CCl",
            "3-[(chloromethyl)(methyl)hydrazinylidene]propanoic acid",
        ),
        (
            "OC(=O)CC=NN(CCl)C",
            "3-[(chloromethyl)(methyl)hydrazinylidene]propanoic acid",
        ),
    ] {
        check(smiles, expected, NameStatus::Preferred);
    }
    for (smiles, expected, reason) in [
        // P-68.3.1.2.5: carbamoylhydrazinylidene and its chalcogen analogues
        // are the preferred prefixes of semicarbazones; the label of these
        // names is left as it was.
        (
            "OC(=O)CCC(=NNC(N)=O)CCC",
            "4-(carbamoylhydrazinylidene)heptanoic acid",
            NameReason::UnretainedParent,
        ),
        (
            "OC(=O)CCC(=NNC(=O)N(C)C)CCC",
            "4-[(dimethylcarbamoyl)hydrazinylidene]heptanoic acid",
            NameReason::UnretainedParent,
        ),
        (
            "OC(=O)CCC(=NNC(=[Se])N(C)C)CCC(O)=O",
            "4-[(dimethylcarbamoselenoyl)hydrazinylidene]heptanedioic acid",
            NameReason::UnretainedParent,
        ),
        // The Blue Book shows no acyl group on hydrazinylidene.
        (
            "OC(=O)CC=NNC(C)=O",
            "3-(acetylhydrazinylidene)propanoic acid",
            NameReason::AminoPrefix,
        ),
        // P-68.3.1.2.2 (2): the functional class name of a hydrazone.
        (
            "NN=CCCC(C)=NN",
            "4-hydrazinylidenepentanal hydrazone",
            NameReason::FunctionalClass,
        ),
    ] {
        check(smiles, expected, NameStatus::NotPreferred(reason));
    }
}
