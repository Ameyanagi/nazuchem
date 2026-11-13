// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | crates/nazuchem/tests/dihydrazides_nitramides.rs | issue #176

#![cfg(all(feature = "n2s", feature = "s2n"))]

//! Issue #176: di- and polyhydrazides cite the parent locant on N and N',
//! N-nitro and N-nitroso amines are named as nitramides and nitrous amides,
//! and a hydrazone with a charged nitrogen is no longer named as the neutral
//! hydrazone (IUPAC 2013 recommendations). Every generated name reads back as
//! the structure it was generated from.

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
fn dihydrazides_cite_the_parent_locant_on_n_and_n_prime() {
    for (smiles, expected) in [
        // P-66.3.1.1: on an alkane the suffix in position 1 is labelled
        // N^1, N'^1 and the one in position x N^x, N'^x. P-66.3.3.2: the
        // parent locant is a superscript on N or N'.
        ("CNNC(=O)CCC(=O)NNC", "N'^1,N'^4-dimethylbutanedihydrazide"),
        ("CNNC(=O)CC(=O)NNC", "N'^1,N'^3-dimethylpropanedihydrazide"),
        ("CNNC(=O)CCC(=O)NN", "N'^1-methylbutanedihydrazide"),
        ("NNC(=O)CCC(=O)N(C)N", "N^1-methylbutanedihydrazide"),
        (
            "CN(C)NC(=O)CCC(=O)NN",
            "N'^1,N'^1-dimethylbutanedihydrazide",
        ),
        (
            "CNNC(=O)CCC(=O)N(C)NC",
            "N^1,N'^1,N'^4-trimethylbutanedihydrazide",
        ),
        // P-66.3.3.2: N^1 is lower than N'^1, so N^1,N'^4 is lower than
        // N'^1,N^4; the Blue Book examples on naphthalene.
        (
            "CN(N)C(=O)c1ccc(C(=O)NNC)c2ccccc12",
            "N^1,N'^4-dimethylnaphthalene-1,4-dicarbohydrazide",
        ),
        (
            "CN(N)C(=O)c1ccc(C(=O)NNC)c2cc(C)ccc12",
            "N^1,N'^4,6-trimethylnaphthalene-1,4-dicarbohydrazide",
        ),
        ("CNNC(=O)CCC(=O)N(C)N", "N^1,N'^4-dimethylbutanedihydrazide"),
        (
            "CNNC(=O)CC(C)CC(=O)N(C)N",
            "N^1,N'^5,3-trimethylpentanedihydrazide",
        ),
        (
            "CNNC(=O)C1CCC(C(=O)N(C)N)CC1",
            "N^1,N'^4-dimethylcyclohexane-1,4-dicarbohydrazide",
        ),
        // P-14.3.5: letter locants are lower than numerals, so N^1 decides
        // before the locant of the methyl group on the parent.
        (
            "CN(N)C(=O)CC(C)C(=O)NNC",
            "N^1,N'^4,3-trimethylbutanedihydrazide",
        ),
        (
            "CN(N)C(=O)c1ccc(C(=O)NNC)c2ccc(C)cc12",
            "N^1,N'^4,7-trimethylnaphthalene-1,4-dicarbohydrazide",
        ),
        // The same order for the N locants of polyamines (P-62.2.4.1.2).
        ("CNCCC(C)CN", "N^1,3-dimethylbutane-1,4-diamine"),
        (
            "CNc1ccc(N)c2ccc(C)cc12",
            "N^1,7-dimethylnaphthalene-1,4-diamine",
        ),
        // P-14.5: the prefix cited first takes the lower locant.
        (
            "CCNNC(=O)CCC(=O)NNC",
            "N'^1-ethyl-N'^4-methylbutanedihydrazide",
        ),
        (
            "CCNNC(=O)CCC(=O)N(C)N",
            "N'^4-ethyl-N^1-methylbutanedihydrazide",
        ),
        // P-66.3.1.2.2: two or more carbohydrazide suffixes.
        (
            "CNNC(=O)c1ccc(C(=O)NNC)cc1",
            "N'^1,N'^4-dimethylbenzene-1,4-dicarbohydrazide",
        ),
        // Ylidene dihydrazides (P-66.3.3.1 with P-68.3.1.3.5.1).
        (
            "CNNC(=O)CCC(=O)NN=C(C)C",
            "N'^1-methyl-N'^4-(propan-2-ylidene)butanedihydrazide",
        ),
        (
            "CC(C)=NNC(=O)CCC(=O)NN=C(C)C",
            "N'^1,N'^4-bis(propan-2-ylidene)butanedihydrazide",
        ),
        // P-66.3.3.3: an acyl group substitutes the senior hydrazide.
        (
            "CNNC(=O)CCC(=O)NNC(C)=O",
            "N'^1-acetyl-N'^4-methylbutanedihydrazide",
        ),
        // The sulfur suffixes and the chalcogen analogues (P-66.3.4).
        (
            "CNNS(=O)(=O)CCS(=O)(=O)NNC",
            "N'^1,N'^2-dimethylethane-1,2-disulfonohydrazide",
        ),
        (
            "CNNC(=S)CCC(=S)NNC",
            "N'^1,N'^4-dimethylbutanedithiohydrazide",
        ),
        // One hydrazide group keeps N and N' (P-66.3.3.1).
        ("CC(=O)N(C)NC", "N,N'-dimethylacetohydrazide"),
        ("NNC(=O)CCC(=O)NN", "butanedihydrazide"),
    ] {
        check(smiles, expected, NameStatus::Preferred);
    }
}

#[test]
fn amides_of_nitric_and_nitrous_acid_substitute_nitramide_and_nitrous_amide() {
    for (smiles, expected) in [
        // P-67.1.2.6.3: preferred names are based on nitramide and nitrous
        // amide, not on nitro and nitroso amines.
        ("CCN[N+](=O)[O-]", "ethylnitramide"),
        ("CCNN(=O)=O", "ethylnitramide"),
        ("N[N+](=O)[O-]", "nitramide"),
        ("NN=O", "nitrous amide"),
        ("CN(C)[N+](=O)[O-]", "dimethylnitramide"),
        ("c1ccccc1N[N+](=O)[O-]", "phenylnitramide"),
        ("C1CCCCC1N(C)[N+](=O)[O-]", "cyclohexyl(methyl)nitramide"),
        // The Blue Book examples of P-67.1.2.6.3.
        ("CCCN(N=O)CCC", "dipropylnitrous amide"),
        ("CCN(N=O)CCCC", "butyl(ethyl)nitrous amide"),
        ("ClCN(C)[N+](=O)[O-]", "(chloromethyl)(methyl)nitramide"),
        ("CN([N+](=O)[O-])[N+](=O)[O-]", "methyl(nitro)nitramide"),
        ("CN(N=O)c1ccccc1", "methyl(phenyl)nitrous amide"),
        ("CC[C@H](C)N[N+](=O)[O-]", "[(2S)-butan-2-yl]nitramide"),
        // P-41: amides (class 11) are senior to nitriles, ketones, hydroxy
        // compounds and amines, which become prefixes.
        ("N#CCCN[N+](=O)[O-]", "(2-cyanoethyl)nitramide"),
        ("OCCN[N+](=O)[O-]", "(2-hydroxyethyl)nitramide"),
        ("OCCN(N=O)CCO", "bis(2-hydroxyethyl)nitrous amide"),
        ("CC(=O)CN(C)N=O", "methyl(2-oxopropyl)nitrous amide"),
        ("NCCN[N+](=O)[O-]", "(2-aminoethyl)nitramide"),
        // A nitro group on carbon is only a prefix.
        (
            "CN(c1c(cc(cc1[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-]",
            "methyl(2,4,6-trinitrophenyl)nitramide",
        ),
        // Substitution by ylidene groups, as in N'-hexylidenenitrous
        // hydrazide; benzylidene is preferred when unsubstituted (P-29.6.2.1).
        ("CC(C)=N[N+](=O)[O-]", "(propan-2-ylidene)nitramide"),
        ("c1ccccc1C=N[N+](=O)[O-]", "benzylidenenitramide"),
        // P-41: acids, esters and amides of senior acids stay principal and
        // cite nitramido or a substituted amino prefix (P-67.1.4.3.2).
        ("OC(=O)CCN[N+](=O)[O-]", "3-nitramidopropanoic acid"),
        ("CCOC(=O)CCN[N+](=O)[O-]", "ethyl 3-nitramidopropanoate"),
        ("NC(=O)CCN[N+](=O)[O-]", "3-nitramidopropanamide"),
        ("CC(=O)N[N+](=O)[O-]", "N-nitroacetamide"),
        ("CN(N=O)C(N)=O", "N-methyl-N-nitrosourea"),
        // A ring nitrogen takes the prefix.
        ("C1CCN(CC1)[N+](=O)[O-]", "1-nitropiperidine"),
        ("O=NN1CCCCC1", "1-nitrosopiperidine"),
    ] {
        check(smiles, expected, NameStatus::Preferred);
    }
}

#[test]
fn a_charged_hydrazone_nitrogen_is_not_named_as_the_neutral_hydrazone() {
    // A hydrazone or azine name describes the neutral N-N, so these names
    // used to read back without the charge. The cation is cited by the
    // suffix iminium or aminium (P-73.1.2.1, Table 7.4). P-73.1.2 takes the
    // largest neutral parent, as in pentamethylhydrazinium (P-73.1.1.2) and
    // N',N',N'-trimethylbenzohydrazid-N'-ium (P-73.1.2.1); that parent is
    // not chosen yet, so the names are not certified.
    for (smiles, expected) in [
        (
            "CC(=O)N[N+](C)=C(C)C",
            "N-acetamido-N-methylpropan-2-iminium",
        ),
        ("CC(=O)N[NH+]=C(C)C", "N-acetamidopropan-2-iminium"),
        ("N[N+](C)=C(C)C", "N-amino-N-methylpropan-2-iminium"),
        (
            "CN[N+](C)=C(C)C",
            "N-methyl-N-(methylamino)propan-2-iminium",
        ),
        (
            "CC(C)=N[N+](C)(C)C",
            "N,N-dimethyl-N-[(propan-2-ylidene)amino]methanaminium",
        ),
        (
            "CC(C)=N[NH+]=C(C)C",
            "N-[(propan-2-ylidene)amino]propan-2-iminium",
        ),
        ("C[N+](C)(C)N", "N-amino-N,N-dimethylmethanaminium"),
        (
            "CC(C)=[N+](C)NO",
            "N-(hydroxyamino)-N-methylpropan-2-iminium",
        ),
    ] {
        check(
            smiles,
            expected,
            NameStatus::Unverified(NameReason::CationParent),
        );
    }
    // The neutral hydrazone and azine names are unchanged.
    check(
        "CC(C)=NN",
        "propan-2-one hydrazone",
        NameStatus::NotPreferred(NameReason::FunctionalClass),
    );
    check(
        "NC(N)=NN=C(C)C",
        "propan-2-one diaminomethylidenehydrazone",
        NameStatus::NotPreferred(NameReason::FunctionalClass),
    );
    // A nitroso group is an acyl group, not a hydrazine nitrogen.
    check(
        "C[N+](C)(C)N=O",
        "N,N-dimethyl-N-nitrosomethanaminium",
        NameStatus::Preferred,
    );
    // An aminium on a ring nitrogen keeps the largest neutral parent.
    check(
        "C[N+](C)(C)N1CCCCC1",
        "N,N,N-trimethylpiperidin-1-aminium",
        NameStatus::Preferred,
    );
}
