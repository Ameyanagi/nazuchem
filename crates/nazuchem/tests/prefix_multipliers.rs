// SPDX-License-Identifier: MIT OR Apache-2.0
// Origin: original
// Source: Nazuchem | crates/nazuchem/tests/prefix_multipliers.rs | issue #134

#![cfg(all(feature = "n2s", feature = "s2n"))]

use nazuchem::{Budget, name_to_smiles, s2n::canonical_smiles, smiles_to_name};

#[test]
fn simple_and_substituted_prefix_multipliers_round_trip() {
    let budget = Budget::default();
    for (smiles, preferred, alternate) in [
        // P-16.3.4(a), (b), (c): attachment locants, unsaturation locants,
        // and initial numerical terms call for a basic multiplier with parentheses.
        (
            "CC(C)c1cc(C(C)C)cc(C(C)C)c1",
            "1,3,5-tri(propan-2-yl)benzene",
            Some("1,3,5-tris(propan-2-yl)benzene"),
        ),
        (
            "CC(=C)c1ccc(C(C)=C)cc1",
            "1,4-di(prop-1-en-2-yl)benzene",
            Some("1,4-bis(prop-1-en-2-yl)benzene"),
        ),
        (
            "CCCCCCCCCCCCCc1ccc(CCCCCCCCCCCCC)cc1",
            "1,4-di(tridecyl)benzene",
            Some("1,4-bis(tridecyl)benzene"),
        ),
        (
            "O=C(O)c1cc(-c2ccc3ccccc3c2)cc(-c2ccc3ccccc3c2)c1",
            "3,5-di(naphthalen-2-yl)benzoic acid",
            Some("3,5-bis(naphthalen-2-yl)benzoic acid"),
        ),
        (
            "C/C=C/SS/C=C/C",
            "1,2-di[(1E)-prop-1-en-1-yl]disulfane",
            Some("1,2-bis[(1E)-prop-1-en-1-yl]disulfane"),
        ),
        // Indicated hydrogen locates the parent hydride in a simple prefix.
        (
            "O=C(O)c1cc(N(c2ccc[nH]2)c2ccc[nH]2)ccc1",
            "3-[di(1H-pyrrol-2-yl)amino]benzoic acid",
            Some("3-[bis(1H-pyrrol-2-yl)amino]benzoic acid"),
        ),
        // P-16.5.1.3.1 explicitly gives this silicon example.
        (
            "CC[SiH](C(C)C)C(C)C",
            "ethyldi(propan-2-yl)silane",
            Some("ethylbis(propan-2-yl)silane"),
        ),
        // P-16.5.1.4 explicitly gives this simple enclosed acyl example.
        (
            "O=C(O)C(S(=O)c1ccccc1)S(=O)c1ccccc1",
            "di(benzenesulfinyl)acetic acid",
            Some("bis(benzenesulfinyl)acetic acid"),
        ),
        (
            "O=C(O)c1cc(C(=O)NN)cc(C(=O)NN)c1",
            "3,5-di(hydrazinecarbonyl)benzoic acid",
            Some("3,5-bis(hydrazinecarbonyl)benzoic acid"),
        ),
        // P-63.2.2.2 classifies retained tert-butoxy as simple. The hyphen
        // sets off the italic descriptor as in di-tert-butyl (P-16.2.4.1(d)).
        (
            "CC(C)(C)Oc1ccc(OC(C)(C)C)cc1",
            "1,4-di-tert-butoxybenzene",
            Some("1,4-bis(tert-butoxy)benzene"),
        ),
        (
            "O=C(O)c1cc(C(C)(C)C)cc(C(C)(C)C)c1",
            "3,5-di-tert-butylbenzoic acid",
            Some("3,5-bis(tert-butyl)benzoic acid"),
        ),
        // The parent hydride in an acyl prefix can itself carry a substituent.
        (
            "O=C(O)c1cc(S(=O)CCl)cc(S(=O)CCl)c1",
            "3,5-bis(chloromethanesulfinyl)benzoic acid",
            None,
        ),
        (
            "O=C(O)c1cc(C(=O)OC)cc(C(=O)OC)c1",
            "3,5-bis(methoxycarbonyl)benzoic acid",
            None,
        ),
        // P-16.3.5(a): chloro substitutes ethyl, so use a derived multiplier.
        (
            "ClCCc1cc(C(=O)O)cc(CCCl)c1",
            "3,5-bis(2-chloroethyl)benzoic acid",
            None,
        ),
        // P-16.3.6(c), (e): these simple prefixes require derived multipliers.
        (
            "O=C(O)c1cc(-c2nocc2)cc(-c2nocc2)c1",
            "3,5-bis(1,2-oxazol-3-yl)benzoic acid",
            None,
        ),
        (
            "[O-]C(=O)c1cc([N+]#N)cc([N+]#N)c1",
            "3,5-bis(diazyn-1-ium-1-yl)benzoate",
            Some("3,5-di(diazyn-1-ium-1-yl)benzoate"),
        ),
        (
            "OC(=O)c1cc(C2=NCCO2)cc(C2=NCCO2)c1",
            "3,5-bis(4,5-dihydro-1,3-oxazol-2-yl)benzoic acid",
            None,
        ),
        ("OC(=O)C(SS)SS", "bis(disulfanyl)acetic acid", None),
        // P-16.3.6(a) explicitly extends the polynuclear ambiguity to Sb.
        ("O=C(O)C([SbH2])[SbH2]", "bis(stibanyl)acetic acid", None),
        // P-16.3.3(a): the neutral suffix is formed by functional replacement.
        (
            "SOc1ccc(OS)cc1",
            "benzene-1,4-di(OS-thioperoxol)",
            Some("benzene-1,4-bis(OS-thioperoxol)"),
        ),
    ] {
        assert_eq!(smiles_to_name(smiles).unwrap().text, preferred, "{smiles}");
        let canonical = canonical_smiles(smiles, &budget).unwrap();
        for name in std::iter::once(preferred).chain(alternate) {
            let parsed = name_to_smiles(name).unwrap_or_else(|error| panic!("{name}: {error:?}"));
            assert_eq!(
                canonical_smiles(&parsed, &budget).unwrap(),
                canonical,
                "{name}"
            );
        }
    }
}
