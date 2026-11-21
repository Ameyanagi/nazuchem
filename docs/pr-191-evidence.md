# Nomenclature and regression evidence for PR #191

The comparison uses base `c3b023bf1551027a93463847b5a27fc17f44976b` and tested code `2f70359adeecfceda4a3f83ec4035d1790b20e71`. It follows the [corrected electronic edition of IUPAC Recommendations and Preferred Names 2013](https://iupac.qmul.ac.uk/BlueBook/), prepared by Henri A. Favre and Warren H. Powell. The issue sections link directly to the applied clauses and identify printed examples where applicable.

The tables record observed API spellings and preference grades. `exact` means the strict parser's output passed the molecular identity check against the input. The check preserves connectivity, formal charge, isotope placement and specified stereochemistry. The thiocarboxylate rows explicitly identify the audited mesomer equivalence. A preferred grade is the API's output. The rule explanation establishes the relationship to the recommendation. Direct parser cases show the actual output SMILES or error.

Already covered base cases provide controls. Remaining cases and conflicting examples are recorded with their outcomes so the evidence can be reviewed independently of issue closure.

The [machine-readable input and outcome table](pr-191-cases.tsv) contains every replayed case, including atom-order variants, deliberate refusals and direct parser controls. An em dash in an identity field marks an outcome without a generated name to compare.

| Measured outcome | Base | Tested code |
| --- | ---: | ---: |
| Generated names with verified graph identity | 334 | 407 |
| Generated names with an identity mismatch | 25 | 0 |
| Generated names rejected by the default strict parser | 7 | 1 |
| Structure-to-name refusals | 52 | 10 |
| Successful direct strict parse probes | 53 | 64 |
| Direct strict parse refusals | 16 | 5 |
| Total comparisons | 487 | 487 |

The one remaining default parse rejection of a generated name is the standalone phosphanyl radical. Its exact read-back uses the parser’s radical option in the [radical round-trip regression](https://github.com/Ameyanagi/nazuchem/blob/311c31b6b330e5cc0589cd96f56ee558ac7f42cc/crates/nazuchem-s2n/src/radical_naming.rs#L603-L625). Direct parser refusals include ambiguity controls and unsupported constructions documented in the issue sections.

The golden gate independently covers the same 1,650 corpus structures on both revisions. Twenty-seven fixture records changed. Strict name parsing successes rose from 1,671 to 1,678, and failures fell from six to four. Of the primary names from those 27 records, the default strict parser accepted 26. A separate replay verified exact graph identity for all 26. The remaining primary name is the standalone phosphanyl radical. [Base fixture](https://github.com/Ameyanagi/nazuchem/blob/c3b023bf1551027a93463847b5a27fc17f44976b/tests/fixtures/golden/s2n_n2s.jsonl) · [Tested fixture](https://github.com/Ameyanagi/nazuchem/blob/2f70359adeecfceda4a3f83ec4035d1790b20e71/tests/fixtures/golden/s2n_n2s.jsonl) · [Golden gate](https://github.com/Ameyanagi/nazuchem/blob/2f70359adeecfceda4a3f83ec4035d1790b20e71/crates/nazuchem/tests/golden.rs#L72)

[#4](#issue-4) · [#6](#issue-6) · [#11](#issue-11) · [#15](#issue-15) · [#24](#issue-24) · [#27](#issue-27) · [#88](#issue-88) · [#89](#issue-89) · [#95](#issue-95) · [#96](#issue-96) · [#100](#issue-100) · [#110](#issue-110) · [#113](#issue-113) · [#121](#issue-121) · [#128](#issue-128) · [#134](#issue-134) · [#138](#issue-138) · [#141](#issue-141) · [#156](#issue-156) · [#160](#issue-160) · [#176](#issue-176) · [#180](#issue-180) · [#182](#issue-182) · [#186](#issue-186)

## Phanes, salts, and prefix multipliers

### Issue 4

**Linear phanes and multiplicative names**

The follow-up scope includes a 10-node phane, an amine node, a branched linker, a methanetetrayl linker and two ring-linker types. The 16-node furan case below is a separate remaining case. [P-26.4.2](https://iupac.qmul.ac.uk/BlueBook/P2.html#260402) · [P-15.3.2.1](https://iupac.qmul.ac.uk/BlueBook/P1.html#15030201)

| Case and input | Base result | Tested result |
| --- | --- | --- |
| `4-ten-node`<br>`c1ccccc1CCCc1ccc(Cc2ccc(CCc3ccccc3)cc2)cc1` | `NotPreferred(LinearPhane)` · `1-{[4-(2-phenylethyl)phenyl]methyl}-4-(3-phenylpropyl)benzene` · `exact` | `Preferred` · `1,10(1),4,6(1,4)-tetrabenzenadecaphane` · `exact` |
| `4-amine`<br>`c1ccccc1SCc1ccc(cc1)Nc1ccc(cc1)CCc1ccccc1` | `NotPreferred(LinearPhane)` · `N-[4-(2-phenylethyl)phenyl]-4-[(phenylsulfanyl)methyl]aniline` · `exact` | `Preferred` · `2-thia-5-aza-1,9(1),4,6(1,4)-tetrabenzenanonaphane` · `exact` |
| `4-branched-linker`<br>`CC(C)(Cc1ccccc1)Cc1ccccc1` | `NotPreferred(Multiplicative)` · `(2,2-dimethyl-3-phenylpropyl)benzene` · `exact` | `Preferred` · `1,1'-(2,2-dimethylpropane-1,3-diyl)dibenzene` · `exact` |
| `4-quaternary-linker`<br>`C(c1ccccc1)(c1ccccc1)(c1ccccc1)c1ccccc1` | `NotPreferred(Multiplicative)` · `(triphenylmethyl)benzene` · `exact` | `Preferred` · `1,1',1'',1'''-methanetetrayltetrabenzene` · `exact` |
| `4-naphthalene-linker`<br>`OC(=O)Cc1c(CC(=O)O)cc2ccccc2c1` | `NotPreferred(Multiplicative)` · `[3-(carboxymethyl)naphthalen-2-yl]acetic acid` · `exact` | `Preferred` · `2,2'-(naphthalene-2,3-diyl)diacetic acid` · `exact` |
| `4-pyridine-linker`<br>`OC(=O)CCCCCCCCCc1ncccc1CCCCCCCCCC(=O)O` | `NotPreferred(Multiplicative)` · `10-[2-(9-carboxynonyl)pyridin-3-yl]decanoic acid` · `exact` | `Preferred` · `10,10'-(pyridine-2,3-diyl)di(decanoic acid)` · `exact` |


Earlier PRs corrected the labels, added single-atom linkers, then added compound and chain linkers, some substituted/heterocyclic units, and unsubstituted linear phanes with 7–9 nodes. The exact baseline already emits preferred multiplicative names for the three-ring ether, the three-arm ether, and the bis(methylene) original issue controls. All three have exact identity in both builds. [Issue #4 and follow-ups](https://github.com/Ameyanagi/nazuchem/issues/4), [#55](https://github.com/Ameyanagi/nazuchem/pull/55), [#132](https://github.com/Ameyanagi/nazuchem/pull/132), [merge batch #142](https://github.com/Ameyanagi/nazuchem/pull/142).

| Case and SMILES | Base | Head |
| --- | --- | --- |
| `4-ether-phenol` · `Oc1ccc(Oc2ccc(Oc3ccc(Oc4ccccc4)cc3)cc2)cc1` | `4-[4-(4-phenoxyphenoxy)phenoxy]phenol`, `NotPreferred(LinearPhane)` | `2,4,6-trioxa-1,7(1),3,5(1,4)-tetrabenzenaheptaphan-1⁴-ol`, `Preferred` |
| `4-methylene-phenol` · `Oc1ccc(Cc2ccc(Cc3ccc(Cc4ccccc4)cc3)cc2)cc1` | `4-({4-[(4-benzylphenyl)methyl]phenyl}methyl)phenol`, `NotPreferred(LinearPhane)` | `1,7(1),3,5(1,4)-tetrabenzenaheptaphan-1⁴-ol`, `Preferred` |
| `4-sulfur-phenol` · `Oc1ccc(Oc2ccc(Sc3ccc(Oc4ccccc4)cc3)cc2)cc1` | `4-{4-[(4-phenoxyphenyl)sulfanyl]phenoxy}phenol`, `NotPreferred(LinearPhane)` | `2,6-dioxa-4-thia-1,7(1),3,5(1,4)-tetrabenzenaheptaphan-1⁴-ol`, `Preferred` |
| `4-carbonate` · `O=C(O)c1ccc(Oc2ccc(OC(=O)Oc3ccc(Oc4ccccc4)cc3)cc2)cc1` | `4-(4-carboxyphenoxy)phenyl 4-phenoxyphenyl carbonate`, `Unverified(UnauditedPath)` | `5-oxo-2,4,6,8-tetraoxa-1,9(1),3,7(1,4)-tetrabenzenanonaphane-1⁴-carboxylic acid`, `Preferred` |
| `4-biphenyl` · `O=C(O)CCCCCCCCCOc1ccc(-c2ccc(OCCCCCCCCCC(=O)O)cc2)cc1` | `10-({4'-[(9-carboxynonyl)oxy][1,1'-biphenyl]-4-yl}oxy)decanoic acid`, `NotPreferred(Multiplicative)` | `10,10'-[[1,1'-biphenyl]-4,4'-diylbis(oxy)]di(decanoic acid)`, `Preferred` |
| `4-alkene-bromo` · `Brc1ccc(C=CCc2ccc(Br)cc2)cc1` | `1-bromo-4-[3-(4-bromophenyl)prop-1-en-1-yl]benzene`, `NotPreferred(Multiplicative)` | `1,1'-(prop-1-ene-1,3-diyl)bis(4-bromobenzene)`, `Preferred` |
| `4-direction-final` · `c1ccccc1Cc1ccc(Cc2ccc(Cc3ccc(Cc4ccccc4)nc3)cn2)cn1` | `2-benzyl-5-({5-[(6-benzylpyridin-3-yl)methyl]pyridin-2-yl}methyl)pyridine`, `NotPreferred(LinearPhane)` | `3,5(2,5),7(5,2)-tripyridina-1,9(1)-dibenzenanonaphane`, `Preferred` |

- All table rows have exact strict identity on both revisions. The first three paths have four rings, two terminal rings, and seven nodes. Example 2 prints the unsuffixed trioxa parent and the smaller multiplicative ether. P-26.4.3.1 supplies the composite OH locant, P-26.5.1 assigns skeletal replacement prefixes, and P-26.6 refers suffix treatment to Phane Part II. The acid in the carbonate row is senior and the carbonate O–C–O segment becomes part of the nine-node path with an oxo prefix. [P-52.2.5.1](https://iupac.qmul.ac.uk/BlueBook/P5.html#52020501), [P-52.2.5.3](https://iupac.qmul.ac.uk/BlueBook/P5.html#52020503), [P-26.4.3.1](https://iupac.qmul.ac.uk/BlueBook/P2.html#26040301), [P-26.5.1](https://iupac.qmul.ac.uk/BlueBook/P2.html#260501), [P-26.6](https://iupac.qmul.ac.uk/BlueBook/P2.html#2606), [PhII-4](https://iupac.qmul.ac.uk/phane2/Ph4.html), [P-41](https://iupac.qmul.ac.uk/BlueBook/P4.html#41), [P-65.6.3.2.3](https://iupac.qmul.ac.uk/BlueBook/P6.html#6506030203), [phane regression fixtures](https://github.com/Ameyanagi/nazuchem/blob/311c31b6b330e5cc0589cd96f56ee558ac7f42cc/crates/nazuchem-s2n/src/phane/tests.rs#L179-L281).
- The biphenyl case is printed in P-16.5.4.1.2. Its assembly bracket is intrinsic and ignored when selecting the outer enclosure. P-15.3.3.1 prints the unsaturated linker joining the substituted bromobenzene units and permits an unsymmetrical center whose later units attach by identical bonds. P-26.4.2.4 prints the last row and selects the lower attachment locant set compared by ascending superatom locants. [P-15.3.2.1](https://iupac.qmul.ac.uk/BlueBook/P1.html#15030201), [P-15.3.2.4.1](https://iupac.qmul.ac.uk/BlueBook/P1.html#1503020401), [P-15.3.3.1](https://iupac.qmul.ac.uk/BlueBook/P1.html#15030301), [P-16.5.4.1.2](https://iupac.qmul.ac.uk/BlueBook/P1.html#1605040102), [P-26.4.2.4](https://iupac.qmul.ac.uk/BlueBook/P2.html#26040204), [multiplicative regression fixtures](https://github.com/Ameyanagi/nazuchem/blob/311c31b6b330e5cc0589cd96f56ee558ac7f42cc/crates/nazuchem-s2n/src/multiplicative.rs#L911-L960).
- Remaining coverage includes the printed P-52.2.5.3 Example 3 `2,5,7,10,12,15-hexathia-1,16(3),6,11(3,4)-tetrafuranahexadecaphane`. Both revisions give a substitutive name with `NotPreferred(LinearPhane)` and exact identity. Strict parsing of the printed PIN reports `APPEARS_AMBIGUOUS` for 12 and 15 on head because the parser also recognizes flattened composite locants. Both revisions also still return exact substitutive fallbacks for the anisole-derived and fully fluorinated identical units in the fixture. [P-52.2.5.3 Example 3](https://iupac.qmul.ac.uk/BlueBook/P5.html#52020503), [remaining multiplicative fixture cases](https://github.com/Ameyanagi/nazuchem/blob/311c31b6b330e5cc0589cd96f56ee558ac7f42cc/crates/nazuchem-s2n/src/multiplicative.rs#L1383-L1418).

### Issue 11

**Citation order across ions in a salt**

The exact base already produces `sodium phenoxide` and `disodium acetate 4-methylbenzoate` as `Preferred`, and `1,2,3,4-tetrahydro-1,4-epoxynaphthalen-6-amine` as `Unverified(BridgedFusedParent)`. These original controls are identical on head and strictly round-trip. Earlier PRs added ordering by descriptor-free citation and ties in narrower cases. [Issue #11 and follow-ups](https://github.com/Ameyanagi/nazuchem/issues/11), [#68](https://github.com/Ameyanagi/nazuchem/pull/68), [#135](https://github.com/Ameyanagi/nazuchem/pull/135), [merge batch #142](https://github.com/Ameyanagi/nazuchem/pull/142).

| Case | Base | Head |
| --- | --- | --- |
| `11-nested` · `[Na+].[Na+].O=C([O-])C12CC(CCCl)C(CC1)C2.O=C([O-])C1C(C(C)Cl)C2CCC1C2` | `disodium 3-(2-chloroethyl)bicyclo[2.2.1]heptane-1-carboxylate 3-(1-chloroethyl)bicyclo[2.2.1]heptane-2-carboxylate`, `Preferred` | `disodium 3-(1-chloroethyl)bicyclo[2.2.1]heptane-2-carboxylate 3-(2-chloroethyl)bicyclo[2.2.1]heptane-1-carboxylate`, `Preferred` |
| `11-marks` · `[Na+].[Na+].O=C([O-])C1CC([SiH2]CN)CC1.O=C([O-])C1CC([SiH](C)N)CC1` | `disodium 3-[(aminomethyl)silyl]cyclopentane-1-carboxylate 3-[amino(methyl)silyl]cyclopentane-1-carboxylate`, `Unverified(SaltIonOrder)` | `disodium 3-[amino(methyl)silyl]cyclopentane-1-carboxylate 3-[(aminomethyl)silyl]cyclopentane-1-carboxylate`, `Preferred` |
| `11-nested-isotope` · `[Na+].[Na+].O=C([O-])C12CC(CCCl)C(CC1)C2.O=C([O-])C1C(C([13CH3])Cl)C2CCC1C2` | `disodium 3-(2-chloroethyl)bicyclo[2.2.1]heptane-1-carboxylate 3-[1-chloro(2-13C)ethyl]bicyclo[2.2.1]heptane-2-carboxylate`, `Unverified(IsotopePlacement)` | `disodium 3-[1-chloro(2-13C)ethyl]bicyclo[2.2.1]heptane-2-carboxylate 3-(2-chloroethyl)bicyclo[2.2.1]heptane-1-carboxylate`, `Unverified(IsotopePlacement)` |
| `11-stereo-repeat` · `[Na+].[Na+].[Na+].[Na+].C/C=C/C(=O)[O-].C/C=C\C(=O)[O-].[O-]C(=O)/C=C/C.CC(=O)[O-]` | `tetrasodium bis[(2E)-but-2-enoate] (2Z)-but-2-enoate acetate`, `Unverified(SaltIonOrder)` | `tetrasodium acetate bis[(2E)-but-2-enoate] (2Z)-but-2-enoate`, `Unverified(SaltIonOrder)` |

- All table rows have exact strict identity on both revisions. The order for a binary salt is cation before anion. Within each sign IR-5.4.1 specifies alphabetical citation. The organic comparison in P-14.5 considers nonitalic Roman letters, italic letters when needed, and locants, with isotope and stereo descriptors omitted. P-14.5.4 ranks the lower internal locant when the corresponding earlier locants tie. This puts `1-chloroethyl` ahead of `2-chloroethyl` in the first and third rows, before the differing parent suffix locants. The isotope remains in the actual cited name. [P-77](https://iupac.qmul.ac.uk/BlueBook/P7.html#77), [P-65.6.2.1](https://iupac.qmul.ac.uk/BlueBook/P6.html#65060201), [Red Book IR-5.4.1](https://iupac.org/wp-content/uploads/2016/07/Red_Book_2005.pdf), [P-14.5](https://iupac.qmul.ac.uk/BlueBook/P1.html#1405), [P-14.5.4](https://iupac.qmul.ac.uk/BlueBook/P1.html#140504).
- P-14.6 compares enclosure and punctuation after the alphanumerical characters tie. Its printed example compares the exact two silyl prefixes of the second row on a cyclopentanol parent and ranks `a` before `(`. In the fourth row the two butenoates still tie once E/Z is omitted. The head applies the provable ranking of acetate before both, preserving the doubled E ion, the Z ion, and four sodium ions. Other unresolved pairs include distinct ring-assembly junction or von Baeyer numerals, superscript composite locants, and stereo-free rendering that regroups prefixes. All four controls keep `Unverified(SaltIonOrder)` and exact identity. [P-14.6](https://iupac.qmul.ac.uk/BlueBook/P1.html#1406), [salt regressions](https://github.com/Ameyanagi/nazuchem/blob/311c31b6b330e5cc0589cd96f56ee558ac7f42cc/crates/nazuchem-s2n/src/species.rs#L1146-L1326).

### Issue 134

**Basic and derived prefix multipliers**

The multiplied whole acyl groups use the basic multiplier with enclosing marks. The rule-derived hydrazinecarbonyl and benzenesulfonyl cases parallel the printed benzenesulfinyl example. The `bis(sulfonyl)` exception applies to the bare prefix. [P-16.5.1.4](https://iupac.qmul.ac.uk/BlueBook/P1.html#16050104) · [P-16.3.6](https://iupac.qmul.ac.uk/BlueBook/P1.html#160306)

| Case and input | Base result | Tested result |
| --- | --- | --- |
| `134-hydrazide`<br>`O=C(O)c1cc(C(=O)NN)cc(C(=O)NN)c1` | `Preferred` · `3,5-bis(hydrazinecarbonyl)benzoic acid` · `exact` | `Preferred` · `3,5-di(hydrazinecarbonyl)benzoic acid` · `exact` |
| `134-benzenesulfonyl-dual`<br>`O=C(O)C(S(=O)(=O)c1ccccc1)S(=O)(=O)c1ccccc1` | `Preferred` · `bis(benzenesulfonyl)acetic acid` · `exact` | `Preferred` · `di(benzenesulfonyl)acetic acid` · `exact` |


P-16.3.4(a)–(d),(f) use a basic multiplier with enclosure for simple locanted prefixes, simple locanted unsaturated prefixes, prefixes beginning with a numerical term or with `dec`, and simple components with brackets. P-16.3.5(a) assigns derived multipliers to substituted compound or complex prefixes. P-16.3.6 lists separate ambiguity exceptions, including mononuclear subsets of polynuclear acyclic structures, retained divalent groups, Hantzsch–Widman and skeletal replacement names, attached fusion components, and names beginning with `di` such as `diazenyl` and `disulfanyl`. The multiplier is independent of locants alone. [Issue #134](https://github.com/Ameyanagi/nazuchem/issues/134), [#146](https://github.com/Ameyanagi/nazuchem/pull/146), [P-16.3.4](https://iupac.qmul.ac.uk/BlueBook/P1.html#160304), [P-16.3.5](https://iupac.qmul.ac.uk/BlueBook/P1.html#160305), [P-16.3.6](https://iupac.qmul.ac.uk/BlueBook/P1.html#160306).

| Case and SMILES | Base | Head |
| --- | --- | --- |
| `134-propyl` · `CC(C)c1cc(C(C)C)cc(C(C)C)c1` | `1,3,5-tris(propan-2-yl)benzene` | `1,3,5-tri(propan-2-yl)benzene` |
| `134-unsaturation` · `CC(=C)c1ccc(C(C)=C)cc1` | `1,4-bis(prop-1-en-2-yl)benzene` | `1,4-di(prop-1-en-2-yl)benzene` |
| `134-tridecyl` · `CCCCCCCCCCCCCc1ccc(CCCCCCCCCCCCC)cc1` | `1,4-bis(tridecyl)benzene` | `1,4-di(tridecyl)benzene` |
| `134-sulfinyl` · `O=C(O)C(S(=O)c1ccccc1)S(=O)c1ccccc1` | `bis(benzenesulfinyl)acetic acid` | `di(benzenesulfinyl)acetic acid` |
| `134-tert-butoxy` · `CC(C)(C)Oc1ccc(OC(C)(C)C)cc1` | `1,4-bis(tert-butoxy)benzene` | `1,4-di-tert-butoxybenzene` |
| `134-peroxol` · `SOc1ccc(OS)cc1` | `benzene-1,4-bis(OS-thioperoxol)` | `benzene-1,4-di(OS-thioperoxol)` |

- Every spelling in the table is emitted as `Preferred` on its respective revision and round-trips exactly. Each old/new pair also strictly parses to identical SMILES on both revisions. P-16.5.1.3.1 prints `ethyldi(propan-2-yl)silane` and P-16.5.1.4 prints `di(benzenesulfinyl)acetic acid`, where an enclosed simple acyl prefix carries a parent hydride name. P-63.2.2.2 classifies retained tert-butoxy as a simple prefix. P-16.2.4.1(d) separates the italic tert descriptor with a hyphen. The thioperoxol suffix is a functional replacement under P-16.3.3(a). [P-16.5.1.3.1](https://iupac.qmul.ac.uk/BlueBook/P1.html#1605010301), [P-16.5.1.4](https://iupac.qmul.ac.uk/BlueBook/P1.html#16050104), [P-63.2.2.2](https://iupac.qmul.ac.uk/BlueBook/P6.html#63020202), [P-16.2.4.1](https://iupac.qmul.ac.uk/BlueBook/P1.html#16020401), [P-16.3.3](https://iupac.qmul.ac.uk/BlueBook/P1.html#160303), [public API fixtures with both spellings](https://github.com/Ameyanagi/nazuchem/blob/311c31b6b330e5cc0589cd96f56ee558ac7f42cc/crates/nazuchem/tests/prefix_multipliers.rs#L10-L129).
- Exact identity and `Preferred` hold on both revisions for `3,5-bis(2-chloroethyl)benzoic acid` under P-16.3.5(a), `3,5-bis(1,2-oxazol-3-yl)benzoic acid` under P-16.3.6(c), and `3,5-bis(diazyn-1-ium-1-yl)benzoate` under P-16.3.6(e). The latter is the charged prefix from the issue report and begins with the numerical term `di` for the two nitrogens. P-16.3.4(f) explicitly lists the bicyclo[3.2.1]octan-3-yl multiplier. The additional `134-printed-bicyclo` probe changes `3,5-bis(bicyclo[3.2.1]octan-3-yl)benzoic acid` to `3,5-di(bicyclo[3.2.1]octan-3-yl)benzoic acid`, with exact identity and `Preferred` on both revisions. The separately tested labelled benzoyl acyl fallback retains `NotPreferred(AcylPrefix)` on both revisions. [P-16.3.5](https://iupac.qmul.ac.uk/BlueBook/P1.html#160305), [P-16.3.6](https://iupac.qmul.ac.uk/BlueBook/P1.html#160306).
- The seven golden rows changed by the #134 commit were reviewed against the exact head fixture and rerun. Lines 398 and 1559 change a simple locanted propan-2-yl prefix to `di(...)`. Lines 1173, 1224, 1225, 1293 and 1294 change multiplied retained tert-butoxy to `di-tert-butoxy` or `tri-tert-butoxy`. Line 1225 also contracts nesting from brackets to parentheses after dropping the inner parentheses. Six generated names pass strict identity on both revisions and both spellings parse identically. Line 1294 is a standalone phosphanyl substituent and strict parsing declines both revisions. The raw table records the naming and parsing results under `134-golden-*`. [Raw cases](pr-191-cases.tsv) · [#134 implementation and fixture diff](https://github.com/Ameyanagi/nazuchem/commit/ddd2ab640f72c9de330a3dae985ef028072103bf) · [Golden fixture](https://github.com/Ameyanagi/nazuchem/blob/2f70359adeecfceda4a3f83ec4035d1790b20e71/tests/fixtures/golden/s2n_n2s.jsonl#L1294)
- Two independent boranyl or silyl groups on acetic acid illustrate a remaining coverage gap. Both revisions mark `diborylacetic acid` and `disilylacetic acid` as preferred. P-29.3.1 gives the preselected prefix `boranyl`, and the note in P-16.3.6(a) extends the ambiguity rule to boron and silicon. The expected `bis(boranyl)acetic acid` and `bis(silyl)acetic acid` strict-parse to the same respective structures. [P-16.3.6(a)](https://iupac.qmul.ac.uk/BlueBook/P1.html#160306), [P-29.3.1](https://iupac.qmul.ac.uk/BlueBook/P2.html#290301), [preferred prefix appendix](https://iupac.qmul.ac.uk/BlueBook/Papp2.html).

## Retained rings and fused cations

### Issue 15

**Heteranthrenes and mixed pheno...ine parents**

The six original issue examples are already preferred and exact in both revisions. They were addressed by #49, and #171 covered an earlier set of mixed chalcogen and pnictogen parents. Fifteen tested examples from the follow-up scope gain the expected preferred names at this head. Both revisions reconstruct the graphs from their own generated names. [Issue #15](https://github.com/Ameyanagi/nazuchem/issues/15) · [PR #49](https://github.com/Ameyanagi/nazuchem/pull/49) · [PR #171](https://github.com/Ameyanagi/nazuchem/pull/171)

| Case | Base → head | IUPAC 2013 rationale |
| --- | --- | --- |
| `c1ccc2c(c1)B=C1C=CC=CC1=B2` | `dibenzo[b,e][1,4]diborinine` (not preferred) → `boranthrene` (preferred, exact) | Both heteroatoms are B. The book explicitly gives the mancude parent as `boranthrene`. [P-25.2.2.2](https://iupac.qmul.ac.uk/BlueBook/P2.html#25020202) |
| `c1ccc2c(c1)Bc1ccccc1B2` | Von Baeyer (not preferred) → `5,10-dihydroboranthrene` (preferred, exact). Direct strict parsing of the latter fails at base locant 10 and succeeds at head. | Saturation introduces hydrogen at both B positions (5 and 10) with the benzenoid rings intact. [P-25.2.2.2](https://iupac.qmul.ac.uk/BlueBook/P2.html#25020202) · [P-31.2.1](https://iupac.qmul.ac.uk/BlueBook/P3.html#310201) |
| `c1ccc2c(c1)Oc1ccc3ccccc3c1[SiH2]2` | `12H-naphtho[2,1-b][1,4]benzoxasiline` (preferred) → `12H-benzo[a]phenoxasiline` (preferred, exact) | The 1,4-O/Si unit with two fused benzenes is `phenoxasiline`. It has more rings than the competing two-ring O/Si component, and the added ring occupies side `a`. Renumbering the complete system moves the indicated H to Si at position 12. [P-25.2.2.3](https://iupac.qmul.ac.uk/BlueBook/P2.html#25020203) · [P-25.3.2.4(b)](https://iupac.qmul.ac.uk/BlueBook/P2.html#25030204) · [P-25.3.3.1.1](https://iupac.qmul.ac.uk/BlueBook/P2.html#2503030101) |
| `c1ccc2c(c1)[Te]c1ccccc1[GeH2]2` | `10H-dibenzo[b,e][1,4]telluragermine` (not preferred) → `10H-phenotelluragermine` (preferred, exact) | The mixed 1,4 pair takes `pheno` and the Hantzsch–Widman name. `tellura` precedes `germa` and Ge takes the 6B stem `ine`. The senior Te gets locant 5, so the indicated H is on Ge at 10. [P-25.2.2.3](https://iupac.qmul.ac.uk/BlueBook/P2.html#25020203) · [Tables 2.4 and 2.5](https://iupac.qmul.ac.uk/BlueBook/P2.html#t204) · [P-25.3.3.1.2(b)](https://iupac.qmul.ac.uk/BlueBook/P2.html#2503030102) |
| `c1ccc2c(c1)P=C1C=CC=CC1=[As]2` | `dibenzo[b,e][1,4]phospharsinine` (not preferred) → `phenophospharsinine` (preferred, exact) | The mixed 1,4 pair takes `pheno`, with P before As and the 6C stem `inine`. The input supplies mancude P=C and As=C bonds. [P-25.2.2.3](https://iupac.qmul.ac.uk/BlueBook/P2.html#25020203) · [Tables 2.4 and 2.5](https://iupac.qmul.ac.uk/BlueBook/P2.html#t204) |
| `c1ccc2c(c1)Pc1ccc3ccccc3c1P2` | Von Baeyer (not preferred) → `7,12-dihydrobenzo[a]phosphanthrene` (preferred, exact) | Use the printed same-element parent `phosphanthrene` as the benzo-fusion component, then cite hydrogenation at the two P atoms at 7 and 12. The As, Si and B analogues show the corresponding result. [P-25.2.2.2](https://iupac.qmul.ac.uk/BlueBook/P2.html#25020202) · [P-25.3.2.1.3](https://iupac.qmul.ac.uk/BlueBook/P2.html#2503020103) · [P-31.2.1](https://iupac.qmul.ac.uk/BlueBook/P3.html#310201) |

The fused N/Si input `c1ccc2c(c1)Nc1ccc3ccccc3c1[SiH2]2` changes from a von Baeyer fallback to exact `7,12-dihydronaphtho[2,1-b][1,4]benzazasiline` with `NotPreferred(UnregisteredRetainedComponent)`. The N/Ge fused control stays on an exact not-preferred von Baeyer fallback and the N/Bi three-ring control stays on the exact not-preferred dibenzo name. The general mixed-pair rule still applies. [P-25.2.2.3](https://iupac.qmul.ac.uk/BlueBook/P2.html#25020203) · [Head regression cases](https://github.com/Ameyanagi/nazuchem/blob/311c31b6b330e5cc0589cd96f56ee558ac7f42cc/crates/nazuchem-s2n/src/lib.rs#L3257-L3346)

### Issue 24

**Derivatives of fused pheno...ine parents**

All seven original issue rows match their expected names and are preferred and exact at head. Six were already preferred and exact at base through #67. Each changed case uses the senior three-ring phenoxazine component under the fusion seniority rules. [Issue #24](https://github.com/Ameyanagi/nazuchem/issues/24) · [PR #67](https://github.com/Ameyanagi/nazuchem/pull/67) · [P-25.2.2.3](https://iupac.qmul.ac.uk/BlueBook/P2.html#25020203) · [P-25.3.2.4](https://iupac.qmul.ac.uk/BlueBook/P2.html#25030204)

| Case | Base → head | IUPAC 2013 rationale |
| --- | --- | --- |
| `O=C1CC=C2Oc3cc4Nc5ccccc5Sc4cc3[N+](C)(C)C2=C1` | Exact not-preferred von Baeyer → `14,14-dimethyl-2-oxo-2,3,7,14-tetrahydro[1,4]benzothiazino[2,3-b]phenoxazin-14-ium` (preferred, exact). The N7 isomer likewise gets the matching `7,7-dimethyl…-7-ium` name. | The dioxoisoindolium example uses prefix `oxo` on the cationic parent. Here the ketone takes `2-oxo` and the charged N takes `-ium`. [P-73.1.2.1(2)](https://iupac.qmul.ac.uk/BlueBook/P7.html#73010201) |
| `OC(=O)c1ccc2nc3c(cc(=O)c4ccccc34)oc2c1` | `5-oxonaphtho[2,1-b][1,4]benzoxazine-9-carboxylic acid` (not preferred, exact) → `5-oxo-5H-benzo[a]phenoxazine-9-carboxylic acid` (preferred, exact) | The acid supplies the suffix. The C5 substituent supplies the `oxo` prefix after the parent has been specified with `5H`. The sulfur and imine analogues give the corresponding `5-sulfanylidene-5H-` and `5-imino-5H-` names. [P-64.5.1](https://iupac.qmul.ac.uk/BlueBook/P6.html#640501) · [P-58.2.4](https://iupac.qmul.ac.uk/BlueBook/P5.html#580204) · [P-64.6.1](https://iupac.qmul.ac.uk/BlueBook/P6.html#640601) |
| `c1ccc2c(c1)Nc1cc3Sc4ccccc4[N-]c3cc1O2` | Declined → `[1,4]benzothiazino[2,3-b]phenoxazin-7(14H)-ide` (preferred, exact). The N14 isomer likewise becomes exact preferred `…-14(7H)-ide`. | The formally dihydro parent loses a hydron at N7. Its surviving hydrogen at N14 is added indicated hydrogen written directly after the anionic locant, as in printed `pyridin-1(2H)-ide`. [P-72.2.2.1.1](https://iupac.qmul.ac.uk/BlueBook/P7.html#7202020101) |
| `N=C1CC=C2Oc3cc4Nc5ccccc5Sc4cc3NC2=C1` | Exact not-preferred von Baeyer → `7,14-dihydro[1,4]benzothiazino[2,3-b]phenoxazin-2(3H)-imine` (preferred, exact) | The imine supplies the suffix at C2. Added indicated hydrogen at C3 accounts for the suffix-induced local saturation while retaining the two hydro N sites. [P-62.3.1.1](https://iupac.qmul.ac.uk/BlueBook/P6.html#62030101) · [P-58.2.2.1](https://iupac.qmul.ac.uk/BlueBook/P5.html#58020201) |

The standalone thione and the mixed N-methyl/quaternary cation controls both still give the same exact `NotPreferred(VonBaeyerFusionFallback)` names. The natural analogous thione spelling `7,14-dihydro[1,4]benzothiazino[2,3-b]phenoxazin-2(3H)-thione` already parses to the same output SMILES in both revisions. The head adds standalone imine generation and thione-prefix generation with a senior acid. [P-64.6.1](https://iupac.qmul.ac.uk/BlueBook/P6.html#640601) · [Head regression cases](https://github.com/Ameyanagi/nazuchem/blob/311c31b6b330e5cc0589cd96f56ee558ac7f42cc/crates/nazuchem-s2n/src/lib.rs#L3420-L3454)

### Issue 96

**Fused cations**

Two additional cases come from the remaining issue report. The dioxolo compound changes from a von Baeyer fallback to a fused cation name. The cyclopentapyrrolinium is a base control. [P-73.1.2.1](https://iupac.qmul.ac.uk/BlueBook/P7.html#73010201) · [Issue #96](https://github.com/Ameyanagi/nazuchem/issues/96)

| Case and input | Base result | Tested result |
| --- | --- | --- |
| `96-dioxolo-isoquinolin`<br>`C[N+]1(C)CCc2cc3OCOc3cc2C1` | `NotPreferred(VonBaeyerFusionFallback)` · `11,11-dimethyl-4,6-dioxa-11-azoniatricyclo[7.4.0.0^{3,7}]trideca-1,3(7),8-triene` · `exact` | `Preferred` · `6,6-dimethyl-5,6,7,8-tetrahydro-2H-[1,3]dioxolo[4,5-g]isoquinolin-6-ium` · `exact` |
| `96-cyclopenta-pyrrol`<br>`C1CC2CC[NH2+]C2C1` | `Preferred` · `octahydrocyclopenta[b]pyrrol-1-ium` · `exact` | `Preferred` · `octahydrocyclopenta[b]pyrrol-1-ium` · `exact` |


Base already gives the expected preferred names and exact graphs for the issue's decahydro-1,6-naphthyridinium, dimethylisoindolium and oxo-naphthyridinium examples after #154. At head, the fourth original example (quinolizidinium) becomes preferred with the same exact spelling, `octahydro-2H-quinolizin-5-ium`. [Issue #96](https://github.com/Ameyanagi/nazuchem/issues/96) · [PR #154](https://github.com/Ameyanagi/nazuchem/pull/154) · [P-25.2.1](https://iupac.qmul.ac.uk/BlueBook/P2.html#250201) · [P-73.1.1.2](https://iupac.qmul.ac.uk/BlueBook/P7.html#73010102)

| Case | Base → head | IUPAC 2013 rationale |
| --- | --- | --- |
| `C[N+]1(C)C(=O)c2ccccc2C1=O` | `2,2-dimethyl-1,3-dioxo-1,3-dihydro-2H-isoindol-2-ium` (preferred, exact) → `2,2-dimethyl-1,3-dioxo-2,3-dihydro-1H-isoindol-2-ium` (preferred, exact) | The printed PIN keeps the two oxo prefixes and gives the indicated hydrogen the lower locant (1). [P-73.1.2.1(2)](https://iupac.qmul.ac.uk/BlueBook/P7.html#73010201) · [P-31.2.2](https://iupac.qmul.ac.uk/BlueBook/P3.html#310202) |
| `O=C1[NH2+]CCC2CCCCC12` | `1-oxo-3,4,4a,5,6,7,8,8a-octahydro-2H-isoquinolin-2-ium` (unverified, exact) → `1-oxodecahydroisoquinolin-2-ium` (preferred, exact) | Complete hydrogenation of the neutral polycyclic parent takes `decahydro` without locants. C1 is oxo on the cationic parent and N2 takes `-ium`. [P-31.2.3.3.2](https://iupac.qmul.ac.uk/BlueBook/P3.html#3102030302) · [P-73.1.2.1(2)](https://iupac.qmul.ac.uk/BlueBook/P7.html#73010201) |
| `c1ccc2c(c1)CC[NH2+]2` | Declined → `2,3-dihydro-1H-indol-1-ium` (preferred, exact) | Adding a hydron to N1 preserves the neutral parent's `1H` and the hydrogenation at C2 and C3. [P-73.1.1.2](https://iupac.qmul.ac.uk/BlueBook/P7.html#73010102) · [P-31.2.2](https://iupac.qmul.ac.uk/BlueBook/P3.html#310202) |
| `C1=Nc2ccccc2[NH2+]1` | Declined → `1H-1,3-benzimidazol-1-ium` (preferred, exact) | The benzo-fused imidazole parent is cited with its N locants. The hydron adds at the indicated-H N1, and the exact graph distinguishes this from protonation at the imine N. [P-25.2.2.4](https://iupac.qmul.ac.uk/BlueBook/P2.html#25020204) · [P-73.1.1.2](https://iupac.qmul.ac.uk/BlueBook/P7.html#73010102) |

The head's `4,5-dihydroquinolizin-5-ium` is graph-exact and `NotPreferred(FusionNitrogenHydron)`, where the base's unverified `4H-quinolizin-5-ium` parsed to a mismatched graph. Likewise, `octahydropyrrolizin-4-ium` is graph-exact and not preferred, replacing a base name that strict parsing rejected. Direct strict parsing of `quinolizinium` and `5lambda^5-quinolizin-5-ylium` yields the same H-free junction-N graph in both revisions and differs from the `4,5-dihydroquinolizin-5-ium` hydron graph. The saturated `octahydroindolizin-4-ium` remains unverified despite exact identity, and the spiro N⁺ example keeps the exact not-preferred fallback. [P-73.1.1.2](https://iupac.qmul.ac.uk/BlueBook/P7.html#73010102) · [P-73.3.1](https://iupac.qmul.ac.uk/BlueBook/P7.html#730301) · [Head regression cases](https://github.com/Ameyanagi/nazuchem/blob/311c31b6b330e5cc0589cd96f56ee558ac7f42cc/crates/nazuchem-s2n/src/fusion/planner.rs#L1785-L1829)

The selected cases did not reveal a new loss of graph identity or an introduced naming regression. Each issue still has tested scope to keep open.

## Bridged and polycyclic naming

### Issue 6

**Oxides, λ heteroatoms, vowel elision and indene suffix locants**

A single amine oxide takes a functional-class oxide name. The Blue Book prints the λ⁴ fused sulfur parent as a PIN. The phosphorus heterone takes a λ⁵ parent and the ketone suffix under method (3). [P-74.2.1.2](https://iupac.qmul.ac.uk/BlueBook/P7.html#74020102) · [P-25.6](https://iupac.qmul.ac.uk/BlueBook/P2.html#2506) · [P-74.2.1.4](https://iupac.qmul.ac.uk/BlueBook/P7.html#74020104)

The multiplying prefix loses its final “a” before “one” or “amine.” For the indene dione, both numberings permit 1H. The suffix locants 1,2 then precede 2,3, and added hydrogen at 3 accommodates the second suffix. [P-16.7.1(c)](https://iupac.qmul.ac.uk/BlueBook/P1.html#160701) · [P-14.4(b,c)](https://iupac.qmul.ac.uk/BlueBook/P1.html#1404) · [P-58.2.3.1.3](https://iupac.qmul.ac.uk/BlueBook/P5.html#5802030103)



| Case and input SMILES | Base c3b023b | Head 2f70359 |
|---|---|---|
| 6-1 · `[O-][n+]1ccccc1` | `Preferred` · `pyridine 1-oxide` · `exact` | `Preferred` · `pyridine 1-oxide` · `exact` |
| 6-2 · `C[N+]1([O-])CCCCC1` | `Preferred` · `1-methylpiperidine 1-oxide` · `exact` | `Preferred` · `1-methylpiperidine 1-oxide` · `exact` |
| 6-3 · `[SH]1=CC=CSC2=C1C=CC=C2` | `NotPreferred(VonBaeyerFusionFallback)` · `2lambda^4,6-dithiabicyclo[5.4.0]undeca-1(11),2,4,7,9-pentaene` · `exact` | `Preferred` · `1lambda^4,5-benzodithiepine` · `exact` |
| 6-4 · `O=[PH]1C=CC=C1` | `Preferred` · `1H-1lambda^5-phosphol-1-one` · `exact` | `Preferred` · `1H-1lambda^5-phosphol-1-one` · `exact` |
| 6-5 · `O=P1(C)C=CC=C1` | `Preferred` · `1-methyl-1H-1lambda^5-phosphol-1-one` · `exact` | `Preferred` · `1-methyl-1H-1lambda^5-phosphol-1-one` · `exact` |
| 6-6 · `O=C1C(=O)C(=O)C(=O)CC1` | `Preferred` · `cyclohexane-1,2,3,4-tetrone` · `exact` | `Preferred` · `cyclohexane-1,2,3,4-tetrone` · `exact` |
| 6-7 · `O=C1C(=O)Cc2ccccc21` | `Preferred` · `1H-indene-1,2(3H)-dione` · `exact` | `Preferred` · `1H-indene-1,2(3H)-dione` · `exact` |
| 6-lambda-order · `S1C=CC=[SH]C2=C1C=CC=C2` | `NotPreferred(VonBaeyerFusionFallback)` · `2lambda^4,6-dithiabicyclo[5.4.0]undeca-1(11),2,4,7,9-pentaene` · `exact` | `Preferred` · `1lambda^4,5-benzodithiepine` · `exact` |
| 6-lambda-isomer · `[SH]1=CSC=CC2=C1C=CC=C2` | `NotPreferred(VonBaeyerFusionFallback)` · `2lambda^4,4-dithiabicyclo[5.4.0]undeca-1(11),2,5,7,9-pentaene` · `exact` | `Preferred` · `1lambda^4,3-benzodithiepine` · `exact` |
| 6-biphenyl · `Nc1ccc(-c2ccc(N)c(N)c2)cc1N` | `Preferred` · `[1,1'-biphenyl]-3,3',4,4'-tetramine` · `exact` | `Preferred` · `[1,1'-biphenyl]-3,3',4,4'-tetramine` · `exact` |
| 6-methyl · `Cc1ccc2c(c1)SC=CC=[SH]2` | `NotPreferred(VonBaeyerFusionFallback)` · `9-methyl-2lambda^4,6-dithiabicyclo[5.4.0]undeca-1(7),2,4,8,10-pentaene` · `exact` | `Preferred` · `7-methyl-1lambda^4,5-benzodithiepine` · `exact` |
| 6-lambda-split · `[SH]1=CC=CC2=C1SC=CC2` | `name-error` · `openclatura: fused lambda-hydrogen parent state not ported` | `NotPreferred(VonBaeyerFusionFallback)` · `2lambda^4,10-dithiabicyclo[4.4.0]deca-1(6),2,4,8-tetraene` · `exact` |

### Issue 88

**Orientation and hydrogen in bridged fused names**

The preferred drawing maximizes a row of ortho-fused rings with vertical shared bonds. Three collinear rings of the triptolide derivatives include just two successive fused rings. Clockwise peripheral numbering starts in the rightmost uppermost ring, from its most counterclockwise nonfusion atom. Applying this orientation changes all locant-dependent parts of both names, including the nine stereodescriptors. [P-25.3.2.3.3(a)](https://iupac.qmul.ac.uk/BlueBook/P2.html#2503020303) · [P-25.3.3.1.1](https://iupac.qmul.ac.uk/BlueBook/P2.html#2503030101)

Insertion of a bridge can change the maximum number of noncumulative double bonds and the need for indicated hydrogen. This explains why the indane example can lose 1H entirely. Hydro prefixes modify the resulting parent as a unit, while any needed indicated hydrogen precedes the bridge prefix. [P-25.4.3.4.1](https://iupac.qmul.ac.uk/BlueBook/P2.html#2504030401) · [P-25.7.1.3.3](https://iupac.qmul.ac.uk/BlueBook/P2.html#2507010303) · [P-31.2.2](https://iupac.qmul.ac.uk/BlueBook/P3.html#310202)



| Case and input SMILES | Base c3b023b | Head 2f70359 |
|---|---|---|
| 88-1 · `C1CC23CCC(CC2C1)CC3` | `Unverified(BridgedFusedParent)` · `hexahydro-3a,6-ethano-2,3-dihydro-1H-indene` · `parse-error: Cannot find atom to add hydrogen to (6 hydrogens requested but only 5 positions that can be hydrogenated)` | `Unverified(BridgedFusedParent)` · `octahydro-3a,6-ethanoindene` · `exact` |
| 88-2 · `CN1CCC23C4C1CC5=C2C(=CC=C5)OC3CC=C4` | `Unverified(OpsinFusionLayout)` · `3-methyl-2,3,4,4a,7,7a-hexahydro-4,12-methano-1H-[1]benzofuro[3,2-e]isoquinoline` · `exact` | `Unverified(OpsinFusionLayout)` · `3-methyl-2,3,4,4a,7,7a-hexahydro-1H-4,12-methano[1]benzofuro[3,2-e]isoquinoline` · `exact` |
| 88-3 · `CCCCCC(=O)O[C@@H]1[C@@]2(C(C)C)O[C@H]2[C@@H]2O[C@]23[C@]12O[C@H]2C[C@H]1C2=C(CC[C@@]13C)C(=O)OC2` | `Unverified(OpsinFusionLayout)` · `(1aR,2R,2aS,3aS,3bS,4aS,4bS,9bS,10aS)-4b-methyl-7-oxo-2a-(propan-2-yl)-2,2a,3a,3b,4b,5,6,9b,10,10a-decahydro-9H-trisoxireno[2',3':4b,5;2'',3'':6,7;2''',3''':8a,9]phenanthro[1,2-c]furan-2-yl hexanoate` · `parse-error: No atom with locant 9b` | `Unverified(OpsinFusionLayout)` · `(3bS,4aS,5aR,6R,6aS,7aS,7bS,8aS,8bS)-8b-methyl-1-oxo-6a-(propan-2-yl)-3b,4,4a,6,6a,7a,7b,8b,9,10-decahydro-3H-trisoxireno[2',3':4b,5;2'',3'':6,7;2''',3''':8a,9]phenanthro[1,2-c]furan-6-yl hexanoate` · `exact` |
| 88-4 · `CC(C)[C@]12O[C@H]1[C@@H]1O[C@]13[C@]1(O[C@H]1C[C@H]1C4=C(C(=O)c5ccccc5)OC(=O)C4CC[C@@]13C)[C@@H]2O` | `Unverified(OpsinFusionLayout)` · `(1aS,2R,2aR,3aS,3bS,4aS,4bS,9bR,10aS)-9-benzoyl-2-hydroxy-4b-methyl-2a-(propan-2-yl)-2,2a,3a,3b,4b,6,6a,9b,10,10a-decahydrotrisoxireno[2',3':4b,5;2'',3'':6,7;2''',3''':8a,9]phenanthro[1,2-c]furan-7(5H)-one` · `parse-error: No atom with locant 9b` | `Unverified(OpsinFusionLayout)` · `(3bR,4aS,5aS,6R,6aR,7aS,7bS,8aS,8bS)-3-benzoyl-6-hydroxy-8b-methyl-6a-(propan-2-yl)-4,4a,6,6a,7a,7b,8b,9,10,10a-decahydrotrisoxireno[2',3':4b,5;2'',3'':6,7;2''',3''':8a,9]phenanthro[1,2-c]furan-1(3bH)-one` · `exact` |

### Issue 95

**Indicated hydrogen in neutral, charged and suffix-bearing fused parents**

The status snapshot for the isoindolobenzodiazepinone records the 12H citation. The historical worker snapshot and bakeoff corpus still contain older name observations for that input. The regressions above report the measured names from the two code revisions.


For a parent hydride, assign the lowest eligible nonfusion peripheral locant to indicated hydrogen before hydro locants. Indicated hydrogen that accommodates a suffix follows the suffix site's eligibility, as in the printed 3aH diol below. Establish the hydrogen state before applying an N-methyl prefix. Proton addition or removal modifies the parent hydride. For a single ketone and a single available indicated hydrogen, place that hydrogen at the ketone-bearing peripheral atom. [P-58.2.1.2](https://iupac.qmul.ac.uk/BlueBook/P5.html#58020102) · [P-31.2.2](https://iupac.qmul.ac.uk/BlueBook/P3.html#310202) · [P-58.2.4](https://iupac.qmul.ac.uk/BlueBook/P5.html#580204) · [P-73.1.1.2](https://iupac.qmul.ac.uk/BlueBook/P7.html#73010102) · [P-72.2.2.1](https://iupac.qmul.ac.uk/BlueBook/P7.html#72020201) · [P-58.2.3.1.1](https://iupac.qmul.ac.uk/BlueBook/P5.html#5802030101)



Five printed reference PINs exercise hydrogen accommodation and component orientation. Read printed rows 1–2 under [P-58.2.3.1.1](https://iupac.qmul.ac.uk/BlueBook/P5.html#5802030101), row 5 under [P-58.2.3.1.2](https://iupac.qmul.ac.uk/BlueBook/P5.html#5802030102), and rows 3–4 under [P-58.2.3.1.4](https://iupac.qmul.ac.uk/BlueBook/P5.html#5802030104). The seven original issue examples already pass on the base.

| Case and input SMILES | Base c3b023b | Head 2f70359 |
|---|---|---|
| 95-1 · `CN1CCc2sccc21` | `Preferred` · `4-methyl-5,6-dihydro-4H-thieno[3,2-b]pyrrole` · `exact` | `Preferred` · `4-methyl-5,6-dihydro-4H-thieno[3,2-b]pyrrole` · `exact` |
| 95-2 · `C[NH+]1Cc2ncccc2C1` | `Preferred` · `6-methyl-6,7-dihydro-5H-pyrrolo[3,4-b]pyridin-6-ium` · `exact` | `Preferred` · `6-methyl-6,7-dihydro-5H-pyrrolo[3,4-b]pyridin-6-ium` · `exact` |
| 95-3 · `C1NCc2ncccc21` | `Preferred` · `6,7-dihydro-5H-pyrrolo[3,4-b]pyridine` · `exact` | `Preferred` · `6,7-dihydro-5H-pyrrolo[3,4-b]pyridine` · `exact` |
| 95-4 · `C1C[NH2+]Cc2c1c1ccccc1[nH]2` | `Preferred` · `2,3,4,9-tetrahydro-1H-pyrido[3,4-b]indol-2-ium` · `exact` | `Preferred` · `2,3,4,9-tetrahydro-1H-pyrido[3,4-b]indol-2-ium` · `exact` |
| 95-5 · `C1C[NH2+]Cc2c1[nH]c1ccccc12` | `Preferred` · `2,3,4,5-tetrahydro-1H-pyrido[4,3-b]indol-2-ium` · `exact` | `Preferred` · `2,3,4,5-tetrahydro-1H-pyrido[4,3-b]indol-2-ium` · `exact` |
| 95-6 · `O=c1ccoc2ncccc12` | `Preferred` · `4H-pyrano[2,3-b]pyridin-4-one` · `exact` | `Preferred` · `4H-pyrano[2,3-b]pyridin-4-one` · `exact` |
| 95-7 · `c1ccc2c(c1)[n-]c1ccncc12` | `Preferred` · `5H-pyrido[4,3-b]indol-5-ide` · `exact` | `Preferred` · `5H-pyrido[4,3-b]indol-5-ide` · `exact` |
| 95-printed1 · `C1C2CC3C(C2(C=C1)O)C(CC3)O` | `Preferred` · `2,3,3a,3b,6,6a,7,7a-octahydro-1H-cyclopenta[a]pentalene-3,3b-diol` · `exact` | `Preferred` · `1,3b,4,5,6,6a,7,7a-octahydro-3aH-cyclopenta[a]pentalene-3a,4-diol` · `exact` |
| 95-printed2 · `O1C2=C(C(CC1)=O)C=C1C(OCCC1=O)=C2` | `Preferred` · `2,3,7,8-tetrahydro-4H,6H-benzo[1,2-b:4,5-e']dipyran-4,6-dione` · `exact` | `Preferred` · `2,3,7,8-tetrahydro-4H,6H-benzo[1,2-b:5,4-b']dipyran-4,6-dione` · `exact` |
| 95-printed3 · `C1=C2C(COC1=O)=CC(C=1OCCCC12)=O` | `NotPreferred(VonBaeyerFusionFallback)` · `5,11-dioxatricyclo[8.4.0.0^{2,7}]tetradeca-1(10),2,7-triene-4,9-dione` · `exact` | `Preferred` · `9,10-dihydro-2H,4H-benzo[1,2-b:4,3-c']dipyran-2,6(8H)-dione` · `exact` |
| 95-printed4 · `C1=C2C(COC1)C(C(C1OC=CC=C12)=O)=O` | `Preferred` · `4,4a-dihydro-2H,5H-benzo[1,2-c:4,3-b']dipyran-5,6(6aH)-dione` · `exact` | `Preferred` · `4,4a-dihydro-2H,5H-benzo[1,2-b:4,3-c']dipyran-5,6(6aH)-dione` · `exact` |
| 95-printed5 · `O1C2=C(C=CC1)C=C1C(OCCC1=O)=C2` | `Preferred` · `7,8-dihydro-2H,6H-benzo[1,2-b:4,5-e']dipyran-6-one` · `exact` | `Preferred` · `7,8-dihydro-2H,6H-benzo[1,2-b:5,4-b']dipyran-6-one` · `exact` |

### Issue 100

**Graph-preserving fallbacks, printed dependent bridges and long series**

The two nonhexagonal systems and four localized π inputs use explicit von Baeyer numbering and double-bond locants. The generated names are marked nonpreferred. The bridged imide is also nonpreferred and reports the acyl-prefix reason. Its reversed traversal preserves the (15R,19R)/(Z) spelling. The hydrazone inversion changes only Z to E, and the bridge epimer changes only 19R to 19S. The Blue Book describes fused-system numbering, bridged-fused parents and selection among alternative ring nomenclatures. It supplies no identity-comparison convention for the four inputs. [P-25.3.3.1.1](https://iupac.qmul.ac.uk/BlueBook/P2.html#2503030101) · [P-25.4.3.4.1](https://iupac.qmul.ac.uk/BlueBook/P2.html#2504030401) · [P-52.2.5.2.1](https://iupac.qmul.ac.uk/BlueBook/P5.html#5202050201)

The Blue Book prints each dependent-bridge example as a PIN. In the tetracycle the longer independent bridge wins. In the pentacycle the independently numbered bridge precedes the dependent one, then the ordered endpoint locants break the tie at 8 versus 13. The software conservatively returns Unverified(UnauditedParent) for each. The long-series checks exercise decacene and undecaphene, both with and without 1,2-dihydro. [P-23.2.6.2.2](https://iupac.qmul.ac.uk/BlueBook/P2.html#2302060202) · [P-23.2.6.3](https://iupac.qmul.ac.uk/BlueBook/P2.html#23020603) · [P-25.1.2.1](https://iupac.qmul.ac.uk/BlueBook/P2.html#25010201) · [P-25.1.2.2](https://iupac.qmul.ac.uk/BlueBook/P2.html#25010202) · [P-31.2.2](https://iupac.qmul.ac.uk/BlueBook/P3.html#310202)

| Case and input SMILES | Base c3b023b | Head 2f70359 |
|---|---|---|
| 100-1 · `CC1=NC=C(S1)C#CC2=C3CN4C(=NN=C4C5=C(N3C=N2)C=CC(=C5)Cl)COC` | `name-error` · `openclatura: nonhexagonal fusion parser numbering not ported` | `NotPreferred(VonBaeyerFusionFallback)` · `5-{2-[15-chloro-9-(methoxymethyl)-2,4,8,10,11-pentaazatetracyclo[11.4.0.0^{2,6}.0^{8,12}]heptadeca-1(13),3,5,9,11,14,16-heptaen-5-yl]ethynyl}-2-methyl-1,3-thiazole` · `exact` |
| 100-1-order · `COCC1=NN=C2C3=C(C=CC(Cl)=C3)N3C=NC(C#CC4=CN=C(C)S4)=C3CN21` | `name-error` · `openclatura: nonhexagonal fusion parser numbering not ported` | `NotPreferred(VonBaeyerFusionFallback)` · `5-{2-[15-chloro-9-(methoxymethyl)-2,4,8,10,11-pentaazatetracyclo[11.4.0.0^{2,6}.0^{8,12}]heptadeca-1(13),3,5,9,11,14,16-heptaen-5-yl]ethynyl}-2-methyl-1,3-thiazole` · `exact` |
| 100-2 · `C1CCCC2=C(CC1)C3=C(S2)N=CC(=NC3=O)C4=CC=CC=C4` | `name-error` · `openclatura: nonhexagonal fusion parser numbering not ported` | `NotPreferred(VonBaeyerFusionFallback)` · `5-phenyl-9-thia-4,7-diazatricyclo[8.6.0.0^{2,8}]hexadeca-1(10),2(8),4,6-tetraen-3-one` · `exact` |
| 100-2-order · `C1=CC=CC=C1C1=NC(=O)C2=C(N=C1)SC1=C2CCCCCC1` | `name-error` · `openclatura: nonhexagonal fusion parser numbering not ported` | `NotPreferred(VonBaeyerFusionFallback)` · `5-phenyl-9-thia-4,7-diazatricyclo[8.6.0.0^{2,8}]hexadeca-1(10),2(8),4,6-tetraen-3-one` · `exact` |
| 100-3 · `N1=C2C3=C(N=NC3=N1)C=C2` | `name-error` · `openclatura: localized nonaromatic fused pi state not ported` | `NotPreferred(VonBaeyerFusionFallback)` · `2,3,5,6-tetraazatricyclo[5.2.1.0^{4,10}]deca-1(10),2,4,6,8-pentaene` · `exact` |
| 100-3-order · `C1=CC2=C3C(=NN=C13)N=N2` | `name-error` · `openclatura: localized nonaromatic fused pi state not ported` | `NotPreferred(VonBaeyerFusionFallback)` · `2,3,5,6-tetraazatricyclo[5.2.1.0^{4,10}]deca-1(10),2,4,6,8-pentaene` · `exact` |
| 100-4 · `N1=C2C3=C(N=NC3=N1)N=N2` | `name-error` · `openclatura: localized nonaromatic fused pi state not ported` | `NotPreferred(VonBaeyerFusionFallback)` · `2,3,5,6,8,9-hexaazatricyclo[5.2.1.0^{4,10}]deca-1(9),2,4(10),5,7-pentaene` · `exact` |
| 100-4-order · `N1=NC2=C3C(=NN=C13)N=N2` | `name-error` · `openclatura: localized nonaromatic fused pi state not ported` | `NotPreferred(VonBaeyerFusionFallback)` · `2,3,5,6,8,9-hexaazatricyclo[5.2.1.0^{4,10}]deca-1(9),2,4(10),5,7-pentaene` · `exact` |
| 100-5 · `CC1=CC=2C3=C1N=NC3=NN2` | `name-error` · `openclatura: localized nonaromatic fused pi state not ported` | `NotPreferred(VonBaeyerFusionFallback)` · `9-methyl-2,3,5,6-tetraazatricyclo[5.2.1.0^{4,10}]deca-1(10),2,4,6,8-pentaene` · `exact` |
| 100-5-order · `N1=C2C3=C(N=NC3=N1)C(C)=C2` | `name-error` · `openclatura: localized nonaromatic fused pi state not ported` | `NotPreferred(VonBaeyerFusionFallback)` · `9-methyl-2,3,5,6-tetraazatricyclo[5.2.1.0^{4,10}]deca-1(10),2,4,6,8-pentaene` · `exact` |
| 100-6 · `FC1=CC=2C3=C1N=NC3=NN2` | `name-error` · `openclatura: localized nonaromatic fused pi state not ported` | `NotPreferred(VonBaeyerFusionFallback)` · `9-fluoro-2,3,5,6-tetraazatricyclo[5.2.1.0^{4,10}]deca-1(10),2,4,6,8-pentaene` · `exact` |
| 100-6-order · `N1=C2C3=C(N=NC3=N1)C(F)=C2` | `name-error` · `openclatura: localized nonaromatic fused pi state not ported` | `NotPreferred(VonBaeyerFusionFallback)` · `9-fluoro-2,3,5,6-tetraazatricyclo[5.2.1.0^{4,10}]deca-1(10),2,4,6,8-pentaene` · `exact` |
| 100-7 · `C1=CC=C(C=C1)CN2C(=O)[C@H]3[C@@H](C2=O)C4(C5=CC=CC=C5C3C6=CC=CC=C64)/C=N\NC(=O)C7=CC=CO7` | `name-error` · `openclatura: polycyclic parent has no audited descriptor` | `NotPreferred(AcylPrefix)` · `(15R,19R)-17-benzyl-1-[(Z)-{[(furan-2-yl)carbonyl]hydrazinylidene}methyl]-17-azapentacyclo[6.6.5.0^{2,7}.0^{9,14}.0^{15,19}]nonadeca-2,4,6,9,11,13-hexaene-16,18-dione` · `exact` |
| 100-7-order · `O1C=CC=C1C(=O)N/N=C\C12C3=CC=CC=C3C(C3=CC=CC=C31)[C@@H]1[C@H]2C(=O)N(CC2=CC=CC=C2)C1=O` | `name-error` · `openclatura: polycyclic parent has no audited descriptor` | `NotPreferred(AcylPrefix)` · `(15R,19R)-17-benzyl-1-[(Z)-{[(furan-2-yl)carbonyl]hydrazinylidene}methyl]-17-azapentacyclo[6.6.5.0^{2,7}.0^{9,14}.0^{15,19}]nonadeca-2,4,6,9,11,13-hexaene-16,18-dione` · `exact` |
| 100-7-E · `C1=CC=C(C=C1)CN2C(=O)[C@H]3[C@@H](C2=O)C4(C5=CC=CC=C5C3C6=CC=CC=C64)/C=N/NC(=O)C7=CC=CO7` | `name-error` · `openclatura: polycyclic parent has no audited descriptor` | `NotPreferred(AcylPrefix)` · `(15R,19R)-17-benzyl-1-[(E)-{[(furan-2-yl)carbonyl]hydrazinylidene}methyl]-17-azapentacyclo[6.6.5.0^{2,7}.0^{9,14}.0^{15,19}]nonadeca-2,4,6,9,11,13-hexaene-16,18-dione` · `exact` |
| 100-7-epimer · `C1=CC=C(C=C1)CN2C(=O)[C@@H]3[C@@H](C2=O)C4(C5=CC=CC=C5C3C6=CC=CC=C64)/C=N\NC(=O)C7=CC=CO7` | `name-error` · `openclatura: polycyclic parent has no audited descriptor` | `NotPreferred(AcylPrefix)` · `(15R,19S)-17-benzyl-1-[(Z)-{[(furan-2-yl)carbonyl]hydrazinylidene}methyl]-17-azapentacyclo[6.6.5.0^{2,7}.0^{9,14}.0^{15,19}]nonadeca-2,4,6,9,11,13-hexaene-16,18-dione` · `exact` |
| 100-bridge1 · `C12C3CCCCCCC(C(CCCCCC1)CCCCCC2)C2CCC(C3)C2` | `name-error` · `openclatura: polycyclic parent has no audited descriptor` | `Unverified(UnauditedParent)` · `tetracyclo[8.6.6.5^{2,9}.1^{23,26}]octacosane` · `exact` |
| 100-bridge2 · `C12CC3CCCCC4CCCCC(CC(CCC5CC5CC1)CCCC2)CC(C4)C3` | `name-error` · `openclatura: polycyclic parent has no audited descriptor` | `Unverified(UnauditedParent)` · `pentacyclo[13.7.4.3^{3,8}.0^{18,20}.1^{13,28}]triacontane` · `exact` |
| 100-decacene · `C1=CC=CC2=CC3=CC4=CC5=CC6=CC7=CC8=CC9=CC%10=CC=CC=C%10C=C9C=C8C=C7C=C6C=C5C=C4C=C3C=C12` | `Preferred` · `decacene` · `exact` | `Preferred` · `decacene` · `exact` |
| 100-1,2-dihydrodecacene · `C1CC=CC2=CC3=CC4=CC5=CC6=CC7=CC8=CC9=CC%10=CC=CC=C%10C=C9C=C8C=C7C=C6C=C5C=C4C=C3C=C12` | `Preferred` · `decacene` · `mismatch: The source structure differs in connectivity, charge, isotope, tautomer or stereochemistry. No name was assigned` | `Preferred` · `1,2-dihydrodecacene` · `exact` |
| 100-undecaphene · `C1=CC=CC2=CC3=CC4=CC5=CC6=CC=C7C=C8C=C9C=C%10C=C%11C=CC=CC%11=CC%10=CC9=CC8=CC7=C6C=C5C=C4C=C3C=C12` | `Preferred` · `undecaphene` · `exact` | `Preferred` · `undecaphene` · `exact` |
| 100-1,2-dihydroundecaphene · `C1CC=CC2=CC3=CC4=CC5=CC6=CC=C7C=C8C=C9C=C%10C=C%11C=CC=CC%11=CC%10=CC9=CC8=CC7=C6C=C5C=C4C=C3C=C12` | `Preferred` · `undecaphene` · `mismatch: The source structure differs in connectivity, charge, isotope, tautomer or stereochemistry. No name was assigned` | `Preferred` · `1,2-dihydroundecaphene` · `exact` |

## Isotope placement and parent selection

### Issue 27

**isotope placement and rare declines**

The original checklist reached this batch through three rounds of fixes. The first covered common cations and radicals, the second simple isotope prefixes, labelled radicals and iodanuide, and the third characteristic-group prefixes, ester alkyls, suffix N-substituents and the iminide. The issue's final comment also identifies thienopyrimidinones, the nitrile-bearing trithio derivative, isotope placement inside more complex names, and two kinds of isotope ties. [Issue and full history](https://github.com/Ameyanagi/nazuchem/issues/27) · [#78](https://github.com/Ameyanagi/nazuchem/pull/78) · [#130](https://github.com/Ameyanagi/nazuchem/pull/130) · [#172](https://github.com/Ameyanagi/nazuchem/pull/172)

For hydrazines the descriptor belongs immediately before the modified unit. The methyl carbon and either backbone nitrogen can be specified independently. Where an entire-name descriptor is used, an element locant can identify a hydrogen bound at sulfur. The graph check must distinguish N1 from N2 even with the same mass label and formula. For competing parents, compare all modified positions, then positions of higher atomic number, then higher mass number. Compare equivalent substituents by their attachment locants using the distinct P-45.4 criteria. [P-82.2.1](https://iupac.qmul.ac.uk/BlueBook/P8.html#820201) · [P-82.5.2](https://iupac.qmul.ac.uk/BlueBook/P8.html#820502) · [P-44.4.1.11.4–.6](https://iupac.qmul.ac.uk/BlueBook/P4.html#4404011104) · [P-45.4](https://iupac.qmul.ac.uk/BlueBook/P4.html#4504)

### Issue 121

**tied parents**

The two diphosphanyl examples are printed as PINs, allowing direct comparison with the program's ASCII lambda notation. In the nitrile, the chain with the directly attached modified diphosphanyl substituent wins. In the acid both chains have one, and the locant 1 for the higher bonding number precedes 2. For the ring-enol pair, an earlier criterion compares ene locants before the number of substituent prefixes. For the shared-atom aminoalcohol, count substituents cited directly on the parent. [P-45.3.1](https://iupac.qmul.ac.uk/BlueBook/P4.html#450301) · [P-45.3.2](https://iupac.qmul.ac.uk/BlueBook/P4.html#450302) · [P-44.4.1.10](https://iupac.qmul.ac.uk/BlueBook/P4.html#44040110) · [P-45.2.1](https://iupac.qmul.ac.uk/BlueBook/P4.html#450201)

P-45.4.2 prints the methoxyethylamine with oxygen-18 in the substituent attached at N. P-45.4.3 assigns the carbon-14-bearing substituent to position 2 of heptan-1-ol. The latter printed example uses square brackets for isotope descriptors. The tested explicit isotopologues use the parentheses of P-82.2.1. P-45.2.2 example (5) prints the 3-amino-2-methyl-3-oxopropyl form. Full multicentre CIP parent ranking was deferred by #161. [P-45.4.2](https://iupac.qmul.ac.uk/BlueBook/P4.html#450402) · [P-45.4.3](https://iupac.qmul.ac.uk/BlueBook/P4.html#450403) · [P-45.2.2](https://iupac.qmul.ac.uk/BlueBook/P4.html#450202) · [#118](https://github.com/Ameyanagi/nazuchem/pull/118) · [#161](https://github.com/Ameyanagi/nazuchem/pull/161) · [Issue](https://github.com/Ameyanagi/nazuchem/issues/121)

### Issue 156

**all four rows from ReShiki**

Row 9 (porphine) compares two compound locants with three in the older form. Row 10 ([2.2]paracyclophane) compares one compound locant with four. P-31.1.4.2(1) minimizes the number of compound locants first. Only ties reach the comparison of numerical locants in (2) and (3). This answers the suggested lower-locant alternatives for the von Baeyer names. [P-31.1.4.2](https://iupac.qmul.ac.uk/BlueBook/P3.html#31010402) · [Source table](https://github.com/Ameyanagi/ReShiki-evidence/blob/bae15056ae7a07a8395aba58a2f12808f16d04f3/pr-reviews/pr1167/name-diff.md) · [#77](https://github.com/Ameyanagi/nazuchem/pull/77)

Row 82 requires enclosing the substituted benzamido prefix, then nesting the larger units in brackets and braces. [P-16.5.1.1](https://iupac.qmul.ac.uk/BlueBook/P1.html#16050101) · [P-16.5.4](https://iupac.qmul.ac.uk/BlueBook/P1.html#160504) · [#58](https://github.com/Ameyanagi/nazuchem/pull/58)

For row 105 the fused parent has one indicated H available for two ketone groups. The rules put indicated H at the lowest permissible nonfusion peripheral site, then use added indicated H for the outstanding principal group. The printed cyclopenta[a]naphthalenedione example makes precisely the analogous choice. This supports the 1H-…-1,4(11H)-dione arrangement without authenticating the full substituted name as a PIN. [P-14.7.2](https://iupac.qmul.ac.uk/BlueBook/P1.html#140702) · [P-58.2.3.1.3](https://iupac.qmul.ac.uk/BlueBook/P5.html#5802030103) · [#125](https://github.com/Ameyanagi/nazuchem/pull/125)

#### Measured cases for issue 27

| Case and input SMILES | Base c3b023b | Head 2f70359 |
|---|---|---|
| 27-hydrazine-N1<br>`C[15NH]N` | `openclatura: an isotope label on a characteristic group or on this kind of substituent prefix is not supported` (name-error) | `1-methyl(1-15N)hydrazine` (Unverified(IsotopePlacement), exact) |
| 27-hydrazine-N2<br>`CN[15NH2]` | `openclatura: an isotope label on a characteristic group or on this kind of substituent prefix is not supported` (name-error) | `1-methyl(2-15N)hydrazine` (Unverified(IsotopePlacement), exact) |
| 27-hydrazine-dual<br>`[13CH3][15NH]N` | `openclatura: an isotope label on a characteristic group or on this kind of substituent prefix is not supported` (name-error) | `1-(13C)methyl(1-15N)hydrazine` (Unverified(IsotopePlacement), exact) |
| 27-hydrazine-sulfanyl<br>`NNS[2H]` | `(S-2H1)-1-sulfanylhydrazine` (Unverified(IsotopePlacement), exact) | `(S-2H1)-1-sulfanylhydrazine` (Unverified(IsotopePlacement), exact) |
| 27-thieno1<br>`CCCNC1=NC(=O)C2=C(N1)CSC2` | `2-(propylamino)-1,5-dihydrothieno[3,4-d]pyrimidin-4(7H)-one` (Preferred, exact) | `2-(propylamino)-5,7-dihydrothieno[3,4-d]pyrimidin-4(1H)-one` (Preferred, exact) |
| 27-thieno2<br>`CC1=NC2=C(C3=C(S2)CCC3)C(=O)N1C[C@H](COCC4=CC=CC=C4Cl)O` | `3-{(2R)-3-[(2-chlorophenyl)methoxy]-2-hydroxypropyl}-2-methyl-3,5,6,7-tetrahydro-4H-cyclopenta[4,5]thieno[2,3-d]pyrimidin-4-one` (Preferred, exact) | `3-{(2R)-3-[(2-chlorophenyl)methoxy]-2-hydroxypropyl}-2-methyl-3,5,6,7-tetrahydro-4H-cyclopenta[4,5]thieno[2,3-d]pyrimidin-4-one` (Preferred, exact) |
| 27-trithio<br>`Cc1ccc(C(=O)SC(=S)SCC#N)cc1` | `openclatura: atom 8 (C) not named` (name-error) | `S-{[(cyanomethyl)sulfanyl]carbonothioyl} 4-methylbenzene-1-carbothioate` (Unverified(FunctionalClass), exact) |
| 27-tbutyl1<br>`CC(C)([13CH3])c1ccccc1C(C)(C)C` | `1-(tert-butyl)-2-[2-methyl(1-13C)propan-2-yl]benzene` (Unverified(IsotopePlacement), exact) | `2-(tert-butyl)-1-[2-methyl(1-13C)propan-2-yl]benzene` (Unverified(IsotopePlacement), exact) |
| 27-tbutyl2<br>`CC(C)(C)c1ccccc1C(C)(C)[13CH3]` | `1-(tert-butyl)-2-[2-methyl(1-13C)propan-2-yl]benzene` (Unverified(IsotopePlacement), exact) | `2-(tert-butyl)-1-[2-methyl(1-13C)propan-2-yl]benzene` (Unverified(IsotopePlacement), exact) |
| 27-isoparent<br>`OC(=O)CC(C[13CH2]C)C(C)C[13CH3]` | `4-methyl-3-(2-13C)propyl(6-13C)hexanoic acid` (Unverified(IsotopePlacement), exact) | `3-[(4-13C)butan-2-yl](5-13C)hexanoic acid` (Unverified(IsotopePlacement), exact) |
| 27-isoparent-Z<br>`OC(=O)CC(CC([2H])[13CH3])C[13CH2]C[2H]` | `3-(2-13C,3-2H1)propyl(6-13C,5-2H1)hexanoic acid` (Unverified(IsotopePlacement), exact) | `3-(3-13C,2-2H1)propyl(5-13C,6-2H1)hexanoic acid` (Unverified(IsotopePlacement), exact) |
| 27-isoparent-A<br>`OC(=O)CC(C[14CH2][13CH3])C[13CH2][14CH3]` | `3-(2-13C,3-14C)propyl(6-13C,5-14C)hexanoic acid` (Unverified(IsotopePlacement), exact) | `3-(2-13C,3-14C)propyl(6-13C,5-14C)hexanoic acid` (Unverified(IsotopePlacement), exact) |
| 27-control<br>`[11CH3]Oc1ccccc1` | `(11C)methoxybenzene` (Unverified(IsotopePlacement), exact) | `(11C)methoxybenzene` (Unverified(IsotopePlacement), exact) |
| 27-hydroxy<br>`[2H]Oc1ccccc1C(=O)O` | `openclatura: an isotope label on a characteristic group or on this kind of substituent prefix is not supported` (name-error) | `2-(2H1)hydroxybenzoic acid` (Unverified(IsotopePlacement), exact) |
| 27-amino<br>`C[15NH]c1ccccc1C(=O)O` | `openclatura: an isotope label on a characteristic group or on this kind of substituent prefix is not supported` (name-error) | `2-[methyl(15N)amino]benzoic acid` (Unverified(IsotopePlacement), exact) |
| 27-carbonyl-rem<br>`Cl[13C](=O)c1ccccc1C(=O)O` | `openclatura: an isotope label on a characteristic group or on this kind of substituent prefix is not supported` (name-error) | `openclatura: an isotope label on a characteristic group or on this kind of substituent prefix is not supported` (name-error) |
| 27-suffix-rem<br>`OC(=O)C1(CCCCC1)[13C](=O)O` | `openclatura: an isotope label on a characteristic group or on this kind of substituent prefix is not supported` (name-error) | `openclatura: an isotope label on a characteristic group or on this kind of substituent prefix is not supported` (name-error) |
| 27-isoparent-A-swapped<br>`OC(=O)CC(C[13CH2][14CH3])C[14CH2][13CH3]` | `3-(3-13C,2-14C)propyl(5-13C,6-14C)hexanoic acid` (Unverified(IsotopePlacement), exact) | `3-(2-13C,3-14C)propyl(6-13C,5-14C)hexanoic acid` (Unverified(IsotopePlacement), exact) |

#### Measured cases for issue 121

| Case and input SMILES | Base c3b023b | Head 2f70359 |
|---|---|---|
| 121-lambda-nitrile-a<br>`P([PH4])CCC(C#N)CCPP` | `openclatura: terminal P atom 1 with 4 hydrogens and bonding number 5 is not supported` (name-error) | `4-(2lambda^5-diphosphan-1-yl)-2-(2-diphosphanylethyl)butanenitrile` (Unverified(UnauditedParent), exact) |
| 121-lambda-nitrile-b<br>`PPCCC(C#N)CCP[PH4]` | `openclatura: terminal P atom 10 with 4 hydrogens and bonding number 5 is not supported` (name-error) | `4-(2lambda^5-diphosphan-1-yl)-2-(2-diphosphanylethyl)butanenitrile` (Unverified(UnauditedParent), exact) |
| 121-lambda-acid-a<br>`[PH3](P)CCCC(C(=O)O)CCCP[PH4]` | `openclatura: terminal P atom 13 with 4 hydrogens and bonding number 5 is not supported` (name-error) | `5-(1lambda^5-diphosphan-1-yl)-2-[3-(2lambda^5-diphosphan-1-yl)propyl]pentanoic acid` (Unverified(UnauditedParent), exact) |
| 121-lambda-acid-b<br>`[PH4]PCCCC(C(=O)O)CCC[PH3]P` | `openclatura: terminal P atom 0 with 4 hydrogens and bonding number 5 is not supported` (name-error) | `5-(1lambda^5-diphosphan-1-yl)-2-[3-(2lambda^5-diphosphan-1-yl)propyl]pentanoic acid` (Unverified(UnauditedParent), exact) |
| 121-ring<br>`OC1C=CC(CC1)OC1CC(O)C(C)C=C1` | `5-[(4-hydroxycyclohex-2-en-1-yl)oxy]-2-methylcyclohex-3-en-1-ol` (Preferred, exact) | `4-[(5-hydroxy-4-methylcyclohex-2-en-1-yl)oxy]cyclohex-2-en-1-ol` (Preferred, exact) |
| 121-amino<br>`CN(C)C1(O)CCC(OC2CCC(O)(N)C(C)C2)CC1` | `1-amino-4-{[4-(dimethylamino)-4-hydroxycyclohexyl]oxy}-2-methylcyclohexan-1-ol` (Preferred, exact) | `1-amino-4-{[4-(dimethylamino)-4-hydroxycyclohexyl]oxy}-2-methylcyclohexan-1-ol` (Preferred, exact) |
| 121-carbamoyl<br>`CNC(=O)CCc1c(C)ccc(CC(C)C(N)=O)c1` | `3-[5-(3-amino-2-methyl-3-oxopropyl)-2-methylphenyl]-N-methylpropanamide` (Preferred, exact) | `3-[5-(3-amino-2-methyl-3-oxopropyl)-2-methylphenyl]-N-methylpropanamide` (Preferred, exact) |
| 121-config-control<br>`CCC(N[C@H](C)[C@@H](C)O)(N[C@@H](C)[C@@H](C)O)CCC` | `(2R,3R)-3-[(3-{[(2S,3R)-3-hydroxybutan-2-yl]amino}hexan-3-yl)amino]butan-2-ol` (Preferred, exact) | `(2R,3R)-3-[(3-{[(2S,3R)-3-hydroxybutan-2-yl]amino}hexan-3-yl)amino]butan-2-ol` (Preferred, exact) |
| 121-CIP-rem<br>`CCC(N[C@H](C)[C@@H](Cl)[C@@H](C)O)(N[C@@H](C)[C@@H](Cl)[C@@H](C)O)CCC` | `(2R,3R,4R)-3-chloro-4-[(3-{[(2S,3R,4R)-3-chloro-4-hydroxypentan-2-yl]amino}hexan-3-yl)amino]pentan-2-ol` (Unverified(ParentConfiguration), exact) | `(2R,3R,4R)-3-chloro-4-[(3-{[(2S,3R,4R)-3-chloro-4-hydroxypentan-2-yl]amino}hexan-3-yl)amino]pentan-2-ol` (Unverified(ParentConfiguration), exact) |
| 121-isoprefix-Z-a<br>`[13CH3]OCCNCC[18O]C` | `2-(13C)methoxy-N-[2-(18O)methoxyethyl]ethan-1-amine` (Unverified(IsotopePlacement), exact) | `2-(13C)methoxy-N-[2-(18O)methoxyethyl]ethan-1-amine` (Unverified(IsotopePlacement), exact) |
| 121-isoprefix-Z-b<br>`C[18O]CCNCCO[13CH3]` | `2-(18O)methoxy-N-[2-(13C)methoxyethyl]ethan-1-amine` (Unverified(IsotopePlacement), exact) | `2-(13C)methoxy-N-[2-(18O)methoxyethyl]ethan-1-amine` (Unverified(IsotopePlacement), exact) |
| 121-isoprefix-A-a<br>`C([13CH3])C(CC(CO)CC(CCC)C[14CH3])CCC` | `4-(2-13C)ethyl-2-[2-(2-14C)ethylpentyl]heptan-1-ol` (Unverified(IsotopePlacement), exact) | `4-(2-13C)ethyl-2-[2-(2-14C)ethylpentyl]heptan-1-ol` (Unverified(IsotopePlacement), exact) |
| 121-isoprefix-A-b<br>`C([14CH3])C(CC(CO)CC(CCC)C[13CH3])CCC` | `4-(2-14C)ethyl-2-[2-(2-13C)ethylpentyl]heptan-1-ol` (Unverified(IsotopePlacement), exact) | `4-(2-13C)ethyl-2-[2-(2-14C)ethylpentyl]heptan-1-ol` (Unverified(IsotopePlacement), exact) |

#### Measured cases for issue 156

| Case and input SMILES | Base c3b023b | Head 2f70359 |
|---|---|---|
| 156-para-other<br>`C1Cc2ccc(cc2)CCc2ccc1cc2` | `tricyclo[8.2.2.2^{4,7}]hexadeca-1(12),4,6,10,13,15-hexaene` (NotPreferred(Phane), exact) | `tricyclo[8.2.2.2^{4,7}]hexadeca-1(12),4,6,10,13,15-hexaene` (NotPreferred(Phane), exact) |
| 156-para-other2<br>`C1CC2=CC=C(C=C2)CCC2=CC=C1C=C2` | `tricyclo[8.2.2.2^{4,7}]hexadeca-1(12),4,6,10,13,15-hexaene` (NotPreferred(Phane), exact) | `tricyclo[8.2.2.2^{4,7}]hexadeca-1(12),4,6,10,13,15-hexaene` (NotPreferred(Phane), exact) |
| 156-porph-other<br>`C1=Cc2cc3ccc(cc4ccc(cc5nc(cc1n2)C=C5)[nH]4)[nH]3` | `21,22,23,24-tetraazapentacyclo[16.2.1.1^{3,6}.1^{8,11}.1^{13,16}]tetracosa-1(21),2,4,6(24),7,9,11,13,15,17,19-undecaene` (Unverified(UnauditedParent), exact) | `21,22,23,24-tetraazapentacyclo[16.2.1.1^{3,6}.1^{8,11}.1^{13,16}]tetracosa-1(21),2,4,6(24),7,9,11,13,15,17,19-undecaene` (Unverified(UnauditedParent), exact) |
| 156-row9<br>`C1=CC2=CC3=CC=C(N3)C=C4C=CC(=N4)C=C5C=CC(=N5)C=C1N2` | `21,22,23,24-tetraazapentacyclo[16.2.1.1^{3,6}.1^{8,11}.1^{13,16}]tetracosa-1(21),2,4,6(24),7,9,11,13,15,17,19-undecaene` (Unverified(UnauditedParent), exact) | `21,22,23,24-tetraazapentacyclo[16.2.1.1^{3,6}.1^{8,11}.1^{13,16}]tetracosa-1(21),2,4,6(24),7,9,11,13,15,17,19-undecaene` (Unverified(UnauditedParent), exact) |
| 156-row10<br>`C1CC2=CC=C(CCC3=CC=C1C=C3)C=C2` | `tricyclo[8.2.2.2^{4,7}]hexadeca-1(12),4,6,10,13,15-hexaene` (NotPreferred(Phane), exact) | `tricyclo[8.2.2.2^{4,7}]hexadeca-1(12),4,6,10,13,15-hexaene` (NotPreferred(Phane), exact) |
| 156-row82<br>`CCCOC1=CC=C(C=C1)C(=O)NC(=S)NC2=[NH+]C=C(C=C2)Br` | `5-bromo-2-[(4-propoxybenzamidocarbonothioyl)amino]pyridin-1-ium` (Preferred, exact) | `5-bromo-2-{[(4-propoxybenzamido)carbonothioyl]amino}pyridin-1-ium` (Preferred, exact) |
| 156-row105<br>`CCc1cccc2c1[nH]c1c3c(c(C(C)=O)cc12)C(=O)C=CC3=O` | `5-acetyl-10-ethyl-1H-benzo[a]carbazole-1,4(11H)-dione` (Unverified(OpsinFusionLayout), exact) | `5-acetyl-10-ethyl-1H-benzo[a]carbazole-1,4(11H)-dione` (Unverified(OpsinFusionLayout), exact) |
| 156-row82-ethoxy<br>`CCOC1=CC=C(C=C1)C(=O)NC(=S)NC2=[NH+]C=C(C=C2)Br` | `5-bromo-2-[(4-ethoxybenzamidocarbonothioyl)amino]pyridin-1-ium` (Preferred, exact) | `5-bromo-2-{[(4-ethoxybenzamido)carbonothioyl]amino}pyridin-1-ium` (Preferred, exact) |
| 156-row82-methoxy<br>`COC1=CC=C(C=C1)C(=O)NC(=S)NC2=[NH+]C=C(C=C2)Br` | `5-bromo-2-[(4-methoxybenzamidocarbonothioyl)amino]pyridin-1-ium` (Preferred, exact) | `5-bromo-2-{[(4-methoxybenzamido)carbonothioyl]amino}pyridin-1-ium` (Preferred, exact) |

## Configuration and pnictogen nomenclature

### Issue 89

[Measured cases for this group](#measured-cases-for-issues-89-128-and-160).

**odd cumulenes**

P-93.4.2.3 applies E/Z to a planar cumulene with an odd number of consecutive double bonds and distinguishable groups at each end. Compare the CIP-senior terminal ligands (same side Z, opposite sides E) and cite the beginning of the cumulated sequence. The book prints (2Z)-hexa-2,3,4-triene. For the issue's fluorinated pair the Z-propenyl arm belongs in the parent under the tied-parent rule, leaving an E-propenyl substituent. Inverting only the F-end slash changes 1Z to 1E and preserves 5Z. [P-93.4.2.3](https://iupac.qmul.ac.uk/BlueBook/P9.html#93040203) · [P-44.4.1.12.1](https://iupac.qmul.ac.uk/BlueBook/P4.html#4404011201)

The public regression distinguishes E, Z and unspecified identity, rejects contradictory alkene directions, and checks reverse traversals. The SMILES naming and strict name-to-SMILES paths carry the descriptor. Molecule IR v1 export raises a serialization error. Unknown ξ and priorities depending on another cumulenyl ligand's E/Z are outside this coverage. Open PR #147 overlaps. [Tests](https://github.com/Ameyanagi/nazuchem/blob/311c31b/crates/nazuchem/tests/cumulene_ez.rs#L183) · [PR #147](https://github.com/Ameyanagi/nazuchem/pull/147)

### Issue 128

[Measured cases for this group](#measured-cases-for-issues-89-128-and-160).

**equivalent amine oxide sites, guanidines and ureas**

The issue comment after PR #140 already records acid/suffix rank corrections, the constitutional threefold symmetry case and the nitrile/hydrazone ties. This PR addresses the remaining stereo-dependent site equivalence. Two representations of the same oxide now give the same locant while retaining 3s in the paired ligands or both E double bonds in the macrocycle. Controls interchange the oxygen between stereochemically inequivalent nitrogens and require different names and different strict identities. The macrocyclic parent carries the program's unaudited-parent status. [Issue update](https://github.com/Ameyanagi/nazuchem/issues/128#issuecomment-6842278088) · [P-62.5](https://iupac.qmul.ac.uk/BlueBook/P6.html#6205) · [Stereo regressions](https://github.com/Ameyanagi/nazuchem/blob/311c31b/crates/nazuchem-s2n/src/amine_oxide_tests.rs#L658)

The guanidine convention assigns N″ to the imino N of a known tautomer; N and N′ distinguish the amino sites. In mixed amides, compare the corresponding acids. The sulfonamide or carbothioamide is senior to the urea, and the shared N bears carbamoyl. The Blue Book explicitly prints N-carbamoylbenzenesulfonamide. [P-66.4.1.2.1.1–2](https://iupac.qmul.ac.uk/BlueBook/P6a.html#660401020101) · [P-66.1.6.1.1.5](https://iupac.qmul.ac.uk/BlueBook/P6a.html#660106010105) · [Regressions](https://github.com/Ameyanagi/nazuchem/blob/311c31b/crates/nazuchem-s2n/src/lib.rs#L1451)

### Issue 160

[Measured cases for this group](#measured-cases-for-issues-89-128-and-160).

**pnictogen acyl groups and H-bearing heterones**

The follow-up generator probes cover additive hydroxy attachment, substitution by chlorine or sulfur and the H-bearing antimony acyl group. The generated names and strict graph results are recorded below. [P-67.1.4.1.1.2–.5](https://iupac.qmul.ac.uk/BlueBook/P6a.html#670104010102)

| Case and input | Base result | Tested result |
| --- | --- | --- |
| `160-hydroxy`<br>`OC(=O)CP(C)(=O)O` | `NotPreferred(AcylPrefix)` · `[hydroxy(methyl)phosphoryl]acetic acid` · `exact` | `Preferred` · `[hydroxy(methylphosphonoyl)]acetic acid` · `exact` |
| `160-ph`<br>`OC(=O)C[PH](=O)O` | `NotPreferred(AcylPrefix)` · `(hydroxyphosphinoyl)acetic acid` · `exact` | `Preferred` · `(hydroxyphosphonoyl)acetic acid` · `exact` |
| `160-chloro`<br>`OC(=O)CP(=O)(O)Cl` | `NotPreferred(AcylPrefix)` · `[chloro(hydroxy)phosphoryl]acetic acid` · `exact` | `Preferred` · `(hydroxyphosphorochloridoyl)acetic acid` · `exact` |
| `160-stib`<br>`OC(=O)C[SbH2]=O` | `name-error` · `openclatura: stibinoyl substituent role not ported` | `Preferred` · `stibinoylacetic acid` · `exact` |
| `160-thio`<br>`OC(=O)CP(=S)(O)O` | `Preferred` · `[dihydroxy(sulfanylidene)phosphanyl]acetic acid` · `exact` | `Preferred` · `(dihydroxyphosphorothioyl)acetic acid` · `exact` |


PR #165 already delivered phosphono/phosphinoyl. For the remaining cases, hydroxy adds to a free valence of methylphosphonoyl, and halide replacement changes the acyl stem. P-67.1.4.1.1.2 assigns stibinoyl to SbH₂(O)–. The mixed P/As example is a printed PIN with phosphonic acid as the senior parent. H-bearing neutral phosphanones use the heterone suffix and the λ⁵ bonding number. Carboxylic-acid examples cover the dioxo-λ⁵ and phosphorothioyl spellings. Configured P can still take a nonpreferred acyl prefix and isotopic heterones requiring further locants can be refused. These examples establish particular bonding numbers and connectivities, without extending to all P/As/Sb valence or charge states. [Issue update](https://github.com/Ameyanagi/nazuchem/issues/160#issuecomment-6851691738) · [Additive operation](https://iupac.qmul.ac.uk/BlueBook/P6a.html#670104010105) · [Halide replacement](https://iupac.qmul.ac.uk/BlueBook/P6a.html#670104010104) · [Stibinoyl](https://iupac.qmul.ac.uk/BlueBook/P6a.html#670104010102) · [Substitutive dioxo group](https://iupac.qmul.ac.uk/BlueBook/P6a.html#670104010106) · [Mixed P/As example](https://iupac.qmul.ac.uk/BlueBook/P6a.html#67010501) · [Heterones](https://iupac.qmul.ac.uk/BlueBook/P6.html#640401) · [Suffix mode](https://iupac.qmul.ac.uk/BlueBook/P6a.html#6803020301) · [Regressions](https://github.com/Ameyanagi/nazuchem/blob/311c31b/crates/nazuchem/tests/phosphorus_acid_prefixes.rs#L98)

#### Measured cases for issues 89, 128, and 160

| Input | Base | PR head |
|---|---|---|
| `C/C=C/C(=C=C=C/F)/C=C\C` | `name-error · Naming does not support cumulene E/Z configuration yet; no less specific name was substituted` | `Preferred · (1Z,5Z)-1-fluoro-4-[(1E)-prop-1-en-1-yl]hepta-1,2,3,5-tetraene · exact` |
| `C/C=C/C(=C=C=C\F)/C=C\C` | `name-error · Naming does not support cumulene E/Z configuration yet; no less specific name was substituted` | `Preferred · (1E,5Z)-1-fluoro-4-[(1E)-prop-1-en-1-yl]hepta-1,2,3,5-tetraene · exact` |
| `(2Z)-hexa-2,3,4-triene` | `parse-error · Could not find bond that: <stereoChemistry locant="2" type="EorZ" value="Z">2Z</stereoChemistry> was referring to` | `parse · C\C=C=C=C/C` |
| `N1([C@H]([C@H](F)C)[C@@H](F)C)CC[N+]([O-])([C@H]([C@H](F)C)[C@@H](F)C)CC1` | `Unverified(AmineChalcogenide) · 1,4-bis[(2R,3s,4S)-2,4-difluoropentan-3-yl]piperazine 4-oxide · exact` | `Preferred · 1,4-bis[(2R,3s,4S)-2,4-difluoropentan-3-yl]piperazine 1-oxide · exact` |
| `CN1CC/C=C/CC[N+](C)([O-])CC/C=C/CC1` | `Unverified(UnauditedParent) · (4E,11E)-1,8-dimethyl-1,8-diazacyclotetradeca-4,11-diene 8-oxide · exact` | `Unverified(UnauditedParent) · (4E,11E)-1,8-dimethyl-1,8-diazacyclotetradeca-4,11-diene 1-oxide · exact` |
| `CN=C(N)N` | `Preferred · 2-methylguanidine · exact` | `Preferred · N''-methylguanidine · exact` |
| `NC(=O)NS(=O)(=O)c1ccccc1` | `NotPreferred(UnretainedParent) · 1-benzenesulfonamidoformamide · exact` | `Preferred · N-carbamoylbenzenesulfonamide · exact` |
| `NC(=O)NC(=S)C` | `NotPreferred(AminoPrefix) · 1-[(methylcarbonothioyl)amino]formamide · exact` | `Preferred · N-carbamoylethanethioamide · exact` |
| `[hydroxy(methylphosphonoyl)]acetic acid` | `parse · OCP(=O)=CC(=O)O` | `parse · OP(=O)(C)CC(=O)O` |
| `stibinoylacetic acid` | `parse · [Sb](=O)CC(=O)O` | `parse · [SbH2](=O)CC(=O)O` |
| `C[PH2]=O` | `NotPreferred(AcylPrefix) · phosphinoylmethane · exact` | `Preferred · methyl-lambda5-phosphanone · exact` |
| `OP(=O)(O)CCCC[As](=O)(O)O` | `Preferred · 1-arsono-4-phosphonobutane · exact` | `Preferred · (4-arsonobutyl)phosphonic acid · exact` |
| `OC(=O)CP(=O)=O` | `Preferred · dioxophosphanylacetic acid · exact` | `Preferred · (dioxo-lambda^5-phosphanyl)acetic acid · exact` |
| `OC(=O)C[P@](=O)(O)C` | `NotPreferred(AcylPrefix) · [(R)-hydroxy(methyl)phosphoryl]acetic acid · parse-error: Could not find atom that: <stereoChemistry type="RorS" value="R" stereoGroup="Abs">R</stereoChemistry> appeared to be referring to` | `NotPreferred(AcylPrefix) · [(R)-hydroxy(methyl)phosphoryl]acetic acid · exact` |

## Anhydrides and acid derivatives

### Issue 110

**Di- and polyanhydrides**

Every generated name in the following table passed the strict round-trip and structure check on both revisions. The base issue-body example already gives `acetic methanedisulfonic propanoic dianhydride` with `Unverified(UnauditedParent)` on both revisions. [Issue history](https://github.com/Ameyanagi/nazuchem/issues/110#issuecomment-6842279045)

| Input SMILES | Base `c3b023b` | Head `2f70359` | Rule applied |
| --- | --- | --- | --- |
| `CC(=O)OC(=O)C(C)CC(=O)OC(=O)CC` | `acetic 2-methyl-4-oxo-4-(propanoyloxy)butanoic anhydride` · Unverified | `1-acetic 2-methylbutanedioic 4-propanoic dianhydride` · Preferred | The unlike acids are located on the substituted dibasic acid. Acetic is on C-1. The base strict parser rejected the locanted name. [P-65.7.6.2](https://iupac.qmul.ac.uk/BlueBook/P6.html#65070602) |
| `CC(=O)OC(=O)CC(C)C(=O)OC(=O)CC` | `acetic 3-[(propanoyloxy)carbonyl]butanoic anhydride` · Unverified | `4-acetic 2-methylbutanedioic 1-propanoic dianhydride` · Preferred | Switching the unlike ends switches their locants. The base strict parser rejected this locanted name too. Both revisions reject the ambiguous unlocanted name. [P-65.7.6.2](https://iupac.qmul.ac.uk/BlueBook/P6.html#65070602) |
| `CC(=O)OC(=O)CCC(=O)OC(=O)CC(C)C(=O)OC(=O)CC` | `acetic 4-oxo-4-({3-[(propanoyloxy)carbonyl]butanoyl}oxy)butanoic anhydride` · Unverified | `4-[4-(acetyloxy)-4-oxobutanoic] 2-methylbutanedioic 1-propanoic dianhydride` · Preferred | Method (1) expresses the third linkage within one adjoining acid's substitutive name. Locants 4 and 1 specify the attachments to 2-methylbutanedioic acid. This is the printed PIN, which the base parser rejected. [P-65.7.6.2](https://iupac.qmul.ac.uk/BlueBook/P6.html#65070602) |
| `CCC(=O)OC(=O)c1cc2ccc(C(=O)OC(C)=O)cc2cc1C(=O)OC(C)=O` | `3,6-bis[(acetyloxy)carbonyl]naphthalene-2-carboxylic propanoic anhydride` · Unverified | `3,6-diacetic 2-propanoic naphthalene-2,3,6-tricarboxylic trianhydride` · Preferred | Cite the three monobasic residues alphabetically with their parent positions and the class multiplier `tri`. This is the printed PIN. [P-65.7.6.3](https://iupac.qmul.ac.uk/BlueBook/P6.html#65070603) |
| `CCC(=O)OP(=O)(OC(C)=O)CCP(=O)(OC(CC)=O)OC(C)=O` | `1-{[(acetyloxy){2-[(acetyloxy)(propanoyloxy)phosphoryl]ethyl}phosphoryl]oxy}propan-1-one` · NotPreferred(AcylPrefix) | `P,P'-diacetic P,P'-dipropanoic (ethane-1,2-diyl)bis(phosphonic) tetraanhydride` · Unverified(UnauditedParent) | Each phosphorus carries one residue of each end acid. The head emits the corrected published PIN spelling and verifies the structure. The base strict parser rejected the P locants. The geminal control `110-phosphorus-geminal` has distinct `P,P` and `P',P'` groups and passes on the head. [P-65.7.6.3](https://iupac.qmul.ac.uk/BlueBook/P6.html#65070603) |
| `CC(=O)OC(=S)CCC(=S)OC(C)=O` | `acetic 4-(acetyloxy)-4-sulfanylidenebutanethioic anhydride` · Unverified | `diacetic butanebis(thioic) dianhydride` · Preferred | The C=S groups modify the central acid, and both bridges remain oxygen. Both revisions could already parse the printed PIN. [P-65.7.6.4.3](https://iupac.qmul.ac.uk/BlueBook/P6.html#6507060403) |
| `CCC(=O)OOC(=O)CCC(=O)OOC(C)=O` | `4-(acetylperoxy)-4-oxobutanoic propanoic peroxyanhydride` · Unverified | `acetic butanedioic propanoic bis(peroxyanhydride)` · Unverified(UnauditedParent) | `peroxyanhydride` is established for a single bridge. `bis(peroxyanhydride)` applies the rule for repeating identical chalcogen linkages by analogy. This is verified structural coverage. [P-65.7.4](https://iupac.qmul.ac.uk/BlueBook/P6.html#650704), [P-65.7.6.4.1](https://iupac.qmul.ac.uk/BlueBook/P6.html#6507060401) |

For `110-higher-linear`, the head generates `acetic butanedioic 4-oxo-4-(propanoyloxy)butanoic dianhydride` as Preferred. P-65.7.6.2 prints `acetic butanedioic 4-(propanoyloxy)-4-oxobutanoic dianhydride` as the PIN. General substitutive citation rules P-14.5 and P-15.1.7.1.5 place `oxo` before `propanoyloxy`, in agreement with the head. The corrected online chapter and PDF retain the printed example. The base already parsed the printed form. The head also parses the rule's method (2) alternative `acetic dibutanedioic propanoic trianhydride` to the same graph. Treat the ordering as a conflict between the worked example and general rules. The phosphoric trianhydride example has the same oxo ordering discrepancy and receives Unverified. [P-65.7.6.2](https://iupac.qmul.ac.uk/BlueBook/P6.html#65070602), [P-65.7.6.3](https://iupac.qmul.ac.uk/BlueBook/P6.html#65070603), [P-14.5](https://iupac.qmul.ac.uk/BlueBook/P1.html#1405), [P-15.1.7.1.5](https://iupac.qmul.ac.uk/BlueBook/P1.html#1501070105), [Chapter PDF](https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf)

The head restores the printed 2-locant in `diacetic 2-methylbutanedioic dianhydride` with `Unverified(UnauditedPath)`. The originally reported S(=O)(=S) end already names and round-trips as Preferred `acetic methanesulfonothioic anhydride` on the base. For the tested sulfur-stereogenic end both revisions emit the same exact graph through the alternative `1-{[(S)-methanesulfinyl]oxy}ethan-1-one`. The isotope-labelled peroxy bridge is declined on both. [P-65.7.6.1](https://iupac.qmul.ac.uk/BlueBook/P6.html#65070601), [Polyanhydride tests](https://github.com/Ameyanagi/nazuchem/blob/311c31b6b330e5cc0589cd96f56ee558ac7f42cc/crates/nazuchem/tests/polyanhydrides.rs#L190)

### Issue 138

**Thiocarboxylates, mixed peroxy acids, and ester/halide citation**

The follow-up audit isolated a heteroatom-bound thioate whose old name reads back with oxygen bearing the negative charge. The tested generator declines the sulfur-localized drawings outside the audited mesomer family. The oxygen-localized control remains graph-exact. This test checks connectivity and charge location. [P-72.2.2.2.1.1](https://iupac.qmul.ac.uk/BlueBook/P7.html#720202020101)

| Case and input | Base result | Tested result |
| --- | --- | --- |
| `138-thioate-hetero`<br>`NC(=O)[S-]` | `Unverified(AnionParent)` · `aminomethanethioate` · `mismatch: The source structure differs in connectivity, charge, isotope, tautomer or stereochemistry. No name was assigned` | `name-error` · `openclatura: heteroatom-bound thioate charge on sulfur is not supported` |
| `138-thioate-hetero-oxygen`<br>`NC(=S)[O-]` | `Unverified(AnionParent)` · `aminomethanethioate` · `exact` | `Unverified(AnionParent)` · `aminomethanethioate` · `exact` |
| `138-thioate-phosphorus-sulfur`<br>`PC(=O)[S-]` | `Unverified(AnionParent)` · `phosphanylmethanethioate` · `mismatch: The source structure differs in connectivity, charge, isotope, tautomer or stereochemistry. No name was assigned` | `name-error` · `openclatura: heteroatom-bound thioate charge on sulfur is not supported` |
| `138-thioate-silicon-sulfur`<br>`[SiH3]C(=O)[S-]` | `Unverified(AnionParent)` · `silylmethanethioate` · `mismatch: The source structure differs in connectivity, charge, isotope, tautomer or stereochemistry. No name was assigned` | `name-error` · `openclatura: heteroatom-bound thioate charge on sulfur is not supported` |


| Input SMILES | Base `c3b023b` | Head `2f70359` and rule |
| --- | --- | --- |
| `CC(=O)[S-]` | `ethanethioate` · Unverified(AnionParent), strict parser returns O⁻ form and base verifier rejects the mismatch | `ethanethioate` · Preferred, strict parse returns `C(C)([O-])=S` as before, and the head verifier accepts the S⁻ input as mesomer-equivalent. P-72.2.2.2.1.1 explicitly prints both depictions with one PIN. The O⁻ input `CC(=S)[O-]` is exact on both revisions and moves from Unverified to Preferred. The sodium salt likewise gains Preferred and mesomer-equivalent verification. [P-72.2.2.2.1.1](https://iupac.qmul.ac.uk/BlueBook/P7.html#720202020101) |
| `OSC(=O)OS` | `(hydroxysulfanyl)(sulfanyloxy)methanone` · Unverified(UnauditedParent), exact. The printed PIN failed strict parsing. | `carbonobis(thioperoxoic) OS,SO-acid` · Preferred, exact. `bis` counts two replacement groups. `OS,SO` identifies their different two-atom sequences from the carbonic centre. [P-65.2.1.2](https://iupac.qmul.ac.uk/BlueBook/P6.html#65020102) |
| `O=C(O)OS` | `hydroxy(sulfanyloxy)methanone` · Unverified(UnauditedParent), exact. The printed PIN failed strict parsing. | `carbono(thioperoxoic) OS-acid` · Preferred, exact. One hydroxy group is retained and the modified arm runs C–O–S. This is a printed example. [P-65.2.1.2](https://iupac.qmul.ac.uk/BlueBook/P6.html#65020102) |
| `O=C(O)SO` | `hydroxy(hydroxysulfanyl)methanone` · Unverified(UnauditedParent), exact. The new name failed strict parsing. | `carbono(thioperoxoic) SO-acid` · Preferred, exact. The same rule gives the opposite element sequence for the distinct C–S–O arm. [P-65.2.1.2](https://iupac.qmul.ac.uk/BlueBook/P6.html#65020102) |
| `CCOC(=O)CCC(=O)O` | `3-(ethoxycarbonyl)propanoic acid` · Preferred, exact | `4-ethoxy-4-oxobutanoic acid` · Preferred, exact. The carbonyl carbon stays in the chain and its two ligands are cited separately. The analogous sodium salt now names the carboxylate substitutively with the ester as a prefix. [P-65.6.3.2.3](https://iupac.qmul.ac.uk/BlueBook/P6.html#6506030203), [P-65.6.3.3.5](https://iupac.qmul.ac.uk/BlueBook/P6.html#6506030305) |
| `CCOC(=O)CC(=O)Oc1ccc(OC(=O)CCC(=O)OC)cc1` | `4-[(3-ethoxy-1,3-dioxopropyl)oxy]phenyl methyl butanedioate` · Preferred, exact | `4-[(3-ethoxy-3-oxopropanoyl)oxy]phenyl methyl butanedioate` · Preferred, exact. The acyl attachment uses `propanoyl` and the other ester end supplies the 3-ethoxy and 3-oxo prefixes. The Blue Book repeats this PIN under both subsections. [P-65.6.3.3.2.1](https://iupac.qmul.ac.uk/BlueBook/P6.html#650603030201), [P-65.6.3.3.4.2](https://iupac.qmul.ac.uk/BlueBook/P6.html#650603030402) |
| `O=C(O)CC(=O)Cl` | `(chlorocarbonyl)acetic acid` · Preferred, exact | `3-chloro-3-oxopropanoic acid` · Preferred, exact. Method (3) cites the chain-end halo and oxo ligands on carbon 3. [P-65.5.4](https://iupac.qmul.ac.uk/BlueBook/P6.html#650504) |
| `CCSC(=O)CCC(=O)SC` | `1-S-ethyl 4-S-methyl butanebis(thioate)` · Preferred, exact | `S-ethyl S-methyl butanebis(thioate)` · Preferred, exact. The acid ends can be exchanged, so position numbers are redundant. The element locants specify sulfur attachment and remain. The chloro-substituted control retains numbers because exchange changes connectivity. [P-65.6.3.3.2.1](https://iupac.qmul.ac.uk/BlueBook/P6.html#650603030201), [P-65.6.3.3.7.1](https://iupac.qmul.ac.uk/BlueBook/P6.html#650603030701) |
| `c1(C(=O)SC)c(C(=O)SC)cc(C(=O)SCC)cc1` | `4-S-ethyl 1,2-bis(S-methyl) benzene-1,2,4-tricarbothioate` · Preferred, strict parse failure | Same Preferred name, now parses exactly. Positions 1,2 specify the two sites, and each repeated methyl is attached through S. [P-65.6.3.3.2.1](https://iupac.qmul.ac.uk/BlueBook/P6.html#650603030201), [P-65.6.3.3.7.1](https://iupac.qmul.ac.uk/BlueBook/P6.html#650603030701) |

The base already passed both original peroxolate cases, both butanediperoxoic parse cases, and the chain-end thiocarboxylic acid case. The head also removes unnecessary organyl position numbers for the symmetric thioperoxy and benzene-1,3,5-tricarboxylate examples in the raw TSV. The tested head declines the sulfur-localized heteroatom-bound thioate `NC(=O)[S-]` because its acid-derived name reads back with charge on oxygen. The thioate equivalence applies to audited carbon- or hydrogen-bound thiocarboxylates. [Issue history](https://github.com/Ameyanagi/nazuchem/issues/138#issuecomment-6871139012), [Thioate regression](https://github.com/Ameyanagi/nazuchem/blob/311c31b6b330e5cc0589cd96f56ee558ac7f42cc/crates/nazuchem/tests/chalcogen_anions.rs#L74), [Chain-end and ester-symmetry regressions](https://github.com/Ameyanagi/nazuchem/blob/311c31b6b330e5cc0589cd96f56ee558ac7f42cc/crates/nazuchem/tests/thio_acids.rs#L356)

### Issue 182

**Anhydride preference follows the cited acid name**

P-65.1.8.1 rules out preferred formic-acid substitution by cyano and azido. P-65.2.1.4 names their preferred acids by replacing a carbonic hydroxy group with cyanido or azido. P-65.7.1 replaces the class word `acid` for two identical halves. P-65.7.2 cites the acid portions alphabetically for two unlike halves, followed by `anhydride`. Applying those rules derives the carbonocyanidic anhydride spellings below. Each generated name strictly round-trips to the exact graph on both revisions. [P-65.1.8.1](https://iupac.qmul.ac.uk/BlueBook/P6.html#65010801), [P-65.2.1.4](https://iupac.qmul.ac.uk/BlueBook/P6.html#65020104), [P-65.7.1](https://iupac.qmul.ac.uk/BlueBook/P6.html#650701), [P-65.7.2](https://iupac.qmul.ac.uk/BlueBook/P6.html#650702)

| Input SMILES | Base `c3b023b` | Head `2f70359` |
| --- | --- | --- |
| `N#CC(=O)O` | `cyanoformic acid` · NotPreferred(UnretainedParent) | `carbonocyanidic acid` · Preferred |
| `N#CC(=O)OC(=O)C#N` | `cyanoformic anhydride` · Preferred | `carbonocyanidic anhydride` · Preferred |
| `N#CC(=O)OC(C)=O` | `acetic cyanoformic anhydride` · Preferred | `acetic carbonocyanidic anhydride` · Preferred with identical result for reverse traversal `CC(=O)OC(=O)C#N` |
| `N#CC(=O)OC(=O)CCC(=O)OC(=O)C#N` | `bis(cyanoformic) butanedioic dianhydride` · Preferred | `bis(carbonocyanidic) butanedioic dianhydride` · Unverified(UnauditedPath), pending audit of the multiplier for this functionally replaced whole-acid citation under P-65.7.6.1 |
| `[N-]=[N+]=NC(=O)OC(=O)N=[N+]=[N-]` | `azidoformic anhydride` · Preferred | Same spelling · NotPreferred(UnretainedParent), matching the acid half. The preferred acid `carbonazidic acid` is printed under P-65.2.1.4. |
| `[N-]=[N+]=NC(=O)OC(C)=O` | `acetic azidoformic anhydride` · Preferred | Same spelling · NotPreferred(UnretainedParent). The same demotion applies to the azidoformic status inherited by a dianhydride. |
| `[O-][N+](=O)C(=O)OC(C)=O` | `acetic nitroformic anhydride` · Preferred | Same spelling and status. The underlying `nitroformic acid` moves from NotPreferred to Preferred in accordance with P-65.1.8.2. |

The strict parser already understood both carbonocyanidic anhydride spellings on the base, so their promotion is a name-generation fix. The azidoformic cases directly verify conservative status inheritance independently of spelling changes. The preferable carbonazidic spelling remains outside the observed implementation. [P-65.1.8.2](https://iupac.qmul.ac.uk/BlueBook/P6.html#65010802), [P-65.7.6.1](https://iupac.qmul.ac.uk/BlueBook/P6.html#65070601), [Nitrile regression](https://github.com/Ameyanagi/nazuchem/blob/311c31b6b330e5cc0589cd96f56ee558ac7f42cc/crates/nazuchem/tests/one_carbon_nitrile.rs#L112), [Preference-inheritance regression](https://github.com/Ameyanagi/nazuchem/blob/311c31b6b330e5cc0589cd96f56ee558ac7f42cc/crates/nazuchem/tests/polyanhydrides.rs#L485)

## Hydrazines and nitramides

### Issue 113

**Hydrazido prefixes, semicarbazones, hydrazine suffixes**

**Audit correction for sulfur attachment.** The first PR head `311c31b` gave Preferred `N-[5-(2-bromobenzohydrazidosulfonyl)-2-methoxyphenyl]acetamide` for the reported structure. That name has an exact read-back. The free valence is on sulfur, which selects the acid-derived hydrazine-1-sulfonyl route. The hydrazido route applies to an N-attached prefix. The corrected name is rule-derived. Substituted acyl prefixes take their own enclosing marks. Both N-methyl cases use graph-exact general names with `NotPreferred(SulfonylPrefix)`. [Issue #113](https://github.com/Ameyanagi/nazuchem/issues/113#issuecomment-6851691239) · [P-66.3.2.1](https://iupac.qmul.ac.uk/BlueBook/P6a.html#66030201) · [P-65.3.2.2.2](https://iupac.qmul.ac.uk/BlueBook/P6.html#6503020202) · [P-66.3.2.3](https://iupac.qmul.ac.uk/BlueBook/P6a.html#66030203) · [P-16.5.1.1](https://iupac.qmul.ac.uk/BlueBook/P1.html#16050101)

| Case and input | Base result | Tested result |
| --- | --- | --- |
| `113-reported-acyl-sulfonyl`<br>`CC(=O)Nc1cc(S(=O)(=O)NNC(=O)c2ccccc2Br)ccc1OC` | `NotPreferred(AcylPrefix)` · `N-(5-{2-[(2-bromophenyl)carbonyl]hydrazin-1-ylsulfonyl}-2-methoxyphenyl)acetamide` · `exact` | `Preferred` · `N-{5-[2-(2-bromobenzoyl)hydrazine-1-sulfonyl]-2-methoxyphenyl}acetamide` · `exact` |
| `113-sulfonyl-audit-reverse`<br>`COc1c(NC(C)=O)cc(S(=O)(=O)NNC(=O)c2ccccc2Br)cc1` | `NotPreferred(AcylPrefix)` · `N-(5-{2-[(2-bromophenyl)carbonyl]hydrazin-1-ylsulfonyl}-2-methoxyphenyl)acetamide` · `exact` | `Preferred` · `N-{5-[2-(2-bromobenzoyl)hydrazine-1-sulfonyl]-2-methoxyphenyl}acetamide` · `exact` |
| `113-minimal-acyl-sulfonyl`<br>`O=C(O)CS(=O)(=O)NNC(C)=O` | `NotPreferred(AminoPrefix)` · `(2-acetylhydrazin-1-ylsulfonyl)acetic acid` · `exact` | `Preferred` · `(2-acetylhydrazine-1-sulfonyl)acetic acid` · `exact` |
| `113-minimal-acyl-sulfonyl-reverse`<br>`CC(=O)NNS(=O)(=O)CC(O)=O` | `NotPreferred(AminoPrefix)` · `(2-acetylhydrazin-1-ylsulfonyl)acetic acid` · `exact` | `Preferred` · `(2-acetylhydrazine-1-sulfonyl)acetic acid` · `exact` |
| `113-sulfonyl-chloro`<br>`OC(=O)CS(=O)(=O)NNC(=O)CCl` | `NotPreferred(AminoPrefix)` · `[2-(2-chloroacetyl)hydrazin-1-ylsulfonyl]acetic acid` · `exact` | `Preferred` · `[2-(chloroacetyl)hydrazine-1-sulfonyl]acetic acid` · `exact` |
| `113-sulfonyl-hydroxy`<br>`OC(=O)CS(=O)(=O)NNC(=O)CO` | `NotPreferred(AminoPrefix)` · `[2-(2-hydroxyacetyl)hydrazin-1-ylsulfonyl]acetic acid` · `exact` | `Preferred` · `[2-(hydroxyacetyl)hydrazine-1-sulfonyl]acetic acid` · `exact` |
| `113-sulfonyl-n1-alkyl`<br>`OC(=O)CS(=O)(=O)N(C)NC(C)=O` | `NotPreferred(AminoPrefix)` · `(2-acetyl-1-methylhydrazin-1-ylsulfonyl)acetic acid` · `exact` | `NotPreferred(SulfonylPrefix)` · `(N'-methylacetohydrazidosulfonyl)acetic acid` · `exact` |
| `113-sulfonyl-n2-alkyl`<br>`OC(=O)CS(=O)(=O)NN(C)C(C)=O` | `NotPreferred(AminoPrefix)` · `(2-acetyl-2-methylhydrazin-1-ylsulfonyl)acetic acid` · `exact` | `NotPreferred(SulfonylPrefix)` · `(N-methylacetohydrazidosulfonyl)acetic acid` · `exact` |


- P-66.3.2.3 method (1) generates the preferred prefix for R–C(=O)–NH–NH– from the corresponding hydrazide. The acetohydrazido examples, including the N-ethyl sulfonic acid, are printed PINs. P-68.3.1.2.6 likewise directly prints `3-(hydrazinecarbohydrazido)propanoic acid`. The bromo and cyclohexane acid cases and the propanamide use the same formation rules. [P-66.3.2.3](https://iupac.qmul.ac.uk/BlueBook/P6a.html#66030203) · [P-66.3.1](https://iupac.qmul.ac.uk/BlueBook/P6a.html#660301) · [P-68.3.1.2.6](https://iupac.qmul.ac.uk/BlueBook/P6a.html#6803010206) · [P-16.5.1.4](https://iupac.qmul.ac.uk/BlueBook/P1.html#16050104)
- P-68.3.1.2.5 uses hydrazinecarboxamide as the functional parent for semicarbazone PINs and demonstrates substitution at 2. P-29.6.2.1 retains unsubstituted benzylidene as a preferred prefix. P-93.4.2.1.3 specifies E/Z for double bonds involving heteroatoms. Combining those rules assigns E or Z to the C=N at hydrazine locant 2. The Blue Book does not give the combined `(2E)` and `(2Z)` semicarbazone examples; these are rule-derived and both configurations read back exactly in separate tests. [P-68.3.1.2.5](https://iupac.qmul.ac.uk/BlueBook/P6a.html#6803010205) · [P-29.6.2.1](https://iupac.qmul.ac.uk/BlueBook/P2.html#29060201) · [P-93.4.2.1.3](https://iupac.qmul.ac.uk/BlueBook/P9.html#9304020103)
- The hydrazine parent takes numeric locants 1 and 2. A polyamide distinguishes suffix nitrogen positions by superscripting the site of attachment to the parent, giving N^1 and N^2 in the dicarboxamide. Both acid groups attach to nitrogen 1 in the geminal diacid, giving 1,1. For the peroxide, the acid suffix for –C(=O)OOH is `carboperoxoic acid`, and the methyl ester applies the acid-derived `carboperoxoate` ending. These are applications of the general rules. [P-68.3.1.2.1](https://iupac.qmul.ac.uk/BlueBook/P6a.html#6803010201) · [P-66.1.1.3.1.1](https://iupac.qmul.ac.uk/BlueBook/P6a.html#660101030101) · [P-65.1.1](https://iupac.qmul.ac.uk/BlueBook/P6.html#650101) · [P-65.1.4.1](https://iupac.qmul.ac.uk/BlueBook/P6.html#65010401) · [P-65.6.1](https://iupac.qmul.ac.uk/BlueBook/P6.html#650601)

| Input SMILES | Base name and status | HEAD name and status | Strict identity, base → HEAD |
| --- | --- | --- | --- |
| `NNC(=O)OOC` | `1-(methylperoxy)formohydrazide` · NotPreferred(UnretainedParent) | `methyl hydrazinecarboperoxoate` · Preferred | exact → exact |
| `COOC(=O)NN` | `1-(methylperoxy)formohydrazide` · NotPreferred(UnretainedParent) | `methyl hydrazinecarboperoxoate` · Preferred | exact → exact |
| `CNC(=O)NNC(=O)NC` | `N-methyl-1-[2-(methylcarbamoyl)hydrazin-1-yl]formamide` · NotPreferred(UnretainedParent) | `N^1,N^2-dimethylhydrazine-1,2-dicarboxamide` · Preferred | exact → exact |
| `NN(C(=O)O)C(=O)O` | `1-hydroxy-N-(hydroxycarbonyl)formohydrazide` · NotPreferred(UnretainedParent) | `hydrazine-1,1-dicarboxylic acid` · Preferred | exact → exact |
| `O=C(O)c1ccc(NNC(C)=O)cc1` | `4-(2-acetylhydrazin-1-yl)benzoic acid` · NotPreferred(AminoPrefix) | `4-acetohydrazidobenzoic acid` · Preferred | exact → exact |
| `O=C(O)c1ccc(NNC(=O)c2ccccc2Br)cc1` | `4-{2-[(2-bromophenyl)carbonyl]hydrazin-1-yl}benzoic acid` · NotPreferred(AcylPrefix) | `4-(2-bromobenzohydrazido)benzoic acid` · Preferred | exact → exact |
| `O=C(O)c1ccc(NNC(=O)C2CCCCC2)cc1` | `4-[2-(cyclohexylcarbonyl)hydrazin-1-yl]benzoic acid` · NotPreferred(AcylPrefix) | `4-(cyclohexanecarbohydrazido)benzoic acid` · Preferred | exact → exact |
| `O=S(=O)(O)c1ccc(NN(CC)C(C)=O)cc1` | `4-(2-acetyl-2-ethylhydrazin-1-yl)benzene-1-sulfonic acid` · NotPreferred(AminoPrefix) | `4-(N-ethylacetohydrazido)benzene-1-sulfonic acid` · Preferred | exact → exact |
| `O=C(O)CCNNC(=O)NN` | `3-[2-(hydrazinylcarbonyl)hydrazin-1-yl]propanoic acid` · NotPreferred(UnretainedParent) | `3-(hydrazinecarbohydrazido)propanoic acid` · Preferred | exact → exact |
| `NC(=O)CCNNC(=O)NN` | `3-[2-(hydrazinylcarbonyl)hydrazin-1-yl]propanamide` · NotPreferred(UnretainedParent) | `3-(hydrazinecarbohydrazido)propanamide` · Preferred | exact → exact |
| `c1ccccc1/C=N/NC(N)=O` | `1-({[(E)-phenylmethylidene]amino}amino)formamide` · NotPreferred(AminoPrefix) | `(2E)-2-benzylidenehydrazine-1-carboxamide` · Preferred | exact → exact |
| `c1ccccc1/C=N\NC(N)=O` | `1-({[(Z)-phenylmethylidene]amino}amino)formamide` · NotPreferred(AminoPrefix) | `(2Z)-2-benzylidenehydrazine-1-carboxamide` · Preferred | exact → exact |

The issue's original mixed anhydride, cyanohydrazide, sodium salt and numbered hydrazinyl prefixes were already Preferred and exact at the base. Configured acylhydrazone and N-linked configured hydrazinyl examples still use exact, NotPreferred fallbacks. Mixed Se/S and Te/S, sulfur-oxoacid neighbours and hydrazinecarboxylic dianhydrides still need coverage. The ring-nitrogen and second-hydrazine probes are controls for the tested structures and do not resolve equal-maximum parent ties. [#113](https://github.com/Ameyanagi/nazuchem/issues/113#issuecomment-6851691239) · [Earlier PR #167](https://github.com/Ameyanagi/nazuchem/pull/167) · [New regressions](https://github.com/Ameyanagi/nazuchem/blob/311c31b6b330e5cc0589cd96f56ee558ac7f42cc/crates/nazuchem/tests/hydrazine_names.rs#L288-L333)

### Issue 141

**O- and N-linked hydrazinecarbonyl prefixes**

The printed functional-class PINs `dicarbonic dihydrazide` and `2-imidodicarbonic dihydrazide` strict-parse to the designated structures on both revisions. The standalone generated forms still have a citation gap. In particular, the engine reports Preferred `hydrazinecarboxylic anhydride` for the dicarbonic structure. The six linked-prefix regressions address structures with a higher parent group. [P-66.3.5.2](https://iupac.qmul.ac.uk/BlueBook/P6a.html#66030502)

| Case and input | Base result | Tested result |
| --- | --- | --- |
| `141-printed-oxy`<br>`NNC(=O)OC(=O)NN` | `Preferred` · `hydrazinecarboxylic anhydride` · `exact` | `Preferred` · `hydrazinecarboxylic anhydride` · `exact` |
| `141-printed-amino`<br>`NNC(=O)NC(=O)NN` | `NotPreferred(UnretainedParent)` · `1-hydrazinyl-N-(hydrazinylcarbonyl)formamide` · `exact` | `NotPreferred(UnretainedParent)` · `N-(hydrazinecarbonyl)-1-hydrazinylformamide` · `exact` |
| `141-dicarbonic-pin-parse`<br>`dicarbonic dihydrazide` | `parse` · `C(=O)(OC(=O)NN)NN` | `parse` · `C(=O)(OC(=O)NN)NN` |
| `141-imidodicarbonic-pin-parse`<br>`2-imidodicarbonic dihydrazide` | `parse` · `C(=O)(NN)NC(=O)NN` | `parse` · `C(=O)(NN)NC(=O)NN` |


P-66.3.5.2 prints `[(hydrazinecarbonyl)oxy]formohydrazide` and `[(hydrazinecarbonyl)amino]formohydrazide` as structural-name examples. The functional-class dicarbonic and 2-imidodicarbonic dihydrazide names receive the PIN labels there. P-66.3.2.1 selects the hydrazine acid’s acyl prefix `hydrazinecarbonyl`. P-16.5.1.4 requires the inner parentheses around the parent-hydride-containing acyl prefix, P-16.5.2.4 brackets the oxy/amino prefix, and P-16.5.3.1 braces the outer ethyl substituent when it already contains both. The six new acid names apply those prefix and enclosure rules. [P-66.3.5.2](https://iupac.qmul.ac.uk/BlueBook/P6a.html#66030502) · [P-66.3.2.1](https://iupac.qmul.ac.uk/BlueBook/P6a.html#66030201) · [P-16.5.1.4](https://iupac.qmul.ac.uk/BlueBook/P1.html#16050104) · [P-16.5.2.4](https://iupac.qmul.ac.uk/BlueBook/P1.html#16050204) · [P-16.5.3.1](https://iupac.qmul.ac.uk/BlueBook/P1.html#16050301)

| Input SMILES | Base name and status | HEAD name and status | Strict identity, base → HEAD |
| --- | --- | --- | --- |
| `OC(=O)CCOC(=O)NN` | `3-[(hydrazinylcarbonyl)oxy]propanoic acid` · NotPreferred(UnretainedParent) | `3-[(hydrazinecarbonyl)oxy]propanoic acid` · Preferred | exact → exact |
| `OC(=O)CCNC(=O)NN` | `3-[(hydrazinylcarbonyl)amino]propanoic acid` · NotPreferred(UnretainedParent) | `3-[(hydrazinecarbonyl)amino]propanoic acid` · Preferred | exact → exact |
| `OC(=O)c1ccc(OC(=O)NN)cc1` | `4-[(hydrazinylcarbonyl)oxy]benzoic acid` · NotPreferred(UnretainedParent) | `4-[(hydrazinecarbonyl)oxy]benzoic acid` · Preferred | exact → exact |
| `OC(=O)c1ccc(NC(=O)NN)cc1` | `4-[(hydrazinylcarbonyl)amino]benzoic acid` · NotPreferred(UnretainedParent) | `4-[(hydrazinecarbonyl)amino]benzoic acid` · Preferred | exact → exact |
| `OC(=O)c1ccc(CCOC(=O)NN)cc1` | `4-{2-[(hydrazinylcarbonyl)oxy]ethyl}benzoic acid` · NotPreferred(UnretainedParent) | `4-{2-[(hydrazinecarbonyl)oxy]ethyl}benzoic acid` · Preferred | exact → exact |
| `OC(=O)c1ccc(CCNC(=O)NN)cc1` | `4-{2-[(hydrazinylcarbonyl)amino]ethyl}benzoic acid` · NotPreferred(UnretainedParent) | `4-{2-[(hydrazinecarbonyl)amino]ethyl}benzoic acid` · Preferred | exact → exact |

Both printed structural names already strict-parse to the matching graph on base and HEAD. So do the earlier repairs `4-(2-hydrazinylethyl)benzoic acid`, `4-(hydrazinecarbonyl)benzene-1-carboximidamide`, and the chain-end `4-hydrazinyl-4-oxobutanoic acid`. The base’s actual six linked spellings use `hydrazinylcarbonyl`. The older issue comment records `aminocarbamoyloxy` and `aminocarbamoylamino`. Methyl on the bridge nitrogen or the terminal hydrazine nitrogen still yields NotPreferred fallback names in the tested examples. [#141 scope](https://github.com/Ameyanagi/nazuchem/issues/141#issuecomment-6842278533) · [Earlier PR #146](https://github.com/Ameyanagi/nazuchem/pull/146) · [Earlier PR #167](https://github.com/Ameyanagi/nazuchem/pull/167) · [New regressions](https://github.com/Ameyanagi/nazuchem/blob/311c31b6b330e5cc0589cd96f56ee558ac7f42cc/crates/nazuchem/tests/hydrazine_names.rs#L31-L86)

### Issue 176

**Hydrazine ions, hydrazide locants, nitric and nitrous derivatives**

- P-73.1.1.2 gives `pentamethylhydrazinium` as the PIN and the general `-ium` hydride suffix. P-73.1.2 and P-73.1.2.1 select a large neutral functional parent and give `N',N',N'-trimethylbenzohydrazid-N'-ium` as the printed hydrazide ion example. The acyl/ylidene cation applies both hydrazide locant and cation rules by analogy. P-72.2.2.1 gives a locanted `-ide` for deprotonation of a parent hydride. The anion uses it at hydrazine N-1. [P-73.1.1.2](https://iupac.qmul.ac.uk/BlueBook/P7.html#73010102) · [P-73.1.2](https://iupac.qmul.ac.uk/BlueBook/P7.html#730102) · [P-73.1.2.1](https://iupac.qmul.ac.uk/BlueBook/P7.html#73010201) · [P-66.3.3.1](https://iupac.qmul.ac.uk/BlueBook/P6a.html#66030301) · [P-72.2.2.1](https://iupac.qmul.ac.uk/BlueBook/P7.html#72020201)
- P-66.3.3.2 assigns superscripted attachment positions to N and N′ for two hydrazides. When both positions are the same, the geminal result combines the distinct priming series with the shared superscript. It follows the numbering analogue in P-66.1.1.3.1.2 for geminal carboxamides. P-66.4.2.1 provides additional support for primed and superscripted locants on a nitrogen-rich suffix. The Blue Book has no directly worked geminal hydrazide example. P-14.3.5 and the locant comparison in P-66.3.3.2 support choosing N^1 before N′^1. For the ylidene multiplier, P-68.3.1.2.3.1 explicitly prints `di(propan-2-ylidene)hydrazine`, an analogue to the corrected dihydrazide multiplier. [P-66.3.3.2](https://iupac.qmul.ac.uk/BlueBook/P6a.html#66030302) · [P-66.1.1.3.1.2](https://iupac.qmul.ac.uk/BlueBook/P6a.html#660101030102) · [P-66.4.2.1](https://iupac.qmul.ac.uk/BlueBook/P6a.html#66040201) · [P-14.3.5](https://iupac.qmul.ac.uk/BlueBook/P1.html#140305) · [P-68.3.1.2.3.1](https://iupac.qmul.ac.uk/BlueBook/P6a.html#680301020301)
- P-67.1.2.6.3 specifies nitric and nitrous hydrazides as parent names and prints `N'-hexylidenenitrous hydrazide` as a PIN. The corresponding methyl substitutions use N′, the hydrazide nitrogen away from the nitric or nitrous group. P-93.4.2.1.3 places E/Z in parentheses at the front if the double bond has no numeric locant. [P-67.1.2.6.3](https://iupac.qmul.ac.uk/BlueBook/P6a.html#6701020603) · [P-66.3.3.1](https://iupac.qmul.ac.uk/BlueBook/P6a.html#66030301) · [P-93.4.2.1.3](https://iupac.qmul.ac.uk/BlueBook/P9.html#9304020103)

| Input SMILES | Base name and status | HEAD name and status | Strict identity, base → HEAD |
| --- | --- | --- | --- |
| `C[N+](C)(C)NC(=O)c1ccccc1` | `N-benzamido-N,N-dimethylmethanaminium` · Unverified(CationParent) | `N',N',N'-trimethylbenzohydrazid-N'-ium` · Preferred | exact → exact |
| `C[N+](C)(C)N(C)C` | `N-(dimethylamino)-N,N-dimethylmethanaminium` · Unverified(CationParent) | `pentamethylhydrazinium` · Preferred | exact → exact |
| `CC(=O)N[N+](C)=C(C)C` | `N-acetamido-N-methylpropan-2-iminium` · Unverified(CationParent) | `N'-methyl-N'-(propan-2-ylidene)acetohydrazid-N'-ium` · Preferred | exact → exact |
| `N[NH3+]` | name-error · `openclatura: no supported parent skeleton` | `hydrazinium` · Preferred | — → exact |
| `N[N+](C)=C(C)C` | `N-amino-N-methylpropan-2-iminium` · Unverified(CationParent) | `1-methyl-1-(propan-2-ylidene)hydrazinium` · Preferred | exact → exact |
| `CC(C)=N[N-]C` | `N-[(propan-2-ylidene)amino]methanaminide` · Preferred | `1-methyl-2-(propan-2-ylidene)hydrazin-1-ide` · Preferred | exact → exact |
| `CC(C)=NNC(=O)CCC(=O)NN=C(C)C` | `N'^1,N'^4-bis(propan-2-ylidene)butanedihydrazide` · Preferred | `N'^1,N'^4-di(propan-2-ylidene)butanedihydrazide` · Preferred | exact → exact |
| `C1CCCCC1(C(=O)NNC)C(=O)NNC` | `N',N'''-dimethylcyclohexane-1,1-dicarbohydrazide` · Preferred | `N'^1,N'''^1-dimethylcyclohexane-1,1-dicarbohydrazide` · Preferred | exact → exact |
| `C1CCCCC1(C(=O)N(C)N)C(=O)NNCC` | `N'''-ethyl-N-methylcyclohexane-1,1-dicarbohydrazide` · Preferred | `N'''^1-ethyl-N^1-methylcyclohexane-1,1-dicarbohydrazide` · Preferred | exact → exact |
| `CCNNC(=O)C1(C(=O)N(C)N)CCCCC1` | `N'-ethyl-N''-methylcyclohexane-1,1-dicarbohydrazide` · Preferred | `N'''^1-ethyl-N^1-methylcyclohexane-1,1-dicarbohydrazide` · Preferred | exact → exact |
| `C1CCCCC1(C(=O)N(C)N)C(=O)NN(C)C` | `N,N''',N'''-trimethylcyclohexane-1,1-dicarbohydrazide` · Preferred | `N^1,N'''^1,N'''^1-trimethylcyclohexane-1,1-dicarbohydrazide` · Preferred | exact → exact |
| `CCCCCC=NNN=O` | `N-nitrosohexanal hydrazone` · NotPreferred(FunctionalClass) | `N'-hexylidenenitrous hydrazide` · Preferred | exact → exact |
| `CCCCC/C=N/NN=O` | `(E)-N-nitrosohexanal hydrazone` · NotPreferred(FunctionalClass) | `(E)-N'-hexylidenenitrous hydrazide` · Preferred | exact → exact |
| `c1ccccc1/C=N/[N+](=O)[O-]` | `(1E)-N-nitro-1-phenylmethanimine` · Preferred | `(E)-benzylidenenitramide` · Preferred | exact → exact |
| `c1ccccc1/C=N\[N+](=O)[O-]` | `(1Z)-N-nitro-1-phenylmethanimine` · Preferred | `(Z)-benzylidenenitramide` · Preferred | exact → exact |
| `CNNN=O` | `(nitrosohydrazinyl)methane` · Preferred | `N'-methylnitrous hydrazide` · Preferred | exact → exact |
| `CNN[N+](=O)[O-]` | `(nitrohydrazinyl)methane` · Preferred | `N'-methylnitric hydrazide` · Preferred | exact → exact |
| `C[NH2+]N` | `N-aminomethanaminium` · Unverified(CationParent) | `1-methylhydrazinium` · Preferred | exact → exact |
| `CN[NH3+]` | `N-azaniumylmethanamine` · Unverified(CationParent) | `2-methylhydrazinium` · Preferred | exact → exact |
| `C[NH2+]NC` | `N-(methylamino)methanaminium` · Unverified(CationParent) | `1,2-dimethylhydrazinium` · Preferred | exact → exact |

The original dimethylbutanedihydrazide, ethylnitramide, and dipropylnitrous amide were already Preferred and exact. The reported hydrazone cation already kept its charge and constitution in the base’s exact iminium fallback, so its change here is parent and preference. Ring and multisite nitramides and geminal diamides are unchanged. `NNC(=O)CCN[N+](=O)[O-]` still gets `3-nitramidopropanehydrazide`, Preferred and exact. P-41 places amides at 11 ahead of hydrazides at 12, and P-67.1.2.6.3 identifies nitramide as an amide. The rule-derived `(3-hydrazinyl-3-oxopropyl)nitramide` strict-parses exactly on both builds. [#176 narrowed scope](https://github.com/Ameyanagi/nazuchem/issues/176#issuecomment-6871139345) · [Earlier PR #185](https://github.com/Ameyanagi/nazuchem/pull/185) · [P-41](https://iupac.qmul.ac.uk/BlueBook/P4.html#41) · [P-67.1.2.6.3](https://iupac.qmul.ac.uk/BlueBook/P6a.html#6701020603) · [P-66.3.2.2](https://iupac.qmul.ac.uk/BlueBook/P6a.html#66030202)

[Geminal hydrazide regressions](https://github.com/Ameyanagi/nazuchem/blob/311c31b6b330e5cc0589cd96f56ee558ac7f42cc/crates/nazuchem/tests/dihydrazides_nitramides.rs#L115-L172) · [Nitramide regressions](https://github.com/Ameyanagi/nazuchem/blob/311c31b6b330e5cc0589cd96f56ee558ac7f42cc/crates/nazuchem/tests/dihydrazides_nitramides.rs#L216-L238) · [Hydrazine ion regressions](https://github.com/Ameyanagi/nazuchem/blob/311c31b6b330e5cc0589cd96f56ee558ac7f42cc/crates/nazuchem/tests/dihydrazides_nitramides.rs#L243-L281)

## Ions and anionic parents

### Issue 180

**parent selection and peroxy-acid anion citation**

[#180's scope after batch 5](https://github.com/Ameyanagi/nazuchem/issues/180#issuecomment-6871139646) lists two remaining parent-selection cases and preferred spelling for a peroxy-acid anion substituent. [PR #187](https://github.com/Ameyanagi/nazuchem/pull/187) had already fixed amide deprotonation and the original polyolate ranking. The base and head both return preferred, exact `acetylazanide`, `1-carboxylatoethane-1,2-bis(olate)`, and `(azanidylmethyl)propanedioate` for the relevant issue examples.

| SMILES | Base name, status, identity | Head name, status, identity | Head charge |
| --- | --- | --- | --- |
| `[O-]c1ccc(C2CC([O-])CC2)cc1CC(=O)[O-]` | `2-(carboxylatomethyl)-4-(3-oxidocyclopentyl)phenoxide`, preferred, exact | `[2-oxido-5-(3-oxidocyclopentyl)phenyl]acetate`, preferred, exact | −3 |
| `[O-]C1CC([O-])CC1[CH-][CH-]C(=O)[O-]` | Declined, anionic parent choice unresolved | `1-(2,4-dioxidocyclopentyl)propane-1,2-diid-3-oate`, preferred, exact | −5 |
| `[O-]C1CC([O-])CC1[CH-]CC(=O)[O-]` | Declined, anionic parent choice unresolved | `1-(2,4-dioxidocyclopentyl)propan-1-id-3-oate`, preferred, exact | −4 |
| `[O-]C1CC([O-])CC1CCC(=O)[O-]` | `4-(2-carboxylatoethyl)cyclopentane-1,3-bis(olate)`, preferred, exact | Same name, preferred, exact | −3 |
| `[O-]C(=O)c1ccc(-c2ccc(C(=O)[O-])cc2)cc1CC([O-])C([O-])C` | `3-(2,3-dioxidobutyl)[1,1'-biphenyl]-4,4'-dicarboxylate`, preferred, exact | Same name, preferred, exact | −4 |
| `C1CC([NH-])CC1C2CC([NH-])CC2CC(C(=O)[O-])C[CH2-]` | Declined, anionic parent choice unresolved | `5-(2-carboxylatobutan-4-id-1-yl)[1,1'-bi(cyclopentane)]-3,3'-bis(aminide)`, preferred, exact | −4 |
| `[NH-]CC(=O)O[O-]` | `1-(dioxidanidylcarbonyl)methanaminide`, unverified (`AnionParent`), exact | `1-carboperoxoatomethanaminide`, preferred, exact | −2 |
| `[O-]OC(=O)CC(=O)[O-]` | `(dioxidanidylcarbonyl)acetate`, unverified (`AnionParent`), exact | `carboperoxoatoacetate`, preferred, exact | −2 |
| `[O-]OC(=O)c1ccccc1C(=O)[O-]` | `2-(dioxidanidylcarbonyl)benzoate`, unverified (`AnionParent`), exact | `2-carboperoxoatobenzoate`, preferred, exact | −2 |
| `[O-]OC(=O)CC(C(=O)O[O-])(C(=O)O[O-])CC(=O)O[O-]` | `3,3-bis(dioxidanidylcarbonyl)pentanediperoxoate`, unverified (`AnionParent`), exact | `3,3-dicarboperoxoatopentanediperoxoate`, preferred, exact | −4 |

- **Different cyclic parents.** The six-membered benzene ring and the saturated five-membered ring in the first row are inequivalent under P-28.1 and P-28.7. Each ring contributes one olate to its own candidacy and the acetate chain has one carboxylate. P-72.7(a) ties all three and the candidates have no 'ide'/'uide' centres. At P-72.7(d), each candidate's anionic atom is oxygen. P-72.7(e) favours carboxylate over olate, as it does in the printed example `3-oxidonaphthalene-2-carboxylate`. The reversed SMILES gives the same before and after names. [P-28.1](https://iupac.qmul.ac.uk/BlueBook/P2.html#2801), [P-28.7](https://iupac.qmul.ac.uk/BlueBook/P2.html#2807), [P-72.7](https://iupac.qmul.ac.uk/BlueBook/P7.html#7207)
- **Skeletal charges.** With two carbanions the chain has three centres and the ring has two, so criterion (a) decides. With one carbanion both candidates have two, and the skeletal 'ide' decides at (b). With zero carbanions the ring has two to the chain's one and wins at (a). P-72.5.3 gives the skeletal anions first claim to low locants and allows the added 'oate' suffix. All five formal charges are present in the first parse-back, including both skeletal charges. [P-72.5.3](https://iupac.qmul.ac.uk/BlueBook/P7.html#720503), [P-72.7](https://iupac.qmul.ac.uk/BlueBook/P7.html#7207)
- **Identical directly joined rings.** Biphenyl and bi(cyclopentane) qualify as parents under P-28.2.1 before the anionic score. Biphenyl's two carboxylates tie the side chain's two olates and win under P-72.7(e). The two cyclopentane aminides tie the carboxylate chain's carbanion and carboxylate at (a), then win by two 'ide' centres to one at (b). The terminal carbanion remains in the `butan-4-id-1-yl` prefix under P-72.6.3. Both SMILES traversals preserve four negative atoms. [P-28.2.1](https://iupac.qmul.ac.uk/BlueBook/P2.html#280201), [P-72.2.2.2.3](https://iupac.qmul.ac.uk/BlueBook/P7.html#7202020203), [P-72.6.3](https://iupac.qmul.ac.uk/BlueBook/P7.html#720603), [P-72.7](https://iupac.qmul.ac.uk/BlueBook/P7.html#7207)
- **Peroxy-acid spelling.** P-65.1.4.1 assigns `carboperoxoic acid` to the external –C(=O)OOH group. P-72.2.2.2.1.1 gives `carboperoxoate` on deprotonation. P-72.6.1 gives `carboperoxoato` for the acid-derived anionic group attached by a single bond. Aminide wins the first row's parent tie under P-72.7(b). Carboxylate wins the acetate and benzoate ties under P-72.7(e) and Table 4.4's suffix order, consistent with the carboxylic/peroxycarboxylic order in Table 4.3. The central double branch in the geminal case produces two pendant prefixes and two chain peroxoate suffixes. Carbon-13- and oxygen-18-labelled peroxy prefixes remain unverified for isotope placement, with exact isotope and graph identity, and the sulfur analogue in the raw rows is declined on both revisions. [P-65.1.4.1](https://iupac.qmul.ac.uk/BlueBook/P6.html#65010401), [P-72.2.2.2.1.1](https://iupac.qmul.ac.uk/BlueBook/P7.html#720202020101), [P-72.6.1](https://iupac.qmul.ac.uk/BlueBook/P7.html#720601), [Table 4.3](https://iupac.qmul.ac.uk/BlueBook/P4.html#t403), [Table 4.4](https://iupac.qmul.ac.uk/BlueBook/P4.html#t404)

[Parent regression](https://github.com/Ameyanagi/nazuchem/blob/311c31b6b330e5cc0589cd96f56ee558ac7f42cc/crates/nazuchem/tests/anionic_prefixes.rs#L319-L421) · [Peroxy prefixes](https://github.com/Ameyanagi/nazuchem/blob/311c31b6b330e5cc0589cd96f56ee558ac7f42cc/crates/nazuchem/tests/anionic_prefixes.rs#L494-L566) · [Aminide assembly](https://github.com/Ameyanagi/nazuchem/blob/311c31b6b330e5cc0589cd96f56ee558ac7f42cc/crates/nazuchem/tests/anionic_prefixes.rs#L569-L593)

### Issue 186

**charged characteristic groups and exact identity**

[#186](https://github.com/Ameyanagi/nazuchem/issues/186) and [PR #185](https://github.com/Ameyanagi/nazuchem/pull/185) reported two neutral parse-backs on the released Python package. Against the actual PR base the isothiocyanate still loses charge and an N-methyl bond, while the hydrazide already declines with `anionic nitrogen in a neutral suffix is not supported`. Strict parsing the older spellings independently on both revisions yields neutral `N(=C=S)CC` and neutral `CNNC(C)=O`, respectively.

| SMILES | Base name, status, identity | Head name, status, identity | Input → base → head net charge |
| --- | --- | --- | --- |
| `CC[N+](C)=C=S` | `isothiocyanatoethane`, unverified (`CationParent`), mismatch | `N-ethyl-N-methyl-1-sulfanylidenemethaniminium`, preferred, exact | +1 → 0 → +1 |
| `CC[N+](C)=C=O` | `isocyanatoethane`, unverified (`CationParent`), mismatch | `N-ethyl-N-methyl-1-oxomethaniminium`, preferred, exact | +1 → 0 → +1 |
| `CN[N-]C(C)=O` | Declined, anionic nitrogen in a neutral suffix | `1-acetyl-2-methylhydrazin-1-ide`, preferred, exact | −1 → error → −1 |
| `C1CCC[N+]1=C=S` | `2-isothiocyanatopyrrolidin-1-ium`, preferred, mismatch | `1-(sulfanylidenemethylidene)pyrrolidin-1-ium`, preferred, exact | +1 → +1 → +1 |
| `[O-]C(=O)CC[N+](C)=C=S` | `3-isothiocyanatopropanoate`, preferred, mismatch | `3-[methyl(sulfanylidenemethylidene)ammonio]propanoate`, unverified (`CationPrefix`), exact | 0 → −1 → 0 |
| `C[NH+](C)[N+](C)=C=S` | `N-isothiocyanato-N-methylmethanaminium`, unverified (`CationParent`), mismatch | `N-methyl-N-[methyl(sulfanylidenemethylidene)ammonio]methanaminium`, unverified (`CationPrefix`), exact | +2 → +1 → +2 |
| `CC[N+](C)=C=S.[Cl-]` | `isothiocyanatoethane chloride`, unverified (`CationParent`), mismatch | `N-ethyl-N-methyl-1-sulfanylidenemethaniminium chloride`, unverified (`UnauditedPath`), exact | 0 → −1 → 0 |
| `C[NH+](C)C=O` | `N,N-dimethylformamide`, preferred, mismatch | `N,N-dimethylformamidium`, preferred, exact | +1 → 0 → +1 |
| `C[NH+](C)C=S` | `N,N-dimethylmethanethioamide`, preferred, mismatch | `N,N-dimethylmethanethioamidium`, preferred, exact | +1 → 0 → +1 |
| `[NH2+]([O-])C(=O)C` | `N-oxidoacetamide`, preferred, mismatch | `N-oxidoacetamidium`, unverified (`CationParent`), exact | 0 → −1 → 0 |
| `[NH3+]C(=O)c1ccc(C(=O)[NH3+])cc1` | `benzene-1,4-dicarboxamide`, preferred, mismatch | `benzene-1,4-bis(carboxamidium)`, preferred, exact | +2 → 0 → +2 |
| `CCN=C=[SH+]` | `isothiocyanatoethane`, unverified (`CationParent`), mismatch | Declined, charged sulfur oxo/imino substituent role unsupported | +1 → 0 → error |

- **Cumulene charge and attachment.** The neutral groups in P-61.8 denote –N=C=O and its chalcogen analogues. The charged N in the first two rows is double-bonded to carbon and bears both alkyl groups. P-73.1.2.1 and Table 7.4 give 'iminium' for a cationic imine characteristic group. For the ring cation the positive N is skeletal and accepts the 'ium' suffix under P-73.1.1.2, and P-65.2.1.8 gives the prefix for the exocyclic =C=S group. The base ring spelling preserves net +1 while moving the C=S unit to a neutral extra N and bonding it through a ring carbon. Strict identity catches the difference. [P-61.8](https://iupac.qmul.ac.uk/BlueBook/P6.html#6108), [P-65.2.1.8](https://iupac.qmul.ac.uk/BlueBook/P6.html#65020108), [P-73.1.1.2](https://iupac.qmul.ac.uk/BlueBook/P7.html#73010102), [P-73.1.2.1 and Table 7.4](https://iupac.qmul.ac.uk/BlueBook/P7.html#730102)
- **Anionic hydrazide.** P-72.2.2.2.4 uses the corresponding anionic parent hydride for anionic amide and hydrazide sites. The book prints `acetylazanide` as a PIN and `acetylhydrazine-1,1-diide` as a hydrazine precedent. The tested monoanion follows this method with hydrazin-1-ide, acetyl on N1, and methyl on N2. Both SMILES traversals preserve charge on the acyl-bound N. [P-72.2.2.2.4](https://iupac.qmul.ac.uk/BlueBook/P7.html#7202020204)
- **Cationic prefixes, zwitterions, and salts.** P-74.1.3 admits the anionic parent with a cationic prefix for ionic centres on different parents. P-73.6 distinguishes the systematic cationic-prefix method from permitted 'onio' alternatives. The public status is `Unverified(CationPrefix)` for the 'ammonio' forms shown. The zwitterion and chloride salt each preserve one positive and one negative atom, and the dication preserves both positive atoms. Neutral isothiocyanate and neutral N-methylacetohydrazide controls remain preferred and exact at both revisions. [P-73.6](https://iupac.qmul.ac.uk/BlueBook/P7.html#7306), [P-74.1.3](https://iupac.qmul.ac.uk/BlueBook/P7.html#740103)
- **Amidium suffix.** P-73.1.2.1 and Table 7.4 give 'amidium' and 'carboxamidium' for cationic amide characteristic groups. P-73.5.1.2 provides the multiplicative treatment of several cationic characteristic groups. The oxygen-protonated amide control retains an exact oxidaniumylidene parse-back. The tested S-protonated cumulenes return a naming error at the head. [P-73.1.2.1 and Table 7.4](https://iupac.qmul.ac.uk/BlueBook/P7.html#730102), [P-73.5.1.2](https://iupac.qmul.ac.uk/BlueBook/P7.html#73050102)

[Charge and attachment regressions through both naming entry points](https://github.com/Ameyanagi/nazuchem/blob/311c31b6b330e5cc0589cd96f56ee558ac7f42cc/crates/nazuchem-s2n/src/identity_validation.rs#L336-L455) · [Multiplied ring amidium regressions](https://github.com/Ameyanagi/nazuchem/blob/311c31b6b330e5cc0589cd96f56ee558ac7f42cc/crates/nazuchem-s2n/src/identity_validation.rs#L557-L628)

## Reproducing the before and after comparison

Save the Rust probe below as `probe.rs`. Run the commands from the repository with an up-to-date `main` and the tested commit available. The two worktrees keep the builds independent. The TSV's first line contains column headings.

```sh
work=$(mktemp -d)
git worktree add --detach "$work/base" c3b023bf1551027a93463847b5a27fc17f44976b
git worktree add --detach "$work/head" 2f70359adeecfceda4a3f83ec4035d1790b20e71
for side in base head; do
    target="$work/$side/target"
    CARGO_TARGET_DIR="$target" cargo build --manifest-path "$work/$side/Cargo.toml" -p nazuchem --lib --locked
    rustc --edition=2024 probe.rs --extern nazuchem="$target/debug/libnazuchem.rlib" \
        -L "dependency=$target/debug/deps" -o "$work/$side-probe"
    awk -F '\t' 'NR > 1 {print $3 "\t" $5}' docs/pr-191-cases.tsv | "$work/$side-probe" > "$work/$side-results.tsv"
done
```

```rust
use std::io::{self, BufRead};
use nazuchem::{Budget, name_to_smiles, s2n::verify_identity, smiles_to_name};

fn clean(text: &str) -> String {
    text.replace(['\t', '\r', '\n'], " ")
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let budget = Budget::default();
    for line in io::stdin().lock().lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let Some((tag, input)) = line.split_once('\t') else {
            continue;
        };
        if tag.ends_with("-parse") {
            match name_to_smiles(input) {
                Ok(smiles) => println!("{tag}\tparse\t{}", clean(&smiles)),
                Err(error) => println!("{tag}\tparse-error\t{}", clean(&error.to_string())),
            }
            continue;
        }
        match smiles_to_name(input) {
            Ok(name) => {
                let identity = match name_to_smiles(&name.text) {
                    Ok(parsed) => match verify_identity(input, &parsed, &budget) {
                        Ok(()) => "exact".to_string(),
                        Err(error) => format!("mismatch: {}", clean(&error)),
                    },
                    Err(error) => format!("parse-error: {}", clean(&error.to_string())),
                };
                println!("{tag}\t{:?}\t{}\t{identity}", name.status, clean(&name.text));
            }
            Err(error) => println!("{tag}\tname-error\t{}", clean(&error.to_string())),
        }
    }
    Ok(())
}
```
