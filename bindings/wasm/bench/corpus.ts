// SPDX-License-Identifier: NOASSERTION
// Origin: data
// Source: Nazuchem | bindings/wasm/bench/corpus.ts | npm/wasm size benchmark

/**
 * Frozen benchmark v1: 300 distinct successful name/SMILES pairs from
 * tests/fixtures/golden/s2n_n2s.jsonl, annotated with the corresponding s2n
 * bakeoff_corpus.tsv, rules_cases.tsv, or openclatura_tests.tsv source.
 * Selection and annotations are original Nazuchem work; the copied fixture
 * values retain their source fixture provenance and unresolved upstream terms.
 *
 * Selection was performed once, before variant measurements: require recorded
 * s2n Ok and n2s SUCCESS; cap names at 110 characters, SMILES at 90 characters,
 * atom tokens at 32, and ring-closure pairs at 3 to exclude large, highly cyclic
 * canonicalization time bombs. Group by source/category; rank each group by
 * SHA256("nazuchem-wasm-bench-v1\0" + smiles + "\0" + name); take round-robin
 * in ascending group order, skipping duplicate names and duplicate SMILES.
 * Source precedence for SMILES present in multiple files is bakeoff, rules,
 * openclatura. These literal inputs never change according to variant success
 * or timing. Golden line numbers and TSV line numbers are one-based.
 */
export interface CorpusCase {
  name: string;
  smiles: string;
  goldenLine: number;
  source: string;
  category: string;
}

export const corpus: readonly CorpusCase[] = [
  {
    "name": "2-aminoethan-1-ol",
    "smiles": "C(CO)N",
    "goldenLine": 48,
    "source": "bakeoff_corpus.tsv:49",
    "category": "alcohol_amine"
  },
  {
    "name": "1-methoxy-4-[(1E)-prop-1-en-1-yl]benzene",
    "smiles": "COc1ccc(/C=C/C)cc1",
    "goldenLine": 65,
    "source": "bakeoff_corpus.tsv:66",
    "category": "alkenyl_benzyl"
  },
  {
    "name": "4-(2-aminoethyl)benzene-1,2-diol",
    "smiles": "C1=CC(=C(C=C1CCN)O)O",
    "goldenLine": 12,
    "source": "bakeoff_corpus.tsv:13",
    "category": "chain_ring_prefix"
  },
  {
    "name": "6-methoxy-2-[(4-methoxy-3,5-dimethylpyridin-2-yl)methanesulfinyl]-1H-1,3-benzimidazole",
    "smiles": "CC1=CN=C(C(=C1OC)C)CS(=O)C2=NC3=C(N2)C=C(C=C3)OC",
    "goldenLine": 100,
    "source": "bakeoff_corpus.tsv:101",
    "category": "drug"
  },
  {
    "name": "3,7-dimethyl-3,7-dihydro-1H-purine-2,6-dione",
    "smiles": "CN1C=NC2=C1C(=O)NC(=O)N2C",
    "goldenLine": 83,
    "source": "bakeoff_corpus.tsv:84",
    "category": "heterocycle"
  },
  {
    "name": "2,2'-bipyridine",
    "smiles": "C1=CC=NC(=C1)C2=CC=CC=N2",
    "goldenLine": 188,
    "source": "bakeoff_corpus.tsv:189",
    "category": "mr_assembly"
  },
  {
    "name": "2,6,6-trimethylbicyclo[3.1.1]hept-2-ene",
    "smiles": "CC1=CCC2CC1C2(C)C",
    "goldenLine": 164,
    "source": "bakeoff_corpus.tsv:165",
    "category": "mr_bridged"
  },
  {
    "name": "9H-carbazole",
    "smiles": "C1=CC=C2C(=C1)C3=CC=CC=C3N2",
    "goldenLine": 180,
    "source": "bakeoff_corpus.tsv:181",
    "category": "mr_fused"
  },
  {
    "name": "cyclododecane",
    "smiles": "C1CCCCCCCCCCC1",
    "goldenLine": 201,
    "source": "bakeoff_corpus.tsv:202",
    "category": "mr_macro_other"
  },
  {
    "name": "spiro[4.5]decan-1-one",
    "smiles": "O=C1CCCC12CCCCC2",
    "goldenLine": 150,
    "source": "bakeoff_corpus.tsv:151",
    "category": "mr_spiro"
  },
  {
    "name": "(5R)-5-[(1S)-1,2-dihydroxyethyl]-3,4-dihydroxyfuran-2(5H)-one",
    "smiles": "C([C@@H]([C@@H]1C(=C(C(=O)O1)O)O)O)O",
    "goldenLine": 126,
    "source": "bakeoff_corpus.tsv:127",
    "category": "natural_product"
  },
  {
    "name": "(2E)-N-[(2-chloro-6-fluorophenyl)methyl]-3-(4-methoxyphenyl)-N-methylprop-2-enamide",
    "smiles": "CN(CC1=C(C=CC=C1Cl)F)C(=O)/C=C/C2=CC=C(C=C2)OC",
    "goldenLine": 624,
    "source": "bakeoff_corpus.tsv:625",
    "category": "pubchem_random"
  },
  {
    "name": "phenyl benzoate",
    "smiles": "C1=CC=C(C=C1)C(=O)OC2=CC=CC=C2",
    "goldenLine": 18,
    "source": "bakeoff_corpus.tsv:19",
    "category": "ring_ON_ester_amide"
  },
  {
    "name": "(2E)-but-2-ene",
    "smiles": "C/C=C/C",
    "goldenLine": 130,
    "source": "bakeoff_corpus.tsv:131",
    "category": "stereo"
  },
  {
    "name": "phenol",
    "smiles": "C1=CC=C(C=C1)O",
    "goldenLine": 32,
    "source": "bakeoff_corpus.tsv:33",
    "category": "suffix_on_ring"
  },
  {
    "name": "furo[2,3-b]thieno[3,2-e]pyridine",
    "smiles": "c1sc2nc3occc3cc2c1",
    "goldenLine": 731,
    "source": "openclatura_tests.tsv:29",
    "category": "data/fusion_pin_cases.json"
  },
  {
    "name": "furo[3,2-b]thieno[2,3-e]pyridine",
    "smiles": "O1C=CC2=NC3=C(C=C21)SC=C3",
    "goldenLine": 764,
    "source": "openclatura_tests.tsv:62",
    "category": "data/systematic_fusion_cases.json"
  },
  {
    "name": "hydrogen iodide",
    "smiles": "[IH]",
    "goldenLine": 873,
    "source": "openclatura_tests.tsv:174",
    "category": "exploratory/test_charge_representation_gaps.py"
  },
  {
    "name": "anisole",
    "smiles": "COc1ccccc1",
    "goldenLine": 850,
    "source": "openclatura_tests.tsv:149",
    "category": "fixtures/diverse_corpus.golden.json"
  },
  {
    "name": "5-methyldodecahydro-1H-cyclopenta[10]annulene",
    "smiles": "CC1CCCCCCC2CCCC2C1",
    "goldenLine": 876,
    "source": "openclatura_tests.tsv:177",
    "category": "integration/test_fusion_annulene_production.py"
  },
  {
    "name": "cyclopenta[c][1,2,5]oxadiazol-4-id-4-aminium",
    "smiles": "[NH3+][C-]1C=CC2=NON=C12",
    "goldenLine": 891,
    "source": "openclatura_tests.tsv:192",
    "category": "integration/test_fusion_charge_numbering.py"
  },
  {
    "name": "1H-imidazo[4,5-e][1,2,4]triazin-2-ium-1-ide",
    "smiles": "[N-]1[NH+]=CN=C2N=CN=C12",
    "goldenLine": 897,
    "source": "openclatura_tests.tsv:198",
    "category": "integration/test_fusion_charged_migration.py"
  },
  {
    "name": "octahydropyrrolo[3,4-b]pyrrol-5-ium",
    "smiles": "N1CCC2C1C[NH2+]C2",
    "goldenLine": 911,
    "source": "openclatura_tests.tsv:212",
    "category": "integration/test_fusion_saturated_composition.py"
  },
  {
    "name": "N,N-dimethylmethanamine N-selenide",
    "smiles": "C[N+]([Se-])(C)C",
    "goldenLine": 927,
    "source": "openclatura_tests.tsv:228",
    "category": "integration/test_human_descriptor.py"
  },
  {
    "name": "azulene",
    "smiles": "C1=CC=C2C=CC=CC=C12",
    "goldenLine": 942,
    "source": "openclatura_tests.tsv:243",
    "category": "integration/test_systematic_fusion.py"
  },
  {
    "name": "methyl acetate",
    "smiles": "COC(=O)C",
    "goldenLine": 957,
    "source": "openclatura_tests.tsv:258",
    "category": "roundtrip/test_acids.py"
  },
  {
    "name": "propan-2-ol",
    "smiles": "CC(O)C",
    "goldenLine": 959,
    "source": "openclatura_tests.tsv:260",
    "category": "roundtrip/test_alcohols.py"
  },
  {
    "name": "cyclohexa-1,3-diene",
    "smiles": "C1=CC=CCC1",
    "goldenLine": 962,
    "source": "openclatura_tests.tsv:263",
    "category": "roundtrip/test_alkenes_alkynes.py"
  },
  {
    "name": "2-cyclohexylethan-1-ol",
    "smiles": "OCCC1CCCCC1",
    "goldenLine": 963,
    "source": "openclatura_tests.tsv:264",
    "category": "roundtrip/test_api.py"
  },
  {
    "name": "1-(2,3-dimethylbutyl)-2-methylcyclopropane",
    "smiles": "CC1CC1CC(C(C)C)C",
    "goldenLine": 966,
    "source": "openclatura_tests.tsv:267",
    "category": "roundtrip/test_bridges.py"
  },
  {
    "name": "1,2-dimethylcyclopropane",
    "smiles": "CC1CC1C",
    "goldenLine": 981,
    "source": "openclatura_tests.tsv:285",
    "category": "roundtrip/test_composite_bridges.py"
  },
  {
    "name": "2,3-dihydro-1H-inden-1-one",
    "smiles": "O=C1CCc2ccccc21",
    "goldenLine": 1023,
    "source": "openclatura_tests.tsv:327",
    "category": "roundtrip/test_fused_hydrocarbon_retained.py"
  },
  {
    "name": "dithieno[2,3-b:3',2'-e]pyridine",
    "smiles": "c1sc2nc3sccc3cc2c1",
    "goldenLine": 1025,
    "source": "openclatura_tests.tsv:329",
    "category": "roundtrip/test_fused_multiplicity.py"
  },
  {
    "name": "1H-inden-4-ol",
    "smiles": "C1=Cc2c(O)cccc2C1",
    "goldenLine": 1074,
    "source": "openclatura_tests.tsv:380",
    "category": "roundtrip/test_fused_parents.py"
  },
  {
    "name": "5,6-dihydroimidazo[2,1-b][1,3]thiazole",
    "smiles": "C1=CSC2=NCCN12",
    "goldenLine": 1103,
    "source": "openclatura_tests.tsv:409",
    "category": "roundtrip/test_fusion_hydrogenation_regressions.py"
  },
  {
    "name": "hydrogen fluoride",
    "smiles": "F",
    "goldenLine": 1110,
    "source": "openclatura_tests.tsv:416",
    "category": "roundtrip/test_fusion_polycyclic_carbon_h.py"
  },
  {
    "name": "2-chloroethan-1-ol",
    "smiles": "ClCCO",
    "goldenLine": 1112,
    "source": "openclatura_tests.tsv:418",
    "category": "roundtrip/test_haloalkanes.py"
  },
  {
    "name": "2-(methylphosphanyl)ethan-1-ol",
    "smiles": "C[PH]CCO",
    "goldenLine": 1116,
    "source": "openclatura_tests.tsv:422",
    "category": "roundtrip/test_heteroatom_substituents.py"
  },
  {
    "name": "1,4-oxathiane",
    "smiles": "O1CCSCC1",
    "goldenLine": 1126,
    "source": "openclatura_tests.tsv:434",
    "category": "roundtrip/test_heterocycles.py"
  },
  {
    "name": "1,2,3,4-tetrahydronaphthalen-1-imine",
    "smiles": "C1CCC2=CC=CC=C2C1=N",
    "goldenLine": 1164,
    "source": "openclatura_tests.tsv:472",
    "category": "roundtrip/test_imines.py"
  },
  {
    "name": "2-hydroxyethanal",
    "smiles": "O=CCO",
    "goldenLine": 1166,
    "source": "openclatura_tests.tsv:474",
    "category": "roundtrip/test_ketones.py"
  },
  {
    "name": "(chlorofluoromethyl)tri(ethyl)silane",
    "smiles": "CC[Si](CC)(CC)C(F)Cl",
    "goldenLine": 1168,
    "source": "openclatura_tests.tsv:476",
    "category": "roundtrip/test_ligand_parentheses.py"
  },
  {
    "name": "methanaminium",
    "smiles": "C[NH3+]",
    "goldenLine": 1188,
    "source": "openclatura_tests.tsv:496",
    "category": "roundtrip/test_nitro.py"
  },
  {
    "name": "2-methylbutane",
    "smiles": "CCC(C)C",
    "goldenLine": 1193,
    "source": "openclatura_tests.tsv:501",
    "category": "roundtrip/test_numbering.py"
  },
  {
    "name": "ethyl dihydrogen phosphate",
    "smiles": "CCOP(=O)(O)O",
    "goldenLine": 1201,
    "source": "openclatura_tests.tsv:509",
    "category": "roundtrip/test_oxyacids.py"
  },
  {
    "name": "fluoro(methoxy)(methyl)(phenyl)-lambda5-phosphane",
    "smiles": "CP(F)(OC)c1ccccc1",
    "goldenLine": 1204,
    "source": "openclatura_tests.tsv:512",
    "category": "roundtrip/test_p13_additive.py"
  },
  {
    "name": "phenylacetic acid",
    "smiles": "c1ccccc1CC(=O)O",
    "goldenLine": 1215,
    "source": "openclatura_tests.tsv:523",
    "category": "roundtrip/test_p13_conjunctive.py"
  },
  {
    "name": "1-[3,5-bis(tert-butoxy)phenyl]ethan-1-ol",
    "smiles": "CC(O)c1cc(OC(C)(C)C)cc(OC(C)(C)C)c1",
    "goldenLine": 1225,
    "source": "openclatura_tests.tsv:533",
    "category": "roundtrip/test_p13_multiplicative.py"
  },
  {
    "name": "piperidine",
    "smiles": "N1CCCCC1",
    "goldenLine": 1226,
    "source": "openclatura_tests.tsv:534",
    "category": "roundtrip/test_p13_operations.py"
  },
  {
    "name": "3,6-diazabicyclo[3.2.1]octane",
    "smiles": "N1CC2CNCC1C2",
    "goldenLine": 1230,
    "source": "openclatura_tests.tsv:538",
    "category": "roundtrip/test_p13_replacement.py"
  },
  {
    "name": "ethanetellurol",
    "smiles": "CC[TeH]",
    "goldenLine": 1236,
    "source": "openclatura_tests.tsv:544",
    "category": "roundtrip/test_p13_substitutive.py"
  },
  {
    "name": "cyclohexene",
    "smiles": "C1CCCC=C1",
    "goldenLine": 1239,
    "source": "openclatura_tests.tsv:547",
    "category": "roundtrip/test_p13_subtractive.py"
  },
  {
    "name": "2-oxabicyclo[2.2.1]heptane",
    "smiles": "O1CC2CCC1C2",
    "goldenLine": 1244,
    "source": "openclatura_tests.tsv:552",
    "category": "roundtrip/test_polycycles.py"
  },
  {
    "name": "1,2,3,4,4a,8a-hexahydronaphthalene",
    "smiles": "C1=CC2CCCCC2C=C1",
    "goldenLine": 1249,
    "source": "openclatura_tests.tsv:560",
    "category": "roundtrip/test_polycyclic_descriptors.py"
  },
  {
    "name": "N-(1-{[(2,4-dimethylphenyl)methyl]amino}-3-methyl-1-oxobutan-2-yl)-3,5-dimethoxybenzamide",
    "smiles": "COc1cc(OC)cc(C(=O)NC(C(=O)NCc2ccc(C)cc2C)C(C)C)c1",
    "goldenLine": 1260,
    "source": "openclatura_tests.tsv:571",
    "category": "roundtrip/test_polyfunctional_fragments.py"
  },
  {
    "name": "(furo[2,3-b]pyridin-2-yl)boronic acid",
    "smiles": "OB(O)c1oc2ncccc2c1",
    "goldenLine": 1281,
    "source": "openclatura_tests.tsv:592",
    "category": "roundtrip/test_principal_group_breadth.py"
  },
  {
    "name": "ammonia",
    "smiles": "[NH3]",
    "goldenLine": 1296,
    "source": "openclatura_tests.tsv:607",
    "category": "roundtrip/test_replacement_parents.py"
  },
  {
    "name": "ethanol",
    "smiles": "C(C)O",
    "goldenLine": 1319,
    "source": "openclatura_tests.tsv:630",
    "category": "roundtrip/test_retained_fused_production.py"
  },
  {
    "name": "dimethylcarbamoyl dihydrogen phosphate",
    "smiles": "CN(C)C(=O)OP(=O)(O)O",
    "goldenLine": 1323,
    "source": "openclatura_tests.tsv:634",
    "category": "roundtrip/test_role_models.py"
  },
  {
    "name": "(R)-bromochlorofluoromethane",
    "smiles": "F[C@H](Cl)Br",
    "goldenLine": 1332,
    "source": "openclatura_tests.tsv:643",
    "category": "roundtrip/test_stereo.py"
  },
  {
    "name": "2-chloropropan-1-ol",
    "smiles": "CC(Cl)CO",
    "goldenLine": 1333,
    "source": "openclatura_tests.tsv:645",
    "category": "roundtrip/test_substituents.py"
  },
  {
    "name": "propanethioic O-acid",
    "smiles": "CCC(=S)O",
    "goldenLine": 1336,
    "source": "openclatura_tests.tsv:648",
    "category": "unit/test_all_chalcogen_acyl_regressions.py"
  },
  {
    "name": "N-[10-(methylamino)phenanthren-9-yl]-2-oxoacetamide",
    "smiles": "CNc1c(NC(=O)C=O)c2ccccc2c2ccccc12",
    "goldenLine": 1349,
    "source": "openclatura_tests.tsv:661",
    "category": "unit/test_analysis.py"
  },
  {
    "name": "1,3-dimethyl-3,7-dihydro-1H-purine-2,6-dione",
    "smiles": "Cn1c(=O)c2[nH]cnc2n(C)c1=O",
    "goldenLine": 1360,
    "source": "openclatura_tests.tsv:674",
    "category": "unit/test_audit_reconstruction.py"
  },
  {
    "name": "methyl 4-(N'-hydroxycarbamimidoyl)benzoate",
    "smiles": "COC(=O)c1ccc(cc1)C(=NO)N",
    "goldenLine": 1367,
    "source": "openclatura_tests.tsv:681",
    "category": "unit/test_carbamimidoyl_stereo.py"
  },
  {
    "name": "Te-methyl ethanetelluroate",
    "smiles": "CC(=O)[Te]C",
    "goldenLine": 1376,
    "source": "openclatura_tests.tsv:690",
    "category": "unit/test_chalcogen_citation_compatibility.py"
  },
  {
    "name": "ethaneselenothioic S-acid",
    "smiles": "CC(=[Se])[SH]",
    "goldenLine": 1400,
    "source": "openclatura_tests.tsv:714",
    "category": "unit/test_chalcogen_functional_groups.py"
  },
  {
    "name": "acetyl isoselenocyanate",
    "smiles": "CC(=O)N=C=[Se]",
    "goldenLine": 1514,
    "source": "openclatura_tests.tsv:828",
    "category": "unit/test_chalcogen_structural_roles.py"
  },
  {
    "name": "6,7-dihydro-5H-cyclopenta[b]pyridin-5-one",
    "smiles": "C1Cc2ncccc2C1=O",
    "goldenLine": 1518,
    "source": "openclatura_tests.tsv:832",
    "category": "unit/test_fused_oxo_added_hydrogen_tier.py"
  },
  {
    "name": "water",
    "smiles": "O",
    "goldenLine": 1027,
    "source": "openclatura_tests.tsv:331",
    "category": "unit/test_fusion_bridge_hydrogen_state.py"
  },
  {
    "name": "2-(methylsulfanyl)benzo[h]quinazolin-4-amine",
    "smiles": "CSc1nc(N)c2ccc3ccccc3c2n1",
    "goldenLine": 1523,
    "source": "openclatura_tests.tsv:837",
    "category": "unit/test_fusion_candidate_budget.py"
  },
  {
    "name": "7-methyl-6,7-dihydro-5H-cyclopenta[b]pyridin-5-one",
    "smiles": "C1C(C)c2ncccc2C1=O",
    "goldenLine": 1527,
    "source": "openclatura_tests.tsv:841",
    "category": "unit/test_fusion_carbon_h_oxo_composition.py"
  },
  {
    "name": "4-methyl-4,7-dihydro[1,3]oxazolo[5,4-d]pyrimidine",
    "smiles": "C1N=CN(C)C2=C1N=CO2",
    "goldenLine": 1530,
    "source": "openclatura_tests.tsv:844",
    "category": "unit/test_fusion_heteroatom_pi_redistribution.py"
  },
  {
    "name": "furo[3,2-d][1,3]oxazol-2(1H)-one",
    "smiles": "O=C1NC2=C(OC=C2)O1",
    "goldenLine": 1533,
    "source": "openclatura_tests.tsv:847",
    "category": "unit/test_fusion_intrinsic_h_suffix.py"
  },
  {
    "name": "4-ethyl-5,6-dihydro-1H-pyrrolo[3,2-c]pyrazole",
    "smiles": "CCN1CCC2=C1C=NN2",
    "goldenLine": 1535,
    "source": "openclatura_tests.tsv:849",
    "category": "unit/test_fusion_mancude_aromatic_h.py"
  },
  {
    "name": "1,2,3,4-tetrahydropyrrolo[3,2-b]pyrrole",
    "smiles": "c1cc2c([nH]1)CCN2",
    "goldenLine": 1538,
    "source": "openclatura_tests.tsv:852",
    "category": "unit/test_fusion_mixed_nh_donors.py"
  },
  {
    "name": "1,4-hexanonaphthalene",
    "smiles": "C12=CC=C(C3=CC=CC=C13)CCCCCC2",
    "goldenLine": 1544,
    "source": "openclatura_tests.tsv:858",
    "category": "unit/test_fusion_wrappers.py"
  },
  {
    "name": "N-(methanesulfonyl)benzaldehyde hydrazone",
    "smiles": "CS(=O)(=O)NN=Cc1ccccc1",
    "goldenLine": 1554,
    "source": "openclatura_tests.tsv:868",
    "category": "unit/test_hydrazone_stereo_scope.py"
  },
  {
    "name": "4-carbamoylbenzoic acid",
    "smiles": "NC(=O)c1ccc(C(=O)O)cc1",
    "goldenLine": 1555,
    "source": "openclatura_tests.tsv:869",
    "category": "unit/test_parent_hydride_handoff.py"
  },
  {
    "name": "1,5-dimethylpentaoxidane",
    "smiles": "COOOOOC",
    "goldenLine": 1558,
    "source": "openclatura_tests.tsv:872",
    "category": "unit/test_polyoxide_chains.py"
  },
  {
    "name": "3-fluoro-7-methoxydibenzofuran",
    "smiles": "COc1ccc2c(c1)oc1cc(F)ccc12",
    "goldenLine": 1576,
    "source": "openclatura_tests.tsv:890",
    "category": "unit/test_retained_parent_precedence_and_order.py"
  },
  {
    "name": "sulfane",
    "smiles": "S",
    "goldenLine": 1028,
    "source": "openclatura_tests.tsv:332",
    "category": "unit/test_spiro_component_scope.py"
  },
  {
    "name": "2,2,5-trimethylhexane",
    "smiles": "C(C)(C)(C)CCC(C)C",
    "goldenLine": 1645,
    "source": "rules_cases.tsv:90",
    "category": "rules"
  },
  {
    "name": "N,N-diethylethanamine",
    "smiles": "CCN(CC)CC",
    "goldenLine": 49,
    "source": "bakeoff_corpus.tsv:50",
    "category": "alcohol_amine"
  },
  {
    "name": "(1E)-2-phenylethenylbenzene",
    "smiles": "C1=CC=C(C=C1)/C=C/C2=CC=CC=C2",
    "goldenLine": 59,
    "source": "bakeoff_corpus.tsv:60",
    "category": "alkenyl_benzyl"
  },
  {
    "name": "(hydroxy)(phenyl)acetic acid",
    "smiles": "C1=CC=C(C=C1)C(C(=O)O)O",
    "goldenLine": 4,
    "source": "bakeoff_corpus.tsv:5",
    "category": "chain_ring_prefix"
  },
  {
    "name": "(2S)-2-(6-methoxynaphthalen-2-yl)propanoic acid",
    "smiles": "C[C@@H](C1=CC2=C(C=C1)C=C(C=C2)OC)C(=O)O",
    "goldenLine": 101,
    "source": "bakeoff_corpus.tsv:102",
    "category": "drug"
  },
  {
    "name": "morpholine",
    "smiles": "C1COCCN1",
    "goldenLine": 67,
    "source": "bakeoff_corpus.tsv:68",
    "category": "heterocycle"
  },
  {
    "name": "1^1,2^1:2^4,3^1-terphenyl",
    "smiles": "C1=CC=C(C=C1)C2=CC=C(C=C2)C3=CC=CC=C3",
    "goldenLine": 190,
    "source": "bakeoff_corpus.tsv:191",
    "category": "mr_assembly"
  },
  {
    "name": "(1R,4S)-bicyclo[2.2.1]hept-2-ene",
    "smiles": "C1C[C@H]2C[C@@H]1C=C2",
    "goldenLine": 154,
    "source": "bakeoff_corpus.tsv:155",
    "category": "mr_bridged"
  },
  {
    "name": "2-methylnaphthalene",
    "smiles": "CC1=CC2=CC=CC=C2C=C1",
    "goldenLine": 168,
    "source": "bakeoff_corpus.tsv:169",
    "category": "mr_fused"
  },
  {
    "name": "1,4(1,4)-dibenzenacyclohexaphane",
    "smiles": "C1CC2=CC=C(CCC3=CC=C1C=C3)C=C2",
    "goldenLine": 198,
    "source": "bakeoff_corpus.tsv:199",
    "category": "mr_macro_other"
  },
  {
    "name": "spiro[4.5]decane",
    "smiles": "C1CCC2(CC1)CCCC2",
    "goldenLine": 146,
    "source": "bakeoff_corpus.tsv:147",
    "category": "mr_spiro"
  },
  {
    "name": "(5R)-2-methyl-5-(prop-1-en-2-yl)cyclohex-2-en-1-one",
    "smiles": "CC1=CC[C@H](CC1=O)C(=C)C",
    "goldenLine": 120,
    "source": "bakeoff_corpus.tsv:121",
    "category": "natural_product"
  },
  {
    "name": "N-(3,4-dimethoxyphenyl)-N-methyl-2-[2-(thiophen-2-yl)-1,3-thiazol-4-yl]acetamide",
    "smiles": "CN(C1=CC(=C(C=C1)OC)OC)C(=O)CC2=CSC(=N2)C3=CC=CS3",
    "goldenLine": 287,
    "source": "bakeoff_corpus.tsv:288",
    "category": "pubchem_random"
  },
  {
    "name": "N,N-dimethylbenzamide",
    "smiles": "CN(C)C(=O)C1=CC=CC=C1",
    "goldenLine": 22,
    "source": "bakeoff_corpus.tsv:23",
    "category": "ring_ON_ester_amide"
  },
  {
    "name": "(2S)-2-aminopropanoic acid",
    "smiles": "C[C@@H](C(=O)O)N",
    "goldenLine": 135,
    "source": "bakeoff_corpus.tsv:136",
    "category": "stereo"
  },
  {
    "name": "aniline",
    "smiles": "C1=CC=C(C=C1)N",
    "goldenLine": 33,
    "source": "bakeoff_corpus.tsv:34",
    "category": "suffix_on_ring"
  },
  {
    "name": "thieno[2,3-c]pyridine",
    "smiles": "S1C=CC=2C1=CN=CC2",
    "goldenLine": 749,
    "source": "openclatura_tests.tsv:47",
    "category": "data/fusion_pin_cases.json"
  },
  {
    "name": "selenopheno[3,2-b]selenophene",
    "smiles": "[Se]1C2=C(C=C1)[Se]C=C2",
    "goldenLine": 763,
    "source": "openclatura_tests.tsv:61",
    "category": "data/systematic_fusion_cases.json"
  },
  {
    "name": "hydrogen bromide",
    "smiles": "[BrH]",
    "goldenLine": 872,
    "source": "openclatura_tests.tsv:173",
    "category": "exploratory/test_charge_representation_gaps.py"
  },
  {
    "name": "propene",
    "smiles": "CC=C",
    "goldenLine": 779,
    "source": "openclatura_tests.tsv:77",
    "category": "fixtures/diverse_corpus.golden.json"
  },
  {
    "name": "dodecahydro-1H-cyclopenta[10]annulen-5-ol",
    "smiles": "OC1CCCCCCC2CCCC2C1",
    "goldenLine": 877,
    "source": "openclatura_tests.tsv:178",
    "category": "integration/test_fusion_annulene_production.py"
  },
  {
    "name": "N-methylcyclopenta[c][1,2,5]oxadiazol-4-id-4-aminium",
    "smiles": "[NH2+](C)[C-]1C=CC2=NON=C12",
    "goldenLine": 892,
    "source": "openclatura_tests.tsv:193",
    "category": "integration/test_fusion_charge_numbering.py"
  },
  {
    "name": "N-methylcyclopenta[d][1,2,3]oxadiazol-4-id-4-aminium",
    "smiles": "[NH2+](C)[C-]1C=CC2=C1N=NO2",
    "goldenLine": 904,
    "source": "openclatura_tests.tsv:205",
    "category": "integration/test_fusion_charged_migration.py"
  },
  {
    "name": "decahydropyrido[4,3-b]pyridine-1,6-diium",
    "smiles": "[NH2+]1CCCC2C1CC[NH2+]C2",
    "goldenLine": 920,
    "source": "openclatura_tests.tsv:221",
    "category": "integration/test_fusion_saturated_composition.py"
  },
  {
    "name": "2-methylthieno[2,3-b]furan",
    "smiles": "CC1=CC2=C(O1)SC=C2",
    "goldenLine": 929,
    "source": "openclatura_tests.tsv:230",
    "category": "integration/test_systematic_fusion.py"
  },
  {
    "name": "phenyl acetate",
    "smiles": "CC(=O)Oc1ccccc1",
    "goldenLine": 958,
    "source": "openclatura_tests.tsv:259",
    "category": "roundtrip/test_acids.py"
  },
  {
    "name": "pent-1-en-4-yne",
    "smiles": "C=CCC#C",
    "goldenLine": 961,
    "source": "openclatura_tests.tsv:262",
    "category": "roundtrip/test_alkenes_alkynes.py"
  },
  {
    "name": "1,3-oxazole",
    "smiles": "c1cocn1",
    "goldenLine": 970,
    "source": "openclatura_tests.tsv:271",
    "category": "roundtrip/test_bridges.py"
  },
  {
    "name": "3-ethyl-2-methylhexane",
    "smiles": "CCCC(C(C)C)CC",
    "goldenLine": 983,
    "source": "openclatura_tests.tsv:287",
    "category": "roundtrip/test_composite_bridges.py"
  },
  {
    "name": "2,3-dihydro-1H-indole",
    "smiles": "c1ccc2c(c1)CCN2",
    "goldenLine": 1017,
    "source": "openclatura_tests.tsv:321",
    "category": "roundtrip/test_fused_hydrocarbon_retained.py"
  },
  {
    "name": "1H-indazol-5-ol",
    "smiles": "Oc1ccc2[nH]ncc2c1",
    "goldenLine": 1078,
    "source": "openclatura_tests.tsv:384",
    "category": "roundtrip/test_fused_parents.py"
  },
  {
    "name": "2-methyl-2,3-dihydropyrrolo[2,1-b]quinazolin-9(1H)-one",
    "smiles": "O=c1c2ccccc2nc2n1CC(C)C2",
    "goldenLine": 1104,
    "source": "openclatura_tests.tsv:410",
    "category": "roundtrip/test_fusion_hydrogenation_regressions.py"
  },
  {
    "name": "chloroethane",
    "smiles": "CCCl",
    "goldenLine": 1111,
    "source": "openclatura_tests.tsv:417",
    "category": "roundtrip/test_haloalkanes.py"
  },
  {
    "name": "2-(dimethylamino)ethan-1-ol",
    "smiles": "CN(C)CCO",
    "goldenLine": 1115,
    "source": "openclatura_tests.tsv:421",
    "category": "roundtrip/test_heteroatom_substituents.py"
  },
  {
    "name": "oxiran-2-one",
    "smiles": "O=C1CO1",
    "goldenLine": 1152,
    "source": "openclatura_tests.tsv:460",
    "category": "roundtrip/test_heterocycles.py"
  },
  {
    "name": "N-methylethanimine",
    "smiles": "CC=NC",
    "goldenLine": 1160,
    "source": "openclatura_tests.tsv:468",
    "category": "roundtrip/test_imines.py"
  },
  {
    "name": "butan-2-one",
    "smiles": "CCC(=O)C",
    "goldenLine": 1165,
    "source": "openclatura_tests.tsv:473",
    "category": "roundtrip/test_ketones.py"
  },
  {
    "name": "bis(tert-butoxy)di(methoxy)silane",
    "smiles": "CO[Si](OC)(OC(C)(C)C)OC(C)(C)C",
    "goldenLine": 1173,
    "source": "openclatura_tests.tsv:481",
    "category": "roundtrip/test_ligand_parentheses.py"
  },
  {
    "name": "nitromethane",
    "smiles": "C[N+](=O)[O-]",
    "goldenLine": 1189,
    "source": "openclatura_tests.tsv:497",
    "category": "roundtrip/test_nitro.py"
  },
  {
    "name": "butan-2-ol",
    "smiles": "CC(O)CC",
    "goldenLine": 1192,
    "source": "openclatura_tests.tsv:500",
    "category": "roundtrip/test_numbering.py"
  },
  {
    "name": "methyl nitrate",
    "smiles": "CO[N+](=O)[O-]",
    "goldenLine": 1196,
    "source": "openclatura_tests.tsv:504",
    "category": "roundtrip/test_oxyacids.py"
  },
  {
    "name": "1,2,3,6-tetrahydropyridine",
    "smiles": "C1=CCNCC1",
    "goldenLine": 1203,
    "source": "openclatura_tests.tsv:511",
    "category": "roundtrip/test_p13_additive.py"
  },
  {
    "name": "2-([1,1'-biphenyl]-2-yl)ethan-1-ol",
    "smiles": "OCCc1ccccc1-c1ccccc1",
    "goldenLine": 1212,
    "source": "openclatura_tests.tsv:520",
    "category": "roundtrip/test_p13_conjunctive.py"
  },
  {
    "name": "2,4-bis(tert-butoxy)-3,5,6-trichloropyridine",
    "smiles": "CC(C)(C)Oc1nc(Cl)c(Cl)c(OC(C)(C)C)c1Cl",
    "goldenLine": 1224,
    "source": "openclatura_tests.tsv:532",
    "category": "roundtrip/test_p13_multiplicative.py"
  },
  {
    "name": "1,2-dichlorodisilene",
    "smiles": "Cl[SiH]=[SiH]Cl",
    "goldenLine": 1228,
    "source": "openclatura_tests.tsv:536",
    "category": "roundtrip/test_p13_replacement.py"
  },
  {
    "name": "2-bromo-4-methyl-N-{(2R)-4-oxo-4-[(propan-2-yl)amino]butan-2-yl}benzamide",
    "smiles": "Cc1ccc(C(=O)N[C@H](C)CC(=O)NC(C)C)c(Br)c1",
    "goldenLine": 1237,
    "source": "openclatura_tests.tsv:545",
    "category": "roundtrip/test_p13_substitutive.py"
  },
  {
    "name": "cyclohexyne",
    "smiles": "C1CC#CCC1",
    "goldenLine": 1240,
    "source": "openclatura_tests.tsv:548",
    "category": "roundtrip/test_p13_subtractive.py"
  },
  {
    "name": "3a,4,5,6,7,7a-hexahydro-1H-indene",
    "smiles": "C1=CC2CCCCC2C1",
    "goldenLine": 1241,
    "source": "openclatura_tests.tsv:549",
    "category": "roundtrip/test_polycycles.py"
  },
  {
    "name": "octahydro-2H-2lambda^4-benzothiophen-2-one",
    "smiles": "O=S1CC2CCCCC2C1",
    "goldenLine": 1247,
    "source": "openclatura_tests.tsv:558",
    "category": "roundtrip/test_polycyclic_descriptors.py"
  },
  {
    "name": "2-(3-aminoazetidin-1-yl)-2-phenylethan-1-ol",
    "smiles": "NC1CN(C(CO)c2ccccc2)C1",
    "goldenLine": 1263,
    "source": "openclatura_tests.tsv:574",
    "category": "roundtrip/test_polyfunctional_fragments.py"
  },
  {
    "name": "ethaneselenonic acid",
    "smiles": "CC[Se](=O)(=O)O",
    "goldenLine": 1268,
    "source": "openclatura_tests.tsv:579",
    "category": "roundtrip/test_principal_group_breadth.py"
  },
  {
    "name": "triaz-1-ene",
    "smiles": "NN=N",
    "goldenLine": 1304,
    "source": "openclatura_tests.tsv:615",
    "category": "roundtrip/test_replacement_parents.py"
  },
  {
    "name": "1,8-phenanthroline-7,9(8H,10H)-dione",
    "smiles": "N1=CC=CC2=CC=C3C(NC(CC3=C12)=O)=O",
    "goldenLine": 1317,
    "source": "openclatura_tests.tsv:628",
    "category": "roundtrip/test_retained_fused_production.py"
  },
  {
    "name": "1-iminoselenophene",
    "smiles": "N=[Se]1C=CC=C1",
    "goldenLine": 1327,
    "source": "openclatura_tests.tsv:638",
    "category": "roundtrip/test_role_models.py"
  },
  {
    "name": "1-ethenylcyclohex-1-ene",
    "smiles": "C=CC1=CCCCC1",
    "goldenLine": 1334,
    "source": "openclatura_tests.tsv:646",
    "category": "roundtrip/test_substituents.py"
  },
  {
    "name": "propanethioic anhydride",
    "smiles": "CCC(=S)OC(=S)CC",
    "goldenLine": 1340,
    "source": "openclatura_tests.tsv:652",
    "category": "unit/test_all_chalcogen_acyl_regressions.py"
  },
  {
    "name": "phosphoric acid",
    "smiles": "OP(=O)(O)O",
    "goldenLine": 1285,
    "source": "openclatura_tests.tsv:596",
    "category": "unit/test_analysis.py"
  },
  {
    "name": "(1r,3r)-3-(methanesulfonyl)cyclobutan-1-amine",
    "smiles": "CS(=O)(=O)[C@H]1C[C@H](N)C1",
    "goldenLine": 1359,
    "source": "openclatura_tests.tsv:673",
    "category": "unit/test_audit_reconstruction.py"
  },
  {
    "name": "(E)-N'-hydroxybenzenecarboximidamide",
    "smiles": "N/C(=N/O)c1ccccc1",
    "goldenLine": 1368,
    "source": "openclatura_tests.tsv:682",
    "category": "unit/test_carbamimidoyl_stereo.py"
  },
  {
    "name": "5-[(propan-2-yl)sulfanylcarbonyl]-1,2-oxazole-3-carboxylic acid",
    "smiles": "CC(C)SC(=O)c1cc(C(=O)O)no1",
    "goldenLine": 1371,
    "source": "openclatura_tests.tsv:685",
    "category": "unit/test_chalcogen_citation_compatibility.py"
  },
  {
    "name": "N,N-dimethylethanamine N-selenide",
    "smiles": "CC[N+](C)(C)[Se-]",
    "goldenLine": 1472,
    "source": "openclatura_tests.tsv:786",
    "category": "unit/test_chalcogen_functional_groups.py"
  },
  {
    "name": "butanoyl cyanide",
    "smiles": "CCCC(=O)C#N",
    "goldenLine": 1510,
    "source": "openclatura_tests.tsv:824",
    "category": "unit/test_chalcogen_structural_roles.py"
  },
  {
    "name": "6-methyl-6,7-dihydro-5H-cyclopenta[b]pyridin-5-one",
    "smiles": "CC1Cc2ncccc2C1=O",
    "goldenLine": 1519,
    "source": "openclatura_tests.tsv:833",
    "category": "unit/test_fused_oxo_added_hydrogen_tier.py"
  },
  {
    "name": "methyl 1-[(4-fluorophenyl)methyl]imidazo[4,5-c]pyridine-6-carboxylate",
    "smiles": "COC(=O)c1cc2c(cn1)ncn2Cc1ccc(F)cc1",
    "goldenLine": 1526,
    "source": "openclatura_tests.tsv:840",
    "category": "unit/test_fusion_candidate_budget.py"
  },
  {
    "name": "6,7-dihydro-5H-cyclopenta[b]pyrazin-5-one",
    "smiles": "C1Cc2nccnc2C1=O",
    "goldenLine": 1529,
    "source": "openclatura_tests.tsv:843",
    "category": "unit/test_fusion_carbon_h_oxo_composition.py"
  },
  {
    "name": "cyclopenta[b]pyrrol-6(1H)-one",
    "smiles": "O=C1C=CC2=C1NC=C2",
    "goldenLine": 1531,
    "source": "openclatura_tests.tsv:845",
    "category": "unit/test_fusion_intrinsic_h_suffix.py"
  },
  {
    "name": "4,6-dimethyl-5,6-dihydro-1H-pyrrolo[3,2-c]pyrazole",
    "smiles": "CN1CC(C)C2=C1C=NN2",
    "goldenLine": 1536,
    "source": "openclatura_tests.tsv:850",
    "category": "unit/test_fusion_mancude_aromatic_h.py"
  },
  {
    "name": "1,2,3,5-tetrahydropyrrolo[3,4-b]pyrrole",
    "smiles": "C1CC2=CNC=C2N1",
    "goldenLine": 1539,
    "source": "openclatura_tests.tsv:853",
    "category": "unit/test_fusion_mixed_nh_donors.py"
  },
  {
    "name": "1,4-prop[1]enonaphthalene",
    "smiles": "C12=CC=C(C3=CC=CC=C13)C=CC2",
    "goldenLine": 1547,
    "source": "openclatura_tests.tsv:861",
    "category": "unit/test_fusion_wrappers.py"
  },
  {
    "name": "(Z)-N-(methanesulfonyl)acetaldehyde hydrazone",
    "smiles": "CS(=O)(=O)N/N=C\\C",
    "goldenLine": 1553,
    "source": "openclatura_tests.tsv:867",
    "category": "unit/test_hydrazone_stereo_scope.py"
  },
  {
    "name": "quinoline",
    "smiles": "c1ccc2ncccc2c1",
    "goldenLine": 972,
    "source": "openclatura_tests.tsv:275",
    "category": "unit/test_parent_hydride_handoff.py"
  },
  {
    "name": "1-hydroperoxyhydrazine",
    "smiles": "OONN",
    "goldenLine": 1567,
    "source": "openclatura_tests.tsv:881",
    "category": "unit/test_polyoxide_chains.py"
  },
  {
    "name": "2-fluoro-7-methoxy-9H-xanthene",
    "smiles": "COc1ccc2c(c1)Cc1cc(F)ccc1O2",
    "goldenLine": 1573,
    "source": "openclatura_tests.tsv:887",
    "category": "unit/test_retained_parent_precedence_and_order.py"
  },
  {
    "name": "N-methylacetamide",
    "smiles": "CC(=O)NC",
    "goldenLine": 813,
    "source": "rules_cases.tsv:76",
    "category": "rules"
  },
  {
    "name": "2-phenylethan-1-ol",
    "smiles": "C1=CC=C(C=C1)CCO",
    "goldenLine": 43,
    "source": "bakeoff_corpus.tsv:44",
    "category": "alcohol_amine"
  },
  {
    "name": "4-ethenylpyridine",
    "smiles": "C=CC1=CC=NC=C1",
    "goldenLine": 62,
    "source": "bakeoff_corpus.tsv:63",
    "category": "alkenyl_benzyl"
  },
  {
    "name": "3-phenylpropanoic acid",
    "smiles": "C1=CC=C(C=C1)CCC(=O)O",
    "goldenLine": 7,
    "source": "bakeoff_corpus.tsv:8",
    "category": "chain_ring_prefix"
  },
  {
    "name": "N-methyl-3-phenyl-3-[4-(trifluoromethyl)phenoxy]propan-1-amine",
    "smiles": "CNCCC(C1=CC=CC=C1)OC2=CC=C(C=C2)C(F)(F)F",
    "goldenLine": 98,
    "source": "bakeoff_corpus.tsv:99",
    "category": "drug"
  },
  {
    "name": "2-amino-1,7-dihydro-6H-purin-6-one",
    "smiles": "C1=NC2=C(N1)C(=O)NC(=N2)N",
    "goldenLine": 85,
    "source": "bakeoff_corpus.tsv:86",
    "category": "heterocycle"
  },
  {
    "name": "1,1'-biphenyl",
    "smiles": "C1=CC=C(C=C1)C2=CC=CC=C2",
    "goldenLine": 187,
    "source": "bakeoff_corpus.tsv:188",
    "category": "mr_assembly"
  },
  {
    "name": "7-oxabicyclo[2.2.1]heptane",
    "smiles": "C1CC2CCC1O2",
    "goldenLine": 161,
    "source": "bakeoff_corpus.tsv:162",
    "category": "mr_bridged"
  },
  {
    "name": "anthracene",
    "smiles": "C1=CC=C2C=C3C=CC=CC3=CC2=C1",
    "goldenLine": 171,
    "source": "bakeoff_corpus.tsv:172",
    "category": "mr_fused"
  },
  {
    "name": "1,4,7,10,13-pentaoxacyclopentadecane",
    "smiles": "C1COCCOCCOCCOCCO1",
    "goldenLine": 200,
    "source": "bakeoff_corpus.tsv:201",
    "category": "mr_macro_other"
  },
  {
    "name": "spiro[5.5]undecane",
    "smiles": "C1CCC2(CC1)CCCCC2",
    "goldenLine": 147,
    "source": "bakeoff_corpus.tsv:148",
    "category": "mr_spiro"
  },
  {
    "name": "2-hydroxypropane-1,2,3-tricarboxylic acid",
    "smiles": "C(C(=O)O)C(CC(=O)O)(C(=O)O)O",
    "goldenLine": 127,
    "source": "bakeoff_corpus.tsv:128",
    "category": "natural_product"
  },
  {
    "name": "2,2,6-trimethyl-6-phenylmorpholine",
    "smiles": "CC1(CNCC(O1)(C)C2=CC=CC=C2)C",
    "goldenLine": 506,
    "source": "bakeoff_corpus.tsv:507",
    "category": "pubchem_random"
  },
  {
    "name": "benzyl acetate",
    "smiles": "CC(=O)OCC1=CC=CC=C1",
    "goldenLine": 13,
    "source": "bakeoff_corpus.tsv:14",
    "category": "ring_ON_ester_amide"
  },
  {
    "name": "(2E)-but-2-enedioic acid",
    "smiles": "C(=C/C(=O)O)\\C(=O)O",
    "goldenLine": 133,
    "source": "bakeoff_corpus.tsv:134",
    "category": "stereo"
  },
  {
    "name": "benzoic acid",
    "smiles": "C1=CC=C(C=C1)C(=O)O",
    "goldenLine": 23,
    "source": "bakeoff_corpus.tsv:24",
    "category": "suffix_on_ring"
  },
  {
    "name": "furo[3,2-g]quinazoline",
    "smiles": "c1oc2cc3ncncc3cc2c1",
    "goldenLine": 719,
    "source": "openclatura_tests.tsv:17",
    "category": "data/fusion_pin_cases.json"
  },
  {
    "name": "2H-cyclopenta[b]furan",
    "smiles": "O1C=2C(=CC1)C=CC2",
    "goldenLine": 762,
    "source": "openclatura_tests.tsv:60",
    "category": "data/systematic_fusion_cases.json"
  },
  {
    "name": "hydrogen chloride",
    "smiles": "[ClH]",
    "goldenLine": 871,
    "source": "openclatura_tests.tsv:172",
    "category": "exploratory/test_charge_representation_gaps.py"
  },
  {
    "name": "pyridine",
    "smiles": "c1ccncc1",
    "goldenLine": 821,
    "source": "openclatura_tests.tsv:120",
    "category": "fixtures/diverse_corpus.golden.json"
  },
  {
    "name": "2,3,3a,6,7,8,9,10,11,11a-decahydro-1H-cyclopenta[10]annulen-5-ol",
    "smiles": "OC1=CC2CCCC2CCCCCC1",
    "goldenLine": 879,
    "source": "openclatura_tests.tsv:180",
    "category": "integration/test_fusion_annulene_production.py"
  },
  {
    "name": "1H-pyrrolo[3,2-c]pyridazin-2-ium-1-ide",
    "smiles": "[N-]1[NH+]=CC=C2N=CC=C12",
    "goldenLine": 894,
    "source": "openclatura_tests.tsv:195",
    "category": "integration/test_fusion_charged_migration.py"
  },
  {
    "name": "decahydropyrido[3,4-b]pyridin-7-ium",
    "smiles": "N1CCCC2C1C[NH2+]CC2",
    "goldenLine": 912,
    "source": "openclatura_tests.tsv:213",
    "category": "integration/test_fusion_saturated_composition.py"
  },
  {
    "name": "1,5-dihydropyrazolo[4,3-c]pyrazole",
    "smiles": "N1N=CC=2C1=CNN2",
    "goldenLine": 939,
    "source": "openclatura_tests.tsv:240",
    "category": "integration/test_systematic_fusion.py"
  },
  {
    "name": "acetaldehyde",
    "smiles": "CC=O",
    "goldenLine": 954,
    "source": "openclatura_tests.tsv:255",
    "category": "roundtrip/test_acids.py"
  },
  {
    "name": "1,2-oxazole",
    "smiles": "c1cnoc1",
    "goldenLine": 971,
    "source": "openclatura_tests.tsv:272",
    "category": "roundtrip/test_bridges.py"
  },
  {
    "name": "pentane",
    "smiles": "CCCCC",
    "goldenLine": 980,
    "source": "openclatura_tests.tsv:284",
    "category": "roundtrip/test_composite_bridges.py"
  },
  {
    "name": "1,3-dihydro-2H-indol-2-one",
    "smiles": "O=C1Cc2ccccc2N1",
    "goldenLine": 1019,
    "source": "openclatura_tests.tsv:323",
    "category": "roundtrip/test_fused_hydrocarbon_retained.py"
  },
  {
    "name": "quinoxaline",
    "smiles": "c1ccc2nccnc2c1",
    "goldenLine": 1051,
    "source": "openclatura_tests.tsv:356",
    "category": "roundtrip/test_fused_parents.py"
  },
  {
    "name": "2,3-dihydropyrrolo[2,1-b]quinazolin-9(1H)-one",
    "smiles": "O=c1c2ccccc2nc2n1CCC2",
    "goldenLine": 1102,
    "source": "openclatura_tests.tsv:408",
    "category": "roundtrip/test_fusion_hydrogenation_regressions.py"
  },
  {
    "name": "2-fluoropropane",
    "smiles": "CC(F)C",
    "goldenLine": 1113,
    "source": "openclatura_tests.tsv:419",
    "category": "roundtrip/test_haloalkanes.py"
  },
  {
    "name": "2-(methylsilyl)ethan-1-ol",
    "smiles": "C[SiH2]CCO",
    "goldenLine": 1118,
    "source": "openclatura_tests.tsv:424",
    "category": "roundtrip/test_heteroatom_substituents.py"
  },
  {
    "name": "imidazolidine",
    "smiles": "N1CNCC1",
    "goldenLine": 1139,
    "source": "openclatura_tests.tsv:447",
    "category": "roundtrip/test_heterocycles.py"
  },
  {
    "name": "N-methylmethanimine",
    "smiles": "C=NC",
    "goldenLine": 1159,
    "source": "openclatura_tests.tsv:467",
    "category": "roundtrip/test_imines.py"
  },
  {
    "name": "bis(chloromethyl)di(methyl)silane",
    "smiles": "C[Si](C)(CCl)CCl",
    "goldenLine": 1170,
    "source": "openclatura_tests.tsv:478",
    "category": "roundtrip/test_ligand_parentheses.py"
  },
  {
    "name": "nitroethane",
    "smiles": "CC[N+](=O)[O-]",
    "goldenLine": 1190,
    "source": "openclatura_tests.tsv:498",
    "category": "roundtrip/test_nitro.py"
  },
  {
    "name": "7-ethyl-2-methylundecan-4-yl hydrogen sulfate",
    "smiles": "CCCCC(CC)CCC(CC(C)C)OS(=O)(=O)O",
    "goldenLine": 1202,
    "source": "openclatura_tests.tsv:510",
    "category": "roundtrip/test_oxyacids.py"
  },
  {
    "name": "2,3-dihydrofuro[2,3-b]pyridine",
    "smiles": "C1COc2ncccc21",
    "goldenLine": 1205,
    "source": "openclatura_tests.tsv:513",
    "category": "roundtrip/test_p13_additive.py"
  },
  {
    "name": "1-phenylpropan-2-ol",
    "smiles": "c1ccccc1CC(C)O",
    "goldenLine": 1219,
    "source": "openclatura_tests.tsv:527",
    "category": "roundtrip/test_p13_conjunctive.py"
  },
  {
    "name": "2-azabicyclo[2.2.1]heptane",
    "smiles": "N1CC2CCC1C2",
    "goldenLine": 1229,
    "source": "openclatura_tests.tsv:537",
    "category": "roundtrip/test_p13_replacement.py"
  },
  {
    "name": "acetaldehyde hydrazone",
    "smiles": "CC=NN",
    "goldenLine": 1234,
    "source": "openclatura_tests.tsv:542",
    "category": "roundtrip/test_p13_substitutive.py"
  },
  {
    "name": "propyne",
    "smiles": "C#CC",
    "goldenLine": 1238,
    "source": "openclatura_tests.tsv:546",
    "category": "roundtrip/test_p13_subtractive.py"
  },
  {
    "name": "3-methylspiro[5.5]undecane",
    "smiles": "CC1CCC2(CC1)CCCCC2",
    "goldenLine": 1245,
    "source": "openclatura_tests.tsv:554",
    "category": "roundtrip/test_polycycles.py"
  },
  {
    "name": "1lambda^4-thiolan-1-one",
    "smiles": "O=S1CCCC1",
    "goldenLine": 1248,
    "source": "openclatura_tests.tsv:559",
    "category": "roundtrip/test_polycyclic_descriptors.py"
  },
  {
    "name": "3-methyl-N'-(4-propoxypyrimidin-2-yl)pentane-1,3-diamine",
    "smiles": "CCCOc1ccnc(NC(C)(CC)CCN)n1",
    "goldenLine": 1262,
    "source": "openclatura_tests.tsv:573",
    "category": "roundtrip/test_polyfunctional_fragments.py"
  },
  {
    "name": "propan-2-one hydrazone",
    "smiles": "CC(=NN)C",
    "goldenLine": 1276,
    "source": "openclatura_tests.tsv:587",
    "category": "roundtrip/test_principal_group_breadth.py"
  },
  {
    "name": "1,2-dimethyldisilyne",
    "smiles": "C[Si]#[Si]C",
    "goldenLine": 1312,
    "source": "openclatura_tests.tsv:623",
    "category": "roundtrip/test_replacement_parents.py"
  },
  {
    "name": "1H-purine-2,6-dione",
    "smiles": "N1C(N=C2N=CN=C2C1=O)=O",
    "goldenLine": 1318,
    "source": "openclatura_tests.tsv:629",
    "category": "roundtrip/test_retained_fused_production.py"
  },
  {
    "name": "4-hydrazinyl-1H-pyrazole",
    "smiles": "NNc1cn[nH]c1",
    "goldenLine": 1321,
    "source": "openclatura_tests.tsv:632",
    "category": "roundtrip/test_role_models.py"
  },
  {
    "name": "propanedithioic acid",
    "smiles": "CCC(=S)S",
    "goldenLine": 1335,
    "source": "openclatura_tests.tsv:647",
    "category": "unit/test_all_chalcogen_acyl_regressions.py"
  },
  {
    "name": "(anthracen-9-yl)acetic acid",
    "smiles": "O=C(O)Cc1c2ccccc2cc2ccccc12",
    "goldenLine": 1348,
    "source": "openclatura_tests.tsv:660",
    "category": "unit/test_analysis.py"
  },
  {
    "name": "(1r,4r)-4-aminocyclohexan-1-ol",
    "smiles": "N[C@H]1CC[C@H](O)CC1",
    "goldenLine": 1357,
    "source": "openclatura_tests.tsv:671",
    "category": "unit/test_audit_reconstruction.py"
  },
  {
    "name": "ethyl 3-(butoxycarbonothioyl)propanoate",
    "smiles": "CCCCOC(=S)CCC(=O)OCC",
    "goldenLine": 1372,
    "source": "openclatura_tests.tsv:686",
    "category": "unit/test_chalcogen_citation_compatibility.py"
  },
  {
    "name": "methanetelluronic acid",
    "smiles": "C[Te](=O)(=O)O",
    "goldenLine": 1448,
    "source": "openclatura_tests.tsv:762",
    "category": "unit/test_chalcogen_functional_groups.py"
  },
  {
    "name": "acetyl isocyanide",
    "smiles": "CC(=O)[N+]#[C-]",
    "goldenLine": 1511,
    "source": "openclatura_tests.tsv:825",
    "category": "unit/test_chalcogen_structural_roles.py"
  },
  {
    "name": "3-chloro-6,7-dihydro-5H-cyclopenta[b]pyridin-5-one",
    "smiles": "C1Cc2ncc(Cl)cc2C1=O",
    "goldenLine": 1520,
    "source": "openclatura_tests.tsv:834",
    "category": "unit/test_fused_oxo_added_hydrogen_tier.py"
  },
  {
    "name": "2-methylnaphtho[2,1-d][1,3]oxazol-5-ol",
    "smiles": "Cc1nc2cc(O)c3ccccc3c2o1",
    "goldenLine": 1524,
    "source": "openclatura_tests.tsv:838",
    "category": "unit/test_fusion_candidate_budget.py"
  },
  {
    "name": "6,7-dihydro-5H-cyclopenta[c]pyridin-5-one",
    "smiles": "C1Cc2cnccc2C1=O",
    "goldenLine": 1528,
    "source": "openclatura_tests.tsv:842",
    "category": "unit/test_fusion_carbon_h_oxo_composition.py"
  },
  {
    "name": "4H-[1,3]dioxolo[4,5-b]pyrrol-2-one",
    "smiles": "O=C1OC2=C(O1)C=CN2",
    "goldenLine": 1534,
    "source": "openclatura_tests.tsv:848",
    "category": "unit/test_fusion_intrinsic_h_suffix.py"
  },
  {
    "name": "1,4-epoxynaphthalen-2(1H)-one",
    "smiles": "C12C(C=C(C3=CC=CC=C13)O2)=O",
    "goldenLine": 1098,
    "source": "openclatura_tests.tsv:404",
    "category": "unit/test_fusion_wrappers.py"
  },
  {
    "name": "(Z)-N-(methanesulfonyl)benzaldehyde hydrazone",
    "smiles": "CS(=O)(=O)N/N=C\\c1ccccc1",
    "goldenLine": 1551,
    "source": "openclatura_tests.tsv:865",
    "category": "unit/test_hydrazone_stereo_scope.py"
  },
  {
    "name": "tetrasulfane",
    "smiles": "SSSS",
    "goldenLine": 1307,
    "source": "openclatura_tests.tsv:618",
    "category": "unit/test_polyoxide_chains.py"
  },
  {
    "name": "2-fluoro-7-methoxy-9H-fluorene",
    "smiles": "COc1ccc2c(c1)Cc1cc(F)ccc12",
    "goldenLine": 1572,
    "source": "openclatura_tests.tsv:886",
    "category": "unit/test_retained_parent_precedence_and_order.py"
  },
  {
    "name": "but-2-ene",
    "smiles": "CC=CC",
    "goldenLine": 863,
    "source": "rules_cases.tsv:9",
    "category": "rules"
  },
  {
    "name": "2-methylpropan-2-ol",
    "smiles": "CC(C)(C)O",
    "goldenLine": 55,
    "source": "bakeoff_corpus.tsv:56",
    "category": "alcohol_amine"
  },
  {
    "name": "2-benzylbutanoic acid",
    "smiles": "CCC(Cc1ccccc1)C(=O)O",
    "goldenLine": 58,
    "source": "bakeoff_corpus.tsv:59",
    "category": "alkenyl_benzyl"
  },
  {
    "name": "2-phenylpropanoic acid",
    "smiles": "CC(C(=O)O)c1ccccc1",
    "goldenLine": 11,
    "source": "bakeoff_corpus.tsv:12",
    "category": "chain_ring_prefix"
  },
  {
    "name": "2-(N,N-dimethylcarbamimidoyl)guanidine",
    "smiles": "CN(C)C(=N)N=C(N)N",
    "goldenLine": 95,
    "source": "bakeoff_corpus.tsv:96",
    "category": "drug"
  },
  {
    "name": "oxane",
    "smiles": "C1CCOCC1",
    "goldenLine": 68,
    "source": "bakeoff_corpus.tsv:69",
    "category": "heterocycle"
  },
  {
    "name": "1^1,2^1:2^3,3^1-terphenyl",
    "smiles": "C1=CC=C(C=C1)C2=CC(=CC=C2)C3=CC=CC=C3",
    "goldenLine": 191,
    "source": "bakeoff_corpus.tsv:192",
    "category": "mr_assembly"
  },
  {
    "name": "bicyclo[1.1.0]butane",
    "smiles": "C1C2C1C2",
    "goldenLine": 166,
    "source": "bakeoff_corpus.tsv:167",
    "category": "mr_bridged"
  },
  {
    "name": "(4as,8as)-decahydronaphthalene",
    "smiles": "C1CC[C@H]2CCCC[C@H]2C1",
    "goldenLine": 179,
    "source": "bakeoff_corpus.tsv:180",
    "category": "mr_fused"
  },
  {
    "name": "cyclohexadecane",
    "smiles": "C1CCCCCCCCCCCCCCC1",
    "goldenLine": 202,
    "source": "bakeoff_corpus.tsv:203",
    "category": "mr_macro_other"
  },
  {
    "name": "1-oxaspiro[4.5]decane",
    "smiles": "C1CCC2(CC1)CCCO2",
    "goldenLine": 151,
    "source": "bakeoff_corpus.tsv:152",
    "category": "mr_spiro"
  },
  {
    "name": "(6E)-N-[(4-hydroxy-3-methoxyphenyl)methyl]-8-methylnon-6-enamide",
    "smiles": "CC(C)/C=C/CCCCC(=O)NCC1=CC(=C(C=C1)O)OC",
    "goldenLine": 121,
    "source": "bakeoff_corpus.tsv:122",
    "category": "natural_product"
  },
  {
    "name": "2-methyloctyl (2R)-2-[(7-chlorodibenzofuran-2-yl)oxy]propanoate",
    "smiles": "CCCCCCC(C)COC(=O)[C@@H](C)OC1=CC2=C(C=C1)OC3=C2C=CC(=C3)Cl",
    "goldenLine": 579,
    "source": "bakeoff_corpus.tsv:580",
    "category": "pubchem_random"
  },
  {
    "name": "N-phenylbenzamide",
    "smiles": "C1=CC=C(C=C1)C(=O)NC2=CC=CC=C2",
    "goldenLine": 20,
    "source": "bakeoff_corpus.tsv:21",
    "category": "ring_ON_ester_amide"
  },
  {
    "name": "(2R,3S)-2,3-dihydroxybutanedioic acid",
    "smiles": "O=C(O)[C@@H](O)[C@@H](O)C(=O)O",
    "goldenLine": 141,
    "source": "bakeoff_corpus.tsv:142",
    "category": "stereo"
  },
  {
    "name": "toluene",
    "smiles": "CC1=CC=CC=C1",
    "goldenLine": 40,
    "source": "bakeoff_corpus.tsv:41",
    "category": "suffix_on_ring"
  },
  {
    "name": "1H-pyrrolo[2,3-c]pyridine",
    "smiles": "N1C=CC=2C1=CN=CC2",
    "goldenLine": 751,
    "source": "openclatura_tests.tsv:49",
    "category": "data/fusion_pin_cases.json"
  },
  {
    "name": "difuro[3,2-b:2',3'-e]pyridine",
    "smiles": "O1C=CC2=C1C=C1C(=N2)C=CO1",
    "goldenLine": 765,
    "source": "openclatura_tests.tsv:63",
    "category": "data/systematic_fusion_cases.json"
  },
  {
    "name": "bromoethane",
    "smiles": "BrCC",
    "goldenLine": 828,
    "source": "openclatura_tests.tsv:127",
    "category": "fixtures/diverse_corpus.golden.json"
  },
  {
    "name": "dodecahydrocyclopenta[9]annulen-5-ol",
    "smiles": "OC1CCCCCC2CCCC2C1",
    "goldenLine": 875,
    "source": "openclatura_tests.tsv:176",
    "category": "integration/test_fusion_annulene_production.py"
  },
  {
    "name": "3-oxo-1,4,5,6-tetrahydropyrrolo[3,4-b]pyrrol-5-ium-2-ide",
    "smiles": "O=C1[CH-]NC2=C1C[NH2+]C2",
    "goldenLine": 900,
    "source": "openclatura_tests.tsv:201",
    "category": "integration/test_fusion_charged_migration.py"
  },
  {
    "name": "1-ethyldecahydropyrido[4,3-b]pyridin-6-ium",
    "smiles": "CCN1CCCC2C1CC[NH2+]C2",
    "goldenLine": 924,
    "source": "openclatura_tests.tsv:225",
    "category": "integration/test_fusion_saturated_composition.py"
  },
  {
    "name": "tetrahydro-1H,3H-thieno[3,4-c]thiophene",
    "smiles": "C1SCC2CSCC12",
    "goldenLine": 937,
    "source": "openclatura_tests.tsv:238",
    "category": "integration/test_systematic_fusion.py"
  },
  {
    "name": "3,6-dihydro-2H-pyran",
    "smiles": "O1CC=CCC1",
    "goldenLine": 977,
    "source": "openclatura_tests.tsv:281",
    "category": "roundtrip/test_bridges.py"
  },
  {
    "name": "1,3-dihydro-2H-inden-2-one",
    "smiles": "O=C1Cc2ccccc2C1",
    "goldenLine": 1024,
    "source": "openclatura_tests.tsv:328",
    "category": "roundtrip/test_fused_hydrocarbon_retained.py"
  },
  {
    "name": "phthalazine",
    "smiles": "c1ccc2cnncc2c1",
    "goldenLine": 1052,
    "source": "openclatura_tests.tsv:357",
    "category": "roundtrip/test_fused_parents.py"
  },
  {
    "name": "6-methyl-6,7-dihydro-5H-[1,3]thiazolo[3,2-a]pyrimidine",
    "smiles": "CC1CN2C(=NC1)SC=C2",
    "goldenLine": 1105,
    "source": "openclatura_tests.tsv:411",
    "category": "roundtrip/test_fusion_hydrogenation_regressions.py"
  },
  {
    "name": "2-methoxyethan-1-ol",
    "smiles": "COCCO",
    "goldenLine": 1114,
    "source": "openclatura_tests.tsv:420",
    "category": "roundtrip/test_heteroatom_substituents.py"
  },
  {
    "name": "pyridine-4-carboxylic acid",
    "smiles": "O=C(O)c1ccncc1",
    "goldenLine": 1150,
    "source": "openclatura_tests.tsv:458",
    "category": "roundtrip/test_heterocycles.py"
  },
  {
    "name": "cyclohexanimine",
    "smiles": "C1CCCCC1=N",
    "goldenLine": 1162,
    "source": "openclatura_tests.tsv:470",
    "category": "roundtrip/test_imines.py"
  },
  {
    "name": "dimethylbis(2-methylpropyl)silane",
    "smiles": "C[Si](C)(CC(C)C)CC(C)C",
    "goldenLine": 1172,
    "source": "openclatura_tests.tsv:480",
    "category": "roundtrip/test_ligand_parentheses.py"
  },
  {
    "name": "nitrobenzene",
    "smiles": "O=[N+]([O-])c1ccccc1",
    "goldenLine": 1191,
    "source": "openclatura_tests.tsv:499",
    "category": "roundtrip/test_nitro.py"
  },
  {
    "name": "ethyl nitrite",
    "smiles": "CCON=O",
    "goldenLine": 1199,
    "source": "openclatura_tests.tsv:507",
    "category": "roundtrip/test_oxyacids.py"
  },
  {
    "name": "1-phenylethan-1-ol",
    "smiles": "c1ccccc1C(C)O",
    "goldenLine": 1218,
    "source": "openclatura_tests.tsv:526",
    "category": "roundtrip/test_p13_conjunctive.py"
  },
  {
    "name": "(aminooxy)methane",
    "smiles": "NOC",
    "goldenLine": 1227,
    "source": "openclatura_tests.tsv:535",
    "category": "roundtrip/test_p13_replacement.py"
  },
  {
    "name": "N-hydroxyethanimine",
    "smiles": "CC=NO",
    "goldenLine": 1233,
    "source": "openclatura_tests.tsv:541",
    "category": "roundtrip/test_p13_substitutive.py"
  },
  {
    "name": "decahydronaphthalene",
    "smiles": "C1CCC2CCCCC2C1",
    "goldenLine": 1242,
    "source": "openclatura_tests.tsv:550",
    "category": "roundtrip/test_polycycles.py"
  },
  {
    "name": "bicyclo[2.2.1]hept-2-ene",
    "smiles": "C1=CC2CCC1C2",
    "goldenLine": 1250,
    "source": "openclatura_tests.tsv:561",
    "category": "roundtrip/test_polycyclic_descriptors.py"
  },
  {
    "name": "1-[(3R)-1-cyclopentylpyrrolidin-3-yl]-3-ethyl-2-[2-(1-methyl-1H-pyrazol-4-yl)ethyl]guanidine",
    "smiles": "CCNC(=NCCc1cnn(C)c1)N[C@@H]1CCN(C2CCCC2)C1",
    "goldenLine": 1256,
    "source": "openclatura_tests.tsv:567",
    "category": "roundtrip/test_polyfunctional_fragments.py"
  },
  {
    "name": "N-hydroxypropan-2-imine",
    "smiles": "CC(=NO)C",
    "goldenLine": 1274,
    "source": "openclatura_tests.tsv:585",
    "category": "roundtrip/test_principal_group_breadth.py"
  },
  {
    "name": "pentafluoro-lambda5-phosphane",
    "smiles": "FP(F)(F)(F)F",
    "goldenLine": 1316,
    "source": "openclatura_tests.tsv:627",
    "category": "roundtrip/test_replacement_parents.py"
  },
  {
    "name": "3-hydroxy-2,2-dimethoxypropyl hydrogen carbonate",
    "smiles": "COC(CO)(COC(=O)O)OC",
    "goldenLine": 1324,
    "source": "openclatura_tests.tsv:635",
    "category": "roundtrip/test_role_models.py"
  },
  {
    "name": "propanethioic S-acid",
    "smiles": "CCC(=O)S",
    "goldenLine": 1337,
    "source": "openclatura_tests.tsv:649",
    "category": "unit/test_all_chalcogen_acyl_regressions.py"
  },
  {
    "name": "anthracen-1-ylacetic acid",
    "smiles": "O=C(O)Cc1cccc2cc3ccccc3cc12",
    "goldenLine": 1347,
    "source": "openclatura_tests.tsv:659",
    "category": "unit/test_analysis.py"
  },
  {
    "name": "3,7-dihydro-1H-purine-2,6-dione",
    "smiles": "O=c1[nH]c(=O)c2[nH]cnc2[nH]1",
    "goldenLine": 1362,
    "source": "openclatura_tests.tsv:676",
    "category": "unit/test_audit_reconstruction.py"
  },
  {
    "name": "O-methyl ethaneselenoate",
    "smiles": "CC(=[Se])OC",
    "goldenLine": 1375,
    "source": "openclatura_tests.tsv:689",
    "category": "unit/test_chalcogen_citation_compatibility.py"
  },
  {
    "name": "methyl ethaneperoxoate",
    "smiles": "CC(=O)OOC",
    "goldenLine": 1432,
    "source": "openclatura_tests.tsv:746",
    "category": "unit/test_chalcogen_functional_groups.py"
  },
  {
    "name": "4-cyano-5-tellanylidenepyrrolidine-2-selone",
    "smiles": "N#CC1CC(=[Se])NC1=[Te]",
    "goldenLine": 1506,
    "source": "openclatura_tests.tsv:820",
    "category": "unit/test_chalcogen_structural_roles.py"
  },
  {
    "name": "3H-imidazo[4,5-b]pyridine",
    "smiles": "c1cnc2[nH]cnc2c1",
    "goldenLine": 1521,
    "source": "openclatura_tests.tsv:835",
    "category": "unit/test_fusion_candidate_budget.py"
  },
  {
    "name": "furo[2,3-d][1,3]oxazol-2(3H)-one",
    "smiles": "O=C1NC2=C(O1)C=CO2",
    "goldenLine": 1532,
    "source": "openclatura_tests.tsv:846",
    "category": "unit/test_fusion_intrinsic_h_suffix.py"
  },
  {
    "name": "1,4-ethenonaphthalene",
    "smiles": "C12=CC=C(C3=CC=CC=C13)C=C2",
    "goldenLine": 1545,
    "source": "openclatura_tests.tsv:859",
    "category": "unit/test_fusion_wrappers.py"
  },
  {
    "name": "(Z)-N-(4-ethoxybenzene-1-sulfonyl)-2,3,4,5,6-pentafluorobenzaldehyde hydrazone",
    "smiles": "CCOc1ccc(S(=O)(=O)N/N=C\\c2c(F)c(F)c(F)c(F)c2F)cc1",
    "goldenLine": 1550,
    "source": "openclatura_tests.tsv:864",
    "category": "unit/test_hydrazone_stereo_scope.py"
  },
  {
    "name": "trisulfane",
    "smiles": "SSS",
    "goldenLine": 1570,
    "source": "openclatura_tests.tsv:884",
    "category": "unit/test_polyoxide_chains.py"
  },
  {
    "name": "2-fluoro-6-methoxyanthracene",
    "smiles": "COc1ccc2cc3cc(F)ccc3cc2c1",
    "goldenLine": 1577,
    "source": "openclatura_tests.tsv:891",
    "category": "unit/test_retained_parent_precedence_and_order.py"
  },
  {
    "name": "propane-1,2-diol",
    "smiles": "CC(O)CO",
    "goldenLine": 1597,
    "source": "rules_cases.tsv:26",
    "category": "rules"
  },
  {
    "name": "ethane-1,2-diol",
    "smiles": "C(CO)O",
    "goldenLine": 54,
    "source": "bakeoff_corpus.tsv:55",
    "category": "alcohol_amine"
  },
  {
    "name": "5-(prop-2-en-1-yl)-2H-1,3-benzodioxole",
    "smiles": "C=CCC1=CC2=C(C=C1)OCO2",
    "goldenLine": 64,
    "source": "bakeoff_corpus.tsv:65",
    "category": "alkenyl_benzyl"
  },
  {
    "name": "(2S)-2-amino-3-(4-hydroxyphenyl)propanoic acid",
    "smiles": "C1=CC(=CC=C1C[C@@H](C(=O)O)N)O",
    "goldenLine": 6,
    "source": "bakeoff_corpus.tsv:7",
    "category": "chain_ring_prefix"
  },
  {
    "name": "4-hydroxy-3-(3-oxo-1-phenylbutyl)-2H-1-benzopyran-2-one",
    "smiles": "CC(=O)CC(C1=CC=CC=C1)C2=C(C3=CC=CC=C3OC2=O)O",
    "goldenLine": 102,
    "source": "bakeoff_corpus.tsv:103",
    "category": "drug"
  },
  {
    "name": "2-methyloxane",
    "smiles": "CC1CCCCO1",
    "goldenLine": 70,
    "source": "bakeoff_corpus.tsv:71",
    "category": "heterocycle"
  },
  {
    "name": "4-phenylpyridine",
    "smiles": "C1=CC=C(C=C1)C2=CC=NC=C2",
    "goldenLine": 193,
    "source": "bakeoff_corpus.tsv:194",
    "category": "mr_assembly"
  },
  {
    "name": "adamantane",
    "smiles": "C1C2CC3CC1CC(C2)C3",
    "goldenLine": 157,
    "source": "bakeoff_corpus.tsv:158",
    "category": "mr_bridged"
  },
  {
    "name": "naphthalene-2-carboxylic acid",
    "smiles": "C1=CC=C2C=C(C=CC2=C1)C(=O)O",
    "goldenLine": 170,
    "source": "bakeoff_corpus.tsv:171",
    "category": "mr_fused"
  },
  {
    "name": "1,4,7,10,13,16-hexaoxacyclooctadecane",
    "smiles": "C1COCCOCCOCCOCCOCCO1",
    "goldenLine": 199,
    "source": "bakeoff_corpus.tsv:200",
    "category": "mr_macro_other"
  },
  {
    "name": "spiro[2.2]pentane",
    "smiles": "C1CC12CC2",
    "goldenLine": 148,
    "source": "bakeoff_corpus.tsv:149",
    "category": "mr_spiro"
  },
  {
    "name": "(4R)-1-methyl-4-(prop-1-en-2-yl)cyclohex-1-ene",
    "smiles": "CC1=CC[C@@H](CC1)C(=C)C",
    "goldenLine": 118,
    "source": "bakeoff_corpus.tsv:119",
    "category": "natural_product"
  },
  {
    "name": "2-(4-bromophenoxy)-N-[1-(1H-pyrazol-1-yl)propan-2-yl]acetamide",
    "smiles": "CC(CN1C=CC=N1)NC(=O)COC2=CC=C(C=C2)Br",
    "goldenLine": 539,
    "source": "bakeoff_corpus.tsv:540",
    "category": "pubchem_random"
  },
  {
    "name": "N-(4-hydroxyphenyl)acetamide",
    "smiles": "CC(=O)NC1=CC=C(C=C1)O",
    "goldenLine": 16,
    "source": "bakeoff_corpus.tsv:17",
    "category": "ring_ON_ester_amide"
  },
  {
    "name": "(2S,3R)-2-amino-3-hydroxybutanoic acid",
    "smiles": "C[C@H]([C@@H](C(=O)O)N)O",
    "goldenLine": 136,
    "source": "bakeoff_corpus.tsv:137",
    "category": "stereo"
  },
  {
    "name": "pyridine-3-carboxamide",
    "smiles": "C1=CC(=CN=C1)C(=O)N",
    "goldenLine": 29,
    "source": "bakeoff_corpus.tsv:30",
    "category": "suffix_on_ring"
  },
  {
    "name": "furo[3,2-g]cinnoline",
    "smiles": "c1oc2cc3nnccc3cc2c1",
    "goldenLine": 724,
    "source": "openclatura_tests.tsv:22",
    "category": "data/fusion_pin_cases.json"
  },
  {
    "name": "thieno[2,3-b]furan",
    "smiles": "O1C2=C(C=C1)C=CS2",
    "goldenLine": 761,
    "source": "openclatura_tests.tsv:59",
    "category": "data/systematic_fusion_cases.json"
  },
  {
    "name": "chlorobenzene",
    "smiles": "C1=CC=CC=C1Cl",
    "goldenLine": 865,
    "source": "openclatura_tests.tsv:166",
    "category": "fixtures/diverse_corpus.golden.json"
  },
  {
    "name": "5-methyldodecahydrocyclopenta[9]annulene",
    "smiles": "CC1CCCCCC2CCCC2C1",
    "goldenLine": 874,
    "source": "openclatura_tests.tsv:175",
    "category": "integration/test_fusion_annulene_production.py"
  },
  {
    "name": "1H-imidazo[4,5-c]pyridazin-2-ium-1-ide",
    "smiles": "[N-]1[NH+]=CC=C2N=CN=C12",
    "goldenLine": 895,
    "source": "openclatura_tests.tsv:196",
    "category": "integration/test_fusion_charged_migration.py"
  },
  {
    "name": "1-methyldecahydropyrido[4,3-b]pyridin-6-ium",
    "smiles": "CN1CCCC2C1CC[NH2+]C2",
    "goldenLine": 923,
    "source": "openclatura_tests.tsv:224",
    "category": "integration/test_fusion_saturated_composition.py"
  },
  {
    "name": "tetrahydro-1H,3H-furo[3,4-c]furan",
    "smiles": "C1OCC2COCC12",
    "goldenLine": 936,
    "source": "openclatura_tests.tsv:237",
    "category": "integration/test_systematic_fusion.py"
  },
  {
    "name": "cycloheptane",
    "smiles": "C1CCCCCC1",
    "goldenLine": 978,
    "source": "openclatura_tests.tsv:282",
    "category": "roundtrip/test_bridges.py"
  },
  {
    "name": "2,3-dihydro-1H-indene",
    "smiles": "c1ccc2c(c1)CCC2",
    "goldenLine": 1021,
    "source": "openclatura_tests.tsv:325",
    "category": "roundtrip/test_fused_hydrocarbon_retained.py"
  },
  {
    "name": "phenanthren-2-ol",
    "smiles": "Oc1ccc2c(c1)ccc1ccccc12",
    "goldenLine": 1072,
    "source": "openclatura_tests.tsv:378",
    "category": "roundtrip/test_fused_parents.py"
  },
  {
    "name": "pyrrolo[2,1-b]quinazolin-9(1H)-one",
    "smiles": "O=c1c2ccccc2nc2n1CC=C2",
    "goldenLine": 1106,
    "source": "openclatura_tests.tsv:412",
    "category": "roundtrip/test_fusion_hydrogenation_regressions.py"
  },
  {
    "name": "oxetane",
    "smiles": "C1COC1",
    "goldenLine": 1123,
    "source": "openclatura_tests.tsv:429",
    "category": "roundtrip/test_heterocycles.py"
  }
];
