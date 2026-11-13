// SPDX-License-Identifier: MIT AND (MIT OR Apache-2.0)
// Origin: mixed
// Source: openclatura | src/openclatura/assembly_utils.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Source: ReShiki | src/naming/openclatura/assembly_utils.rs | 5f559fd482e237307a589dfcb7ea333cee0a6e13
// Source: openclatura | tests/unit/test_analysis.py | 90dc5bcb4910c196122ec0405d06a32929a8f308
// Copyright (c) 2026 lamalab-org
// Copyright (c) 2026 Ameyanagi and ReShiki contributors
// Changes: Openclatura grammar helpers with local compound-locant and stereo-hyphen fixes.
// Local contributions: MIT OR Apache-2.0; see NOTICE and PROVENANCE-nazuchem.toml.

//! Small grammar utilities shared by the plain-name formatters.
// Provenance: fn parse_locant | adapted
// Source: openclatura | src/openclatura/assembly_utils.py::parse_locant
pub fn parse_locant(text: &str) -> (u8, usize, String) {
    let before = text.split('(').next().unwrap_or(text);
    let digit_count = before.bytes().take_while(u8::is_ascii_digit).count();
    if let Some((primary, secondary)) = before.split_once('^') {
        // An italic element locant such as N^4 precedes numerals, like N; its
        // superscript is a parent locant, compared as one (P-62.2.4.1.1).
        if primary.starts_with(|c: char| c.is_ascii_uppercase()) {
            // A primed nitrogen follows its unprimed twin, N^3 then N'^3, and
            // both precede the nitrogens on the primed parent locant, N^3'
            // (P-14.3.5). The space sorts the nitrogen's primes after the
            // parent locant's own primes and letters.
            let letter = primary.trim_end_matches('\'');
            let primes = &primary[letter.len()..];
            let (_, number, tail) = parse_locant(secondary);
            return (0, 0, format!("{letter}^{number:020}{tail} {primes}"));
        }
        let n = secondary.bytes().take_while(u8::is_ascii_digit).count();
        if let (Ok(primary), Ok(secondary_number)) =
            (primary.parse::<usize>(), secondary[..n].parse::<usize>())
        {
            // P-14.3.1: compare composite locants by primary locant, then superscript.
            return (
                1,
                primary,
                format!("^{secondary_number:020}{}", &secondary[n..]),
            );
        }
    }
    if digit_count > 0
        && before[digit_count..]
            .chars()
            .all(|c| c.is_ascii_alphabetic() || c == '\'')
    {
        return (
            1,
            before[..digit_count].parse().unwrap_or(0),
            format!(
                "{}{}",
                before[digit_count..].replace('\'', ""),
                "'".repeat(before[digit_count..].matches('\'').count())
            ),
        );
    }
    let digits = text.trim_start_matches(|c: char| !c.is_ascii_digit());
    let digit_count = digits.bytes().take_while(u8::is_ascii_digit).count();
    if digit_count > 0 {
        return (
            1,
            digits[..digit_count].parse().unwrap_or(0),
            text.to_owned(),
        );
    }
    (0, 0, text.to_owned())
}
// Provenance: fn needs_hyphen | adapted
// Source: openclatura | src/openclatura/assembly_utils.py::needs_hyphen
pub fn needs_hyphen(left: &str, right: &str) -> bool {
    if left.is_empty() || right.is_empty() {
        return false;
    }
    // P-16.2.4.1(d): a hyphen sets off the italic s or as that opens the
    // retained parents s- and as-indacene (P-25.1.1): 3-chloro-s-indacene.
    if right.starts_with(|c: char| c.is_ascii_digit())
        || left.ends_with(|c: char| c.is_ascii_digit())
        || right.starts_with("s-")
        || right.starts_with("as-")
    {
        return true;
    }
    if right
        .strip_prefix(['N', 'O', 'S', 'P'])
        .is_some_and(|tail| {
            // N', N^4 or N^4a: a superscript parent locant may carry a fusion letter.
            tail.trim_start_matches(|c: char| {
                matches!(c, '\'' | '^') || c.is_ascii_digit() || c.is_ascii_lowercase()
            })
            .starts_with([',', '-'])
        })
    {
        return true;
    }
    // A stereo descriptor such as `(1R,2E)` opening the right-hand side.
    right
        .strip_prefix('(')
        .and_then(|s| s.split_once(')'))
        .is_some_and(|(inner, _)| {
            inner.split(',').all(|term| {
                let letters = term.trim_start_matches(|c: char| c.is_ascii_digit());
                letters.len() < term.len()
                    && letters.ends_with(['R', 'S', 'E', 'Z', 'M', 'P'])
                    && letters
                        .chars()
                        .all(|c| c.is_ascii_alphabetic() || c == '\'')
            })
        })
}
// Provenance: fn is_fully_enclosed | port
// Source: openclatura | src/openclatura/assembly_utils.py::is_fully_enclosed
pub fn is_fully_enclosed(text: &str) -> bool {
    if !text.starts_with('(') || !text.ends_with(')') {
        return false;
    }
    let mut depth = 0_i32;
    for (i, c) in text.char_indices() {
        match c {
            '(' => depth += 1,
            ')' => depth -= 1,
            _ => (),
        }
        if depth == 0 && i < text.len() - 1 {
            return false;
        }
        if depth < 0 {
            return false;
        }
    }
    depth == 0
}

#[cfg(test)]
// Provenance: mod tests | original
// Source: ReShiki | src/naming/openclatura/assembly_utils.rs::tests
mod tests {
    use super::{needs_hyphen, parse_locant};

    #[test]
    fn assembly_locants_compare_primes_and_superscripts_numerically() {
        let mut primed = ["2'", "1'", "2", "1"];
        primed.sort_by_key(|locant| parse_locant(locant));
        assert_eq!(primed, ["1", "1'", "2", "2'"]);
        let mut fusion = ["3'a", "3a", "3'", "3"];
        fusion.sort_by_key(|locant| parse_locant(locant));
        assert_eq!(fusion, ["3", "3'", "3a", "3'a"]);
        let mut composite = ["2^1", "1^10", "1^2", "2^2"];
        composite.sort_by_key(|locant| parse_locant(locant));
        assert_eq!(composite, ["1^2", "1^10", "2^1", "2^2"]);
    }

    #[test]
    fn superscript_n_locants_precede_numerals_in_parent_locant_order() {
        // P-62.2.4.1.1, P-62.2.4.1.2: N^4,2-dimethylpentane-2,4-diamine.
        let mut locants = ["2", "N^12", "N^1^10", "N^4a", "N^4", "N^1^2"];
        locants.sort_by_key(|locant| parse_locant(locant));
        assert_eq!(locants, ["N^1^2", "N^1^10", "N^4", "N^4a", "N^12", "2"]);
        assert!(needs_hyphen("ethyl", "N^4a-methyl"));
        assert!(needs_hyphen("ethyl", "N^1^10,N^1^2-dimethyl"));
        // P-14.3.5: a primed locant follows its unprimed twin, so N^1,N'^3 as in
        // N^3-ethyl-N^1,N'^3-dimethylhexane-1,3,3,6-tetramine (P-62.2.4.1.2).
        let mut locants = ["N'^3", "N^4a", "N^3", "N''^3", "N^1", "N'^4"];
        locants.sort_by_key(|locant| parse_locant(locant));
        assert_eq!(locants, ["N^1", "N^3", "N'^3", "N''^3", "N'^4", "N^4a"]);
        assert!(needs_hyphen("ethyl", "N'^3-methyl"));
        // P-14.3.5: 4 is lower than 4', and 4' than 4a, so every nitrogen on
        // locant 4 precedes those on 4' (issue #107).
        let mut locants = ["N^4a", "N'^4'", "N^4'", "N''^4", "N'^4", "N^4"];
        locants.sort_by_key(|locant| parse_locant(locant));
        assert_eq!(locants, ["N^4", "N'^4", "N''^4", "N^4'", "N'^4'", "N^4a"]);
    }

    #[test]
    // Provenance: fn tests::only_a_leading_stereo_descriptor_takes_a_hyphen | adapted
    // Source: ReShiki | src/naming/openclatura/assembly_utils.rs::tests::only_a_leading_stereo_descriptor_takes_a_hyphen
    // Source: openclatura | tests/unit/test_analysis.py::test_distinct_center_ligands_only_group_ligands_after_the_first
    // Test provenance: upstream inputs/expected names with local assertions and extensions;
    // exact matches are conservatively adapted; independent ancestry is not established.
    fn only_a_leading_stereo_descriptor_takes_a_hyphen() {
        assert!(needs_hyphen("methyl", "(1R,2E)-but-2-enyl"));
        assert!(needs_hyphen("methyl", "(5'E)-1'-butyl"));
        assert!(needs_hyphen("methyl", "(1M)-biphenyl-2-yl"));
        assert!(needs_hyphen("methyl", "(1'P,2R)-biphenyl-2-yl"));
        assert!(!needs_hyphen(
            "(cyclopentyl)",
            "(4-(((3R)-thiolan-3-yl)sulfanyl)pyridin-3-yl)"
        ));
        assert!(!needs_hyphen("methyl", "SO-x"));
    }

    #[test]
    fn italic_indacene_descriptors_keep_their_hyphen() {
        // P-16.2.4.1(d): a hyphen separates the italic s or as of s- and
        // as-indacene (P-25.1.1) from a substituent or hydro prefix.
        for (smiles, expected) in [
            ("ClC1=c2cc3c(cc2C=C1)=CC=C3", "3-chloro-s-indacene"),
            ("ClC=1C=CC2=C3C=CC=C3C=CC12", "3-chloro-as-indacene"),
            ("ClCC1=CC=C2C=C3C=CC=C3C=C12", "1-(chloromethyl)-s-indacene"),
            (
                "C1CC2=CC3=C(CCC3)C=C2C1",
                "1,2,3,5,6,7-hexahydro-s-indacene",
            ),
            ("C1CCC2=CC=C3CCCC3=C12", "1,2,3,6,7,8-hexahydro-as-indacene"),
            ("C1CCC2CC3CCCC3CC12", "dodecahydro-s-indacene"),
            // Added indicated hydrogen cited after the suffix (P-14.7.2).
            (
                "O=C1CC2=CC3=C(CCC3)C=C2C1",
                "3,5,6,7-tetrahydro-s-indacen-2(1H)-one",
            ),
            (
                "O=C1CCC2=CC=C3CCCC3=C12",
                "2,3,7,8-tetrahydro-as-indacen-1(6H)-one",
            ),
        ] {
            assert_eq!(
                super::super::name_smiles(smiles, &Default::default())
                    .map(|name| name.text)
                    .as_deref(),
                Ok(expected),
                "{smiles}"
            );
            // The parser reads the hyphenated spelling back to the structure.
            let parsed = nazuchem_n2s::Parser::new()
                .unwrap()
                .parse(expected, &nazuchem_n2s::ParseOptions::strict());
            let round_trip = parsed.structure.unwrap().semantic_cxsmiles().unwrap();
            crate::test_support::verify_identity(smiles, &round_trip).unwrap();
        }
    }
}
