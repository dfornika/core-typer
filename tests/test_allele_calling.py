import csv

from core_typer import allele_calling


def make_kma_result(allele_id, score, template_identity=100.0, query_identity=100.0,
                     template_coverage=100.0, depth=20.0, allele_hash=None):
    template = f"locusA_{allele_id}" if allele_hash is None else f"locusA_{allele_id}_{allele_hash}"
    return {
        "template": template,
        "locus_id": "locusA",
        "allele_id": allele_id,
        "allele_hash": allele_hash,
        "score": score,
        "template_length": 500,
        "template_identity": template_identity,
        "query_identity": query_identity,
        "template_coverage": template_coverage,
        "depth": depth,
    }


def test_choose_best_allele_picks_highest_score():
    kma_results = [
        make_kma_result("1", score=480),
        make_kma_result("2", score=500),
    ]

    best = allele_calling.choose_best_allele(kma_results)

    assert best["allele_id"] == "2"
    assert best["closest_allele_id"] == "2"
    assert best["call_status"] == "called"
    assert best["num_hits"] == 2


def test_choose_best_allele_below_identity_threshold_is_uncalled():
    kma_results = [make_kma_result("1", score=500, query_identity=95.0)]

    best = allele_calling.choose_best_allele(kma_results, min_identity=100.0)

    assert best["allele_id"] == "-"
    assert best["closest_allele_id"] == "1"
    assert best["call_status"] == "divergent"


def test_choose_best_allele_below_coverage_threshold_is_uncalled():
    kma_results = [make_kma_result("1", score=500, template_coverage=90.0)]

    best = allele_calling.choose_best_allele(kma_results, min_coverage=100.0)

    assert best["allele_id"] == "-"
    assert best["closest_allele_id"] == "1"
    assert best["call_status"] == "divergent"


def test_choose_best_allele_uses_query_identity_not_template_identity():
    # Regression test for the CGE-cgMLSTFinder-derived fix: template_identity
    # conflates mismatches with missing coverage into one number, so a locus
    # with partial coverage but a perfect match within the aligned region can
    # have a low template_identity despite query_identity being 100. Calling
    # must gate on query_identity (+ template_coverage separately), not
    # template_identity.
    kma_results = [make_kma_result("1", score=500, template_identity=95.0, query_identity=100.0, template_coverage=100.0)]

    best = allele_calling.choose_best_allele(kma_results, min_identity=100.0, min_coverage=100.0)

    assert best["call_status"] == "called"
    assert best["allele_id"] == "1"


def test_choose_best_allele_highest_score_wins_even_if_uncalled():
    # Matches CGE-style allele calling: the single best-scoring KMA hit is
    # chosen for the locus; if it fails thresholds, the locus is uncalled
    # even if a lower-scoring hit would have passed.
    kma_results = [
        make_kma_result("1", score=500, query_identity=90.0),
        make_kma_result("2", score=400, query_identity=100.0),
    ]

    best = allele_calling.choose_best_allele(kma_results, min_identity=100.0)

    assert best["allele_id"] == "-"
    assert best["closest_allele_id"] == "1"
    assert best["call_status"] == "divergent"
    assert best["score"] == 500


def test_choose_best_allele_propagates_closest_allele_hash():
    kma_results = [make_kma_result("1", score=500, allele_hash="5d41402abc4b2a76b9719d911017c592")]

    best = allele_calling.choose_best_allele(kma_results)

    assert best["closest_allele_hash"] == "5d41402abc4b2a76b9719d911017c592"
    assert best["closest_allele_template"] == "locusA_1_5d41402abc4b2a76b9719d911017c592"
    assert best["novel_allele_hash"] is None


def test_build_complete_allele_calls_fills_missing_loci():
    parsed_kma_result = {
        "locusA": [make_kma_result("1", score=500)],
    }
    locus_ids = ["locusA", "locusB"]

    allele_calls = allele_calling.build_complete_allele_calls(locus_ids, parsed_kma_result)

    assert [c["locus_id"] for c in allele_calls] == ["locusA", "locusB"]
    assert allele_calls[0]["allele_id"] == "1"
    assert allele_calls[0]["call_status"] == "called"
    assert allele_calls[0]["num_hits"] == 1
    assert allele_calls[1]["allele_id"] == "-"
    assert allele_calls[1]["closest_allele_id"] is None
    assert allele_calls[1]["closest_allele_hash"] is None
    assert allele_calls[1]["call_status"] == "no_hit"
    assert allele_calls[1]["num_hits"] == 0
    assert allele_calls[1]["depth"] is None


def test_find_possible_multicopy_loci_flags_loci_with_multiple_hits():
    parsed_kma_result = {
        "locusA": [make_kma_result("1", score=500), make_kma_result("2", score=300, query_identity=97.0)],
        "locusB": [make_kma_result("1", score=400)],
    }

    multicopy_records = allele_calling.find_possible_multicopy_loci(parsed_kma_result)

    assert [r["locus_id"] for r in multicopy_records] == ["locusA", "locusA"]
    assert [r["allele_id"] for r in multicopy_records] == ["1", "2"]
    assert multicopy_records[0]["score"] == 500
    assert multicopy_records[1]["query_identity"] == 97.0


def test_find_possible_multicopy_loci_empty_when_no_locus_has_multiple_hits():
    parsed_kma_result = {
        "locusA": [make_kma_result("1", score=500)],
        "locusB": [make_kma_result("1", score=400)],
    }

    multicopy_records = allele_calling.find_possible_multicopy_loci(parsed_kma_result)

    assert multicopy_records == []


def test_write_possible_multicopy_loci(tmp_path):
    multicopy_records = [
        {"locus_id": "locusA", "allele_id": "1", "query_identity": 100.0, "template_coverage": 100.0, "depth": 30.0, "score": 500},
        {"locus_id": "locusA", "allele_id": "2", "query_identity": 97.0, "template_coverage": 100.0, "depth": 28.0, "score": 300},
    ]
    multicopy_path = tmp_path / "possible_multicopy_loci.csv"

    allele_calling.write_possible_multicopy_loci(multicopy_records, str(multicopy_path))

    with open(multicopy_path) as f:
        rows = list(csv.DictReader(f))

    assert len(rows) == 2
    assert rows[0]["locus_id"] == "locusA"
    assert rows[0]["allele_id"] == "1"
    assert rows[1]["allele_id"] == "2"
    assert rows[1]["query_identity"] == "97.0"


def test_is_clean_sequence():
    assert allele_calling.is_clean_sequence("ACGTACGT") is True
    assert allele_calling.is_clean_sequence("ACGTNACGT") is False
    assert allele_calling.is_clean_sequence("acgtacgt") is False


def test_hash_sequence_matches_known_md5():
    assert allele_calling.hash_sequence("hello") == "5d41402abc4b2a76b9719d911017c592"


def test_extract_novel_alleles_finds_divergent_full_coverage_clean_consensus():
    allele_calls = [
        allele_calling.choose_best_allele(
            [make_kma_result("1", score=500, query_identity=98.0, template_coverage=100.0)],
            min_identity=100.0,
        ),
    ]
    consensus_by_template = {"locusA_1": "ACGTACGT"}

    novel_alleles = allele_calling.extract_novel_alleles(allele_calls, consensus_by_template)

    assert len(novel_alleles) == 1
    assert novel_alleles[0]["locus_id"] == "locusA"
    assert novel_alleles[0]["sequence"] == "ACGTACGT"
    assert novel_alleles[0]["hash"] == allele_calling.hash_sequence("ACGTACGT")
    assert allele_calls[0]["novel_allele_hash"] == novel_alleles[0]["hash"]


def test_extract_novel_alleles_looks_up_consensus_by_full_template_name():
    # Regression test: for schemes prepared with allele hashing, kma's -ef
    # fsa output is keyed by the *full* template name including the hash
    # suffix (e.g. "locusA_1_<md5>"), not a "<locus_id>_<allele_id>"
    # reconstruction - the lookup must use closest_allele_template as-is.
    allele_calls = [
        allele_calling.choose_best_allele(
            [make_kma_result("1", score=500, query_identity=98.0, template_coverage=100.0,
                              allele_hash="5d41402abc4b2a76b9719d911017c592")],
            min_identity=100.0,
        ),
    ]
    consensus_by_template = {"locusA_1_5d41402abc4b2a76b9719d911017c592": "ACGTACGT"}

    novel_alleles = allele_calling.extract_novel_alleles(allele_calls, consensus_by_template)

    assert len(novel_alleles) == 1
    assert novel_alleles[0]["sequence"] == "ACGTACGT"


def test_extract_novel_alleles_skips_called_loci():
    allele_calls = [allele_calling.choose_best_allele([make_kma_result("1", score=500)])]
    consensus_by_template = {"locusA_1": "ACGTACGT"}

    novel_alleles = allele_calling.extract_novel_alleles(allele_calls, consensus_by_template)

    assert novel_alleles == []
    assert allele_calls[0]["novel_allele_hash"] is None


def test_extract_novel_alleles_skips_partial_coverage():
    allele_calls = [
        allele_calling.choose_best_allele(
            [make_kma_result("1", score=500, query_identity=98.0, template_coverage=90.0)],
            min_identity=100.0,
        ),
    ]
    consensus_by_template = {"locusA_1": "ACGTACGT"}

    novel_alleles = allele_calling.extract_novel_alleles(allele_calls, consensus_by_template)

    assert novel_alleles == []


def test_extract_novel_alleles_skips_ambiguous_consensus():
    allele_calls = [
        allele_calling.choose_best_allele(
            [make_kma_result("1", score=500, query_identity=98.0, template_coverage=100.0)],
            min_identity=100.0,
        ),
    ]
    consensus_by_template = {"locusA_1": "ACGTNACGT"}

    novel_alleles = allele_calling.extract_novel_alleles(allele_calls, consensus_by_template)

    assert novel_alleles == []


def test_extract_novel_alleles_skips_missing_consensus():
    allele_calls = [
        allele_calling.choose_best_allele(
            [make_kma_result("1", score=500, query_identity=98.0, template_coverage=100.0)],
            min_identity=100.0,
        ),
    ]

    novel_alleles = allele_calling.extract_novel_alleles(allele_calls, consensus_by_template={})

    assert novel_alleles == []


def test_write_novel_alleles(tmp_path):
    novel_alleles = [{"locus_id": "locusA", "hash": "abc123", "sequence": "ACGT"}]
    novel_alleles_path = tmp_path / "novel_alleles.fsa"

    allele_calling.write_novel_alleles(novel_alleles, str(novel_alleles_path))

    content = novel_alleles_path.read_text()
    assert content == ">locusA novel_allele_hash=abc123\nACGT\n"


def test_write_allele_profile_two_rows_matching_order(tmp_path):
    allele_calls = [
        {"locus_id": "locusA", "allele_id": "1"},
        {"locus_id": "locusB", "allele_id": "-"},
    ]
    profile_path = tmp_path / "allele_profile.csv"

    allele_calling.write_allele_profile(allele_calls, str(profile_path))

    with open(profile_path) as f:
        rows = list(csv.reader(f))

    assert rows == [["locusA", "locusB"], ["1", "-"]]


def test_write_allele_calls_includes_all_records(tmp_path):
    allele_calls = [
        allele_calling.choose_best_allele([make_kma_result("1", score=500, allele_hash="5d41402abc4b2a76b9719d911017c592")]),
        {
            "locus_id": "locusB", "allele_id": "-", "closest_allele_id": None, "closest_allele_hash": None,
            "novel_allele_hash": None, "num_hits": 0, "call_status": "no_hit", "score": None, "template_length": None,
            "template_identity": None, "query_identity": None, "template_coverage": None, "depth": None,
        },
    ]
    calls_path = tmp_path / "allele_calls.csv"

    allele_calling.write_allele_calls(str(calls_path), allele_calls)

    with open(calls_path) as f:
        rows = list(csv.DictReader(f))

    assert len(rows) == 2
    assert rows[0]["locus_id"] == "locusA"
    assert rows[0]["call_status"] == "called"
    assert rows[0]["closest_allele_id"] == "1"
    assert rows[0]["closest_allele_hash"] == "5d41402abc4b2a76b9719d911017c592"
    assert rows[0]["query_identity"] == "100.0"
    assert rows[0]["num_hits"] == "1"
    assert rows[1]["locus_id"] == "locusB"
    assert rows[1]["allele_id"] == "-"
    assert rows[1]["call_status"] == "no_hit"
    assert rows[1]["closest_allele_id"] == ""
    assert rows[1]["num_hits"] == "0"
