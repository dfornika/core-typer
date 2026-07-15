import os

from core_typer import parsers

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


def test_parse_kma_result_groups_by_locus():
    result = parsers.parse_kma_result(os.path.join(DATA_DIR, "kma-out.res"))

    assert set(result.keys()) == {"locusA", "locusB"}
    assert len(result["locusA"]) == 2
    assert len(result["locusB"]) == 1


def test_parse_kma_result_sorted_by_score_descending():
    result = parsers.parse_kma_result(os.path.join(DATA_DIR, "kma-out.res"))

    scores = [row["score"] for row in result["locusA"]]
    assert scores == sorted(scores, reverse=True)
    assert result["locusA"][0]["allele_id"] == "1"


def test_parse_kma_result_field_types():
    result = parsers.parse_kma_result(os.path.join(DATA_DIR, "kma-out.res"))

    row = result["locusA"][0]
    assert row["score"] == 500
    assert row["template_length"] == 500
    assert row["template_identity"] == 100.00
    assert row["template_coverage"] == 100.00
    assert row["depth"] == 20.5
    assert row["locus_id"] == "locusA"
    assert row["allele_id"] == "1"
    assert row["allele_hash"] is None


def test_parse_template_name_without_hash():
    locus_id, allele_id, allele_hash = parsers.parse_template_name("locusA_1")

    assert locus_id == "locusA"
    assert allele_id == "1"
    assert allele_hash is None


def test_parse_template_name_with_hash():
    locus_id, allele_id, allele_hash = parsers.parse_template_name(
        "locusA_1_5d41402abc4b2a76b9719d911017c592"
    )

    assert locus_id == "locusA"
    assert allele_id == "1"
    assert allele_hash == "5d41402abc4b2a76b9719d911017c592"


def test_parse_template_name_locus_id_with_underscore():
    # locus ids can themselves contain underscores (e.g. Enterobase-style
    # "STMMW_00001") - the trailing allele id/hash must still parse correctly.
    locus_id, allele_id, allele_hash = parsers.parse_template_name(
        "STMMW_00001_5_098f6bcd4621d373cade4e832627b4f6"
    )

    assert locus_id == "STMMW_00001"
    assert allele_id == "5"
    assert allele_hash == "098f6bcd4621d373cade4e832627b4f6"


def test_parse_kma_result_extracts_allele_hash_when_present():
    result = parsers.parse_kma_result(os.path.join(DATA_DIR, "kma-out-hashed.res"))

    assert result["locusA"][0]["allele_hash"] == "5d41402abc4b2a76b9719d911017c592"
    assert result["STMMW_00001"][0]["locus_id"] == "STMMW_00001"
    assert result["STMMW_00001"][0]["allele_id"] == "5"
    assert result["STMMW_00001"][0]["allele_hash"] == "098f6bcd4621d373cade4e832627b4f6"


def test_parse_kma_consensus_fasta():
    consensus = parsers.parse_kma_consensus_fasta(os.path.join(DATA_DIR, "kma-out.fsa"))

    assert consensus == {
        "locusA_1": "ACGTACGTACGTACGT",
        "locusB_1": "TTTTGGGGCCCC",
    }


def test_parse_kma_mapstat_groups_by_locus():
    result = parsers.parse_kma_mapstat(os.path.join(DATA_DIR, "kma-out.mapstat"))

    assert set(result.keys()) == {"locusA", "locusB"}
    row = result["locusA"][0]
    assert row["read_count"] == 1000
    assert row["depth_variance"] == 1.2
    assert row["locus_id"] == "locusA"
    assert row["allele_id"] == "1"


def test_parse_locus_names_dedups_and_preserves_order():
    locus_names = parsers.parse_locus_names(os.path.join(DATA_DIR, "scheme.name"))

    assert locus_names == ["locusA", "locusB", "locusC"]


def test_parse_blast_result_groups_by_locus():
    result = parsers.parse_blast_result(os.path.join(DATA_DIR, "blast-out.tsv"))

    assert set(result.keys()) == {"locusA", "locusB"}


def test_parse_blast_result_dedupes_overlapping_hits_to_best_scoring():
    # locusA_1 and locusA_2 both match the *same* query region (contig1:1-300
    # vs 1-297) - two different catalogued alleles matching one genomic
    # copy, not two copies. Only the better-scoring one should survive.
    result = parsers.parse_blast_result(os.path.join(DATA_DIR, "blast-out.tsv"))

    assert len(result["locusA"]) == 1
    assert result["locusA"][0]["allele_id"] == "1"


def test_parse_blast_result_keeps_hits_at_distinct_query_regions():
    # locusB_1 (401-700) and locusB_2 (701-1000) are non-overlapping query
    # regions - genuine evidence of two separate copies - both should survive.
    result = parsers.parse_blast_result(os.path.join(DATA_DIR, "blast-out.tsv"))

    assert len(result["locusB"]) == 2
    assert {r["allele_id"] for r in result["locusB"]} == {"1", "2"}


def test_parse_blast_result_field_values():
    result = parsers.parse_blast_result(os.path.join(DATA_DIR, "blast-out.tsv"))

    row = result["locusA"][0]
    assert row["locus_id"] == "locusA"
    assert row["allele_id"] == "1"
    assert row["template_length"] == 300
    assert row["template_identity"] == 100.0
    assert row["template_coverage"] == 100.0
    assert row["query_identity"] == 100.0
    assert row["score"] == 555.0
    assert row["depth"] is None


def test_parse_locus_names_from_fasta_dedups_and_preserves_order():
    locus_names = parsers.parse_locus_names_from_fasta(os.path.join(DATA_DIR, "scheme.fasta"))

    assert locus_names == ["locusA", "locusB"]


def test_parse_allele_calls_treats_blank_numeric_fields_as_none(tmp_path):
    calls_path = tmp_path / "allele_calls.csv"
    calls_path.write_text(
        "locus_id,call_status,allele_id,closest_allele_id,closest_allele_hash,novel_allele_hash,"
        "query_identity,template_identity,template_coverage,depth,score,template_length\n"
        "locusA,called,1,1,,,100.0,100.0,100.0,20.5,500,500\n"
        "locusB,no_hit,-,,,,,,,,,\n"
    )

    rows = parsers.parse_allele_calls(str(calls_path))

    assert rows[0]["score"] == 500
    assert rows[0]["depth"] == 20.5
    assert rows[1]["score"] is None
    assert rows[1]["depth"] is None
    assert rows[1]["template_length"] is None
