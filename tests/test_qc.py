import csv

from core_typer import qc


def make_call(locus_id, allele_id, depth, call_status, num_hits=1):
    return {"locus_id": locus_id, "allele_id": allele_id, "depth": depth, "call_status": call_status, "num_hits": num_hits}


def test_calculate_qc_stats_all_called():
    allele_calls = [
        make_call("locusA", "1", 20.0, "called"),
        make_call("locusB", "1", 18.0, "called"),
    ]

    stats = qc.calculate_qc_stats(allele_calls)

    assert stats["num_loci"] == 2
    assert stats["num_called_alleles"] == 2
    assert stats["num_divergent"] == 0
    assert stats["num_no_hit"] == 0
    assert stats["percent_called"] == 100.0
    assert stats["mean_depth"] == 19.0
    assert stats["stdev_depth"] is not None


def test_calculate_qc_stats_counts_missing_loci_in_denominator():
    # Regression test: a locus with zero KMA hits must still count toward
    # num_loci / percent_called, not just loci that got at least one hit.
    allele_calls = [
        make_call("locusA", "1", 20.0, "called"),
        make_call("locusB", "-", None, "no_hit"),
        make_call("locusC", "-", None, "no_hit"),
    ]

    stats = qc.calculate_qc_stats(allele_calls)

    assert stats["num_loci"] == 3
    assert stats["num_called_alleles"] == 1
    assert stats["percent_called"] == round(1 / 3 * 100, 3)
    assert stats["mean_depth"] == 20.0


def test_calculate_qc_stats_distinguishes_divergent_from_no_hit():
    allele_calls = [
        make_call("locusA", "1", 20.0, "called"),
        make_call("locusB", "-", 5.0, "divergent"),
        make_call("locusC", "-", None, "no_hit"),
    ]

    stats = qc.calculate_qc_stats(allele_calls)

    assert stats["num_called_alleles"] == 1
    assert stats["num_divergent"] == 1
    assert stats["num_no_hit"] == 1
    # mean_depth only counts loci with an actual depth measurement (called or
    # divergent), not "no_hit" loci where nothing mapped at all.
    assert stats["mean_depth"] == round((20.0 + 5.0) / 2, 3)


def test_calculate_qc_stats_counts_possible_multicopy_loci():
    allele_calls = [
        make_call("locusA", "1", 20.0, "called", num_hits=1),
        make_call("locusB", "1", 18.0, "called", num_hits=2),
        make_call("locusC", "-", 5.0, "divergent", num_hits=3),
    ]

    stats = qc.calculate_qc_stats(allele_calls)

    assert stats["num_possible_multicopy_loci"] == 2


def test_calculate_qc_stats_single_locus_no_stdev_error():
    allele_calls = [make_call("locusA", "1", 20.0, "called")]

    stats = qc.calculate_qc_stats(allele_calls)

    assert stats["mean_depth"] == 20.0
    assert stats["stdev_depth"] is None


def test_calculate_qc_stats_empty_scheme():
    stats = qc.calculate_qc_stats([])

    assert stats["num_loci"] == 0
    assert stats["percent_called"] == 0.0
    assert stats["mean_depth"] is None
    assert stats["stdev_depth"] is None


def test_write_qc_stats(tmp_path):
    stats = qc.calculate_qc_stats([make_call("locusA", "1", 20.0, "called"), make_call("locusB", "1", 22.0, "called")])
    qc_path = tmp_path / "qc.csv"

    qc.write_qc_stats(stats, str(qc_path))

    with open(qc_path) as f:
        rows = list(csv.DictReader(f))

    assert len(rows) == 1
    assert rows[0]["num_loci"] == "2"
    assert rows[0]["num_divergent"] == "0"
    assert rows[0]["num_no_hit"] == "0"
    assert rows[0]["num_possible_multicopy_loci"] == "0"
    assert rows[0]["percent_called"] == "100.0"
